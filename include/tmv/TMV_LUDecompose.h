///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef TMV_LUDecompose_H
#define TMV_LUDecompose_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_MultMV.h"
#include "TMV_MultMM.h"
#include "TMV_DivVU.h"
#include "TMV_DivMU.h"
#include "TMV_SwapV.h"
#include "TMV_Permutation.h"
#include "TMV_PermuteM.h"

#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#include "TMV_TriMatrixIO.h"

#include "TMV_Matrix.h"
#include "TMV_SmallMatrix.h"

#ifdef PRINTALGO_LU
#include <iostream>
#endif

// TMV_LU_BLOCKSIZE is the block size to use in algo 21
#define TMV_LU_BLOCKSIZE 64

// Q2 is the maximum size to stop recursing in algo 27
#define TMV_Q2 1

namespace tmv {

    // Defined below:
    template <class M> 
    void LU_Decompose(
        BaseMatrix_Rec_Mutable<M>& m, int* P, int& detp);
    template <class M> 
    void InlineLU_Decompose(
        BaseMatrix_Rec_Mutable<M>& m, int* P, int& detp);

    // Defined in TMV_LUDecompose.cpp
    template <class T> 
    void InstLU_Decompose(MatrixView<T,1> m, int* P, int& detp);

    template <int algo, int cs, int rs, class M>
    struct LUDecompose_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or 1, or N == 0)
    template <int cs, int rs, class M>
    struct LUDecompose_Helper<0,cs,rs,M>
    { static void call(M& A, int* P, int& detp) {} };

    // algo 1: N == 1
    template <int cs, class M1>
    struct LUDecompose_Helper<1,cs,1,M1>
    {
        static void call(M1& A, int* P, int& detp)
        {
            const int M = cs==UNKNOWN ? int(A.colsize()) : cs;
#ifdef PRINTALGO_LU
            std::cout<<"LUDecompose algo 1: M,N,cs,rs = "<<M<<','<<1<<
                ','<<cs<<','<<1<<std::endl;
#endif
            // Same as NonBlock version, but with R==1 hard coded
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            typename M1::col_type A0 = A.get_col(0);
            RT piv = A0.maxAbsElement(P);

            if (TMV_Underflow(piv)) {
                *P = 0;
                A0.setZero();
            } else {
                if (*P != 0) {
                    A0.swap(*P,0);
                    detp = -detp;
                }
                //A0.cSubVector(1,M) /= A0.cref(0);
                typename M1::col_sub_type A0b = A0.cSubVector(1,M);
                InlineScale(Scaling<0,T>(RT(1)/A0.cref(0)),A0b);
            }
        }
    };

    // algo 2: N == 2
    template <int cs, class M1>
    struct LUDecompose_Helper<2,cs,2,M1>
    {
        static void call(M1& A, int* P, int& detp)
        {
            const int M = cs==UNKNOWN ? int(A.colsize()) : cs;
#ifdef PRINTALGO_LU
            std::cout<<"LUDecompose algo 2: M,N,cs,rs = "<<M<<','<<2<<
                ','<<cs<<','<<2<<std::endl;
#endif
            // Same as NonBlock version, but with N==2 hard coded
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            typename M1::col_type A0 = A.get_col(0);
            typename M1::col_type A1 = A.get_col(1);
            typename M1::col_type::iterator it0 = A0.begin();
            typename M1::col_type::iterator it1 = A1.begin();

            int ip0,ip1;
            RT piv = InlineMaxAbsElement(A0,&ip0);

            if (TMV_Underflow(piv)) {
                A0.setZero();
                ip0 = 0;
                piv = InlineMaxAbsElement(A1.cSubVector(1,M),&ip1);
                ip1++;
            } else {
                if (ip0 != 0) {
                    A0.swap(ip0,0);
                    A1.swap(ip0,0);
                    detp = -detp;
                } 

                // A0.subVector(1,M) /= A00;
                // A1.subVector(1,M) -= A0.subVector(1,M) * A01;
                const T invA00 = RT(1)/(*it0++);
                const T A01 = *it1++;
                piv = RT(0); // next pivot element
                ip1 = 1;
                for(int i=1;i<M;++i,++it0,++it1) {
                    *it0 *= invA00;
                    *it1 -= *it0 * A01;
                    RT absAi1 = TMV_ABS(*it1);
                    if (absAi1 > piv) { piv = absAi1; ip1=i; }
                }
            }

            if (TMV_Underflow(piv)) {
                A1.cSubVector(1,M).setZero();
                ip1 = 1;
            } else {
                if (M > 2) {
                    if (ip1 != 1) {
                        A1.swap(ip1,1);
                        A0.swap(ip1,1);
                        detp = -detp;
                    } 

                    //A1.cSubVector(2,M) /= A1.cref(1);
                    typename M1::col_sub_type A1b = A1.cSubVector(2,M);
                    InlineScale(Scaling<0,T>(RT(1)/A1.cref(1)),A1b);
                }
            }

            P[0] = ip0;
            P[1] = ip1;
        }
    };

    // algo 3: M == 2
    template <int rs, class M1>
    struct LUDecompose_Helper<3,2,rs,M1>
    {
        static void call(M1& A, int* P, int& detp)
        {
            const int N = rs==UNKNOWN ? int(A.rowsize()) : rs;
#ifdef PRINTALGO_LU
            std::cout<<"LUDecompose algo 3: M,N,cs,rs = "<<2<<','<<N<<
                ','<<2<<','<<rs<<std::endl;
#endif
            // Same as NonBlock version, but with M==2 hard coded
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            typename M1::col_type A0 = A.get_col(0);
            typename M1::col_type A1 = A.get_col(1);

            int ip0;
            RT piv = InlineMaxAbsElement(A0,&ip0);

            if (TMV_Underflow(piv)) {
                A0.setZero();
                ip0 = 0;
            } else {
                if (ip0 != 0) {
                    A0.swap(ip0,0);
                    A1.swap(ip0,0);
                    detp = -detp;
                } 

                // A0.subVector(1,M) /= A00;
                // A1.subVector(1,M) -= A0.subVector(1,M) * A01;
                A0.ref(1) /= A0.cref(0);
                A1.ref(1) -= A0.cref(1) * A1.cref(0);
            }

            P[0] = ip0;
            P[1] = 1;

            if (N <= 2) return;
            else {
                typename M1::row_sub_type Aa = A.get_row(0,2,N);
                typename M1::row_sub_type Ab = A.get_row(1,2,N);
                // M=2, N>2, so solve for U(0:2,2:N))
                // A.colRange(2,N).permuteRows(P);
                if (ip0 == 1) {
                    //A.cColRange(2,N).swapRows(0,1);
                    InlineSwap(Aa,Ab);
                }

                // A.colRange(2,N) /= A.colRange(0,2).unitLowerTri();
                //A.get_row(1,2,N) -= A1.cref(0) * A.get_row(0,2,N);
                InlineMultXV<true>(Scaling<0,T>(-A1.cref(0)),Aa,Ab);
            }
        }
    };

    // algo 11: Non-block algorithm, loop over n
    template <int cs, int rs, class M1>
    struct LUDecompose_Helper<11,cs,rs,M1>
    {
        static void call(M1& A, int* P, int& detp)
        {
            // LU Decompostion with partial pivoting.
            //
            // We want to decompose the matrix (input as A) into P * L * U
            // where P is a permutation, 
            // L is a unit-diag lower triangle matrix,
            // and U is an upper triangle matrix.  
            //
            // After doing j cols of the calculation, we have calculated
            // U(0:j,0:j) and L(0:N,0:j) 
            // (where my a:b notation does not include the index b)
            //
            // The equation A = LU gives for the j col:
            //
            // A(0:N,j) = L(0:N,0:N) U(0:N,j)
            //
            // which breaks up into:
            //
            // (1) A(0:j,j) = L(0:j,0:N) U(0:N,j)
            // (2) A(j:N,j) = L(j:N,0:N) U(0:N,j)
            //
            // The first of these (1) simplifies to:
            // 
            // (1*) A(0:j,j) = L(0:j,0:j) U(0:j,j)
            //
            // since L is lower triangular, so L(0:j,j:N) = 0.
            // L(0:j,0:j) is already known, so this equation can be solved for
            // U(0:j,j) by forward substitution.
            // 
            // The second equation (2) simplifies to:
            //
            //      A(j:N,j) = L(j:N,0:j+1) U(0:j+1,j)
            // (2*)          = L(j:N,0:j) U(0:j,j) + L(j:N,j) U(j,j)
            //
            // since U is upper triangular so U(j+1:N,j) = 0.
            // Since we now know U(0:j,j) from (1*) above, this equation can
            // be solved for the product L(j:N,j) U(j,j)
            // 
            // This means we have some leeway on the values for L(j,j) and 
            // U(j,j), as only their product is specified.
            //
            // If we take U to have unit diagonal, then L(j,j) is set here along
            // with the rest of the L(j:N,j) column.  However, this will mean 
            // that the forward substutions in the (1*) steps will require 
            // divisions by the non-unit-diagonal elements of L.
            // It is faster to take L to have unit-diagonal elements,
            // and do the division by U(j,j) here, since then we can calculate 
            // 1/U(j,j) and multiply.  So 1 division and N-j multiplies 
            // which is generally faster than N-j divisions.
            //
            // However, another potential problem is that U(j,j) could be 0,
            // or close to 0.  
            // This would lead to either an error or inaccurate results.
            // Thus we add a step in the middle of the (2*) calculation:
            //
            // Define v(j:N) = A(j:N,j) - L(j:N,0:j) U(0:j,j)
            // 
            // We search v for the element with the largest absolute value and 
            // apply a permutation to swap it into the j spot.
            // This element then becomes U(j,j), which is then the divisor for 
            // the rest of the vector.  This will minimize the possibility of 
            // roundoff errors due to small U(j,j)'s.
            // 
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            const int N = rs==UNKNOWN ? int(A.rowsize()) : rs;
            const int M = cs==UNKNOWN ? int(A.colsize()) : cs;
            const int R = TMV_MIN(N,M);
            const int xx = UNKNOWN;
#ifdef PRINTALGO_LU
            std::cout<<"LUDecompose algo 11: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif

            typedef typename M1::col_sub_type M1c;
            typedef typename M1::row_type M1r;
            typedef typename M1::const_col_sub_type M1cc;
            typedef typename M1::colrange_type M1cr;
            typedef typename M1::const_submatrix_type M1s;
            typedef typename M1::submatrix_type::const_unit_lowertri_type M1l;

            for (int j=0; j<R; ++j) {
                M1c Ajb = A.get_col(j,j,M);
                M1l L = A.cSubMatrix(0,j,0,j).unitLowerTri();

                if (j > 0) {
                    M1c Aja = A.get_col(j,0,j);
                    // Solve for U(0:j,j))
                    //A.get_col(j,0,j) /= A.cSubMatrix(0,j,0,j).unitLowerTri();
                    LDivEqVU_Helper<-4,xx,M1c,M1l>::call(Aja,L);

                    // Solve for v = L(j:M,j) U(j,j)
                    //A.get_col(j,j,M) -= 
                        //A.cSubMatrix(j,M,0,j) * A.get_col(j,0,j);
                    MultMV_Helper<-4,xx,xx,true,-1,RT,M1s,M1cc,M1c>::call(
                        Scaling<-1,RT>(),A.cSubMatrix(j,M,0,j),Aja,Ajb);
                }

                // Find the pivot element
                int ip;
                RT piv = InlineMaxAbsElement(Ajb,&ip);
                // ip is relative to j index, not absolute.

                // Check for underflow:
                if (TMV_Underflow(piv)) {
                    P[j] = j;
                    Ajb.setZero();
                    continue;
                }

                // Swap the pivot row with j if necessary
                if (ip != 0) {
                    ip += j;
                    TMVAssert(ip < int(A.colsize()));
                    TMVAssert(j < int(A.colsize()));
                    //A.cSwapRows(ip,j);  // This does both Lkb and A'
                    M1r Aip = A.get_row(ip);
                    M1r Aj = A.get_row(j);
                    InlineSwap(Aip,Aj);
                    P[j] = ip;
                    detp = -detp;
                } else P[j] = j;

                // Solve for L(j+1:M,j)
                // If Ujj is 0, then all of the L's are 0.
                // ie. Ujj Lij = 0 for all i>j
                // Any value for Lij is valid, so leave them 0.
                if (A.cref(j,j) != T(0)) {
                    //A.col(j,j+1,M) /= *Ujj;
                    M1c Ajc = A.get_col(j,j+1,M);
                    InlineScale(Scaling<0,T>(RT(1)/A.cref(j,j)),Ajc);
                }
            }
            if (N > M) {
                // Solve for U(0:M,M:N))
                //A.cColRange(M,N) /= A.cColRange(0,M).unitLowerTri();
                M1cr Acr = A.cColRange(M,N);
                M1l L = A.cColRange(0,M).unitLowerTri();
                LDivEqMU_Helper<-2,cs,xx,M1cr,M1l>::call(Acr,L);
            }
        }
    };

    // algo 21: Block algorithm
    template <int cs, int rs, class M1>
    struct LUDecompose_Helper<21,cs,rs,M1>
    {
        typedef typename M1::value_type T;
        static void call(M1& A, int* P, int& detp)
        {
            // If A is large, we can take advantage of Blas Level 3 speed
            // by partitioning the matrix by columns.
            //
            // We do this calculation one block at a time down the diagonal. 
            // Each block contains k columns (expect for possibly the last 
            // block which may be fewer).
            //
            // For the first block, we decompose A into:
            //
            // ( A00 A01 ) = ( L00  0  ) ( U00 U01 )
            // ( A10 A11 )   ( L10  A~ ) (  0   1  )
            //
            // From this we obtain:
            //
            // (1) A00 = L00 U00
            // (2) A10 = L10 U00
            // (3) A01 = L00 U01
            // (4) A11 = L10 U01 + A~
            //
            // For (1) we decompose A00 in place using the non-blocked 
            // algorithm.
            // (2) and (3) then give us L10 and U01.
            // Finally, (4) lets us solve for A~.
            // Repeat until done.
            //
            // With pivoting, the only real change is to combine equations 
            // (1),(2) and solve both together with the non-blocked algorithm.
            //

            typedef typename M1::real_type RT;

            const int N = rs==UNKNOWN ? int(A.rowsize()) : rs;
            const int M = cs==UNKNOWN ? int(A.colsize()) : cs;
            const int Nx = TMV_LU_BLOCKSIZE;
            const int R = TMV_MIN(N,M);
            const int xx = UNKNOWN;
#ifdef PRINTALGO_LU
            std::cout<<"LUDecompose algo 21: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif

            int jk=0;
            for (int jkpk=jk+Nx; jkpk<=R; jk=jkpk, jkpk+=Nx) {
                typedef typename M1::submatrix_type M1s;
                typedef typename M1::const_submatrix_type M1sc;
                typedef typename M1s::const_unit_lowertri_type M1l;
                M1s A0 = A.cSubMatrix(jk,M,jk,jkpk); // Both A00 and A10
                M1s A1 = A.cSubMatrix(jk,M,jkpk,N); // Both A01 and A11
                M1s A2 = A.cSubMatrix(jk,M,0,jk); // Previous blocks
                M1s A00 = A.cSubMatrix(jk,jkpk,jk,jkpk);
                M1s A10 = A.cSubMatrix(jkpk,M,jk,jkpk);
                M1s A01 = A.cSubMatrix(jk,jkpk,jkpk,N);
                M1s A11 = A.cSubMatrix(jkpk,M,jkpk,N);

                // Solve for L00, U00
                LUDecompose_Helper<11,xx,Nx,M1s>::call(A0,P+jk,detp);

                // Apply the permutation to the rest of the matrix
                if (jk > 0) {
                    //A2.cPermuteRows(P+jk,0,Nx);
                    InlinePermuteRows(A2,P+jk,0,Nx);
                }
                if (jkpk < N) {
                    //A1.cPermuteRows(P+jk,0,Nx);
                    InlinePermuteRows(A1,P+jk,0,Nx);

                    // Solve for U01
                    //A01 /= A00.unitLowerTri();
                    M1l L00 = A00.unitLowerTri();
                    LDivEqMU_Helper<-2,Nx,xx,M1s,M1l>::call(A01,L00);

                    // Solve for A~
                    if (jkpk < M) {
                        //A11 -= A10 * A01;
                        MultMM_Helper<-2,xx,xx,Nx,true,-1,RT,M1sc,M1sc,M1s>::call(
                            Scaling<-1,RT>(),A10,A01,A11);
                    }
                }
                for(int i=jk;i<jkpk;++i) P[i]+=jk;
            }

            if (jk < R) { // Last block is not full size

                const int Ny = R-jk;
                typedef typename M1::submatrix_type M1s;
                typedef typename M1s::const_unit_lowertri_type M1l;
                M1s A0 = A.cSubMatrix(jk,M,jk,R); // Both A00 and A10
                M1s A1 = A.cSubMatrix(jk,M,R,N); // Both A01 and A11
                M1s A2 = A.cSubMatrix(jk,M,0,jk); // Previous blocks
                M1s A00 = A.cSubMatrix(jk,R,jk,R);
                M1s A01 = A.cSubMatrix(jk,R,R,N);

                // Solve for L00, U00
                LUDecompose_Helper<11,xx,xx,M1s>::call(A0,P+jk,detp);

                // Apply the permutation to the rest of the matrix
                if (jk > 0) {
                    //A2.cPermuteRows(P+jk,0,Nx);
                    InlinePermuteRows(A2,P+jk,0,Ny);
                }
                if (R < N) {
                    //A1.cPermuteRows(P+jk,0,Nx);
                    InlinePermuteRows(A1,P+jk,0,Ny);

                    // Solve for U01
                    //A01 /= A00.unitLowerTri();
                    M1l L00 = A00.unitLowerTri();
                    LDivEqMU_Helper<-2,xx,xx,M1s,M1l>::call(A01,L00);
                }
                for(int i=jk;i<R;++i) P[i]+=jk;
            }
        }
    };

    // algo 27: Recursive algorithm
    template <int cs, int rs, class M1>
    struct LUDecompose_Helper<27,cs,rs,M1>
    {
        static void call(M1& A, int* P, int& detp)
        {
            // The recursive LU algorithm is similar to the block algorithm,
            // except that the block is roughly half the size of the whole 
            // matrix.
            // We keep dividing the matrix in half (column-wise) until we 
            // get down to an Mx2 or Mx1 matrix.
            typedef typename M1::real_type RT;

            const int N = rs==UNKNOWN ? int(A.rowsize()) : rs;
            const int M = cs==UNKNOWN ? int(A.colsize()) : cs;
            const int R = TMV_MIN(N,M);
#ifdef PRINTALGO_LU
            std::cout<<"LUDecompose algo 27: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const int rx = IntTraits2<rs,cs>::min;
            const int algo2a = 
                (rs==UNKNOWN || rs==1) ? 1 : 0;
            const int algo2b = 
                TMV_Q2 == 1 ? 0 :
                (rs==UNKNOWN || rs==2) && (cs==UNKNOWN || cs>2) ? 2 : 0;
            const int algo2c = 
                TMV_Q2 == 1 ? 0 :
                (cs==UNKNOWN || cs==2) ? 3 : 0;
            const int algo2d = 
                TMV_Q2 <= 2 ? 0 :
                (rx != UNKNOWN && rx > TMV_Q2) ? 0 :
                11;

            const int algo3 =  // The algorithm for R > 32
                (rx == UNKNOWN || rx > 32) ? 27 : 0;
            const int algo4 =  // The algorithm for MultMM, LDivEqMU
                rx == UNKNOWN ? -2 : rx > 32 ? -3 : 0;
#if TMV_Q2 < 32
            const int algo3b =  // The algorithm for R > Q2
                (rx == UNKNOWN || rx > TMV_Q2) ? 27 : 0;
            const int algo4b =  // The algorithm for MultMM, LDivEqMU
                (rx == UNKNOWN || rx > TMV_Q2) ? -4 : 0;
#endif
            typedef typename M1::colrange_type M1c;
            typedef typename M1c::rowrange_type M1s;
            typedef typename M1c::const_rowrange_type M1sc;
            typedef typename M1s::const_unit_lowertri_type M1l;

            if (R > 32) { // For large R, make sure to use good MultMM algo

                const int Nx = (((R-1)>>5)+1)<<4;
                // (Round R/2 up to a multiple of 16.)
                const int rsx = IntTraits<rx>::half_roundup;
                const int rsy = IntTraits2<rs,rsx>::diff;
                const int csy = IntTraits2<cs,rsx>::diff;

                M1c A0 = A.cColRange(0,Nx);
                M1s A00 = A0.cRowRange(0,Nx);
                M1s A10 = A0.cRowRange(Nx,M);
                M1c A1 = A.cColRange(Nx,N);
                M1s A01 = A1.cRowRange(0,Nx);
                M1s A11 = A1.cRowRange(Nx,M);
                M1l L00 = A00.unitLowerTri();

                // Decompose left half into PLU
                LUDecompose_Helper<algo3,cs,rsx,M1c>::call(A0,P,detp);

                // Apply the permutation to the right half of the matrix
                //A1.cPermuteRows(P,0,Nx);
                InlinePermuteRows(A1,P,0,Nx);

                // Solve for U01
                //A01 /= A00.unitLowerTri();
                LDivEqMU_Helper<algo4,rsx,rsy,M1s,M1l>::call(A01,L00);

                // Solve for A~
                //A11 -= A10 * A01;
                MultMM_Helper<algo4,csy,rsy,rsx,true,-1,RT,M1sc,M1sc,M1s>::call(
                    Scaling<-1,RT>(),A10,A01,A11);

                // Decompose A~ into PLU
                LUDecompose_Helper<algo3,csy,rsy,M1s>::call(A11,P+Nx,detp);
                for(int i=Nx;i<R;++i) P[i]+=Nx;

                // Apply the new permutations to the left half
                //A0.cPermuteRows(P,Nx,R);
                InlinePermuteRows(A0,P,Nx,R);

#if TMV_Q2 < 32
            } else if (R > TMV_Q2) { // For small R, use simpler algo inline

                const int Nx = R > 16 ? ((((R-1)>>5)+1)<<4) : (R>>1);
                // (If R > 16, round R/2 up to a multiple of 16.)
                const int rsx = IntTraits<rx>::half_roundup;
                const int rsy = IntTraits2<rs,rsx>::diff;
                const int csy = IntTraits2<cs,rsx>::diff;

                M1c A0 = A.cColRange(0,Nx);
                M1s A00 = A0.cRowRange(0,Nx);
                M1s A10 = A0.cRowRange(Nx,M);
                M1c A1 = A.cColRange(Nx,N);
                M1s A01 = A1.cRowRange(0,Nx);
                M1s A11 = A1.cRowRange(Nx,M);
                M1l L00 = A00.unitLowerTri();

                // Decompose left half into PLU
                LUDecompose_Helper<algo3b,cs,rsx,M1c>::call(A0,P,detp);

                // Apply the permutation to the right half of the matrix
                //A1.cPermuteRows(P,0,Nx);
                InlinePermuteRows(A1,P,0,Nx);

                // Solve for U01
                //A01 /= A00.unitLowerTri();
                LDivEqMU_Helper<algo4b,rsx,rsy,M1s,M1l>::call(A01,L00);

                // Solve for A~
                //A11 -= A10 * A01;
                MultMM_Helper<algo4b,csy,rsy,rsx,true,-1,RT,M1sc,M1sc,M1s>::call(
                    Scaling<-1,RT>(),A10,A01,A11);

                // Decompose A~ into PLU
                LUDecompose_Helper<algo3b,csy,rsy,M1s>::call(A11,P+Nx,detp);
                for(int i=Nx;i<R;++i) P[i]+=Nx;

                // Apply the new permutations to the left half
                //A0.cPermuteRows(P,Nx,R);
                InlinePermuteRows(A0,P,Nx,R);
#endif
            }
            else if (N == 1) 
                LUDecompose_Helper<algo2a,cs,1,M1>::call(A,P,detp);
            else if (N == 2 && M > 2) 
                LUDecompose_Helper<algo2b,cs,2,M1>::call(A,P,detp);
            else if (M == 2) 
                LUDecompose_Helper<algo2c,2,rs,M1>::call(A,P,detp);
            else if (R > 2) 
                LUDecompose_Helper<algo2d,cs,rs,M1>::call(A,P,detp);
            else 
                TMVAssert(N==0 || M==0 || M==1); // and thus nothing to do.
        }
    };

    // algo 81: Copy to colmajor
    template <int cs, int rs, class M>
    struct LUDecompose_Helper<81,cs,rs,M>
    {
        static void call(M& m, int* P, int& detp)
        {
#ifdef PRINTALGO_LU
            std::cout<<"LUDecompose algo 81: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M::value_type T;
            typedef typename MCopyHelper<T,Rec,cs,rs,false,false>::type Mcm;
            Mcm mcm = m;
            LUDecompose_Helper<-1,cs,rs,Mcm>::call(mcm,P,detp);
            NoAliasCopy(mcm,m);
        }
    };

    // algo -4: No copies or branches
    template <int cs, int rs, class M1>
    struct LUDecompose_Helper<-4,cs,rs,M1>
    {
        static void call(M1& m, int* P, int& detp)
        {
#if 0
            const int algo = 27;
#else
            const int algo = 
                cs == 0 || rs == 0 || cs == 1 ? 0 :
                rs == 1 ? 1 : 
                cs == 2 ? 3 :
                rs == 2 ? 2 :
                27;
#endif
#ifdef PRINTALGO_LU
            std::cout<<"Inline LUDecompose: \n";
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"sizes = "<<m.colsize()<<"  "<<m.rowsize()<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            //std::cout<<"m = "<<m<<std::endl;
            //std::cout<<"detp = "<<detp<<std::endl;
#endif
            LUDecompose_Helper<algo,cs,rs,M1>::call(m,P,detp);
#ifdef PRINTALGO_LU
            //std::cout<<"m => "<<m<<std::endl;
            //std::cout<<"P => "<<ConstVectorView<int>(P,m.rowsize(),1)<<std::endl;
            //std::cout<<"detp -> "<<detp<<std::endl;
#endif
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1>
    struct LUDecompose_Helper<-3,cs,rs,M1>
    {
        static void call(M1& m, int* P, int& detp)
        {
            const int algo = (
                ( cs != UNKNOWN && rs != UNKNOWN &&
                  cs <= 16 && rs <= 16 ) ? -4 :
                !M1::_colmajor ? 81 :
                -4 );
#ifdef PRINTALGO_LU
            const int M = cs==UNKNOWN ? int(m.colsize()) : cs;
            const int N = rs==UNKNOWN ? int(m.rowsize()) : rs;
            std::cout<<"LUDecompose algo -3: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            LUDecompose_Helper<algo,cs,rs,M1>::call(m,P,detp);
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M>
    struct LUDecompose_Helper<97,cs,rs,M>
    {
        static void call(M& m, int* P, int& detp)
        {
            typedef typename M::conjugate_type Mc;
            Mc mc = m.conjugate();
            LUDecompose_Helper<-1,cs,rs,Mc>::call(mc,P,detp);
        }
    };

    // algo 98: call InstLU_Decompose
    template <int cs, int rs, class M>
    struct LUDecompose_Helper<98,cs,rs,M>
    {
        static void call(M& m, int* P, int& detp)
        { InstLU_Decompose(m.xView().cmView(),P,detp); }
    };

    // algo -1: Check for inst
    template <int cs, int rs, class M>
    struct LUDecompose_Helper<-1,cs,rs,M>
    {
        static void call(M& m, int* P, int& detp)
        {
            typedef typename M::value_type T;
            const bool inst = 
                (cs == UNKNOWN || cs > 16) &&
                (rs == UNKNOWN || rs > 16) &&
                Traits<T>::isinst;
            const int algo = 
                cs == 0 || rs == 0 || cs == 1 ? 0 :
                M::_conj ? 97 :
                inst ? (M::_colmajor ? 98 : 81) :
                -3;
            LUDecompose_Helper<algo,cs,rs,M>::call(m,P,detp);
        }
    };

    template <class M> 
    static void InlineLU_Decompose(
        BaseMatrix_Rec_Mutable<M>& m, int* P, int& detp)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::cview_type Mv;
        Mv mv = m.cView();
        LUDecompose_Helper<-3,cs,rs,Mv>::call(mv,P,detp);
    }

    template <class M> 
    static void LU_Decompose(
        BaseMatrix_Rec_Mutable<M>& m, int* P, int& detp)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::cview_type Mv;
        Mv mv = m.cView();
        LUDecompose_Helper<-1,cs,rs,Mv>::call(mv,P,detp);
    }

    // This function is a friend of Permutation class.
    template <class M>
    static void LU_Decompose(BaseMatrix_Rec_Mutable<M>& m, Permutation& P)
    {
        TMVAssert(P.size() == m.colsize());
        P.allocateMem();
        LU_Decompose(m,P.getMem(),P.itsdet=1);
        P.isinv = true;
    }

    // Allow views as an argument by value (for convenience)
    template <class T, int Si, int Sj, bool C, IndexStyle I>
    static void LU_Decompose(MatrixView<T,Si,Sj,C,I> m, Permutation& P)
    { 
        typedef MatrixView<T,Si,Sj,C,I> M;
        LU_Decompose(static_cast<BaseMatrix_Rec_Mutable<M>&>(m),P); 
    }
    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    static void LU_Decompose(
        SmallMatrixView<T,M,N,Si,Sj,C,I> m, Permutation& P)
    { 
        typedef SmallMatrixView<T,M,N,Si,Sj,C,I> MM;
        LU_Decompose(static_cast<BaseMatrix_Rec_Mutable<MM>&>(m),P); 
    }

#undef TMV_Q2
#undef TMV_LU_BLOCKSIZE

} // namespace tmv

#endif

