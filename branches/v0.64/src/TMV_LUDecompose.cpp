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


#define XDEBUG


#include "TMV_Blas.h"
#include "tmv/TMV_LUD.h"
#include "TMV_LUDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_TriMatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_TriMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define LU_BLOCKSIZE TMV_BLOCKSIZE
#define LU_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define LU_BLOCKSIZE 64
#define LU_BLOCKSIZE2 2
#endif

    //
    // Decompose
    //

    template <class T> 
    static void NonBlockLUDecompose(
        const MatrixView<T>& A, int* P, int& detp)
    {
        // LU Decompostion with partial pivoting.
        //
        // We want to decompose the matrix (input as A) into P * L * U
        // where P is a permutation, L is a lower triangle matrix with 1's for 
        // the diagonal, and U is an upper triangle matrix.  
        //
        // We do this calculation one column at a time.  There are other versions
        // of this algorithm which use Rank 1 updates.  This version uses mostly
        // martix-vector products and forward substitutions.  The different
        // algorithms all require the same number of operation, but this
        // one requires fewer vector touches which often means it will be 
        // slightly faster.
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
        // This means we have some leeway on the values for L(j,j) and U(j,j),
        // as only their product is specified.
        //
        // If we take U to have unit diagonal, then L(j,j) is set here along
        // with the rest of the L(j:N,j) column.  However, this will mean that the
        // forward substutions in the (1*) steps will require divisions by the 
        // non-unit-diagonal elements of L.  It is faster to take L to have
        // unit-diagonal elements, and do the division by U(j,j) here, since then
        // we can calculate 1/U(j,j) and multiply.  So 1 division and N-j
        // multiplies which is generally faster than N-j divisions.
        //
        // However, another potential problem is that U(j,j) could be 0, or close
        // to 0.  This would lead to either an error or inaccurate results.
        // Thus we add a step in the middle of the (2*) calculation:
        //
        // Define v(j:N) = A(j:N,j) - L(j:N,0:j) U(0:j,j)
        // 
        // We search v for the element with the largest absolute value and apply
        // a permutation to swap it into the j spot.  This element then becomes 
        // U(j,j), which is then the divisor for the rest of the vector.  This 
        // will minimize the possibility of roundoff errors due to small U(j,j)'s.
         
        
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.iscm());
        typedef TMV_RealType(T) RT;

        // I use x*halfeps to check for underflow.  Normally x*eps should
        // be sufficient.  But if x is complex, then 1/x can still overflow
        // if x*eps != 0, but x*halfeps == 0.
        const RT halfeps = TMV_Epsilon<RT>()/RT(2);

        const int N = A.rowsize();
        const int M = A.colsize();
        const int R = TMV_MIN(N,M);
#ifdef XDEBUG
        Matrix<T> A0(A);
#endif

        const T* Ujj = A.cptr();
        const int Ads = A.stepj()+1;
        int* Pj = P;

        for (int j=0; j<R; ++j,Ujj+=Ads,++Pj) {
            if (j > 0) {
                // Solve for U(0:j,j))
                A.col(j,0,j) /= A.subMatrix(0,j,0,j).lowerTri(UnitDiag);

                // Solve for v = L(j:M,j) U(j,j)
                A.col(j,j,M) -= A.subMatrix(j,M,0,j) * A.col(j,0,j);
            }

            // Find the pivot element
            int ip;
            RT piv = A.col(j,j,M).maxAbsElement(&ip);
            // ip is relative to j index, not absolute.
            
            // Check for underflow:
            if (piv * halfeps == RT(0)) {
                *Pj = j;
                A.col(j,j,M).setZero();
                continue;
            }

            // Swap the pivot row with j if necessary
            if (ip != 0) {
                ip += j;
                TMVAssert(ip < int(A.colsize()));
                TMVAssert(j < int(A.colsize()));
                A.swapRows(ip,j);  // This does both Lkb and A'
                *Pj = ip;
                detp = -detp;
            } else *Pj = j;

            // Solve for L(j+1:M,j)
            // If Ujj is 0, then all of the L's are 0.
            // ie. Ujj Lij = 0 for all i>j
            // Any value for Lij is valid, so leave them 0.
            if (*Ujj != T(0)) A.col(j,j+1,M) /= *Ujj;
        }
        if (N > M) {
            // Solve for U(0:M,M:N))
            A.colRange(M,N) /= A.colRange(0,M).lowerTri(UnitDiag);
        }
#ifdef XDEBUG
        Matrix<T> L = A;
        if (R > 1)
            L.upperTri().offDiag().setZero();
        L.diag().setAllTo(T(1));
        Matrix<T> U = A.upperTri();
        Matrix<T> AA = L * U;
        AA.reversePermuteRows(P,0,R);
        if (Norm(AA-A0) > 0.001*Norm(A0) && Norm(A0)*halfeps > RT(0)) {
            cerr<<"Done NonBlock LU: \n";
            cerr<<"A0 = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"LU = "<<A<<endl;
            cerr<<"P = (";
            for(int i=0;i<R;i++) cerr<<P[i]<<" ";
            cerr<<")\n";
            cerr<<"AA = "<<AA<<endl;
            cerr<<"Norm(A0-AA) = "<<Norm(AA-A0)<<endl;
            cerr<<"Norm(A0) = "<<Norm(A0)<<endl;
            abort();
        }
#endif
    }

    template <class T> 
    static void RecursiveLUDecompose(
        const MatrixView<T>& A, int* P, int& detp)
    {
        //cout<<"Start recursive LU: A = "<<A<<endl;
        // The recursive LU algorithm is similar to the block algorithm, except 
        // that the block is roughly half the size of the whole matrix.
        // We keep dividing the matrix in half (column-wise) until we get down
        // to an Mx2 or Mx1 matrix.

#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<T,ColMajor> A2 = A;
        std::vector<int> P2(A.colsize());
        int detp2=1;
        NonBlockLUDecompose(A2.view(),&P2[0],detp2);
#endif

        typedef TMV_RealType(T) RT;
        const RT halfeps = TMV_Epsilon<RT>()/RT(2);

        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.iscm());
        const int N = A.rowsize();
        const int M = A.colsize();
        const int R = TMV_MIN(N,M);

        if (R > LU_BLOCKSIZE2) {
            // Split N in half, with N1 being rounded to multiple of BLOCKSIZE
            // if appropriate.
            int N1 = R/2;
            if (N1 > LU_BLOCKSIZE) N1 = (N1/LU_BLOCKSIZE)*LU_BLOCKSIZE;

            MatrixView<T> A0 = A.colRange(0,N1);
            MatrixView<T> A00 = A0.rowRange(0,N1);
            MatrixView<T> A10 = A0.rowRange(N1,M);
            MatrixView<T> A1 = A.colRange(N1,N);
            MatrixView<T> A01 = A1.rowRange(0,N1);
            MatrixView<T> A11 = A1.rowRange(N1,M);

            // Decompose left half into PLU
            RecursiveLUDecompose(A0,P,detp);

            // Apply the permutation to the right half of the matrix
            A1.permuteRows(P,0,N1);

            // Solve for U01
            A01 /= A00.lowerTri(UnitDiag);

            // Solve for A~
            A11 -= A10 * A01;

            // Decompose A~ into PLU
            RecursiveLUDecompose(A11,P+N1,detp);
            for(int i=N1;i<R;++i) P[i]+=N1;

            // Apply the new permutations to the left half
            A0.permuteRows(P,N1,R);
        } else if (LU_BLOCKSIZE2 > 2 && R > 2) {
            NonBlockLUDecompose(A,P,detp);
        } else if (R == 2) {
            // Same as NonBlock version, but with R==2 hard coded
            VectorView<T> A0 = A.col(0);
            VectorView<T> A1 = A.col(1);

            int ip0,ip1;
            RT piv = A0.maxAbsElement(&ip0);

            // Check for underflow:
            if (piv * halfeps == RT(0)) {
                ip0 = 0;
                piv = RT(0);
                A0.setZero();
            }

            if (piv != RT(0)) {
                if (ip0 != 0) {
                    A0.swap(ip0,0);
                    A1.swap(ip0,0);
                    detp = -detp;
                } 

                // A0.subVector(1,M) /= A00;
                // A1.subVector(1,M) -= A0.subVector(1,M) * A01;
                const T invA00 = RT(1)/(*A0.cptr());
                const T A01 = (*A1.cptr());
                piv = RT(0); // next pivot element
                ip1 = 1;
                T* Ai0 = A0.ptr()+1;
                T* Ai1 = A1.ptr()+1;
                for(int i=1;i<M;++i,++Ai0,++Ai1) {
#ifdef TMVFLDEBUG
                    TMVAssert(Ai0 >= A.first);
                    TMVAssert(Ai0 < A.last);
                    TMVAssert(Ai1 >= A.first);
                    TMVAssert(Ai1 < A.last);
#endif
                    *Ai0 *= invA00;
                    *Ai1 -= *Ai0 * A01;
                    RT absAi1 = TMV_ABS(*Ai1);
                    if (absAi1 > piv) { piv = absAi1; ip1=i; }
                }
            } else {
                piv = A1.subVector(1,M).maxAbsElement(&ip1); 
                ++ip1;
            }

            // Check for underflow:
            if (piv * halfeps == RT(0)) {
                piv = RT(0);
                ip1 = 1;
                A1.subVector(1,M).setZero();
            }

            if (piv != RT(0) && M>2) {
                if (ip1 != 1) {
                    A1.swap(ip1,1);
                    A0.swap(ip1,1);
                    detp = -detp;
                } 

                //A1.subVector(2,M) /= A1(1);
                const T A11 = (*(A1.cptr()+1));
                A1.subVector(2,M) /= A11;
            }

            if (N > 2) {
                // M=2, N>2, so solve for U(0:2,2:N))
                // A.colRange(2,N).permuteRows(P);
                if (*P == 1) A.colRange(2,N).swapRows(0,1);
                // A.colRange(2,N) /= A.colRange(0,2).lowerTri(UnitDiag);
                const T A01 = (*A1.cptr());
                A.row(1,2,N) -= A01 * A.row(0,2,N);
            }
            P[0] = ip0;
            P[1] = ip1;
        } else if (R == 1) {
            // Same as NonBlock version, but with R==1 hard coded
            VectorView<T> A0 = A.col(0);

            RT piv = A0.maxAbsElement(P);

            // Check for underflow:
            if (piv * halfeps == RT(0)) {
                piv = RT(0);
                *P = 0;
                A0.setZero();
            }

            if (piv != RT(0)) {
                if (*P != 0) {
                    A0.swap(*P,0);
                    detp = -detp;
                }
                A0.subVector(1,M) /= (*A0.cptr());
            }
        }
#ifdef XDEBUG
        Matrix<T> L = A;
        if (R > 1) L.upperTri().offDiag().setZero();
        L.diag().setAllTo(T(1));
        Matrix<T> U = A.upperTri();
        Matrix<T> AA = L * U;
        AA.reversePermuteRows(P,0,R);
        if (Norm(AA-A0) > 0.001*Norm(A0) && Norm(A0)*halfeps > RT(0)) {
            cerr<<"Done Recursive LU: \n";
            cerr<<"A0 = "<<A0<<endl;
            cerr<<"LU = "<<A<<endl;
            cerr<<"P = (";
            for(int i=0;i<R;i++) cerr<<P[i]<<" ";
            cerr<<")\n";
            cerr<<"A2 = "<<A2<<endl;
            cerr<<"Norm(A-A2) = "<<Norm(A-A2)<<endl;
            cerr<<"Norm(AA-A0) = "<<Norm(AA-A0)<<endl;
            cerr<<"correct P = (";
            for(int i=0;i<R;i++) cerr<<P2[i]<<" ";
            cerr<<")\n";
            abort();
        }
#endif
    }

    template <class T> 
    static inline void NonLapLUDecompose(
        const MatrixView<T>& A, int* P, int& detp)
    {
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.iscm());

        RecursiveLUDecompose(A,P,detp);
    }

#ifdef ALAP
    template <class T> 
    static void LapLUDecompose(const MatrixView<T>& A, int* P, int& detp)
    { NonLapLUDecompose(A,P,detp); }
#ifdef INST_DOUBLE
    template <> 
    void LapLUDecompose(const MatrixView<double>& A, int* P, int& detp)
    {
        TMVAssert(A.iscm());
        TMVAssert(A.ct()==NonConj);

        int m = A.colsize();
        int n = A.rowsize();
        int lda = A.stepj();
        auto_array<int> lap_p(new int[n]);
        LAPNAME(dgetrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
                         LAPP(lap_p.get()) LAPINFO);
        LAP_Results("dgetrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = (lap_p.get())[i] LAPMINUS1;
            if (P[i]!=i) detp = -detp;
        }
    }
    template <> 
    void LapLUDecompose(
        const MatrixView<std::complex<double> >& A, int* P, int& detp)
    {
        TMVAssert(A.iscm());
        TMVAssert(A.ct()==NonConj);

        int m = A.colsize();
        int n = A.rowsize();
        int lda = A.stepj();
        auto_array<int> lap_p(new int[n]);
        LAPNAME(zgetrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
                         LAPP(lap_p.get()) LAPINFO);
        LAP_Results("zgetrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = (lap_p.get())[i] LAPMINUS1;
            if (P[i]!=i) detp = -detp;
        }
    }
#endif
#ifdef INST_FLOAT
#ifndef MKL
    // This is giving me a weird runtime error sometimes with MKL:
    // OMP abort: Unable to set worker thread stack size to 2098176 bytes
    // Try reducing KMP_STACKSIZE or increasing the shell stack limit.
    // So I'm cutting it out for MKL compilations
    template <> 
    void LapLUDecompose(const MatrixView<float>& A, int* P, int& detp)
    {
        TMVAssert(A.iscm());
        TMVAssert(A.ct()==NonConj);

        int m = A.colsize();
        int n = A.rowsize();
        int lda = A.stepj();
        auto_array<int> lap_p(new int[n]);
        LAPNAME(sgetrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
                         LAPP(lap_p.get()) LAPINFO);
        LAP_Results("sgetrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = (lap_p.get())[i] LAPMINUS1;
            if (P[i]!=i) detp = -detp;
        }
    }
    template <> 
    void LapLUDecompose(
        const MatrixView<std::complex<float> >& A, int* P, int& detp)
    {
        TMVAssert(A.iscm());
        TMVAssert(A.ct()==NonConj);

        int m = A.colsize();
        int n = A.rowsize();
        int lda = A.stepj();
        auto_array<int> lap_p(new int[n]);
        LAPNAME(cgetrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
                         LAPP(lap_p.get()) LAPINFO);
        LAP_Results("cgetrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = (lap_p.get())[i] LAPMINUS1;
            if (P[i]!=i) detp = -detp;
        }
    }
#endif // MKL
#endif // FLOAT
#endif // ALAP

    template <class T> 
    void LU_Decompose(const MatrixView<T>& A, int* P, int& detp)
    {
#ifdef XDEBUG
        std::cout<<"Start LUDecompose:"<<std::endl;
        std::cout<<"A = "<<A<<std::endl;
        Matrix<T> A0 = A;
#endif
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.iscm());

        if (A.colsize() > 0 && A.rowsize() > 0) {
#ifdef ALAP
            LapLUDecompose(A,P,detp);
#else
            NonLapLUDecompose(A,P,detp);
#endif
        }
#ifdef XDEBUG
        std::cout<<"A => "<<A<<std::endl;
        Matrix<T> A2 = A.unitLowerTri() * A.upperTri();
        A2.reversePermuteRows(P);
        std::cout<<"Norm(A0-A2) = "<<Norm(A0-A2)<<std::endl;
        typedef TMV_RealType(T) RT;
        const RT halfeps = TMV_Epsilon<RT>()/RT(2);
        if (Norm(A2-A0) > 0.001*Norm(A0) && Norm(A0)*halfeps > RT(0)) {
            cerr<<"Done LU_Decompose: \n";
            cerr<<"A0 = "<<A0<<endl;
            cerr<<"LU = "<<A<<endl;
            cerr<<"P = (";
            for(int i=0;i<int(A.rowsize());i++) cerr<<P[i]<<" ";
            cerr<<")\n";
            cerr<<"A2 = "<<A2<<endl;
            cerr<<"Norm(A0-A2) = "<<Norm(A0-A2)<<endl;
            abort();
        }
#endif
    }

    template <class T> 
    void LU_Decompose(const MatrixView<T>& A, int* P)
    { 
        if (A.isconj()) LU_Decompose(A.conjugate(),P);
        else {
            int d=0;
            if (A.isrm()) {
                A.transposeSelf();
                LU_Decompose(A.transpose(),P,d);
                A.transposeSelf();
            } else {
                LU_Decompose(A,P,d); 
            }
        }
    }

#define InstFile "TMV_LUDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


