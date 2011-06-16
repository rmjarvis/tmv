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


#ifndef TMV_QRDecompose_H
#define TMV_QRDecompose_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Householder.h"
//#include "TMV_MultMV.h"
//#include "TMV_MultMM.h"
//#include "TMV_DivVU.h"
//#include "TMV_DivMU.h"

#ifdef PRINTALGO_QR
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#include "TMV_TriMatrixIO.h"
#include "TMV_ProdMM.h"
#include "TMV_MultUL.h"
#include "TMV_PermuteM.h"
#endif

// BLOCKSIZE is the block size to use in algo 21
#define TMV_QR_BLOCKSIZE 64

// RECURSE = the maximum size to stop recursing in algo 27
#if TMV_OPT >= 3
#define TMV_QR_RECURSE 2
#elif TMV_OPT == 2
#define TMV_QR_RECURSE 32
#else
#define TMV_QR_RECURSE 1
#endif

// INLINE_MV = Inline the MV (MultMV, LDivEqVU) calls.
#if TMV_OPT >= 1
#define TMV_QR_INLINE_MV
#endif

// INLINE_MM = Inline the small-sized MM (MultMM, LDivEqMU) calls.
#if TMV_OPT >= 3
#define TMV_QR_INLINE_MM
#endif

namespace tmv {

    // Defined in TMV_QRDecompose.cpp
    template <class T, class RT>
    void InstQR_Decompose(MatrixView<T> m, VectorView<RT> beta);

    template <int algo, int cs, int rs, class M, class V>
    struct QRDecompose_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or 1, or N == 0)
    template <int cs, int rs, class M, class V>
    struct QRDecompose_Helper<0,cs,rs,M,V>
    { static TMV_INLINE void call(M& A, V& beta) {} };

    // algo 11: Non-block algorithm, loop over n
    template <int cs, int rs, class M1, class V>
    struct QRDecompose_Helper<11,cs,rs,M1,V>
    {
        static void call(M1& A, V& beta)
        {
            // QR Decompostion
            //
            // We want to decompose the matrix (input as A) into Q * R
            // where Q is unitary and R is upper triangular.
            // Q and R are stored in the same matrix (output of A),
            // with the beta values for the Householder matrices returned
            // in beta.
            
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            const int M = cs==UNKNOWN ? int(A.colsize()) : cs;
            const int N = rs==UNKNOWN ? int(A.rowsize()) : rs;
#ifdef PRINTALGO_QR
            std::cout<<"QRDecompose algo 11: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename V::iterator IT;
            typedef typename M1::col_sub_type Mc;
            typedef typename M1::submatrix_type Ms;

            IT bj = beta.begin();
            for (int j=0;j<N;++j,++bj) {
                // Work on the lower part of column j
                Mc Acolj = A.col(j,j,M);
                // Compute the Householder reflection of this column
                // (and perform the reflection on this column).
                Householder<Mc> H(Acolj);
                *bj = H.getBeta();
                // Reflect the rest of the matrix to the right of this column.
                Ms Asub = A.subMatrix(j,M,j+1,N);
                H.multEq(Asub);
            }
        }
    };

#if 0
    // algo 21: Block algorithm
    template <int cs, int rs, class M1, class V>
    struct QRDecompose_Helper<21,cs,rs,M1,V>
    {
        typedef typename M1::value_type T;
        static void call(M1& A, V& beta)
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
            const int Nx = TMV_QR_BLOCKSIZE;
            const int R = TMV_MIN(N,M);
            const int xx = UNKNOWN;
#ifdef PRINTALGO_QR
            std::cout<<"QRDecompose algo 21: M,N,cs,rs = "<<M<<','<<N<<
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
                QRDecompose_Helper<11,xx,Nx,M1s>::call(A0,P+jk);

                // Apply the permutation to the rest of the matrix
                if (jk > 0) {
                    A2.cPermuteRows(P+jk,0,Nx);
                }
                if (jkpk < N) {
                    A1.cPermuteRows(P+jk,0,Nx);

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
                QRDecompose_Helper<11,xx,xx,M1s>::call(A0,P+jk);

                // Apply the permutation to the rest of the matrix
                if (jk > 0) {
                    A2.cPermuteRows(P+jk,0,Nx);
                }
                if (R < N) {
                    A1.cPermuteRows(P+jk,0,Nx);

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
    template <int cs, int rs, class M1, class V>
    struct QRDecompose_Helper<27,cs,rs,M1,V>
    {
        static void call(M1& A, V& beta)
        {
            // The recursive QR algorithm is similar to the block algorithm,
            // except that the block is roughly half the size of the whole 
            // matrix.
            // We keep dividing the matrix in half (column-wise) until we 
            // get down to an Mx2 or Mx1 matrix.
            typedef typename M1::real_type RT;

            const int N = rs==UNKNOWN ? int(A.rowsize()) : rs;
            const int M = cs==UNKNOWN ? int(A.colsize()) : cs;
            const int R = TMV_MIN(N,M);
#ifdef PRINTALGO_QR
            std::cout<<"QRDecompose algo 27: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const int rx = IntTraits2<rs,cs>::min;
            const int algo2a = 
                (rs==UNKNOWN || rs==1) ? 1 : 0;
#if TMV_QR_RECURSE > 1
            const int algo2b = 
                (rs==UNKNOWN || rs==2) && (cs==UNKNOWN || cs>2) ? 2 : 0;
            const int algo2c = 
                (cs==UNKNOWN || cs==2) ? 3 : 0;
#endif
#if TMV_QR_RECURSE > 2
            const int algo2d = 
                (rx != UNKNOWN && rx > TMV_QR_RECURSE) ? 0 :
                11;
#endif

            const int algo3 =  // The algorithm for R > 32
                (rx == UNKNOWN || rx > 32) ? 27 : 0;
            const int algo4 =  // The algorithm for MultMM, LDivEqMU
                rx == UNKNOWN ? -2 : rx > 32 ? -3 : 0;
#if TMV_QR_RECURSE < 32
            const int algo3b =  // The algorithm for R > QR_RECURSE
                (rx == UNKNOWN || rx > TMV_QR_RECURSE) ? 27 : 0;
            const int algo4b =  // The algorithm for R<32 MultMM, LDivEqMU
#ifdef TMV_QR_INLINE_MM
                rx == UNKNOWN ? -4 :
#else
                rx == UNKNOWN ? -402 :
#endif
                rx > TMV_QR_RECURSE ? -4 : 0;
#endif
            typedef typename M1::colrange_type M1c;
            typedef typename M1c::rowrange_type M1s;
            typedef typename M1c::const_rowrange_type M1sc;
            typedef typename M1s::const_unit_lowertri_type M1l;

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

            if (R > 32) { 
                // For large R, make sure to use good MultMM algo

                // Decompose left half into PQR
                QRDecompose_Helper<algo3,cs,rsx,M1c>::call(A0,P);

                // Apply the permutation to the right half of the matrix
                A1.cPermuteRows(P,0,Nx);

                // Solve for U01
                //A01 /= A00.unitLowerTri();
                LDivEqMU_Helper<algo4,rsx,rsy,M1s,M1l>::call(A01,L00);

                // Solve for A~
                //A11 -= A10 * A01;
                MultMM_Helper<algo4,csy,rsy,rsx,true,-1,RT,M1sc,M1sc,M1s>::call(
                    Scaling<-1,RT>(),A10,A01,A11);

                // Decompose A~ into PQR
                QRDecompose_Helper<algo3,csy,rsy,M1s>::call(A11,P+Nx);
                for(int i=Nx;i<R;++i) P[i]+=Nx;

                // Apply the new permutations to the left half
                A0.cPermuteRows(P,Nx,R);

#if TMV_QR_RECURSE < 32
            } else if (R > TMV_QR_RECURSE) {
                // For small R, use simpler inline algorithms

                // Decompose left half into PQR
                QRDecompose_Helper<algo3b,cs,rsx,M1c>::call(A0,P);

                // Apply the permutation to the right half of the matrix
                A1.cPermuteRows(P,0,Nx);

                // Solve for U01
                //A01 /= A00.unitLowerTri();
                LDivEqMU_Helper<algo4b,rsx,rsy,M1s,M1l>::call(A01,L00);

                // Solve for A~
                //A11 -= A10 * A01;
                MultMM_Helper<algo4b,csy,rsy,rsx,true,-1,RT,M1sc,M1sc,M1s>::call(
                    Scaling<-1,RT>(),A10,A01,A11);

                // Decompose A~ into PQR
                QRDecompose_Helper<algo3b,csy,rsy,M1s>::call(A11,P+Nx);
                for(int i=Nx;i<R;++i) P[i]+=Nx;

                // Apply the new permutations to the left half
                A0.cPermuteRows(P,Nx,R);
#endif
            } else if (N == 1) 
                QRDecompose_Helper<algo2a,cs,rs,M1>::call(A,P);
#if TMV_QR_RECURSE > 1
            else if (N == 2 && M > 2) 
                QRDecompose_Helper<algo2b,cs,rs,M1>::call(A,P);
            else if (M == 2) 
                QRDecompose_Helper<algo2c,cs,rs,M1>::call(A,P);
#endif
#if TMV_QR_RECURSE > 2
            else if (R > 2) 
                QRDecompose_Helper<algo2d,cs,rs,M1>::call(A,P);
#endif
            else 
                TMVAssert(N==0 || M==0 || M==1); // and thus nothing to do.
        }
    };
#endif

    // algo 81: Copy to colmajor
    template <int cs, int rs, class M, class V>
    struct QRDecompose_Helper<81,cs,rs,M,V>
    {
        static inline void call(M& m, V& beta)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRDecompose algo 81: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M::value_type T;
            typedef typename MCopyHelper<T,Rec,cs,rs,false,false>::type Mcm;
            Mcm mcm = m;
            QRDecompose_Helper<-2,cs,rs,Mcm,V>::call(mcm,beta);
            NoAliasCopy(mcm,m);
        }
    };

    // algo 90: call InstQR_Decompose
    template <int cs, int rs, class M, class V>
    struct QRDecompose_Helper<90,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta)
        { InstQR_Decompose(m.xView(),beta.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M, class V>
    struct QRDecompose_Helper<97,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta)
        {
            typedef typename M::conjugate_type Mc;
            Mc mc = m.conjugate();
            QRDecompose_Helper<-2,cs,rs,Mc,V>::call(mc,beta);
        }
    };

    // algo -4: No copies or branches
    template <int cs, int rs, class M, class V>
    struct QRDecompose_Helper<-4,cs,rs,M,V>
    {
        typedef typename M::value_type T;
        static TMV_INLINE void call(M& m, V& beta)
        {
            const int algo = 
                cs == 0 || rs == 0 || cs == 1 ? 0 :
                TMV_OPT == 0 ? 11 :
                11;
                //27;
#ifdef PRINTALGO_QR
            std::cout<<"Inline QRDecompose: \n";
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"sizes = "<<m.colsize()<<"  "<<m.rowsize()<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            QRDecompose_Helper<algo,cs,rs,M,V>::call(m,beta);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1, class V>
    struct QRDecompose_Helper<-3,cs,rs,M1,V>
    {
        static TMV_INLINE void call(M1& m, V& beta)
        {
            const int algo = (
                ( cs != UNKNOWN && rs != UNKNOWN &&
                  cs <= 16 && rs <= 16 ) ? -4 :
                !M1::_colmajor ? 81 :
                -4 );
#ifdef PRINTALGO_QR
            const int M = cs==UNKNOWN ? int(m.colsize()) : cs;
            const int N = rs==UNKNOWN ? int(m.rowsize()) : rs;
            std::cout<<"QRDecompose algo -3: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            QRDecompose_Helper<algo,cs,rs,M1,V>::call(m,beta);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M, class V>
    struct QRDecompose_Helper<-2,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta)
        {
            typedef typename M::value_type T;
            const bool inst = 
                (cs == UNKNOWN || cs > 16) &&
                (rs == UNKNOWN || rs > 16) &&
                Traits<T>::isinst;
            const int algo = 
                cs == 0 || rs == 0 || cs == 1 ? 0 :
                M::_conj ? 97 :
                inst ? 90 :
                -3;
            QRDecompose_Helper<algo,cs,rs,M,V>::call(m,beta);
        }
    };

    template <int cs, int rs, class M, class V>
    struct QRDecompose_Helper<-1,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& m, V& beta)
        { QRDecompose_Helper<-2,cs,rs,M,V>::call(m,beta); }
    };

    template <class M, class V>
    static inline void InlineQR_Decompose(
        BaseMatrix_Rec_Mutable<M>& m, BaseVector_Mutable<V>& beta)
    {
        const int cs = M::_colsize;
        const int rs = Sizes<M::_rowsize,V::_size>::size;
        typedef typename M::cview_type Mv;
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        TMV_MAYBE_REF(V,Vv) betav = beta.cView();
        QRDecompose_Helper<-3,cs,rs,Mv,Vv>::call(mv,betav);
    }

    // This is the basic functionality
    template <class M, class V>
    static inline void QR_Decompose(
        BaseMatrix_Rec_Mutable<M>& m, BaseVector_Mutable<V>& beta)
    {
        TMVStaticAssert(V::isreal);
        TMVStaticAssert((Traits2<
                         typename M::value_type,
                         typename V::value_type>::samebase));
        TMVAssert(m.colsize() >= m.rowsize());
        TMVAssert(beta.size() == m.rowsize());
        TMVStaticAssert((Sizes<M::_rowsize,V::_size>::same));
        TMVAssert(m.rowsize() == beta.size());

        const int cs = M::_colsize;
        const int rs = Sizes<M::_rowsize,V::_size>::size;
        typedef typename M::cview_type Mv;
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        TMV_MAYBE_REF(V,Vv) betav = beta.cView();
        QRDecompose_Helper<-2,cs,rs,Mv,Vv>::call(mv,betav);
    }

    // The rest of these below are basically convenience functions
    // to allow the user to provide fewer or different arguments.
    template <class M1, class M2>
    static inline void QR_Decompose(
        BaseMatrix_Rec_Mutable<M1>& Q, BaseMatrix_Tri_Mutable<M2>& R)
    {
        TMVStaticAssert((Traits2<
                         typename M1::value_type,
                         typename M2::value_type>::sametype));
        TMVStaticAssert(M2::_upper);
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(R.colsize() == Q.rowsize());

        typedef typename M1::real_type RT;
        Vector<RT> beta(Q.rowsize());
        QR_Decompose(Q,beta);
        NoAliasCopy(Q.upperTri(),R);
        UnpackQ(Q,beta);
    }

    template <class M>
    static inline void QR_Decompose(BaseMatrix_Rec_Mutable<M>& A)
    {
        TMVAssert(A.colsize() >= A.rowsize());
        typedef typename M::real_type RT;
        Vector<RT> beta(A.rowsize());
        QR_Decompose(A,beta);
    }


    // Allow views as an argument by value (for convenience)
    template <class T, int A, int A2>
    static inline void QR_Decompose(
        MatrixView<T,A> Q, UpperTriMatrixView<T,A2> R)
    {
        typedef MatrixView<T,A> M1;
        typedef UpperTriMatrixView<T,A2> M2;
        QR_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),
            static_cast<BaseMatrix_Tri_Mutable<M2>&>(R)); 
    }

    template <class T, int M, int N, int Si, int Sj, int A, int Si2, int Sj2, int A2>
    static inline void QR_Decompose(
        SmallMatrixView<T,M,N,Si,Sj,A> Q,
        SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> R)
    {
        typedef SmallMatrixView<T,M,N,Si,Sj,A> M1;
        typedef SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> M2;
        QR_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<M1>&>(Q),
            static_cast<BaseMatrix_Tri_Mutable<M2>&>(R));
    }

    template <class T, int A>
    static inline void QR_Decompose(MatrixView<T,A> m)
    {
        typedef MatrixView<T,A> M1;
        QR_Decompose(static_cast<BaseMatrix_Rec_Mutable<M1>&>(m));
    }

    template <class T, int M, int N, int Si, int Sj, int A>
    static inline void QR_Decompose(SmallMatrixView<T,M,N,Si,Sj,A> m)
    {
        typedef SmallMatrixView<T,M,N,Si,Sj,A> M1;
        QR_Decompose(static_cast<BaseMatrix_Rec_Mutable<M1>&>(m));
    }

#undef TMV_QR_RECURSE
#undef TMV_QR_BLOCKSIZE
#undef TMV_QR_INLINE_MV
#undef TMV_QR_INLINE_MM

} // namespace tmv

#endif

