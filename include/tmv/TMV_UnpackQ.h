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


#ifndef TMV_UnpackQ_H
#define TMV_UnpackQ_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_MultMV.h"
#include "TMV_MultMM.h"
#include "TMV_DivVU.h"
#include "TMV_DivMU.h"
#include "TMV_Householder.h"

#ifdef PRINTALGO_QR
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#include "TMV_TriMatrixIO.h"
#include "TMV_ProdMM.h"
#include "TMV_MultUL.h"
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

    // Defined in TMV_UnpackQ.cpp
    template <class T, class RT>
    void InstUnpackQ(MatrixView<T> Q, const ConstVectorView<RT>& beta);

    template <int algo, int cs, int rs, class M, class V>
    struct UnpackQ_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or 1, or N == 0)
    template <int cs, int rs, class M, class V>
    struct UnpackQ_Helper<0,cs,rs,M,V>
    { static TMV_INLINE void call(M& A, const V& beta) {} };

    // algo 11: Non-block algorithm, loop over n
    template <int cs, int rs, class M1, class V>
    struct UnpackQ_Helper<11,cs,rs,M1,V>
    {
        static void call(M1& Q, const V& beta)
        {
            // This is essentially the reverse of the QR Decomposition
            // algorithm where R is taken to be an identity matrix.
            
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            const int M = cs==UNKNOWN ? int(Q.colsize()) : cs;
            const int N = rs==UNKNOWN ? int(Q.rowsize()) : rs;
#ifdef PRINTALGO_QR
            std::cout<<"UnpackQ algo 11: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename V::const_reverse_iterator IT;
            typedef typename M1::reference Mr;
            typedef typename M1::col_sub_type Mc;
            typedef typename M1::submatrix_type Ms;

            Q.upperTri().offDiag().setZero();

            IT bj = beta.rbegin(); // rbegin, so iterate from end.
            for (int j=N-1;j>=0;--j,++bj) {
                // Work on the lower part of column j
                Mr Qjj = Q.ref(j,j);
                Mc Qcolj = Q.col(j,j+1,M);
                Householder<Mc> H(Qcolj,*bj);
                // Reflect the rest of the matrix to the right of this column.
                Ms Qsub = Q.subMatrix(j,M,j+1,N);
                H.multEq(Qsub);
                H.unpack(Qjj);
            }
        }
    };

#if 0
    // algo 21: Block algorithm
    template <int cs, int rs, class M1, class V>
    struct UnpackQ_Helper<21,cs,rs,M1,V>
    {
        typedef typename M1::value_type T;
        static void call(M1& Q, const V& beta)
    };

    // algo 27: Recursive algorithm
    template <int cs, int rs, class M1, class V>
    struct UnpackQ_Helper<27,cs,rs,M1,V>
    {
        static void call(M1& Q, const V& beta)
    };
#endif

    // algo 81: Copy to colmajor
    template <int cs, int rs, class M, class V>
    struct UnpackQ_Helper<81,cs,rs,M,V>
    {
        static inline void call(M& Q, const V& beta)
        {
#ifdef PRINTALGO_QR
            std::cout<<"UnpackQ algo 81: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M::value_type T;
            typedef typename MCopyHelper<T,Rec,cs,rs,false,false>::type Mcm;
            Mcm Qcm = Q;
            UnpackQ_Helper<-2,cs,rs,Mcm,V>::call(Qcm,beta);
            NoAliasCopy(Qcm,Q);
        }
    };

    // algo 90: call InstUnpackQ
    template <int cs, int rs, class M, class V>
    struct UnpackQ_Helper<90,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& Q, const V& beta)
        { InstUnpackQ(Q.xView(),beta.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M, class V>
    struct UnpackQ_Helper<97,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& Q, const V& beta)
        {
            typedef typename M::conjugate_type Mc;
            Mc Qc = Q.conjugate();
            UnpackQ_Helper<-2,cs,rs,Mc,V>::call(Qc,beta);
        }
    };

    // algo -4: No copies or branches
    template <int cs, int rs, class M1, class V>
    struct UnpackQ_Helper<-4,cs,rs,M1,V>
    {
        typedef typename M1::value_type T;
        static TMV_INLINE void call(M1& Q, const V& beta)
        {
            const int algo = 
                cs == 0 || rs == 0 || cs == 1 ? 0 :
                TMV_OPT == 0 ? 11 :
                11;
                //27;
#ifdef PRINTALGO_QR
            std::cout<<"Inline UnpackQ: \n";
            std::cout<<"Q = "<<TMV_Text(Q)<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"sizes = "<<Q.colsize()<<"  "<<Q.rowsize()<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            UnpackQ_Helper<algo,cs,rs,M1,V>::call(Q,beta);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1, class V>
    struct UnpackQ_Helper<-3,cs,rs,M1,V>
    {
        static TMV_INLINE void call(M1& Q, const V& beta)
        {
            const int algo = (
                ( cs != UNKNOWN && rs != UNKNOWN &&
                  cs <= 16 && rs <= 16 ) ? -4 :
                !M1::_colmajor ? 81 :
                -4 );
#ifdef PRINTALGO_QR
            const int M = cs==UNKNOWN ? int(Q.colsize()) : cs;
            const int N = rs==UNKNOWN ? int(Q.rowsize()) : rs;
            std::cout<<"UnpackQ algo -3: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            UnpackQ_Helper<algo,cs,rs,M1,V>::call(Q,beta);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M, class V>
    struct UnpackQ_Helper<-2,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& Q, const V& beta)
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
            UnpackQ_Helper<algo,cs,rs,M,V>::call(Q,beta);
        }
    };

    template <int cs, int rs, class M, class V>
    struct UnpackQ_Helper<-1,cs,rs,M,V>
    {
        static TMV_INLINE void call(M& Q, const V& beta)
        { UnpackQ_Helper<-2,cs,rs,M,V>::call(Q,beta); }
    };

    template <class M, class V>
    static inline void InlineUnpackQ(
        BaseMatrix_Rec_Mutable<M>& Q, const BaseVector_Calc<V>& beta)
    {
        TMVStaticAssert(V::isreal);
        TMVStaticAssert((Traits2<
                         typename M::value_type,
                         typename V::value_type>::samebase));
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVStaticAssert((Sizes<M::_rowsize,V::_size>::same));
        TMVAssert(Q.rowsize() == beta.size());

        const int cs = M::_colsize;
        const int rs = Sizes<M::_rowsize,V::_size>::size;
        typedef typename M::cview_type Mv;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_REF(M,Mv) Qv = Q.cView();
        TMV_MAYBE_CREF(V,Vv) betav = beta.cView();
        UnpackQ_Helper<-3,cs,rs,Mv,Vv>::call(Qv,betav);
    }

    // This is the basic functionality
    template <class M, class V>
    static inline void UnpackQ(
        BaseMatrix_Rec_Mutable<M>& Q, const BaseVector_Calc<V>& beta)
    {
        TMVStaticAssert(V::isreal);
        TMVStaticAssert((Traits2<
                         typename M::value_type,
                         typename V::value_type>::samebase));
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVStaticAssert((Sizes<M::_rowsize,V::_size>::same));
        TMVAssert(Q.rowsize() == beta.size());

        const int cs = M::_colsize;
        const int rs = Sizes<M::_rowsize,V::_size>::size;
        typedef typename M::cview_type Mv;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_REF(M,Mv) Qv = Q.cView();
        TMV_MAYBE_CREF(V,Vv) betav = beta.cView();
        UnpackQ_Helper<-2,cs,rs,Mv,Vv>::call(Qv,betav);
    }

#undef TMV_QR_RECURSE
#undef TMV_QR_BLOCKSIZE
#undef TMV_QR_INLINE_MV
#undef TMV_QR_INLINE_MM

} // namespace tmv

#endif

