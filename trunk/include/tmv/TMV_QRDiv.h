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


#ifndef TMV_QRDiv_H
#define TMV_QRDiv_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_DivVM_Funcs.h"
#include "TMV_DivMM_Funcs.h"
#include "TMV_MultMM_Funcs.h"
#include "TMV_Permutation.h"
#include "TMV_PackedQ.h"

#ifdef PRINTALGO_QR
#include <iostream>
#include "TMV_ProdMV.h"
#include "TMV_ProdMM.h"
#endif

namespace tmv {

    // Defined in TMV_QRDiv.cpp
    template <class T1, int C1, class RT1, class T2>
    void InstQR_SolveInPlace(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, MatrixView<T2> m2);
    template <class T1, int C1, class RT1, class T2>
    void InstQR_SolveTransposeInPlace(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, MatrixView<T2> m2);
    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstQR_Solve(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);
    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstQR_SolveTranspose(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);
    template <class T1, int C1, class RT1, class T2>
    void InstQR_SolveInPlace(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, VectorView<T2> v2);
    template <class T1, int C1, class RT1, class T2>
    void InstQR_SolveTransposeInPlace(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, VectorView<T2> v2);
    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstQR_Solve(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);
    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstQR_SolveTranspose(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    template <int algo, bool trans, int cs, int rs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0)
    // Also used for invalid real/complex combination from the virtual calls.
    template <bool trans, int cs, int rs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<0,trans,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& , const V1& , const Permutation* , int, M2& ) {} 
    };

    // The only difference between the vector and matrix versions of this
    // are the rowRange part.  So here is a little helper class to 
    // do the right thing for vector and matrix objects.
    template <bool vec, class M>
    struct RowRangeHelper;
    template <class V>
    struct RowRangeHelper<true,V>
    {
        typedef typename V::subvector_type rrtype;
        typedef typename V::const_subvector_type const_rrtype;
        static TMV_INLINE rrtype getRowRange(V& v, int i1, int i2) 
        { return v.subVector(i1,i2); }
        template <class V2>
        static TMV_INLINE typename V2::const_subvector_type getRowRange(
            const V2& v, int i1, int i2) 
        { return v.subVector(i1,i2); }
    };
    template <class M>
    struct RowRangeHelper<false,M>
    {
        typedef typename M::rowrange_type rrtype;
        typedef typename M::const_rowrange_type const_rrtype;
        static TMV_INLINE rrtype getRowRange(M& m, int i1, int i2) 
        { return m.rowRange(i1,i2); }
        template <class M2>
        static TMV_INLINE typename M2::const_rowrange_type getRowRange(
            const M2& m, int i1, int i2) 
        { return m.rowRange(i1,i2); }
    };

    // algo 11: Normal case
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<11,false,cs,rs,M1,V1,M2>
    {
        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRSolveInPlace algo 11: trans,cs,rs = "<<
                false<<','<<cs<<','<<rs<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            // m2 = (QRP)^-1 m2
            //    = Pt R^-1 Qt m2
            PackedQ_LDivEq(QR,beta,m2);
            const bool isvec = ShapeTraits<M2::_shape>::vector;
            typedef typename RowRangeHelper<isvec,M2>::rrtype rrtype;
            const int M = QR.colsize();
            rrtype m2a = RowRangeHelper<isvec,M2>::getRowRange(m2,0,N1);
            rrtype m2b = RowRangeHelper<isvec,M2>::getRowRange(m2,N1,M);
            m2b.setZero();
            NoAliasTriLDivEq(m2a,QR.upperTri().subTriMatrix(0,N1));
            if (P) P->inverse().applyOnLeft(m2);
        }
    };
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<11,true,cs,rs,M1,V1,M2>
    {
        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRSolveInPlace algo 11: trans,cs,rs = "<<
                true<<','<<cs<<','<<rs<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            // m2 = (QRP)^-1T m2
            //    = (Pt R^-1 Qt)T m2
            //    = Q* R^-1T P m2
            if (P) P->applyOnLeft(m2);
            const bool isvec = ShapeTraits<M2::_shape>::vector;
            typedef typename RowRangeHelper<isvec,M2>::rrtype rrtype;
            const int M = QR.colsize();
            rrtype m2a = RowRangeHelper<isvec,M2>::getRowRange(m2,0,N1);
            rrtype m2b = RowRangeHelper<isvec,M2>::getRowRange(m2,N1,M);
            m2b.setZero();
            NoAliasTriLDivEq(m2a,QR.upperTri().subTriMatrix(0,N1).transpose());
            PackedQ_MultEq(QR.conjugate(),beta,m2);
        }
    };

    // algo 90: call InstQR_SolveInPlace
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<90,false,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        { InstQR_SolveInPlace(QR.xView(),beta.xView(),P,N1,m2.xView()); }
    };
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<90,true,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
            InstQR_SolveTransposeInPlace(
                QR.xView(),beta.xView(),P,N1,m2.xView()); 
        }
    };

    // algo 95: Turn m2 into vector
    template <bool trans, int cs, int rs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<95,trans,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
            TMVStaticAssert(!ShapeTraits<M2::_shape>::vector);
            typedef typename M2::col_type M2c;
            M2c m2c = m2.col(0);
            QR_SolveInPlace_Helper<-2,trans,cs,rs,M1,V1,M2c>::call(
                QR,beta,P,N1,m2c);
        }
    };

    // algo 97: Conjugate
    template <bool trans, int cs, int rs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<97,trans,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c QRc = QR.conjugate();
            M2c m2c = m2.conjugate();
            QR_SolveInPlace_Helper<-2,trans,cs,rs,M1c,V1,M2c>::call(
                QRc,beta,P,N1,m2c);
        }
    };

    // algo -3: Determine which algorithm to use
    template <bool trans, int cs, int rs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<-3,trans,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
            const bool invalid =
                M1::iscomplex && M2::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                11;
#ifdef PRINTALGO_QR
            std::cout<<"Inline QRSolveInPlace\n";
            std::cout<<"trans,cs,rs = "<<trans<<','<<cs<<','<<rs<<std::endl;
            std::cout<<"QR = "<<TMV_Text(QR)<<std::endl;
            std::cout<<"beta = "<<TMV_Text(beta)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            QR_SolveInPlace_Helper<algo,trans,cs,rs,M1,V1,M2>::call(
                QR,beta,P,N1,m2);
#ifdef PRINTALGO_QR
            //std::cout<<"m2 => "<<m2<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <bool trans, int cs, int rs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<-2,trans,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == UNKNOWN || cs > 16) &&
                (rs == UNKNOWN || rs > 16 || rs == 1) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                M1::_colmajor && 
                Traits<T1>::isinst;
            const bool makevector =
                rs == 1 && !ShapeTraits<M2::_shape>::vector;
            const bool invalid =
                M1::iscomplex && M2::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                makevector ? 95 :
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            QR_SolveInPlace_Helper<algo,trans,cs,rs,M1,V1,M2>::call(
                QR,beta,P,N1,m2);
        }
    };

    template <bool trans, int cs, int rs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<-1,trans,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        { 
            QR_SolveInPlace_Helper<-2,trans,cs,rs,M1,V1,M2>::call(
                QR,beta,P,N1,m2); 
        }
    };

    template <class M1, class V1, class M2>
    static inline void InlineQR_SolveInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const int cs = M2::_colsize;
        const int rs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        QR_SolveInPlace_Helper<-3,false,cs,rs,M1v,V1,M2v>::call(
            QRv,betav,P,N1,m2v);
    }

    template <class M1, class V1, class M2>
    static inline void QR_SolveInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const int cs = M2::_colsize;
        const int rs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        QR_SolveInPlace_Helper<-2,false,cs,rs,M1v,V1,M2v>::call(
            QRv,betav,P,N1,m2v);
    }

    template <class M1, class V1, class V2>
    static inline void InlineQR_SolveInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseVector_Mutable<V2>& v2)
    {
        const int cs = V2::_size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        QR_SolveInPlace_Helper<-3,false,cs,1,M1v,V1,V2v>::call(
            QRv,betav,P,N1,v2v);
    }

    template <class M1, class V1, class V2>
    static inline void QR_SolveInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseVector_Mutable<V2>& v2)
    {
        const int cs = V2::_size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        QR_SolveInPlace_Helper<-1,false,cs,1,M1v,V1,V2v>::call(
            QRv,betav,P,N1,v2v);
    }

    template <class M1, class V1, class M2>
    static inline void InlineQR_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const int cs = M2::_colsize;
        const int rs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        QR_SolveInPlace_Helper<-3,true,cs,rs,M1v,V1,M2v>::call(
            QRv,betav,P,N1,m2v);
    }

    template <class M1, class V1, class M2>
    static inline void QR_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const int cs = M2::_colsize;
        const int rs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        QR_SolveInPlace_Helper<-1,true,cs,rs,M1v,V1,M2v>::call(
            QRv,betav,P,N1,m2v);
    }

    template <class M1, class V1, class V2>
    static inline void InlineQR_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseVector_Mutable<V2>& v2)
    {
        const int cs = V2::_size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        QR_SolveInPlace_Helper<-3,true,cs,1,M1v,V1,V2v>::call(
            QRv,betav,P,N1,v2v);
    }

    template <class M1, class V1, class V2>
    static inline void QR_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseVector_Mutable<V2>& v2)
    {
        const int cs = V2::_size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        QR_SolveInPlace_Helper<-1,true,cs,1,M1v,V1,V2v>::call(
            QRv,betav,P,N1,v2v);
    }

    template <int algo, bool trans, int cs, int rs, int xs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0)
    // Also used for invalid real/complex combination from the virtual calls.
    template <bool trans, int cs, int rs, int xs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<0,trans,cs,rs,xs,M1,V1,M2,M3>
    { 
        static TMV_INLINE void call(
            const M1& , const V1& , const Permutation* , int , 
            const M2& , M3& ) {}
    };

    // Now we need another helper for the difference between matrix
    // and vectors for M2,M3.  Specifically, what object do we copy 
    // m2 to for the temporary.
    template <bool isvec, class M1, class M2, int cs, int rs>
    struct QR_Solve_Aux2;
    template <class M1, class V2, int cs, int rs>
    struct QR_Solve_Aux2<true,M1,V2,cs,rs>
    {
        typedef typename M1::value_type T1;
        typedef typename V2::value_type T2;
        typedef typename Traits2<T1,T2>::type T12;
        typedef typename VCopyHelper<T12,cs,false>::type M2c;
    };
    template <class M1, class M2, int cs, int rs>
    struct QR_Solve_Aux2<false,M1,M2,cs,rs>
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        typedef typename Traits2<T1,T2>::type T12;
        enum { rm = M2::_rowmajor };
        typedef typename MCopyHelper<T12,Rec,cs,rs,rm,false>::type M2c;
    };

    // algo 11: Normal case
    template <int cs, int rs, int xs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<11,false,cs,rs,xs,M1,V1,M2,M3>
    {
        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRSolve algo 11: trans,cs,rs,xs = "<<
                false<<','<<cs<<','<<rs<<','<<xs<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
            //std::cout<<"m3 = "<<m3<<std::endl;
#endif
            const bool isvec = ShapeTraits<M2::_shape>::vector;
            typedef typename RowRangeHelper<isvec,M3>::rrtype rrtype3;
            typedef typename QR_Solve_Aux2<isvec,M1,M2,xs,rs>::M2c M2c;
            const int N = QR.rowsize();
            rrtype3 m3a = RowRangeHelper<isvec,M3>::getRowRange(m3,0,N1);
            rrtype3 m3b = RowRangeHelper<isvec,M3>::getRowRange(m3,N1,N);

            // m3 = (QRP)^-1 m2
            //    = Pt R^-1 Qt m2
            if (QR.isSquare()) {
                m3a = RowRangeHelper<isvec,M3>::getRowRange(m2,0,N1);
                PackedQ_LDivEq(QR,beta,m3);
            } else {
                M2c m2c = m2;
                PackedQ_LDivEq(QR,beta,m2c);
                m3a = RowRangeHelper<isvec,M3>::getRowRange(m2c,0,N1);
            }
            m3b.setZero();
            NoAliasTriLDivEq(m3a,QR.upperTri().subTriMatrix(0,N1));
            if (P) P->inverse().applyOnLeft(m3);
        }
    };
    template <int cs, int rs, int xs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<11,true,cs,rs,xs,M1,V1,M2,M3>
    {
        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRSolve algo 11: trans,cs,rs,xs = "<<
                true<<','<<cs<<','<<rs<<','<<xs<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
            //std::cout<<"m3 = "<<m3<<std::endl;
#endif
            const bool isvec = ShapeTraits<M2::_shape>::vector;
            typedef typename RowRangeHelper<isvec,M3>::rrtype rrtype;
            const int M = QR.colsize();
            const int N = QR.rowsize();
            rrtype m3a = RowRangeHelper<isvec,M3>::getRowRange(m3,0,N1);
            rrtype m3ax = RowRangeHelper<isvec,M3>::getRowRange(m3,0,N);
            rrtype m3b = RowRangeHelper<isvec,M3>::getRowRange(m3,N1,M);

            // m3 = (QRP)^-1T m2
            //    = (Pt R^-1 Qt)T m2
            //    = Q* R^-1T P m2
            m3ax = m2;
            if (P) P->applyOnLeft(m3ax);

            m3b.setZero();
            NoAliasTriLDivEq(m3a,QR.upperTri().subTriMatrix(0,N1).transpose());
            PackedQ_MultEq(QR.conjugate(),beta,m3);
        }
    };

    // algo 90: call InstQR_Solve
    template <int cs, int rs, int xs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<90,false,cs,rs,xs,M1,V1,M2,M3>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        { InstQR_Solve(QR.xView(),beta.xView(),P,N1,m2.xView(),m3.xView()); }
    };
    template <int cs, int rs, int xs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<90,true,cs,rs,xs,M1,V1,M2,M3>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
            InstQR_SolveTranspose(
                QR.xView(),beta.xView(),P,N1,m2.xView(),m3.xView()); 
        }
    };

    // algo 95: Turn QR,m3 into vector
    template <bool trans, int cs, int rs, int xs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<95,trans,cs,rs,xs,M1,V1,M2,M3>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
            TMVStaticAssert(!ShapeTraits<M2::_shape>::vector);
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;
            M2c m2c = m2.col(0);
            M3c m3c = m3.col(0);
            QR_Solve_Helper<-2,trans,cs,rs,xs,M1,V1,M2c,M3c>::call(
                QR,beta,P,N1,m2c,m3c);
        }
    };

    // algo 97: Conjugate
    template <bool trans, int cs, int rs, int xs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<97,trans,cs,rs,xs,M1,V1,M2,M3>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c QRc = QR.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            QR_Solve_Helper<-2,trans,cs,rs,xs,M1c,V1,M2c,M3c>::call(
                QRc,beta,P,N1,m2c,m3c);
        }
    };

    // algo -3: Determine which algorithm to use
    template <bool trans, int cs, int rs, int xs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<-3,trans,cs,rs,xs,M1,V1,M2,M3>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
            const bool invalid =
                (M1::iscomplex || M2::iscomplex) && M3::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                11;
#ifdef PRINTALGO_QR
            std::cout<<"Inline QRSolve\n";
            std::cout<<"trans,cs,rs,xs = "<<trans<<','<<cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"QR = "<<TMV_Text(QR)<<std::endl;
            std::cout<<"beta = "<<TMV_Text(beta)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
            //std::cout<<"m3 = "<<m3<<std::endl;
#endif
            QR_Solve_Helper<algo,trans,cs,rs,xs,M1,V1,M2,M3>::call(
                QR,beta,P,N1,m2,m3);
#ifdef PRINTALGO_QR
            //std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <bool trans, int cs, int rs, int xs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<-2,trans,cs,rs,xs,M1,V1,M2,M3>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == UNKNOWN || cs > 16) &&
                (rs == UNKNOWN || rs > 16 || rs == 1) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                M1::_colmajor && 
                Traits<T1>::isinst;
            const bool makevector =
                rs == 1 && !ShapeTraits<M2::_shape>::vector;
            const bool invalid =
                (M1::iscomplex || M2::iscomplex) && M3::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                makevector ? 95 :
                M3::_conj ? 97 :
                inst ? 90 :
                -3;
            QR_Solve_Helper<algo,trans,cs,rs,xs,M1,V1,M2,M3>::call(
                QR,beta,P,N1,m2,m3);
        }
    };

    template <bool trans, int cs, int rs, int xs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<-1,trans,cs,rs,xs,M1,V1,M2,M3>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
            QR_Solve_Helper<-2,trans,cs,rs,xs,M1,V1,M2,M3>::call(
                QR,beta,P,N1,m2,m3); 
        }
    };

    template <class M1, class V1, class M2, class M3>
    static inline void InlineQR_Solve(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(QR.rowsize() == beta.size());
        TMVAssert(QR.colsize() == m2.colsize());
        TMVAssert(QR.rowsize() == m3.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());
        if (P) TMVAssert(QR.rowsize() == P->size());
        TMVAssert(N1 <= int(QR.rowsize()));

        // cs = m3.colsize
        // rs = m3.rowsize
        // xs = "extra size" = QR.colsize, m2.colsize
        const int cs = Sizes<M3::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_colsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        QR_Solve_Helper<-3,false,cs,rs,xs,M1v,V1,M2v,M3v>::call(
            QRv,betav,P,N1,m2v,m3v);
    }

    template <class M1, class V1, class M2, class M3>
    static inline void QR_Solve(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(QR.rowsize() == beta.size());
        TMVAssert(QR.colsize() == m2.colsize());
        TMVAssert(QR.rowsize() == m3.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());
        if (P) TMVAssert(QR.rowsize() == P->size());
        TMVAssert(N1 <= int(QR.rowsize()));

        const int cs = Sizes<M3::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_colsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        QR_Solve_Helper<-2,false,cs,rs,xs,M1v,V1,M2v,M3v>::call(
            QRv,betav,P,N1,m2v,m3v);
    }

    template <class M1, class V1, class V2, class V3>
    static inline void InlineQR_Solve(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, 
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V2::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V3::_size>::same));
        TMVAssert(QR.rowsize() == beta.size());
        TMVAssert(QR.colsize() == v2.size());
        TMVAssert(QR.rowsize() == v3.size());
        if (P) TMVAssert(QR.rowsize() == P->size());
        TMVAssert(N1 <= int(QR.rowsize()));

        const int cs = Sizes<V3::_size,M1::_rowsize>::size;
        const int xs = Sizes<M1::_colsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        QR_Solve_Helper<-3,false,cs,1,xs,M1v,V1v,V2v,V3v>::call(
            QRv,betav,P,N1,v2v,v3v);
    }

    template <class M1, class V1, class V2, class V3>
    static inline void QR_Solve(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, 
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V2::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V3::_size>::same));
        TMVAssert(QR.rowsize() == beta.size());
        TMVAssert(QR.colsize() == v2.size());
        TMVAssert(QR.rowsize() == v3.size());
        if (P) TMVAssert(QR.rowsize() == P->size());
        TMVAssert(N1 <= int(QR.rowsize()));

        const int cs = Sizes<V3::_size,M1::_rowsize>::size;
        const int xs = Sizes<M1::_colsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        QR_Solve_Helper<-1,false,cs,1,xs,M1v,V1v,V2v,V3v>::call(
            QRv,betav,P,N1,v2v,v3v);
    }

    template <class M1, class V1, class M2, class M3>
    static inline void InlineQR_SolveTranspose(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(QR.rowsize() == beta.size());
        TMVAssert(QR.colsize() == m3.colsize());
        TMVAssert(QR.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());
        if (P) TMVAssert(QR.rowsize() == P->size());
        TMVAssert(N1 <= int(QR.rowsize()));

        const int cs = Sizes<M3::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_colsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        QR_Solve_Helper<-3,true,cs,rs,xs,M1v,V1v,M2v,M3v>::call(
            QRv,betav,P,N1,m2v,m3v);
    }

    template <class M1, class V1, class M2, class M3>
    static inline void QR_SolveTranspose(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(QR.rowsize() == beta.size());
        TMVAssert(QR.colsize() == m3.colsize());
        TMVAssert(QR.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());
        if (P) TMVAssert(QR.rowsize() == P->size());
        TMVAssert(N1 <= int(QR.rowsize()));

        const int cs = Sizes<M3::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_colsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        QR_Solve_Helper<-1,true,cs,rs,xs,M1v,V1v,M2v,M3v>::call(
            QRv,betav,P,N1,m2v,m3v);
    }

    template <class M1, class V1, class V2, class V3>
    static inline void InlineQR_SolveTranspose(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, 
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V3::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V2::_size>::same));
        TMVAssert(QR.rowsize() == beta.size());
        TMVAssert(QR.colsize() == v3.size());
        TMVAssert(QR.rowsize() == v2.size());
        if (P) TMVAssert(QR.rowsize() == P->size());
        TMVAssert(N1 <= int(QR.rowsize()));

        const int cs = Sizes<V3::_size,M1::_rowsize>::size;
        const int xs = Sizes<M1::_colsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        QR_Solve_Helper<-3,true,cs,1,xs,M1v,V1v,V2v,V3v>::call(
            QRv,betav,P,N1,v2v,v3v);
    }

    template <class M1, class V1, class V2, class V3>
    static inline void QR_SolveTranspose(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, 
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V3::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V2::_size>::same));
        TMVAssert(QR.rowsize() == beta.size());
        TMVAssert(QR.colsize() == v3.size());
        TMVAssert(QR.rowsize() == v2.size());
        if (P) TMVAssert(QR.rowsize() == P->size());
        TMVAssert(N1 <= int(QR.rowsize()));

        const int cs = Sizes<V3::_size,M1::_rowsize>::size;
        const int xs = Sizes<M1::_colsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        QR_Solve_Helper<-1,true,cs,1,xs,M1v,V1v,V2v,V3v>::call(
            QRv,betav,P,N1,v2v,v3v);
    }

} // namespace tmv

#endif

