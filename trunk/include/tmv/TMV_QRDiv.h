

#ifndef TMV_QRDiv_H
#define TMV_QRDiv_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_DivVM_Funcs.h"
#include "TMV_DivMM_Funcs.h"
#include "TMV_MultMM_Funcs.h"
#include "TMV_Permutation.h"
#include "TMV_PackedQ.h"

//#define PRINTALGO_QR

#ifdef PRINTALGO_QR
#include <iostream>
#include "TMV_ProdMV.h"
#include "TMV_ProdMM.h"
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#include "TMV_AddMM.h"
#include "TMV_AddVV.h"
#include "TMV_MultMM.h"
#include "TMV_MultMV.h"
#include "TMV_MultUV.h"
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

    template <int algo, bool trans, int cs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper;

    // algo 0: Trivial, nothing to do (M == 0)
    // Also used for invalid real/complex combination from the virtual calls.
    template <bool trans, int cs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<0,trans,cs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& , const V1& , const Permutation* , int, M2& ) {} 
    };

    // algo 11: Normal case
    template <int cs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<11,false,cs,M1,V1,M2>
    {
        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRSolveInPlace algo 11: trans,cs = "<<
                false<<','<<cs<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            // m2 = (QRP)^-1 m2
            //    = Pt R^-1 Qt m2
            PackedQ_LDivEq(QR,beta,m2);
            typedef typename M2::rowrange_type::noalias_type M2r;
            const int M = QR.colsize();
            M2r m2a = m2.rowRange(0,N1).noAlias();
            M2r m2b = m2.rowRange(N1,M).noAlias();
            m2b.setZero();
            TriLDivEq(m2a,QR.upperTri().subTriMatrix(0,N1));
            if (P) P->inverse().applyOnLeft(m2);
        }
    };
    template <int cs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<11,true,cs,M1,V1,M2>
    {
        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRSolveInPlace algo 11: trans,cs = "<<
                true<<','<<cs<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            // m2 = (QRP)^-1T m2
            //    = (Pt R^-1 Qt)T m2
            //    = Q* R^-1T P m2
            if (P) P->applyOnLeft(m2);
            typedef typename M2::rowrange_type::noalias_type M2r;
            const int M = QR.colsize();
            M2r m2a = m2.rowRange(0,N1).noAlias();
            M2r m2b = m2.rowRange(N1,M).noAlias();
            m2b.setZero();
            TriLDivEq(m2a,QR.upperTri().subTriMatrix(0,N1).transpose());
            PackedQ_MultEq(QR.conjugate(),beta,m2);
        }
    };

    // algo 12: Normal case: M2 is a vector
    template <int cs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<12,false,cs,M1,V1,M2>
    {
        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& v2)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRSolveInPlace algo 12: trans,cs = "<<
                false<<','<<cs<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"v2 = "<<v2<<std::endl;
#endif
            // v2 = (QRP)^-1 v2
            //    = Pt R^-1 Qt v2
            PackedQ_LDivEq(QR,beta,v2);
            typedef typename M2::subvector_type::noalias_type M2s;
            const int M = QR.colsize();
            M2s v2a = v2.subVector(0,N1).noAlias();
            M2s v2b = v2.subVector(N1,M).noAlias();
            v2b.setZero();
            TriLDivEq(v2a,QR.upperTri().subTriMatrix(0,N1));
            if (P) P->inverse().applyOnLeft(v2);
        }
    };
    template <int cs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<12,true,cs,M1,V1,M2>
    {
        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRSolveInPlace algo 12: trans,cs = "<<
                true<<','<<cs<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            // m2 = (QRP)^-1T m2
            //    = (Pt R^-1 Qt)T m2
            //    = Q* R^-1T P m2
            if (P) P->applyOnLeft(m2);
            typedef typename M2::subvector_type::noalias_type M2s;
            const int M = QR.colsize();
            M2s m2a = m2.subVector(0,N1).noAlias();
            M2s m2b = m2.subVector(N1,M).noAlias();
            m2b.setZero();
            TriLDivEq(m2a,QR.upperTri().subTriMatrix(0,N1).transpose());
            PackedQ_MultEq(QR.conjugate(),beta,m2);
        }
    };

    // algo 90: call InstQR_SolveInPlace
    template <int cs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<90,false,cs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        { InstQR_SolveInPlace(QR.xView(),beta.xView(),P,N1,m2.xView()); }
    };
    template <int cs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<90,true,cs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
            InstQR_SolveTransposeInPlace(
                QR.xView(),beta.xView(),P,N1,m2.xView()); 
        }
    };

    // algo 97: Conjugate
    template <bool trans, int cs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<97,trans,cs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c QRc = QR.conjugate();
            M2c m2c = m2.conjugate();
            QR_SolveInPlace_Helper<-2,trans,cs,M1c,V1,M2c>::call(
                QRc,beta,P,N1,m2c);
        }
    };

    // algo -3: Determine which algorithm to use
    template <bool trans, int cs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<-3,trans,cs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
            const bool invalid =
                M1::iscomplex && M2::isreal;
            const int algo = 
                cs == 0 || invalid ? 0 : 
                ShapeTraits<M2::_shape>::vector ? 12 : 11;
#ifdef PRINTALGO_QR
            std::cout<<"Inline QRSolveInPlace\n";
            std::cout<<"trans,cs = "<<trans<<','<<cs<<std::endl;
            std::cout<<"QR = "<<TMV_Text(QR)<<std::endl;
            std::cout<<"beta = "<<TMV_Text(beta)<<std::endl;
            std::cout<<"N1 = "<<N1<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //if (P) std::cout<<"P = "<<*P<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            QR_SolveInPlace_Helper<algo,trans,cs,M1,V1,M2>::call(
                QR,beta,P,N1,m2);
#ifdef PRINTALGO_QR
            //std::cout<<"m2 => "<<m2<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <bool trans, int cs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<-2,trans,cs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const bool invalid =
                M1::iscomplex && M2::isreal;
            const int algo = 
                cs == 0 || invalid ? 0 : 
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            QR_SolveInPlace_Helper<algo,trans,cs,M1,V1,M2>::call(
                QR,beta,P,N1,m2);
        }
    };

    template <bool trans, int cs, class M1, class V1, class M2>
    struct QR_SolveInPlace_Helper<-1,trans,cs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
        { 
            QR_SolveInPlace_Helper<-2,trans,cs,M1,V1,M2>::call(
                QR,beta,P,N1,m2); 
        }
    };

    template <class M1, class V1, class M2>
    inline void InlineQR_SolveInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        const int cs1 = Sizes<M1::_colsize,M1::_rowsize>::size;
        const int cs2 = Sizes<V1::_size,M2::_colsize>::size;
        const int cs = Sizes<cs1,cs2>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        QR_SolveInPlace_Helper<-3,false,cs,M1v,V1,M2v>::call(
            QRv,betav,P,N1,m2v);
    }

    template <class M1, class V1, class M2>
    inline void QR_SolveInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        const int cs1 = Sizes<M1::_colsize,M1::_rowsize>::size;
        const int cs2 = Sizes<V1::_size,M2::_colsize>::size;
        const int cs = Sizes<cs1,cs2>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        QR_SolveInPlace_Helper<-2,false,cs,M1v,V1,M2v>::call(
            QRv,betav,P,N1,m2v);
    }

    template <class M1, class V1, class V2>
    inline void InlineQR_SolveInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V1::_size>::same));
        const int cs1 = Sizes<M1::_colsize,M1::_rowsize>::size;
        const int cs2 = Sizes<V1::_size,V2::_size>::size;
        const int cs = Sizes<cs1,cs2>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        QR_SolveInPlace_Helper<-3,false,cs,M1v,V1,V2v>::call(
            QRv,betav,P,N1,v2v);
    }

    template <class M1, class V1, class V2>
    inline void QR_SolveInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V1::_size>::same));
        const int cs1 = Sizes<M1::_colsize,M1::_rowsize>::size;
        const int cs2 = Sizes<V1::_size,V2::_size>::size;
        const int cs = Sizes<cs1,cs2>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        QR_SolveInPlace_Helper<-2,false,cs,M1v,V1,V2v>::call(
            QRv,betav,P,N1,v2v);
    }

    template <class M1, class V1, class M2>
    inline void InlineQR_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        const int cs1 = Sizes<M1::_colsize,M1::_rowsize>::size;
        const int cs2 = Sizes<V1::_size,M2::_colsize>::size;
        const int cs = Sizes<cs1,cs2>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        QR_SolveInPlace_Helper<-3,true,cs,M1v,V1,M2v>::call(
            QRv,betav,P,N1,m2v);
    }

    template <class M1, class V1, class M2>
    inline void QR_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        const int cs1 = Sizes<M1::_colsize,M1::_rowsize>::size;
        const int cs2 = Sizes<V1::_size,M2::_colsize>::size;
        const int cs = Sizes<cs1,cs2>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        QR_SolveInPlace_Helper<-2,true,cs,M1v,V1,M2v>::call(
            QRv,betav,P,N1,m2v);
    }

    template <class M1, class V1, class V2>
    inline void InlineQR_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V1::_size>::same));
        const int cs1 = Sizes<M1::_colsize,M1::_rowsize>::size;
        const int cs2 = Sizes<V1::_size,V2::_size>::size;
        const int cs = Sizes<cs1,cs2>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        QR_SolveInPlace_Helper<-3,true,cs,M1v,V1,V2v>::call(
            QRv,betav,P,N1,v2v);
    }

    template <class M1, class V1, class V2>
    inline void QR_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,V1::_size>::same));
        const int cs1 = Sizes<M1::_colsize,M1::_rowsize>::size;
        const int cs2 = Sizes<V1::_size,V2::_size>::size;
        const int cs = Sizes<cs1,cs2>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        QR_SolveInPlace_Helper<-2,true,cs,M1v,V1,V2v>::call(
            QRv,betav,P,N1,v2v);
    }

    template <int algo, bool trans, int cs, int rs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0)
    // Also used for invalid real/complex combination from the virtual calls.
    template <bool trans, int cs, int rs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<0,trans,cs,rs,M1,V1,M2,M3>
    { 
        static TMV_INLINE void call(
            const M1& , const V1& , const Permutation* , int , 
            const M2& , M3& ) {}
    };

    // algo 11: Normal case
    template <int cs, int rs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<11,false,cs,rs,M1,V1,M2,M3>
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        typedef typename Traits2<T1,T2>::type T12;
        enum { xs = M2::_rowsize };
        enum { A = M2::_rowmajor ? RowMajor : ColMajor };
        typedef typename MCopyHelper<T12,Rec,rs,xs,A>::type M2c;
        typedef typename M3::rowrange_type::noalias_type M3r;

        template <int square, int dummy>
        struct Helper;

        template <int dummy>
        struct Helper<0,dummy> // not square
        {
            static void call(
                const M1& QR, const V1& beta, const M2& m2, M3& m3)
            {
                M2c m2c = m2;
                PackedQ_LDivEq(QR,beta,m2c);
                m3 = m2c.rowRange(0,m3.colsize());
            }
        };
        template <int dummy>
        struct Helper<1,dummy> // square
        {
            static void call(
                const M1& QR, const V1& beta, const M2& m2, M3& m3)
            {
                m3 = m2;
                PackedQ_LDivEq(QR,beta,m3);
            }
        };
        template <int dummy>
        struct Helper<2,dummy> // maybe square
        {
            static void call(
                const M1& QR, const V1& beta, const M2& m2, M3& m3)
            {
                if (QR.isSquare()) {
                    m3 = m2;
                    PackedQ_LDivEq(QR,beta,m3);
                } else {
                    M2c m2c = m2;
                    PackedQ_LDivEq(QR,beta,m2c);
                    m3 = m2c.rowRange(0,m3.colsize());
                }
            }
        };


        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRSolve algo 11: trans,cs,rs = "<<
                false<<','<<cs<<','<<rs<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
            //std::cout<<"m3 = "<<m3<<std::endl;
#endif
            const int N = QR.rowsize();

            M3r m3a = m3.rowRange(0,N1).noAlias();
            M3r m3b = m3.rowRange(N1,N).noAlias();

            // m3 = (QRP)^-1 m2
            //    = Pt R^-1 Qt m2
            const int square =
                cs == Unknown || rs == Unknown ? 2 :
                cs == rs ? 1 : 0;
            //std::cout<<"square = "<<square<<std::endl;
            Helper<square,1>::call(QR,beta,m2,m3);
            //std::cout<<"m3 => "<<m3<<std::endl;
            m3b.setZero();
            //std::cout<<"m3 => "<<m3<<std::endl;
            TriLDivEq(m3a,QR.upperTri().subTriMatrix(0,N1));
            //std::cout<<"m3 => "<<m3<<std::endl;
            if (P) P->inverse().applyOnLeft(m3);
            //std::cout<<"m3 => "<<m3<<std::endl;
        }
    };
    template <int cs, int rs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<11,true,cs,rs,M1,V1,M2,M3>
    {
        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRSolve algo 11: trans,cs,rs = "<<
                true<<','<<cs<<','<<rs<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
            //std::cout<<"m3 = "<<m3<<std::endl;
#endif
            typedef typename M3::rowrange_type::noalias_type M3r;
            const int M = QR.colsize();
            const int N = QR.rowsize();
            M3r m3a = m3.rowRange(0,N1).noAlias();
            M3r m3ax = m3.rowRange(0,N).noAlias();
            M3r m3b = m3.rowRange(N1,M).noAlias();

            // m3 = (QRP)^-1T m2
            //    = (Pt R^-1 Qt)T m2
            //    = Q* R^-1T P m2
            m3ax = m2;
            //std::cout<<"m3 => "<<m3<<std::endl;
            if (P) P->applyOnLeft(m3ax);
            //std::cout<<"m3 => "<<m3<<std::endl;

            m3b.setZero();
            //std::cout<<"m3 => "<<m3<<std::endl;
            TriLDivEq(m3a,QR.upperTri().subTriMatrix(0,N1).transpose());
            //std::cout<<"m3 => "<<m3<<std::endl;
            PackedQ_MultEq(QR.conjugate(),beta,m3);
            //std::cout<<"m3 => "<<m3<<std::endl;
        }
    };

    // algo 12: Normal case - M2, M3 are vectors
    template <int cs, int rs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<12,false,cs,rs,M1,V1,M2,M3>
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        typedef typename Traits2<T1,T2>::type T12;
        typedef typename VCopyHelper<T12,rs>::type M2c;
        typedef typename M3::subvector_type::noalias_type M3s;

        template <int square, int dummy>
        struct Helper;

        template <int dummy>
        struct Helper<0,dummy> // not square
        {
            static void call(
                const M1& QR, const V1& beta, const M2& v2, M3& v3)
            {
                //std::cout<<"Helper 0\n";
                M2c v2c = v2;
                //std::cout<<"v2c = "<<v2c<<std::endl;
                PackedQ_LDivEq(QR,beta,v2c);
                //std::cout<<"v2c => "<<v2c<<std::endl;
                v3 = v2c.subVector(0,v3.size());
                //std::cout<<"v3 => "<<v3<<std::endl;
            }
        };
        template <int dummy>
        struct Helper<1,dummy> // square
        {
            static void call(
                const M1& QR, const V1& beta, const M2& v2, M3& v3)
            {
                //std::cout<<"Helper 1\n";
                v3 = v2;
                //std::cout<<"v3 = "<<v3<<std::endl;
                PackedQ_LDivEq(QR,beta,v3);
                //std::cout<<"v3 => "<<v3<<std::endl;
            }
        };
        template <int dummy>
        struct Helper<2,dummy> // maybe square
        {
            static void call(
                const M1& QR, const V1& beta, const M2& v2, M3& v3)
            {
                //std::cout<<"Helper 2\n";
                if (QR.isSquare()) {
                    v3 = v2;
                    //std::cout<<"v3 = "<<v3<<std::endl;
                    PackedQ_LDivEq(QR,beta,v3);
                    //std::cout<<"v3 => "<<v3<<std::endl;
                } else {
                    M2c v2c = v2;
                    //std::cout<<"v2c = "<<v2c<<std::endl;
                    PackedQ_LDivEq(QR,beta,v2c);
                    //std::cout<<"v2c => "<<v2c<<std::endl;
                    v3 = v2c.subVector(0,v3.size());
                    //std::cout<<"v3 => "<<v3<<std::endl;
                }
            }
        };


        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& v2, M3& v3)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRSolve algo 12: trans,cs,rs = "<<
                false<<','<<cs<<','<<rs<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"v2 = "<<v2<<std::endl;
            //std::cout<<"v3 = "<<v3<<std::endl;
#endif
            const int N = QR.rowsize();

            M3s v3a = v3.subVector(0,N1).noAlias();
            M3s v3b = v3.subVector(N1,N).noAlias();

            // v3 = (QRP)^-1 v2
            //    = Pt R^-1 Qt v2
            const int square =
                cs == Unknown || rs == Unknown ? 2 :
                cs == rs ? 1 : 0;
            //std::cout<<"square = "<<square<<std::endl;
            Helper<square,1>::call(QR,beta,v2,v3);
            //std::cout<<"After /= Q: v3 => "<<v3<<std::endl;
            //typename M1::copy_type Q = QR;
            //UnpackQ(Q,beta);
            //std::cout<<"Q = "<<Q<<std::endl;
            //std::cout<<"Q*v3 = "<<Q*v3<<std::endl;
            v3b.setZero();
            //std::cout<<"v3 => "<<v3<<std::endl;
            //std::cout<<"Q*v3 = "<<Q*v3<<std::endl;
            //std::cout<<"N1 = "<<N1<<std::endl;
            //std::cout<<"v3a = "<<v3a<<std::endl;
            TriLDivEq(v3a,QR.upperTri().subTriMatrix(0,N1));
            //std::cout<<"v3a => "<<v3a<<std::endl;
            //std::cout<<"after /= R v3 => "<<v3<<std::endl;
            //std::cout<<"R = "<<QR.upperTri()<<std::endl;
            //std::cout<<"R * v3 = "<<QR.upperTri()*v3<<std::endl;
            if (P) P->inverse().applyOnLeft(v3);
            //std::cout<<"after /= P v3 => "<<v3<<std::endl;
        }
    };
    template <int cs, int rs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<12,true,cs,rs,M1,V1,M2,M3>
    {
        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& v2, M3& v3)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRSolve algo 12: trans,cs,rs = "<<
                true<<','<<cs<<','<<rs<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //std::cout<<"v2 = "<<v2<<std::endl;
            //std::cout<<"v3 = "<<v3<<std::endl;
#endif
            typedef typename M3::subvector_type::noalias_type M3s;
            const int M = QR.colsize();
            const int N = QR.rowsize();
            M3s v3a = v3.subVector(0,N1).noAlias();
            M3s v3ax = v3.subVector(0,N).noAlias();
            M3s v3b = v3.subVector(N1,M).noAlias();

            // v3 = (QRP)^-1T v2
            //    = (Pt R^-1 Qt)T v2
            //    = Q* R^-1T P v2
            v3ax = v2;
            //std::cout<<"v3 => "<<v3<<std::endl;
            if (P) P->applyOnLeft(v3ax);
            //std::cout<<"after /= PT: v3 => "<<v3<<std::endl;

            v3b.setZero();
            //std::cout<<"v3 => "<<v3<<std::endl;
            TriLDivEq(v3a,QR.upperTri().subTriMatrix(0,N1).transpose());
            //std::cout<<"After /= RT: v3 => "<<v3<<std::endl;
            PackedQ_MultEq(QR.conjugate(),beta,v3);
            //std::cout<<"After /= QT: v3 => "<<v3<<std::endl;
        }
    };

    // algo 90: call InstQR_Solve
    template <int cs, int rs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<90,false,cs,rs,M1,V1,M2,M3>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        { InstQR_Solve(QR.xView(),beta.xView(),P,N1,m2.xView(),m3.xView()); }
    };
    template <int cs, int rs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<90,true,cs,rs,M1,V1,M2,M3>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
            InstQR_SolveTranspose(
                QR.xView(),beta.xView(),P,N1,m2.xView(),m3.xView()); 
        }
    };

    // algo 97: Conjugate
    template <bool trans, int cs, int rs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<97,trans,cs,rs,M1,V1,M2,M3>
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
            QR_Solve_Helper<-2,trans,cs,rs,M1c,V1,M2c,M3c>::call(
                QRc,beta,P,N1,m2c,m3c);
        }
    };

    // algo -3: Determine which algorithm to use
    template <bool trans, int cs, int rs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<-3,trans,cs,rs,M1,V1,M2,M3>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
            TMVStaticAssert((
                    ShapeTraits<M2::_shape>::vector == 
                    ShapeTraits<M3::_shape>::vector));

            const bool invalid =
                (M1::iscomplex || M2::iscomplex) && M3::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                ShapeTraits<M3::_shape>::vector ? 12 : 11;
#ifdef PRINTALGO_QR
            std::cout<<"Inline QRSolve\n";
            std::cout<<"trans,cs,rs = "<<trans<<','<<cs<<','<<rs<<std::endl;
            std::cout<<"QR = "<<TMV_Text(QR)<<std::endl;
            std::cout<<"beta = "<<TMV_Text(beta)<<std::endl;
            std::cout<<"N1 = "<<N1<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            //std::cout<<"QR = "<<QR<<std::endl;
            //std::cout<<"beta = "<<beta<<std::endl;
            //if (P) std::cout<<"P = "<<*P<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
            //std::cout<<"m3 = "<<m3<<std::endl;
#endif
            QR_Solve_Helper<algo,trans,cs,rs,M1,V1,M2,M3>::call(
                QR,beta,P,N1,m2,m3);
#ifdef PRINTALGO_QR
            //std::cout<<"m3 => "<<m3<<std::endl;
            //typename M1::copy_type A = QR;
            //UnpackQ(A,beta);
            //A *= QR.upperTri();
            //if (P) P->applyOnRight(A);
            //std::cout<<"A = "<<A<<std::endl;
            //std::cout<<"A*m3 = "<<A*m3<<std::endl;
            //std::cout<<"A*m3-m2 = "<<A*m3-m2<<std::endl;
            //std::cout<<"AtA*m3-Atm2 = "<<A.adjoint()*A*m3-A.adjoint()*m2<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <bool trans, int cs, int rs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<-2,trans,cs,rs,M1,V1,M2,M3>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const bool invalid =
                (M1::iscomplex || M2::iscomplex) && M3::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                M3::_conj ? 97 :
                inst ? 90 :
                -3;
            QR_Solve_Helper<algo,trans,cs,rs,M1,V1,M2,M3>::call(
                QR,beta,P,N1,m2,m3);
        }
    };

    template <bool trans, int cs, int rs, class M1, class V1, class M2, class M3>
    struct QR_Solve_Helper<-1,trans,cs,rs,M1,V1,M2,M3>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, 
            const M2& m2, M3& m3)
        {
            QR_Solve_Helper<-2,trans,cs,rs,M1,V1,M2,M3>::call(
                QR,beta,P,N1,m2,m3); 
        }
    };

    template <class M1, class V1, class M2, class M3>
    inline void InlineQR_Solve(
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
        TMVAssert(N1 <= QR.rowsize());

        // cs = QR.colsize
        // rs = QR.rowsize
        const int cs = Sizes<M3::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M1::_colsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        QR_Solve_Helper<-3,false,cs,rs,M1v,V1,M2v,M3v>::call(
            QRv,betav,P,N1,m2v,m3v);
    }

    template <class M1, class V1, class M2, class M3>
    inline void QR_Solve(
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
        TMVAssert(N1 <= QR.rowsize());

        const int cs = Sizes<M3::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M1::_colsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        QR_Solve_Helper<-2,false,cs,rs,M1v,V1,M2v,M3v>::call(
            QRv,betav,P,N1,m2v,m3v);
    }

    template <class M1, class V1, class V2, class V3>
    inline void InlineQR_Solve(
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
        TMVAssert(N1 <= QR.rowsize());

        const int cs = Sizes<V3::_size,M1::_rowsize>::size;
        const int rs = Sizes<M1::_colsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        QR_Solve_Helper<-3,false,cs,rs,M1v,V1,V2v,V3v>::call(
            QRv,betav,P,N1,v2v,v3v);
    }

    template <class M1, class V1, class V2, class V3>
    inline void QR_Solve(
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
        TMVAssert(N1 <= QR.rowsize());

        const int cs = Sizes<V3::_size,M1::_rowsize>::size;
        const int rs = Sizes<M1::_colsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        QR_Solve_Helper<-2,false,cs,rs,M1v,V1,V2v,V3v>::call(
            QRv,betav,P,N1,v2v,v3v);
    }

    template <class M1, class V1, class M2, class M3>
    inline void InlineQR_SolveTranspose(
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
        TMVAssert(N1 <= QR.rowsize());

        const int cs = Sizes<M3::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M1::_colsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        QR_Solve_Helper<-3,true,cs,rs,M1v,V1v,M2v,M3v>::call(
            QRv,betav,P,N1,m2v,m3v);
    }

    template <class M1, class V1, class M2, class M3>
    inline void QR_SolveTranspose(
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
        TMVAssert(N1 <= QR.rowsize());

        const int cs = Sizes<M3::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M1::_colsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        QR_Solve_Helper<-2,true,cs,rs,M1v,V1v,M2v,M3v>::call(
            QRv,betav,P,N1,m2v,m3v);
    }

    template <class M1, class V1, class V2, class V3>
    inline void InlineQR_SolveTranspose(
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
        TMVAssert(N1 <= QR.rowsize());

        const int cs = Sizes<V3::_size,M1::_rowsize>::size;
        const int rs = Sizes<M1::_colsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        QR_Solve_Helper<-3,true,cs,rs,M1v,V1v,V2v,V3v>::call(
            QRv,betav,P,N1,v2v,v3v);
    }

    template <class M1, class V1, class V2, class V3>
    inline void QR_SolveTranspose(
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
        TMVAssert(N1 <= QR.rowsize());

        const int cs = Sizes<V3::_size,M1::_rowsize>::size;
        const int rs = Sizes<M1::_colsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        QR_Solve_Helper<-2,true,cs,rs,M1v,V1v,V2v,V3v>::call(
            QRv,betav,P,N1,v2v,v3v);
    }

} // namespace tmv

#endif

