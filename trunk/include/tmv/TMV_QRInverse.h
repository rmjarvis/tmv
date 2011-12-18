

#ifndef TMV_QRInverse_H
#define TMV_QRInverse_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Permutation.h"
#include "TMV_MultMM_Funcs.h"
#include "TMV_PackedQ.h"

#ifdef PRINTALGO_QR
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_TriMatrixIO.h"
#endif

namespace tmv {

    // Defined in TMV_QRInverse.cpp
    template <class T1, int C1, class RT1, class T2>
    void InstQR_Inverse(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, MatrixView<T2> minv);

    template <class T1, int C1, class RT1, class T2>
    void InstQR_InverseATA(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, MatrixView<T2> ata);



    // Note: cs,rs refer to M1, not M2 (which is reverse of M1)
    template <int algo, int cs, int rs, class M1, class V1, class M2>
    struct QR_Inverse_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0)
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_Inverse_Helper<0,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& , const V1& , const Permutation* , int, M2& ) {} 
    };

    // algo 11: Normal case
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_Inverse_Helper<11,cs,rs,M1,V1,M2>
    {
        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& minv)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRInverse algo 11: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            // m2 = (QRP)^-1
            //    = Pt R^-1 Qt
#ifdef XDEBUG_QR
            typedef typename M1::value_type T;
            Matrix<T> Q1 = QR;
            UnpackQ(Q1,beta);
            Matrix<T> R1 = QR.upperTri();
            Matrix<T> A1 = Q1*R1;
            if (P) P->applyOnRight(A1);
#endif
            minv.setZero();
            const int N = rs == TMV_UNKNOWN ? QR.rowsize() : rs;
            typename M2::colrange_type::uppertri_type R = 
                minv.colRange(0,N1).upperTri();
            R = QR.upperTri().subTriMatrix(0,N1);
            R.invertSelf();
            if (P) minv.colRange(0,N) /= *P;
            minv %= PackedQ<M1,V1>(QR,beta);

#ifdef XDEBUG_QR
            if (Norm(minv*A1-T(1)) > 1.e-3*Norm(A1)) {
                std::cout<<"QR = "<<QR<<std::endl;
                std::cout<<"Q = "<<Q1<<std::endl;
                std::cout<<"R = "<<R1<<std::endl;
                if (P) std::cout<<"P = "<<*P<<std::endl;
                std::cout<<"A = QRP = "<<A1<<std::endl;
                std::cout<<"minv = "<<minv<<std::endl;
                std::cout<<"minv * A = "<<(minv*A1)<<std::endl;
                std::cout<<"A * minv = "<<(A1*minv)<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo 12: Alternate method: Unpack Q first.
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_Inverse_Helper<12,cs,rs,M1,V1,M2>
    {
        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& minv)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRInverse algo 12: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            // m2 = (QRP)^-1
            //    = Pt R^-1 Qt
            // m2T = Q* RT^-1 P
#ifdef XDEBUG_QR
            typedef typename M1::value_type T;
            Matrix<T> Q1 = QR;
            UnpackQ(Q1,beta);
            Matrix<T> R1 = QR.upperTri();
            Matrix<T> A1 = Q1*R1;
            if (P) P->applyOnRight(A1);
#endif
            const int N = rs == TMV_UNKNOWN ? QR.rowsize() : rs;
            typename M2::transpose_type m2t = minv.transpose();
            m2t = QR.conjugate();
            UnpackQ(m2t,beta);
            m2t.colRange(0,N1) %= QR.upperTri().subTriMatrix(0,N1).transpose();
            m2t.colRange(N1,N).setZero();
            if (P) m2t *= *P;

#ifdef XDEBUG_QR
            if (Norm(minv*A1-T(1)) > 1.e-3*Norm(A1)) {
                std::cout<<"QR = "<<QR<<std::endl;
                std::cout<<"Q = "<<Q1<<std::endl;
                std::cout<<"R = "<<R1<<std::endl;
                if (P) std::cout<<"P = "<<*P<<std::endl;
                std::cout<<"A = QRP = "<<A1<<std::endl;
                std::cout<<"minv = "<<minv<<std::endl;
                std::cout<<"minv * A = "<<(minv*A1)<<std::endl;
                std::cout<<"A * minv = "<<(A1*minv)<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo 90: call InstQR_Inverse
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_Inverse_Helper<90,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& minv)
        { InstQR_Inverse(QR.xView(),beta.xView(),P,N1,minv.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_Inverse_Helper<97,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& minv)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c QRc = QR.conjugate();
            M2c minvc = minv.conjugate();
            QR_Inverse_Helper<-2,cs,rs,M1c,V1,M2c>::call(QRc,beta,P,N1,minvc);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_Inverse_Helper<-3,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& minv)
        {
            const int algo = 
                cs == 0 || rs == 0 ? 0 : 
                // algo 11 is a bit better for smaller matrices.
                // Not a big enough difference to warrant a runtime 
                // check if rs == TMV_UNKNOWN though.
                (rs != TMV_UNKNOWN && rs < 16) ? 11 :
                12;
#ifdef PRINTALGO_QR
            std::cout<<"Inline QRInverse\n";
            std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"QR = "<<TMV_Text(QR)<<std::endl;
            std::cout<<"beta = "<<TMV_Text(beta)<<std::endl;
            std::cout<<"minv = "<<TMV_Text(minv)<<std::endl;
            std::cout<<"P = "<<bool(P)<<std::endl;
            std::cout<<"N1 = "<<N1<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            QR_Inverse_Helper<algo,cs,rs,M1,V1,M2>::call(QR,beta,P,N1,minv);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_Inverse_Helper<-2,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& minv)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T2>::isinst;
            const int algo = 
                cs == 0 || rs == 0 ? 0 : 
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            QR_Inverse_Helper<algo,cs,rs,M1,V1,M2>::call(QR,beta,P,N1,minv);
        }
    };

    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_Inverse_Helper<-1,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& minv)
        { QR_Inverse_Helper<-2,cs,rs,M1,V1,M2>::call(QR,beta,P,N1,minv); }
    };

    template <class M1, class V1, class M2>
    inline void InlineQR_Inverse(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& minv)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVAssert(QR.rowsize() == beta.size());
        TMVAssert(QR.colsize() == minv.rowsize());
        TMVAssert(QR.rowsize() == minv.colsize());
        if (P) TMVAssert(QR.rowsize() == P->size());
        TMVAssert(N1 <= int(QR.rowsize()));

        const int cs = Sizes<M2::_rowsize,M1::_colsize>::size;
        const int rs = Sizes<M2::_colsize,M1::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) minvv = minv.cView();
        QR_Inverse_Helper<-3,cs,rs,M1v,V1v,M2v>::call(QRv,betav,P,N1,minvv);
    }

    template <class M1, class V1, class M2>
    inline void QR_Inverse(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& minv)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVAssert(QR.rowsize() == beta.size());
        TMVAssert(QR.colsize() == minv.rowsize());
        TMVAssert(QR.rowsize() == minv.colsize());
        if (P) TMVAssert(QR.rowsize() == P->size());
        TMVAssert(N1 <= int(QR.rowsize()));

        const int cs = Sizes<M2::_rowsize,M1::_colsize>::size;
        const int rs = Sizes<M2::_colsize,M1::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) minvv = minv.cView();
        QR_Inverse_Helper<-2,cs,rs,M1v,V1v,M2v>::call(QRv,betav,P,N1,minvv);
    }

    template <int algo, int cs, int rs, class M1, class V1, class M2>
    struct QR_InverseATA_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0)
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_InverseATA_Helper<0,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& , const V1& , const Permutation* , int, M2& ) {} 
    };

    // algo 11: Normal case
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_InverseATA_Helper<11,cs,rs,M1,V1,M2>
    {
        static void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& ata)
        {
#ifdef PRINTALGO_QR
            std::cout<<"QRInverseATA algo 11: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            // A = QRP
            // AtA = PtRtQtQRP
            //     = Pt Rt R P
            // (AtA)^-1 = Pt R^-1 Rt^-1 P
            ata.setZero();
            typename M2::colrange_type::uppertri_type rinv = 
                ata.colRange(0,N1).upperTri();
            rinv = QR.upperTri().subTriMatrix(0,N1);
            rinv.invertSelf();
            ata.subMatrix(0,N1,0,N1) = rinv * rinv.adjoint();
            if (P) {
                ata /= *P;
                ata *= *P;
            }
        }
    };

    // algo 90: call InstQR_InverseATA
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_InverseATA_Helper<90,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& ata)
        { InstQR_InverseATA(QR.xView(),beta.xView(),P,N1,ata.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_InverseATA_Helper<97,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& ata)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c QRc = QR.conjugate();
            M2c atac = ata.conjugate();
            QR_InverseATA_Helper<-2,cs,rs,M1c,V1,M2c>::call(QRc,beta,P,N1,atac);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_InverseATA_Helper<-3,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& ata)
        {
            const int algo = 
                cs == 0 || rs == 0 ? 0 : 
                11;
#ifdef PRINTALGO_QR
            std::cout<<"Inline QRInverseATA\n";
            std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            QR_InverseATA_Helper<algo,cs,rs,M1,V1,M2>::call(QR,beta,P,N1,ata);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_InverseATA_Helper<-2,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& ata)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T2>::isinst;
            const int algo = 
                cs == 0 || rs == 0 ? 0 : 
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            QR_InverseATA_Helper<algo,cs,rs,M1,V1,M2>::call(QR,beta,P,N1,ata);
        }
    };

    template <int cs, int rs, class M1, class V1, class M2>
    struct QR_InverseATA_Helper<-1,cs,rs,M1,V1,M2>
    {
        static TMV_INLINE void call(
            const M1& QR, const V1& beta, const Permutation* P, int N1, M2& ata)
        { QR_InverseATA_Helper<-2,cs,rs,M1,V1,M2>::call(QR,beta,P,N1,ata); }
    }; 

    template <class M1, class V1, class M2>
    inline void InlineQR_InverseATA(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& ata)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVAssert(QR.rowsize() == beta.size());
        TMVAssert(QR.rowsize() == ata.rowsize());
        TMVAssert(QR.rowsize() == ata.colsize());
        if (P) TMVAssert(QR.rowsize() == P->size());
        TMVAssert(N1 <= int(QR.rowsize()));

        const int cs = M1::_colsize;
        const int rs = Sizes<M2::_colsize,
              Sizes<M2::_rowsize,M1::_rowsize>::size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) atav = ata.cView();
        QR_InverseATA_Helper<-3,cs,rs,M1v,V1v,M2v>::call(QRv,betav,P,N1,atav);
    }

    template <class M1, class V1, class M2>
    inline void QR_InverseATA(
        const BaseMatrix_Rec<M1>& QR, const BaseVector_Calc<V1>& beta,
        const Permutation* P, int N1, BaseMatrix_Rec_Mutable<M2>& ata)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVAssert(QR.rowsize() == beta.size());
        TMVAssert(QR.rowsize() == ata.rowsize());
        TMVAssert(QR.rowsize() == ata.colsize());
        if (P) TMVAssert(QR.rowsize() == P->size());
        TMVAssert(N1 <= int(QR.rowsize()));

        const int cs = M1::_colsize;
        const int rs = Sizes<M2::_colsize,
              Sizes<M2::_rowsize,M1::_rowsize>::size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) QRv = QR.cView();
        TMV_MAYBE_CREF(V1,V1v) betav = beta.cView();
        TMV_MAYBE_REF(M2,M2v) atav = ata.cView();
        QR_InverseATA_Helper<-2,cs,rs,M1v,V1v,M2v>::call(QRv,betav,P,N1,atav);
    }

} // namespace tmv

#endif

