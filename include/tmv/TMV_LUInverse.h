
#ifndef TMV_LUInverse_H
#define TMV_LUInverse_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Permutation.h"
#include "TMV_MultMM_Funcs.h"

#ifdef PRINTALGO_LU
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_TriMatrixIO.h"
#endif

namespace tmv {

    // Defined in TMV_LUInverse.cpp
    template <class T1>
    void InstLU_Inverse(MatrixView<T1> m1, const Permutation& P);

    template <class T1, class T2, int C1>
    void InstLU_InverseATA(
        const ConstMatrixView<T1,C1>& m1, const Permutation& P, 
        const bool trans, MatrixView<T2> m2);


    template <int algo, int cs, int rs, class M2>
    struct LU_Inverse_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0)
    template <int cs, int rs, class M1>
    struct LU_Inverse_Helper<0,cs,rs,M1>
    { static TMV_INLINE void call(M1& , const Permutation& ) {} };

    // algo 11: Normal case
    template <int cs, int rs, class M1>
    struct LU_Inverse_Helper<11,cs,rs,M1>
    {
        static void call(M1& m1, const Permutation& P)
        {
#ifdef PRINTALGO_LU
            std::cout<<"LUInverse algo 11: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            // m1 = (PLU)^-1
            //    = U^-1 L^-1 Pt
#ifdef PRINTALGO_LU
            std::cout<<"Start LU_Inverse\n";
            typedef typename M1::value_type T;
#endif
            typename M1::uppertri_type U = m1.upperTri();
            typename M1::unit_lowertri_type L = m1.unitLowerTri();
#ifdef XDEBUG_LU
            Matrix<T> L0 = L;
            Matrix<T> U0 = U;
            Matrix<T> A = L*U;
            P.applyOnLeft(A);
#endif
            U.invertSelf();
            L.invertSelf();
#ifdef XDEBUG_LU
            Matrix<T> Linv = L;
            Matrix<T> Uinv = U;
#endif
            const Scaling<1,typename M1::real_type> one;
            typename M1::noalias_type m1na = m1.noAlias();
            MultMM<false>(one,U,L,m1na);
            P.inverse().applyOnRight(m1);
#ifdef XDEBUG_LU
            if (Norm(m1*A-T(1)) > 1.e-3*Norm(A)) {
                std::cout<<"L = "<<L0<<std::endl;
                std::cout<<"U = "<<U0<<std::endl;
                std::cout<<"A = PLU = "<<A<<std::endl;
                std::cout<<"U^-1 = "<<Uinv<<std::endl;
                std::cout<<"L^-1 = "<<Linv<<std::endl;
                std::cout<<"U^-1 * U = "<<(Uinv*U0)<<std::endl;
                std::cout<<"U * U^-1 = "<<(U0*Uinv)<<std::endl;
                std::cout<<"L^-1 * L = "<<(Linv*L0)<<std::endl;
                std::cout<<"L * L^-1 = "<<(L0*Linv)<<std::endl;
                std::cout<<"Direct U^-1 L^-1 = "<<(Uinv*Linv)<<std::endl;
                std::cout<<"U^-1 L^-1 P^-1 = "<<m1<<std::endl;
                std::cout<<"m1 * A = "<<(m1*A)<<std::endl;
                std::cout<<"A * m1 = "<<(A*m1)<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo 90: call InstLU_Inverse
    template <int cs, int rs, class M1>
    struct LU_Inverse_Helper<90,cs,rs,M1>
    {
        static TMV_INLINE void call(M1& m1, const Permutation& P)
        { InstLU_Inverse(m1.xView(),P); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1>
    struct LU_Inverse_Helper<97,cs,rs,M1>
    {
        static TMV_INLINE void call(M1& m1, const Permutation& P)
        {
            typedef typename M1::conjugate_type M1c;
            M1c m1c = m1.conjugate();
            LU_Inverse_Helper<-2,cs,rs,M1c>::call(m1c,P);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1>
    struct LU_Inverse_Helper<-3,cs,rs,M1>
    {
        static TMV_INLINE void call(M1& m1, const Permutation& P)
        {
            const int algo = 
                cs == 0 || rs == 0 ? 0 : 
                11;
#ifdef PRINTALGO_LU
            std::cout<<"Inline LUInverse\n";
            std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
#ifdef XDEBUG_LU
            typedef typename M1::value_type T;
            typename M1::uppertri_type U = m1.upperTri();
            typename M1::unit_lowertri_type L = m1.unitLowerTri();
            Matrix<T> L0 = L;
            Matrix<T> U0 = U;
            Matrix<T> A = L*U;
            P.applyOnLeft(A);
#endif
            LU_Inverse_Helper<algo,cs,rs,M1>::call(m1,P);
#ifdef XDEBUG_LU
            if (Norm(m1*A-T(1)) > 1.e-3*Norm(A)) {
                std::cout<<"L = "<<L0<<std::endl;
                std::cout<<"U = "<<U0<<std::endl;
                std::cout<<"A = PLU = "<<A<<std::endl;
                std::cout<<"U^-1 L^-1 P^-1 = "<<m1<<std::endl;
                std::cout<<"m1 * A = "<<(m1*A)<<std::endl;
                std::cout<<"A * m1 = "<<(A*m1)<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M1>
    struct LU_Inverse_Helper<-2,cs,rs,M1>
    {
        static TMV_INLINE void call(M1& m1, const Permutation& P)
        {
            typedef typename M1::value_type T1;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
                Traits<T1>::isinst;
            const int algo = 
                cs == 0 || rs == 0 ? 0 : 
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
            LU_Inverse_Helper<algo,cs,rs,M1>::call(m1,P);
        }
    };

    template <int cs, int rs, class M1>
    struct LU_Inverse_Helper<-1,cs,rs,M1>
    {
        static TMV_INLINE void call(M1& m1, const Permutation& P)
        { LU_Inverse_Helper<-2,cs,rs,M1>::call(m1,P); }
    };

    template <class M1>
    inline void InlineLU_Inverse(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& P)
    {
        const int cs = M1::_colsize;
        const int rs = M1::_rowsize;
        typedef typename M1::cview_type M1v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        LU_Inverse_Helper<-3,cs,rs,M1v>::call(m1v,P);
    }

    template <class M1>
    inline void LU_Inverse(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& P)
    {
        const int cs = M1::_colsize;
        const int rs = M1::_rowsize;
        typedef typename M1::cview_type M1v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        LU_Inverse_Helper<-2,cs,rs,M1v>::call(m1v,P);
    }

    template <int algo, int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0)
    // Also used for invalid real/complex combination from the virtual calls.
    template <int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper<0,cs,rs,M1,M2>
    { 
        static TMV_INLINE void call(
            const M1& , const Permutation& , const bool, M2& ) {} 
    };

    // algo 11: Normal case
    template <int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper<11,cs,rs,M1,M2>
    {
        static void call(
            const M1& m1, const Permutation& P, const bool trans, M2& m2)
        {
#ifdef PRINTALGO_LU
            std::cout<<"LUInverseATA algo 11: cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"trans = "<<trans<<std::endl;
            std::cout<<"m1 (=LU) = "<<m1<<std::endl;
            std::cout<<"m2 = "<<m2<<std::endl;
#endif
            // (At A)^-1 = A^-1 (A^-1)t
            // = (U^-1 L^-1 Pt) (P L^-1t U^-1t)
            // = U^-1 L^-1 L^-1t U^-1t
            //
            // if PLU is really AT, then
            // A^-1 = P L^-1T U^-1T
            // (At A)^-1 = P L^-1T U^-1T U^-1* L^-1* Pt

            typename M1::const_unit_lowertri_type L = m1.unitLowerTri();
            typename M1::const_uppertri_type U = m1.upperTri();
            Scaling<1,typename M2::real_type> one;
            typedef typename M2::noalias_type M2na;
            M2na m2na = m2.noAlias();

            if (trans) {
                typename M2na::uppertri_type uinv = m2na.upperTri();
                Copy(U,uinv);
                uinv.invertSelf();
                // m2 = uinv.transpose() * uinv.conjugate();
                // -> m2.transpose() = uinv.adjoint() * uinv;
                // Doing it this way means that the NoAlias call works,
                // even though there actually are aliases.
                typename M2na::transpose_type m2T = m2na.transpose();
                MultMM<false>(one,uinv.adjoint(),uinv,m2T);
                // m2 /= L.transpose();
                TriLDivEq(m2na,L.transpose());
                // m2 %= L.conjugate();
                // m2 = m2 * L.conjugate()^-1
                // -> m2.adjoint() = L.transpose()^-1 * m2.adjoint()
                typename M2na::adjoint_type m2t = m2na.adjoint();
                TriLDivEq(m2t,L.transpose());
                P.inverse().applyOnRight(m2);
                P.applyOnLeft(m2);
            } else {
                typename M2na::unit_lowertri_type linv = m2na.unitLowerTri();
                Copy(L,linv);
                linv.invertSelf();
                // m2 = linv * linv.adjoint();
                MultMM<false>(one,linv,linv.adjoint(),m2na);
                // m2 /= U;
                TriLDivEq(m2na,U);
                // m2 %= U.adjoint();
                typename M2na::adjoint_type m2t = m2na.adjoint();
                TriLDivEq(m2t,U);
            }
        }
    };

    // algo 90: call InstLU_InverseATA
    template <int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper<90,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(
            const M1& m1, const Permutation& P, const bool trans, M2& m2)
        {
            TMVAssert(m1.iscm());
            InstLU_InverseATA(m1.xView(),P,trans,m2.xView()); 
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper<97,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(
            const M1& m1, const Permutation& P, const bool trans, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            LU_InverseATA_Helper<-2,cs,rs,M1c,M2c>::call(m1c,P,trans,m2c);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper<-3,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(
            const M1& m1, const Permutation& P, const bool trans, M2& m2)
        {
            const bool invalid =
                M1::iscomplex && M2::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                11;
#ifdef PRINTALGO_LU
            std::cout<<"Inline LUInverseATA\n";
            std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
#ifdef XDEBUG_LU
            typedef typename M1::value_type T;
            typename M1::const_uppertri_type U = m1.upperTri();
            typename M1::const_unit_lowertri_type L = m1.unitLowerTri();
            Matrix<T> L0 = L;
            Matrix<T> U0 = U;
            Matrix<T> PLU = L*U;
            P.applyOnLeft(PLU);
            Matrix<T> A = PLU;
            if (trans) A.transposeSelf();
            Matrix<T> Ainv = m1;
            LU_Inverse(Ainv,P);
            if (trans) Ainv.transposeSelf();
            Matrix<T> invAtA = Ainv * Ainv.adjoint();
#endif
            LU_InverseATA_Helper<algo,cs,rs,M1,M2>::call(m1,P,trans,m2);
#ifdef XDEBUG_LU
            if (Norm(m2-invAtA) > 1.e-3*Norm(invAtA)) {
                std::cout<<"L = "<<L0<<std::endl;
                std::cout<<"U = "<<U0<<std::endl;
                std::cout<<"PLU = "<<PLU<<std::endl;
                std::cout<<"A = "<<A<<std::endl;
                std::cout<<"Ainv = "<<Ainv<<std::endl;
                std::cout<<"A*Ainv = "<<A*Ainv<<std::endl;
                std::cout<<"Ainv*A = "<<Ainv*A<<std::endl;
                std::cout<<"Ainv * Ainvt = "<<invAtA<<std::endl;
                std::cout<<"At * A = "<<A.adjoint()*A<<std::endl;
                std::cout<<"m2 = "<<m2<<std::endl;
                std::cout<<"(At * A) * m2 = "<<A.adjoint()*A*m2<<std::endl;
                std::cout<<"m2 * (At * A) = "<<m2*A.adjoint()*A<<std::endl;
                std::cout<<"(At * A) * Ainv*Ainvt = "<<A.adjoint()*A*invAtA<<std::endl;
                std::cout<<"Ainv*Ainvt * (At * A) = "<<invAtA*A.adjoint()*A<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper<-2,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(
            const M1& m1, const Permutation& P, const bool trans, M2& m2)
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
                Traits<T2>::isinst;
            const bool invalid =
                M1::iscomplex && M2::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            LU_InverseATA_Helper<algo,cs,rs,M1,M2>::call(m1,P,trans,m2);
        }
    };

    // algo -1: Check for aliases?  Not now.  Do we need this?
    template <int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper<-1,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(
            const M1& m1, const Permutation& P, const bool trans, M2& m2)
        { LU_InverseATA_Helper<-2,cs,rs,M1,M2>::call(m1,P,trans,m2); }
    };

    template <class M1, class M2>
    inline void InlineLU_InverseATA(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        const bool trans, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const int cs = M2::_colsize;
        const int rs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        LU_InverseATA_Helper<-3,cs,rs,M1v,M2v>::call(m1v,P,trans,m2v);
    }

    template <class M1, class M2>
    inline void LU_InverseATA(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        const bool trans, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const int cs = M2::_colsize;
        const int rs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        LU_InverseATA_Helper<-2,cs,rs,M1v,M2v>::call(m1v,P,trans,m2v);
    }

} // namespace tmv

#endif

