
#include "tmv/TMV_MultUD.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultXU.h"
#include "tmv/TMV_ScaleU.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_ConjugateV.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_MultXD.h"

namespace tmv {

    // 
    //
    // MultUD
    //


    template <class M1, class M2>
    static void DoMultEq(M1& m1, const M2& m2)
    {
        typedef typename M1::real_type RT;
        Scaling<1,RT> one;
        if (m2.step() == 1) {
            InlineMultMM<false>(one,m1.constView(),m2.unitView(),m1);
        } else {
            InlineMultMM<false>(one,m1.constView(),m2.copy().view(),m1);
        }
    }

    template <class M1, class M2, class M3, class T3>
    static inline void DoInstMultMM(
        const T3 x, const M1& m1, const M2& m2, M3& m3)
    {
        if (SameStorage(m1,m3)) {
            // Must be exact same storage.  Just check conj and unitdiag.
            Maybe<M1::_conj>::conjself(m3);
            if (m1.isunit()) m3.diag().setAllTo(T3(1));
            InstScale(x,m3);
        } else {
            InstMultXM(x,m1,m3);
        }
        if (m3.iscm()) {
            typename M3::cmview_type m3cm = m3.cmView();
            DoMultEq(m3cm,m2);
        } else if (m3.isrm()) {
            typename M3::rmview_type m3rm = m3.rmView();
            DoMultEq(m3rm,m2);
        } else {
            typedef typename TypeSelect<M3::_upper,
                    UpperTriMatrix<T3,NonUnitDiag|ColMajor>,
                    LowerTriMatrix<T3,NonUnitDiag|ColMajor> >::type M3c;
            M3c m3c = m3;
            typename M3c::cmview_type m3cm = m3c.view();
            DoMultEq(m3cm,m2);
            InstCopy(m3c.constView().xView(),m3.xView());
        }
    }

    template <class M1, class M2, class M3, class T3>
    static inline void DoInstAddMultMM(
        const T3 x, const M1& m1, const M2& m2, M3& m3)
    {
        typedef typename Traits2<
            typename M1::value_type,typename M2::value_type>::type T12;
        Matrix<T12,ColMajor|NoDivider> m1c(m3.size(),m3.size());
        typedef typename TypeSelect<M3::_upper,
                UpperTriMatrixView<T12,NonUnitDiag|ColMajor>,
                LowerTriMatrixView<T12,NonUnitDiag|ColMajor> >::type M1t;
        M1t m1ct = Maybe<M3::_upper>::uppertri(m1c);
        InstCopy(m1,m1ct.xView());
        DoMultEq(m1ct,m2);
        InstAddMultXM(x,m1ct.constView().xView(),m3);
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3)
    { DoInstMultMM(x,m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3)
    { DoInstAddMultMM(x,m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3)
    { DoInstMultMM(x,m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3)
    { DoInstAddMultMM(x,m1,m2,m3); }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3)
    { InlineAliasMultMM<false>(Scaling<0,T3>(x),m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3)
    { InlineAliasMultMM<true>(Scaling<0,T3>(x),m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3)
    { InlineAliasMultMM<false>(Scaling<0,T3>(x),m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3)
    { InlineAliasMultMM<true>(Scaling<0,T3>(x),m1,m2,m3); }

#define InstFile "TMV_MultUD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


