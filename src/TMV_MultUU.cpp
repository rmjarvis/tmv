
//#define PRINTALGO_UU

#include "tmv/TMV_MultUU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultXU.h"
#include "tmv/TMV_ScaleU.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_ConjugateV.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_ProdMM.h"

namespace tmv {

    template <class M2, class M3>
    static inline void DoMultEqUU(const M2& m2, M3& m3)
    {
        typedef typename M2::value_type T2;
        typedef typename M3::real_type RT;
        Scaling<1,RT> one;
        if (m2.iscm()) 
            InlineMultMM<false>(one,m3.constView(),m2.cmView(),m3);
        else if (m2.isrm())
            InlineMultMM<false>(one,m3.constView(),m2.rmView(),m3);
        else {
            DoMultEqUU(m2.copy().xView(),m3);
        }
    }

    template <class M1, class M2, class M3>
    static inline void DoMultUU(const M1& m1, const M2& m2, M3& m3)
    {
        typedef typename M2::value_type T2;
        typedef typename M3::real_type RT;

        // Can save some compiled code by copying m1 to m3.
        // This means the InlineMultMM call will always have m1 and m3
        // with the same majority.  And m1 will never be conjugated.
        //
        // Given the alias checking that has been done already,
        // there are some cases that were allowed that won't work
        // with the copy.
        //
        // Specifically, if diagonal of m2 has the same storage as m3,
        // then we need to make sure the copy won't change it.
        // (Remember that the only way m2 can have the same storage as m3
        // is if it is in the opposite triangle, so the diagonal is all 
        // we have to worry about.)
        if ( SameStorage(m2,m3) && !m2.isunit() &&
             (M1::_conj || m1.isunit() || !SameStorage(m1,m3)) ) {
            // Then the copy will clobber m2.diag().  
            // Need to use a temporary for m2.
            DoMultUU(m1,m2.copy().xView(),m3);
        } else {
            // OK to copy m1 to m3.  However, the copy might be trivial.
            // Check if m1 is same storage as m3.
            if ( SameStorage(m1,m3) ) {
                if (ExactSameStorage(m1,m3)) {
                    // Just check diagonal and conj.
                    Maybe<M1::_conj>::conjself(m3);
                    if (m1.isunit() && !m3.isunit()) m3.diag().setAllTo(RT(1));
                } else { // Opposite storage
                    if (m1.size() > 1)
                        InstCopy(m1.offDiag().xView(),m3.offDiag().xView());
                    if (!m3.isunit()) {
                        if (m1.isunit()) m3.diag().setAllTo(RT(1));
                        else Maybe<M1::_conj>::conjself2(m3.diag());
                    }
                }
            } else {
                InstCopy(m1,m3.xView());
            }
            DoMultEqUU(m2,m3);
        }
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstUpperTriMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3)
    {
        if (m3.iscm()) {
            UpperTriMatrixView<T3,ColMajor> m3cm = m3.cmView();
            DoMultUU(m1,m2,m3cm);
            if (x != T3(1)) InstScale(x,m3);
        } else if (m3.isrm()) {
            UpperTriMatrixView<T3,RowMajor> m3rm = m3.rmView();
            DoMultUU(m1,m2,m3rm);
            if (x != T3(1)) InstScale(x,m3);
        } else {
            Matrix<T3,ColMajor|NoDivider> m3c(m3.size(),m3.size(),T3(0));
            UpperTriMatrixView<T3,ColMajor> m3ct = m3c.upperTri(m3.dt());
            DoMultUU(m1,m2,m3ct);
            if (m3.isunit())
                InstCopy(m3ct.constView().xView(),m3);
            else
                InstMultXM(x,m3ct.constView().xView(),m3);
        }
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstUpperTriMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3)
    {
        typedef typename Traits2<T1,T2>::type T12;
        Matrix<T12,ColMajor|NoDivider> m1c(m3.size(),m3.size());
        UpperTriMatrixView<T12,ColMajor> m1ct = m1c.upperTri(
            (m1.isunit() && m2.isunit()) ? UnitDiag : NonUnitDiag);
        InstCopy(m1,m1ct.xView());
        DoMultUU(m1,m2,m1ct);
        InstAddMultXM(x,m1ct.constView().xView(),m3);
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstLowerTriMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3)
    {
        if (m3.iscm()) {
            LowerTriMatrixView<T3,ColMajor> m3cm = m3.cmView();
            DoMultUU(m1,m2,m3cm);
            if (x != T3(1)) InstScale(x,m3);
        } else if (m3.isrm()) {
            LowerTriMatrixView<T3,RowMajor> m3rm = m3.rmView();
            DoMultUU(m1,m2,m3rm);
            if (x != T3(1)) InstScale(x,m3);
        } else {
            Matrix<T3,ColMajor|NoDivider> m3c(m3.size(),m3.size(),T3(0));
            LowerTriMatrixView<T3,ColMajor> m3ct = m3c.lowerTri(m3.dt());
            DoMultUU(m1,m2,m3ct);
            if (m3.isunit())
                InstCopy(m3ct.constView().xView(),m3);
            else
                InstMultXM(x,m3ct.constView().xView(),m3);
        }
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstLowerTriMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3)
    {
        typedef typename Traits2<T1,T2>::type T12;
        Matrix<T12,ColMajor|NoDivider> m1c(m3.size(),m3.size());
        LowerTriMatrixView<T12,ColMajor> m1ct = m1c.lowerTri(
            (m1.isunit() && m2.isunit()) ? UnitDiag : NonUnitDiag);
        InstCopy(m1,m1ct.xView());
        DoMultUU(m1,m2,m1ct);
        InstAddMultXM(x,m1ct.constView().xView(),m3);
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstUpperTriMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3)
    { InlineAliasMultMM<false>(Scaling<0,T3>(x),m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstUpperTriMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3)
    { InlineAliasMultMM<true>(Scaling<0,T3>(x),m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstLowerTriMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3)
    { InlineAliasMultMM<false>(Scaling<0,T3>(x),m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstLowerTriMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3)
    { InlineAliasMultMM<true>(Scaling<0,T3>(x),m1,m2,m3); }

#define InstFile "TMV_MultUU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


