
//#define PRINTALGO_DivU

#include "tmv/TMV_DivUU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallTriMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_Det.h"

namespace tmv {

    template <class M1, class M2>
    static void DoTriLDivEq2(M1& m1, const M2& m2)
    {
        TMVAssert(m1.isrm() || m1.iscm());
        TMVAssert(m2.isrm() || m2.iscm());
        if (m1.iscm()) {
            typename M1::cmview_type m1cm = m1.cmView();
            if (m2.iscm())
                InlineTriLDivEq(m1cm,m2.cmView());
            else
                InlineTriLDivEq(m1cm,m2.rmView());
        } else {
            typename M1::rmview_type m1rm = m1.rmView();
            if (m2.iscm())
                InlineTriLDivEq(m1rm,m2.cmView());
            else
                InlineTriLDivEq(m1rm,m2.rmView());
        }
    }

    template <class M1, class M2>
    static void DoTriLDivEq(M1& m1, const M2& m2)
    {
        if (m1.iscm() || m1.isrm()) {
            if (m2.iscm() || m2.isrm()) {
                DoTriLDivEq2(m1,m2);
            } else {
                if (m2.isunit()) 
                    DoTriLDivEq2(
                        m1,m2.copy().viewAsUnitDiag().constView().xView());
                else 
                    DoTriLDivEq2(m1,m2.copy().constView().xView());
            }
        } else {
            if (m1.isunit()) {
                typename M1::copy_type m1c(m1);
                typename M1::copy_type::unitdiag_type::xview_type m1cv = 
                    m1c.viewAsUnitDiag().xView();
                DoTriLDivEq(m1cv,m2);
                InstCopy(m1cv.constView(),m1);
            } else {
                typename M1::copy_type m1c(m1);
                typename M1::copy_type::xview_type m1cv = m1c.xView();
                DoTriLDivEq(m1cv,m2);
                InstCopy(m1cv.constView(),m1);
            }
        }
    }

    template <class T1, class T2, int C2>
    void InstTriLDivEq(
        UpperTriMatrixView<T1> m1, const ConstUpperTriMatrixView<T2,C2>& m2)
    { DoTriLDivEq(m1,m2); }
    template <class T1, class T2, int C2>
    void InstTriLDivEq(
        LowerTriMatrixView<T1> m1, const ConstLowerTriMatrixView<T2,C2>& m2)
    { DoTriLDivEq(m1,m2); }

    template <class T1, class T2, int C2>
    void InstAliasTriLDivEq(
        UpperTriMatrixView<T1> m1, const ConstUpperTriMatrixView<T2,C2>& m2)
    { InlineAliasTriLDivEq(m1,m2); }
    template <class T1, class T2, int C2>
    void InstAliasTriLDivEq(
        LowerTriMatrixView<T1> m1, const ConstLowerTriMatrixView<T2,C2>& m2)
    { InlineAliasTriLDivEq(m1,m2); }

#define InstFile "TMV_DivUU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


