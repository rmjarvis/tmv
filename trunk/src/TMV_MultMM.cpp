
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_ProdMM.h"

namespace tmv {

    // Defined in TMV_MultMM_CCC.cpp
    // Defined in TMV_MultMM_CRC.cpp
    // Defined in TMV_MultMM_RCC.cpp
    // Defined in TMV_MultMM_RRC.cpp
    template <bool add, class T1, int C1, class T2, int C2, class T3>
    void DoInstMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3,ColMajor> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
#if TMV_OPT <= 2
        m3.setZero();
        InstAddMultMM(x,m1,m2,m3);
#else
        if (m3.iscm()) {
            MatrixView<T3,ColMajor> m3cm = m3.cmView();
            if (m1.iscm()) {
                if (m2.iscm())
                    DoInstMultMM<false>(x,m1.cmView(),m2.cmView(),m3cm);
                else if (m2.isrm())
                    DoInstMultMM<false>(x,m1.cmView(),m2.rmView(),m3cm);
                else {
                    Matrix<T2,ColMajor|NoDivider> m2c = m2;
                    DoInstMultMM<false>(x,m1.cmView(),m2c.constView(),m3cm);
                }
            } else if (m1.isrm()) {
                if (m2.iscm())
                    DoInstMultMM<false>(x,m1.rmView(),m2.cmView(),m3cm);
                else if (m2.isrm())
                    DoInstMultMM<false>(x,m1.rmView(),m2.rmView(),m3cm);
                else {
                    Matrix<T2,ColMajor|NoDivider> m2c = m2;
                    DoInstMultMM<false>(x,m1.rmView(),m2c.constView(),m3cm);
                }
            } else {
                Matrix<T1,RowMajor|NoDivider> m1c = m1;
                InstMultMM(x,m1c.constView().xView(),m2,m3);
            }
        } else if (m3.isrm()) {
            InstMultMM(x,m2.transpose(),m1.transpose(),m3.transpose());
        } else  {
            Matrix<T3,ColMajor|NoDivider> m3c(m3.colsize(),m3.rowsize());
            InstMultMM(x,m1,m2,m3c.xView());
            InstCopy(m3c.constView().xView(),m3);
        }
#endif
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
        if (m3.iscm()) {
            MatrixView<T3,ColMajor> m3cm = m3.cmView();
            if (m1.iscm()) {
                if (m2.iscm())
                    DoInstMultMM<true>(x,m1.cmView(),m2.cmView(),m3cm);
                else if (m2.isrm())
                    DoInstMultMM<true>(x,m1.cmView(),m2.rmView(),m3cm);
                else {
                    Matrix<T2,ColMajor|NoDivider> m2c = m2;
                    DoInstMultMM<true>(x,m1.cmView(),m2c.constView(),m3cm);
                }
            } else if (m1.isrm()) {
                if (m2.iscm())
                    DoInstMultMM<true>(x,m1.rmView(),m2.cmView(),m3cm);
                else if (m2.isrm())
                    DoInstMultMM<true>(x,m1.rmView(),m2.rmView(),m3cm);
                else {
                    Matrix<T2,ColMajor|NoDivider> m2c = m2;
                    DoInstMultMM<true>(x,m1.rmView(),m2c.constView(),m3cm);
                }
            } else {
                Matrix<T1,RowMajor|NoDivider> m1c = m1;
                InstAddMultMM(x,m1c.constView().xView(),m2,m3);
            }
        } else if (m3.isrm()) {
            InstAddMultMM(x,m2.transpose(),m1.transpose(),m3.transpose());
        } else  {
            Matrix<T3,ColMajor|NoDivider> m3c(m3.colsize(),m3.rowsize());
            InstMultMM(T3(1),m1,m2,m3c.xView());
            InstAddMultXM(x,m3c.constView().xView(),m3);
        }
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<false>(Scaling<0,T3>(x),m1,m2,m3); }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<true>(Scaling<0,T3>(x),m1,m2,m3); }

#define InstFile "TMV_MultMM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


