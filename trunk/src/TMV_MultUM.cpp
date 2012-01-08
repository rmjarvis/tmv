
//#define PRINTALGO_UM

#include "TMV_Blas.h"
#include "tmv/TMV_MultUM.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_MultXU.h"
#include "tmv/TMV_ScaleM.h"

namespace tmv {

    template <class M1, class M2>
    static void DoMultEq(const M1& m1, M2& m2)
    {
        TMVAssert(m1.iscm() || m1.isrm());
        TMVAssert(m2.iscm() || m2.isrm());
        typedef typename M1::real_type RT;
        Scaling<1,RT> one;
        if (m2.iscm()) {
            typename M2::cmview_type m2cm = m2.cmView();
            if (m1.iscm())
                InlineMultMM<false>(one,m1.cmView(),m2cm,m2cm);
            else
                InlineMultMM<false>(one,m1.rmView(),m2cm,m2cm);
        } else {
            typename M2::rmview_type m2rm = m2.rmView();
            if (m1.iscm())
                InlineMultMM<false>(one,m1.cmView(),m2rm,m2rm);
            else
                InlineMultMM<false>(one,m1.rmView(),m2rm,m2rm);
        }
    }

#ifdef BLAS
    template <class M1, class M2, class T1>
    static void BlasMultEq(const M1& m1, M2& m2, T1)
    { DoMultEq(m1,m2); }
#ifdef TMV_INST_DOUBLE
    template <class M1>
    static void BlasMultEq(const M1& A, MatrixView<double> B, double)
    {
        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.iscm()?A.stepj():A.stepi();
        int ldb = B.iscm()?B.stepj():B.stepi();
        const double alpha(1);

        BLASNAME(dtrmm) (
            BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
            A.iscm()==M1::_upper?BLASCH_UP:BLASCH_LO, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
            BLASV(alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
            BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <class M1>
    static void BlasMultEq(
        const M1& A, MatrixView<std::complex<double> > B, std::complex<double>)
    {
        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.iscm()?A.stepj():A.stepi();
        int ldb = B.iscm()?B.stepj():B.stepi();
        const std::complex<double> alpha(1);

        if (A.iscm()==B.iscm() && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ztrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()==M1::_upper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ztrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()==M1::_upper?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <class M1>
    static void BlasMultEq(
        const M1& A, MatrixView<std::complex<double> > B, double t1)
    {
        if (B.iscm()) {
            Matrix<double,ColMajor> B1 = B.realPart();
            BlasMultEq(A,B1.xView(),t1);
            B.realPart() = B1;
            B1 = B.imagPart();
            BlasMultEq(A,B1.xView(),t1);
            B.imagPart() = B1;
        } else {
            int m = 2*B.rowsize();
            int n = B.colsize();
            int lda = A.iscm()?A.stepj():A.stepi();
            int ldb = 2*B.stepi();
            double alpha(1);
            BLASNAME(dtrmm) (
                BLASCM BLASCH_R, 
                A.iscm()==M1::_upper?BLASCH_UP:BLASCH_LO, 
                A.iscm()?BLASCH_T:BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASV(alpha),BLASP(A.cptr()),BLASV(lda), 
                BLASP((double*)B.ptr()), BLASV(ldb)
                BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
#endif
#ifdef TMV_INST_FLOAT
    template <class M1>
    static void BlasMultEq(const M1& A, MatrixView<float> B, float)
    {
        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.iscm()?A.stepj():A.stepi();
        int ldb = B.iscm()?B.stepj():B.stepi();
        const float alpha(1);

        BLASNAME(strmm) (
            BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
            A.iscm()==M1::_upper?BLASCH_UP:BLASCH_LO, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
            BLASV(alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
            BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <class M1>
    static void BlasMultEq(
        const M1& A, MatrixView<std::complex<float> > B, std::complex<float>)
    {
        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.iscm()?A.stepj():A.stepi();
        int ldb = B.iscm()?B.stepj():B.stepi();
        const std::complex<float> alpha(1);

        if (A.iscm()==B.iscm() && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ctrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()==M1::_upper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ctrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()==M1::_upper?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <class M1>
    static void BlasMultEq(
        const M1& A, MatrixView<std::complex<float> > B, float t1)
    {
        if (B.iscm()) {
            Matrix<float,ColMajor> B1 = B.realPart();
            BlasMultEq(A,B1.xView(),t1);
            B.realPart() = B1;
            B1 = B.imagPart();
            BlasMultEq(A,B1.xView(),t1);
            B.imagPart() = B1;
        } else {
            int m = 2*B.rowsize();
            int n = B.colsize();
            int lda = A.iscm()?A.stepj():A.stepi();
            int ldb = 2*B.stepi();
            float alpha(1);
            BLASNAME(strmm) (
                BLASCM BLASCH_R, 
                A.iscm()==M1::_upper?BLASCH_UP:BLASCH_LO, 
                A.iscm()?BLASCH_T:BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASV(alpha),BLASP(A.cptr()),BLASV(lda), 
                BLASP((float*)B.ptr()), BLASV(ldb)
                BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
#endif
#endif // BLAS

    template <class M1, class T2>
    static void CallMultEq(const M1& m1, MatrixView<T2> m2)
    {
        typedef typename M1::value_type T1;
#ifdef BLAS
        if (m2.colsize() > 0 && m2.rowsize() > 0) {
            const T1 t1(0);
            if ((m2.iscm() && m2.stepj()>0) || (m2.isrm() && m2.stepi()>0)) {
                if ((m1.iscm() && m1.stepj()>0) || (m1.isrm() && m1.stepi()>0)) {
                    BlasMultEq(m1,m2,t1);
                } else {
                    Matrix<T1,ColMajor|NoDivider> m1c(m1.size(),m1.size());
                    typedef typename TypeSelect<M1::_upper,
                            UpperTriMatrixView<T1,ColMajor>,
                            LowerTriMatrixView<T1,ColMajor> >::type M1t;
                    M1t m1ct = Maybe<M1::_upper>::uppertri(m1c,m1.dt());
                    InstCopy(m1,m1ct.xView());
                    BlasMultEq(m1ct.constView(),m2,t1);
                }
            } else {
                Matrix<T2,ColMajor|NoDivider> m2c(m2);
                CallMultEq(m1,m2c.xView());
                InstCopy(m2c.constView().xView(),m2);
            }
        }
#else
        if (m2.iscm() || m2.isrm()) {
            if (m1.iscm() || m1.isrm()) {
                DoMultEq(m1,m2);
            } else {
                Matrix<T1,ColMajor|NoDivider> m1c(m1.size(),m1.size());
                typedef typename TypeSelect<M1::_upper,
                        UpperTriMatrixView<T1,ColMajor>,
                        LowerTriMatrixView<T1,ColMajor> >::type M1t;
                M1t m1ct = Maybe<M1::_upper>::uppertri(m1c,m1.dt());
                InstCopy(m1,m1ct.xView());
                DoMultEq(m1ct.constView(),m2);
            }
        } else {
            Matrix<T2,ColMajor|NoDivider> m2c(m2);
            CallMultEq(m1,m2c.xView());
            InstCopy(m2c.constView().xView(),m2);
        }
#endif
    }

    template <class T, class M1, class M2>
    static void DoInstMultMM(
        const T x, const M1& m1, const M2& m2, MatrixView<T> m3)
    {
        InstMultXM(x,m2,m3);
        CallMultEq(m1,m3);
    }

    template <class T, class M1, class M2>
    static void DoInstAddMultMM(
        const T x, const M1& m1, const M2& m2, MatrixView<T> m3)
    {
        Matrix<T,ColMajor|NoDivider> m3c(m3.colsize(),m3.rowsize());
        InstMultXM(x,m2,m3c.xView());
        CallMultEq(m1,m3c.xView());
        InstAddMultXM(T(1),m3c.xView().constView(),m3);
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { DoInstMultMM(x,m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { DoInstAddMultMM(x,m1,m2,m3); }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { DoInstMultMM(x,m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { DoInstAddMultMM(x,m1,m2,m3); }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<false>(Scaling<0,T3>(x),m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<true>(Scaling<0,T3>(x),m1,m2,m3); }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<false>(Scaling<0,T3>(x),m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<true>(Scaling<0,T3>(x),m1,m2,m3); }

#define InstFile "TMV_MultUM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


