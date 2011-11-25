
//#define PRINTALGO_DivU
//#define PRINTALGO_UM
//#define XDEBUG_DivU

#include "TMV_Blas.h"
#include "tmv/TMV_DivMU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SmallTriMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_Det.h"
#include "tmv/TMV_ScaleM.h"

namespace tmv {

    template <class M1, class M2>
    static void DoTriLDivEq(M1& m1, const M2& m2)
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

#ifdef BLAS
    template <class T1, class M2, class T2> 
    static inline void BlasTriLDivEq(MatrixView<T1> A, const M2& B, T2)
    { DoTriLDivEq(A,B); }
#ifdef TMV_INST_DOUBLE
    template <class M2> 
    static void BlasTriLDivEq(MatrixView<double> A, const M2& B, double )
    {
        TMVAssert(B.size() == A.colsize());
        TMVAssert(A.colsize()>0);
        TMVAssert(A.rowsize()>0);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());

        int m=A.iscm()?A.colsize():A.rowsize();
        int n=A.iscm()?A.rowsize():A.colsize();
        double alpha(1);
        int lda = A.iscm()?A.stepj():A.stepi();
        int ldb = B.iscm()?B.stepj():B.stepi();

        BLASNAME(dtrsm) (
            BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
            B.iscm()==M2::_upper?BLASCH_UP:BLASCH_LO, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            B.isunit()?BLASCH_U:BLASCH_NU, 
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(B.cptr()),BLASV(ldb),
            BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <class M2>
    static void BlasTriLDivEq(
        MatrixView<std::complex<double> > A, const M2& B, std::complex<double> )
    {
        TMVAssert(B.size() == A.colsize());
        TMVAssert(A.colsize()>0);
        TMVAssert(A.rowsize()>0);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());

        int m=A.iscm()?A.colsize():A.rowsize();
        int n=A.iscm()?A.rowsize():A.colsize();
        std::complex<double> alpha(1);
        int lda = A.iscm()?A.stepj():A.stepi();
        int ldb = B.iscm()?B.stepj():B.stepi();
        if (A.iscm()==B.iscm() && B.isconj()) {
            A.conjugateSelf();
            BLASNAME(ztrsm) (
                BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
                B.iscm()==M2::_upper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                B.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(B.cptr()),BLASV(ldb),
                BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
            A.conjugateSelf();
        } else {
            BLASNAME(ztrsm) (
                BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
                B.iscm()==M2::_upper?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
                B.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(B.cptr()),BLASV(ldb),
                BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <class M2>
    static void BlasTriLDivEq(
        MatrixView<std::complex<double> > A, const M2& B, double )
    {
        TMVAssert(B.size() == A.colsize());
        TMVAssert(A.colsize()>0);
        TMVAssert(A.rowsize()>0);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());

        if (A.iscm()) {
            Matrix<double,ColMajor> A1 = A.realPart();
            BlasTriLDivEq(A1.xView(),B,double(0));
            A.realPart() = A1;
            A1 = A.imagPart();
            BlasTriLDivEq(A1.xView(),B,double(0));
            A.imagPart() = A1;
        } else {
            int m=2*A.rowsize();
            int n=A.colsize();
            double alpha(1);
            int lda = 2*A.stepi();
            int ldb = B.iscm()?B.stepj():B.stepi();

            BLASNAME(dtrsm) (
                BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
                B.iscm()==M2::_upper?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
                B.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASV(alpha),BLASP(B.cptr()),BLASV(ldb),
                BLASP((double*)A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
#endif
#ifdef TMV_INST_FLOAT
    template <class M2> 
    static void BlasTriLDivEq(MatrixView<float> A, const M2& B, float )
    {
        TMVAssert(B.size() == A.colsize());
        TMVAssert(A.colsize()>0);
        TMVAssert(A.rowsize()>0);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());

        int m=A.iscm()?A.colsize():A.rowsize();
        int n=A.iscm()?A.rowsize():A.colsize();
        float alpha(1);
        int lda = A.iscm()?A.stepj():A.stepi();
        int ldb = B.iscm()?B.stepj():B.stepi();

        BLASNAME(strsm) (
            BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
            B.iscm()==M2::_upper?BLASCH_UP:BLASCH_LO, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            B.isunit()?BLASCH_U:BLASCH_NU, 
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(B.cptr()),BLASV(ldb),
            BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <class M2>
    static void BlasTriLDivEq(
        MatrixView<std::complex<float> > A, const M2& B, std::complex<float> )
    {
        TMVAssert(B.size() == A.colsize());
        TMVAssert(A.colsize()>0);
        TMVAssert(A.rowsize()>0);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());

        int m=A.iscm()?A.colsize():A.rowsize();
        int n=A.iscm()?A.rowsize():A.colsize();
        std::complex<float> alpha(1);
        int lda = A.iscm()?A.stepj():A.stepi();
        int ldb = B.iscm()?B.stepj():B.stepi();
        if (A.iscm()==B.iscm() && B.isconj()) {
            A.conjugateSelf();
            BLASNAME(ctrsm) (
                BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
                B.iscm()==M2::_upper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                B.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(B.cptr()),BLASV(ldb),
                BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
            A.conjugateSelf();
        } else {
            BLASNAME(ctrsm) (
                BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
                B.iscm()==M2::_upper?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
                B.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(B.cptr()),BLASV(ldb),
                BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <class M2>
    static void BlasTriLDivEq(
        MatrixView<std::complex<float> > A, const M2& B, float )
    {
        TMVAssert(B.size() == A.colsize());
        TMVAssert(A.colsize()>0);
        TMVAssert(A.rowsize()>0);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());

        if (A.iscm()) {
            Matrix<float,ColMajor> A1 = A.realPart();
            BlasTriLDivEq(A1.xView(),B,float(0));
            A.realPart() = A1;
            A1 = A.imagPart();
            BlasTriLDivEq(A1.xView(),B,float(0));
            A.imagPart() = A1;
        } else {
            int m=2*A.rowsize();
            int n=A.colsize();
            float alpha(1);
            int lda = 2*A.stepi();
            int ldb = B.iscm()?B.stepj():B.stepi();

            BLASNAME(strsm) (
                BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
                B.iscm()==M2::_upper?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
                B.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASV(alpha),BLASP(B.cptr()),BLASV(ldb),
                BLASP((float*)A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
#endif
#endif // BLAS

    template <class T1, class M2>
    static void CallTriLDivEq(MatrixView<T1> m1, const M2& m2)
    {
        typedef typename M2::value_type T2;
        if (m1.colsize() > 0 && m1.rowsize() > 0) {
#ifdef BLAS
            const T2 t2(0);
            if ((m1.isrm() && m1.stepi()>0) || (m1.iscm() && m1.stepj()>0) ) {
                if ((m2.isrm() && m2.stepi()>0) || (m2.iscm() && m2.stepj()>0)) {
                    BlasTriLDivEq(m1,m2,t2);
                } else {
                    Matrix<T2,ColMajor|NoDivider> m2c(m2.size(),m2.size());
                    typedef typename TypeSelect<M2::_upper,
                            UpperTriMatrixView<T2,ColMajor>,
                            LowerTriMatrixView<T2,ColMajor> >::type M2t;
                    M2t m2ct = Maybe<M2::_upper>::uppertri(m2c,m2.dt());
                    InstCopy(m2,m2ct.xView());
                    DoTriLDivEq(m1,m2ct.constView());
                }
            } else {
                Matrix<T1,ColMajor|NoDivider> m1c(m1);
                CallTriLDivEq(m1c.xView(),m2);
                InstCopy(m1c.constView().xView(),m1);
            }
#else
            if (m1.iscm() || m1.isrm()) {
                if (m2.iscm() || m2.isrm()) {
                    DoTriLDivEq(m1,m2);
                } else {
                    Matrix<T2,ColMajor|NoDivider> m2c(m2.size(),m2.size());
                    typedef typename TypeSelect<M2::_upper,
                            UpperTriMatrixView<T2,ColMajor>,
                            LowerTriMatrixView<T2,ColMajor> >::type M2t;
                    M2t m2ct = Maybe<M2::_upper>::uppertri(m2c,m2.dt());
                    InstCopy(m2,m2ct.xView());
                    DoTriLDivEq(m1,m2ct.constView());
                }
            } else {
                Matrix<T1,ColMajor|NoDivider> m1c(m1);
                CallTriLDivEq(m1c.xView(),m2);
                InstCopy(m1c.constView().xView(),m1);
            }
#endif
        }
    }

    template <class T1, class T2, int C2>
    void InstTriLDivEq(
        MatrixView<T1> m1, const ConstUpperTriMatrixView<T2,C2>& m2)
    { CallTriLDivEq(m1,m2); }

    template <class T1, class T2, int C2>
    void InstTriLDivEq(
        MatrixView<T1> m1, const ConstLowerTriMatrixView<T2,C2>& m2)
    { CallTriLDivEq(m1,m2); }

    template <class T1, class T2, int C2>
    void InstAliasTriLDivEq(
        MatrixView<T1> m1, const ConstUpperTriMatrixView<T2,C2>& m2)
    { InlineAliasTriLDivEq(m1,m2); }

    template <class T1, class T2, int C2>
    void InstAliasTriLDivEq(
        MatrixView<T1> m1, const ConstLowerTriMatrixView<T2,C2>& m2)
    { InlineAliasTriLDivEq(m1,m2); }


#define InstFile "TMV_DivMU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


