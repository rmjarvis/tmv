
//#undef NDEBUG
//#include "TMV.h"

#include "TMV_Blas.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_ProdMM.h"

#ifdef BLAS
#include "tmv/TMV_AddMM.h"
#include "tmv/TMV_SumMM.h"
#endif

namespace tmv {

    // Defined in TMV_MultMM_CCC.cpp
    // Defined in TMV_MultMM_CRC.cpp
    // Defined in TMV_MultMM_RCC.cpp
    // Defined in TMV_MultMM_RRC.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void DoInstMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3,ColMajor> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void DoInstAddMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3,ColMajor> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    static void DoMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
        TMVAssert(m1.iscm() || m1.isrm());
        TMVAssert(m2.iscm() || m2.isrm());
        TMVAssert(m3.iscm());
        if (m1.iscm())
            if (m2.iscm())
                DoInstMultMM(x,m1.cmView(),m2.cmView(),m3.cmView());
            else
                DoInstMultMM(x,m1.cmView(),m2.rmView(),m3.cmView());
        else 
            if (m2.iscm())
                DoInstMultMM(x,m1.rmView(),m2.cmView(),m3.cmView());
            else
                DoInstMultMM(x,m1.rmView(),m2.rmView(),m3.cmView());
    }

    template <class T1, int C1, class T2, int C2, class T3>
    static void DoAddMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
        TMVAssert(m1.iscm() || m1.isrm());
        TMVAssert(m2.iscm() || m2.isrm());
        TMVAssert(m3.iscm());
        if (m1.iscm())
            if (m2.iscm())
                DoInstAddMultMM(x,m1.cmView(),m2.cmView(),m3.cmView());
            else
                DoInstAddMultMM(x,m1.cmView(),m2.rmView(),m3.cmView());
        else 
            if (m2.iscm())
                DoInstAddMultMM(x,m1.rmView(),m2.cmView(),m3.cmView());
            else
                DoInstAddMultMM(x,m1.rmView(),m2.rmView(),m3.cmView());
    }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
    template <int C1, int C2>
    static void DoAddMultMM(
        double alpha, const ConstMatrixView<double,C1>& A,
        const ConstMatrixView<double,C2>& B, MatrixView<double> C)
    {
        int m = C.colsize();
        int n = C.rowsize();
        int k = A.rowsize();
        int lda = BlasIsCM(A)?A.stepj():A.stepi();
        int ldb = BlasIsCM(B)?B.stepj():B.stepi();
        int ldc = C.stepj();
        double xbeta(1.);
        BLASNAME(dgemm) (
            BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
            BlasIsCM(B)?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(k),BLASV(alpha),
            BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
            BLASV(xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
    template <int C1, int C2>
    static void DoAddMultMM(
        std::complex<double> alpha,
        const ConstMatrixView<std::complex<double>,C1>& A,
        const ConstMatrixView<std::complex<double>,C2>& B,
        MatrixView<std::complex<double> > C)
    {
        if (BlasIsCM(A) && A.isconj()) {
            const Matrix<std::complex<double>,ColMajor> AA = alpha*A;
            DoAddMultMM(std::complex<double>(1),AA.xView(),B,C);
        } else if (BlasIsCM(B) && B.isconj()) {
            const Matrix<std::complex<double>,ColMajor> BB = alpha*B;
            DoAddMultMM(std::complex<double>(1),A,BB.xView(),C);
        } else {
            int m = C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = BlasIsCM(A)?A.stepj():A.stepi();
            int ldb = BlasIsCM(B)?B.stepj():B.stepi();
            int ldc = C.stepj();
            std::complex<double> xbeta(1.);
            BLASNAME(zgemm) (
                BLASCM BlasIsCM(A)?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                BlasIsCM(B)?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(k),BLASP(&alpha),
                BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
                BLASP(&xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
        }
    }
#ifdef TMV_INST_MIX
    template <int C1, int C2>
    static void DoAddMultMM(
        std::complex<double> alpha,
        const ConstMatrixView<std::complex<double>,C1>& A,
        const ConstMatrixView<double,C2>& B,
        MatrixView<std::complex<double> > C)
    {
        if (BlasIsCM(A) && !A.isconj() && TMV_IMAG(alpha)==0.) {
            int m = 2*C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = 2*A.stepj();
            int ldb = BlasIsCM(B)?B.stepj():B.stepi();
            int ldc = 2*C.stepj();
            double xalpha(TMV_REAL(alpha));
            double xbeta(1.);
            BLASNAME(dgemm) (
                BLASCM BLASCH_NT, BlasIsCM(B)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
                BLASP((const double*)(A.cptr())),BLASV(lda),
                BLASP(B.cptr()),BLASV(ldb),
                BLASV(xbeta),BLASP((double*)(C.ptr())),BLASV(ldc) BLAS1 BLAS1);
        } else {
            if (TMV_IMAG(alpha) == 0.) {
                Matrix<double,ColMajor> Ax = A.realPart();
                Matrix<double,ColMajor> Cx = TMV_REAL(alpha)*Ax*B;
                C.realPart() += Cx;
                Ax = A.imagPart();
                if (A.isconj()) Cx = -TMV_REAL(alpha)*Ax*B;
                else Cx = TMV_REAL(alpha)*Ax*B;
                C.imagPart() += Cx;
            } else {
                Matrix<double,ColMajor> Ar = A.realPart();
                Matrix<double,ColMajor> Ai = A.imagPart();
                Matrix<double,ColMajor> Cx = TMV_REAL(alpha)*Ar*B;
                if (A.isconj()) Cx += TMV_IMAG(alpha)*Ai*B;
                else Cx -= TMV_IMAG(alpha)*Ai*B;
                C.realPart() += Cx;

                if (A.isconj()) Cx = -TMV_REAL(alpha)*Ai*B;
                else Cx = TMV_REAL(alpha)*Ai*B;
                Cx += TMV_IMAG(alpha)*Ar*B;
                C.imagPart() += Cx;
            }
        }
    }
    template <int C1, int C2>
    static void DoAddMultMM(
        std::complex<double> alpha,
        const ConstMatrixView<double,C1>& A,
        const ConstMatrixView<std::complex<double>,C2>& B,
        MatrixView<std::complex<double> > C)
    {
        if (TMV_IMAG(alpha) == 0.) {
            Matrix<double,ColMajor> Bx = B.realPart();
            Matrix<double,ColMajor> Cx = TMV_REAL(alpha)*A*Bx;
            C.realPart() += Cx;
            Bx = B.imagPart();
            if (B.isconj()) Cx = -TMV_REAL(alpha)*A*Bx;
            else Cx = TMV_REAL(alpha)*A*Bx;
            C.imagPart() += Cx;
        } else {
            Matrix<double,ColMajor> Br = B.realPart();
            Matrix<double,ColMajor> Bi = B.imagPart();
            Matrix<double,ColMajor> Cx = TMV_REAL(alpha)*A*Br;
            if (B.isconj()) Cx += TMV_IMAG(alpha)*A*Bi;
            else Cx -= TMV_IMAG(alpha)*A*Bi;
            C.realPart() += Cx;

            if (B.isconj()) Cx = -TMV_REAL(alpha)*A*Bi;
            else Cx = TMV_REAL(alpha)*A*Bi;
            Cx += TMV_IMAG(alpha)*A*Br;
            C.imagPart() += Cx;
        }
    }
    template <int C1, int C2>
    static void DoAddMultMM(
        std::complex<double> alpha,
        const ConstMatrixView<double,C1>& A,
        const ConstMatrixView<double,C2>& B,
        MatrixView<std::complex<double> > C)
    {
        Matrix<double,ColMajor> Cx = A*B;
        C += alpha*Cx;
    }
#endif
#endif // INST_DOUBLE
#ifdef TMV_INST_FLOAT
    template <int C1, int C2>
    static void DoAddMultMM(
        float alpha, const ConstMatrixView<float,C1>& A,
        const ConstMatrixView<float,C2>& B, MatrixView<float> C)
    {
        int m = C.colsize();
        int n = C.rowsize();
        int k = A.rowsize();
        int lda = BlasIsCM(A)?A.stepj():A.stepi();
        int ldb = BlasIsCM(B)?B.stepj():B.stepi();
        int ldc = C.stepj();
        float xbeta(1.F);
        BLASNAME(sgemm) (
            BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
            BlasIsCM(B)?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(k),BLASV(alpha),
            BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
            BLASV(xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
    template <int C1, int C2>
    static void DoAddMultMM(
        std::complex<float> alpha,
        const ConstMatrixView<std::complex<float>,C1>& A,
        const ConstMatrixView<std::complex<float>,C2>& B,
        MatrixView<std::complex<float> > C)
    {
        if (BlasIsCM(A) && A.isconj()) {
            const Matrix<std::complex<float>,ColMajor> AA = alpha*A;
            DoAddMultMM(std::complex<float>(1),AA.xView(),B,C);
        } else if (BlasIsCM(B) && B.isconj()) {
            const Matrix<std::complex<float>,ColMajor> BB = alpha*B;
            DoAddMultMM(std::complex<float>(1),A,BB.xView(),C);
        } else {
            int m = C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = BlasIsCM(A)?A.stepj():A.stepi();
            int ldb = BlasIsCM(B)?B.stepj():B.stepi();
            int ldc = C.stepj();
            std::complex<float> xbeta(1.F);
            BLASNAME(cgemm) (
                BLASCM BlasIsCM(A)?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                BlasIsCM(B)?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(k),BLASP(&alpha),
                BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
                BLASP(&xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
        }
    }
#ifdef TMV_INST_MIX
    template <int C1, int C2>
    static void DoAddMultMM(
        std::complex<float> alpha,
        const ConstMatrixView<std::complex<float>,C1>& A,
        const ConstMatrixView<float,C2>& B,
        MatrixView<std::complex<float> > C)
    {
        if (BlasIsCM(A) && !A.isconj() && TMV_IMAG(alpha)==0.) {
            int m = 2*C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = 2*A.stepj();
            int ldb = BlasIsCM(B)?B.stepj():B.stepi();
            int ldc = 2*C.stepj();
            float xalpha(TMV_REAL(alpha));
            float xbeta(1.F);
            BLASNAME(sgemm) (
                BLASCM BLASCH_NT, BlasIsCM(B)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
                BLASP((const float*)(A.cptr())),BLASV(lda),
                BLASP(B.cptr()),BLASV(ldb),
                BLASV(xbeta),BLASP((float*)(C.ptr())),BLASV(ldc) BLAS1 BLAS1);
        } else {
            if (TMV_IMAG(alpha) == 0.F) {
                Matrix<float,ColMajor> Ax = A.realPart();
                Matrix<float,ColMajor> Cx = TMV_REAL(alpha)*Ax*B;
                C.realPart() += Cx;
                Ax = A.imagPart();
                if (A.isconj()) Cx = -TMV_REAL(alpha)*Ax*B;
                else Cx = TMV_REAL(alpha)*Ax*B;
                C.imagPart() += Cx;
            } else {
                Matrix<float,ColMajor> Ar = A.realPart();
                Matrix<float,ColMajor> Ai = A.imagPart();
                Matrix<float,ColMajor> Cx = TMV_REAL(alpha)*Ar*B;
                if (A.isconj()) Cx += TMV_IMAG(alpha)*Ai*B;
                else Cx -= TMV_IMAG(alpha)*Ai*B;
                C.realPart() += Cx;

                if (A.isconj()) Cx = -TMV_REAL(alpha)*Ai*B;
                else Cx = TMV_REAL(alpha)*Ai*B;
                Cx += TMV_IMAG(alpha)*Ar*B;
                C.imagPart() += Cx;
            }
        }
    }
    template <int C1, int C2>
    static void DoAddMultMM(
        std::complex<float> alpha,
        const ConstMatrixView<float,C1>& A,
        const ConstMatrixView<std::complex<float>,C2>& B,
        MatrixView<std::complex<float> > C)
    {
        if (TMV_IMAG(alpha) == 0.F) {
            Matrix<float,ColMajor> Bx = B.realPart();
            Matrix<float,ColMajor> Cx = TMV_REAL(alpha)*A*Bx;
            C.realPart() += Cx;
            Bx = B.imagPart();
            if (B.isconj()) Cx = -TMV_REAL(alpha)*A*Bx;
            else Cx = TMV_REAL(alpha)*A*Bx;
            C.imagPart() += Cx;
        } else {
            Matrix<float,ColMajor> Br = B.realPart();
            Matrix<float,ColMajor> Bi = B.imagPart();
            Matrix<float,ColMajor> Cx = TMV_REAL(alpha)*A*Br;
            if (B.isconj()) Cx += TMV_IMAG(alpha)*A*Bi;
            else Cx -= TMV_IMAG(alpha)*A*Bi;
            C.realPart() += Cx;

            if (B.isconj()) Cx = -TMV_REAL(alpha)*A*Bi;
            else Cx = TMV_REAL(alpha)*A*Bi;
            Cx += TMV_IMAG(alpha)*A*Br;
            C.imagPart() += Cx;
        }
    }
    template <int C1, int C2>
    static void DoAddMultMM(
        std::complex<float> alpha,
        const ConstMatrixView<float,C1>& A,
        const ConstMatrixView<float,C2>& B,
        MatrixView<std::complex<float> > C)
    {
        Matrix<float,ColMajor> Cx = A*B;
        C += alpha*Cx;
    }
#endif
#endif // INST_DOUBLE
#endif

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
#if TMV_OPT <= 2 || defined(BLAS)
        m3.setZero();
        InstAddMultMM(x,m1,m2,m3);
#else
        if (m3.colsize() > 0 && m3.rowsize() > 0) {
            if (m1.rowsize() > 0) {
                if (m3.iscm()) {
                    if (m1.iscm() || m1.isrm()) {
                        if (m2.iscm() || m2.isrm()) {
                            DoMultMM(x,m1,m2,m3);
                        } else {
                            Matrix<T2,ColMajor|NoDivider> m2c = m2;
                            DoMultMM(x,m1,m2c.constView().xView(),m3);
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
            } else {
                m3.setZero();
            }
        }
#endif
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
#if 0
        std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
        std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
        std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
        std::cout<<"Norm(m1) = "<<Norm(m1)<<std::endl;
        std::cout<<"Norm(m2) = "<<Norm(m2)<<std::endl;
        std::cout<<"Norm(m3) = "<<Norm(m3)<<std::endl;
#endif
        if (m3.colsize() > 0 && m3.rowsize() > 0 && m1.rowsize() > 0) {
            if (BlasIsCM(m3)) {
                if (BlasIsCM(m1) || BlasIsRM(m1)) {
                    if (BlasIsCM(m2) || BlasIsRM(m2)) {
                        //std::cout<<"All OK\n";
                        DoAddMultMM(x,m1,m2,m3);
                    } else {
                        Matrix<T2,ColMajor|NoDivider> m2c = m2;
                        //std::cout<<"Copy m2 to ColMajor\n";
                        DoAddMultMM(x,m1,m2c.constView().xView(),m3);
                    }
                } else {
                    Matrix<T1,RowMajor|NoDivider> m1c = m1;
                    //std::cout<<"Copy m1 to RowMajor\n";
                    InstAddMultMM(x,m1c.constView().xView(),m2,m3);
                }
            } else if (BlasIsRM(m3)) {
                //std::cout<<"Transpose\n";
                InstAddMultMM(x,m2.transpose(),m1.transpose(),m3.transpose());
            } else  {
                //std::cout<<"Copy m3 to ColMajor\n";
                Matrix<T3,ColMajor|NoDivider> m3c(m3.colsize(),m3.rowsize());
                InstMultMM(T3(1),m1,m2,m3c.xView());
                InstAddMultXM(x,m3c.constView().xView(),m3);
            }
        }
#if 0
        std::cout<<"Norm(m3) => "<<Norm(m3)<<std::endl;
#endif
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


