///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef TMV_MultMM_Blas_H
#define TMV_MultMM_Blas_H

#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_SumMM.h"
#include "tmv/TMV_SwapM.h"
#include "tmv/TMV_TransposeM.h"

namespace tmv {

#ifdef TMV_INST_DOUBLE
    template <int A1, int A2>
    static void BlasMultMM(
        double alpha, const ConstMatrixView<double,A1>& A,
        const ConstMatrixView<double,A2>& B,
        double beta, MatrixView<double,ColMajor> C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);

        int m = C.colsize();
        int n = C.rowsize();
        int k = A.rowsize();
        int lda = A.iscm()?A.stepj():A.stepi();
        int ldb = B.iscm()?B.stepj():B.stepi();
        int ldc = C.stepj();
        double xbeta(beta);
        BLASNAME(dgemm) (
            BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
            B.iscm()?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(k),BLASV(alpha),
            BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
            BLASV(xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
    template <int A1, int A2>
    static void BlasMultMM(
        std::complex<double> alpha,
        const ConstMatrixView<std::complex<double>,A1>& A,
        const ConstMatrixView<std::complex<double>,A2>& B,
        double beta, MatrixView<std::complex<double>,ColMajor> C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);

        if (A.iscm() && A.isconj()) {
            const Matrix<std::complex<double>,ColMajor> AA = alpha*A;
            BlasMultMM(std::complex<double>(1),AA.xView(),B,beta,C);
        } else if (B.iscm() && B.isconj()) {
            const Matrix<std::complex<double>,ColMajor> BB = alpha*B;
            BlasMultMM(std::complex<double>(1),A,BB.xView(),beta,C);
        } else {
            int m = C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = A.iscm()?A.stepj():A.stepi();
            int ldb = B.iscm()?B.stepj():B.stepi();
            int ldc = C.stepj();
            std::complex<double> xbeta(beta);
            BLASNAME(zgemm) (
                BLASCM A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                B.iscm()?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(k),BLASP(&alpha),
                BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
                BLASP(&xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
        }
    }
#ifdef TMV_INST_MIX
    template <int A1, int A2>
    static void BlasMultMM(
        std::complex<double> alpha,
        const ConstMatrixView<std::complex<double>,A1>& A,
        const ConstMatrixView<double,A2>& B,
        double beta, MatrixView<std::complex<double>,ColMajor> C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);

        if (A.iscm() && ((!A.isconj() && TMV_IMAG(alpha)==0.) || beta == 0)) {
            int m = 2*C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = 2*A.stepj();
            int ldb = B.iscm()?B.stepj():B.stepi();
            int ldc = 2*C.stepj();
            if (beta == 0) {
                double xalpha(1);
                double xbeta(0);
                BLASNAME(dgemm) (
                    BLASCM BLASCH_NT, B.iscm()?BLASCH_NT:BLASCH_T,
                    BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
                    BLASP((const double*)(A.cptr())),BLASV(lda),
                    BLASP(B.cptr()),BLASV(ldb),
                    BLASV(xbeta),BLASP((double*)(C.ptr())),BLASV(ldc) BLAS1 BLAS1);
                if (A.isconj()) C.conjugateSelf();
                C *= alpha;
            } else {
                double xalpha(TMV_REAL(alpha));
                double xbeta(beta);
                BLASNAME(dgemm) (
                    BLASCM BLASCH_NT, B.iscm()?BLASCH_NT:BLASCH_T,
                    BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
                    BLASP((const double*)(A.cptr())),BLASV(lda),
                    BLASP(B.cptr()),BLASV(ldb),
                    BLASV(xbeta),BLASP((double*)(C.ptr())),BLASV(ldc) BLAS1 BLAS1);
            } 
        } else {
            if (TMV_IMAG(alpha) == 0.) {
                Matrix<double,ColMajor> Ax = A.realPart();
                Matrix<double,ColMajor> Cx = TMV_REAL(alpha)*Ax*B;
                if (beta == 0) C.realPart() = Cx;
                else C.realPart() += Cx;
                Ax = A.imagPart();
                if (A.isconj()) Cx = -TMV_REAL(alpha)*Ax*B;
                else Cx = TMV_REAL(alpha)*Ax*B;
                if (beta == 0) C.imagPart() = Cx;
                else C.imagPart() += Cx;
            } else {
                Matrix<double,ColMajor> Ar = A.realPart();
                Matrix<double,ColMajor> Ai = A.imagPart();
                Matrix<double,ColMajor> Cx = TMV_REAL(alpha)*Ar*B;
                if (A.isconj()) Cx += TMV_IMAG(alpha)*Ai*B;
                else Cx -= TMV_IMAG(alpha)*Ai*B;
                if (beta == 0) C.realPart() = Cx;
                else C.realPart() += Cx;

                if (A.isconj()) Cx = -TMV_REAL(alpha)*Ai*B;
                else Cx = TMV_REAL(alpha)*Ai*B;
                Cx += TMV_IMAG(alpha)*Ar*B;
                if (beta == 0) C.imagPart() = Cx;
                else C.imagPart() += Cx;
            }
        }
    }
    template <int A1, int A2>
    static void BlasMultMM(
        std::complex<double> alpha,
        const ConstMatrixView<double,A1>& A,
        const ConstMatrixView<std::complex<double>,A2>& B,
        double beta, MatrixView<std::complex<double>,ColMajor> C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);

        if (TMV_IMAG(alpha) == 0.) {
            Matrix<double,ColMajor> Bx = B.realPart();
            Matrix<double,ColMajor> Cx = TMV_REAL(alpha)*A*Bx;
            if (beta == 0) C.realPart() = Cx;
            else C.realPart() += Cx;
            Bx = B.imagPart();
            if (B.isconj()) Cx = -TMV_REAL(alpha)*A*Bx;
            else Cx = TMV_REAL(alpha)*A*Bx;
            if (beta == 0) C.imagPart() = Cx;
            else C.imagPart() += Cx;
        } else {
            Matrix<double,ColMajor> Br = B.realPart();
            Matrix<double,ColMajor> Bi = B.imagPart();
            Matrix<double,ColMajor> Cx = TMV_REAL(alpha)*A*Br;
            if (B.isconj()) Cx += TMV_IMAG(alpha)*A*Bi;
            else Cx -= TMV_IMAG(alpha)*A*Bi;
            if (beta == 0) C.realPart() = Cx;
            else C.realPart() += Cx;

            if (B.isconj()) Cx = -TMV_REAL(alpha)*A*Bi;
            else Cx = TMV_REAL(alpha)*A*Bi;
            Cx += TMV_IMAG(alpha)*A*Br;
            if (beta == 0) C.imagPart() = Cx;
            else C.imagPart() += Cx;
        }
    }
    template <int A1, int A2>
    static void BlasMultMM(
        std::complex<double> alpha,
        const ConstMatrixView<double,A1>& A,
        const ConstMatrixView<double,A2>& B,
        double beta, MatrixView<std::complex<double>,ColMajor> C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);

        Matrix<double,ColMajor> Cx = A*B;
        if (beta == 0) C = alpha*Cx;
        else C += alpha*Cx;
    }
#endif
#endif // INST_DOUBLE
#ifdef TMV_INST_DOUBLE
    template <int A1, int A2>
    static void BlasMultMM(
        float alpha, const ConstMatrixView<float,A1>& A,
        const ConstMatrixView<float,A2>& B,
        float beta, MatrixView<float,ColMajor> C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);

        int m = C.colsize();
        int n = C.rowsize();
        int k = A.rowsize();
        int lda = A.iscm()?A.stepj():A.stepi();
        int ldb = B.iscm()?B.stepj():B.stepi();
        int ldc = C.stepj();
        float xbeta(beta);
        BLASNAME(sgemm) (
            BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
            B.iscm()?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(k),BLASV(alpha),
            BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
            BLASV(xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
    template <int A1, int A2>
    static void BlasMultMM(
        std::complex<float> alpha,
        const ConstMatrixView<std::complex<float>,A1>& A,
        const ConstMatrixView<std::complex<float>,A2>& B,
        float beta, MatrixView<std::complex<float>,ColMajor> C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);

        if (A.iscm() && A.isconj()) {
            const Matrix<std::complex<float>,ColMajor> AA = alpha*A;
            BlasMultMM(std::complex<float>(1),AA.xView(),B,beta,C);
        } else if (B.iscm() && B.isconj()) {
            const Matrix<std::complex<float>,ColMajor> BB = alpha*B;
            BlasMultMM(std::complex<float>(1),A,BB.xView(),beta,C);
        } else {
            int m = C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = A.iscm()?A.stepj():A.stepi();
            int ldb = B.iscm()?B.stepj():B.stepi();
            int ldc = C.stepj();
            std::complex<float> xbeta(beta);
            BLASNAME(cgemm) (
                BLASCM A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                B.iscm()?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(k),BLASP(&alpha),
                BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
                BLASP(&xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
        }
    }
#ifdef TMV_INST_MIX
    template <int A1, int A2>
    static void BlasMultMM(
        std::complex<float> alpha,
        const ConstMatrixView<std::complex<float>,A1>& A,
        const ConstMatrixView<float,A2>& B,
        float beta, MatrixView<std::complex<float>,ColMajor> C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);

        if (A.iscm() && ((!A.isconj() && TMV_IMAG(alpha)==0.) || beta == 0)) {
            int m = 2*C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = 2*A.stepj();
            int ldb = B.iscm()?B.stepj():B.stepi();
            int ldc = 2*C.stepj();
            if (beta == 0) {
                float xalpha(1);
                float xbeta(0);
                BLASNAME(sgemm) (
                    BLASCM BLASCH_NT, B.iscm()?BLASCH_NT:BLASCH_T,
                    BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
                    BLASP((const float*)(A.cptr())),BLASV(lda),
                    BLASP(B.cptr()),BLASV(ldb),
                    BLASV(xbeta),BLASP((float*)(C.ptr())),BLASV(ldc) BLAS1 BLAS1);
                if (A.isconj()) C.conjugateSelf();
                C *= alpha;
            } else {
                float xalpha(TMV_REAL(alpha));
                float xbeta(beta);
                BLASNAME(sgemm) (
                    BLASCM BLASCH_NT, B.iscm()?BLASCH_NT:BLASCH_T,
                    BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
                    BLASP((const float*)(A.cptr())),BLASV(lda),
                    BLASP(B.cptr()),BLASV(ldb),
                    BLASV(xbeta),BLASP((float*)(C.ptr())),BLASV(ldc) BLAS1 BLAS1);
            } 
        } else {
            if (TMV_IMAG(alpha) == 0.) {
                Matrix<float,ColMajor> Ax = A.realPart();
                Matrix<float,ColMajor> Cx = TMV_REAL(alpha)*Ax*B;
                if (beta == 0) C.realPart() = Cx;
                else C.realPart() += Cx;
                Ax = A.imagPart();
                if (A.isconj()) Cx = -TMV_REAL(alpha)*Ax*B;
                else Cx = TMV_REAL(alpha)*Ax*B;
                if (beta == 0) C.imagPart() = Cx;
                else C.imagPart() += Cx;
            } else {
                Matrix<float,ColMajor> Ar = A.realPart();
                Matrix<float,ColMajor> Ai = A.imagPart();
                Matrix<float,ColMajor> Cx = TMV_REAL(alpha)*Ar*B;
                if (A.isconj()) Cx += TMV_IMAG(alpha)*Ai*B;
                else Cx -= TMV_IMAG(alpha)*Ai*B;
                if (beta == 0) C.realPart() = Cx;
                else C.realPart() += Cx;

                if (A.isconj()) Cx = -TMV_REAL(alpha)*Ai*B;
                else Cx = TMV_REAL(alpha)*Ai*B;
                Cx += TMV_IMAG(alpha)*Ar*B;
                if (beta == 0) C.imagPart() = Cx;
                else C.imagPart() += Cx;
            }
        }
    }
    template <int A1, int A2>
    static void BlasMultMM(
        std::complex<float> alpha,
        const ConstMatrixView<float,A1>& A,
        const ConstMatrixView<std::complex<float>,A2>& B,
        float beta, MatrixView<std::complex<float>,ColMajor> C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);

        if (TMV_IMAG(alpha) == 0.) {
            Matrix<float,ColMajor> Bx = B.realPart();
            Matrix<float,ColMajor> Cx = TMV_REAL(alpha)*A*Bx;
            if (beta == 0) C.realPart() = Cx;
            else C.realPart() += Cx;
            Bx = B.imagPart();
            if (B.isconj()) Cx = -TMV_REAL(alpha)*A*Bx;
            else Cx = TMV_REAL(alpha)*A*Bx;
            if (beta == 0) C.imagPart() = Cx;
            else C.imagPart() += Cx;
        } else {
            Matrix<float,ColMajor> Br = B.realPart();
            Matrix<float,ColMajor> Bi = B.imagPart();
            Matrix<float,ColMajor> Cx = TMV_REAL(alpha)*A*Br;
            if (B.isconj()) Cx += TMV_IMAG(alpha)*A*Bi;
            else Cx -= TMV_IMAG(alpha)*A*Bi;
            if (beta == 0) C.realPart() = Cx;
            else C.realPart() += Cx;

            if (B.isconj()) Cx = -TMV_REAL(alpha)*A*Bi;
            else Cx = TMV_REAL(alpha)*A*Bi;
            Cx += TMV_IMAG(alpha)*A*Br;
            if (beta == 0) C.imagPart() = Cx;
            else C.imagPart() += Cx;
        }
    }
    template <int A1, int A2>
    static void BlasMultMM(
        std::complex<float> alpha,
        const ConstMatrixView<float,A1>& A,
        const ConstMatrixView<float,A2>& B,
        float beta, MatrixView<std::complex<float>,ColMajor> C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);

        Matrix<float,ColMajor> Cx = A*B;
        if (beta == 0) C = alpha*Cx;
        else C += alpha*Cx;
    }
#endif
#endif // INST_DOUBLE


} // namespace tmv

#endif


