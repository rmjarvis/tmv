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



//#define XDEBUG

#include "TMV_Blas.h"
#include "tmv/TMV_MatrixArithFunc.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MatrixArith.h"
#include "TMV_MultMM.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef XDEBUG
#include "tmv/TMV_VectorArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    template <bool add, class T, class Ta, class Tb> 
    static void NonBlasMultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());
        TMVAssert(C.iscm());

        const int M = C.colsize();
        const int N = C.rowsize();
        const int K = A.rowsize();
        const int Mb = (M>>6); // = M/64
        const int Nb = (N>>6); // = N/64
        const int Kb = (K>>6); // = K/64
        const int Mc = M < 16 ? 1 : (M>>4); // = M/16
        const int Nc = N < 16 ? 1 : (N>>4); // = N/16
        const int Kc = K < 16 ? 1 : (K>>4); // = K/16
        const bool twobig = (Mb&&Nb) || (Mb&&Kb) || (Nb&&Kb);

        if ( (M < 16 && N < 16 && K < 16) ||
             (M <= 3 || N <= 3 || K <= 3) ||
             ( ( M < 16 || N < 16 || K < 16 ) &&
               ( !twobig || (Mc * Nc * Kc < 4) ) ) ) {
            // Use a simple algorithm
            if (A.iscm()) {
                if (B.iscm()) 
                    CCCMultMM<add>(alpha,A,B,C);
                else 
                    CRCMultMM<add>(alpha,A,B,C);
            } else {
                if (B.iscm()) 
                    RCCMultMM<add>(alpha,A,B,C);
                else {
                    Matrix<T,ColMajor> B1 = B;
                    RCCMultMM<add>(alpha,A,B1,C);
                }
            }
#ifdef _OPENMP
        } else if (!omp_in_parallel() && (Mb || Nb) &&
                 ( Mc * Nc * Kc >= 64 ) ) {
            OpenMPMultMM<add>(alpha,A,B,C);
#endif
        } else {
            BlockMultMM<add>(alpha,A,B,C);
        }
    }

#ifdef BLAS
    template <class T, class Ta, class Tb> 
    static inline void BlasMultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const int beta, const MatrixView<T>& C)
    {
        if (beta == 0) NonBlasMultMM<false>(alpha,A,B,C); 
        else NonBlasMultMM<true>(alpha,A,B,C); 
    }
#ifdef INST_DOUBLE
    template <> 
    void BlasMultMM(
        const double alpha, const GenMatrix<double>& A,
        const GenMatrix<double>& B, const int beta, const MatrixView<double>& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != 0.);
        TMVAssert(A.ct()==NonConj);
        TMVAssert(B.ct()==NonConj);
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());
        TMVAssert(C.iscm());

        int m = C.colsize();
        int n = C.rowsize();
        int k = A.rowsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();
        int ldc = C.stepj();
        double xbeta(beta);
        BLASNAME(dgemm) (
            BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
            B.iscm()?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(k),BLASV(alpha),
            BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
            BLASV(xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
    template <> 
    void BlasMultMM(
        const std::complex<double> alpha,
        const GenMatrix<std::complex<double> >& A,
        const GenMatrix<std::complex<double> >& B,
        const int beta, const MatrixView<std::complex<double> >& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != 0.);
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());
        TMVAssert(C.iscm());

        if (A.iscm() && A.isconj()) {
            Matrix<std::complex<double>,ColMajor> AA = alpha*A;
            return BlasMultMM(std::complex<double>(1),AA,B,beta,C);
        } else if (B.iscm() && B.isconj()) {
            Matrix<std::complex<double>,ColMajor> BB = alpha*B;
            return BlasMultMM(std::complex<double>(1),A,BB,beta,C);
        } else {
            int m = C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = A.isrm()?A.stepi():A.stepj();
            int ldb = B.isrm()?B.stepi():B.stepj();
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
    template <> 
    void BlasMultMM(
        const std::complex<double> alpha,
        const GenMatrix<std::complex<double> >& A,
        const GenMatrix<double>& B,
        const int beta, const MatrixView<std::complex<double> >& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != 0.);
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());
        TMVAssert(C.iscm());

        if (A.iscm() && ((!A.isconj() && TMV_IMAG(alpha)==0.) || beta == 0)) {
            int m = 2*C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = 2*A.stepj();
            int ldb = B.isrm()?B.stepi():B.stepj();
            int ldc = 2*C.stepj();
            if (beta == 0) {
                double xalpha(1);
                double xbeta(0);
                BLASNAME(dgemm) (
                    BLASCM BLASCH_NT, B.iscm()?BLASCH_NT:BLASCH_T,
                    BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
                    BLASP((const double*)(A.cptr())),BLASV(lda),
                    BLASP(B.cptr()),BLASV(ldb),
                    BLASV(xbeta),BLASP((double*)(C.ptr())),BLASV(ldc) 
                    BLAS1 BLAS1);
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
                    BLASV(xbeta),BLASP((double*)(C.ptr())),BLASV(ldc) 
                    BLAS1 BLAS1);
            } 
        } else {
            if (TMV_IMAG(alpha) == 0.) {
                Matrix<double,ColMajor> A1 = A.realPart();
                Matrix<double,ColMajor> C1 = TMV_REAL(alpha)*A1*B;
                if (beta == 0) C.realPart() = C1;
                else C.realPart() += C1;
                A1 = A.imagPart();
                if (A.isconj()) C1 = -TMV_REAL(alpha)*A1*B;
                else C1 = TMV_REAL(alpha)*A1*B;
                if (beta == 0) C.imagPart() = C1;
                else C.imagPart() += C1;
            } else {
                Matrix<double,ColMajor> Ar = A.realPart();
                Matrix<double,ColMajor> Ai = A.imagPart();
                Matrix<double,ColMajor> C1 = TMV_REAL(alpha)*Ar*B;
                if (A.isconj()) C1 += TMV_IMAG(alpha)*Ai*B;
                else C1 -= TMV_IMAG(alpha)*Ai*B;
                if (beta == 0) C.realPart() = C1;
                else C.realPart() += C1;

                if (A.isconj()) C1 = -TMV_REAL(alpha)*Ai*B;
                else C1 = TMV_REAL(alpha)*Ai*B;
                C1 += TMV_IMAG(alpha)*Ar*B;
                if (beta == 0) C.imagPart() = C1;
                else C.imagPart() += C1;
            }
        }
    }
    template <> 
    void BlasMultMM(
        const std::complex<double> alpha,
        const GenMatrix<double>& A,
        const GenMatrix<std::complex<double> >& B,
        const int beta, const MatrixView<std::complex<double> >& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != 0.);
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());
        TMVAssert(C.iscm());

        if (TMV_IMAG(alpha) == 0.) {
            Matrix<double,ColMajor> B1 = B.realPart();
            Matrix<double,ColMajor> C1 = TMV_REAL(alpha)*A*B1;
            if (beta == 0) C.realPart() = C1;
            else C.realPart() += C1;
            B1 = B.imagPart();
            if (B.isconj()) C1 = -TMV_REAL(alpha)*A*B1;
            else C1 = TMV_REAL(alpha)*A*B1;
            if (beta == 0) C.imagPart() = C1;
            else C.imagPart() += C1;
        } else {
            Matrix<double,ColMajor> Br = B.realPart();
            Matrix<double,ColMajor> Bi = B.imagPart();
            Matrix<double,ColMajor> C1 = TMV_REAL(alpha)*A*Br;
            if (B.isconj()) C1 += TMV_IMAG(alpha)*A*Bi;
            else C1 -= TMV_IMAG(alpha)*A*Bi;
            if (beta == 0) C.realPart() = C1;
            else C.realPart() += C1;

            if (B.isconj()) C1 = -TMV_REAL(alpha)*A*Bi;
            else C1 = TMV_REAL(alpha)*A*Bi;
            C1 += TMV_IMAG(alpha)*A*Br;
            if (beta == 0) C.imagPart() = C1;
            else C.imagPart() += C1;
        }
    }
    template <> 
    void BlasMultMM(
        const std::complex<double> alpha,
        const GenMatrix<double>& A,
        const GenMatrix<double>& B,
        const int beta, const MatrixView<std::complex<double> >& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != 0.);
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());
        TMVAssert(C.iscm());

        Matrix<double,ColMajor> C1 = A*B;
        if (beta == 0) C = alpha*C1;
        else C += alpha*C1;
    }
#endif // INST_DOUBLE
#ifdef INST_FLOAT
    template <> 
    void BlasMultMM(
        const float alpha, const GenMatrix<float>& A,
        const GenMatrix<float>& B, const int beta, const MatrixView<float>& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != 0.F);
        TMVAssert(A.ct()==NonConj);
        TMVAssert(B.ct()==NonConj);
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());
        TMVAssert(C.iscm());

        int m = C.colsize();
        int n = C.rowsize();
        int k = A.rowsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();
        int ldc = C.stepj();
        float xbeta(beta);
        //std::cout<<"Before sgemm"<<std::endl;
        //std::cout<<"A = "<<TMV_Text(A)<<std::endl;
        //std::cout<<"B = "<<TMV_Text(B)<<std::endl;
        //std::cout<<"C = "<<TMV_Text(C)<<std::endl;
        //std::cout<<"A.ptr = "<<A.cptr()<<" .. "<<A.cptr()+(A.colsize()-1)*A.stepi()+(A.rowsize()-1)*A.stepj()<<std::endl;
        //std::cout<<"B.ptr = "<<B.cptr()<<" .. "<<B.cptr()+(B.colsize()-1)*B.stepi()+(B.rowsize()-1)*B.stepj()<<std::endl;
        //std::cout<<"C.ptr = "<<C.cptr()<<" .. "<<C.cptr()+(C.colsize()-1)*C.stepi()+(C.rowsize()-1)*C.stepj()<<std::endl;
        //std::cout<<"m n k = "<<m<<"  "<<n<<"  "<<k<<std::endl;
        //std::cout<<"lda b c = "<<lda<<"  "<<ldb<<"  "<<ldc<<std::endl;
        BLASNAME(sgemm) (
            BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
            B.iscm()?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(k),BLASV(alpha),
            BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
            BLASV(xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
        //std::cout<<"After sgemm"<<std::endl;
    }
    template <> 
    void BlasMultMM(
        const std::complex<float> alpha,
        const GenMatrix<std::complex<float> >& A,
        const GenMatrix<std::complex<float> >& B,
        const int beta, const MatrixView<std::complex<float> >& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != 0.F);
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());
        TMVAssert(C.iscm());

        if (A.iscm() && A.isconj()) {
            Matrix<std::complex<float> > AA = alpha*A;
            return BlasMultMM(std::complex<float>(1),AA,B,beta,C);
        } else if (B.iscm() && B.isconj()) {
            Matrix<std::complex<float> > BB = alpha*B;
            return BlasMultMM(std::complex<float>(1),A,BB,beta,C);
        } else {
            int m = C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = A.isrm()?A.stepi():A.stepj();
            int ldb = B.isrm()?B.stepi():B.stepj();
            int ldc = C.stepj();
            std::complex<float> xbeta(beta);
            //std::cout<<"Before cgemm"<<std::endl;
            BLASNAME(cgemm) (
                BLASCM A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                B.iscm()?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(k),BLASP(&alpha),
                BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
                BLASP(&xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
            //std::cout<<"After cgemm"<<std::endl;
        }
    }
    template <> 
    void BlasMultMM(
        const std::complex<float> alpha,
        const GenMatrix<std::complex<float> >& A,
        const GenMatrix<float>& B,
        const int beta, const MatrixView<std::complex<float> >& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != 0.F);
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());
        TMVAssert(C.iscm());

        if (A.iscm() && !A.isconj() && (TMV_IMAG(alpha)==0.F || beta == 0)) {
            int m = 2*C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = 2*A.stepj();
            int ldb = B.isrm()?B.stepi():B.stepj();
            int ldc = 2*C.stepj();
            if (TMV_IMAG(alpha)==0.F) {
                float xalpha(TMV_REAL(alpha));
                float xbeta(beta);
                BLASNAME(sgemm) (
                    BLASCM BLASCH_NT, B.iscm()?BLASCH_NT:BLASCH_T,
                    BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
                    BLASP((const float*)(A.cptr())),BLASV(lda),
                    BLASP(B.cptr()),BLASV(ldb),
                    BLASV(xbeta),BLASP((float*)(C.ptr())),BLASV(ldc) 
                    BLAS1 BLAS1);
            } else { // beta == 0
                float xalpha(1);
                float xbeta(0);
                BLASNAME(sgemm) (
                    BLASCM BLASCH_NT, B.iscm()?BLASCH_NT:BLASCH_T,
                    BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
                    BLASP((const float*)(A.cptr())),BLASV(lda),
                    BLASP(B.cptr()),BLASV(ldb),
                    BLASV(xbeta),BLASP((float*)(C.ptr())),BLASV(ldc) 
                    BLAS1 BLAS1);
                C *= alpha;
            } 
        } else {
            if (TMV_IMAG(alpha) == 0.F) {
                Matrix<float,ColMajor> A1 = A.realPart();
                Matrix<float,ColMajor> C1 = TMV_REAL(alpha)*A1*B;
                if (beta == 0) C.realPart() = C1;
                else C.realPart() += C1;
                A1 = A.imagPart();
                if (A.isconj()) C1 = -TMV_REAL(alpha)*A.imagPart()*B;
                else C1 = TMV_REAL(alpha)*A1*B;
                if (beta == 0) C.imagPart() = C1;
                else C.imagPart() += C1;
            } else {
                Matrix<float,ColMajor> Ar = A.realPart();
                Matrix<float,ColMajor> Ai = A.imagPart();
                Matrix<float,ColMajor> C1 = TMV_REAL(alpha)*Ar*B;
                if (A.isconj()) C1 += TMV_IMAG(alpha)*Ai*B;
                else C1 -= TMV_IMAG(alpha)*Ai*B;
                if (beta == 0) C.realPart() = C1;
                else C.realPart() += C1;

                if (A.isconj()) C1 = -TMV_REAL(alpha)*Ai*B;
                else C1 = TMV_REAL(alpha)*Ai*B;
                C1 += TMV_IMAG(alpha)*Ar*B;
                if (beta == 0) C.imagPart() = C1;
                else C.imagPart() += C1;
            }
        }
    }
    template <> 
    void BlasMultMM(
        const std::complex<float> alpha,
        const GenMatrix<float>& A,
        const GenMatrix<std::complex<float> >& B,
        const int beta, const MatrixView<std::complex<float> >& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != 0.F);
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());
        TMVAssert(C.iscm());

        if (TMV_IMAG(alpha) == 0.F) {
            Matrix<float,ColMajor> B1 = B.realPart();
            Matrix<float,ColMajor> C1 = TMV_REAL(alpha)*A*B1;
            if (beta == 0) C.realPart() = C1;
            else C.realPart() += C1;
            B1 = B.imagPart();
            if (B.isconj()) C1 = -TMV_REAL(alpha)*A*B1;
            else C1 = TMV_REAL(alpha)*A*B1;
            if (beta == 0) C.imagPart() = C1;
            else C.imagPart() += C1;
        } else {
            Matrix<float,ColMajor> Br = B.realPart();
            Matrix<float,ColMajor> Bi = B.imagPart();
            Matrix<float,ColMajor> C1 = TMV_REAL(alpha)*A*Br;
            if (B.isconj()) C1 += TMV_IMAG(alpha)*A*Bi;
            else C1 -= TMV_IMAG(alpha)*A*Bi;
            if (beta == 0) C.realPart() = C1;
            else C.realPart() += C1;

            if (B.isconj()) C1 = -TMV_REAL(alpha)*A*Bi;
            else C1 = TMV_REAL(alpha)*A*Bi;
            C1 += TMV_IMAG(alpha)*A*Br;
            if (beta == 0) C.imagPart() = C1;
            else C.imagPart() += C1;
        }
    }
    template <> 
    void BlasMultMM(
        const std::complex<float> alpha,
        const GenMatrix<float>& A,
        const GenMatrix<float>& B,
        const int beta, const MatrixView<std::complex<float> >& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != 0.F);
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());
        TMVAssert(C.iscm());

        Matrix<float,ColMajor> C1 = A*B;
        if (beta == 0) C = alpha*C1;
        else C += alpha*C1;
    }
#endif // INST_FLOAT
#endif // BLAS

    template <bool add, class T, class Ta, class Tb> 
    static void DoMultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(C.ct() == NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());
        TMVAssert(C.iscm());

#ifdef BLAS
        if (isComplex(T()) && (isReal(Ta()) || isReal(Tb())))
            BlasMultMM(alpha,A,B,add?1:0,C);
        else if (!((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0))) {
            Matrix<T,ColMajor> A2 = alpha*A;
            DoMultMM<add>(T(1),A2,B,C);
        } else if (!((B.isrm() && B.stepi()>0) || (B.iscm() && B.stepj()>0))) {
            Matrix<T,ColMajor> B2 = alpha*B;
            DoMultMM<add>(T(1),A,B2,C);
        } else {
            BlasMultMM(alpha,A,B,add?1:0,C);
        }
#else
        NonBlasMultMM<add>(alpha,A,B,C);
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    static void FullTempMultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
        // C (+)= alpha * A * B  via a temporary
    {
        Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
        DoMultMM<false>(T(1),A,B,C2.view());

        if (add) C += alpha*C2;
        else C = alpha*C2;
    }

    // Block Temp allows B to be the same storage as C
    template <bool add, class T, class Ta, class Tb> 
    static void BlockTempMultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
        // C (+)= alpha * A * C
    {
        const int N = C.rowsize();
        for (int j=0,j2;j<N;j=j2) {
            j2 = TMV_MIN(N,j+16);
            if (C.isrm()) {
                Matrix<T,ColMajor> B2 = alpha * B.colRange(j,j2);
                DoMultMM<add>(
                    T(1),B2.transpose(),A.transpose(),
                    C.colRange(j,j2).transpose());
            } else {
                Matrix<T,ColMajor> B2 = alpha * B.colRange(j,j2);
                DoMultMM<add>(T(1),A,B2,C.colRange(j,j2));
            }
        }
    }


    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
        // C (+)= alpha * A * B
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<T> B0 = B;
        Matrix<T> C0 = C;
        Matrix<T> C2 = C;
        for(size_t i=0;i<C.colsize();i++)
            for(size_t j=0;j<C.rowsize();j++)
                C2(i,j) = A0.row(i) * B0.col(j);
        C2 *= alpha;
        if (add) C2 += C0;
        cout<<"Start MultMM: add = "<<add<<", alpha = "<<alpha<<endl;
        cout<<"A = "<<TMV_Text(A)<<" "<<A0<<endl;
        cout<<"B = "<<TMV_Text(B)<<" "<<B0<<endl;
        cout<<"C = "<<TMV_Text(C)<<" "<<C0<<endl;
#endif

        if (C.colsize() > 0 && C.rowsize() > 0) {
            if (A.rowsize() == 0 || alpha == T(0))  {
                if (!add) C.setZero();
            } else if (C.isconj()) 
                MultMM<add>(
                    TMV_CONJ(alpha),A.conjugate(),B.conjugate(),C.conjugate());
            else if (C.isrm()) 
                MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose());
            else if (!(B.isrm() || B.iscm())) 
                MultMM<add>(alpha,A,Matrix<Tb,ColMajor>(B),C);
            else if (!(A.isrm() || A.iscm())) 
                if (B.iscm())
                    MultMM<add>(alpha,Matrix<Ta,RowMajor>(A),B,C);
                else
                    MultMM<add>(alpha,Matrix<Ta,ColMajor>(A),B,C);
            else if (!C.iscm()) 
                FullTempMultMM<add>(alpha,A,B,C);
            else if (SameStorage(A,C)) 
                if (SameStorage(B,C)) 
                    FullTempMultMM<add>(alpha,A,B,C);
                else if (C.stepi() == A.stepi() && C.stepj() == A.stepj())
                    BlockTempMultMM<add>(alpha,B.transpose(),A.transpose(),
                                         C.transpose());
                else
                    FullTempMultMM<add>(alpha,A,B,C);
            else if (SameStorage(B,C))
                if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
                    BlockTempMultMM<add>(alpha,A,B,C);
                else
                    FullTempMultMM<add>(alpha,A,B,C);
            else
                DoMultMM<add>(alpha,A,B,C);
        }

#ifdef XDEBUG
        cout<<"Done MultMM\n";
        cout<<"C = "<<C<<std::endl;
        cout<<"C2 = "<<C2<<std::endl;
        cout<<"Norm(C-C2) = "<<Norm(C-C2)<<std::endl;
        cout<<"Norm(A0) = "<<Norm(A0)<<std::endl;
        cout<<"Norm(B0) = "<<Norm(B0)<<std::endl;
        cout<<"Norm(C0) = "<<Norm(C0)<<std::endl;
        if (Norm(C2-C) > 0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                                (add?Norm(C0):TMV_RealType(T)(0)))) {
            cerr<<"MultMM: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C0<<endl;
            cerr<<"--> C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_MultMM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


