///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
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
        MatrixView<T> C)
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

        const ptrdiff_t M = C.colsize();
        const ptrdiff_t N = C.rowsize();
        const ptrdiff_t K = A.rowsize();
        const ptrdiff_t Mb = (M>>6); // = M/64
        const ptrdiff_t Nb = (N>>6); // = N/64
        const ptrdiff_t Kb = (K>>6); // = K/64
        const ptrdiff_t Mc = M < 16 ? 1 : (M>>4); // = M/16
        const ptrdiff_t Nc = N < 16 ? 1 : (N>>4); // = N/16
        const ptrdiff_t Kc = K < 16 ? 1 : (K>>4); // = K/16
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
        MatrixView<T> C)
    { NonBlasMultMM<true>(alpha,A,B,C); }
#ifdef INST_DOUBLE
    template <> 
    void BlasMultMM(
        const double alpha, const GenMatrix<double>& A,
        const GenMatrix<double>& B, MatrixView<double> C)
    {
        int m = C.colsize();
        int n = C.rowsize();
        int k = A.rowsize();
        int lda = BlasIsCM(A)?A.stepj():A.stepi();
        int ldb = BlasIsCM(B)?B.stepj():B.stepi();
        int ldc = C.stepj();
        double xbeta(1);
        BLASNAME(dgemm) (
            BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
            BlasIsCM(B)?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(k),BLASV(alpha),
            BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
            BLASV(xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
    template <> 
    void BlasMultMM(
        const std::complex<double> alpha,
        const GenMatrix<std::complex<double> >& A,
        const GenMatrix<std::complex<double> >& B,
        MatrixView<std::complex<double> > C)
    {
        if (BlasIsCM(A) && A.isconj()) {
            Matrix<std::complex<double>,ColMajor> AA = alpha*A;
            return BlasMultMM(std::complex<double>(1),AA,B,C);
        } else if (BlasIsCM(B) && B.isconj()) {
            Matrix<std::complex<double>,ColMajor> BB = alpha*B;
            return BlasMultMM(std::complex<double>(1),A,BB,C);
        } else {
            int m = C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = BlasIsCM(A)?A.stepj():A.stepi();
            int ldb = BlasIsCM(B)?B.stepj():B.stepi();
            int ldc = C.stepj();
            std::complex<double> xbeta(1);
            BLASNAME(zgemm) (
                BLASCM BlasIsCM(A)?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                BlasIsCM(B)?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
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
        MatrixView<std::complex<double> > C)
    {
        if (BlasIsCM(A) && ((!A.isconj() && TMV_IMAG(alpha)==0.))) {
            int m = 2*C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = 2*A.stepj();
            int ldb = BlasIsCM(B)?B.stepj():B.stepi();
            int ldc = 2*C.stepj();
            double xalpha(TMV_REAL(alpha));
            double xbeta(1);
            BLASNAME(dgemm) (
                BLASCM BLASCH_NT, BlasIsCM(B)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
                BLASP((const double*)(A.cptr())),BLASV(lda),
                BLASP(B.cptr()),BLASV(ldb),
                BLASV(xbeta),BLASP((double*)(C.ptr())),BLASV(ldc) 
                BLAS1 BLAS1);
        } else {
            if (TMV_IMAG(alpha) == 0.) {
                Matrix<double,ColMajor> A1 = A.realPart();
                Matrix<double,ColMajor> C1 = TMV_REAL(alpha)*A1*B;
                C.realPart() += C1;
                A1 = A.imagPart();
                if (A.isconj()) C1 = -TMV_REAL(alpha)*A1*B;
                else C1 = TMV_REAL(alpha)*A1*B;
                C.imagPart() += C1;
            } else {
                Matrix<double,ColMajor> Ar = A.realPart();
                Matrix<double,ColMajor> Ai = A.imagPart();
                Matrix<double,ColMajor> C1 = TMV_REAL(alpha)*Ar*B;
                if (A.isconj()) C1 += TMV_IMAG(alpha)*Ai*B;
                else C1 -= TMV_IMAG(alpha)*Ai*B;
                C.realPart() += C1;

                if (A.isconj()) C1 = -TMV_REAL(alpha)*Ai*B;
                else C1 = TMV_REAL(alpha)*Ai*B;
                C1 += TMV_IMAG(alpha)*Ar*B;
                C.imagPart() += C1;
            }
        }
    }
    template <> 
    void BlasMultMM(
        const std::complex<double> alpha,
        const GenMatrix<double>& A,
        const GenMatrix<std::complex<double> >& B,
        MatrixView<std::complex<double> > C)
    {
        if (TMV_IMAG(alpha) == 0.) {
            Matrix<double,ColMajor> B1 = B.realPart();
            Matrix<double,ColMajor> C1 = TMV_REAL(alpha)*A*B1;
            C.realPart() += C1;
            B1 = B.imagPart();
            if (B.isconj()) C1 = -TMV_REAL(alpha)*A*B1;
            else C1 = TMV_REAL(alpha)*A*B1;
            C.imagPart() += C1;
        } else {
            Matrix<double,ColMajor> Br = B.realPart();
            Matrix<double,ColMajor> Bi = B.imagPart();
            Matrix<double,ColMajor> C1 = TMV_REAL(alpha)*A*Br;
            if (B.isconj()) C1 += TMV_IMAG(alpha)*A*Bi;
            else C1 -= TMV_IMAG(alpha)*A*Bi;
            C.realPart() += C1;

            if (B.isconj()) C1 = -TMV_REAL(alpha)*A*Bi;
            else C1 = TMV_REAL(alpha)*A*Bi;
            C1 += TMV_IMAG(alpha)*A*Br;
            C.imagPart() += C1;
        }
    }
    template <> 
    void BlasMultMM(
        const std::complex<double> alpha,
        const GenMatrix<double>& A, const GenMatrix<double>& B,
        MatrixView<std::complex<double> > C)
    {
        Matrix<double,ColMajor> C1 = A*B;
        C += alpha*C1;
    }
#endif // INST_DOUBLE
#ifdef INST_FLOAT
    template <> 
    void BlasMultMM(
        const float alpha, const GenMatrix<float>& A,
        const GenMatrix<float>& B, MatrixView<float> C)
    {
        int m = C.colsize();
        int n = C.rowsize();
        int k = A.rowsize();
        int lda = BlasIsCM(A)?A.stepj():A.stepi();
        int ldb = BlasIsCM(B)?B.stepj():B.stepi();
        int ldc = C.stepj();
        float xbeta(1);
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
            BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
            BlasIsCM(B)?BLASCH_NT:BLASCH_T,
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
        MatrixView<std::complex<float> > C)
    {
        if (BlasIsCM(A) && A.isconj()) {
            Matrix<std::complex<float> > AA = alpha*A;
            return BlasMultMM(std::complex<float>(1),AA,B,C);
        } else if (BlasIsCM(B) && B.isconj()) {
            Matrix<std::complex<float> > BB = alpha*B;
            return BlasMultMM(std::complex<float>(1),A,BB,C);
        } else {
            int m = C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = BlasIsCM(A)?A.stepj():A.stepi();
            int ldb = BlasIsCM(B)?B.stepj():B.stepi();
            int ldc = C.stepj();
            std::complex<float> xbeta(1);
            //std::cout<<"Before cgemm"<<std::endl;
            BLASNAME(cgemm) (
                BLASCM BlasIsCM(A)?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                BlasIsCM(B)?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(k),BLASP(&alpha),
                BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
                BLASP(&xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
            //std::cout<<"After cgemm"<<std::endl;
        }
    }
    template <> 
    void BlasMultMM(
        const std::complex<float> alpha,
        const GenMatrix<std::complex<float> >& A, const GenMatrix<float>& B,
        MatrixView<std::complex<float> > C)
    {
        if (BlasIsCM(A) && !A.isconj() && (TMV_IMAG(alpha)==0.F)) {
            int m = 2*C.colsize();
            int n = C.rowsize();
            int k = A.rowsize();
            int lda = 2*A.stepj();
            int ldb = BlasIsCM(B)?B.stepj():B.stepi();
            int ldc = 2*C.stepj();
            float xalpha(TMV_REAL(alpha));
            float xbeta(1);
            BLASNAME(sgemm) (
                BLASCM BLASCH_NT, BlasIsCM(B)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
                BLASP((const float*)(A.cptr())),BLASV(lda),
                BLASP(B.cptr()),BLASV(ldb),
                BLASV(xbeta),BLASP((float*)(C.ptr())),BLASV(ldc) 
                BLAS1 BLAS1);
        } else {
            if (TMV_IMAG(alpha) == 0.F) {
                Matrix<float,ColMajor> A1 = A.realPart();
                Matrix<float,ColMajor> C1 = TMV_REAL(alpha)*A1*B;
                C.realPart() += C1;
                A1 = A.imagPart();
                if (A.isconj()) C1 = -TMV_REAL(alpha)*A.imagPart()*B;
                else C1 = TMV_REAL(alpha)*A1*B;
                C.imagPart() += C1;
            } else {
                Matrix<float,ColMajor> Ar = A.realPart();
                Matrix<float,ColMajor> Ai = A.imagPart();
                Matrix<float,ColMajor> C1 = TMV_REAL(alpha)*Ar*B;
                if (A.isconj()) C1 += TMV_IMAG(alpha)*Ai*B;
                else C1 -= TMV_IMAG(alpha)*Ai*B;
                C.realPart() += C1;

                if (A.isconj()) C1 = -TMV_REAL(alpha)*Ai*B;
                else C1 = TMV_REAL(alpha)*Ai*B;
                C1 += TMV_IMAG(alpha)*Ar*B;
                C.imagPart() += C1;
            }
        }
    }
    template <> 
    void BlasMultMM(
        const std::complex<float> alpha,
        const GenMatrix<float>& A,
        const GenMatrix<std::complex<float> >& B,
        MatrixView<std::complex<float> > C)
    {
        if (TMV_IMAG(alpha) == 0.F) {
            Matrix<float,ColMajor> B1 = B.realPart();
            Matrix<float,ColMajor> C1 = TMV_REAL(alpha)*A*B1;
            C.realPart() += C1;
            B1 = B.imagPart();
            if (B.isconj()) C1 = -TMV_REAL(alpha)*A*B1;
            else C1 = TMV_REAL(alpha)*A*B1;
            C.imagPart() += C1;
        } else {
            Matrix<float,ColMajor> Br = B.realPart();
            Matrix<float,ColMajor> Bi = B.imagPart();
            Matrix<float,ColMajor> C1 = TMV_REAL(alpha)*A*Br;
            if (B.isconj()) C1 += TMV_IMAG(alpha)*A*Bi;
            else C1 -= TMV_IMAG(alpha)*A*Bi;
            C.realPart() += C1;

            if (B.isconj()) C1 = -TMV_REAL(alpha)*A*Bi;
            else C1 = TMV_REAL(alpha)*A*Bi;
            C1 += TMV_IMAG(alpha)*A*Br;
            C.imagPart() += C1;
        }
    }
    template <> 
    void BlasMultMM(
        const std::complex<float> alpha,
        const GenMatrix<float>& A, const GenMatrix<float>& B,
        MatrixView<std::complex<float> > C)
    {
        Matrix<float,ColMajor> C1 = A*B;
        C += alpha*C1;
    }
#endif // INST_FLOAT
#endif // BLAS

    template <bool add, class T, class Ta, class Tb> 
    static void DoMultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
#ifdef BLAS
        if (!add) C.setZero();
        BlasMultMM(alpha,A,B,C);
#else
        NonBlasMultMM<add>(alpha,A,B,C);
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
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
        for(ptrdiff_t i=0;i<C.colsize();i++)
            for(ptrdiff_t j=0;j<C.rowsize();j++)
                C2(i,j) = A0.row(i) * B0.col(j);
        C2 *= alpha;
        if (add) C2 += C0;
        //cout<<"Start MultMM: add = "<<add<<", alpha = "<<alpha<<endl;
        //cout<<"A = "<<TMV_Text(A)<<" "<<A0<<endl;
        //cout<<"B = "<<TMV_Text(B)<<" "<<B0<<endl;
        //cout<<"C = "<<TMV_Text(C)<<" "<<C0<<endl;
#endif

        if (C.colsize() > 0 && C.rowsize() > 0) {
            if (A.rowsize() == 0 || alpha == T(0))  {
                if (!add) C.setZero();
            } else if (C.isconj()) {
                MultMM<add>(
                    TMV_CONJ(alpha),A.conjugate(),B.conjugate(),C.conjugate());
            } else if (BlasIsCM(C)) {
                if (!SameStorage(A,C) && (BlasIsCM(A) || BlasIsRM(A))) {
                    if (!SameStorage(B,C) && (BlasIsCM(B) || BlasIsRM(B))) {
                        DoMultMM<add>(alpha,A,B,C);
                    } else {
                        Matrix<T,ColMajor> B2 = alpha*B;
                        DoMultMM<add>(T(1),A,B2,C);
                    }
                } else {
                    Matrix<T,ColMajor> A2 = alpha*A;
                    MultMM<add>(T(1),A2,B,C);
                }
            } else if (BlasIsRM(C)) {
                MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose());
            } else {
                Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
                MultMM<false>(T(1),A,B,C2.view());
                if (add) C += alpha*C2;
                else C = alpha*C2;
            }
        }
#ifdef XDEBUG
        //cout<<"Done MultMM\n";
        //cout<<"C = "<<C<<std::endl;
        //cout<<"C2 = "<<C2<<std::endl;
        //cout<<"Norm(C-C2) = "<<Norm(C-C2)<<std::endl;
        //cout<<"Norm(A0) = "<<Norm(A0)<<std::endl;
        //cout<<"Norm(B0) = "<<Norm(B0)<<std::endl;
        //cout<<"Norm(C0) = "<<Norm(C0)<<std::endl;
        if (!(Norm(C2-C) <= 
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                     (add?Norm(C0):TMV_RealType(T)(0))))) {
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


