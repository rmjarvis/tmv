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
#include "tmv/TMV_SymMatrixArithFunc.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"
#ifdef BLAS
#include "tmv/TMV_TriMatrixArith.h"
#endif

#ifdef XDEBUG
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_MM_BLOCKSIZE TMV_BLOCKSIZE
#define SYM_MM_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define SYM_MM_BLOCKSIZE 64
#define SYM_MM_BLOCKSIZE2 32
#endif

    //
    // MultMM
    //

    template <bool add, class T, class Ta, class Tb> 
    static void RRowMultMM(
        const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.uplo() == Lower);

        const ptrdiff_t N = A.size();
        for(ptrdiff_t j=0;j<N;++j) {
            if (add) C.row(j) += alpha * A.row(j,0,j+1) * B.rowRange(0,j+1);
            else C.row(j) = alpha * A.row(j,0,j+1) * B.rowRange(0,j+1);
            C.rowRange(0,j) += alpha * A.col(j,0,j) ^ B.row(j);
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static void CRowMultMM(
        const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.uplo() == Lower);

        const ptrdiff_t N = A.size();
        for(ptrdiff_t j=N-1;j>=0;--j) {
            if (add) C.row(j) += alpha * A.row(j,j,N) * B.rowRange(j,N);
            else C.row(j) = alpha * A.row(j,j,N) * B.rowRange(j,N);
            C.rowRange(j+1,N) += alpha * A.col(j,j+1,N) ^ B.row(j);
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static inline void RowMultMM(
        const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        if (A.iscm()) CRowMultMM<add>(alpha,A,B,C);
        else RRowMultMM<add>(alpha,A,B,C);
    }

    template <bool add, class T, class Ta, class Tb> 
    static void ColMultMM(
        const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.uplo() == Lower);

        const ptrdiff_t N = C.rowsize();
        for(ptrdiff_t j=0;j<N;++j) 
            if (add) C.col(j) += alpha * A * B.col(j);
            else C.col(j) = alpha * A * B.col(j);
    }

    template <bool add, class T, class Ta, class Tb> 
    static void RecursiveMultMM(
        const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.uplo() == Lower);

        const ptrdiff_t N = A.size();
        if (N <= SYM_MM_BLOCKSIZE2) {
            if (B.isrm() && C.isrm()) RowMultMM<add>(alpha,A,B,C);
            else if (B.iscm() && C.iscm()) ColMultMM<add>(alpha,A,B,C);
            else if (C.colsize() < C.rowsize()) RowMultMM<add>(alpha,A,B,C);
            else ColMultMM<add>(alpha,A,B,C);
        } else {
            ptrdiff_t k = N/2;
            const ptrdiff_t nb = SYM_MM_BLOCKSIZE;
            if (k > nb) k = k/nb*nb;

            // [ A00 A10t ] [ B0 ] = [ A00 B0 + A10t B1 ]
            // [ A10 A11  ] [ B1 ]   [ A10 B0 + A11 B1  ]

            ConstSymMatrixView<Ta> A00 = A.subSymMatrix(0,k);
            ConstSymMatrixView<Ta> A11 = A.subSymMatrix(k,N);
            ConstMatrixView<Ta> A10 = A.subMatrix(k,N,0,k);
            ConstMatrixView<Tb> B0 = B.rowRange(0,k);
            ConstMatrixView<Tb> B1 = B.rowRange(k,N);
            MatrixView<T> C0 = C.rowRange(0,k);
            MatrixView<T> C1 = C.rowRange(k,N);

            RecursiveMultMM<add>(alpha,A00,B0,C0);
            RecursiveMultMM<add>(alpha,A11,B1,C1);
            C1 += alpha * A10 * B0;
            if (A.issym())
                C0 += alpha * A10.transpose() * B1;
            else
                C0 += alpha * A10.adjoint() * B1;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static void NonBlasMultMM(
        const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != T(0));

        if (A.uplo() == Upper)
            if (A.isherm()) NonBlasMultMM<add>(alpha,A.adjoint(),B,C);
            else NonBlasMultMM<add>(alpha,A.transpose(),B,C);
        else if (C.isconj())
            NonBlasMultMM<add>(
                TMV_CONJ(alpha),A.conjugate(),B.conjugate(),C.conjugate());
        else RecursiveMultMM<add>(alpha,A,B,C);
    }

#ifdef BLAS
    template <class T, class Ta, class Tb> 
    static inline void BlasMultMM(
        const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const int beta, MatrixView<T> C)
    { 
        if (beta == 1) NonBlasMultMM<true>(alpha,A,B,C); 
        else NonBlasMultMM<false>(alpha,A,B,C); 
    }
#ifdef INST_DOUBLE
    template <> 
    void BlasMultMM(
        const double alpha, const GenSymMatrix<double>& A,
        const GenMatrix<double>& B, const int beta, MatrixView<double> C)
    {
        int m = C.iscm() ? C.colsize() : C.rowsize();
        int n = C.iscm() ? C.rowsize() : C.colsize();
        int lda = A.stepj();
        int ldb = B.iscm()?B.stepj():B.stepi();
        int ldc = C.iscm()?C.stepj():C.stepi();
        if (beta == 0) C.setZero();
        double xbeta(1);
        BLASNAME(dsymm) (
            BLASCM C.iscm()?BLASCH_L:BLASCH_R,
            A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(B.cptr()),BLASV(ldb),BLASV(xbeta),
            BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
    template <> 
    void BlasMultMM(
        std::complex<double> alpha,
        const GenSymMatrix<std::complex<double> >& A,
        const GenMatrix<std::complex<double> >& B,
        const int beta, MatrixView<std::complex<double> > C)
    {
        int m = C.iscm() ? C.colsize() : C.rowsize();
        int n = C.iscm() ? C.rowsize() : C.colsize();
        int lda = A.stepj();
        int ldb = B.iscm()?B.stepj():B.stepi();
        int ldc = C.iscm()?C.stepj():C.stepi();
        if (beta == 0) C.setZero();
        std::complex<double> xbeta(1);
        if (A.issym())
            BLASNAME(zsymm) (
                BLASCM C.iscm()?BLASCH_L:BLASCH_R,
                A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(B.cptr()),BLASV(ldb),BLASP(&xbeta),
                BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
        else {
            if (!C.iscm()) alpha = TMV_CONJ(alpha);
            BLASNAME(zhemm) (
                BLASCM C.iscm()?BLASCH_L:BLASCH_R,
                A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(B.cptr()),BLASV(ldb),BLASP(&xbeta),
                BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void BlasMultMM(
        const float alpha, const GenSymMatrix<float>& A,
        const GenMatrix<float>& B, const int beta, MatrixView<float> C)
    {
        int m = C.iscm() ? C.colsize() : C.rowsize();
        int n = C.iscm() ? C.rowsize() : C.colsize();
        int lda = A.stepj();
        int ldb = B.iscm()?B.stepj():B.stepi();
        int ldc = C.iscm()?C.stepj():C.stepi();
        if (beta == 0) C.setZero();
        float xbeta(1);
        BLASNAME(ssymm) (
            BLASCM C.iscm()?BLASCH_L:BLASCH_R,
            A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(B.cptr()),BLASV(ldb),BLASV(xbeta),
            BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
    template <> 
    void BlasMultMM(
        std::complex<float> alpha,
        const GenSymMatrix<std::complex<float> >& A,
        const GenMatrix<std::complex<float> >& B,
        const int beta, MatrixView<std::complex<float> > C)
    {
        int m = C.iscm() ? C.colsize() : C.rowsize();
        int n = C.iscm() ? C.rowsize() : C.colsize();
        int lda = A.stepj();
        int ldb = B.iscm()?B.stepj():B.stepi();
        int ldc = C.iscm()?C.stepj():C.stepi();
        if (beta == 0) C.setZero();
        std::complex<float> xbeta(1);
        if (A.issym())
            BLASNAME(csymm) (
                BLASCM C.iscm()?BLASCH_L:BLASCH_R,
                A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(B.cptr()),BLASV(ldb),BLASP(&xbeta),
                BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
        else {
            if (!C.iscm()) alpha = TMV_CONJ(alpha);
            BLASNAME(chemm) (
                BLASCM C.iscm()?BLASCH_L:BLASCH_R,
                A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(B.cptr()),BLASV(ldb),BLASP(&xbeta),
                BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
        }
    }
#endif
    template <class T> 
    static void BlasMultMM(
        const std::complex<T> alpha,
        const GenSymMatrix<std::complex<T> >& A, const GenMatrix<T>& B,
        const int beta, MatrixView<std::complex<T> > C)
    {
        if (TMV_IMAG(alpha) == T(0)) {
            SymMatrix<T,Lower|ColMajor> A1 = A.realPart();
            Matrix<T,ColMajor> C1 = TMV_REAL(alpha)*A1*B;
            if (beta == 0) C.realPart() = C1;
            else C.realPart() += C1;
            if (A.issym()) {
                A1 = A.imagPart();
                if (C.isconj()) C1 = -TMV_REAL(alpha)*A1*B;
                else C1 = TMV_REAL(alpha)*A1*B;
            } else {
                LowerTriMatrixView<T> L = A1.lowerTri();
                L = A.lowerTri().imagPart();
                // A.imagPart() = L - LT
                if (A.lowerTri().isconj() != C.isconj()) {
                    C1 = -TMV_REAL(alpha)*L*B;
                    C1 += TMV_REAL(alpha)*L.transpose()*B;
                } else {
                    C1 = TMV_REAL(alpha)*L*B;
                    C1 -= TMV_REAL(alpha)*L.transpose()*B;
                }
            }
            if (beta == 0) C.imagPart() = C1;
            else C.imagPart() += C1;
        } else {
            SymMatrix<T,Lower|ColMajor> Ar = A.realPart();
            SymMatrix<T,Lower|ColMajor> Ai(A.size());
            LowerTriMatrixView<T> L = Ai.lowerTri();
            Matrix<T,ColMajor> C1 = TMV_REAL(alpha)*Ar*B;
            if (A.issym()) {
                Ai = A.imagPart();
                C1 -= TMV_IMAG(alpha)*Ai*B;
            } else {
                L = A.lowerTri().imagPart();
                if (A.lowerTri().isconj()) L *= T(-1);
                C1 -= TMV_IMAG(alpha)*L*B;
                C1 += TMV_IMAG(alpha)*L.transpose()*B;
            }
            if (beta == 0) C.realPart() = C1;
            else C.realPart() += C1;
            C1 = TMV_IMAG(alpha)*Ar*B;
            if (A.issym()) {
                C1 += TMV_REAL(alpha)*Ai*B;
            } else {
                C1 += TMV_REAL(alpha)*L*B;
                C1 -= TMV_REAL(alpha)*L.transpose()*B;
            }
            if (C.isconj()) C1 *= T(-1);
            if (beta == 0) C.imagPart() = C1;
            else C.imagPart() += C1;
        }
    }
    template <class T> 
    static void BlasMultMM(
        const std::complex<T> alpha,
        const GenSymMatrix<T>& A, const GenMatrix<std::complex<T> >& B,
        const int beta, MatrixView<std::complex<T> > C)
    {
        if (TMV_IMAG(alpha) == T(0)) {
            Matrix<T,ColMajor> B1 = B.realPart();
            Matrix<T,ColMajor> C1 = TMV_REAL(alpha)*A*B1;
            if (beta == 0) C.realPart() = C1;
            else C.realPart() += C1;
            B1 = B.imagPart();
            if (B.isconj()) C1 = -TMV_REAL(alpha)*A*B1;
            else C1 = TMV_REAL(alpha)*A*B1;
            if (beta == 0) C.imagPart() = C1;
            else C.imagPart() += C1;
        } else {
            Matrix<T,ColMajor> Br = B.realPart();
            Matrix<T,ColMajor> Bi = B.imagPart();
            Matrix<T,ColMajor> C1 = TMV_REAL(alpha)*A*Br;
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
    template <class T> 
    static void BlasMultMM(
        const std::complex<T> alpha,
        const GenSymMatrix<T>& A, const GenMatrix<T>& B,
        const int beta, MatrixView<std::complex<T> > C)
    {
        Matrix<T,ColMajor> C1 = A*B;
        if (beta == 0) C = alpha*C1;
        else C += alpha*C1;
    }
#endif // BLAS

    template <bool add, class T, class Ta, class Tb> 
    static void DoMultMM(
        const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != T(0));

#ifdef BLAS
        if (A.isrm())
            DoMultMM<add>(alpha,A.issym()?A.transpose():A.adjoint(),B,C);
        else if (A.isconj())
            DoMultMM<add>(
                TMV_CONJ(alpha),A.conjugate(),B.conjugate(),C.conjugate());
        else if ( !((C.isrm() && C.stepi()>0) || (C.iscm() && C.stepj()>0)) ||
                  (C.iscm() && C.isconj()) || 
                  (C.isrm() && C.isconj()==A.issym()) ) {
            Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
            DoMultMM<false>(T(1),A,B,C2.view());
            if (add) C += alpha*C2;
            else C = alpha*C2;
        } else if (!(A.iscm() && A.stepj()>0)) {
            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                if (A.isherm()) {
                    if (A.uplo() == Upper) {
                        HermMatrix<Ta,Upper|ColMajor> A2 = TMV_REAL(alpha)*A;
                        DoMultMM<add>(T(1),A2,B,C);
                    } else {
                        HermMatrix<Ta,Lower|ColMajor> A2 = TMV_REAL(alpha)*A;
                        DoMultMM<add>(T(1),A2,B,C);
                    }
                } else {
                    if (A.uplo() == Upper) {
                        SymMatrix<Ta,Upper|ColMajor> A2 = TMV_REAL(alpha)*A;
                        DoMultMM<add>(T(1),A2,B,C);
                    } else {
                        SymMatrix<Ta,Lower|ColMajor> A2 = TMV_REAL(alpha)*A;
                        DoMultMM<add>(T(1),A2,B,C);
                    }
                }
            } else {
                if (!A.issym()) {
                    if (A.uplo() == Upper) {
                        // alpha * A is not Hermitian, so can't do 
                        // A2 = alpha * A
                        HermMatrix<Ta,Upper|ColMajor> A2 = A;
                        DoMultMM<add>(alpha,A2,B,C);
                    } else {
                        HermMatrix<Ta,Lower|ColMajor> A2 = A;
                        DoMultMM<add>(alpha,A2,B,C);
                    }
                } else {
                    if (A.uplo() == Upper) {
                        SymMatrix<T,Upper|ColMajor> A2 = alpha*A;
                        DoMultMM<add>(T(1),A2,B,C);
                    } else {
                        SymMatrix<T,Lower|ColMajor> A2 = alpha*A;
                        DoMultMM<add>(T(1),A2,B,C);
                    }
                }
            }
        } else if (!(B.isrm()==C.isrm() && B.iscm()==C.iscm()) ||
                   (isComplex(Tb()) && B.isconj() != C.isconj()) || 
                   !((B.isrm() && B.stepi()>0) || (B.iscm() && B.stepj()>0))) {
            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                if (C.isconj()) {
                    if (C.iscm()) {
                        Matrix<Tb,ColMajor> B2 = TMV_REAL(alpha)*B.conjugate();
                        DoMultMM<add>(T(1),A,B2.conjugate(),C);
                    } else {
                        Matrix<Tb,RowMajor> B2 = TMV_REAL(alpha)*B.conjugate();
                        DoMultMM<add>(T(1),A,B2.conjugate(),C);
                    }
                } else {
                    if (C.iscm()) {
                        Matrix<Tb,ColMajor> B2 = TMV_REAL(alpha)*B;
                        DoMultMM<add>(T(1),A,B2,C);
                    } else {
                        Matrix<Tb,RowMajor> B2 = TMV_REAL(alpha)*B;
                        DoMultMM<add>(T(1),A,B2,C);
                    }
                }
            } else {
                if (C.isconj()) {
                    if (C.iscm()) {
                        Matrix<T,ColMajor> B2 = TMV_CONJ(alpha)*B.conjugate();
                        DoMultMM<add>(T(1),A,B2.conjugate(),C);
                    } else {
                        Matrix<T,RowMajor> B2 = TMV_CONJ(alpha)*B.conjugate();
                        DoMultMM<add>(T(1),A,B2.conjugate(),C);
                    }
                } else {
                    if (C.iscm()) {
                        Matrix<T,ColMajor> B2 = alpha*B;
                        DoMultMM<add>(T(1),A,B2,C);
                    } else {
                        Matrix<T,RowMajor> B2 = alpha*B;
                        DoMultMM<add>(T(1),A,B2,C);
                    }
                }
            }
        } else {
            BlasMultMM(alpha,A,B,add?1:0,C);
        }
#else
        NonBlasMultMM<add>(alpha,A,B,C);
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    static void FullTempMultMM(
        const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        if (C.isrm()) {
            Matrix<T,RowMajor> C2(C.colsize(),C.rowsize());
            DoMultMM<false>(T(1),A,B,C2.view());
            if (add) C += alpha*C2;
            else C = alpha*C2;
        } else {
            Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
            DoMultMM<false>(T(1),A,B,C2.view());
            if (add) C += alpha*C2;
            else C = alpha*C2;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static void BlockTempMultMM(
        const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        const ptrdiff_t N = C.rowsize();
        for(ptrdiff_t j=0;j<N;) {
            ptrdiff_t j2 = TMV_MIN(N,j+SYM_MM_BLOCKSIZE);
            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                if (C.isrm()) {
                    Matrix<Tb,RowMajor> B2 = TMV_REAL(alpha) * B.colRange(j,j2);
                    DoMultMM<add>(T(1),A,B2,C.colRange(j,j2));
                } else {
                    Matrix<Tb,ColMajor> B2 = TMV_REAL(alpha) * B.colRange(j,j2);
                    DoMultMM<add>(T(1),A,B2,C.colRange(j,j2));
                }
            } else {
                if (C.isrm()) {
                    Matrix<T,RowMajor> B2 = alpha * B.colRange(j,j2);
                    DoMultMM<add>(T(1),A,B2,C.colRange(j,j2));
                } else {
                    Matrix<T,ColMajor> B2 = alpha * B.colRange(j,j2);
                    DoMultMM<add>(T(1),A,B2,C.colRange(j,j2));
                }
            }
            j = j2;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    // C (+)= alpha * A * B
    {
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
#ifdef XDEBUG
        //cout<<"Start MultMM: alpha = "<<alpha<<endl;
        //cout<<"A = "<<A.cptr()<<"  "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<"B = "<<B.cptr()<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        //cout<<"C = "<<C.cptr()<<"  "<<TMV_Text(C)<<"  "<<C<<endl;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C0 = C;
        Matrix<T> C2 = alpha*A0*B0;
        if (add) C2 += C0;
#endif

        if (C.colsize() > 0 && C.rowsize() > 0) {
            if (alpha == T(0)) {
                if (!add) C.setZero();
            }
            else if (SameStorage(A,C)) 
                FullTempMultMM<add>(alpha,A,B,C);
            else if (SameStorage(B,C)) 
                if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
                    BlockTempMultMM<add>(alpha,A,B,C);
                else
                    FullTempMultMM<add>(alpha,A,B,C);
            else DoMultMM<add>(alpha, A, B, C);
        }

#ifdef XDEBUG
        //cout<<"Done: C = "<<C<<endl;
        if (!(Norm(C-C2) <= 
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

    template <bool add, class T, class Ta, class Tb> 
    static void BlockTempMultMM(
        const T alpha, const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
        MatrixView<T> C)
    {
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == C.rowsize());
        TMVAssert(A.size() > 0);
        TMVAssert(alpha != T(0));

        const ptrdiff_t N = A.size();

        for(ptrdiff_t j=0;j<N;) {
            ptrdiff_t j2 = TMV_MIN(N,j+SYM_MM_BLOCKSIZE);
            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                if (C.isrm()) {
                    Matrix<Tb,RowMajor> B2(N,j2-j);
                    B2.rowRange(0,j) = TMV_REAL(alpha) * B.subMatrix(0,j,j,j2);
                    B2.rowRange(j,j2) = TMV_REAL(alpha) * B.subSymMatrix(j,j2);
                    B2.rowRange(j2,N) = TMV_REAL(alpha) * B.subMatrix(j2,N,j,j2);
                    DoMultMM<add>(T(1),A,B2.view(),C.colRange(j,j2));
                } else {
                    Matrix<Tb,ColMajor> B2(N,j2-j);
                    B2.rowRange(0,j) = TMV_REAL(alpha) * B.subMatrix(0,j,j,j2);
                    B2.rowRange(j,j2) = TMV_REAL(alpha) * B.subSymMatrix(j,j2);
                    B2.rowRange(j2,N) = TMV_REAL(alpha) * B.subMatrix(j2,N,j,j2);
                    DoMultMM<add>(T(1),A,B2.view(),C.colRange(j,j2));
                }
            } else {
                if (C.isrm()) {
                    Matrix<T,RowMajor> B2(N,j2-j);
                    B2.rowRange(0,j) = alpha * B.subMatrix(0,j,j,j2);
                    B2.rowRange(j,j2) = alpha * B.subSymMatrix(j,j2);
                    B2.rowRange(j2,N) = alpha * B.subMatrix(j2,N,j,j2);
                    DoMultMM<add>(T(1),A,B2.view(),C.colRange(j,j2));
                } else {
                    Matrix<T,ColMajor> B2(N,j2-j);
                    B2.rowRange(0,j) = alpha * B.subMatrix(0,j,j,j2);
                    B2.rowRange(j,j2) = alpha * B.subSymMatrix(j,j2);
                    B2.rowRange(j2,N) = alpha * B.subMatrix(j2,N,j,j2);
                    DoMultMM<add>(T(1),A,B2.view(),C.colRange(j,j2));
                }
            }
            j = j2;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static void FullTempMultMM(
        const T alpha, const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
        MatrixView<T> C)
    {
        if (C.isrm()) {
            Matrix<T,RowMajor> C2(C.colsize(),C.rowsize());
            BlockTempMultMM<false>(T(1),A,B,C2.view());
            if (add) C += alpha*C2;
            else C = alpha*C2;
        } else {
            Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
            BlockTempMultMM<false>(T(1),A,B,C2.view());
            if (add) C += alpha*C2;
            else C = alpha*C2;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenSymMatrix<Ta>& A,
        const GenSymMatrix<Tb>& B, MatrixView<T> C)
    // C (+)= alpha * A * B
    {
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == C.rowsize());
#ifdef XDEBUG
        //cout<<"Start MultMM: alpha = "<<alpha<<endl;
        //cout<<"A = "<<A.cptr()<<"  "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<"B = "<<B.cptr()<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        //cout<<"C = "<<C.cptr()<<"  "<<TMV_Text(C)<<"  "<<C<<endl;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C0 = C;
        Matrix<T> C2 = alpha*A0*B0;
        if (add) C2 += C0;
#endif

        if (A.size() > 0) {
            if (SameStorage(A,C) || SameStorage(B,C))
                FullTempMultMM<add>(alpha,A,B,C);
            else BlockTempMultMM<add>(alpha, A, B, C);
        }

#ifdef XDEBUG
        //cout<<"done: C = "<<C<<endl;
        if (!(Norm(C-C2) <= 
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

#define InstFile "TMV_MultSM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


