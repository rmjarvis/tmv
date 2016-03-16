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
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"
#ifdef BLAS
#include "tmv/TMV_SymMatrixArith.h"
#endif

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_SymMatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_R2K_BLOCKSIZE TMV_BLOCKSIZE
#define SYM_R2K_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define SYM_R2K_BLOCKSIZE 64
#define SYM_R2K_BLOCKSIZE2 1
#endif

    // 
    // Rank2KUpdate
    //

    template <bool ha, bool a1, bool add, class T, class Tx, class Ty> 
    static void RecursiveRank2KUpdate(
        const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
        SymMatrixView<T> A)
    {
        TMVAssert(A.size() == x.colsize());
        TMVAssert(A.size() == y.colsize());
        TMVAssert(x.rowsize() == y.rowsize());
        TMVAssert(alpha != T(0));
        TMVAssert(x.colsize() > 0);
        TMVAssert(x.rowsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.uplo() == Lower);
        TMVAssert(ha == A.isherm());
        TMVAssert(a1 == (alpha == T(1)));

        const ptrdiff_t nb = SYM_R2K_BLOCKSIZE;
        ptrdiff_t N = A.size();

        if (N <= SYM_R2K_BLOCKSIZE2) {
            if (N == 1) {
                T temp = x.row(0) * (ha ? y.row(0).conjugate() : y.row(0));
                if (!a1) temp *= alpha;
#ifdef TMVFLDEBUG
                TMVAssert(A.ptr() >= A._first);
                TMVAssert(A.ptr() < A._last);
#endif
                if (ha)
                    if (add) *(A.ptr()) += TMV_RealType(T)(2) * TMV_REAL(temp);
                    else *(A.ptr()) = TMV_RealType(T)(2) * TMV_REAL(temp);
                else
                    if (add) *(A.ptr()) += TMV_RealType(T)(2) * temp;
                    else *(A.ptr()) = TMV_RealType(T)(2) * temp;
            } else {
                if (x.isrm() && y.isrm()) {
                    if (A.isrm()) {
                        for (ptrdiff_t i=0;i<N;++i) {
                            if (add) 
                                A.row(i,0,i+1) += alpha * x.row(i) * 
                                    (ha ?
                                     y.rowRange(0,i+1).adjoint() :
                                     y.rowRange(0,i+1).transpose());
                            else 
                                A.row(i,0,i+1) = alpha * x.row(i) * 
                                    (ha ?
                                     y.rowRange(0,i+1).adjoint() :
                                     y.rowRange(0,i+1).transpose());
                            A.row(i,0,i+1) += y.row(i) * 
                                (ha ?
                                 (TMV_CONJ(alpha)*x.rowRange(0,i+1).adjoint()) :
                                 (alpha*x.rowRange(0,i+1).transpose()));
                        }
                    }
                    else {
                        for (ptrdiff_t j=0;j<N;++j) {
                            if (add) 
                                A.col(j,j,N) += alpha * x.rowRange(j,N) * 
                                    (ha ? y.row(j).conjugate() : y.row(j));
                            else 
                                A.col(j,j,N) = alpha * x.rowRange(j,N) * 
                                    (ha ? y.row(j).conjugate() : y.row(j));
                            A.col(j,j,N) += y.rowRange(j,N) * 
                                (ha ?
                                 (TMV_CONJ(alpha)*x.row(j).conjugate()) :
                                 (alpha*x.row(j)));
                        }
                    }
                } else { // x,y not row major
                    for (ptrdiff_t i=0;i<N;++i) {
                        Rank2Update<add>(alpha,x.col(i),y.col(i),A);
                    }
                }
            }
        } else { // Not <= BLOCKSIZE2, so do recurse...
            ptrdiff_t k = N/2;
            if (k > nb) k = k/nb*nb;

            RecursiveRank2KUpdate<ha,a1,add>(
                alpha,x.rowRange(0,k),y.rowRange(0,k), A.subSymMatrix(0,k));

            if (add) 
                A.subMatrix(k,N,0,k) += alpha * x.rowRange(k,N) * 
                    (ha ? y.rowRange(0,k).adjoint() :
                     y.rowRange(0,k).transpose());
            else 
                A.subMatrix(k,N,0,k) = alpha * x.rowRange(k,N) * 
                    (ha ? y.rowRange(0,k).adjoint() :
                     y.rowRange(0,k).transpose());

            A.subMatrix(k,N,0,k) += y.rowRange(k,N) * 
                (ha ? (TMV_CONJ(alpha)*x.rowRange(0,k).adjoint()) : 
                 (alpha * x.rowRange(0,k).transpose()));

            RecursiveRank2KUpdate<ha,a1,add>(
                alpha,x.rowRange(k,N),y.rowRange(k,N), A.subSymMatrix(k,N));
        }
    }

    template <bool ha, bool a1, bool add, class T, class Tx, class Ty> 
    static void RecursiveInPlaceRank2KUpdate(
        const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
        SymMatrixView<T> A)
    {
        TMVAssert(SameStorage(x,A) || SameStorage(y,A));
        TMVAssert(A.size() > 0);
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.ct() == NonConj);

        ptrdiff_t N = A.size();
        if (N == 1) {
            Tx x00 = x(0,0);
            Ty y00 = y(0,0);
#ifdef TMVFLDEBUG
            TMVAssert(A.ptr() >= A._first);
            TMVAssert(A.ptr() < A._last);
#endif
            if (ha) {
                TMV_RealType(T) temp = TMV_RealType(T)(2) *
                    (a1 ? TMV_REAL(x00*TMV_CONJ(y00)) :
                     TMV_REAL(alpha*(x00*(TMV_CONJ(y00)))));
                if (add) *A.ptr() += temp;
                else *A.ptr() = temp;
            } else {
                T temp = TMV_RealType(T)(2) * 
                    (a1 ? (x00 * y00) : (alpha * (x00 * y00)));
                if (add) *A.ptr() += temp;
                else *A.ptr() = temp;
            }
        } else {
            const ptrdiff_t k = N/2;
            const ConstMatrixView<Tx> x00 = x.subMatrix(0,k,0,k);
            const ConstMatrixView<Tx> x10 = x.subMatrix(k,N,0,k);
            const ConstMatrixView<Tx> x01 = x.subMatrix(0,k,k,N);
            const ConstMatrixView<Tx> x11 = x.subMatrix(k,N,k,N);
            const ConstMatrixView<Ty> y00 = y.subMatrix(0,k,0,k);
            const ConstMatrixView<Ty> y10 = y.subMatrix(k,N,0,k);
            const ConstMatrixView<Ty> y01 = y.subMatrix(0,k,k,N);
            const ConstMatrixView<Ty> y11 = y.subMatrix(k,N,k,N);
            SymMatrixView<T> A00 = A.subSymMatrix(0,k);
            SymMatrixView<T> A11 = A.subSymMatrix(k,N);
            MatrixView<T> A10 = A.subMatrix(k,N,0,k);

            Matrix<T> tempA10 = x10 * (ha ? y00.adjoint() : y00.transpose());
            tempA10 += x11 * (ha ? y01.adjoint() : y01.transpose());
            if (!a1) tempA10 *= alpha;
            T ca = ha ? TMV_CONJ(alpha) : alpha;
            tempA10 += ca * y10 * (ha ? x00.adjoint() : x00.transpose());
            tempA10 += ca * y11 * (ha ? x01.adjoint() : x01.transpose());
            RecursiveInPlaceRank2KUpdate<ha,a1,add>(alpha,x11,y11,A11);
            RecursiveRank2KUpdate<ha,a1,true>(alpha,x10,y10,A11);
            RecursiveInPlaceRank2KUpdate<ha,a1,add>(alpha,x00,y00,A00);
            RecursiveRank2KUpdate<ha,a1,true>(alpha,x01,y01,A00);

            if (add) A10 += tempA10;
            else A10 = tempA10;

        }
    }


    template <bool add, class T, class Tx, class Ty> 
    static inline void InPlaceRank2KUpdate(
        const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
        SymMatrixView<T> A)
    {
        TMVAssert(SameStorage(x,A) || SameStorage(y,A));
        TMVAssert(A.size() > 0);
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.ct() == NonConj);

        if (A.isherm())
            if (alpha == T(1))
                RecursiveInPlaceRank2KUpdate<true,true,add>(alpha,x,y,A); 
            else
                RecursiveInPlaceRank2KUpdate<true,false,add>(alpha,x,y,A); 
        else
            if (alpha == T(1))
                RecursiveInPlaceRank2KUpdate<false,true,add>(alpha,x,y,A); 
            else
                RecursiveInPlaceRank2KUpdate<false,false,add>(alpha,x,y,A); 
    } 
    template <bool add, class T, class Tx, class Ty> 
    static void NonBlasRank2KUpdate(
        const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
        SymMatrixView<T> A)
    { 
        if (A.uplo() == Upper)
            return NonBlasRank2KUpdate<add>(
                alpha,x,y,A.issym()?A.transpose():A.adjoint());
        else if (A.isconj())
            return NonBlasRank2KUpdate<add>(
                TMV_CONJ(alpha),x.conjugate(),y.conjugate(),A.conjugate());
        else if (SameStorage(x,A) || SameStorage(y,A)) {
            const ptrdiff_t N = A.size();
            TMVAssert(x.colsize() == N);
            TMVAssert(y.colsize() == N);
            const ptrdiff_t K = x.rowsize();
            TMVAssert(y.rowsize() == K);

            if (K >= N) {
                InPlaceRank2KUpdate<add>(
                    alpha,x.colRange(0,N),y.colRange(0,N),A);
                if (K > N) {
                    NonBlasRank2KUpdate<add>(
                        alpha,x.colRange(N,K),y.colRange(N,K),A);
                }
            } else { 
                NonBlasRank2KUpdate<add>(
                    alpha,x.rowRange(K,N),y.rowRange(K,N), A.subSymMatrix(K,N));
                const bool ha = A.isherm();
                if ((x.stepi() < x.stepj()) == (A.stepi() < A.stepj())) {
                    if ((y.stepi() < y.stepj()) == (A.stepi() < A.stepj())) {
                        // Then both x and y overlap with A(K:N,0:K)
                        // Need a temporary
                        if (A.iscm()) {
                            if (ha) {
                                Matrix<T,ColMajor> temp = 
                                    alpha * x.rowRange(K,N) *
                                    y.rowRange(0,K).adjoint();
                                temp += TMV_CONJ(alpha) * y.rowRange(K,N) *
                                    x.rowRange(0,K).adjoint();
                                if (add) A.subMatrix(K,N,0,K) += temp;
                                else A.subMatrix(K,N,0,K) = temp;
                            } else {
                                Matrix<T,ColMajor> temp = 
                                    alpha * x.rowRange(K,N) * 
                                    y.rowRange(0,K).transpose();
                                temp += alpha * y.rowRange(K,N) *
                                    x.rowRange(0,K).transpose();
                                if (add) A.subMatrix(K,N,0,K) += temp;
                                else A.subMatrix(K,N,0,K) = temp;
                            }
                        } else {
                            if (ha) {
                                Matrix<T,RowMajor> temp = 
                                    alpha * x.rowRange(K,N) *
                                    y.rowRange(0,K).adjoint();
                                temp += TMV_CONJ(alpha) * y.rowRange(K,N) *
                                    x.rowRange(0,K).adjoint();
                                if (add) A.subMatrix(K,N,0,K) += temp;
                                else A.subMatrix(K,N,0,K) = temp;
                            } else {
                                Matrix<T,RowMajor> temp = 
                                    alpha * x.rowRange(K,N) *
                                    y.rowRange(0,K).transpose();
                                temp += alpha * y.rowRange(K,N) *
                                    x.rowRange(0,K).transpose();
                                if (add) A.subMatrix(K,N,0,K) += temp;
                                else A.subMatrix(K,N,0,K) = temp;
                            }
                        }
                    } else {
                        MultMM<true>(
                            alpha, x.rowRange(K,N), 
                            ha ? y.rowRange(0,K).adjoint() : 
                            y.rowRange(0,K).transpose(),
                            A.subMatrix(K,N,0,K) );
                        MultMM<add>(
                            ha ? TMV_CONJ(alpha) : alpha, y.rowRange(K,N),
                            ha ? x.rowRange(0,K).adjoint() :
                            x.rowRange(0,K).transpose(),
                            A.subMatrix(K,N,0,K) );
                    }
                }
                else {
                    MultMM<true>(
                        ha ? TMV_CONJ(alpha) : alpha, y.rowRange(K,N),
                        ha ? x.rowRange(0,K).adjoint() :
                        x.rowRange(0,K).transpose(),
                        A.subMatrix(K,N,0,K) );
                    MultMM<add>(
                        alpha, x.rowRange(K,N),
                        ha ? y.rowRange(0,K).adjoint() :
                        y.rowRange(0,K).transpose(),
                        A.subMatrix(K,N,0,K) );
                }
                InPlaceRank2KUpdate<add>(
                    alpha,x.rowRange(0,K),y.rowRange(0,K),A.subSymMatrix(0,K));
            }
        }
        else 
            if (A.isherm())
                if (alpha == T(1))
                    RecursiveRank2KUpdate<true,true,add>(alpha,x,y,A); 
                else
                    RecursiveRank2KUpdate<true,false,add>(alpha,x,y,A); 
            else
                if (alpha == T(1))
                    RecursiveRank2KUpdate<false,true,add>(alpha,x,y,A); 
                else
                    RecursiveRank2KUpdate<false,false,add>(alpha,x,y,A); 
    }

#ifdef BLAS
    template <class T, class Tx, class Ty> 
    static inline void BlasRank2KUpdate(
        const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
        SymMatrixView<T> A)
    { NonBlasRank2KUpdate<true>(alpha,x,y,A); }
#ifdef INST_DOUBLE
    template <> 
    void BlasRank2KUpdate(
        const double alpha, const GenMatrix<double>& x,
        const GenMatrix<double>& y, 
        SymMatrixView<double> A)
    {
        int n=A.size();
        int k=x.rowsize();
        int ldx=x.iscm()?x.stepj():x.stepi();
        int ldy=y.iscm()?y.stepj():y.stepi();
        double beta(1);
        int lda=A.stepj();
        BLASNAME(dsyr2k) (
            BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, 
            x.iscm()?BLASCH_NT:BLASCH_T,BLASV(n),BLASV(k),BLASV(alpha),
            BLASP(x.cptr()),BLASV(ldx),BLASP(y.cptr()),BLASV(ldy),
            BLASV(beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
    }
    template <> 
    void BlasRank2KUpdate(
        const std::complex<double> alpha,
        const GenMatrix<std::complex<double> >& x, 
        const GenMatrix<std::complex<double> >& y,
        SymMatrixView<std::complex<double> > A)
    {
        int n=A.size();
        int k=x.rowsize();
        int ldx=x.iscm()?x.stepj():x.stepi();
        int ldy=y.iscm()?y.stepj():y.stepi();
        int lda=A.stepj();
        if (A.isherm()) {
            double beta(1);
            BLASNAME(zher2k) (
                BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, 
                x.iscm()?BLASCH_NT:BLASCH_CT,BLASV(n),BLASV(k),BLASP(&alpha),
                BLASP(x.cptr()),BLASV(ldx),BLASP(y.cptr()),BLASV(ldy),
                BLASV(beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
        }
        else {
            std::complex<double> beta(1);
            BLASNAME(zsyr2k) (
                BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, 
                x.iscm()?BLASCH_NT:BLASCH_T,BLASV(n),BLASV(k),BLASP(&alpha),
                BLASP(x.cptr()),BLASV(ldx),BLASP(y.cptr()),BLASV(ldy),
                BLASP(&beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasRank2KUpdate(
        const std::complex<double> alpha, 
        const GenMatrix<std::complex<double> >& x,
        const GenMatrix<double>& y, 
        SymMatrixView<std::complex<double> > A)
    {
        Matrix<double,ColMajor> x1 = x.realPart();
        SymMatrix<double,Lower|ColMajor> A1(A.size(),0.);
        BlasRank2KUpdate(1.,x1,y,A1.view());
        A += alpha*A1;
        x1 = x.imagPart();
        BlasRank2KUpdate(1.,x1,y,A1.view());
        A += std::complex<double>(0,1)*alpha*A1;
    }
    template <> 
    void BlasRank2KUpdate(
        const std::complex<double> alpha, const GenMatrix<double>& x,
        const GenMatrix<std::complex<double> >& y, 
        SymMatrixView<std::complex<double> > A)
    {
        Matrix<double,RowMajor> y1 = y.realPart();
        SymMatrix<double,Lower|ColMajor> A1(A.size(),0.);
        BlasRank2KUpdate(1.,x,y1,A1.view());
        A += alpha*A1;
        y1 = y.imagPart();
        BlasRank2KUpdate(1.,x,y1,A1.view());
        A += std::complex<double>(0,1)*alpha*A1;
    }
    template <> 
    void BlasRank2KUpdate(
        const std::complex<double> alpha, const GenMatrix<double>& x,
        const GenMatrix<double>& y, 
        SymMatrixView<std::complex<double> > A)
    {
        SymMatrix<double,Lower|ColMajor> A1(A.size(),0.);
        BlasRank2KUpdate(1.,x,y,A1.view());
        A += alpha*A1;
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void BlasRank2KUpdate(
        const float alpha, const GenMatrix<float>& x,
        const GenMatrix<float>& y, 
        SymMatrixView<float> A)
    {
        int n=A.size();
        int k=x.rowsize();
        int ldx=x.iscm()?x.stepj():x.stepi();
        int ldy=y.iscm()?y.stepj():y.stepi();
        float beta(1);
        int lda=A.stepj();
        BLASNAME(ssyr2k) (
            BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, 
            x.iscm()?BLASCH_NT:BLASCH_T,BLASV(n),BLASV(k),BLASV(alpha),
            BLASP(x.cptr()),BLASV(ldx),BLASP(y.cptr()),BLASV(ldy),
            BLASV(beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
    }
    template <> 
    void BlasRank2KUpdate(
        const std::complex<float> alpha,
        const GenMatrix<std::complex<float> >& x, 
        const GenMatrix<std::complex<float> >& y,
        SymMatrixView<std::complex<float> > A)
    {
        int n=A.size();
        int k=x.rowsize();
        int ldx=x.iscm()?x.stepj():x.stepi();
        int ldy=y.iscm()?y.stepj():y.stepi();
        int lda=A.stepj();
        if (A.isherm()) {
            float beta(1);
            BLASNAME(cher2k) (
                BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, 
                x.iscm()?BLASCH_NT:BLASCH_CT,BLASV(n),BLASV(k),BLASP(&alpha),
                BLASP(x.cptr()),BLASV(ldx),BLASP(y.cptr()),BLASV(ldy),
                BLASV(beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
        }
        else {
            std::complex<float> beta(1);
            BLASNAME(csyr2k) (
                BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, 
                x.iscm()?BLASCH_NT:BLASCH_T,BLASV(n),BLASV(k),BLASP(&alpha),
                BLASP(x.cptr()),BLASV(ldx),BLASP(y.cptr()),BLASV(ldy),
                BLASP(&beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasRank2KUpdate(
        const std::complex<float> alpha, 
        const GenMatrix<std::complex<float> >& x,
        const GenMatrix<float>& y, 
        SymMatrixView<std::complex<float> > A)
    {
        Matrix<float,ColMajor> x1 = x.realPart();
        SymMatrix<float,Lower|ColMajor> A1(A.size(),0.F);
        BlasRank2KUpdate(1.F,x1,y,A1.view());
        A += alpha*A1;
        x1 = x.imagPart();
        BlasRank2KUpdate(1.F,x1,y,A1.view());
        A += std::complex<float>(0,1)*alpha*A1;
    }
    template <> 
    void BlasRank2KUpdate(
        const std::complex<float> alpha, const GenMatrix<float>& x,
        const GenMatrix<std::complex<float> >& y, 
        SymMatrixView<std::complex<float> > A)
    {
        Matrix<float,RowMajor> y1 = y.realPart();
        SymMatrix<float,Lower|ColMajor> A1(A.size(),0.F);
        BlasRank2KUpdate(1.F,x,y1,A1.view());
        A += alpha*A1;
        y1 = y.imagPart();
        BlasRank2KUpdate(1.F,x,y1,A1.view());
        A += std::complex<float>(0,1)*alpha*A1;
    }
    template <> 
    void BlasRank2KUpdate(
        const std::complex<float> alpha, const GenMatrix<float>& x,
        const GenMatrix<float>& y, 
        SymMatrixView<std::complex<float> > A)
    {
        SymMatrix<float,Lower|ColMajor> A1(A.size(),0.F);
        BlasRank2KUpdate(1.F,x,y,A1.view());
        A += alpha*A1;
    }
#endif
#endif // BLAS

    template <bool add, class T, class Tx, class Ty> 
    void Rank2KUpdate(
        const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
        SymMatrixView<T> A)
    // if A is sym:  A (+)= alpha * (x ^ y + y ^ x)
    // if A is herm: A (+)= alpha * x ^ y* + conj(alpha) * y ^ x*
    {
        TMVAssert(A.size() == x.colsize());
        TMVAssert(A.size() == y.colsize());
        TMVAssert(x.rowsize() == y.rowsize());

#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<Tx> x0 = x;
        Matrix<Ty> y0 = y;
        Matrix<T> A2 = A;
        if (A.isherm()) {
            if (add) A2 += alpha*x*y.adjoint();
            else A2 = alpha*x*y.adjoint();
            A2 += TMV_CONJ(alpha)*y*x.adjoint();
        }
        else {
            if (add) A2 += alpha*x*y.transpose();
            else A2 = alpha*x*y.transpose();
            A2 += alpha*y*x.transpose();
        }
#endif

        if (alpha != T(0) && A.size() > 0) {
            if (x.rowsize() == 1)
                Rank2Update<add>(alpha,x.col(0),y.col(0),A);
#ifdef BLAS
            else if (!A.iscm() && A.isrm())
                Rank2KUpdate<add>(
                    alpha,x,y,A.issym()?A.transpose():A.adjoint());
            else if (A.isconj()) 
                Rank2KUpdate<add>(
                    TMV_CONJ(alpha),x.conjugate(),y.conjugate(),A.conjugate());
            else if (A.iscm() && A.stepj()>0) {
                if (!add) A.setZero();
                if (!((x.isrm() && x.stepi()>0) || (x.iscm() && x.stepj()>0)) ||
                    (!A.issym() && x.iscm() == x.isconj()) ||
                    (!A.isherm() && x.isconj()) || SameStorage(x,A)) {
                    if (!((y.isrm() && y.stepi()>0) ||
                          (y.iscm() && y.stepj()>0)) ||
                        y.isconj() || SameStorage(y,A)) {
                        if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                            Matrix<Tx,ColMajor> xx = TMV_REAL(alpha)*x;
                            Matrix<Ty,ColMajor> yy = y;
                            BlasRank2KUpdate(T(1),xx,yy,A);
                        } else {
                            Matrix<T,ColMajor> xx = alpha*x;
                            Matrix<Ty,ColMajor> yy = y;
                            BlasRank2KUpdate(T(1),xx,yy,A);
                        }
                    } else {
                        if (y.iscm()) {
                            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                                Matrix<Tx,ColMajor> xx = TMV_REAL(alpha)*
                                    (y.isconj() ? x.conjugate() : x.view());
                                BlasRank2KUpdate(
                                    T(1),
                                    (y.isconj() ? xx.conjugate() : xx.view()),
                                    y,A);
                            } else {
                                Matrix<T,ColMajor> xx = alpha*
                                    (y.isconj() ? x.conjugate() : x.view());
                                BlasRank2KUpdate(
                                    T(1),
                                    (y.isconj() ? xx.conjugate() : xx.view()),
                                    y,A);
                            }
                        } else {
                            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                                Matrix<Tx,RowMajor> xx = TMV_REAL(alpha)*
                                    (y.isconj() ? x.conjugate() : x.view());
                                BlasRank2KUpdate(
                                    T(1),
                                    (y.isconj() ? xx.conjugate() : xx.view()),
                                    y,A);
                            } else {
                                Matrix<T,RowMajor> xx = alpha*
                                    (y.isconj() ? x.conjugate() : x.view());
                                BlasRank2KUpdate(
                                    T(1),
                                    (y.isconj() ? xx.conjugate() : xx.view()),
                                    y,A);
                            }
                        }
                    }
                } else {
                    if (!((y.isrm() && y.stepi()>0) || 
                          (y.iscm() && y.stepj()>0)) ||
                        !((x.isrm()==y.isrm()) && (x.iscm()==y.iscm())) ||
                        y.isconj() || SameStorage(y,A)) {
                        if (x.iscm()) {
                            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                                Matrix<Ty,ColMajor> yy = TMV_REAL(alpha)*
                                    (x.isconj() ? y.conjugate() : y.view());
                                BlasRank2KUpdate(
                                    T(1),x,
                                    (x.isconj() ? yy.conjugate() : yy.view()),
                                    A);
                            } else {
                                Matrix<T,ColMajor> yy = alpha*
                                    (x.isconj() ? y.conjugate() : y.view());
                                BlasRank2KUpdate(
                                    T(1),x,
                                    (x.isconj() ? yy.conjugate() : yy.view()),
                                    A);
                            }
                        } else {
                            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                                Matrix<Ty,RowMajor> yy = TMV_REAL(alpha)*
                                    (x.isconj() ? y.conjugate() : y.view());
                                BlasRank2KUpdate(
                                    T(1),x,
                                    (x.isconj() ? yy.conjugate() : yy.view()),
                                    A);
                            } else {
                                Matrix<T,RowMajor> yy = alpha*
                                    (x.isconj() ? y.conjugate() : y.view());
                                BlasRank2KUpdate(
                                    T(1),x,
                                    (x.isconj() ? yy.conjugate() : yy.view()),
                                    A);
                            }
                        }
                    } else {
                        BlasRank2KUpdate(alpha,x,y,A);
                    }
                }
            } else {
                if (A.isherm()) {
                    HermMatrix<T,Lower|ColMajor> AA(A.size());
                    Rank2KUpdate<false>(alpha,x,y,AA.view());
                    if (add) A += AA;
                    else A = AA;
                } else {
                    SymMatrix<T,Lower|ColMajor> AA(A.size());
                    Rank2KUpdate<false>(alpha,x,y,AA.view());
                    if (add) A += AA;
                    else A = AA;
                }
            }
#else
            else NonBlasRank2KUpdate<add>(alpha,x,y,A);
#endif
        }

#ifdef XDEBUG
        TMVAssert(A.isHermOK());
        if (!(Norm(A-A2) < 0.001*(TMV_ABS(alpha)*Norm(x0)*Norm(y0)+
                                  (add?Norm(A0):TMV_RealType(T)(0))))) {
            cerr<<"Rank2KUpdate: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"x = "<<TMV_Text(x)<<endl;
            //cerr<<"  "<<x0<<endl;
            cerr<<"y = "<<TMV_Text(y)<<endl;
            //cerr<<"  "<<y0<<endl;
            cerr<<"A = "<<TMV_Text(A)<<endl;
            //cerr<<"  "<<A0<<endl;
            //if (A.isherm()) cerr<<"x*yt = "<<x0*y0.adjoint()<<endl;
            //else cerr<<"x*yT = "<<x0*y0.transpose()<<endl;
            //cerr<<"-> A = "<<A<<endl;
            //cerr<<"A2 = "<<A2<<endl;
            cerr<<"Norm(A-A2) = "<<Norm(A-A2)<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_Rank2K_MMS.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


