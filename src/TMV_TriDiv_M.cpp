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
#include "tmv/TMV_TriDiv.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef NOTHROW
#include <iostream>
#endif

#ifdef XDEBUG
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define TRI_DIV_BLOCKSIZE TMV_BLOCKSIZE
#define TRI_DIV_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define TRI_DIV_BLOCKSIZE 64
#define TRI_DIV_BLOCKSIZE2 32
#endif

    //
    // TriLDivEq M
    //

    template <class T, class Ta> 
    static void RowTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, MatrixView<T> B)
    {
        // Solve A X = B  where A is an upper triangle matrix
        const ptrdiff_t N = A.size();
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.colsize()>0);
        TMVAssert(B.rowsize()>0);
        TMVAssert(B.ct() == NonConj);

        if (A.isunit()) {
            for(ptrdiff_t i=N-1; i>=0; --i) 
                B.row(i) -= A.row(i,i+1,N) * B.rowRange(i+1,N);
        } else {
            const ptrdiff_t Ads = A.stepi() + A.stepj();
            const Ta* Aii = A.cptr() + (N-1) * Ads;
            for(ptrdiff_t i=N-1; i>=0; --i,Aii-=Ads) {
                B.row(i) -= A.row(i,i+1,N) * B.rowRange(i+1,N);
                if (*Aii==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular UpperTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularUpperTriMatrix<Ta>(A);
#endif
                }
                B.row(i) /= (A.isconj() ? TMV_CONJ(*Aii) : *Aii);
            }
        }
    }

    template <class T, class Ta> 
    static void ColTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, MatrixView<T> B)
    {
        // Solve A X = B  where A is an upper triangle matrix
        const ptrdiff_t N = A.size();
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.colsize()>0);
        TMVAssert(B.rowsize()>0);
        TMVAssert(B.ct() == NonConj);

        if (A.isunit()) {
            for(ptrdiff_t j=N-1; j>0; --j) 
                B.rowRange(0,j) -= A.col(j,0,j) ^ B.row(j);
        } else {
            const ptrdiff_t Ads = A.stepi()+A.stepj();
            const Ta* Ajj = A.cptr() + (N-1)*Ads;
            for(ptrdiff_t j=N-1; j>=0; --j,Ajj-=Ads) {
                if (*Ajj==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular UpperTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularUpperTriMatrix<Ta>(A);
#endif
                }
                B.row(j) /= (A.isconj() ? TMV_CONJ(*Ajj) : *Ajj);
                B.rowRange(0,j) -= A.col(j,0,j) ^ B.row(j);
            }
        } 
    }

    template <class T, class Ta> 
    static void RowTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, MatrixView<T> B)
    {
        // Solve A X = B  where A is a lower triangle matrix
        const ptrdiff_t N = A.size();
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.colsize()>0);
        TMVAssert(B.rowsize()>0);
        TMVAssert(B.ct() == NonConj);

        if (A.isunit()) {
            for(ptrdiff_t i=0;i<N;++i) 
                B.row(i) -= A.row(i,0,i) * B.rowRange(0,i);
        } else {
            const ptrdiff_t Ads = A.stepi()+A.stepj();
            const Ta* Aii = A.cptr();
            for(ptrdiff_t i=0;i<N;++i,Aii+=Ads) {
                B.row(i) -= A.row(i,0,i) * B.rowRange(0,i);
                if (*Aii==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular LowerTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularLowerTriMatrix<Ta>(A);
#endif
                }
                B.row(i) /= (A.isconj() ? TMV_CONJ(*Aii) : *Aii);
            }
        }
    }

    template <class T, class Ta> 
    static void ColTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, MatrixView<T> B)
    {
        // Solve A X = B  where A is a lower triangle matrix
        const ptrdiff_t N = A.size();
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.colsize()>0);
        TMVAssert(B.rowsize()>0);
        TMVAssert(B.ct() == NonConj);

        if (A.isunit()) {
            for(ptrdiff_t j=0;j<N;++j) 
                B.rowRange(j+1,N) -= A.col(j,j+1,N) ^ B.row(j);
        } else {
            const ptrdiff_t Ads = A.stepi()+A.stepj();
            const Ta* Ajj = A.cptr();
            for(ptrdiff_t j=0;j<N;++j,Ajj+=Ads) {
                if (*Ajj==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular LowerTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularLowerTriMatrix<Ta>(A);
#endif
                }
                B.row(j) /= (A.isconj() ? TMV_CONJ(*Ajj) : *Ajj);
                B.rowRange(j+1,N) -= A.col(j,j+1,N) ^ B.row(j);
            }
        } 
    }

    template <class T, class Ta> 
    static void NonBlasTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, MatrixView<T> B)
    // B = A^-1 * B
    // where A is a triangle matrix
    {
        //cout<<"Upper Tri LDivEq Matrix\n";
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.colsize()>0);
        TMVAssert(B.rowsize()>0);
        TMVAssert(B.ct() == NonConj);

        const ptrdiff_t nb = TRI_DIV_BLOCKSIZE;
        const ptrdiff_t N = A.size();

        if (N <= TRI_DIV_BLOCKSIZE2) {
            if (B.isrm()) {
                if (A.isrm()) RowTriLDivEq(A,B);
                else ColTriLDivEq(A,B);
            } else {
                const ptrdiff_t K = B.rowsize();
                for(ptrdiff_t j=0;j<K;++j) TriLDivEq(A,B.col(j));
            }
        } else {
            ptrdiff_t k = N/2;
            if (k > nb) k = k/nb*nb;

            ConstUpperTriMatrixView<Ta> A00 = A.subTriMatrix(0,k);
            ConstUpperTriMatrixView<Ta> A11 = A.subTriMatrix(k,N);
            ConstMatrixView<Ta> A01 = A.subMatrix(0,k,k,N);
            MatrixView<T> B0 = B.rowRange(0,k);
            MatrixView<T> B1 = B.rowRange(k,N);

            NonBlasTriLDivEq(A11,B1);
            B0 -= A01*B1;
            NonBlasTriLDivEq(A00,B0);
        }
    }

    template <class T, class Ta> 
    static void NonBlasTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, MatrixView<T> B)
    // B = A^-1 * B
    // where A is a triangle matrix
    {
        //cout<<"Lower Tri LDivEq Matrix\n";
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.colsize()>0);
        TMVAssert(B.rowsize()>0);
        TMVAssert(B.ct() == NonConj);

        const ptrdiff_t nb = TRI_DIV_BLOCKSIZE;
        const ptrdiff_t N = A.size();

        if (N <= TRI_DIV_BLOCKSIZE2) {
            if (B.isrm()) {
                if (A.isrm()) RowTriLDivEq(A,B);
                else ColTriLDivEq(A,B);
            } else {
                const ptrdiff_t K = B.rowsize();
                for(ptrdiff_t j=0;j<K;++j) TriLDivEq(A,B.col(j));
            }
        } else {
            ptrdiff_t k = N/2;
            if (k > nb) k = k/nb*nb;

            ConstLowerTriMatrixView<Ta> A00 = A.subTriMatrix(0,k);
            ConstLowerTriMatrixView<Ta> A11 = A.subTriMatrix(k,N);
            ConstMatrixView<Ta> A10 = A.subMatrix(k,N,0,k);
            MatrixView<T> B0 = B.rowRange(0,k);
            MatrixView<T> B1 = B.rowRange(k,N);

            NonBlasTriLDivEq(A00,B0);
            B1 -= A10*B0;
            return NonBlasTriLDivEq(A11,B1);
        }
    }

#ifdef BLAS
    template <class T, class Ta> 
    static inline void BlasTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, MatrixView<T> B)
    { NonBlasTriLDivEq(A,B); }
    template <class T, class Ta> 
    static inline void BlasTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, MatrixView<T> B)
    { NonBlasTriLDivEq(A,B); }
#ifdef INST_DOUBLE
    template <> 
    void BlasTriLDivEq(
        const GenUpperTriMatrix<double>& A, MatrixView<double> B)
    {
        int m=BlasIsCM(B)?B.colsize():B.rowsize();
        int n=BlasIsCM(B)?B.rowsize():B.colsize();
        double alpha = 1.;
        int lda = BlasIsCM(A)?A.stepj():A.stepi();
        int ldb = BlasIsCM(B)?B.stepj():B.stepi();

        BLASNAME(dtrsm) (
            BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
            BlasIsCM(A)?BLASCH_UP:BLASCH_LO, 
            BlasIsCM(A)==BlasIsCM(B)?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, 
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasTriLDivEq(
        const GenLowerTriMatrix<double>& A, MatrixView<double> B)
    {
        int m=BlasIsCM(B)?B.colsize():B.rowsize();
        int n=BlasIsCM(B)?B.rowsize():B.colsize();
        double alpha = 1.;
        int lda = BlasIsCM(A)?A.stepj():A.stepi();
        int ldb = BlasIsCM(B)?B.stepj():B.stepi();

        BLASNAME(dtrsm) (
            BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
            BlasIsCM(A)?BLASCH_LO:BLASCH_UP, 
            BlasIsCM(A)==BlasIsCM(B)?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, 
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasTriLDivEq(
        const GenUpperTriMatrix<std::complex<double> >& A,
        MatrixView<std::complex<double> > B)
    {
        int m=BlasIsCM(B)?B.colsize():B.rowsize();
        int n=BlasIsCM(B)?B.rowsize():B.colsize();
        std::complex<double> alpha = 1.;
        int lda = BlasIsCM(A)?A.stepj():A.stepi();
        int ldb = BlasIsCM(B)?B.stepj():B.stepi();
        if (BlasIsCM(A)==BlasIsCM(B) && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ztrsm) (
                BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
                BlasIsCM(A)?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ztrsm) (
                BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
                BlasIsCM(A)?BLASCH_UP:BLASCH_LO, 
                BlasIsCM(A)==BlasIsCM(B)?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasTriLDivEq(
        const GenLowerTriMatrix<std::complex<double> >& A,
        MatrixView<std::complex<double> > B)
    {
        int m=BlasIsCM(B)?B.colsize():B.rowsize();
        int n=BlasIsCM(B)?B.rowsize():B.colsize();
        std::complex<double> alpha = 1.;
        int lda = BlasIsCM(A)?A.stepj():A.stepi();
        int ldb = BlasIsCM(B)?B.stepj():B.stepi();
        if (BlasIsCM(A)==BlasIsCM(B) && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ztrsm) (
                BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
                BlasIsCM(A)?BLASCH_LO:BLASCH_UP, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ztrsm) (
                BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
                BlasIsCM(A)?BLASCH_LO:BLASCH_UP, 
                BlasIsCM(A)==BlasIsCM(B)?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasTriLDivEq(
        const GenUpperTriMatrix<double>& A,
        MatrixView<std::complex<double> > B)
    {
        if (BlasIsRM(B)) {
            int m=2*B.rowsize();
            int n=B.colsize();
            double alpha = 1.;
            int lda = A.isrm()?A.stepi():A.stepj();
            int ldb = 2*B.stepi();

            BLASNAME(dtrsm) (
                BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
                BlasIsCM(A)?BLASCH_UP:BLASCH_LO, 
                BlasIsCM(A)==BlasIsCM(B)?BLASCH_NT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP((double*)B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        } else {
            Matrix<double,ColMajor> B1 = B.realPart();
            BlasTriLDivEq(A,B1.view());
            B.realPart() = B1;
            B1 = B.imagPart();
            BlasTriLDivEq(A,B1.view());
            B.imagPart() = B1;
        }
    }
    template <> 
    void BlasTriLDivEq(
        const GenLowerTriMatrix<double>& A,
        MatrixView<std::complex<double> > B)
    {
        if (BlasIsRM(B)) {
            int m=2*B.rowsize();
            int n=B.colsize();
            double alpha = 1.;
            int lda = BlasIsCM(A)?A.stepj():A.stepi();
            int ldb = 2*B.stepi();

            BLASNAME(dtrsm) (
                BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
                BlasIsCM(A)?BLASCH_LO:BLASCH_UP, 
                BlasIsCM(A)==BlasIsCM(B)?BLASCH_NT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP((double*)B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        } else {
            Matrix<double,ColMajor> B1 = B.realPart();
            BlasTriLDivEq(A,B1.view());
            B.realPart() = B1;
            B1 = B.imagPart();
            BlasTriLDivEq(A,B1.view());
            B.imagPart() = B1;
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void BlasTriLDivEq(
        const GenUpperTriMatrix<float>& A, MatrixView<float> B)
    {
        int m=BlasIsCM(B)?B.colsize():B.rowsize();
        int n=BlasIsCM(B)?B.rowsize():B.colsize();
        float alpha = 1.;
        int lda = BlasIsCM(A)?A.stepj():A.stepi();
        int ldb = BlasIsCM(B)?B.stepj():B.stepi();

        BLASNAME(strsm) (
            BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
            BlasIsCM(A)?BLASCH_UP:BLASCH_LO, 
            BlasIsCM(A)==BlasIsCM(B)?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, 
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasTriLDivEq(
        const GenLowerTriMatrix<float>& A, MatrixView<float> B)
    {
        int m=BlasIsCM(B)?B.colsize():B.rowsize();
        int n=BlasIsCM(B)?B.rowsize():B.colsize();
        float alpha = 1.;
        int lda = BlasIsCM(A)?A.stepj():A.stepi();
        int ldb = BlasIsCM(B)?B.stepj():B.stepi();

        BLASNAME(strsm) (
            BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
            BlasIsCM(A)?BLASCH_LO:BLASCH_UP, 
            BlasIsCM(A)==BlasIsCM(B)?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, 
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasTriLDivEq(
        const GenUpperTriMatrix<std::complex<float> >& A,
        MatrixView<std::complex<float> > B)
    {
        int m=BlasIsCM(B)?B.colsize():B.rowsize();
        int n=BlasIsCM(B)?B.rowsize():B.colsize();
        std::complex<float> alpha = 1.;
        int lda = BlasIsCM(A)?A.stepj():A.stepi();
        int ldb = BlasIsCM(B)?B.stepj():B.stepi();
        if (BlasIsCM(A)==BlasIsCM(B) && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ctrsm) (
                BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
                BlasIsCM(A)?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ctrsm) (
                BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
                BlasIsCM(A)?BLASCH_UP:BLASCH_LO, 
                BlasIsCM(A)==BlasIsCM(B)?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasTriLDivEq(
        const GenLowerTriMatrix<std::complex<float> >& A,
        MatrixView<std::complex<float> > B)
    {
        int m=BlasIsCM(B)?B.colsize():B.rowsize();
        int n=BlasIsCM(B)?B.rowsize():B.colsize();
        std::complex<float> alpha = 1.;
        int lda = BlasIsCM(A)?A.stepj():A.stepi();
        int ldb = BlasIsCM(B)?B.stepj():B.stepi();
        if (BlasIsCM(A)==BlasIsCM(B) && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ctrsm) (
                BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
                BlasIsCM(A)?BLASCH_LO:BLASCH_UP, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ctrsm) (
                BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
                BlasIsCM(A)?BLASCH_LO:BLASCH_UP, 
                BlasIsCM(A)==BlasIsCM(B)?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasTriLDivEq(
        const GenUpperTriMatrix<float>& A,
        MatrixView<std::complex<float> > B)
    {
        if (BlasIsRM(B)) {
            int m=2*B.rowsize();
            int n=B.colsize();
            float alpha = 1.;
            int lda = BlasIsCM(A)?A.stepj():A.stepi();
            int ldb = 2*B.stepi();

            BLASNAME(strsm) (
                BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
                BlasIsCM(A)?BLASCH_UP:BLASCH_LO, 
                BlasIsCM(A)==BlasIsCM(B)?BLASCH_NT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP((float*)B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        } else {
            Matrix<float,ColMajor> B1 = B.realPart();
            BlasTriLDivEq(A,B1.view());
            B.realPart() = B1;
            B1 = B.imagPart();
            BlasTriLDivEq(A,B1.view());
            B.imagPart() = B1;
        }
    }
    template <> 
    void BlasTriLDivEq(
        const GenLowerTriMatrix<float>& A,
        MatrixView<std::complex<float> > B)
    {
        if (BlasIsRM(B)) {
            int m=2*B.rowsize();
            int n=B.colsize();
            float alpha = 1.;
            int lda = BlasIsCM(A)?A.stepj():A.stepi();
            int ldb = 2*B.stepi();

            BLASNAME(strsm) (
                BLASCM BlasIsCM(B)?BLASCH_L:BLASCH_R, 
                BlasIsCM(A)?BLASCH_LO:BLASCH_UP, 
                BlasIsCM(A)==BlasIsCM(B)?BLASCH_NT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP((float*)B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        } else {
            Matrix<float,ColMajor> B1 = B.realPart();
            BlasTriLDivEq(A,B1.view());
            B.realPart() = B1;
            B1 = B.imagPart();
            BlasTriLDivEq(A,B1.view());
            B.imagPart() = B1;
        }
    }
#endif // FLOAT
#endif // BLAS

    template <class T, class Ta> 
    void TriLDivEq(
        const GenUpperTriMatrix<Ta>& A, MatrixView<T> B)
    {
#ifdef XDEBUG
        Matrix<Ta> A0(A);
        Matrix<T> B0(B);
#endif

        TMVAssert(A.size() == B.colsize());
        if (B.colsize() > 0 && B.rowsize() > 0) {
            if (B.isconj()) TriLDivEq(A.conjugate(),B.conjugate());
            else if (B.rowsize() == 1) TriLDivEq(A,B.col(0));
            else if (SameStorage(A,B)) {
                if (A.dt() == NonUnitDiag) {
                    if (A.isrm()) {
                        UpperTriMatrix<Ta,NonUnitDiag|RowMajor> tempA = A;
                        TriLDivEq(tempA,B);
                    } else {
                        UpperTriMatrix<Ta,NonUnitDiag|ColMajor> tempA = A;
                        TriLDivEq(tempA,B);
                    }
                } else {
                    if (A.isrm()) {
                        UpperTriMatrix<Ta,UnitDiag|RowMajor> tempA = A;
                        TriLDivEq(tempA,B);
                    } else {
                        UpperTriMatrix<Ta,UnitDiag|ColMajor> tempA = A;
                        TriLDivEq(tempA,B);
                    }
                }
            } else {
#ifdef BLAS
                if (!(BlasIsCM(A) || BlasIsRM(A))) {
                    if (A.isunit()) {
                        UpperTriMatrix<Ta,UnitDiag|ColMajor> AA(A);
                        BlasTriLDivEq(AA,B);
                    } else {
                        UpperTriMatrix<Ta,NonUnitDiag|ColMajor> AA(A);
                        TriLDivEq(AA,B);
                    }
                } else if (!(BlasIsCM(B) || BlasIsRM(B))) {
                    Matrix<T,ColMajor> BB(B);
                    BlasTriLDivEq(A,BB.view());
                    B = BB;
                } else 
                    BlasTriLDivEq(A,B);
#else
                NonBlasTriLDivEq(A,B);
#endif
            }
        }
#ifdef XDEBUG
        Matrix<T> BB = A0*B;
        if (Norm(BB-B0) > 0.001*Norm(A0)*Norm(B0)) {
            cerr<<"TriLDivEq: M/Upper\n";
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"Done: B = "<<B<<endl;
            cerr<<"A*B = "<<BB<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta> 
    void TriLDivEq(const GenLowerTriMatrix<Ta>& A, MatrixView<T> B)
    {
#ifdef XDEBUG
        Matrix<Ta> A0(A);
        Matrix<T> B0(B);
#endif

        TMVAssert(A.size() == B.colsize());
        if (B.colsize() > 0 && B.rowsize() > 0) {
            if (B.isconj()) TriLDivEq(A.conjugate(),B.conjugate());
            else if (B.rowsize() == 1) TriLDivEq(A,B.col(0));
            else if (SameStorage(A,B)) {
                if (A.dt() == NonUnitDiag) {
                    if (A.isrm()) {
                        LowerTriMatrix<Ta,NonUnitDiag|RowMajor> tempA = A;
                        TriLDivEq(tempA,B);
                    } else {
                        LowerTriMatrix<Ta,NonUnitDiag|ColMajor> tempA = A;
                        TriLDivEq(tempA,B);
                    }
                } else {
                    if (A.isrm()) {
                        LowerTriMatrix<Ta,UnitDiag|RowMajor> tempA = A;
                        TriLDivEq(tempA,B);
                    } else {
                        LowerTriMatrix<Ta,UnitDiag|ColMajor> tempA = A;
                        TriLDivEq(tempA,B);
                    }
                }
            } else {
#ifdef BLAS
                if (!(BlasIsCM(A) || BlasIsRM(A))) {
                    if (A.isunit()) {
                        LowerTriMatrix<Ta,UnitDiag|ColMajor> AA(A);
                        TriLDivEq(AA,B);
                    } else {
                        LowerTriMatrix<Ta,NonUnitDiag|ColMajor> AA(A);
                        TriLDivEq(AA,B);
                    }
                } else if (!(BlasIsCM(B) || BlasIsRM(B))) {
                    Matrix<T,ColMajor> BB(B);
                    BlasTriLDivEq(A,BB.view());
                    B = BB;
                } else 
                    BlasTriLDivEq(A,B);
#else
                NonBlasTriLDivEq(A,B);
#endif
            }
        }
#ifdef XDEBUG
        Matrix<T> BB = A0*B;
        if (Norm(BB-B0) > 0.001*Norm(A0)*Norm(B0)) {
            cerr<<"TriLDivEq: M/Lower\n";
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"Done: B = "<<B<<endl;
            cerr<<"A*B = "<<BB<<endl;
            abort();
        }
#endif
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_TriDiv_M.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


