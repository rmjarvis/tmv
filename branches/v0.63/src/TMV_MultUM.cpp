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
#include "tmv/TMV_TriMatrixArithFunc.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_TriMatrixArith.h"

#ifdef XDEBUG
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define TRI_MM_BLOCKSIZE TMV_BLOCKSIZE
#define TRI_MM_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define TRI_MM_BLOCKSIZE 64
#define TRI_MM_BLOCKSIZE2 32
#endif

    //
    // MultEqMM: M = U * M
    //

    template <class T, class Ta> 
    static void RRMultEqMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
    {
        TMVAssert(A.isrm());
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(alpha != T(0));
        TMVAssert(B.rowsize()>0);
        TMVAssert(B.colsize()>0);
        // A,B Same storage is ok
        const int N = B.colsize();

        if (A.isunit()) {
            for(int i=0; i<N; ++i) {
                B.row(i) += A.row(i,i+1,N) * B.rowRange(i+1,N);
            }
        } else {
            const int Ads = A.stepi()+1;
            const Ta* Aii = A.cptr();
            for(int i=0; i<N; ++i,Aii+=Ads) {
                B.row(i) *= A.isconj() ? TMV_CONJ(*Aii) : *Aii;
                B.row(i) += A.row(i,i+1,N) * B.rowRange(i+1,N);
            }
        }
        if (alpha != T(1)) B *= alpha;
    }

    template <class T, class Ta> 
    static void CRMultEqMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
    {
        TMVAssert(A.iscm());
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(alpha != T(0));
        TMVAssert(B.rowsize()>0);
        TMVAssert(B.colsize()>0);
        // A,B Same storage is ok
        const int N = B.colsize();

        if (A.isunit()) {
            for(int j=0; j<N; ++j) 
                B.rowRange(0,j) += A.col(j,0,j) ^ B.row(j);
        } else {
            const int Ads = A.stepi()+A.stepj();
            const Ta* Ajj = A.cptr();
            for(int j=0; j<N; ++j,Ajj+=Ads) {
                B.rowRange(0,j) += A.col(j,0,j) ^ B.row(j);
                B.row(j) *= (A.isconj()?TMV_CONJ(*Ajj):*Ajj);
            }
        }
        if (alpha != T(1)) B *= alpha;
    }

    template <class T, class Ta> 
    static void CMultEqMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize()>0);
        TMVAssert(B.colsize()>0);
        TMVAssert(alpha != T(0));
        TMVAssert(B.ct() == NonConj);
        TMVAssert(!SameStorage(A,B));

        const int N = B.rowsize();
        for(int j=0;j<N;++j) 
            B.col(j) = alpha * A * B.col(j);
    }

    template <class T, class Ta> 
    static void NonBlasMultEqMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
    // B = alpha * A * B
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(B.ct() == NonConj);
#ifdef XDEBUG
        Matrix<T> BB0 = B;
        Matrix<Ta> A0 = A;
        Matrix<T> B2 = alpha * A0 * BB0;
#endif

        const int nb = TRI_MM_BLOCKSIZE;
        const int N = A.size();

        //cout<<"A,TriMMB = "<<N<<","<<TRI_MM_BLOCKSIZE2<<endl;
        if (N <= TRI_MM_BLOCKSIZE2) {
            if (A.isrm() && B.isrm()) {
                RRMultEqMM(alpha,A,B);
            } else if (A.iscm() && B.isrm()) {
                CRMultEqMM(alpha,A,B);
            } else if (!(A.iscm() || A.isrm()) || SameStorage(A,B)) {
                UpperTriMatrix<T,NonUnitDiag,ColMajor> AA = alpha*A;
                NonBlasMultEqMM(T(1),AA,B);
            } else {
                CMultEqMM(alpha,A,B);
            }
        } else {
            int k = N/2;
            if (k > nb) k = k/nb*nb;

            // [ A00 A01 ] [ B0 ] = [ A00 B0 + A01 B1 ]
            // [  0  A11 ] [ B1 ]   [      A11 B1     ]

            ConstUpperTriMatrixView<Ta> A00 = A.subTriMatrix(0,k);
            ConstUpperTriMatrixView<Ta> A11 = A.subTriMatrix(k,N);
            ConstMatrixView<Ta> A01 = A.subMatrix(0,k,k,N);
            MatrixView<T> B0 = B.rowRange(0,k);
            MatrixView<T> B1 = B.rowRange(k,N);

            NonBlasMultEqMM(alpha,A00,B0);
            B0 += alpha * A01 * B1;
            NonBlasMultEqMM(alpha,A11,B1);
        }
#ifdef XDEBUG
        if (Norm(B-B2) > 0.001*(Norm(A0)*Norm(BB0))) {
            cerr<<"NonBlas MultEqMM alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<BB0<<endl;
            cerr<<"B = "<<B<<endl;
            cerr<<"B2 = "<<B2<<endl;
            cerr<<"Norm(B2-B) = "<<Norm(B2-B)<<endl;
            abort();
        }
#endif
    }


    //
    // MultEqMM: M = L * M
    //

    template <class T, class Ta> 
    static void RRMultEqMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
    {
        TMVAssert(A.isrm());
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(alpha != T(0));
        TMVAssert(B.rowsize()>0);
        TMVAssert(B.colsize()>0);
        // A,B Same storage is ok
        const int N = B.colsize();

        if (A.isunit()) {
            for(int i=N-1; i>=0; --i) {
                B.row(i) += A.row(i,0,i) * B.rowRange(0,i);
            }
        } else {
            const int Ads = A.stepi()+A.stepj();
            const Ta* Aii = A.cptr()+(N-1)*Ads;
            for(int i=N-1; i>=0; --i,Aii-=Ads) {
                B.row(i) *= A.isconj() ? TMV_CONJ(*Aii) : *Aii;
                B.row(i) += A.row(i,0,i) * B.rowRange(0,i);
            }
        }
        B *= alpha;
    }

    template <class T, class Ta> 
    static void CRMultEqMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
    {
        TMVAssert(A.iscm());
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize()>0);
        TMVAssert(B.colsize()>0);
        TMVAssert(alpha != T(0));
        TMVAssert(B.ct() == NonConj);
        // A,B Same storage is ok
        const int N = B.colsize();

        if (A.isunit()) {
            for(int jj=N,j=jj-1; jj>0; --j,--jj) { // jj = j+1
                B.rowRange(jj,N) += A.col(j,jj,N) ^ B.row(j);
            }
        } else {
            const int Ads = A.stepi()+A.stepj();
            const Ta* Ajj = A.cptr()+(N-1)*Ads;
            for(int jj=N,j=jj-1; jj>0; --j,--jj,Ajj-=Ads) { // jj = j+1
                B.rowRange(jj,N) += A.col(j,jj,N) ^ B.row(j);
                B.row(j) *= A.isconj()?TMV_CONJ(*Ajj):*Ajj;
            }
        }
        B *= alpha;
    }

    template <class T, class Ta> 
    static void CMultEqMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize()>0);
        TMVAssert(B.colsize()>0);
        TMVAssert(alpha != T(0));
        TMVAssert(B.ct() == NonConj);
        TMVAssert(!SameStorage(A,B));

        const int N = B.rowsize();
        for(int j=0;j<N;++j) 
            B.col(j) = alpha * A * B.col(j);
    }

    template <class T, class Ta> 
    static void NonBlasMultEqMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
    // B = alpha * A * B
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(B.ct() == NonConj);

        const int nb = TRI_MM_BLOCKSIZE;
        const int N = A.size();

        if (N <= TRI_MM_BLOCKSIZE2) {
            if (A.isrm() && B.isrm()) {
                RRMultEqMM(alpha,A,B);
            } else if (A.iscm() && B.isrm()) {
                CRMultEqMM(alpha,A,B);
            } else if (!(A.iscm() || A.isrm()) || SameStorage(A,B)) {
                LowerTriMatrix<T,NonUnitDiag,ColMajor> AA = alpha*A;
                NonBlasMultEqMM(T(1),AA,B);
            } else {
                CMultEqMM(alpha,A,B);
            }
        } else {
            int k = N/2;
            if (k > nb) k = k/nb*nb;

            // [ A00  0  ] [ B0 ] = [      A00 B0     ]
            // [ A10 A11 ] [ B1 ]   [ A10 B0 + A11 B1 ]

            ConstLowerTriMatrixView<Ta> A00 = A.subTriMatrix(0,k);
            ConstLowerTriMatrixView<Ta> A11 = A.subTriMatrix(k,N);
            ConstMatrixView<Ta> A10 = A.subMatrix(k,N,0,k);
            MatrixView<T> B0 = B.rowRange(0,k);
            MatrixView<T> B1 = B.rowRange(k,N);

            NonBlasMultEqMM(alpha,A11,B1);
            B1 += alpha * A10 * B0;
            NonBlasMultEqMM(alpha,A00,B0);
        }
    }

#ifdef BLAS
    template <class T, class Ta> 
    static inline void BlasMultEqMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
    { NonBlasMultEqMM(alpha,A,B); }
    template <class T, class Ta> 
    static inline void BlasMultEqMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
    { NonBlasMultEqMM(alpha,A,B); }
#ifdef INST_DOUBLE
    template <> 
    void BlasMultEqMM(
        double alpha, const GenUpperTriMatrix<double>& A, 
        const MatrixView<double>& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != 0.);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();

        BLASNAME(dtrmm) (
            BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
            A.iscm()?BLASCH_UP:BLASCH_LO, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
            BLASV(alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
            BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasMultEqMM(
        double alpha, const GenLowerTriMatrix<double>& A, 
        const MatrixView<double>& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != 0.);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();

        BLASNAME(dtrmm) (
            BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
            A.iscm()?BLASCH_LO:BLASCH_UP, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
            BLASV(alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
            BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasMultEqMM(
        std::complex<double> alpha,
        const GenUpperTriMatrix<std::complex<double> >& A,
        const MatrixView<std::complex<double> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != 0.);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();

        if (A.iscm()==B.iscm() && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ztrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ztrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasMultEqMM(
        std::complex<double> alpha,
        const GenLowerTriMatrix<std::complex<double> >& A,
        const MatrixView<std::complex<double> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != 0.);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();

        if (A.iscm()==B.iscm() && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ztrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ztrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, 
                A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasMultEqMM(
        std::complex<double> alpha, const GenUpperTriMatrix<double>& A,
        const MatrixView<std::complex<double> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != 0.);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);

        if (B.isrm()) {
            int m = 2*B.rowsize();
            int n = B.colsize();
            int lda = A.isrm()?A.stepi():A.stepj();
            int ldb = 2*B.stepi();
            double xalpha(1);
            BLASNAME(dtrmm) (
                BLASCM BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, 
                A.isrm()?BLASCH_NT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASV(xalpha),BLASP(A.cptr()),BLASV(lda), 
                BLASP((double*)B.ptr()), BLASV(ldb)
                BLAS1 BLAS1 BLAS1 BLAS1);
            B *= alpha;
        } else {
            Matrix<double,ColMajor> B1 = B.real();
            BlasMultEqMM(1.,A,B1.view());
            B.real() = B1;
            B1 = B.imag();
            BlasMultEqMM(1.,A,B1.view());
            B.imag() = B1;
            B *= alpha;
        }
    }
    template <> 
    void BlasMultEqMM(
        std::complex<double> alpha, const GenLowerTriMatrix<double>& A,
        const MatrixView<std::complex<double> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != 0.);
        TMVAssert(B.ct() == NonConj);

        if (B.isrm()) {
            int m=2*B.rowsize();
            int n=B.colsize();
            int lda = A.isrm()?A.stepi():A.stepj();
            int ldb = B.stepi();
            double xalpha(1);
            BLASNAME(dtrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, 
                A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASV(xalpha),BLASP(A.cptr()),BLASV(lda), 
                BLASP((double*)B.ptr()), BLASV(ldb)
                BLAS1 BLAS1 BLAS1 BLAS1);
            B *= alpha;
        } else {
            Matrix<double,ColMajor> B1 = B.real();
            BlasMultEqMM(1.,A,B1.view());
            B.real() = B1;
            B1 = B.imag();
            BlasMultEqMM(1.,A,B1.view());
            B.imag() = B1;
            B *= alpha;
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void BlasMultEqMM(
        float alpha, const GenUpperTriMatrix<float>& A, 
        const MatrixView<float>& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != 0.F);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();

        BLASNAME(strmm) (
            BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
            A.iscm()?BLASCH_UP:BLASCH_LO, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
            BLASV(alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
            BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasMultEqMM(
        float alpha, const GenLowerTriMatrix<float>& A, 
        const MatrixView<float>& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != 0.F);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();

        BLASNAME(strmm) (
            BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
            A.iscm()?BLASCH_LO:BLASCH_UP, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
            BLASV(alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
            BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasMultEqMM(
        std::complex<float> alpha,
        const GenUpperTriMatrix<std::complex<float> >& A,
        const MatrixView<std::complex<float> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != 0.F);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();

        if (A.iscm()==B.iscm() && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ctrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ctrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasMultEqMM(
        std::complex<float> alpha,
        const GenLowerTriMatrix<std::complex<float> >& A,
        const MatrixView<std::complex<float> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != 0.F);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();

        if (A.iscm()==B.iscm() && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ctrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ctrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, 
                A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasMultEqMM(
        std::complex<float> alpha, const GenUpperTriMatrix<float>& A,
        const MatrixView<std::complex<float> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != 0.F);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);

        if (B.isrm()) {
            int m = 2*B.rowsize();
            int n = B.colsize();
            int lda = A.isrm()?A.stepi():A.stepj();
            int ldb = 2*B.stepi();
            float xalpha(1);
            BLASNAME(strmm) (
                BLASCM BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, 
                A.isrm()?BLASCH_NT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASV(xalpha),BLASP(A.cptr()),BLASV(lda), 
                BLASP((float*)B.ptr()), BLASV(ldb)
                BLAS1 BLAS1 BLAS1 BLAS1);
            B *= alpha;
        } else {
            Matrix<float,ColMajor> B1 = B.real();
            BlasMultEqMM(1.F,A,B1.view());
            B.real() = B1;
            B1 = B.imag();
            BlasMultEqMM(1.F,A,B1.view());
            B.imag() = B1;
            B *= alpha;
        }
    }
    template <> 
    void BlasMultEqMM(
        std::complex<float> alpha, const GenLowerTriMatrix<float>& A,
        const MatrixView<std::complex<float> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != 0.F);
        TMVAssert(B.ct() == NonConj);

        if (B.isrm()) {
            int m=2*B.rowsize();
            int n=B.colsize();
            int lda = A.isrm()?A.stepi():A.stepj();
            int ldb = B.stepi();
            float xalpha(1);
            BLASNAME(strmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, 
                A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP((float*)B.ptr()), BLASV(ldb)
                BLAS1 BLAS1 BLAS1 BLAS1);
            B *= alpha;
        } else {
            Matrix<float,ColMajor> B1 = B.real();
            BlasMultEqMM(1.F,A,B1.view());
            B.real() = B1;
            B1 = B.imag();
            BlasMultEqMM(1.F,A,B1.view());
            B.imag() = B1;
            B *= alpha;
        }
    }
#endif
#endif // BLAS

    template <class T, class Ta> 
    static void MultEqMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
    {
#ifdef XDEBUG
        //cout<<"MultEqMM: "<<alpha<<"  "<<TMV_Text(A)<<"  "<<A<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        Matrix<T> B0 = B;
        Matrix<Ta> A0 = A;
        Matrix<T> B2 = alpha * A0 * B0;
#endif
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(B.ct() == NonConj);
#ifdef BLAS
        if ( !((B.isrm() && B.stepi()>0) || (B.iscm() && B.stepj()>0)) ||
             !SameStorage(A,B)) {
            Matrix<T> BB = alpha*B;
            BlasMultEqMM(T(1),A,BB.view());
            B = BB;
        } else if (!((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0))) {
            if (A.isunit()) {
                UpperTriMatrix<Ta,UnitDiag,RowMajor> AA = A;
                BlasMultEqMM(alpha,AA,B);
            } else if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                UpperTriMatrix<Ta,NonUnitDiag,RowMajor> AA = TMV_REAL(alpha)*A;
                BlasMultEqMM(T(1),AA,B);
            } else {
                UpperTriMatrix<T,NonUnitDiag,RowMajor> AA = alpha*A;
                BlasMultEqMM(T(1),AA,B);
            }
        } else {
            BlasMultEqMM(alpha,A,B);
        }
#else
        NonBlasMultEqMM(alpha,A,B);
#endif

#ifdef XDEBUG
        if (Norm(B-B2) > 0.001*(Norm(A0)*Norm(B0))) {
            cerr<<"MultEqMM alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"B = "<<B<<endl;
            cerr<<"B2 = "<<B2<<endl;
            cerr<<"Norm(B2-B) = "<<Norm(B2-B)<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta> 
    static void MultEqMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
    {
#ifdef XDEBUG
        //cout<<"MultEqMM: "<<alpha<<"  "<<TMV_Text(A)<<"  "<<A<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        Matrix<T> B0 = B;
        Matrix<Ta> A0 = A;
        Matrix<T> B2 = alpha * A0 * B0;
#endif
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(B.ct() == NonConj);
#ifdef BLAS
        if ( !((B.isrm() && B.stepi()>0) || (B.iscm() && B.stepj()>0)) ||
             !SameStorage(A,B)) {
            Matrix<T> BB = alpha*B;
            BlasMultEqMM(T(1),A,BB.view());
            B = BB;
        } else if (!((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0))) {
            if (A.isunit()) {
                LowerTriMatrix<Ta,UnitDiag,RowMajor> AA = A;
                BlasMultEqMM(alpha,AA,B);
            } else if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                LowerTriMatrix<Ta,NonUnitDiag,RowMajor> AA = TMV_REAL(alpha)*A;
                BlasMultEqMM(T(1),AA,B);
            } else {
                LowerTriMatrix<T,NonUnitDiag,RowMajor> AA = alpha*A;
                BlasMultEqMM(T(1),AA,B);
            }
        } else {
            BlasMultEqMM(alpha,A,B);
        }
#else
        NonBlasMultEqMM(alpha,A,B);
#endif

#ifdef XDEBUG
        if (Norm(B-B2) > 0.001*(Norm(A0)*Norm(B0))) {
            cerr<<"MultEqMM alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"B = "<<B<<endl;
            cerr<<"B2 = "<<B2<<endl;
            cerr<<"Norm(B2-B) = "<<Norm(B2-B)<<endl;
            abort();
        }
#endif
    }

    //
    // addMultMM: M += U * M
    //

    template <class T, class Ta, class Tb> 
    static void RRAddMultMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        TMVAssert(A.isrm());
        TMVAssert(C.isrm());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(alpha != T(0));
        TMVAssert(C.rowsize()>0);
        TMVAssert(C.colsize()>0);
        TMVAssert(!A.isunit());
        TMVAssert(B.isrm() || !SameStorage(B,C));
        // A,C Same storage is ok

        const int N = C.colsize();

        for(int i=0; i<N; ++i) 
            C.row(i) += alpha * A.row(i,i,N) * B.rowRange(i,N);
    }

    template <class T, class Ta, class Tb> 
    static void CRAddMultMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        TMVAssert(A.iscm());
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(alpha != T(0));
        TMVAssert(C.rowsize()>0);
        TMVAssert(C.colsize()>0);
        TMVAssert(!A.isunit());
        TMVAssert(C.isrm() || !SameStorage(B,C));
        TMVAssert(!SameStorage(A,C));

        const int N = C.colsize();

        if (SameStorage(B,C)) {
            for(int j=0; j<N; ++j) {
                C.rowRange(0,j) += alpha * A.col(j,0,j) ^ B.row(j);
                C.row(j) += alpha * A(j,j) * B.row(j);
            }
        } else {
            for(int j=0; j<N; ++j) 
                C.rowRange(0,j+1) += alpha * A.col(j,0,j+1) ^ B.row(j);
        }
    }

    template <class T, class Ta, class Tb> 
    static void CAddMultMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(alpha != T(0));
        TMVAssert(C.rowsize()>0);
        TMVAssert(C.colsize()>0);
        TMVAssert(!SameStorage(A,B) || (A.stepi() == B.stepi()));
        TMVAssert(!SameStorage(A,C));

        const int N = C.rowsize();
        for(int j=0;j<N;++j) 
            C.col(j) += alpha * A * B.col(j);
    }

    template <class T, class Ta, class Tb> 
    static void addMultMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    // C += alpha * A * B
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(alpha != T(0));
        TMVAssert(C.rowsize() > 0);
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.ct() == NonConj);

#ifdef XDEBUG
        Matrix<Ta> A0 = A;
        Matrix<Tb> BB0 = B;
        Matrix<T> CC0 = C;
        Matrix<T> C2 = CC0 + alpha*A0*BB0;
#endif

        const int nb = TRI_MM_BLOCKSIZE;
        const int N = A.size();

        if (N <= TRI_MM_BLOCKSIZE2) {
            if (A.isunit() && ((A.isrm()&&C.isrm()) || (A.iscm()&&B.isrm())) ) {
                if (SameStorage(B,C) && N > 1) {
                    Matrix<Tb> BB = alpha*B;
                    addMultMM(T(1),A.offDiag(),BB.rowRange(1,N),
                              C.rowRange(0,N-1));
                    C += BB;
                } else {
                    if (N > 1) addMultMM(alpha,A.offDiag(),B.rowRange(1,N),
                                         C.rowRange(0,N-1));
                    C += alpha * B;
                }
            } else if (A.isrm() && C.isrm()) {
                RRAddMultMM(alpha,A,B,C);
            } else if (A.iscm() && B.isrm()) {
                CRAddMultMM(alpha,A,B,C);
            } else if (!(A.isrm() || A.iscm())) {
                UpperTriMatrix<T,NonUnitDiag,ColMajor> AA = A;
                addMultMM(alpha,AA,B,C);
            }
            else CAddMultMM(alpha,A,B,C);
        } else {
            int k = N/2;
            if (k > nb) k = k/nb*nb;

            ConstUpperTriMatrixView<Ta> A00 = A.subTriMatrix(0,k);
            ConstUpperTriMatrixView<Ta> A11 = A.subTriMatrix(k,N);
            ConstMatrixView<Ta> A01 = A.subMatrix(0,k,k,N);
            ConstMatrixView<Tb> B0 = B.rowRange(0,k);
            ConstMatrixView<Tb> B1 = B.rowRange(k,N);
            MatrixView<T> C0 = C.rowRange(0,k);
            MatrixView<T> C1 = C.rowRange(k,N);

            addMultMM(alpha,A00,B0,C0);
            C0 += alpha * A01 * B1;
            addMultMM(alpha,A11,B1,C1);
        }
#ifdef XDEBUG
        if (Norm(C-C2) > 0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(BB0))) {
            cerr<<"addMultMM alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<BB0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<CC0<<endl;
            cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr()<<", Cptr = "<<C.cptr()<<endl;
            cerr<<"C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            cerr<<"Norm(C2-C) = "<<Norm(C2-C)<<endl;
            abort();
        }
#endif
    }

    //
    // addMultMM: M += L * M
    //

    template <class T, class Ta, class Tb> 
    static void RRAddMultMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        TMVAssert(A.isrm());
        TMVAssert(C.isrm());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(alpha != T(0));
        TMVAssert(C.rowsize()>0);
        TMVAssert(C.colsize()>0);

        const int N = C.colsize();

        for(int i=0; i<N; ++i) 
            C.row(i) += alpha * A.row(i,0,i+1) * B.rowRange(0,i+1);
    }

    template <class T, class Ta, class Tb> 
    static void CRAddMultMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        TMVAssert(A.iscm());
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(alpha != T(0));
        TMVAssert(C.rowsize()>0);
        TMVAssert(C.colsize()>0);
        TMVAssert(!A.isunit());

        const int N = C.colsize();
        for(int j=0; j<N; ++j) 
            C.rowRange(j,N) += alpha * A.col(j,j,N) ^ B.row(j);
    }

    template <class T, class Ta, class Tb> 
    static void CAddMultMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(alpha != T(0));
        TMVAssert(C.rowsize()>0);
        TMVAssert(C.colsize()>0);
        TMVAssert(!SameStorage(A,B) || (A.stepi() == B.stepi()));
        TMVAssert(!SameStorage(A,C));

        const int N = C.rowsize();
        for(int j=0;j<N;++j) 
            C.col(j) += alpha * A * B.col(j);
    }

    template <class T, class Ta, class Tb> 
    static void addMultMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    // C += alpha * A * B
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(alpha != T(0));
        TMVAssert(C.rowsize() > 0);
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.ct() == NonConj);

#ifdef XDEBUG
        Matrix<Ta> A0 = A;
        Matrix<Tb> BB0 = B;
        Matrix<T> CC0 = C;
        Matrix<T> C2 = CC0 + alpha*A0*BB0;
#endif

        const int nb = TRI_MM_BLOCKSIZE;
        const int N = A.size();

        if (N <= TRI_MM_BLOCKSIZE2) {
            if (A.isunit() && ((A.isrm()&&C.isrm()) || (A.iscm()&&B.isrm())) ) {
                if (SameStorage(B,C) && N > 1) {
                    Matrix<Tb> BB = alpha*B;
                    addMultMM(T(1),A.offDiag(),BB.rowRange(0,N-1),
                              C.rowRange(1,N));
                    C += BB;
                } else {
                    if (N > 1) addMultMM(alpha,A.offDiag(),B.rowRange(0,N-1),
                                         C.rowRange(1,N));
                    C += alpha * B;
                }
            } else if (A.isrm() && C.isrm()) {
                RRAddMultMM(alpha,A,B,C);
            } else if (A.iscm() && B.isrm()) {
                CRAddMultMM(alpha,A,B,C);
            } else if (!(A.isrm() || A.iscm())) {
                LowerTriMatrix<T,NonUnitDiag,ColMajor> AA = A;
                addMultMM(alpha,AA,B,C);
            } else {
                CAddMultMM(alpha,A,B,C); 
            }
        } else {
            int k = N/2;
            if (k > nb) k = k/nb*nb;

            ConstLowerTriMatrixView<Ta> A00 = A.subTriMatrix(0,k);
            ConstMatrixView<Ta> A10 = A.subMatrix(k,N,0,k);
            ConstLowerTriMatrixView<Ta> A11 = A.subTriMatrix(k,N);
            ConstMatrixView<Tb> B0 = B.rowRange(0,k);
            ConstMatrixView<Tb> B1 = B.rowRange(k,N);
            MatrixView<T> C0 = C.rowRange(0,k);
            MatrixView<T> C1 = C.rowRange(k,N);

            addMultMM(alpha,A11,B1,C1);
            C1 += alpha * A10 * B0;
            addMultMM(alpha,A00,B0,C0);
        }
#ifdef XDEBUG
        if (Norm(C-C2) > 0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(BB0))) {
            cerr<<"addMultMM alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<BB0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<CC0<<endl;
            cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr()<<", Cptr = "<<C.cptr()<<endl;
            cerr<<"C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            cerr<<"Norm(C2-C) = "<<Norm(C2-C)<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    static void BlockTempMultMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    // C (+)= alpha * A * B
    { 
        const int N = C.rowsize();
        for (int j=0;j<N;) {
            int j2 = TMV_MIN(N,j+TRI_MM_BLOCKSIZE);
            if (B.isrm()) {
                Matrix<T,RowMajor> B2 = alpha * B.colRange(j,j2);
                MultEqMM(T(1),A,B2.view());
                if (add) C.colRange(j,j2) += B2;
                else C.colRange(j,j2) = B2;
            } else {
                Matrix<T,ColMajor> B2 = alpha * B.colRange(j,j2);
                MultEqMM(T(1),A,B2.view());
                if (add) C.colRange(j,j2) += B2;
                else C.colRange(j,j2) = B2;
            }
            j = j2;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static void FullTempMultMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    // C (+)= alpha * A * B
    { 
        if (B.isrm()) {
            Matrix<T,RowMajor> B2 = alpha * B;
            MultEqMM(T(1),A,B2.view());
            if (add) C += B2;
            else C = B2;
        } else {
            Matrix<T,ColMajor> B2 = alpha * B;
            MultEqMM(T(1),A,B2.view());
            if (add) C += B2;
            else C = B2;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
        // C (+)= alpha * A * B
    { 
#ifdef XDEBUG
        //cout<<"MultMM: "<<alpha<<"  "<<TMV_Text(A)<<"  "<<A<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        //cout<<TMV_Text(C)<<"  "<<C<<endl;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C0 = C;
        Matrix<T> C2 = alpha*A0*B0;
        if (add) C2 += C0;
#endif
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());

        if (C.colsize() > 0 && C.rowsize() > 0) {
            if (C.isconj()) 
                MultMM<add>(
                    TMV_CONJ(alpha),A.conjugate(),B.conjugate(),C.conjugate());
            else if (alpha==T(0)) {
                if (!add) C.zero();
            }
            else if (SameStorage(A,C)) 
                FullTempMultMM<add>(alpha,A,B,C);
            else if (!add) {
                C = alpha * B;
                MultEqMM(T(1),A,C);
            } 
            else if (SameStorage(B,C)) 
                if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
                    BlockTempMultMM<add>(alpha,A,B,C);
                else
                    FullTempMultMM<add>(alpha,A,B,C);
            else 
                addMultMM(alpha,A,B,C);
        }
#ifdef XDEBUG
        if (Norm(C-C2) > 0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                                (add?Norm(C0):TMV_RealType(T)(0)))) {
            cerr<<"MultMM alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C0<<endl;
            cerr<<"C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            cerr<<"Norm(C2-C) = "<<Norm(C2-C)<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    static void FullTempMultMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    // C (+)= alpha * A * B
    { 
        if (B.isrm()) {
            Matrix<T,RowMajor> B2 = alpha * B;
            MultEqMM(T(1),A,B2.view());
            if (add) C += B2;
            else C = B2;
        } else {
            Matrix<T,ColMajor> B2 = alpha * B;
            MultEqMM(T(1),A,B2.view());
            if (add) C += B2;
            else C = B2;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static void BlockTempMultMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    // C (+)= alpha * A * B
    { 
        const int N = C.rowsize();
        for (int j=0;j<N;) {
            int j2 = TMV_MIN(N,j+TRI_MM_BLOCKSIZE);
            if (B.isrm()) {
                Matrix<T,RowMajor> B2 = alpha * B.colRange(j,j2);
                MultEqMM(T(1),A,B2.view());
                if (add) C.colRange(j,j2) += B2;
                else C.colRange(j,j2) = B2;
            } else {
                Matrix<T,ColMajor> B2 = alpha * B.colRange(j,j2);
                MultEqMM(T(1),A,B2.view());
                if (add) C.colRange(j,j2) += B2;
                else C.colRange(j,j2) = B2;
            }
            j = j2;
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    // C (+)= alpha * A * B
    { 
#ifdef XDEBUG
        //cout<<"MultMM: "<<alpha<<"  "<<TMV_Text(A)<<"  "<<A<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        //cout<<TMV_Text(C)<<"  "<<C<<endl;
        Matrix<T> C0 = C;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C2 = alpha*A0*B0;
        if (add) C2 += C0;
#endif
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());

        if (C.colsize() > 0 && C.rowsize() > 0) {
            if (C.isconj()) 
                MultMM<add>(
                    TMV_CONJ(alpha),A.conjugate(),B.conjugate(),C.conjugate());
            else if (alpha==T(0)) {
                if (!add) C.zero();
            }
            else if (SameStorage(A,C)) 
                FullTempMultMM<add>(alpha,A,B,C);
            else if (!add) {
                C = alpha * B;
                MultEqMM(T(1),A,C);
            } 
            else if (SameStorage(B,C)) 
                if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
                    BlockTempMultMM<add>(alpha,A,B,C);
                else
                    FullTempMultMM<add>(alpha,A,B,C);
            else 
                addMultMM(alpha,A,B,C);
        }
#ifdef XDEBUG
        if (Norm(C-C2) > 0.001*(
                TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                (add?Norm(C0):TMV_RealType(T)(0)))) {
            cerr<<"MultMM alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C0<<endl;
            cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
            cerr<<", Cptr = "<<C.cptr()<<endl;
            cerr<<"C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            cerr<<"Norm(C2-C) = "<<Norm(C2-C)<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_MultUM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


