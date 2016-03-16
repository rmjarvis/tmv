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


#include "tmv/TMV_TriMatrixArithFunc.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"

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
    // MultMM: U = U * U
    //

    template <class T, class Ta> 
    static void RRMultEqMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, UpperTriMatrixView<T> B)
    {
        TMVAssert(A.isrm());
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
        TMVAssert(B.size() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(B.ct()==NonConj);

        const ptrdiff_t N = B.size();

        if (A.isunit()) {
            for(ptrdiff_t i=0; i<N; ++i) {
                B.row(i,i+1,N) += A.row(i,i+1,N) * B.subTriMatrix(i+1,N);
                if (!B.isunit()) 
                    B.row(i,i,N) *= alpha;
                else TMVAssert(alpha == T(1)); 
            }
        } else {
            TMVAssert(!B.isunit());
            const ptrdiff_t Ads = A.stepi()+A.stepj();
            const Ta* Aii = A.cptr();
            const ptrdiff_t Bds = B.stepi()+B.stepj();
            T* Bii = B.ptr();
            for(ptrdiff_t i=0; i<N; ++i,Aii+=Ads,Bii+=Bds) {
                T aa = A.isconj()?TMV_CONJ(*Aii):*Aii;
                if (alpha != T(1)) aa *= alpha;
                B.row(i,i+1,N) = aa * B.row(i,i+1,N) +
                    alpha * A.row(i,i+1,N) * B.subTriMatrix(i+1,N);
#ifdef TMVFLDEBUG
                TMVAssert(Bii >= B._first);
                TMVAssert(Bii < B._last);
#endif
                *Bii *= aa;
            }
        }
    }

    template <class T, class Ta> 
    static void CRMultEqMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, UpperTriMatrixView<T> B)
    {
        TMVAssert(A.iscm());
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
        TMVAssert(B.size() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(B.ct()==NonConj);

        const ptrdiff_t N = B.size();

        if (A.isunit()) {
            if (B.isunit()) {
                TMVAssert(alpha == T(1));
                for(ptrdiff_t j=1; j<N; ++j) {
                    B.subMatrix(0,j,j+1,N) += A.col(j,0,j) ^ B.row(j,j+1,N);
                    B.col(j,0,j) += A.col(j,0,j);
                }
            } else {
                for(ptrdiff_t j=0; j<N; ++j) 
                    B.subMatrix(0,j,j,N) += A.col(j,0,j) ^ B.row(j,j,N);
                B *= alpha;
            }
        } else {
            TMVAssert(!B.isunit());
            const ptrdiff_t Ads = A.stepi()+A.stepj();
            const Ta* Ajj = A.cptr();
            for(ptrdiff_t j=0; j<N; ++j,Ajj+=Ads) {
                B.subMatrix(0,j,j,N) += A.col(j,0,j) ^ B.row(j,j,N);
                B.row(j,j,N) *= A.isconj()?TMV_CONJ(*Ajj):*Ajj;
            }
            B *= alpha;
        } 
    }

    template <class T, class Ta> 
    static void CMultEqMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, UpperTriMatrixView<T> B)
    {
        TMVAssert(B.iscm());
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
        TMVAssert(B.size() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(B.ct()==NonConj);
#ifdef XDEBUG
        Matrix<Ta> A0 = A;
        Matrix<T> B0 = B;
        Matrix<T> B2 = alpha*A0*B0;
        //cout<<"CMultEqMM Upper\n";
        //cout<<"alpha = "<<alpha<<endl;
        //cout<<"A = "<<TMV_Text(A)<<" = "<<A<<endl;
        //cout<<"B = "<<TMV_Text(B)<<" = "<<B<<endl;
        //cout<<"alpha*A*B = "<<B2<<std::endl;
#endif

        const ptrdiff_t N = B.size();

        if (B.isunit()) {
            // Then alpha = 1 and A.isunit
            if (SameStorage(A,B)) {
                Vector<T> Acolj(N-1);
                for(ptrdiff_t j=N-1;j>0;--j) {
                    Acolj.subVector(0,j) = A.col(j,0,j);
                    B.col(j,0,j) = A.subTriMatrix(0,j) * B.col(j,0,j);
                    B.col(j,0,j) += Acolj.subVector(0,j);
                }
            } else {
                for(ptrdiff_t j=N-1;j>0;--j) {
                    B.col(j,0,j) = A.subTriMatrix(0,j) * B.col(j,0,j);
                    B.col(j,0,j) += A.col(j,0,j);
                }
            }
        } else {
            if (SameStorage(A,B)) {
                Vector<T> Btemp(N);
                for(ptrdiff_t jj=N,j=jj-1;jj>0;--jj,--j) { // jj = j+1
                    Btemp.subVector(0,jj) = 
                        alpha * A.subTriMatrix(0,jj) * B.col(j,0,jj);
                    B.col(j,0,jj) = Btemp.subVector(0,jj);
                }
            } else {
                for(ptrdiff_t j=0;j<N;++j) 
                    B.col(j,0,j+1) = alpha * A.subTriMatrix(0,j+1) * B.col(j,0,j+1);
            }
        }
#ifdef XDEBUG
        if (!(Norm(B-B2) <= 0.001*(Norm(A)*Norm(B)))) {
            cerr<<"CMultEqMM alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B.cptr()<<"  "<<B0<<endl;
            cerr<<"--> B = "<<B<<endl;
            cerr<<"B2 = "<<B2<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta> 
    static void MultEqMM(T alpha,
              const GenUpperTriMatrix<Ta>& A, UpperTriMatrixView<T> B)
        // B = alpha * A * B
    {
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
        TMVAssert(alpha != T(0));
        TMVAssert(B.ct()==NonConj);
        TMVAssert(B.size() > 0);
#ifdef XDEBUG
        //cout<<"MultEqMM: "<<alpha<<"  "<<TMV_Text(A)<<"  "<<A<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        Matrix<Ta> A0 = A;
        Matrix<T> B0 = B;
        Matrix<T> B2 = alpha*A0*B0;
#endif

        const ptrdiff_t nb = TRI_MM_BLOCKSIZE;
        const ptrdiff_t N = A.size();
        const bool samestorage = ( 
            SameStorage(A,B) &&
            ((A.stepi()>A.stepj()) == (B.stepi()>B.stepj()) )  );

        if (N <= TRI_MM_BLOCKSIZE2) {
            if (A.isrm() && B.isrm()) {
                RRMultEqMM(alpha,A,B);
            } else if (A.iscm() && B.isrm()) {
                CRMultEqMM(alpha,A,B);
            } else if (B.iscm()) {
                CMultEqMM(alpha,A,B);
            } else {
                if (B.isunit()) {
                    UpperTriMatrix<T,UnitDiag|ColMajor> BB = B;
                    if (!(A.isrm() || A.iscm())) {
                        if (A.isunit()) {
                            UpperTriMatrix<T,UnitDiag|ColMajor> AA = A;
                            CMultEqMM(alpha,AA,BB.view());
                        } else {
                            UpperTriMatrix<T,NonUnitDiag|ColMajor> AA = A;
                            CMultEqMM(alpha,AA,BB.view());
                        }
                    } else {
                        CMultEqMM(alpha,A,BB.view());
                    }
                    B = BB;
                } else {
                    UpperTriMatrix<T,NonUnitDiag|ColMajor> BB = B;
                    if (!(A.isrm() || A.iscm())) {
                        if (A.isunit()) {
                            UpperTriMatrix<T,UnitDiag|ColMajor> AA = A;
                            CMultEqMM(alpha,AA,BB.view());
                        } else {
                            UpperTriMatrix<T,NonUnitDiag|ColMajor> AA = A;
                            CMultEqMM(alpha,AA,BB.view());
                        }
                    } else {
                        CMultEqMM(alpha,A,BB.view());
                    }
                    B = BB;
                }
            }
        } else {
            ptrdiff_t k = N/2;
            if (k > nb) k = samestorage ? nb :  k/nb*nb;

            // [ A00 A01 ] [ B00 B01 ] = [ A00 B00   A00 B01 + A01 B11 ]
            // [  0  A11 ] [  0  B11 ]   [    0           A11 B11      ]

            ConstUpperTriMatrixView<Ta> A00 = A.subTriMatrix(0,k);
            ConstMatrixView<Ta> A01 = A.subMatrix(0,k,k,N);
            ConstUpperTriMatrixView<Ta> A11 = A.subTriMatrix(k,N);
            UpperTriMatrixView<T> B00 = B.subTriMatrix(0,k);
            MatrixView<T> B01 = B.subMatrix(0,k,k,N);
            UpperTriMatrixView<T> B11 = B.subTriMatrix(k,N);

            if (samestorage) {
                Matrix<T> xB01 = alpha * A00 * B01;
                xB01 += alpha * A01 * B11;
                B01 = xB01;
            } else {
                B01 = alpha * A00 * B01;
                B01 += alpha * A01 * B11;
            }
            MultEqMM(alpha,A00,B00);
            MultEqMM(alpha,A11,B11);
        }
#ifdef XDEBUG
        if (!(Norm(B-B2) <= 0.001*(Norm(A)*Norm(B)))) {
            cerr<<"MultEqMM alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B.cptr()<<"  "<<B0<<endl;
            cerr<<"--> B = "<<B<<endl;
            cerr<<"B2 = "<<B2<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta> 
    static void RRMultEqMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, LowerTriMatrixView<T> B)
    {
        //cout<<"RRMultEqMM Lower\n";
        TMVAssert(A.isrm());
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
        TMVAssert(B.size() > 0);
        TMVAssert(alpha != T(0));

        const ptrdiff_t N = B.size();

        if (A.isunit()) {
            for(ptrdiff_t i=N-1; i>=0; --i) {
                B.row(i,0,i) += A.row(i,0,i) * B.subTriMatrix(0,i);
                B.row(i,0,i) *= alpha;
            }
            if (!B.isunit()) 
                B.diag() *= alpha;
        } else {
            TMVAssert(!B.isunit());
            const ptrdiff_t Ads = A.stepi()+A.stepj();
            const Ta* Aii = A.cptr()+(N-1)*Ads;
            const ptrdiff_t Bds = B.stepi()+B.stepj();
            T* Bii = B.ptr()+(N-1)*Bds;
            for(ptrdiff_t i=N-1; i>=0; --i,Aii-=Ads,Bii-=Bds) {
                T aa = A.isconj()?TMV_CONJ(*Aii):*Aii;
                if (alpha != T(1)) aa *= alpha;
                B.row(i,0,i) = aa * B.row(i,0,i) +
                    alpha * A.row(i,0,i) * B.subTriMatrix(0,i);
#ifdef TMVFLDEBUG
                TMVAssert(Bii >= B._first);
                TMVAssert(Bii < B._last);
#endif
                *Bii *= aa;
            }
        }
    }

    template <class T, class Ta> 
    static void CRMultEqMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, LowerTriMatrixView<T> B)
    {
        TMVAssert(A.iscm());
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
        TMVAssert(B.size() > 0);
        TMVAssert(alpha != T(0));

        const ptrdiff_t N = B.size();

        if (A.isunit()) {
            if (B.isunit()) {
                for(ptrdiff_t j=N-1,jj=N; jj>0; --j,--jj) {
                    B.subMatrix(jj,N,0,j) += A.col(j,jj,N) ^ B.row(j,0,j);
                    B.col(j,jj,N) += A.col(j,jj,N);
                }
            } else {
                for(ptrdiff_t j=N-1,jj=N; jj>0; --j,--jj) 
                    B.subMatrix(jj,N,0,jj) += A.col(j,jj,N) ^ B.row(j,0,jj);
                B *= alpha;
            }
        } else {
            const ptrdiff_t Ads = A.stepi()+A.stepj();
            const Ta* Ajj = A.cptr()+(N-1)*Ads;
            TMVAssert(!B.isunit());
            for(ptrdiff_t j=N-1,jj=N; jj>0; --j,--jj,Ajj-=Ads) {
                B.subMatrix(jj,N,0,jj) += A.col(j,jj,N) ^ B.row(j,0,jj);
                B.row(j,0,jj) *= A.isconj()?TMV_CONJ(*Ajj):*Ajj;
            }
            B *= alpha;
        } 
    }

    template <class T, class Ta> 
    static void CMultEqMM(
        T alpha, const GenLowerTriMatrix<Ta>& A, LowerTriMatrixView<T> B)
    {
        TMVAssert(B.iscm());
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
        TMVAssert(B.size() > 0);
        TMVAssert(alpha != T(0));

        const ptrdiff_t N = B.size();
        if (B.isunit()) {
            // Then alpha = 1 and A.isunit
            for(ptrdiff_t j=0;j<N-1;++j) {
                B.col(j,j+1,N) = A.subTriMatrix(j+1,N) * B.col(j,j+1,N);
                B.col(j,j+1,N) += A.col(j,j+1,N);
            }
        } else {
            for(ptrdiff_t j=0;j<N;++j) 
                B.col(j,j,N) = alpha * A.subTriMatrix(j,N) * B.col(j,j,N);
        }
    }

    template <class T, class Ta> 
    static void MultEqMM(T alpha,
              const GenLowerTriMatrix<Ta>& A, LowerTriMatrixView<T> B)
        // B = alpha * A * B
    {
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
        TMVAssert(alpha != T(0));
        TMVAssert(B.ct()==NonConj);
        TMVAssert(B.size() > 0);
#ifdef XDEBUG
        //cout<<"MultEqMM: "<<alpha<<"  "<<TMV_Text(A)<<"  "<<A<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        Matrix<Ta> A0 = A;
        Matrix<T> B0 = B;
        Matrix<T> B2 = alpha*A0*B0;
#endif

        const ptrdiff_t nb = TRI_MM_BLOCKSIZE;
        const ptrdiff_t N = A.size();
        const bool samestorage = ( 
            SameStorage(A,B) &&
            ((A.stepi()>A.stepj()) == (B.stepi()>B.stepj()) )  );

        if (N <= TRI_MM_BLOCKSIZE2) {
            if (A.isrm() && B.isrm()) {
                RRMultEqMM(alpha,A,B);
            } else if (A.iscm() && B.isrm()) {
                CRMultEqMM(alpha,A,B);
            } else if (B.iscm()) {
                CMultEqMM(alpha,A,B);
            } else {
                if (B.isunit()) {
                    LowerTriMatrix<T,UnitDiag|ColMajor> BB = B;
                    if (!(A.isrm() || A.iscm())) {
                        if (A.isunit()) {
                            LowerTriMatrix<T,UnitDiag|ColMajor> AA = A;
                            CMultEqMM(alpha,AA,BB.view());
                        } else {
                            LowerTriMatrix<T,NonUnitDiag|ColMajor> AA = A;
                            CMultEqMM(alpha,AA,BB.view());
                        }
                    } else {
                        CMultEqMM(alpha,A,BB.view());
                    }
                    B = BB;
                } else {
                    LowerTriMatrix<T,NonUnitDiag|ColMajor> BB = B;
                    if (!(A.isrm() || A.iscm())) {
                        if (A.isunit()) {
                            LowerTriMatrix<T,UnitDiag|ColMajor> AA = A;
                            CMultEqMM(alpha,AA,BB.view());
                        } else {
                            LowerTriMatrix<T,NonUnitDiag|ColMajor> AA = A;
                            CMultEqMM(alpha,AA,BB.view());
                        }
                    } else {
                        CMultEqMM(alpha,A,BB.view());
                    }
                    B = BB;
                }
            }
        } else {
            ptrdiff_t k = N/2;
            if (k > nb) k = samestorage ? nb : k/nb*nb;

            // [ A00  0  ] [ B00  0  ] = [      A00 B00           0    ]
            // [ A10 A11 ] [ B10 B11 ]   [ A10 B00 + A11 B10   A11 B11 ]

            ConstLowerTriMatrixView<Ta> A00 = A.subTriMatrix(0,k);
            ConstMatrixView<Ta> A10 = A.subMatrix(k,N,0,k);
            ConstLowerTriMatrixView<Ta> A11 = A.subTriMatrix(k,N);
            LowerTriMatrixView<T> B00 = B.subTriMatrix(0,k);
            MatrixView<T> B10 = B.subMatrix(k,N,0,k);
            LowerTriMatrixView<T> B11 = B.subTriMatrix(k,N);

            if (samestorage) {
                Matrix<T> xB10 = alpha * A11 * B10;
                xB10 += alpha * A10 * B00;
                B10 = xB10;
            } else {
                B10 = alpha * A11 * B10;
                B10 += alpha * A10 * B00;
            }
            MultEqMM(alpha,A00,B00);
            MultEqMM(alpha,A11,B11);
        }
#ifdef XDEBUG
        if (!(Norm(B-B2) <= 0.001*(Norm(A)*Norm(B)))) {
            cerr<<"MultEqMM alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B.cptr()<<"  "<<B0<<endl;
            cerr<<"--> B = "<<B<<endl;
            cerr<<"B2 = "<<B2<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta, class Tb> 
    static void RRAddMultMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, 
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<T> C)
    {
        TMVAssert(A.isrm());
        TMVAssert(B.isrm());
        TMVAssert(C.isrm());
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() == C.size());
        TMVAssert(!C.isunit());
        TMVAssert(C.size() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(C.ct()==NonConj);

        const ptrdiff_t N = C.size();

        if (A.isunit()) {
            if (B.isunit()) {
                const ptrdiff_t Cds = C.stepi()+C.stepj();
                T* Cii = C.ptr();
                for(ptrdiff_t i=0; i<N; ++i,Cii+=Cds) {
                    C.row(i,i+1,N) += 
                        alpha * A.row(i,i+1,N) * B.subTriMatrix(i+1,N);
                    C.row(i,i+1,N) += alpha * B.row(i,i+1,N);
#ifdef TMVFLDEBUG
                    TMVAssert(Cii >= C._first);
                    TMVAssert(Cii < C._last);
#endif
                    *Cii += alpha;
                }
            } else {
                for(ptrdiff_t i=0; i<N; ++i) {
                    C.row(i,i+1,N) += 
                        alpha * A.row(i,i+1,N) * B.subTriMatrix(i+1,N);
                    C.row(i,i,N) += alpha * B.row(i,i,N);
                }
            }
        } else {
            for(ptrdiff_t i=0; i<N; ++i) 
                C.row(i,i,N) += alpha * A.row(i,i,N) * B.subTriMatrix(i,N);
        }
    }

    template <class T, class Ta, class Tb> 
    static void CRAddMultMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
        UpperTriMatrixView<T> C)
    {
        TMVAssert(A.iscm());
        TMVAssert(B.isrm());
        TMVAssert(C.isrm());
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() == C.size());
        TMVAssert(!C.isunit());
        TMVAssert(C.size() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(C.ct()==NonConj);

        const ptrdiff_t N = C.size();

        if (A.isunit()) {
            if (B.isunit()) {
                const ptrdiff_t Cds = C.stepi()+C.stepj();
                T* Cjj = C.ptr();
                for(ptrdiff_t j=0; j<N; ++j,Cjj+=Cds) {
                    C.subMatrix(0,j,j+1,N) += alpha * A.col(j,0,j) ^ B.row(j,j+1,N);
                    C.col(j,0,j) += alpha * A.col(j,0,j);
                    C.row(j,j+1,N) += alpha * B.row(j,j+1,N);
#ifdef TMVFLDEBUG
                    TMVAssert(Cjj >= C._first);
                    TMVAssert(Cjj < C._last);
#endif
                    *Cjj += alpha;
                }
            } else {
                for(ptrdiff_t j=0; j<N; ++j) {
                    C.subMatrix(0,j,j,N) += alpha * A.col(j,0,j) ^ B.row(j,j,N);
                    C.row(j,j,N) += alpha * B.row(j,j,N);
                }
            }
        } else {
            if (B.isunit()) {
                for(ptrdiff_t j=0; j<N; ++j) {
                    C.subMatrix(0,j+1,j+1,N) += alpha * A.col(j,0,j+1)^B.row(j,j+1,N);
                    C.col(j,0,j+1) += alpha * A.col(j,0,j+1);
                }
            } else {
                for(ptrdiff_t j=0; j<N; ++j) 
                    C.subMatrix(0,j+1,j,N) += alpha * A.col(j,0,j+1) ^ B.row(j,j,N);
            }
        } 
    }

    template <class T, class Ta, class Tb> 
    static void CAddMultMM(
        T alpha, const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
        UpperTriMatrixView<T> C)
    {
        TMVAssert(B.iscm());
        TMVAssert(C.iscm());
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() == C.size());
        TMVAssert(!C.isunit());
        TMVAssert(C.size() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(C.ct()==NonConj);

        const ptrdiff_t N = C.size();

        if (B.isunit()) {
            if (A.isunit()) {
                const ptrdiff_t Cds = C.stepi()+C.stepj();
                T* Cjj = C.ptr();
                for(ptrdiff_t j=0;j<N;++j,Cjj+=Cds) {
                    C.col(j,0,j) += alpha * A.subTriMatrix(0,j) * B.col(j,0,j);
                    C.col(j,0,j) += alpha * A.col(j,0,j);
#ifdef TMVFLDEBUG
                    TMVAssert(Cjj >= C._first);
                    TMVAssert(Cjj < C._last);
#endif
                    *Cjj += alpha;
                }
            } else {
                for(ptrdiff_t j=0;j<N;++j) {
                    C.col(j,0,j) += alpha * A.subTriMatrix(0,j) * B.col(j,0,j);
                    C.col(j,0,j+1) += alpha * A.col(j,0,j+1);
                }
            }
        } else {
            for(ptrdiff_t j=0;j<N;++j)
                C.col(j,0,j+1) += alpha * A.subTriMatrix(0,j+1) * B.col(j,0,j+1);
        }
    }

    template <class T, class Ta, class Tb> 
    static void AddMultMM(
        T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<T> C)
        // C += alpha * A * B
    {
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() == C.size());
        TMVAssert(!C.isunit());
        TMVAssert(C.size() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(C.ct()==NonConj);

        const ptrdiff_t nb = TRI_MM_BLOCKSIZE;
        const ptrdiff_t N = A.size();

        if (N <= TRI_MM_BLOCKSIZE2) {
            if (A.isrm() && B.isrm() && C.isrm()) {
                RRAddMultMM(alpha,A,B,C);
            } else if (A.iscm() && B.isrm() && C.isrm()) {
                CRAddMultMM(alpha,A,B,C);
            } else if (!A.isrm() && !A.iscm()) {
                if (A.isunit()) {
                    UpperTriMatrix<T,UnitDiag|ColMajor> AA = A;
                    AddMultMM(alpha,AA,B,C);
                } else {
                    UpperTriMatrix<T,NonUnitDiag|ColMajor> AA = A;
                    AddMultMM(alpha,AA,B,C);
                }
            } else if (B.iscm() && C.iscm()) {
                CAddMultMM(alpha,A,B,C);
            } else if (C.isrm()) { 
                TMVAssert(!B.isrm());
                if (B.isunit()) {
                    UpperTriMatrix<T,UnitDiag|RowMajor> BB = B;
                    AddMultMM(alpha,A,BB,C);
                } else {
                    UpperTriMatrix<T,NonUnitDiag|RowMajor> BB = B;
                    AddMultMM(alpha,A,BB,C);
                }
            } else if (C.iscm()) { 
                TMVAssert(!B.iscm());
                if (B.isunit()) {
                    UpperTriMatrix<T,UnitDiag|ColMajor> BB = B;
                    AddMultMM(alpha,A,BB,C);
                } else {
                    UpperTriMatrix<T,NonUnitDiag|ColMajor> BB = B;
                    AddMultMM(alpha,A,BB,C);
                }
            } else {
                UpperTriMatrix<T,NonUnitDiag|ColMajor> CC = C;
                AddMultMM(alpha,A,B,CC.view());
                C = CC;
            }
        } else {
            ptrdiff_t k = N/2;
            if (k > nb) k = k/nb*nb;

            // [ A00 A01 ] [ B00 B01 ] = [ A00 B00   A00 B01 + A01 B11 ]
            // [  0  A11 ] [  0  B11 ]   [    0           A11 B11      ]

            ConstUpperTriMatrixView<Ta> A00 = A.subTriMatrix(0,k);
            ConstMatrixView<Ta> A01 = A.subMatrix(0,k,k,N);
            ConstUpperTriMatrixView<Ta> A11 = A.subTriMatrix(k,N);
            ConstUpperTriMatrixView<Tb> B00 = B.subTriMatrix(0,k);
            ConstMatrixView<Tb> B01 = B.subMatrix(0,k,k,N);
            ConstUpperTriMatrixView<Tb> B11 = B.subTriMatrix(k,N);
            UpperTriMatrixView<T> C00 = C.subTriMatrix(0,k);
            MatrixView<T> C01 = C.subMatrix(0,k,k,N);
            UpperTriMatrixView<T> C11 = C.subTriMatrix(k,N);

            AddMultMM(alpha,A00,B00,C00);
            C01 += alpha * A00 * B01;
            C01 += alpha * A01 * B11;
            AddMultMM(alpha,A11,B11,C11);
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    static void TempMultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<T> C)
        // C (+)= alpha * A * B
    { 
        if (B.isrm()) {
            if (C.isunit() && B.isunit()) {
                UpperTriMatrix<T,UnitDiag|RowMajor> tempB(B);
                MultEqMM(alpha,A,tempB.view());
                if (add) C += tempB;
                else C = tempB;
            } else {
                UpperTriMatrix<T,NonUnitDiag|RowMajor> tempB(B);
                MultEqMM(alpha,A,tempB.view());
                if (add) C += tempB;
                else C = tempB;
            }
        } else {
            if (C.isunit() && B.isunit()) {
                UpperTriMatrix<T,UnitDiag|ColMajor> tempB(B);
                MultEqMM(alpha,A,tempB.view());
                if (add) C += tempB;
                else C = tempB;
            } else {
                UpperTriMatrix<T,NonUnitDiag|ColMajor> tempB(B);
                MultEqMM(alpha,A,tempB.view());
                if (add) C += tempB;
                else C = tempB;
            }
        }
    }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<T> C)
        // C (+)= alpha * A * B
    { 
#ifdef XDEBUG
        cout<<"MultMM:  alpha = "<<alpha<<endl;
        cout<<"A = "<<A.cptr()<<"  "<<TMV_Text(A)<<"  "<<A<<endl;
        cout<<"B = "<<B.cptr()<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        cout<<"C = "<<C.cptr()<<"  "<<TMV_Text(C);
        if (add) cout<<"  "<<C;
        cout<<endl;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C0 = C;
        Matrix<T> C2 = alpha*A0*B0;
        if (add) C2 += C0;
#endif
        TMVAssert(A.size() == C.size());
        TMVAssert(A.size() == B.size());
        TMVAssert(!C.isunit() || 
                  (A.isunit() && B.isunit() && alpha == T(1) && !add) );

        if (C.size() > 0) {
            if (C.isconj()) {
                MultMM<add>(TMV_CONJ(alpha),A.conjugate(),B.conjugate(),
                            C.conjugate());
            } else if (alpha==T(0)) {
                if (!add) C.setZero();
            } else if (add) {
                if (SameStorage(A,C) || SameStorage(B,C)) {
                    TempMultMM<add>(alpha,A,B,C);
                }
                else {
                    AddMultMM(alpha,A,B,C);
                }
            } else {
                if (SameStorage(A,C)) {
                    if (SameStorage(B,C)) {
                        if (C.isSameAs(B)) MultEqMM(alpha,A,C);
                        else TempMultMM<add>(alpha,A,B,C);
                    } else {
                        C = A;
                        MultEqMM(alpha,B.transpose(),C.transpose());
                    }
                } else {
                    C = B;
                    MultEqMM(alpha,A,C);
                }
            }
        }
#ifdef XDEBUG
        cout<<"Done: C = "<<C<<endl;
        if (!(Norm(C-C2) <= 
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                     (add?Norm(C0):TMV_RealType(T)(0))))) {
            cerr<<"MultMM alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B.cptr()<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C.cptr()<<C0<<endl;
            cerr<<"C => "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            cerr<<"Norm(C2-C) = "<<Norm(C2-C)<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_MultUU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


