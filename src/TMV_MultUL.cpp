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
    // MultMM: M = U * L
    //

    template <bool add, class T, class Ta, class Tb> 
    static void ColMultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, MatrixView<T> C)
    {
#ifdef XDEBUG
        //cout<<"MultMM: "<<alpha<<"  "<<TMV_Text(A)<<"  "<<A<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        //cout<<TMV_Text(C)<<"  "<<C<<endl;
        Matrix<T> C2 = C;
        Matrix<T> C0 = C;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        if (add) C2 += alpha*A0*B0;
        else C2 = alpha*A0*B0;
        //cout<<"ColMultMM UL\n";
        //cout<<"A = "<<A<<endl;
        //cout<<"B = "<<B<<endl;
        //cout<<"C = "<<C<<endl;
        //cout<<"alpha = "<<alpha<<endl;
        //cout<<"Correct result = "<<C2<<endl;
#endif
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.size() == C.rowsize());
        TMVAssert(A.size() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(C.ct() == NonConj);
        TMVAssert(!C.isrm() || C.iscm());
        TMVAssert(!C.isconj());
        const ptrdiff_t N = A.size();

        if (SameStorage(A,C) && A.stepj() == C.stepi()) {
            // Then need temporary (see below)
            if (A.isrm()) {
                UpperTriMatrix<Ta,NonUnitDiag|RowMajor> A2 = A;
                ColMultMM<add>(alpha,A2,B,C);
            }
            else {
                UpperTriMatrix<Ta,NonUnitDiag|ColMajor> A2 = A;
                ColMultMM<add>(alpha,A2,B,C);
            }
        } else {
            if (A.isunit()) {
                if (B.isunit()) {
                    T* Cjj = C.ptr();
                    const ptrdiff_t Cds = C.stepi()+C.stepj();
                    for(ptrdiff_t j=0,jj=1;j<N;++j,++jj,Cjj+=Cds) {
                        // C.col(j) (+=) alpha*A*B.col(j)
                        //
                        // C.col(j) (+=) alpha*A.colRange(j,N)*B.col(j,j,N)
                        //
                        // C(j,j) (+=) alpha*A.row(j,j,N) * B.col(j,j,N)
                        // C.col(j,0,j) (+=) 
                        //     alpha*A.subMatrix(0,j,j,N)*B.col(j,j,N)
                        // C.col(j,j+1,N) (+=) 
                        //     alpha*A.subTriMatrix(j,N)*B.col(j,j+1,N)
                        //
                        // Requirements on storage: 
                        //   B can be stored in either triangle
                        //   A cannot be stored in C's lower triangle

                        T newcjj = A.row(j,jj,N)*B.col(j,jj,N) + T(1);
                        if (alpha != T(1)) newcjj *= alpha;
                        if (add) newcjj += *Cjj;

                        if (add) C.col(j,0,j) += alpha * A.col(j,0,j);
                        else C.col(j,0,j) = alpha * A.col(j,0,j);

                        C.col(j,0,j) += 
                            alpha * A.subMatrix(0,j,jj,N) * B.col(j,jj,N);

                        MultMV<add>(alpha,A.subTriMatrix(jj,N),B.col(j,jj,N),
                                    C.col(j,jj,N));

#ifdef TMVFLDEBUG
                        TMVAssert(Cjj >= C._first);
                        TMVAssert(Cjj < C._last);
#endif
                        *Cjj = newcjj;
                    }
                } else {
                    const Tb* Bjj = B.cptr();
                    T* Cjj = C.ptr();
                    const ptrdiff_t Bds = B.stepi()+B.stepj();
                    const ptrdiff_t Cds = C.stepi()+C.stepj();
                    for(ptrdiff_t j=0,jj=1;j<N;++j,++jj,Bjj+=Bds,Cjj+=Cds) {
                        T xBjj = B.isconj()?TMV_CONJ(*Bjj):*Bjj;
                        T newcjj = A.row(j,jj,N)*B.col(j,jj,N) + xBjj;
                        if (alpha != T(1)) {
                            newcjj *= alpha;
                            xBjj *= alpha; // xBjj is now alpha*B(j,j)
                        }
                        if (add) newcjj += *Cjj;

                        if (add) C.col(j,0,j) += xBjj * A.col(j,0,j);
                        else C.col(j,0,j) = xBjj * A.col(j,0,j);

                        C.col(j,0,j) += 
                            alpha * A.subMatrix(0,j,jj,N) * B.col(j,jj,N);

                        MultMV<add>(alpha,A.subTriMatrix(jj,N),B.col(j,jj,N),
                                    C.col(j,jj,N));

#ifdef TMVFLDEBUG
                        TMVAssert(Cjj >= C._first);
                        TMVAssert(Cjj < C._last);
#endif
                        *Cjj = newcjj;
                    }
                }
            } else {
                if (B.isunit()) {
                    const Ta* Ajj = A.cptr();
                    T* Cjj = C.ptr();
                    const ptrdiff_t Ads = A.stepi()+A.stepj();
                    const ptrdiff_t Cds = C.stepi()+C.stepj();
                    for(ptrdiff_t j=0,jj=1;j<N;++j,++jj,Ajj+=Ads,Cjj+=Cds) {
                        Ta xAjj = A.isconj() ? TMV_CONJ(*Ajj) : *Ajj;
                        T newcjj = A.row(j,jj,N)*B.col(j,jj,N) + xAjj;
                        if (alpha != T(1)) newcjj *= alpha;
                        if (add) newcjj += *Cjj;

                        if (add) C.col(j,0,j) += alpha * A.col(j,0,j);
                        else C.col(j,0,j) = alpha * A.col(j,0,j);

                        C.col(j,0,j) += 
                            alpha * A.subMatrix(0,j,jj,N) * B.col(j,jj,N);

                        MultMV<add>(alpha,A.subTriMatrix(jj,N),B.col(j,jj,N),
                                    C.col(j,jj,N));

#ifdef TMVFLDEBUG
                        TMVAssert(Cjj >= C._first);
                        TMVAssert(Cjj < C._last);
#endif
                        *Cjj = newcjj;
                    }
                } else {
                    const Tb* Bjj = B.cptr();
                    T* Cjj = C.ptr();
                    const ptrdiff_t Bds = B.stepi()+B.stepj();
                    const ptrdiff_t Cds = C.stepi()+C.stepj();
                    for(ptrdiff_t j=0,jj=1;j<N;++j,++jj,Bjj+=Bds,Cjj+=Cds) {
                        T xBjj = B.isconj() ? TMV_CONJ(*Bjj) : *Bjj;
                        T newcjj = A.row(j,j,N)*B.col(j,j,N);
                        if (alpha != T(1)) {
                            newcjj *= alpha;
                            xBjj *= alpha;
                        }
                        if (add) newcjj += *Cjj;

                        if (add) C.col(j,0,j) += xBjj * A.col(j,0,j);
                        else C.col(j,0,j) = xBjj * A.col(j,0,j);

                        C.col(j,0,j) += 
                            alpha * A.subMatrix(0,j,jj,N) * B.col(j,jj,N);

                        MultMV<add>(alpha,A.subTriMatrix(jj,N),B.col(j,jj,N),
                                    C.col(j,jj,N));

#ifdef TMVFLDEBUG
                        TMVAssert(Cjj >= C._first);
                        TMVAssert(Cjj < C._last);
#endif
                        *Cjj = newcjj;
                    }
                }
            }
        }
#ifdef XDEBUG
        if (!(Norm(C2-C) <=
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                     (add?Norm(C0):TMV_RealType(T)(0))))) {
            cerr<<"ColMultMM: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C0<<endl;
            cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
            cerr<<", Cptr = "<<C.cptr()<<endl;
            cerr<<"--> C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    static void DoMultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, MatrixView<T> C)
    // C (+)= alpha * A * B
    // This is designed to work even if A,B are in same storage as C
    {
#ifdef XDEBUG
        //cout<<"MultMM: "<<alpha<<"  "<<TMV_Text(A)<<"  "<<A<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        //cout<<TMV_Text(C)<<"  "<<C<<endl;
        Matrix<T> C2 = C;
        Matrix<T> C0 = C;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        if (add) C2 += alpha*A0*B0;
        else C2 = alpha*A0*B0;
#endif
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.size() == C.rowsize());
        TMVAssert(A.size() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(C.ct() == NonConj);

        const ptrdiff_t nb = TRI_MM_BLOCKSIZE;
        const ptrdiff_t N = A.size();

        if (N <= TRI_MM_BLOCKSIZE2) {
            if (C.isrm()) 
                ColMultMM<add>(
                    alpha,B.transpose(),A.transpose(),C.transpose());
            else ColMultMM<add>(alpha,A,B,C);
        } else {
            ptrdiff_t k = N/2;
            if (k > nb) k = k/nb*nb;

            // [ A00 A01 ] [ B00  0  ] = [ A00 B00 + A01 B10   A01 B11 ]
            // [  0  A11 ] [ B10 B11 ]   [      A11 B10        A11 B11 ]

            ConstUpperTriMatrixView<Ta> A00 = A.subTriMatrix(0,k);
            ConstMatrixView<Ta> A01 = A.subMatrix(0,k,k,N);
            ConstUpperTriMatrixView<Ta> A11 = A.subTriMatrix(k,N);
            ConstLowerTriMatrixView<Tb> B00 = B.subTriMatrix(0,k);
            ConstMatrixView<Tb> B10 = B.subMatrix(k,N,0,k);
            ConstLowerTriMatrixView<Tb> B11 = B.subTriMatrix(k,N);
            MatrixView<T> C00 = C.subMatrix(0,k,0,k);
            MatrixView<T> C01 = C.subMatrix(0,k,k,N);
            MatrixView<T> C10 = C.subMatrix(k,N,0,k);
            MatrixView<T> C11 = C.subMatrix(k,N,k,N);

            DoMultMM<add>(alpha,A00,B00,C00);
            C00 += alpha*A01*B10;
            if (SameStorage(A01,C10)) {
                if (SameStorage(B10,C01)) {
                    // This is the only case where we need temporary storage,
                    // and I don't image that it is often needed, but it's 
                    // worth checking.
                    Matrix<T> A01x = A01;
                    MultMM<add>(alpha,A11,B10,C10);
                    MultMM<add>(alpha,B11.transpose(),A01x.transpose(),
                                C01.transpose());
                } else {
                    MultMM<add>(alpha,B11.transpose(),A01.transpose(),
                                C01.transpose());
                    MultMM<add>(alpha,A11,B10,C10);
                }
            } else {
                MultMM<add>(alpha,A11,B10,C10);
                MultMM<add>(alpha,B11.transpose(),A01.transpose(),
                            C01.transpose());
            }
            DoMultMM<add>(alpha,A11,B11,C11);
        }
#ifdef XDEBUG
        if (!(Norm(C2-C) <=
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                     (add?Norm(C0):TMV_RealType(T)(0))))) {
            cerr<<"DoMultMM: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C0<<endl;
            cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
            cerr<<", Cptr = "<<C.cptr()<<endl;
            cerr<<"--> C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, MatrixView<T> C)
    {
#ifdef XDEBUG
        //cout<<"MultMM: "<<alpha<<"  "<<TMV_Text(A)<<"  "<<A<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        //cout<<TMV_Text(C)<<"  "<<C<<endl;
        Matrix<T> C2 = C;
        Matrix<T> C0 = C;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        if (add) C2 += alpha*A0*B0;
        else C2 = alpha*A0*B0;
#endif
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.size() == C.rowsize());

        const ptrdiff_t N = A.size();

        if (N==0) return;
        else if (alpha == T(0)) {
            if (!add) C.setZero();
        }
        else if (C.isconj()) 
            DoMultMM<add>(
                TMV_CONJ(alpha),A.conjugate(),B.conjugate(),C.conjugate());
        else
            DoMultMM<add>(alpha,A,B,C);

#ifdef XDEBUG
        if (!(Norm(C2-C) <=
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                     (add?Norm(C0):TMV_RealType(T)(0))))) {
            cerr<<"MultMM: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C0<<endl;
            cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
            cerr<<", Cptr = "<<C.cptr()<<endl;
            cerr<<"--> C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif

    }

    //
    // MultMM: M = L * U
    //

    template <bool add, class T, class Ta, class Tb> 
    static void ColMultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, MatrixView<T> C)
    {
#ifdef XDEBUG
        //cout<<"MultMM: "<<alpha<<"  "<<TMV_Text(A)<<"  "<<A<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        //cout<<TMV_Text(C)<<"  "<<C<<endl;
        Matrix<T> C0 = C;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C2 = C;
        if (add) C2 += alpha*A0*B0;
        else C2 = alpha*A0*B0;
        //cout<<"Start ColMultMM: L * U\n";
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<"B = "<<TMV_Text(B)<<"  "<<B<<endl;
        //cout<<"C = "<<TMV_Text(C)<<"  "<<C<<endl;
        //cout<<"Correct result = "<<C2<<endl;
#endif

        if (SameStorage(A,C) && A.stepi() == C.stepj()) {
            // Then need temporary (see below)
            if (A.isrm()) {
                LowerTriMatrix<Ta,NonUnitDiag|RowMajor> A2 = A;
                ColMultMM<add>(alpha,A2,B,C);
            } else {
                LowerTriMatrix<Ta,NonUnitDiag|ColMajor> A2 = A;
                ColMultMM<add>(alpha,A2,B,C);
            }
        } else {
            if (A.isunit()) {
                if (B.isunit()) {
                    const ptrdiff_t Cds = C.stepi()+C.stepj();
                    T* Cjj = C.ptr()+(C.rowsize()-1)*Cds;
                    for(ptrdiff_t jj=C.rowsize(),j=jj-1;jj>0;--jj,--j, Cjj-=Cds) {
                        // jj = j+1
                        // C.col(j) (+=) alpha*A*B.col(j)
                        //
                        // C.col(j) (+=) alpha*A.colRange(0,j+1)*B.col(j,0,j+1)
                        //
                        // C(j,j) (+=) alpha*A.row(j,0,j+1)*B.col(j,0,j+1)
                        // C.col(j,j+1,N) (+=) 
                        //     alpha*A.subMatrix(j+1,N,0,j+1)*B.col(j,0,j+1)
                        // C.col(j,0,j) (+=) 
                        //     alpha*A.subTriMatrix(0,j)*B.col(j,0,j)
                        //
                        // Requirements on storage: 
                        //   B can be stored in either triangle
                        //   A cannot be stored in C's upper triangle

                        const ptrdiff_t N = A.size();

                        T newcjj = A.row(j,0,j)*B.col(j,0,j) + T(1);
                        if (alpha != T(1)) newcjj *= alpha;
                        if (add) newcjj += *Cjj;

                        if (add) C.col(j,jj,N) += alpha * A.col(j,jj,N);
                        else C.col(j,jj,N) = alpha * A.col(j,jj,N);

                        C.col(j,jj,N) += 
                            alpha * A.subMatrix(jj,N,0,j) * B.col(j,0,j);

                        MultMV<add>(alpha,A.subTriMatrix(0,j),B.col(j,0,j),
                                    C.col(j,0,j));

#ifdef TMVFLDEBUG
                        TMVAssert(Cjj >= C._first);
                        TMVAssert(Cjj < C._last);
#endif
                        *Cjj = newcjj;
                    }
                } else {
                    const ptrdiff_t Nm1 = C.rowsize()-1;
                    const ptrdiff_t Bds = B.stepi()+B.stepj();
                    const ptrdiff_t Cds = C.stepi()+C.stepj();
                    const Tb* Bjj = B.cptr()+Nm1*Bds;
                    T* Cjj = C.ptr()+Nm1*Cds;
                    for(ptrdiff_t jj=C.rowsize(),j=jj-1;
                        jj>0;
                        --jj,--j, Bjj-=Bds,Cjj-=Cds) {

                        const ptrdiff_t N = A.size();

                        T xBjj = B.isconj()?TMV_CONJ(*Bjj):*Bjj;
                        T newcjj = A.row(j,0,j)*B.col(j,0,j) + xBjj;
                        if (alpha != T(1)) {
                            newcjj *= alpha;
                            xBjj *= alpha;
                        }
                        if (add) newcjj += *Cjj;

                        if (add) C.col(j,jj,N) += xBjj * A.col(j,jj,N);
                        else C.col(j,jj,N) = xBjj * A.col(j,jj,N);

                        C.col(j,jj,N) += 
                            alpha * A.subMatrix(jj,N,0,j) * B.col(j,0,j);

                        MultMV<add>(alpha,A.subTriMatrix(0,j),B.col(j,0,j),
                                    C.col(j,0,j));

#ifdef TMVFLDEBUG
                        TMVAssert(Cjj >= C._first);
                        TMVAssert(Cjj < C._last);
#endif
                        *Cjj = newcjj;
                    }
                }
            } else {
                if (B.isunit()) {
                    const ptrdiff_t Nm1 = C.rowsize()-1;
                    const ptrdiff_t Ads = A.stepi()+A.stepj();
                    const ptrdiff_t Cds = C.stepi()+C.stepj();
                    const Ta* Ajj = A.cptr()+Nm1*Ads;
                    T* Cjj = C.ptr()+Nm1*Cds;
                    for(ptrdiff_t jj=C.rowsize(),j=jj-1;
                        jj>0;
                        --jj,--j, Ajj-=Ads,Cjj-=Cds) {

                        const ptrdiff_t N = A.size();

                        T newcjj = 
                            A.row(j,0,j)*B.col(j,0,j) +
                            (A.isconj()?TMV_CONJ(*Ajj):*Ajj);
                        if (alpha != T(1)) newcjj *= alpha;
                        if (add) newcjj += *Cjj;

                        if (add) C.col(j,jj,N) += alpha * A.col(j,jj,N);
                        else C.col(j,jj,N) = alpha * A.col(j,jj,N);

                        C.col(j,jj,N) += 
                            alpha * A.subMatrix(jj,N,0,j) * B.col(j,0,j);

                        MultMV<add>(alpha,A.subTriMatrix(0,j),B.col(j,0,j),
                                    C.col(j,0,j));

#ifdef TMVFLDEBUG
                        TMVAssert(Cjj >= C._first);
                        TMVAssert(Cjj < C._last);
#endif
                        *Cjj = newcjj;
                    }
                } else {
                    const ptrdiff_t Nm1 = C.rowsize()-1;
                    const ptrdiff_t Bds = B.stepi()+B.stepj();
                    const ptrdiff_t Cds = C.stepi()+C.stepj();
                    const Tb* Bjj = B.cptr()+Nm1*Bds;
                    T* Cjj = C.ptr()+Nm1*Cds;
                    for(ptrdiff_t jj=C.rowsize(),j=jj-1;
                        jj>0;
                        --jj,--j, Bjj-=Bds,Cjj-=Cds) {

                        const ptrdiff_t N = A.size();

                        T xBjj = B.isconj() ? TMV_CONJ(*Bjj):*Bjj;
                        T newcjj = A.row(j,0,j+1)*B.col(j,0,j+1);
                        if (alpha != T(1)) {
                            newcjj *= alpha;
                            xBjj *= alpha;
                        }
                        if (add) newcjj += *Cjj;

                        if (add) C.col(j,jj,N) += xBjj * A.col(j,jj,N);
                        else C.col(j,jj,N) = xBjj * A.col(j,jj,N);
                        C.col(j,jj,N) += 
                            alpha * A.subMatrix(jj,N,0,j) * B.col(j,0,j);

                        MultMV<add>(alpha,A.subTriMatrix(0,j),B.col(j,0,j),
                                    C.col(j,0,j));

#ifdef TMVFLDEBUG
                        TMVAssert(Cjj >= C._first);
                        TMVAssert(Cjj < C._last);
#endif
                        *Cjj = newcjj;
                    }
                }
            }
        }
#ifdef XDEBUG
        if (!(Norm(C2-C) <=
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                     (add?Norm(C0):TMV_RealType(T)(0))))) {
            cerr<<"MultMM: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C0<<endl;
            cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
            cerr<<", Cptr = "<<C.cptr()<<endl;
            cerr<<"--> C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    static void DoMultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, MatrixView<T> C)
    // C (+)= alpha * A * B
    // This is designed to work even if A,B are in same storage as C
    {
#ifdef XDEBUG
        //cout<<"MultMM: "<<alpha<<"  "<<TMV_Text(A)<<"  "<<A<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        //cout<<TMV_Text(C)<<"  "<<C<<endl;
        Matrix<T> C0 = C;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C2 = C;
        if (add) C2 += alpha*A0*B0;
        else C2 = alpha*A0*B0;
#endif
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.size() == C.rowsize());
        TMVAssert(A.size() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(C.ct() == NonConj);

        const ptrdiff_t nb = TRI_MM_BLOCKSIZE;
        const ptrdiff_t N = A.size();

        if (N <= TRI_MM_BLOCKSIZE2) {
            if (C.isrm()) 
                ColMultMM<add>(
                    alpha,B.transpose(),A.transpose(),C.transpose());
            else
                ColMultMM<add>(alpha,A,B,C);
        } else {
            ptrdiff_t k = N/2;
            if (k > nb) k = k/nb*nb;

            // [ A00  0  ] [ B00 B01 ] = [ A00 B00       A00 B01      ]
            // [ A10 A11 ] [  0  B11 ]   [ A10 B00  A10 B01 + A11 B11 ]

            ConstLowerTriMatrixView<Ta> A00 = A.subTriMatrix(0,k);
            ConstMatrixView<Ta> A10 = A.subMatrix(k,N,0,k);
            ConstLowerTriMatrixView<Ta> A11 = A.subTriMatrix(k,N);
            ConstUpperTriMatrixView<Tb> B00 = B.subTriMatrix(0,k);
            ConstMatrixView<Tb> B01 = B.subMatrix(0,k,k,N);
            ConstUpperTriMatrixView<Tb> B11 = B.subTriMatrix(k,N);
            MatrixView<T> C00 = C.subMatrix(0,k,0,k);
            MatrixView<T> C01 = C.subMatrix(0,k,k,N);
            MatrixView<T> C10 = C.subMatrix(k,N,0,k);
            MatrixView<T> C11 = C.subMatrix(k,N,k,N);

            DoMultMM<add>(alpha,A11,B11,C11);
            C11 += alpha*A10*B01;
            if (SameStorage(A10,C01)) {
                if (SameStorage(B01,C10)) {
                    // This is the only case where we need temporary storage,
                    // and I don't image that it is often needed, but it's 
                    // worth checking.
                    Matrix<T> A10x = A10;
                    MultMM<add>(alpha,A00,B01,C01);
                    MultMM<add>(alpha,B00.transpose(),A10x.transpose(),
                                C10.transpose());
                } else {
                    MultMM<add>(alpha,B00.transpose(),A10.transpose(),
                                C10.transpose());
                    MultMM<add>(alpha,A00,B01,C01);
                }
            } else {
                MultMM<add>(alpha,A00,B01,C01);
                MultMM<add>(alpha,B00.transpose(),A10.transpose(),
                            C10.transpose());
            }
            DoMultMM<add>(alpha,A00,B00,C00);
        }
#ifdef XDEBUG
        if (!(Norm(C2-C) <= 
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                     (add?Norm(C0):TMV_RealType(T)(0))))) {
            cerr<<"MultMM: alpha= "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C0<<endl;
            cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
            cerr<<", Cptr = "<<C.cptr()<<endl;
            cerr<<"--> C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, MatrixView<T> C)
    {
#ifdef XDEBUG
        //cout<<"MultMM: "<<alpha<<"  "<<TMV_Text(A)<<"  "<<A<<"  "<<TMV_Text(B)<<"  "<<B<<endl;
        //cout<<TMV_Text(C)<<"  "<<C<<endl;
        Matrix<T> C0 = C;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C2 = C;
        if (add) C2 += alpha*A0*B0;
        else C2 = alpha*A0*B0;
#endif
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.size() == C.rowsize());

        const ptrdiff_t N = A.size();

        if (N==0) return;
        else if (alpha == T(0)) {
            if (!add) C.setZero();
        } else if (C.isconj()) {
            DoMultMM<add>(
                TMV_CONJ(alpha),A.conjugate(),B.conjugate(),C.conjugate());
        } else {
            DoMultMM<add>(alpha,A,B,C);
        }

#ifdef XDEBUG
        if (!(Norm(C2-C) <= 
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                     (add?Norm(C0):TMV_RealType(T)(0))))) {
            cerr<<"MultMM: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C0<<endl;
            cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
            cerr<<", Cptr = "<<C.cptr()<<endl;
            cerr<<"--> C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_MultUL.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


