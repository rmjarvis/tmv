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


#include "tmv/TMV_SymMatrixArithFunc.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"

#ifdef XDEBUG
#include <iostream>
#include "tmv/TMV_SymMatrixArith.h"
using std::cout;
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
    // SymMultMM
    // A += alpha * x * y
    // where x,y are Matrices, and the product x*y is assumed to be symmetric.
    //

    template <bool ha, bool a1, bool add, class T, class Tx, class Ty> 
    static void RecursiveSymMultMM(
        const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
        SymMatrixView<T> A)
    {
        TMVAssert(A.size() == x.colsize());
        TMVAssert(A.size() == y.rowsize());
        TMVAssert(x.rowsize() == y.colsize());
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
                T temp = x.row(0) * y.col(0);
                if (!a1) temp *= alpha;
#ifdef TMVFLDEBUG
                TMVAssert(A.ptr() >= A._first);
                TMVAssert(A.ptr() < A._last);
#endif
                if (ha)
                    if (add) *(A.ptr()) += TMV_REAL(temp);
                    else *(A.ptr()) = TMV_REAL(temp);
                else
                    if (add) *(A.ptr()) += temp;
                    else *(A.ptr()) = temp;
            } else {
                if (A.isrm()) {
                    for (ptrdiff_t i=0;i<N;++i) {
                        if (add) 
                            A.row(i,0,i+1) += 
                                alpha * x.row(i) * y.colRange(0,i+1);
                        else 
                            A.row(i,0,i+1) = 
                                alpha * x.row(i) * y.colRange(0,i+1);
                    }
                } else {
                    for (ptrdiff_t j=0;j<N;++j) {
                        if (add) 
                            A.col(j,j,N) += alpha * x.rowRange(j,N) * y.col(j);
                        else 
                            A.col(j,j,N) = alpha * x.rowRange(j,N) * y.col(j);
                    }
                }
                if (ha && isComplex(T())) {
#ifdef XDEBUG
                    TMVAssert(NormInf(A.diag().imagPart()) <
                              A.size()*TMV_Epsilon<T>()*
                              (Norm(A)+Norm(x)*Norm(y)));
#endif
                    A.diag().imagPart().setZero();
                }
            }
        } else { // Not <= BLOCKSIZE2, so do recurse...
            ptrdiff_t k = N/2;
            if (k > nb) k = k/nb*nb;

            RecursiveSymMultMM<ha,a1,add>(
                alpha,x.rowRange(0,k),y.colRange(0,k),A.subSymMatrix(0,k));

            if (add) 
                A.subMatrix(k,N,0,k) += 
                    alpha * x.rowRange(k,N) * y.colRange(0,k);
            else 
                A.subMatrix(k,N,0,k) = 
                    alpha * x.rowRange(k,N) * y.colRange(0,k);

            RecursiveSymMultMM<ha,a1,add>(
                alpha,x.rowRange(k,N),y.colRange(k,N), A.subSymMatrix(k,N));
        }
    }

    template <bool ha, bool a1, bool add, class T, class Tx, class Ty> 
    static void RecursiveInPlaceSymMultMM(
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
            if (ha) {
                TMV_RealType(T) temp = 
                    a1 ? TMV_REAL(x00*y00) : TMV_REAL(alpha*x00*y00);
#ifdef TMVFLDEBUG
                TMVAssert(A.ptr() >= A._first);
                TMVAssert(A.ptr() < A._last);
#endif
                if (add) *A.ptr() += temp;
                else *A.ptr() = temp;
            } else {
                T temp = a1 ? (x00 * y00) : (alpha * (x00 * y00));
#ifdef TMVFLDEBUG
                TMVAssert(A.ptr() >= A._first);
                TMVAssert(A.ptr() < A._last);
#endif
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

            Matrix<T> tempA10 = x10 * y00;
            tempA10 += x11 * y10;

            RecursiveInPlaceSymMultMM<ha,a1,add>(alpha,x11,y11,A11);
            RecursiveSymMultMM<ha,a1,true>(alpha,x10,y01,A11);
            RecursiveInPlaceSymMultMM<ha,a1,add>(alpha,x00,y00,A00);
            RecursiveSymMultMM<ha,a1,true>(alpha,x01,y10,A00);

            if (add) A10 += alpha * tempA10;
            else A10 = alpha * tempA10;
        }
    }

    template <bool add, class T, class Tx, class Ty> 
    static inline void InPlaceSymMultMM(
        const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
        SymMatrixView<T> A)
    {
        if (A.isherm())
            if (alpha == T(1))
                RecursiveInPlaceSymMultMM<true,true,add>(alpha,x,y,A); 
            else
                RecursiveInPlaceSymMultMM<true,false,add>(alpha,x,y,A); 
        else
            if (alpha == T(1))
                RecursiveInPlaceSymMultMM<false,true,add>(alpha,x,y,A); 
            else
                RecursiveInPlaceSymMultMM<false,false,add>(alpha,x,y,A); 
    }

    template <bool add, class T, class Tx, class Ty> 
    static void DoSymMultMM(
        const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
        SymMatrixView<T> A)
    { 
        if (SameStorage(x,A) || SameStorage(y,A)) {
            const ptrdiff_t N = A.size();
            TMVAssert(x.colsize() == N);
            TMVAssert(y.rowsize() == N);
            const ptrdiff_t K = x.rowsize();
            TMVAssert(y.colsize() == K);

            if (K >= N) {
                InPlaceSymMultMM<add>(alpha,x.colRange(0,N),y.rowRange(0,N),A);
                if (K > N) 
                    DoSymMultMM<add>(alpha,x.colRange(N,K),y.rowRange(N,K),A);
            } else { 
                DoSymMultMM<add>(alpha,x.rowRange(K,N),y.colRange(K,N),
                                 A.subSymMatrix(K,N));
                if ((y.stepi() < y.stepj()) != (A.stepi() < A.stepj())) {
                    if ((x.stepi() < x.stepj()) == (A.stepi() < A.stepj())) {
                        // Then both x and y overlap with A(K:N,0:K)
                        // Need a temporary
                        if (A.iscm()) {
                            Matrix<T,ColMajor> temp =
                                x.rowRange(K,N) * y.colRange(0,K);
                            if (add) A.subMatrix(K,N,0,K) += alpha * temp;
                            else A.subMatrix(K,N,0,K) = alpha * temp;
                        } else {
                            Matrix<T,RowMajor> temp =
                                x.rowRange(K,N) * y.colRange(0,K);
                            if (add) A.subMatrix(K,N,0,K) += alpha * temp;
                            else A.subMatrix(K,N,0,K) = alpha * temp;
                        }
                    } else {
                        MultMM<add>(
                            alpha,x.subMatrix(K,N,K,N),y.subMatrix(K,N,0,K),
                            A.subMatrix(K,N,0,K) );
                        MultMM<true>(
                            alpha,x.subMatrix(K,N,0,K),y.subMatrix(0,K,0,K),
                            A.subMatrix(K,N,0,K) );
                    }
                }
                else {
                    MultMM<add>(
                        alpha,x.rowRange(K,N),y.colRange(0,K),
                        A.subMatrix(K,N,0,K) );
                }
                InPlaceSymMultMM<add>(
                    alpha,x.rowRange(0,K),y.rowRange(0,K),A.subSymMatrix(0,K));
            }
        } else {
            if (A.isherm())
                if (alpha == T(1))
                    RecursiveSymMultMM<true,true,add>(alpha,x,y,A); 
                else
                    RecursiveSymMultMM<true,false,add>(alpha,x,y,A); 
            else
                if (alpha == T(1))
                    RecursiveSymMultMM<false,true,add>(alpha,x,y,A); 
                else
                    RecursiveSymMultMM<false,false,add>(alpha,x,y,A); 
        }
    }

    template <bool add, class T, class Tx, class Ty> 
    void SymMultMM(
        const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
        SymMatrixView<T> A)
    // A (+)= alpha * x * y
    {
        TMVAssert(A.size() == x.colsize());
        TMVAssert(A.size() == y.rowsize());
        TMVAssert(x.rowsize() == y.colsize());

#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<Tx> x0 = x;
        Matrix<Ty> y0 = y;
        Matrix<T> A2 = A;
        if (add) A2 += alpha*x*y;
        else A2 = alpha*x*y;
        //cout<<"Start SymMultMM: alpha = "<<alpha<<endl;
        //cout<<"add = "<<add<<endl;
        //cout<<"x = "<<TMV_Text(x)<<"  "<<x0<<endl;
        //cout<<"y = "<<TMV_Text(y)<<"  "<<y0<<endl;
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
#endif

        if (alpha != T(0) && A.size() > 0) {
            if (A.uplo() == Upper) {
                if (A.isherm()) SymMultMM<add>(alpha,x,y,A.adjoint());
                else SymMultMM<add>(alpha,x,y,A.transpose());
            }
            else if (A.isconj()) 
                SymMultMM<add>(
                    TMV_CONJ(alpha),x.conjugate(),y.conjugate(),A.conjugate());
            else DoSymMultMM<add>(alpha,x,y,A);
        }

#ifdef XDEBUG
        if (Norm(A-A2) > 0.001*(TMV_ABS(alpha)*Norm(x0)*Norm(y0)+
                                (add?Norm(A0):TMV_RealType(T)(0)))) {
            cerr<<"SymMultMM: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"x = "<<TMV_Text(x)<<"  "<<x0<<endl;
            cerr<<"y = "<<TMV_Text(y)<<"  "<<y0<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"x*y = "<<x0*y0<<endl;
            cerr<<"-> A = "<<A<<endl;
            cerr<<"A2 = "<<A2<<endl;
            abort();
        }
#endif
        TMVAssert(A.isHermOK());
    }

#define InstFile "TMV_SymMultMMS.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


