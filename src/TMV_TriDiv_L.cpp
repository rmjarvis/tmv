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


#include "tmv/TMV_TriDiv.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef NOTHROW
#include <iostream>
#endif

#ifdef XDEBUG
#include <iostream>
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

    template <class T, class Ta> 
    static void RowTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, UpperTriMatrixView<T> B)
    // B = A^-1 * B
    // where A is a triangle matrix
    {
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || A.isunit());
        TMVAssert(B.size()>0);
        TMVAssert(B.ct() == NonConj);
        const ptrdiff_t N = B.size();

        if (A.isunit()) {
            for(ptrdiff_t i=N-1; i>=0; --i) 
                B.row(i,i+1,N) -= A.row(i,i+1,N) * B.subTriMatrix(i+1,N);
        } else {
            TMVAssert(!B.isunit());
            const ptrdiff_t Ads = A.stepi()+A.stepj();
            const Ta* Aii = A.cptr() + (N-1)*Ads;
            ptrdiff_t len=1;
            for(ptrdiff_t i=N-1; i>=0; --i,Aii-=Ads,++len) {
                B.row(i,i+1,N) -= A.row(i,i+1,N) * B.subTriMatrix(i+1,N);
                if (*Aii==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular UpperTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularUpperTriMatrix<Ta>(A);
#endif
                }
                if (*Aii != Ta(1)) 
                    B.row(i,i,N) /= (A.isconj()?TMV_CONJ(*Aii):*Aii);
            }
        }
    }

    template <class T, class Ta> 
    static void ColTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, UpperTriMatrixView<T> B)
    // B = A^-1 * B
    // where A is a triangle matrix
    {
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || A.isunit());
        TMVAssert(B.size()>0);
        TMVAssert(B.ct() == NonConj);
        const ptrdiff_t N = B.size();

        if (A.isunit()) {
            if (B.isunit()) 
                for(ptrdiff_t j=N-1; j>=0; --j) {
                    B.subMatrix(0,j,j+1,N) -= A.col(j,0,j) ^ B.row(j,j+1,N);
                    B.col(j,0,j) -= A.col(j,0,j);
                }
            else
                for(ptrdiff_t j=N-1; j>=0; --j) 
                    B.subMatrix(0,j,j,N) -= A.col(j,0,j) ^ B.row(j,j,N);
        } else {
            TMVAssert(!B.isunit());
            const ptrdiff_t Ads = A.stepi()+A.stepj();
            const Ta* Ajj = A.cptr() + (N-1)*Ads;
            ptrdiff_t len=1;
            for(ptrdiff_t j=N-1; j>=0; --j,Ajj-=Ads,++len) {
                if (*Ajj==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular UpperTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularUpperTriMatrix<Ta>(A);
#endif
                }
                if (*Ajj != Ta(1)) 
                    B.row(j,j,N) /= (A.isconj()?TMV_CONJ(*Ajj):*Ajj);
                B.subMatrix(0,j,j,N) -= A.col(j,0,j) ^ B.row(j,j,N);
            }
        } 
    }

    template <class T, class Ta> 
    static void RowTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, LowerTriMatrixView<T> B)
    // B = A^-1 * B
    // where A is a triangle matrix
    {
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || A.isunit());
        TMVAssert(B.size()>0);
        TMVAssert(B.ct() == NonConj);
        const ptrdiff_t N = B.size();

        if (A.isunit()) {
            for(ptrdiff_t i=0;i<N;++i) 
                B.row(i,0,i) -= A.row(i,0,i) * B.subTriMatrix(0,i);
        } else {
            TMVAssert(!B.isunit());
            const ptrdiff_t Ads = A.stepi()+A.stepj();
            const Ta* Aii = A.cptr();
            for(ptrdiff_t i=0;i<N;++i,Aii+=Ads) {
                B.row(i,0,i) -= A.row(i,0,i) * B.subTriMatrix(0,i);
                if (*Aii==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular LowerTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularLowerTriMatrix<Ta>(A);
#endif
                }
                if (*Aii != Ta(1)) 
                    B.row(i,0,i+1) /= (A.isconj()?TMV_CONJ(*Aii):*Aii);
            }
        }
    }

    template <class T, class Ta> 
    static void ColTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, LowerTriMatrixView<T> B)
    // B = A^-1 * B
    // where A is a triangle matrix
    {
        TMVAssert(B.isrm());
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || A.isunit());
        TMVAssert(B.size()>0);
        TMVAssert(B.ct() == NonConj);
        const ptrdiff_t N = B.size();

        if (A.isunit()) {
            if (B.isunit())
                for(ptrdiff_t j=0;j<N;++j) {
                    B.col(j,j+1,N) -= A.col(j,j+1,N);
                    B.subMatrix(j+1,N,0,j) -= A.col(j,j+1,N) ^ B.row(j,0,j);
                }
            else
                for(ptrdiff_t j=0;j<N;++j) 
                    B.subMatrix(j+1,N,0,j+1) -= A.col(j,j+1,N) ^ B.row(j,0,j+1);
        } else {
            TMVAssert(!B.isunit());
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
                if (*Ajj != Ta(1)) 
                    B.row(j,0,j+1) /= (A.isconj()?TMV_CONJ(*Ajj):*Ajj);
                B.subMatrix(j+1,N,0,j+1) -= A.col(j,j+1,N) ^ B.row(j,0,j+1);
            }
        }
    }

    template <class T, class Ta> 
    static void DoTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, UpperTriMatrixView<T> B)
    {
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || A.isunit());
        TMVAssert(B.size()>0);
        TMVAssert(B.ct() == NonConj);

        const ptrdiff_t nb = TRI_DIV_BLOCKSIZE;
        const ptrdiff_t N = B.size();

        if (N <= TRI_DIV_BLOCKSIZE2) {
            if (B.isrm()) {
                if (A.isrm()) RowTriLDivEq(A,B);
                else ColTriLDivEq(A,B);
            } else {
                if (B.isunit())
                    for(ptrdiff_t j=0;j<N;++j) {
                        B.col(j,0,j) -= A.col(j,0,j);
                        TriLDivEq(A.subTriMatrix(0,j),B.col(j,0,j));
                    }
                else // B is NonUnitDiag
                    for(ptrdiff_t j=0;j<N;++j)
                        TriLDivEq(A.subTriMatrix(0,j+1),B.col(j,0,j+1));
            }
        } else {
            ptrdiff_t k = N/2;
            if (k > nb) k = k/nb*nb;

            ConstUpperTriMatrixView<Ta> A00 = A.subTriMatrix(0,k);
            ConstMatrixView<Ta> A01 = A.subMatrix(0,k,k,N);
            ConstUpperTriMatrixView<Ta> A11 = A.subTriMatrix(k,N);
            UpperTriMatrixView<T> B00 = B.subTriMatrix(0,k);
            MatrixView<T> B01 = B.subMatrix(0,k,k,N);
            UpperTriMatrixView<T> B11 = B.subTriMatrix(k,N);

            DoTriLDivEq(A11,B11);
            B01 -= A01 * B11;
            TriLDivEq(A00,B01);
            DoTriLDivEq(A00,B00);
        }
    }

    template <class T, class Ta> 
    void TriLDivEq(
        const GenUpperTriMatrix<Ta>& A, UpperTriMatrixView<T> B)
    // B = A^-1 * B
    // where A is a triangle matrix
    {
#ifdef XDEBUG
        Matrix<Ta> A0(A);
        Matrix<T> B0(B);
#endif
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || A.isunit());

        if (B.isconj()) TriLDivEq(A.conjugate(),B.conjugate());
        else DoTriLDivEq(A,B);

#ifdef XDEBUG
        Matrix<T> BB = A0*B;
        if (Norm(BB-B0) > 0.001*Norm(A0)*Norm(B0)) {
            cerr<<"TriLDivEq: Upper/Upper\n";
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"Done: B = "<<B<<endl;
            cerr<<"A*B = "<<BB<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta> 
    static void DoTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, LowerTriMatrixView<T> B)
    {
        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || A.isunit());
        TMVAssert(B.size()>0);
        TMVAssert(B.ct() == NonConj);

        const ptrdiff_t nb = TRI_DIV_BLOCKSIZE;
        const ptrdiff_t N = B.size();

        if (N <= TRI_DIV_BLOCKSIZE2) {
            if (B.isrm()) {
                if (A.isrm()) RowTriLDivEq(A,B);
                else ColTriLDivEq(A,B);
            } else {
                if (B.isunit())
                    for(ptrdiff_t j=0;j<N;++j) {
                        B.col(j,j+1,N) -= A.col(j,j+1,N);
                        TriLDivEq(A.subTriMatrix(j+1,N),B.col(j,j+1,N));
                    }
                else // B is NonUnitDiag
                    for(ptrdiff_t j=0;j<N;++j)
                        TriLDivEq(A.subTriMatrix(j,N),B.col(j,j,N));
            }
        } else {
            ptrdiff_t k = N/2;
            if (k > nb) k = k/nb*nb;

            ConstLowerTriMatrixView<Ta> A00 = A.subTriMatrix(0,k);
            ConstMatrixView<Ta> A10 = A.subMatrix(k,N,0,k);
            ConstLowerTriMatrixView<Ta> A11 = A.subTriMatrix(k,N);
            LowerTriMatrixView<T> B00 = B.subTriMatrix(0,k);
            MatrixView<T> B10 = B.subMatrix(k,N,0,k);
            LowerTriMatrixView<T> B11 = B.subTriMatrix(k,N);

            DoTriLDivEq(A00,B00);
            B10 -= A10 * B00;
            TriLDivEq(A11,B10);
            TriLDivEq(A11,B11);
        }
    }

    template <class T, class Ta> 
    void TriLDivEq(
        const GenLowerTriMatrix<Ta>& A, LowerTriMatrixView<T> B)
    // B = A^-1 * B
    // where A is a triangle matrix
    {
#ifdef XDEBUG
        Matrix<Ta> A0(A);
        Matrix<T> B0(B);
#endif

        TMVAssert(A.size() == B.size());
        TMVAssert(!B.isunit() || A.isunit());

        if (B.isconj()) TriLDivEq(A.conjugate(),B.conjugate());
        else if (B.size() > 0) DoTriLDivEq(A,B);

#ifdef XDEBUG
        Matrix<T> BB = A0*B;
        if (Norm(BB-B0) > 0.001*Norm(A0)*Norm(B0)) {
            cerr<<"TriLDivEq: Lower/Lower\n";
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
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

#define InstFile "TMV_TriDiv_L.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


