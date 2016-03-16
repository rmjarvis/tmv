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

//#define TMV_DEBUG
//#include <iostream>

#include "tmv/TMV_TriMatrixArithFunc.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_VectorArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    //
    // MultXM
    //

    template <class T, class T1> 
    static void RowMajorMultXM(const T1 alpha, UpperTriMatrixView<T> A)
    {
        TMVAssert(A.isrm());
        TMVAssert(!A.isunit());
        TMVAssert(alpha != T1(1));
        TMVAssert(A.size() > 0);
        TMVAssert(A.ct() == NonConj);

        T* Aii = A.ptr();
        const ptrdiff_t ds = A.stepi()+1;
        const ptrdiff_t N = A.size();

        for(ptrdiff_t len=N;len>0;--len,Aii+=ds) {
            // A.row(i,i,N) *= alpha;
            T* Aij = Aii;
            for(ptrdiff_t j=len;j>0;--j,++Aij) {
#ifdef TMVFLDEBUG
                TMVAssert(Aij >= A._first);
                TMVAssert(Aij < A._last);
#endif
                *Aij *= alpha;
            }
        }
    }

    template <class T, class T1> 
    static void ColMajorMultXM(const T1 alpha, UpperTriMatrixView<T> A)
    {
        TMVAssert(A.iscm());
        TMVAssert(!A.isunit());
        TMVAssert(alpha != T1(1));
        TMVAssert(A.size() > 0);
        TMVAssert(A.ct() == NonConj);

        T* A0j = A.ptr();
        const ptrdiff_t Astepj = A.stepj();
        const ptrdiff_t N = A.size();

        for(ptrdiff_t j=N,len=1;j>0;--j,++len,A0j+=Astepj) {
            // A.col(j,0,j+1) *= alpha;
            T* Aij = A0j;
            for(ptrdiff_t i=len;i>0;--i,++Aij) {
#ifdef TMVFLDEBUG
                TMVAssert(Aij >= A._first);
                TMVAssert(Aij < A._last);
#endif
                *Aij *= alpha;
            }
        }
    }

    template <class T> 
    void MultXM(const T alpha, UpperTriMatrixView<T> A)
    // A = alpha * A
    {
#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<T> A2 = alpha * A0;
        //cout<<"MultXM: alpha = "<<alpha<<", A = "<<TMV_Text(A)<<" "<<A<<endl;
#endif

        if (A.size() > 0 && alpha != T(1)) {
            TMVAssert(!A.isunit());
            if (A.isconj()) {
                MultXM(TMV_CONJ(alpha),A.conjugate());
            } else if (alpha == T(0)) {
                A.setZero();
            } else if (A.isrm()) {
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0))
                    RowMajorMultXM(TMV_REAL(alpha),A);
                else
                    RowMajorMultXM(alpha,A);
            } else if (A.iscm()) {
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0))
                    ColMajorMultXM(TMV_REAL(alpha),A);
                else
                    ColMajorMultXM(alpha,A);
            } else {
                const ptrdiff_t M = A.colsize();
                const ptrdiff_t N = A.rowsize();
                for(ptrdiff_t i=0;i<M;++i) 
                    A.row(i,i,N) *= alpha;
            }
        }
#ifdef XDEBUG
        //cout<<"Done MultXM: A = "<<A<<endl;
        if (!(Norm(Matrix<T>(A)-A2) <= 0.001*TMV_ABS(alpha)*Norm(A))) {
            cerr<<"MultXM: alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"-> A = "<<A<<endl;
            cerr<<"A2 = "<<A2<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    void ElemMultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<T> C)
    {
        //std::cout<<"Start ElemMultMM:\n";
        //std::cout<<"add = "<<add<<std::endl;
        //std::cout<<"alpha = "<<alpha<<std::endl;
        //std::cout<<"A = "<<TMV_Text(A)<<"  "<<A<<std::endl;
        //std::cout<<"B = "<<TMV_Text(B)<<"  "<<B<<std::endl;
        //std::cout<<"C = "<<TMV_Text(C)<<"  "<<C<<std::endl;
        TMVAssert(A.size() == C.size());
        TMVAssert(B.size() == C.size());
        if (C.isunit()) {
            TMVAssert(alpha == T(1));
            TMVAssert(A.isunit());
            TMVAssert(B.isunit());
            TMVAssert(add == false);
            if (C.size() > 1) 
                ElemMultMM<add>(alpha,A.offDiag(),B.offDiag(),C.offDiag());
        } else if (A.isunit()) {
            if (B.isunit()) {
                if (add) C.diag().addToAll(alpha);
                else C.diag().setAllTo(alpha);
            } else {
                if (add) C.diag() += alpha * B.diag();
                else C.diag() = alpha * B.diag();
            }
            if (C.size() > 1) 
                ElemMultMM<add>(alpha,A.offDiag(),B.offDiag(),C.offDiag());
        } else if (B.isunit()) {
            if (add) C.diag() += alpha * A.diag();
            else C.diag() = alpha * A.diag();
            if (C.size() > 1) 
                ElemMultMM<add>(alpha,A.offDiag(),B.offDiag(),C.offDiag());
        } else {
            const ptrdiff_t N = C.size();
            if (C.isrm()) {
                for(ptrdiff_t i=0;i<N;i++)
                    ElemMultVV<add>(
                        alpha,A.row(i,i,N),B.row(i,i,N),C.row(i,i,N));
            } else {
                for(ptrdiff_t j=0;j<N;j++)
                    ElemMultVV<add>(
                        alpha,A.col(j,0,j+1),B.col(j,0,j+1),C.col(j,0,j+1));
            }
        }
        //std::cout<<"C => "<<TMV_Text(C)<<"  "<<C<<std::endl;
    }

#define InstFile "TMV_MultXU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv

