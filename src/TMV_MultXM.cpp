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


#include "tmv/TMV_MatrixArithFunc.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_VectorArith.h"

#ifdef XDEBUG
#include <iostream>
#include "tmv/TMV_VIt.h"
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    //
    // MultXM
    //

    template <class T, class Ta> 
    static void rowMajorMultXM(
        const Ta alpha, MatrixView<T> A)
    {
        TMVAssert(A.isrm());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.colsize() > 0 && A.rowsize() > 0);
        TMVAssert(alpha != Ta(1));
        TMVAssert(alpha != Ta(0));

        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();

        T* Ai0 = A.ptr();

        for(ptrdiff_t i=M;i>0;--i,Ai0+=A.stepi()) {
            // A.row(i) *= alpha;
            T* Aij = Ai0;
            for(ptrdiff_t j=N;j>0;--j,++Aij) {
#ifdef TMVFLDEBUG
                TMVAssert(Aij >= A._first);
                TMVAssert(Aij < A._last);
#endif
                *Aij *= alpha;
            }
        }
    }

    template <class T> 
    void MultXM(const T alpha, MatrixView<T> A)
    // A = alpha * A
    {
#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<T> A2 = A;
        for(ptrdiff_t i=0;i<A.colsize();i++)
            for(ptrdiff_t j=0;j<A.rowsize();j++)
                A2(i,j) *= alpha;
        //cout<<"MultXM: alpha = "<<alpha<<", A = "<<TMV_Text(A)<<"  "<<A<<endl;
#endif

        if (A.colsize() > 0 && A.rowsize() > 0 && alpha != T(1)) {
            if (A.isconj()) MultXM(TMV_CONJ(alpha),A.conjugate());
            else if (alpha == T(0)) A.setZero();
            else if (A.canLinearize()) A.linearView() *= alpha;
            else if (A.isrm())
                if (isComplex(T()) && TMV_IMAG(alpha) == TMV_RealType(T)(0))
                    rowMajorMultXM(TMV_REAL(alpha),A);
                else
                    rowMajorMultXM(alpha,A);
            else if (A.iscm())
                if (isComplex(T()) && TMV_IMAG(alpha) == TMV_RealType(T)(0))
                    rowMajorMultXM(TMV_REAL(alpha),A.transpose());
                else 
                    rowMajorMultXM(alpha,A.transpose());
            else {
                const ptrdiff_t M = A.colsize();
                const ptrdiff_t N = A.rowsize();
                if (M < N)
                    for(ptrdiff_t i=0;i<M;++i) A.row(i) *= alpha;
                else 
                    for(ptrdiff_t j=0;j<N;++j) A.col(j) *= alpha;
            }
        }
#ifdef XDEBUG
        //cout<<"Done: A = "<<A<<endl;
        // Doing A-A2 becomes recursive call to MultXM
        Matrix<T> diff(A.colsize(),A.rowsize());
        for(ptrdiff_t i=0;i<A.colsize();i++)
            for(ptrdiff_t j=0;j<A.rowsize();j++)
                diff(i,j) = A(i,j) - A2(i,j);
        if (!(Norm(diff) <= 0.001*TMV_ABS(alpha)*Norm(A))) {
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
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == C.rowsize());
        TMVAssert(B.colsize() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        if (A.canLinearize() && B.canLinearize() && C.canLinearize() &&
            A.stepi() == C.stepi() && A.stepj() == C.stepj() &&
            B.stepi() == C.stepi() && B.stepj() == C.stepj()) {
            TMVAssert(A.stepi() == C.stepi() && A.stepj() == C.stepj());
            TMVAssert(B.stepi() == C.stepi() && B.stepj() == C.stepj());
            ElemMultVV<add>(
                alpha,A.constLinearView(),B.constLinearView(),C.linearView());
        } else if (C.isrm()) {
            const ptrdiff_t M = C.colsize();
            for(ptrdiff_t i=0;i<M;i++)
                ElemMultVV<add>(alpha,A.row(i),B.row(i),C.row(i));
        } else {
            const ptrdiff_t N = C.rowsize();
            for(ptrdiff_t j=0;j<N;j++)
                ElemMultVV<add>(alpha,A.col(j),B.col(j),C.col(j));
        }
    }

#define InstFile "TMV_MultXM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


