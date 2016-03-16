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

#include "tmv/TMV_BandMatrixArithFunc.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_VectorArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_BandMatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

    //
    // MultXM
    //

    template <class T> 
    void MultXM(const T alpha, BandMatrixView<T> A)
    {
#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<T> A2 = alpha*A0;
#endif
        if (A.rowsize() > 0 && A.colsize() > 0 && alpha != T(1)) {
            if (A.isconj()) MultXM(TMV_CONJ(alpha),A.conjugate());
            else if (alpha == T(0)) A.setZero();
            else if (A.canLinearize()) A.linearView() *= alpha;
            else {
                for(ptrdiff_t i=-A.nlo();i<=A.nhi();++i) A.diag(i) *= alpha;
            }
        }
#ifdef XDEBUG
        if (!(Norm(A2-A) <= 0.001*TMV_ABS(alpha)*Norm(A0))) {
            cerr<<"MultXM: alpha = "<<alpha;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"A2 = "<<A2<<endl;
            cerr<<"Norm(diff) = "<<Norm(A-A2)<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    void ElemMultMM(
        const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
        BandMatrixView<T> C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == C.rowsize());
        TMVAssert(B.colsize() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        const ptrdiff_t lo = TMV_MIN(A.nlo(),B.nlo());
        const ptrdiff_t hi = TMV_MIN(A.nhi(),B.nhi());
        if (A.nlo() == lo && A.nhi() == hi && 
            B.nlo() == lo && B.nhi() == hi &&
            C.nlo() == lo && C.nhi() == hi) {
            if (A.canLinearize() && B.canLinearize() && C.canLinearize() &&
                A.stepi() == C.stepi() && A.stepj() == C.stepj() &&
                B.stepi() == C.stepi() && B.stepj() == C.stepj()) {
                ElemMultVV<add>(
                    alpha,A.constLinearView(),B.constLinearView(),
                    C.linearView());
            } else {
                for(ptrdiff_t i=-lo;i<=hi;++i) 
                    ElemMultVV<add>(alpha,A.diag(i),B.diag(i),C.diag(i));
            }
        } else {
            if (!add) {
                if (C.nlo() > lo) {
                    C.diagRange(-C.nlo(),-lo).setZero();
                }
                if (C.nhi() > hi) {
                    C.diagRange(hi+1,C.nhi()+1).setZero();
                }
            }
            ElemMultMM<add>(
                alpha,A.diagRange(-lo,hi+1),B.diagRange(-lo,hi+1),
                C.diagRange(-lo,hi+1));
        }
    }

#define InstFile "TMV_MultXB.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


