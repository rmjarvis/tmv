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
#include "tmv/TMV_BandMatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    //
    // AddMM
    //

    template <class T, class Ta> 
    static void DoAddMM(
        const T alpha, const GenBandMatrix<Ta>& A, BandMatrixView<T> B)
    { 
        TMVAssert(A.colsize() == B.colsize());
        TMVAssert(A.rowsize() == B.rowsize());
        TMVAssert(B.nlo() >= A.nlo());
        TMVAssert(B.nhi() >= A.nhi());
        TMVAssert(alpha != T(0));
        TMVAssert(B.colsize() > 0);
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.ct() == NonConj);
        TMVAssert(!SameStorage(A,B));

        if (A.nlo() == B.nlo() && A.nhi() == B.nhi() && 
            A.canLinearize() && B.canLinearize() &&
            A.stepi() == B.stepi() && A.stepj() == B.stepj()) {
            B.linearView() += alpha*A.constLinearView();
        } else {
            for(ptrdiff_t i=-A.nlo();i<=A.nhi();++i) {
                B.diag(i) += alpha * A.diag(i);
            }
        }
    }

    template <class T, class Ta> 
    void AddMM(
        const T alpha, const GenBandMatrix<Ta>& A, BandMatrixView<T> B)
        // B += alpha * A
    {
#ifdef XDEBUG
        //cout<<"Band AddMM: alpha = "<<alpha<<endl;
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<"B = "<<TMV_Text(B)<<"  "<<B<<endl;
        Matrix<Ta> A0 = A;
        Matrix<T> B0 = B;
        Matrix<T> B2 = B0 + alpha*A0;
#endif
        TMVAssert(A.colsize() == B.colsize());
        TMVAssert(A.rowsize() == B.rowsize());
        TMVAssert(B.nlo() >= A.nlo());
        TMVAssert(B.nhi() >= A.nhi());

        if (B.colsize() > 0 && B.rowsize() > 0) {
            if (B.isconj())
                AddMM(TMV_CONJ(alpha),A.conjugate(),B.conjugate());
            else {
                if (SameStorage(A,B)) {
                    if (B.isrm()) {
                        BandMatrix<Ta,RowMajor> A2 = A;
                        DoAddMM(alpha,A2,B);
                    } else if (B.iscm()) {
                        BandMatrix<Ta,ColMajor> A2 = A;
                        DoAddMM(alpha,A2,B);
                    } else {
                        BandMatrix<Ta,DiagMajor> A2 = A;
                        DoAddMM(alpha,A2,B);
                    }
                } 
                else DoAddMM(alpha,A,B);
            }
        }
#ifdef XDEBUG
        if (Norm(B2-Matrix<T>(B)) > 0.001*TMV_ABS(alpha)*Norm(A0)*Norm(B0)) {
            cerr<<"Band AddMM\n";
            cerr<<"alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"   "<<A<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B.cptr()<<"   "<<B0<<endl;
            cerr<<"B => "<<B<<endl;
            cerr<<"B2 = "<<B2<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta, class Tb> 
    void AddMM(
        const T alpha, const GenBandMatrix<Ta>& A,
        const T beta, const GenBandMatrix<Tb>& B, BandMatrixView<T> C)
    {
#ifdef XDEBUG
        //cout<<"Band AddMM: alpha = "<<alpha<<", beta = "<<beta<<endl;
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<"B = "<<TMV_Text(B)<<"  "<<B<<endl;
        //cout<<"C = "<<TMV_Text(C)<<endl;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C2 = alpha*A0+beta*B0;
#endif
        TMVAssert(A.colsize() == B.colsize());
        TMVAssert(A.rowsize() == B.rowsize());
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == C.rowsize());
        TMVAssert(C.nlo() >= A.nlo());
        TMVAssert(C.nhi() >= A.nhi());
        TMVAssert(C.nlo() >= B.nlo());
        TMVAssert(C.nhi() >= B.nhi());

        if (C.isconj()) {
            AddMM(TMV_CONJ(alpha),A.conjugate(),TMV_CONJ(beta),B.conjugate(),C.conjugate());
        }
        else {
            if (B.colsize() > 0 && B.rowsize() > 0) {
                if (SameStorage(A,C)) {
                    if (SameStorage(B,C)) {
                        if (B.isrm()) {
                            BandMatrix<T,RowMajor> tempB = B;
                            C = alpha*A;
                            DoAddMM(beta,tempB,C);
                        } else if (C.iscm()) {
                            BandMatrix<T,ColMajor> tempB = B;
                            C = alpha*A;
                            DoAddMM(beta,tempB,C);
                        } else {
                            BandMatrix<T,DiagMajor> tempB = B;
                            C = alpha*A;
                            DoAddMM(beta,tempB,C);
                        }
                    } else {
                        C = alpha*A;
                        DoAddMM(beta,B,C);
                    }
                } else {
                    C = beta*B;
                    DoAddMM(alpha,A,C);
                }
            }
        }
        //cout<<"Done C => "<<C<<endl;
#ifdef XDEBUG
        if (Norm(C2-Matrix<T>(C)) > 
            0.001*(TMV_ABS(alpha)*Norm(A0)+TMV_ABS(beta)*Norm(B0))) {
            cerr<<"Band AddMM\n";
            cerr<<"alpha,beta = "<<alpha<<","<<beta<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"   "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B.cptr()<<"   "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C.cptr()<<" ->  "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta, class Tb> 
    void AddMM(
        const T alpha, const GenBandMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C)
    {
#ifdef XDEBUG
        //cout<<"Band AddMM: alpha = "<<alpha<<", beta = "<<beta<<endl;
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"  "<<A<<endl;
        //cout<<"B = "<<TMV_Text(B)<<"  "<<B.cptr()<<"  "<<B<<endl;
        //cout<<"C = "<<TMV_Text(C)<<"  "<<C.cptr()<<endl;
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C2 = alpha*A0+beta*B0;
#endif
        TMVAssert(A.colsize() == B.colsize());
        TMVAssert(A.rowsize() == B.rowsize());
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == C.rowsize());

        if (C.isconj()) {
            AddMM(TMV_CONJ(alpha),A.conjugate(),beta,B.conjugate(),
                  C.conjugate());
        } else {
            if (C.colsize() > 0 && C.rowsize() > 0) {
                if (SameStorage(A,C)) {
                    if (SameStorage(B,C)) {
                        if (A.isrm()) {
                            BandMatrix<Ta,RowMajor> tempA = A;
                            C = beta*B;
                            DoAddMM(alpha,tempA,
                                    BandMatrixView<T>(C,A.nlo(),A.nhi()));
                        } else if (A.iscm()) {
                            BandMatrix<Ta,ColMajor> tempA = A;
                            C = beta*B;
                            DoAddMM(alpha,tempA,
                                    BandMatrixView<T>(C,A.nlo(),A.nhi()));
                        } else {
                            BandMatrix<Ta,DiagMajor> tempA = A;
                            C = beta*B;
                            DoAddMM(alpha,tempA,
                                    BandMatrixView<T>(C,A.nlo(),A.nhi()));
                        }
                    } else {
                        C = alpha*A;
                        AddMM(beta,B,C);
                    }
                } else {
                    C = beta*B;
                    DoAddMM(alpha,A,BandMatrixView<T>(C,A.nlo(),A.nhi()));
                }
            }
        }
#ifdef XDEBUG
        //cout<<"Done C => "<<C<<endl;
        if (Norm(C2-Matrix<T>(C)) > 
            0.001*(TMV_ABS(alpha)*Norm(A0)+TMV_ABS(beta)*Norm(B0))) {
            cerr<<"Band AddMM\n";
            cerr<<"alpha = "<<alpha<<","<<beta<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"   "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B.cptr()<<"   "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C.cptr()<<" ->  "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_AddBB.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


