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
#include "tmv/TMV_MatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_VIt.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    //
    // AddMM
    //

    template <bool rm, bool a1, bool ca, class T, class T1, class T2> 
    static void DoRowAddMM(
        const T1 alpha, const GenMatrix<T2>& A, MatrixView<T> B)
    {
        TMVAssert(A.colsize() == B.colsize());
        TMVAssert(A.rowsize() == B.rowsize());
        TMVAssert(alpha != T1(0));
        TMVAssert(B.colsize() > 0);
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.ct() == NonConj);
        TMVAssert(!SameStorage(A,B));
        TMVAssert(rm == (A.isrm() && B.isrm()));
        TMVAssert(a1 == (alpha == T1(1)));
        TMVAssert(ca == A.isconj());

        const T2* Arowi = A.cptr();
        T* Browi = B.ptr();
        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();
        const ptrdiff_t Asi = A.stepi();
        const ptrdiff_t Asj = (rm ? 1 : A.stepj());
        const ptrdiff_t Bsi = B.stepi();
        const ptrdiff_t Bsj = (rm ? 1 : B.stepj());

        for(ptrdiff_t i=M;i>0;--i,Arowi+=Asi,Browi+=Bsi) {
            const T2* Aij = Arowi;
            T* Bij = Browi;
            for(ptrdiff_t j=N;j>0;--j,(rm?++Aij:Aij+=Asj),(rm?++Bij:Bij+=Bsj)) {
#ifdef TMVFLDEBUG
                TMVAssert(Bij >= B._first);
                TMVAssert(Bij < B._last);
#endif
                if (a1) *Bij += (ca ? TMV_CONJ(*Aij) : *Aij);
                else *Bij += alpha * (ca ? TMV_CONJ(*Aij) : *Aij);
            }
        }
    }

    template <bool rm, class T, class Ta> 
    static void rowAddMM(
        const T alpha, const GenMatrix<Ta>& A, MatrixView<T> B)
    { 
        if (TMV_IMAG(alpha) == TMV_RealType(T)(0))
            if (TMV_REAL(alpha) == TMV_RealType(T)(1))
                if (A.isconj()) DoRowAddMM<rm,true,true>(TMV_REAL(alpha),A,B);
                else DoRowAddMM<rm,true,false>(TMV_REAL(alpha),A,B);
            else
                if (A.isconj()) DoRowAddMM<rm,false,true>(TMV_REAL(alpha),A,B);
                else DoRowAddMM<rm,false,false>(TMV_REAL(alpha),A,B);
        else
            if (A.isconj()) DoRowAddMM<rm,false,true>(alpha,A,B);
            else DoRowAddMM<rm,false,false>(alpha,A,B);
    }

    template <class T, class Ta> 
    static void DoAddMM(
        const T alpha, const GenMatrix<Ta>& A, MatrixView<T> B)
    { 
        TMVAssert(A.colsize() == B.colsize());
        TMVAssert(A.rowsize() == B.rowsize());
        TMVAssert(alpha != T(0));
        TMVAssert(B.colsize() > 0);
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.ct() == NonConj);
        TMVAssert(!SameStorage(A,B));

        //cout<<"!SameStorage\n";
        //cout<<"A = "<<A.cptr()<<" "<<A.stepi()<<"  "<<A.stepj()<<std::endl;;
        //cout<<"B = "<<B.cptr()<<" "<<B.stepi()<<"  "<<B.stepj()<<std::endl;;
        //cout<<"A.isrm = "<<A.isrm()<<std::endl;
        //cout<<"B.isrm = "<<A.isrm()<<std::endl;
        //cout<<"A.iscm = "<<A.iscm()<<std::endl;
        //cout<<"B.iscm = "<<A.iscm()<<std::endl;
        if (A.canLinearize() && B.canLinearize() &&
            A.stepi() == B.stepi() && A.stepj() == B.stepj()) {
            //cout<<"linearize\n";
            B.linearView() += alpha * A.constLinearView();
        } else {
            if (A.isrm() && B.isrm())
                rowAddMM<true>(alpha,A,B); 
            else if (A.iscm() && B.iscm())
                rowAddMM<true>(alpha,A.transpose(),B.transpose()); 
            else if (A.rowsize() > A.colsize())
                rowAddMM<false>(alpha,A,B); 
            else
                rowAddMM<false>(alpha,A.transpose(),B.transpose()); 
        }
    }

    template <class T, class Ta> 
    void AddMM(const T alpha, const GenMatrix<Ta>& A, MatrixView<T> B)
    {
        TMVAssert(A.colsize() == B.colsize());
        TMVAssert(A.rowsize() == B.rowsize());
#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<T> B0 = B;
        Matrix<T> B2 = B;
        for(ptrdiff_t i=0;i<A.colsize();i++)
            for(ptrdiff_t j=0;j<A.rowsize();j++)
                B2(i,j) += alpha*A(i,j);
        cout<<"AddMM: alpha = "<<alpha<<", A = "<<TMV_Text(A)<<"  "<<A<<endl;
        cout<<", B = "<<TMV_Text(B)<<"  "<<B<<endl;
#endif

        if (alpha != T(0) && B.colsize() > 0 && B.rowsize() > 0) {
            if (B.isconj()) 
                AddMM(TMV_CONJ(alpha),A.conjugate(),B.conjugate());
            else {
                if (SameStorage(A,B)) {
                    if (B.isrm()) {
                        //cout<<"SameStorage, B isrm\n";
                        Matrix<T,RowMajor> A2 = A;
                        DoAddMM(alpha,A2,B);
                    } else {
                        //cout<<"SameStorage, B iscm\n";
                        Matrix<T,ColMajor> A2 = A;
                        DoAddMM(alpha,A2,B);
                    }
                } 
                else DoAddMM(alpha,A,B);
            }
        }
#ifdef XDEBUG
        cout<<"done: B = "<<B<<endl;
        Matrix<T> diff(B.colsize(),B.rowsize());
        for(ptrdiff_t i=0;i<B.colsize();i++)
            for(ptrdiff_t j=0;j<B.rowsize();j++)
                diff(i,j) = B(i,j) - B2(i,j);
        if (Norm(diff) > 0.001*(TMV_ABS(alpha)*Norm(A0)+Norm(B0))) {
            cerr<<"AddMM: alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"->B = "<<B<<endl;
            cerr<<"B2 = "<<B2<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta, class Tb> 
    void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C)
    {
        TMVAssert(A.colsize() == B.colsize());
        TMVAssert(A.rowsize() == B.rowsize());
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == C.rowsize());
#ifdef XDEBUG
        Matrix<T> A0 = A;
        Matrix<T> B0 = B;
        Matrix<T> C2(A.colsize(),A.rowsize());
        for(ptrdiff_t i=0;i<A.colsize();i++)
            for(ptrdiff_t j=0;j<A.rowsize();j++)
                C2(i,j) = alpha*A(i,j) + beta*B(i,j);
        cout<<"AddMM: alpha = "<<alpha<<", A = "<<TMV_Text(A)<<"  "<<A<<endl;
        cout<<"beta = "<<beta<<", B = "<<TMV_Text(B)<<"  "<<B;
        cout<<", C = "<<TMV_Text(C)<<"  "<<C<<endl;
#endif

        if (C.colsize() > 0 && C.rowsize() > 0) {
            if (SameStorage(A,C)) {
                if (SameStorage(B,C)) {
                    if (A.isrm()) {
                        Matrix<Ta,RowMajor> tempA = A;
                        C = B;
                        C *= beta;
                        AddMM(alpha,tempA,C);
                    } else {
                        Matrix<Ta,ColMajor> tempA = A;
                        C = B;
                        C *= beta;
                        AddMM(alpha,tempA,C);
                    }
                } else {
                    C = A;
                    C *= alpha;
                    AddMM(beta,B,C);
                }
            } else {
                C = B;
                C *= beta;
                AddMM(alpha,A,C);
            }
        }

#ifdef XDEBUG
        cout<<"Done: C = "<<C<<endl;
        Matrix<T> diff(C.colsize(),C.rowsize());
        for(ptrdiff_t i=0;i<C.colsize();i++)
            for(ptrdiff_t j=0;j<C.rowsize();j++)
                diff(i,j) = C(i,j) - C2(i,j);
        if (Norm(diff) >
            0.001*(1.+TMV_ABS(alpha)*Norm(A0)+TMV_ABS(beta)*Norm(B0))) {
            cerr<<"AddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"->C = "<<TMV_Text(C)<<"  "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_AddMM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


