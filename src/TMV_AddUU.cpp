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

    //
    // AddMM
    //

    template <bool rm, bool a1, bool ca, class T1, class T2, class T3> 
    static void DoRowAddMM(
        const T1 alpha, const GenUpperTriMatrix<T2>& A, 
        UpperTriMatrixView<T3> B)
    {
        TMVAssert(!A.isunit());
        TMVAssert(!B.isunit());
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() > 0);
        TMVAssert(alpha != T1(0));
        TMVAssert(B.ct() == NonConj);

        const ptrdiff_t N = A.size();
        const ptrdiff_t Astepj = rm ? 1 : A.stepj();
        const ptrdiff_t Ads = A.stepi() + Astepj;
        const ptrdiff_t Bstepj = rm ? 1 : B.stepj();
        const ptrdiff_t Bds = B.stepi() + Bstepj;
        const T2* Aii = A.cptr();
        T3* Bii = B.ptr();

        for(ptrdiff_t len=N;len>0;--len,Aii+=Ads,Bii+=Bds) {
            const T2* Aij = Aii;
            T3* Bij = Bii;
            for(ptrdiff_t j=len;j>0;--j,(rm?++Aij:Aij+=Astepj),
                (rm?++Bij:Bij+=Bstepj)) {
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
    static inline void RowAddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        UpperTriMatrixView<T> B)
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

    template <bool cm, bool a1, bool ca, class T1, class T2, class T3> 
    static void DoColAddMM(
        const T1 alpha, const GenUpperTriMatrix<T2>& A, 
        UpperTriMatrixView<T3> B)
    {
        TMVAssert(!A.isunit());
        TMVAssert(!B.isunit());
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() > 0);
        TMVAssert(alpha != T1(0));
        TMVAssert(B.ct() == NonConj);

        const ptrdiff_t N = A.size();
        const ptrdiff_t Astepi = (cm ? 1 : A.stepi());
        const ptrdiff_t Astepj = A.stepj();
        const ptrdiff_t Bstepi = (cm ? 1 : B.stepi());
        const ptrdiff_t Bstepj = B.stepj();
        const T2* Acolj = A.cptr()+(N-1)*Astepj;
        T3* Bcolj = B.ptr()+(N-1)*Bstepj;

        for(ptrdiff_t j=N;j>0;--j,Acolj-=Astepj,Bcolj-=Bstepj) {
            const T2* Aij = Acolj;
            T3* Bij = Bcolj;
            for(ptrdiff_t i=j;i>0;--i,(cm?++Aij:Aij+=Astepi),(cm?++Bij:Bij+=Bstepi)) {
                if (a1) *Bij += (ca ? TMV_CONJ(*Aij) : *Aij);
                else *Bij += alpha * (ca ? TMV_CONJ(*Aij) : *Aij);
            }
        }
    }

    template <bool cm, class T, class Ta> 
    static inline void ColAddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        UpperTriMatrixView<T> B)
    {
        if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) 
            if (TMV_REAL(alpha) == TMV_RealType(T)(1))
                if (A.isconj()) DoColAddMM<cm,true,true>(TMV_REAL(alpha),A,B); 
                else DoColAddMM<cm,true,false>(TMV_REAL(alpha),A,B);
            else
                if (A.isconj()) DoColAddMM<cm,false,true>(TMV_REAL(alpha),A,B); 
                else DoColAddMM<cm,false,false>(TMV_REAL(alpha),A,B);
        else
            if (A.isconj()) DoColAddMM<cm,false,true>(alpha,A,B); 
            else DoColAddMM<cm,false,false>(alpha,A,B);
    }

    template <class T, class Ta> 
    static inline void DoAddMM(const T alpha,
              const GenUpperTriMatrix<Ta>& A, UpperTriMatrixView<T> B)
    { 
        TMVAssert(!B.isunit());
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(B.ct() == NonConj);

        if (A.isunit()) {
            B.diag().addToAll(alpha);
            if (A.size() > 1)
                DoAddMM(alpha,A.offDiag(),B.offDiag());
        } else {
            if (A.isrm() && B.isrm()) RowAddMM<true>(alpha,A,B); 
            else if (A.iscm() && B.iscm()) ColAddMM<true>(alpha,A,B); 
            else RowAddMM<false>(alpha,A,B); 
        }
    }

    template <class T, class Ta> 
    void AddMM(const T alpha,
              const GenUpperTriMatrix<Ta>& A, UpperTriMatrixView<T> B)
        // B += alpha * A
    {
        TMVAssert(!B.isunit());
        TMVAssert(A.size() == B.size());
#ifdef XDEBUG
        Matrix<Ta> A0 = A;
        Matrix<T> B0 = B;
        Matrix<T> B2 = alpha*A0 + B0;
        //cout<<"AddMM: alpha = "<<alpha<<", A = "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<", B = "<<TMV_Text(B)<<"  "<<B<<endl;
#endif

        if (A.size() > 0) {
            if (B.isconj()) AddMM(TMV_CONJ(alpha),A.conjugate(),B.conjugate());
            else {
                if (SameStorage(A,B)) {
                    if (B.isrm()) {
                        UpperTriMatrix<Ta,NonUnitDiag|RowMajor> tempA = alpha*A;
                        DoAddMM(T(1),tempA,B);
                    } else {
                        UpperTriMatrix<Ta,NonUnitDiag|ColMajor> tempA = alpha*A;
                        DoAddMM(T(1),tempA,B);
                    }
                } 
                else {
                    DoAddMM(alpha,A,B);
                }
            }
        }
#ifdef XDEBUG
        //cout<<"done: B = "<<B<<endl;
        if (Norm(Matrix<T>(B)-B2) > 0.001*Norm(B)) {
            cerr<<"AddMM: alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"->B = "<<B<<endl;
            cerr<<"B2 = "<<B2<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta, class Tb> 
    void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        const T beta, const GenUpperTriMatrix<Tb>& B,
        UpperTriMatrixView<T> C)
    {
        TMVAssert(A.size() == B.size());
        TMVAssert(A.size() == C.size());
#ifdef XDEBUG
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C2 = alpha*A0 + beta*B0;
        //cout<<"AddMM: alpha = "<<alpha<<", A = "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<"beta = "<<beta<<", B = "<<TMV_Text(B)<<"  "<<B;
        //cout<<", C = "<<TMV_Text(C)<<"  "<<C<<endl;
#endif

        if (C.size() == 0) return;
        if (A.isunit()) {
            if (B.isunit()) {
                if (A.size() > 1)
                    AddMM(alpha,A.offDiag(),beta,B.offDiag(),C.offDiag());
                C.diag().setAllTo(alpha+beta);
            } else {
                if (A.size() > 1)
                    AddMM(alpha,A.offDiag(),beta,B.offDiag(),C.offDiag());
                C.diag() = beta * B.diag();
                C.diag().addToAll(alpha);
            }
        } else {
            if (B.isunit()) {
                if (A.size() > 1)
                    AddMM(alpha,A.offDiag(),beta,B.offDiag(),C.offDiag());
                C.diag() = alpha * A.diag();
                C.diag().addToAll(beta);
            } else {
                if (SameStorage(A,C)) {
                    if (SameStorage(B,C)) {
                        if (A.isunit()) {
                            if (A.isrm()) {
                                UpperTriMatrix<Ta,UnitDiag|RowMajor> tempA 
                                    = alpha*A;
                                C = beta*B;
                                AddMM(T(1),tempA,C);
                            } else {
                                UpperTriMatrix<Ta,UnitDiag|ColMajor> tempA 
                                    = alpha*A;
                                C = beta*B;
                                AddMM(T(1),tempA,C);
                            }
                        } else {
                            if (A.isrm()) {
                                UpperTriMatrix<Ta,NonUnitDiag|RowMajor> tempA 
                                    = alpha*A;
                                C = beta*B;
                                AddMM(T(1),tempA,C);
                            } else {
                                UpperTriMatrix<Ta,NonUnitDiag|ColMajor> tempA 
                                    = alpha*A;
                                C = beta*B;
                                AddMM(T(1),tempA,C);
                            }
                        }
                    } else {
                        C = alpha*A;
                        AddMM(beta,B,C);
                    }
                } else {
                    C = beta*B;
                    AddMM(alpha,A,C);
                }
            }
        }

#ifdef XDEBUG
        if (Norm(Matrix<T>(C)-C2) > 0.001*
            (TMV_ABS(alpha)*Norm(A0)+TMV_ABS(beta)*Norm(B0))) {
            cerr<<"AddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"->C = "<<TMV_Text(C)<<"  "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta, class Tb> 
    void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == B.rowsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == C.rowsize());
#ifdef XDEBUG
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C2 = alpha*A0 + beta*B0;
        //cout<<"AddMM: alpha = "<<alpha<<", A = "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<"beta = "<<beta<<", B = "<<TMV_Text(B)<<"  "<<B;
        //cout<<", C = "<<TMV_Text(C)<<"  "<<C<<endl;
#endif

        if (A.size() > 0) {
            if (SameStorage(A,C)) {
                if (SameStorage(B,C)) {
                    if (C.isrm()) {
                        Matrix<T,RowMajor> temp = beta*B;
                        AddMM(alpha,A,temp.view());
                        C = temp;
                    } else {
                        Matrix<T,ColMajor> temp = beta*B;
                        AddMM(alpha,A,temp.view());
                        C = temp;
                    }
                } else {
                    C.upperTri() = alpha*A;
                    AddMM(beta,B,C);
                }
            } else {
                C = beta*B;
                AddMM(alpha,A,C.upperTri());
            }
        }

#ifdef XDEBUG
        if (Norm(Matrix<T>(C)-C2) > 0.001*
            (TMV_ABS(alpha)*Norm(A0)+TMV_ABS(beta)*Norm(B0))) {
            cerr<<"AddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"->C = "<<TMV_Text(C)<<"  "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta, class Tb> 
    void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        const T beta, const GenLowerTriMatrix<Tb>& B, MatrixView<T> C)
        // C = alpha * A + beta * B
    {
        TMVAssert(A.size() == B.size());
        TMVAssert(C.rowsize() == A.size());
        TMVAssert(C.colsize() == A.size());
#ifdef XDEBUG
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C2 = alpha*A0 + beta*B0;
        //cout<<"AddMM: alpha = "<<alpha<<", A = "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<"beta = "<<beta<<", B = "<<TMV_Text(B)<<"  "<<B;
        //cout<<", C = "<<TMV_Text(C)<<"  "<<C<<endl;
#endif

        if (A.size() > 0) {
            if (SameStorage(A,C)) {
                if (SameStorage(B,C)) {
                    if (C.isrm()) {
                        Matrix<T,RowMajor> temp = beta*B;
                        AddMM(alpha,A,temp.view());
                        C = temp;
                    } else {
                        Matrix<T,ColMajor> temp = beta*B;
                        AddMM(alpha,A,temp.view());
                        C = temp;
                    }
                } else {
                    C = alpha*A;
                    AddMM(beta,B,C);
                }
            } else {
                C = beta*B;
                AddMM(alpha,A,C);
            }
        }

#ifdef XDEBUG
        if (Norm(Matrix<T>(C)-C2) > 0.001*
            (TMV_ABS(alpha)*Norm(A0)+TMV_ABS(beta)*Norm(B0))) {
            cerr<<"AddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"->C = "<<TMV_Text(C)<<"  "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_AddUU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


