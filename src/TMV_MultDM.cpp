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


#include "tmv/TMV_DiagMatrixArithFunc.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_DiagMatrixArith.h"

#ifdef XDEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

    template <bool rm, bool ca, class T, class Ta> 
    static void RowMultEqMM(
        const GenDiagMatrix<Ta>& A, MatrixView<T> B)
    {
        // B = A * B
        // Bij = Ai * Bij
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(B.ct() == NonConj);
        TMVAssert(rm == B.isrm());
        TMVAssert(ca == A.diag().isconj());

        const ptrdiff_t M = B.colsize();
        const ptrdiff_t N = B.rowsize();

        const Ta* Ai = A.diag().cptr();
        T* Browi = B.ptr();
        const ptrdiff_t Astep = A.diag().step();
        const ptrdiff_t stepj = B.stepj();
        const ptrdiff_t stepi = B.stepi();

        for(ptrdiff_t i=M;i>0;--i,Ai+=Astep,Browi+=stepi) {
            T* Bij = Browi;
            if (*Ai == Ta(0)) {
                for(ptrdiff_t j=N;j>0;--j,(rm?++Bij:Bij+=stepj)) {
#ifdef TMVFLDEBUG
                    TMVAssert(Bij >= B._first);
                    TMVAssert(Bij < B._last);
#endif
                    *Bij = T(0);
                }
            } else if (TMV_IMAG(*Ai) == TMV_RealType(Ta)(0)) {
                for(ptrdiff_t j=N;j>0;--j,(rm?++Bij:Bij+=stepj)) {
#ifdef TMVFLDEBUG
                    TMVAssert(Bij >= B._first);
                    TMVAssert(Bij < B._last);
#endif
                    *Bij *= TMV_REAL(*Ai);
                }
            } else {
                for(ptrdiff_t j=N;j>0;--j,(rm?++Bij:Bij+=stepj)) {
#ifdef TMVFLDEBUG
                    TMVAssert(Bij >= B._first);
                    TMVAssert(Bij < B._last);
#endif
                    *Bij *= (ca?TMV_CONJ(*Ai):*Ai);
                }
            }
        }
    }

    template <bool cm, bool ca, class T, class Ta> 
    static void DoColMultEqMM(
        const GenDiagMatrix<Ta>& A, MatrixView<T> B)
    {
        // B = A * B 
        // Bij = Ai * Bij
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(B.ct()==NonConj);
        TMVAssert(A.diag().step() == 1);
        TMVAssert(cm == B.iscm());
        TMVAssert(ca == A.diag().isconj());

        const Ta*const Aptr = A.diag().cptr();
        T*const Bptr = B.ptr();
        const ptrdiff_t stepj = B.stepj();
        const ptrdiff_t stepi = B.stepi();
        const ptrdiff_t M = B.colsize();
        ptrdiff_t N = B.rowsize();

        T* Bcolj = Bptr;
        for(;N>0;--N,Bcolj+=stepj) {
            T* Bij = Bcolj;
            const Ta* Ai = Aptr;
            for(ptrdiff_t i=M;i>0;--i,++Ai,(cm?++Bij:Bij+=stepi)) {
#ifdef TMVFLDEBUG
                TMVAssert(Bij >= B._first);
                TMVAssert(Bij < B._last);
#endif
                *Bij *= (ca ? TMV_CONJ(*Ai) : *Ai);
            }
        }
    }

    template <bool cm, bool ca, class T, class Ta> 
    static inline void ColMultEqMM(
        const GenDiagMatrix<Ta>& A, MatrixView<T> B)
    {
        if (A.diag().step() == 1)
            DoColMultEqMM<cm,ca>(A,B);
        else {
            DiagMatrix<Ta> AA = A;
            DoColMultEqMM<cm,false>(AA,B);
        }
    }

    template <class T, class Ta> 
    void MultEqMM(
        const T alpha,
        const GenDiagMatrix<Ta>& A, MatrixView<T> B)
    // B = alpha * A * B
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(alpha != T(0));
#ifdef XDEBUG
        Matrix<T> B0 = B;
        Matrix<Ta> A0 = A;
        Matrix<T> B2 = alpha*A0*B0;
#endif

        if (B.isconj()) MultEqMM(TMV_CONJ(alpha),A.conjugate(),B.conjugate());
        else if (B.colsize() > 0 && B.rowsize() > 0) {
            if (alpha != T(1)) {
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                    DiagMatrix<Ta> AA = TMV_REAL(alpha) * A;
                    if (B.isrm()) RowMultEqMM<true,false>(AA,B);
                    else if (B.iscm()) DoColMultEqMM<true,false>(AA,B);
                    else if (B.colsize() > B.rowsize()) 
                        DoColMultEqMM<false,false>(AA,B);
                    else RowMultEqMM<false,false>(AA,B);
                }
                else {
                    // AA = alpha * A;
                    DiagMatrix<T> AA = alpha * A;
                    if (B.isrm()) RowMultEqMM<true,false>(AA,B);
                    else if (B.iscm()) DoColMultEqMM<true,false>(AA,B);
                    else if (B.colsize() > B.rowsize()) 
                        DoColMultEqMM<false,false>(AA,B);
                    else RowMultEqMM<false,false>(AA,B);
                }
            }
            else if (A.diag().isconj()) {
                if (B.isrm()) RowMultEqMM<true,true>(A,B);
                else if (B.iscm()) ColMultEqMM<true,true>(A,B);
                else if (B.colsize() > B.rowsize()) 
                    ColMultEqMM<false,true>(A,B);
                else RowMultEqMM<false,true>(A,B);
            }
            else {
                if (B.isrm()) RowMultEqMM<true,false>(A,B);
                else if (B.iscm()) ColMultEqMM<true,false>(A,B);
                else if (B.colsize() > B.rowsize()) 
                    ColMultEqMM<false,false>(A,B);
                else RowMultEqMM<false,false>(A,B);
            }
        }

#ifdef XDEBUG
        if (!(Norm(Matrix<T>(B)-B2) <= 
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)))) {
            cerr<<"MultEqMM: alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<" step "<<A.diag().step()<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"-> B = "<<B<<endl;
            cerr<<"B2 = "<<B2<<endl;
            abort();
        }
#endif
    }

    template <bool rm, bool ca, bool cb, class T, class Ta, class Tb>
    static void DoRowAddMultMM(
        const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        // C += A * B
        // Cij += Ai * Bij
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.rowsize() > 0);
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.ct() == NonConj);
        TMVAssert(rm == (B.isrm() && C.isrm()));
        TMVAssert(ca == A.diag().isconj());
        TMVAssert(cb == B.isconj());

        const Ta* Ai = A.diag().cptr();
        const Tb* Browi = B.cptr();
        T* Crowi = C.ptr();
        const ptrdiff_t Astep = A.diag().step();
        const ptrdiff_t Bstepj = B.stepj();
        const ptrdiff_t Bstepi = B.stepi();
        const ptrdiff_t Cstepj = C.stepj();
        const ptrdiff_t Cstepi = C.stepi();
        const ptrdiff_t M = C.colsize();
        const ptrdiff_t N = C.rowsize();

        for(ptrdiff_t i=M;i>0;--i,Ai+=Astep,Browi+=Bstepi,Crowi+=Cstepi) {
            const Tb* Bij = Browi;
            T* Cij = Crowi;
            if (TMV_IMAG(*Ai) == TMV_RealType(Ta)(0)) {
                if (TMV_REAL(*Ai) != TMV_RealType(Ta)(0))
                    for(ptrdiff_t j=N;j>0;--j,(rm?++Bij:Bij+=Bstepj),
                        (rm?++Cij:Cij+=Cstepj)) {
#ifdef TMVFLDEBUG
                        TMVAssert(Cij >= C._first);
                        TMVAssert(Cij < C._last);
#endif
                        *Cij += TMV_REAL(*Ai)*(cb?TMV_CONJ(*Bij):*Bij);
                    }
            }
            else
                for(ptrdiff_t j=N;j>0;--j,(rm?++Bij:Bij+=Bstepj),
                    (rm?++Cij:Cij+=Cstepj)) {
#ifdef TMVFLDEBUG
                    TMVAssert(Cij >= C._first);
                    TMVAssert(Cij < C._last);
#endif
                    *Cij += (ca?TMV_CONJ(*Ai):*Ai)*(cb?TMV_CONJ(*Bij):*Bij);
                }
        }
    }

    template <bool rm, bool ca, class T, class Ta, class Tb>
    static inline void RowAddMultMM(
        const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        if (B.isconj()) DoRowAddMultMM<rm,ca,true>(A,B,C);
        else DoRowAddMultMM<rm,ca,false>(A,B,C);
    }

    template <bool cm, bool ca, bool cb, class T, class Ta, class Tb> 
    static void DoColAddMultMM(
        const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        // C += A * B 
        // Cij = Ai * Bij
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.rowsize() > 0);
        TMVAssert(C.colsize() > 0);
        TMVAssert(A.diag().step() == 1);
        TMVAssert(cm == (B.iscm() && C.iscm()));
        TMVAssert(ca == A.diag().isconj());
        TMVAssert(cb == B.isconj());

        const Ta*const Aptr = A.diag().cptr();
        const Tb* Bcolj = B.cptr();
        T* Ccolj = C.ptr();
        const ptrdiff_t Cstepj = C.stepj();
        const ptrdiff_t Cstepi = C.stepi();
        const ptrdiff_t Bstepj = B.stepj();
        const ptrdiff_t Bstepi = B.stepi();
        const ptrdiff_t M = C.colsize();
        const ptrdiff_t N = C.rowsize();

        for(ptrdiff_t j=N;j>0;--j,Bcolj+=Bstepj,Ccolj+=Cstepj) {
            const Tb* Bij = Bcolj;
            T* Cij = Ccolj;
            const Ta* Ai = Aptr;
            for(ptrdiff_t i=M;i>0;--i,++Ai,(cm?++Bij:Bij+=Bstepi),
                (cm?++Cij:Cij+=Cstepi)) {
#ifdef TMVFLDEBUG
                TMVAssert(Cij >= C._first);
                TMVAssert(Cij < C._last);
#endif
                *Cij += (ca ? TMV_CONJ(*Ai) : *Ai) * (cb ? TMV_CONJ(*Bij) : *Bij);
            }
        }
    }

    template <bool cm, bool ca, class T, class Ta, class Tb> 
    static void ColAddMultMM(
        const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    { 
        if (A.diag().step() == 1)
            if (B.isconj())
                DoColAddMultMM<cm,ca,true>(A,B,C);
            else
                DoColAddMultMM<cm,ca,false>(A,B,C);
        else {
            DiagMatrix<Ta> AA = A;
            if (B.isconj())
                DoColAddMultMM<cm,ca,true>(AA,B,C);
            else
                DoColAddMultMM<cm,ca,false>(AA,B,C);
        }
    }

    template <class T, class Ta, class Tb> 
    static void AddMultMM(const T alpha,
                          const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
                          MatrixView<T> C)
    // C += alpha * A * B
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(alpha != T(0));
#ifdef XDEBUG
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C0 = C;
        Matrix<T> C2 = C0+alpha*A0*B0;
#endif

        if (C.isconj()) {
            AddMultMM(
                TMV_CONJ(alpha),A.conjugate(),B.conjugate(),C.conjugate());
        } else if (C.colsize() > 0 && C.rowsize() > 0) {
            if (alpha != T(1)) {
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                    DiagMatrix<Ta> AA = TMV_REAL(alpha) * A;
                    if (B.isrm() && C.isrm()) 
                        RowAddMultMM<true,false>(AA,B,C);
                    else if (B.iscm() && C.iscm()) 
                        ColAddMultMM<true,false>(AA,B,C);
                    else if (B.colsize() > B.rowsize()) 
                        ColAddMultMM<false,false>(AA,B,C);
                    else 
                        RowAddMultMM<false,false>(AA,B,C);
                }
                else {
                    DiagMatrix<T> AA = alpha * A;
                    if (B.isrm() && C.isrm()) 
                        RowAddMultMM<true,false>(AA,B,C);
                    else if (B.iscm() && C.iscm()) 
                        ColAddMultMM<true,false>(AA,B,C);
                    else if (B.colsize() > B.rowsize()) 
                        ColAddMultMM<false,false>(AA,B,C);
                    else 
                        RowAddMultMM<false,false>(AA,B,C);
                }
            }
            else if (A.diag().isconj()) {
                if (B.isrm() && C.isrm()) 
                    RowAddMultMM<true,true>(A,B,C);
                else if (B.iscm() && C.iscm()) 
                    ColAddMultMM<true,true>(A,B,C);
                else if (B.colsize() > B.rowsize()) 
                    ColAddMultMM<false,true>(A,B,C);
                else 
                    RowAddMultMM<false,true>(A,B,C);
            }
            else {
                if (B.isrm() && C.isrm()) 
                    RowAddMultMM<true,false>(A,B,C);
                else if (B.iscm() && C.iscm()) 
                    ColAddMultMM<true,false>(A,B,C);
                else if (B.colsize() > B.rowsize()) 
                    ColAddMultMM<false,false>(A,B,C);
                else 
                    RowAddMultMM<false,false>(A,B,C);
            }
        }

#ifdef XDEBUG
        if (!(Norm(Matrix<T>(C)-C2) <= 
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+Norm(C0)))) {
            cerr<<"AddMultMM: alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<" step "<<
                A.diag().step()<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C0<<endl;
            cerr<<"-> C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    // C (+)= alpha * A * B
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
#ifdef XDEBUG
        Matrix<Ta> A0 = A;
        Matrix<Tb> B0 = B;
        Matrix<T> C0 = C;
        Matrix<T> C2 = alpha*A0*B0;
        if (add) C2 += C0;
#endif

        if (C.colsize() > 0 && C.rowsize() > 0) {
            if (alpha==T(0)) {
                if (!add) C.setZero();
            } else if (SameStorage(A.diag(),C)) {
                DiagMatrix<T> tempA = A;
                MultMM<add>(alpha,tempA,B,C);
            } else if (!add) {
                C = B;
                MultEqMM(alpha,A,C);
            } else if (SameStorage(B,C)) {
                if (B.isrm()) {
                    Matrix<T,RowMajor> tempB = B;
                    MultEqMM(alpha,A,tempB.view());
                    C += tempB;
                } else {
                    Matrix<T,ColMajor> tempB = B;
                    MultEqMM(alpha,A,tempB.view());
                    C += tempB;
                }
            } else {
                AddMultMM(alpha,A,B,C);
            }
        }
#ifdef XDEBUG
        if (1(Norm(Matrix<T>(C)-C2) <= 
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(B0)+
                     add?Norm(C0):TMV_RealType(T)(0)))) {
            cerr<<"MultMM: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<" step "<<A.diag().step()<<"  "<<A0<<endl;
            cerr<<"B = "<<TMV_Text(B)<<"  "<<B0<<endl;
            cerr<<"C = "<<TMV_Text(C)<<"  "<<C0<<endl;
            cerr<<"-> C = "<<C<<endl;
            cerr<<"C2 = "<<C2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_MultDM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


