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


//#include <iostream>

#include "TMV_Blas.h" // Sets BLAS if appropriate
#include "TMV_MultMM.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MatrixArith.h"

namespace tmv {

    template <bool ca, class T, class Ta, class Tb> 
    static void ColMult(ptrdiff_t M, const Tb b00, const Ta* A, T* C)
    {
        for(ptrdiff_t i=M/8;i;--i) {
            C[0] += (ca ? TMV_CONJ(A[0]) : A[0]) * b00;
            C[1] += (ca ? TMV_CONJ(A[1]) : A[1]) * b00;
            C[2] += (ca ? TMV_CONJ(A[2]) : A[2]) * b00;
            C[3] += (ca ? TMV_CONJ(A[3]) : A[3]) * b00;
            C[4] += (ca ? TMV_CONJ(A[4]) : A[4]) * b00;
            C[5] += (ca ? TMV_CONJ(A[5]) : A[5]) * b00;
            C[6] += (ca ? TMV_CONJ(A[6]) : A[6]) * b00;
            C[7] += (ca ? TMV_CONJ(A[7]) : A[7]) * b00;
            A += 8; C += 8;
        }
        M %= 8;
        if (M) {
            if (M >= 4) {
                C[0] += (ca ? TMV_CONJ(A[0]) : A[0]) * b00;
                C[1] += (ca ? TMV_CONJ(A[1]) : A[1]) * b00;
                C[2] += (ca ? TMV_CONJ(A[2]) : A[2]) * b00;
                C[3] += (ca ? TMV_CONJ(A[3]) : A[3]) * b00;
                M -= 4; A += 4; C += 4;
            }
            if (M >= 2) {
                C[0] += (ca ? TMV_CONJ(A[0]) : A[0]) * b00;
                C[1] += (ca ? TMV_CONJ(A[1]) : A[1]) * b00;
                M -= 2; A += 2; C += 2;
            }
            if (M) {
                C[0] += (ca ? TMV_CONJ(A[0]) : A[0]) * b00;
            }
        }
    }

    template <bool ca, bool cb, class T, class Ta, class Tb> 
    static void RecursiveCRCMultMM(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Ta* A, const Tb* B, T* C,
        const ptrdiff_t Ask, const ptrdiff_t Bsk, const ptrdiff_t Csj)
    {
        if (K > N) {
            ptrdiff_t K1 = K/2;
            RecursiveCRCMultMM<ca,cb>(M,N,K1,A,B,C,Ask,Bsk,Csj);
            RecursiveCRCMultMM<ca,cb>(M,N,K-K1,A+K1*Ask,B+K1*Bsk,C,Ask,Bsk,Csj);
        } else if (N > 1) {
            ptrdiff_t N1 = N/2;
            RecursiveCRCMultMM<ca,cb>(M,N1,K,A,B,C,Ask,Bsk,Csj);
            RecursiveCRCMultMM<ca,cb>(M,N-N1,K,A,B+N1,C+N1*Csj,Ask,Bsk,Csj);
        } else {
            ColMult<ca>(M,(cb ? TMV_CONJ(B[0]) : B[0]),A,C);
        }
    }

    template <bool add, bool ca, bool cb, class T, class Ta, class Tb> 
    static void DoCRCMultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        TMVAssert(A.colsize() == C.colsize());
        TMVAssert(A.rowsize() == B.colsize());
        TMVAssert(B.rowsize() == C.rowsize());
        TMVAssert(C.colsize() > 0);
        TMVAssert(C.rowsize() > 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(alpha != T(0));
        TMVAssert(ca == A.isconj());
        TMVAssert(cb == B.isconj());
        TMVAssert(C.ct()==NonConj);
        TMVAssert(A.iscm());
        TMVAssert(B.isrm());
        TMVAssert(C.iscm());

        const ptrdiff_t M = C.colsize();
        const ptrdiff_t N = C.rowsize();
        const ptrdiff_t K = A.rowsize();
        const ptrdiff_t Ask = A.stepj();
        const ptrdiff_t Bsk = B.stepi();

        const Ta* Ap = A.cptr();
        const Tb* Bp = B.cptr();

        Matrix<T,ColMajor> Ctemp(M,N,T(0));
        T* Ct = Ctemp.ptr();
        const ptrdiff_t Ctsj = Ctemp.stepj();
        RecursiveCRCMultMM<ca,cb>(M,N,K,Ap,Bp,Ct, Ask,Bsk,Ctsj);
        if (alpha != T(1)) Ctemp *= alpha;
        if (add) C += Ctemp;
        else C = Ctemp;
    }

    template <bool add, class T, class Ta, class Tb> void CRCMultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        if (A.isconj())
            if (B.isconj())
                DoCRCMultMM<add,true,true>(alpha,A,B,C);
            else
                DoCRCMultMM<add,true,false>(alpha,A,B,C);
        else
            if (B.isconj())
                DoCRCMultMM<add,false,true>(alpha,A,B,C);
            else
                DoCRCMultMM<add,false,false>(alpha,A,B,C);
    }

#ifdef BLAS
#define INST_SKIP_BLAS
#endif

#define InstFile "TMV_MultMM_CRC.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


