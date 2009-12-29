///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



//#include <iostream>

#include "TMV_Blas.h" // Sets BLAS if appropriate
#include "TMV_MultMM.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MatrixArith.h"

namespace tmv {

    template <bool ca, class T, class Ta, class Tb> 
    static void ColMult(
        int M, const Tb b00, const Ta* A, T* C)
    {
        for(int i=M/8;i;--i) {
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
    extern void RecursiveCCCMultMM(
        const int M, const int N, const int K,
        const Ta* A, const Tb* B, T* C,
        const int Ask, const int Bsj, const int Csj)
    {
        if (K > N) {
            int K1 = K/2;
            RecursiveCCCMultMM<ca,cb>(M,N,K1,A,B,C,Ask,Bsj,Csj);
            RecursiveCCCMultMM<ca,cb>(M,N,K-K1,A+K1*Ask,B+K1,C,Ask,Bsj,Csj);
        } else if (N > 1) {
            int N1 = N/2;
            RecursiveCCCMultMM<ca,cb>(M,N1,K,A,B,C,Ask,Bsj,Csj);
            RecursiveCCCMultMM<ca,cb>(M,N-N1,K,A,B+N1*Bsj,C+N1*Csj,Ask,Bsj,Csj);
        } else {
            ColMult<ca>(M,(cb ? TMV_CONJ(B[0]) : B[0]),A,C);
        }
    }

    static void makeTaskList(
        int M,int N, int i, int j, 
        std::vector<int>& ilist, std::vector<int>& jlist, int& index)
    {
        if (M > N) {
            int M1 = M/2;
            makeTaskList(M1,N,i,j,ilist,jlist,index);
            makeTaskList(M-M1,N,i+M1,j,ilist,jlist,index);
        } else if (N > 1) {
            int N1 = N/2;
            makeTaskList(M,N1,i,j,ilist,jlist,index);
            makeTaskList(M,N-N1,i,j+N1,ilist,jlist,index);
        } else {
            ilist[index] = i;
            jlist[index] = j;
            index++;
        }
    }

    template <bool add, bool ca, bool cb, class T, class Ta, class Tb> 
    extern void DoCCCMultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
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
        TMVAssert(B.iscm());
        TMVAssert(C.iscm());

        const int M = C.colsize();
        const int N = C.rowsize();
        const int K = A.rowsize();
        const int Ask = A.stepj();
        const int Bsj = B.stepj();

        const Ta* Ap = A.cptr();
        const Tb* Bp = B.cptr();

        const int MB = MMCCC_BLOCKSIZE_M;
        const int NB = MMCCC_BLOCKSIZE_N*sizeof(double)/sizeof(T);
        const int Mb = M/MB;
        const int Nb = N/NB;
        if (Mb || Nb) {
            const int M2 = M-Mb*MB;
            const int N2 = N-Nb*NB;
            const int M1 = M-M2;
            const int N1 = N-N2;

            // This is best done recursively to avoid cache misses as much 
            // as possible, but openmp doesn't do recursion very well (yet?), 
            // so I just figure out the recursive order first and do a 
            // simple loop that openmp can divvy up appropriately.
            const int MbNb = Mb*Nb;
            std::vector<int> ilist(MbNb);
            std::vector<int> jlist(MbNb);
            if (MbNb) {
                int listindex=0;
                makeTaskList(Mb,Nb,0,0,ilist,jlist,listindex);
            }
#ifdef _OPENMP
#pragma omp parallel
            {
                // OpenMP seems to have trouble with new [], so just use a 
                // static array.  
                // However, the compiler needs to know the size at compile time
                // so a little wasted memory here, always using the full MB*NB.
                T Ct[ MB * NB ];
                MatrixView<T> Ctemp = 
                    MatrixViewOf(Ct,TMV_MIN(M,MB),TMV_MIN(N,NB),ColMajor);
#else
                Matrix<T,ColMajor> Ctemp(TMV_MIN(M,MB),TMV_MIN(N,NB));
                T* Ct = Ctemp.ptr();
#endif
                const int Ctsj = Ctemp.stepj();

                // Full MBxNB blocks in C matrix
                if (Mb && Nb) {
#ifdef _OPENMP
#pragma omp for nowait schedule(guided)
#endif
                    for(int ij=0;ij<int(MbNb);ij++) {
                        const int i=ilist[ij];
                        const int j=jlist[ij];
                        const int ii = i*MB;
                        const int jj = j*NB;

                        Ctemp.zero();
                        RecursiveCCCMultMM<ca,cb>(
                            MB,NB,K,Ap+ii,Bp+jj*Bsj,Ct,Ask,Bsj,Ctsj);
                        if (alpha != T(1)) Ctemp *= alpha;
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            if (add) C.subMatrix(ii,ii+MB,jj,jj+NB) += Ctemp;
                            else C.subMatrix(ii,ii+MB,jj,jj+NB) = Ctemp;
                        }
                    }
                } // Mb && Nb

                // MB x N2 partial blocks:
                if (Mb && N2) {
#ifdef _OPENMP
#pragma omp for nowait schedule(guided)
#endif
                    for(int i=0;i<int(Mb);i++) {
                        const int ii = i*MB;
                        Ctemp.zero();
                        RecursiveCCCMultMM<ca,cb>(
                            MB,N2,K,Ap+ii,Bp+N1*Bsj,Ct,Ask,Bsj,Ctsj);
                        if (alpha != T(1)) Ctemp.colRange(0,N2) *= alpha;
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            if (add) C.subMatrix(ii,ii+MB,N1,N) += 
                                Ctemp.colRange(0,N2);
                            else C.subMatrix(ii,ii+MB,N1,N) = 
                                Ctemp.colRange(0,N2);
                        }
                    }
                } // Mb && N2

                // M2 x NB partial blocks:
                if (M2 && Nb) {
#ifdef _OPENMP
#pragma omp for nowait schedule(guided)
#endif
                    for(int j=0;j<int(Nb);j++) {
                        const int jj = j*NB;
                        Ctemp.zero();
                        RecursiveCCCMultMM<ca,cb>(
                            M2,NB,K,Ap+M1,Bp+jj*Bsj,Ct,Ask,Bsj,Ctsj);
                        if (alpha != T(1)) Ctemp.rowRange(0,M2) *= alpha;
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            if (add) C.subMatrix(M1,M,jj,jj+NB) += 
                                Ctemp.rowRange(0,M2);
                            else C.subMatrix(M1,M,jj,jj+NB) = 
                                Ctemp.rowRange(0,M2);
                        }
                    }
                } // M2 && Nb

                // Final M2 x N2 partial block
                if (M2 && N2) {
#ifdef _OPENMP
#pragma omp single
#endif
                    {
                        Ctemp.zero();
                        RecursiveCCCMultMM<ca,cb>(
                            M2,N2,K,Ap+M1,Bp+N1*Bsj,Ct,Ask,Bsj,Ctsj);
                        if (alpha != T(1)) Ctemp.subMatrix(0,M2,0,N2) *= alpha;
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            if (add) C.subMatrix(M1,M,N1,N) += 
                                Ctemp.subMatrix(0,M2,0,N2);
                            else C.subMatrix(M1,M,N1,N) = 
                                Ctemp.subMatrix(0,M2,0,N2);
                        }
                    } 
                } // M2 && N2
#ifdef _OPENMP
            } // end parallel
#endif
        } else {
            Matrix<T,ColMajor> Ctemp(M,N,T(0));
            T* Ct = Ctemp.ptr();
            const int Ctsj = Ctemp.stepj();
            RecursiveCCCMultMM<ca,cb>(M,N,K,Ap,Bp,Ct, Ask,Bsj,Ctsj);
            if (alpha != T(1)) Ctemp *= alpha;
            if (add) C += Ctemp;
            else C = Ctemp;
        }
    }

    template <bool add, class T, class Ta, class Tb> void CCCMultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        if (A.isconj())
            if (B.isconj())
                DoCCCMultMM<add,true,true>(alpha,A,B,C);
            else
                DoCCCMultMM<add,true,false>(alpha,A,B,C);
        else
            if (B.isconj())
                DoCCCMultMM<add,false,true>(alpha,A,B,C);
            else
                DoCCCMultMM<add,false,false>(alpha,A,B,C);
    }

#ifdef BLAS
#define INST_SKIP_BLAS
#endif

#define InstFile "TMV_MultMM_CCC.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


