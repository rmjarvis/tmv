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

//#define PRINTALGO_LU

#include "TMV_Blas.h"
#include "tmv/TMV_LUDecompose.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_SwapV.h"
#include "tmv/TMV_MinMax.h"
#include "tmv/TMV_CopyM.h"

namespace tmv {

#ifdef ALAP
    template <class T> 
    static inline void LapLU_Decompose(
        MatrixView<T,ColMajor>& A, int* P, int& signdet)
    { InlineLU_Decompose(A,P,signdet); }
#ifdef TMV_INST_DOUBLE
    static void LapLU_Decompose(
        MatrixView<double,ColMajor>& A, int* P, int& signdet)
    {
        TMVAssert(A.iscm());
        TMVAssert(A.ct()==NonConj);

        int m = A.colsize();
        int n = A.rowsize();
        int lda = A.stepj();
        int* lap_p = new int[n];
        LAPNAME(dgetrf) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p) LAPINFO);
        LAP_Results("dgetrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = lap_p[i] LAPMINUS1;
            if (P[i]!=i) signdet = -signdet;
        }
        delete [] lap_p;
    }
    static void LapLU_Decompose(
        MatrixView<std::complex<double>,ColMajor>& A, int* P, int& signdet)
    {
        TMVAssert(A.iscm());
        TMVAssert(A.ct()==NonConj);

        int m = A.colsize();
        int n = A.rowsize();
        int lda = A.stepj();
        int* lap_p = new int[n];
        LAPNAME(zgetrf) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p) LAPINFO);
        LAP_Results("zgetrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = lap_p[i] LAPMINUS1;
            if (P[i]!=i) signdet = -signdet;
        }
        delete [] lap_p;
    }
#endif
#ifdef TMV_INST_FLOAT
#ifndef MKL
    // This is giving me a weird runtime error sometimes with MKL:
    //   OMP abort: Unable to set worker thread stack size to 2098176 bytes
    //   Try reducing KMP_STACKSIZE or increasing the shell stack limit.
    // So I'm cutting it out for MKL compilations
    static void LapLU_Decompose(
        MatrixView<float,ColMajor>& A, int* P, int& signdet)
    {
        int m = A.colsize();
        int n = A.rowsize();
        int lda = A.stepj();
        int* lap_p = new int[n];
        LAPNAME(sgetrf) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p) LAPINFO);
        LAP_Results("sgetrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = lap_p[i] LAPMINUS1;
            if (P[i]!=i) signdet = -signdet;
        }
        delete [] lap_p;
    }
    static void LapLU_Decompose(
        MatrixView<std::complex<float>,ColMajor>& A, int* P, int& signdet)
    {
        int m = A.colsize();
        int n = A.rowsize();
        int lda = A.stepj();
        int* lap_p = new int[n];
        LAPNAME(cgetrf) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p) LAPINFO);
        LAP_Results("cgetrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = lap_p[i] LAPMINUS1;
            if (P[i]!=i) signdet = -signdet;
        }
        delete [] lap_p;
    }
#endif // MKL
#endif // FLOAT
#endif // ALAP

    template <class T> 
    void InstLU_Decompose(MatrixView<T> A, int* P, int& signdet)
    {
        //std::cout<<"Start LUDecompose:\n";
        //std::cout<<"A = "<<A<<std::endl;
        //Matrix<T> A0 = A;
        if (A.colsize() > 0 && A.rowsize() > 0) {
            if (A.iscm()) {
                MatrixView<T,ColMajor> Acm = A;
#ifdef ALAP
                LapLU_Decompose(Acm,P,signdet);
#else
                InlineLU_Decompose(Acm,P,signdet);
#endif
            } else {
                Matrix<T,ColMajor|NoDivider> Ac = A;
                MatrixView<T,ColMajor> Acm = Ac.view();
                InlineLU_Decompose(Acm,P,signdet);
                InstCopy(Ac.constView().xView(),A);
            }
        }
        //std::cout<<"A -> "<<A<<std::endl;
        //Matrix<T> lu = A.unitLowerTri() * A.upperTri();
        //std::cout<<"LU = "<<lu<<std::endl;
        //lu.reversePermuteRows(P);
        //std::cout<<"PLU = "<<lu<<std::endl;
        //std::cout<<"diff = ";
        //(lu-A0).write(std::cout,1.e-8);
        //std::cout<<"\nNorm(diff) = "<<Norm(lu-A0)<<std::endl;
    }

#define InstFile "TMV_LUDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv

