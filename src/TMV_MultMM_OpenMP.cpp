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


#ifdef _OPENMP

#include "TMV_Blas.h"
#include "TMV_MultMM.h"
#include <omp.h>

namespace tmv {

    template <bool add, class T, class Ta, class Tb> 
    void OpenMPMultMM(
        const T x, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C)
    {
        const ptrdiff_t M = C.colsize();
        const ptrdiff_t N = C.rowsize();

#pragma omp parallel
        {
            ptrdiff_t num_threads = omp_get_num_threads();
            ptrdiff_t mythread = omp_get_thread_num();
            if (num_threads == 1) {
                BlockMultMM<add>(x,A,B,C);
            } else if (M > N) {
                ptrdiff_t Mx = M / num_threads;
                Mx = ((((Mx-1)>>4)+1)<<4); // round up to mult of 16
                ptrdiff_t i1 = mythread * Mx;
                ptrdiff_t i2 = (mythread+1) * Mx;
                if (i2 > M || mythread == num_threads-1) i2 = M;
                if (i1 < M) {
                    // Need to make sure, since we rounded up Mx!
                    BlockMultMM<add>(x,A.rowRange(i1,i2),B,C.rowRange(i1,i2));
                }
            } else {
                ptrdiff_t Nx = N / num_threads;
                Nx = ((((Nx-1)>>4)+1)<<4); 
                ptrdiff_t j1 = mythread * Nx;
                ptrdiff_t j2 = (mythread+1) * Nx;
                if (j2 > N || mythread == num_threads-1) j2 = N;
                if (j1 < N) {
                    BlockMultMM<add>(x,A,B.colRange(j1,j2),C.colRange(j1,j2));
                }
            }
        }
    }

#ifdef BLAS
#define INST_SKIP_BLAS
#endif

#define InstFile "TMV_MultMM_OpenMP.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv

#endif

