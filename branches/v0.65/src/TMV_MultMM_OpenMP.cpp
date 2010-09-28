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


#ifdef _OPENMP

#include "TMV_Blas.h"
#include "TMV_MultMM.h"
#include <omp.h>

namespace tmv {

    template <bool add, class T, class Ta, class Tb> 
    void OpenMPMultMM(
        const T x, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C)
    {
        const int M = C.colsize();
        const int N = C.rowsize();

#pragma omp parallel
        {
            int num_threads = omp_get_num_threads();
            int mythread = omp_get_thread_num();
            if (num_threads == 1) {
                BlockMultMM<add>(x,A,B,C);
            } else if (M > N) {
                int Mx = M / num_threads;
                Mx = ((((Mx-1)>>4)+1)<<4); // round up to mult of 16
                int i1 = mythread * Mx;
                int i2 = (mythread+1) * Mx;
                if (i2 > M || mythread == num_threads-1) i2 = M;
                if (i1 < M) {
                    // Need to make sure, since we rounded up Mx!
                    BlockMultMM<add>(x,A.rowRange(i1,i2),B,C.rowRange(i1,i2));
                }
            } else {
                int Nx = N / num_threads;
                Nx = ((((Nx-1)>>4)+1)<<4); 
                int j1 = mythread * Nx;
                int j2 = (mythread+1) * Nx;
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

