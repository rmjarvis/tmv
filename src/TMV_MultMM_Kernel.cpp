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

#define TMV_CompilingLibrary

#include "TMV_Blas.h"
#include "tmv/TMV_MultMM_Kernel.h"

namespace tmv {

#ifndef BLAS
#ifdef TMV_INST_FLOAT
#ifdef __SSE__
    template void multmm_16_16_32(
        const int M, const int N, const int K,
        const Scaling<1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_16_16_32(
        const int M, const int N, const int K,
        const Scaling<-1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_M_16_32(
        const int M, const int N, const int K,
        const Scaling<1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_M_16_32(
        const int M, const int N, const int K,
        const Scaling<-1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_16_N_32(
        const int M, const int N, const int K,
        const Scaling<1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_16_N_32(
        const int M, const int N, const int K,
        const Scaling<-1,float>& x,
        const float* A, const float* B, float* C0);



    template void multmm_16_16_64(
        const int M, const int N, const int K,
        const Scaling<1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_16_16_64(
        const int M, const int N, const int K,
        const Scaling<-1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_M_16_64(
        const int M, const int N, const int K,
        const Scaling<1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_M_16_64(
        const int M, const int N, const int K,
        const Scaling<-1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_16_N_64(
        const int M, const int N, const int K,
        const Scaling<1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_16_N_64(
        const int M, const int N, const int K,
        const Scaling<-1,float>& x,
        const float* A, const float* B, float* C0);
#endif
#endif

#ifdef TMV_INST_DOUBLE
#ifdef __SSE2__
    template void multmm_16_16_32(
        const int M, const int N, const int K,
        const Scaling<1,double>& x,
        const double* A, const double* B, double* C0);
    template void multmm_16_16_32(
        const int M, const int N, const int K,
        const Scaling<-1,double>& x,
        const double* A, const double* B, double* C0);
    template void multmm_M_16_32(
        const int M, const int N, const int K,
        const Scaling<1,double>& x,
        const double* A, const double* B, double* C0);
    template void multmm_M_16_32(
        const int M, const int N, const int K,
        const Scaling<-1,double>& x,
        const double* A, const double* B, double* C0);
    template void multmm_16_N_32(
        const int M, const int N, const int K,
        const Scaling<1,double>& x,
        const double* A, const double* B, double* C0);
    template void multmm_16_N_32(
        const int M, const int N, const int K,
        const Scaling<-1,double>& x,
        const double* A, const double* B, double* C0);

#endif
#endif
#endif

#define TMV_INST_SKIP_BLAS

#define InstFile "TMV_MultMM_Kernel.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


