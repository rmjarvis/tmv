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

#ifndef TMV_MultMM_Kernel_H
#define TMV_MultMM_Kernel_H

// This file has the kernel routines for doing a single block matrix product.
// These are called from the Block algorithm in TMV_MultMM_Block.cpp
//
// There are different kernels depending on what is known at compile time.
// e.g. multmm_16_N_64 has M = 16, K = 64, and N is a runtime value.

#if defined(__SSE2__) || defined(__SSE__)
#include "xmmintrin.h"
#endif

namespace tmv {

    // First the generic implementations for any type T.
    // If we have SSE commands, therea are specialized double and float
    // versions below.
    
    template <typename T>
    static void generic_multmm_16_16_16(
        const T* A0, const T* B0, T* C0)
    {
        T a0, a1;
        T b0, b1;
        T C00, C01, C10, C11;
        T c00, c01, c10, c11;

        ptrdiff_t i2,j2;

        j2 = 8; do {
            b0 = B0[0]; 
            b1 = B0[16]; 
            i2 = 4; do {
                // 1
                a0 = A0[0];  c00 = a0 * b0; C00 = c00;
                a1 = A0[16]; c10 = a1 * b0; C10 = c10;
                c01 = a0 * b1; C01 = c01;
                b0 = B0[1];  c11 = a1 * b1; C11 = c11;

                // 2
                a0 = A0[1];  c00 = a0 * b0; C00 += c00;
                a1 = A0[17]; c10 = a1 * b0; C10 += c10;
                b1 = B0[17]; c01 = a0 * b1; C01 += c01;
                b0 = B0[2];  c11 = a1 * b1; C11 += c11;

                // 3
                a0 = A0[2];  c00 = a0 * b0; C00 += c00;
                a1 = A0[18]; c10 = a1 * b0; C10 += c10;
                b1 = B0[18]; c01 = a0 * b1; C01 += c01;
                b0 = B0[3];  c11 = a1 * b1; C11 += c11;

                // 4
                a0 = A0[3];  c00 = a0 * b0; C00 += c00;
                a1 = A0[19]; c10 = a1 * b0; C10 += c10;
                b1 = B0[19]; c01 = a0 * b1; C01 += c01;
                b0 = B0[4];  c11 = a1 * b1; C11 += c11;

                // 5
                a0 = A0[4];  c00 = a0 * b0; C00 += c00;
                a1 = A0[20]; c10 = a1 * b0; C10 += c10;
                b1 = B0[20]; c01 = a0 * b1; C01 += c01;
                b0 = B0[5];  c11 = a1 * b1; C11 += c11;

                // 6
                a0 = A0[5];  c00 = a0 * b0; C00 += c00;
                a1 = A0[21]; c10 = a1 * b0; C10 += c10;
                b1 = B0[21]; c01 = a0 * b1; C01 += c01;
                b0 = B0[6];  c11 = a1 * b1; C11 += c11;

                // 7
                a0 = A0[6];  c00 = a0 * b0; C00 += c00;
                a1 = A0[22]; c10 = a1 * b0; C10 += c10;
                b1 = B0[22]; c01 = a0 * b1; C01 += c01;
                b0 = B0[7];  c11 = a1 * b1; C11 += c11;

                // 8
                a0 = A0[7];  c00 = a0 * b0; C00 += c00;
                a1 = A0[23]; c10 = a1 * b0; C10 += c10;
                b1 = B0[23]; c01 = a0 * b1; C01 += c01;
                b0 = B0[8];  c11 = a1 * b1; C11 += c11;

                // 9
                a0 = A0[8];  c00 = a0 * b0; C00 += c00;
                a1 = A0[24]; c10 = a1 * b0; C10 += c10;
                b1 = B0[24]; c01 = a0 * b1; C01 += c01;
                b0 = B0[9];  c11 = a1 * b1; C11 += c11;

                // 10
                a0 = A0[9];  c00 = a0 * b0; C00 += c00;
                a1 = A0[25]; c10 = a1 * b0; C10 += c10;
                b1 = B0[25]; c01 = a0 * b1; C01 += c01;
                b0 = B0[10]; c11 = a1 * b1; C11 += c11;

                // 11
                a0 = A0[10]; c00 = a0 * b0; C00 += c00;
                a1 = A0[26]; c10 = a1 * b0; C10 += c10;
                b1 = B0[26]; c01 = a0 * b1; C01 += c01;
                b0 = B0[11]; c11 = a1 * b1; C11 += c11;

                // 12
                a0 = A0[11]; c00 = a0 * b0; C00 += c00;
                a1 = A0[27]; c10 = a1 * b0; C10 += c10;
                b1 = B0[27]; c01 = a0 * b1; C01 += c01;
                b0 = B0[12]; c11 = a1 * b1; C11 += c11;

                // 13
                a0 = A0[12]; c00 = a0 * b0; C00 += c00;
                a1 = A0[28]; c10 = a1 * b0; C10 += c10;
                b1 = B0[28]; c01 = a0 * b1; C01 += c01;
                b0 = B0[13]; c11 = a1 * b1; C11 += c11;

                // 14
                a0 = A0[13]; c00 = a0 * b0; C00 += c00;
                a1 = A0[29]; c10 = a1 * b0; C10 += c10;
                b1 = B0[29]; c01 = a0 * b1; C01 += c01;
                b0 = B0[14]; c11 = a1 * b1; C11 += c11;

                // 15
                a0 = A0[14]; c00 = a0 * b0; C00 += c00;
                a1 = A0[30]; c10 = a1 * b0; C10 += c10;
                b1 = B0[30]; c01 = a0 * b1; C01 += c01;
                b0 = B0[15]; c11 = a1 * b1; C11 += c11;

                // 16
                a0 = A0[15]; c00 = a0 * b0; C00 += c00;
                C0[0] += C00;
                a1 = A0[31]; c10 = a1 * b0; C10 += c10;
                C0[1] += C10;
                b1 = B0[31]; c01 = a0 * b1; C01 += c01;
                C0[16] += C01;
                c11 = a1 * b1; C11 += c11;
                C0[17] += C11;

                // 16
                a0 = A0[47]; c00 = a0 * b0; C00 = c00;
                a1 = A0[63]; c10 = a1 * b0; C10 = c10;
                c01 = a0 * b1; C01 = c01;
                b0 = B0[14]; c11 = a1 * b1; C11 = c11;

                // 15
                a0 = A0[46]; c00 = a0 * b0; C00 += c00;
                a1 = A0[62]; c10 = a1 * b0; C10 += c10;
                b1 = B0[30]; c01 = a0 * b1; C01 += c01;
                b0 = B0[13]; c11 = a1 * b1; C11 += c11;

                // 14
                a0 = A0[45]; c00 = a0 * b0; C00 += c00;
                a1 = A0[61]; c10 = a1 * b0; C10 += c10;
                b1 = B0[29]; c01 = a0 * b1; C01 += c01;
                b0 = B0[12]; c11 = a1 * b1; C11 += c11;

                // 13
                a0 = A0[44]; c00 = a0 * b0; C00 += c00;
                a1 = A0[60]; c10 = a1 * b0; C10 += c10;
                b1 = B0[28]; c01 = a0 * b1; C01 += c01;
                b0 = B0[11]; c11 = a1 * b1; C11 += c11;

                // 12
                a0 = A0[43]; c00 = a0 * b0; C00 += c00;
                a1 = A0[59]; c10 = a1 * b0; C10 += c10;
                b1 = B0[27]; c01 = a0 * b1; C01 += c01;
                b0 = B0[10]; c11 = a1 * b1; C11 += c11;

                // 11
                a0 = A0[42]; c00 = a0 * b0; C00 += c00;
                a1 = A0[58]; c10 = a1 * b0; C10 += c10;
                b1 = B0[26]; c01 = a0 * b1; C01 += c01;
                b0 = B0[9];  c11 = a1 * b1; C11 += c11;

                // 10
                a0 = A0[41]; c00 = a0 * b0; C00 += c00;
                a1 = A0[57]; c10 = a1 * b0; C10 += c10;
                b1 = B0[25]; c01 = a0 * b1; C01 += c01;
                b0 = B0[8];  c11 = a1 * b1; C11 += c11;

                // 9
                a0 = A0[40]; c00 = a0 * b0; C00 += c00;
                a1 = A0[56]; c10 = a1 * b0; C10 += c10;
                b1 = B0[24]; c01 = a0 * b1; C01 += c01;
                b0 = B0[7];  c11 = a1 * b1; C11 += c11;

                // 8
                a0 = A0[39]; c00 = a0 * b0; C00 += c00;
                a1 = A0[55]; c10 = a1 * b0; C10 += c10;
                b1 = B0[23]; c01 = a0 * b1; C01 += c01;
                b0 = B0[6];  c11 = a1 * b1; C11 += c11;

                // 7
                a0 = A0[38]; c00 = a0 * b0; C00 += c00;
                a1 = A0[54]; c10 = a1 * b0; C10 += c10;
                b1 = B0[22]; c01 = a0 * b1; C01 += c01;
                b0 = B0[5];  c11 = a1 * b1; C11 += c11;

                // 6
                a0 = A0[37]; c00 = a0 * b0; C00 += c00;
                a1 = A0[53]; c10 = a1 * b0; C10 += c10;
                b1 = B0[21]; c01 = a0 * b1; C01 += c01;
                b0 = B0[4];  c11 = a1 * b1; C11 += c11;

                // 5
                a0 = A0[36]; c00 = a0 * b0; C00 += c00;
                a1 = A0[52]; c10 = a1 * b0; C10 += c10;
                b1 = B0[20]; c01 = a0 * b1; C01 += c01;
                b0 = B0[3];  c11 = a1 * b1; C11 += c11;

                // 4
                a0 = A0[35]; c00 = a0 * b0; C00 += c00;
                a1 = A0[51]; c10 = a1 * b0; C10 += c10;
                b1 = B0[19]; c01 = a0 * b1; C01 += c01;
                b0 = B0[2];  c11 = a1 * b1; C11 += c11;

                // 3
                a0 = A0[34]; c00 = a0 * b0; C00 += c00;
                a1 = A0[50]; c10 = a1 * b0; C10 += c10;
                b1 = B0[18]; c01 = a0 * b1; C01 += c01;
                b0 = B0[1];  c11 = a1 * b1; C11 += c11;

                // 2
                a0 = A0[33]; c00 = a0 * b0; C00 += c00;
                a1 = A0[49]; c10 = a1 * b0; C10 += c10;
                b1 = B0[17]; c01 = a0 * b1; C01 += c01;
                b0 = B0[0];  c11 = a1 * b1; C11 += c11;

                // 1
                a0 = A0[32]; c00 = a0 * b0; C00 += c00;
                C0[ 2] += C00;
                a1 = A0[48]; c10 = a1 * b0; C10 += c10;
                C0[ 3] += C10;
                b1 = B0[16]; c01 = a0 * b1; C01 += c01;
                C0[18] += C01;
                A0 += 64;    c11 = a1 * b1; C11 += c11;
                C0[19] += C11;
                C0 += 4;
            } while (--i2);
            A0 -= 256;
            B0 += 32;
            C0 += 16; 
        } while (--j2);
    }
    template <typename T>
    static void generic_multmm_M_16_16(
        const ptrdiff_t M, const T* A0, const T* B0, T* C0)
    {
        const ptrdiff_t M_4 = M>>2; // M_4 = M/4
        const ptrdiff_t Mc = M-(M_4<<2); // Mc = M%4
        const ptrdiff_t Mx16 = M<<4; // = M*16

        T a0, a1;
        T b0, b1;
        T C00, C01, C10, C11;
        T c00, c01, c10, c11;

        T* C1 = C0 + M;

        ptrdiff_t j,i;

        j = 8; do {
            b0 = B0[0]; 
            b1 = B0[16]; 
            i = M_4; if (i) do {
                // 1
                a0 = A0[0];  c00 = a0 * b0; C00 = c00;
                a1 = A0[16]; c10 = a1 * b0; C10 = c10;
                c01 = a0 * b1; C01 = c01;
                b0 = B0[1];  c11 = a1 * b1; C11 = c11;

                // 2
                a0 = A0[1];  c00 = a0 * b0; C00 += c00;
                a1 = A0[17]; c10 = a1 * b0; C10 += c10;
                b1 = B0[17]; c01 = a0 * b1; C01 += c01;
                b0 = B0[2];  c11 = a1 * b1; C11 += c11;

                // 3
                a0 = A0[2];  c00 = a0 * b0; C00 += c00;
                a1 = A0[18]; c10 = a1 * b0; C10 += c10;
                b1 = B0[18]; c01 = a0 * b1; C01 += c01;
                b0 = B0[3];  c11 = a1 * b1; C11 += c11;

                // 4
                a0 = A0[3];  c00 = a0 * b0; C00 += c00;
                a1 = A0[19]; c10 = a1 * b0; C10 += c10;
                b1 = B0[19]; c01 = a0 * b1; C01 += c01;
                b0 = B0[4];  c11 = a1 * b1; C11 += c11;

                // 5
                a0 = A0[4];  c00 = a0 * b0; C00 += c00;
                a1 = A0[20]; c10 = a1 * b0; C10 += c10;
                b1 = B0[20]; c01 = a0 * b1; C01 += c01;
                b0 = B0[5];  c11 = a1 * b1; C11 += c11;

                // 6
                a0 = A0[5];  c00 = a0 * b0; C00 += c00;
                a1 = A0[21]; c10 = a1 * b0; C10 += c10;
                b1 = B0[21]; c01 = a0 * b1; C01 += c01;
                b0 = B0[6];  c11 = a1 * b1; C11 += c11;

                // 7
                a0 = A0[6];  c00 = a0 * b0; C00 += c00;
                a1 = A0[22]; c10 = a1 * b0; C10 += c10;
                b1 = B0[22]; c01 = a0 * b1; C01 += c01;
                b0 = B0[7];  c11 = a1 * b1; C11 += c11;

                // 8
                a0 = A0[7];  c00 = a0 * b0; C00 += c00;
                a1 = A0[23]; c10 = a1 * b0; C10 += c10;
                b1 = B0[23]; c01 = a0 * b1; C01 += c01;
                b0 = B0[8];  c11 = a1 * b1; C11 += c11;

                // 9
                a0 = A0[8];  c00 = a0 * b0; C00 += c00;
                a1 = A0[24]; c10 = a1 * b0; C10 += c10;
                b1 = B0[24]; c01 = a0 * b1; C01 += c01;
                b0 = B0[9];  c11 = a1 * b1; C11 += c11;

                // 10
                a0 = A0[9];  c00 = a0 * b0; C00 += c00;
                a1 = A0[25]; c10 = a1 * b0; C10 += c10;
                b1 = B0[25]; c01 = a0 * b1; C01 += c01;
                b0 = B0[10]; c11 = a1 * b1; C11 += c11;

                // 11
                a0 = A0[10]; c00 = a0 * b0; C00 += c00;
                a1 = A0[26]; c10 = a1 * b0; C10 += c10;
                b1 = B0[26]; c01 = a0 * b1; C01 += c01;
                b0 = B0[11]; c11 = a1 * b1; C11 += c11;

                // 12
                a0 = A0[11]; c00 = a0 * b0; C00 += c00;
                a1 = A0[27]; c10 = a1 * b0; C10 += c10;
                b1 = B0[27]; c01 = a0 * b1; C01 += c01;
                b0 = B0[12]; c11 = a1 * b1; C11 += c11;

                // 13
                a0 = A0[12]; c00 = a0 * b0; C00 += c00;
                a1 = A0[28]; c10 = a1 * b0; C10 += c10;
                b1 = B0[28]; c01 = a0 * b1; C01 += c01;
                b0 = B0[13]; c11 = a1 * b1; C11 += c11;

                // 14
                a0 = A0[13]; c00 = a0 * b0; C00 += c00;
                a1 = A0[29]; c10 = a1 * b0; C10 += c10;
                b1 = B0[29]; c01 = a0 * b1; C01 += c01;
                b0 = B0[14]; c11 = a1 * b1; C11 += c11;

                // 15
                a0 = A0[14]; c00 = a0 * b0; C00 += c00;
                a1 = A0[30]; c10 = a1 * b0; C10 += c10;
                b1 = B0[30]; c01 = a0 * b1; C01 += c01;
                b0 = B0[15]; c11 = a1 * b1; C11 += c11;

                // 16
                a0 = A0[15]; c00 = a0 * b0; C00 += c00;
                C0[0] += C00;
                a1 = A0[31]; c10 = a1 * b0; C10 += c10;
                C0[1] += C10;
                b1 = B0[31]; c01 = a0 * b1; C01 += c01;
                C1[0] += C01;
                c11 = a1 * b1; C11 += c11;
                C1[1] += C11;

                // 16
                a0 = A0[47]; c00 = a0 * b0; C00 = c00;
                a1 = A0[63]; c10 = a1 * b0; C10 = c10;
                c01 = a0 * b1; C01 = c01;
                b0 = B0[14]; c11 = a1 * b1; C11 = c11;

                // 15
                a0 = A0[46]; c00 = a0 * b0; C00 += c00;
                a1 = A0[62]; c10 = a1 * b0; C10 += c10;
                b1 = B0[30]; c01 = a0 * b1; C01 += c01;
                b0 = B0[13]; c11 = a1 * b1; C11 += c11;

                // 14
                a0 = A0[45]; c00 = a0 * b0; C00 += c00;
                a1 = A0[61]; c10 = a1 * b0; C10 += c10;
                b1 = B0[29]; c01 = a0 * b1; C01 += c01;
                b0 = B0[12]; c11 = a1 * b1; C11 += c11;

                // 13
                a0 = A0[44]; c00 = a0 * b0; C00 += c00;
                a1 = A0[60]; c10 = a1 * b0; C10 += c10;
                b1 = B0[28]; c01 = a0 * b1; C01 += c01;
                b0 = B0[11]; c11 = a1 * b1; C11 += c11;

                // 12
                a0 = A0[43]; c00 = a0 * b0; C00 += c00;
                a1 = A0[59]; c10 = a1 * b0; C10 += c10;
                b1 = B0[27]; c01 = a0 * b1; C01 += c01;
                b0 = B0[10]; c11 = a1 * b1; C11 += c11;

                // 11
                a0 = A0[42]; c00 = a0 * b0; C00 += c00;
                a1 = A0[58]; c10 = a1 * b0; C10 += c10;
                b1 = B0[26]; c01 = a0 * b1; C01 += c01;
                b0 = B0[9];  c11 = a1 * b1; C11 += c11;

                // 10
                a0 = A0[41]; c00 = a0 * b0; C00 += c00;
                a1 = A0[57]; c10 = a1 * b0; C10 += c10;
                b1 = B0[25]; c01 = a0 * b1; C01 += c01;
                b0 = B0[8];  c11 = a1 * b1; C11 += c11;

                // 9
                a0 = A0[40]; c00 = a0 * b0; C00 += c00;
                a1 = A0[56]; c10 = a1 * b0; C10 += c10;
                b1 = B0[24]; c01 = a0 * b1; C01 += c01;
                b0 = B0[7];  c11 = a1 * b1; C11 += c11;

                // 8
                a0 = A0[39]; c00 = a0 * b0; C00 += c00;
                a1 = A0[55]; c10 = a1 * b0; C10 += c10;
                b1 = B0[23]; c01 = a0 * b1; C01 += c01;
                b0 = B0[6];  c11 = a1 * b1; C11 += c11;

                // 7
                a0 = A0[38]; c00 = a0 * b0; C00 += c00;
                a1 = A0[54]; c10 = a1 * b0; C10 += c10;
                b1 = B0[22]; c01 = a0 * b1; C01 += c01;
                b0 = B0[5];  c11 = a1 * b1; C11 += c11;

                // 6
                a0 = A0[37]; c00 = a0 * b0; C00 += c00;
                a1 = A0[53]; c10 = a1 * b0; C10 += c10;
                b1 = B0[21]; c01 = a0 * b1; C01 += c01;
                b0 = B0[4];  c11 = a1 * b1; C11 += c11;

                // 5
                a0 = A0[36]; c00 = a0 * b0; C00 += c00;
                a1 = A0[52]; c10 = a1 * b0; C10 += c10;
                b1 = B0[20]; c01 = a0 * b1; C01 += c01;
                b0 = B0[3];  c11 = a1 * b1; C11 += c11;

                // 4
                a0 = A0[35]; c00 = a0 * b0; C00 += c00;
                a1 = A0[51]; c10 = a1 * b0; C10 += c10;
                b1 = B0[19]; c01 = a0 * b1; C01 += c01;
                b0 = B0[2];  c11 = a1 * b1; C11 += c11;

                // 3
                a0 = A0[34]; c00 = a0 * b0; C00 += c00;
                a1 = A0[50]; c10 = a1 * b0; C10 += c10;
                b1 = B0[18]; c01 = a0 * b1; C01 += c01;
                b0 = B0[1];  c11 = a1 * b1; C11 += c11;

                // 2
                a0 = A0[33]; c00 = a0 * b0; C00 += c00;
                a1 = A0[49]; c10 = a1 * b0; C10 += c10;
                b1 = B0[17]; c01 = a0 * b1; C01 += c01;
                b0 = B0[0];  c11 = a1 * b1; C11 += c11;

                // 1
                a0 = A0[32]; c00 = a0 * b0; C00 += c00;
                C0[ 2] += C00;
                a1 = A0[48]; c10 = a1 * b0; C10 += c10;
                C0[ 3] += C10;
                b1 = B0[16]; c01 = a0 * b1; C01 += c01;
                C1[ 2] += C01;
                A0 += 64;    c11 = a1 * b1; C11 += c11;
                C1[ 3] += C11;
                C0 += 4;
                C1 += 4;
            } while (--i);
            i = Mc; if (i) do {
                // 1
                a0 = A0[0];  c00 = a0 * b0; C00 = c00;
                b1 = B0[16]; c01 = a0 * b1; C01 = c01;
                b0 = B0[1];

                // 2
                a0 = A0[1];  c00 = a0 * b0; C00 += c00;
                b1 = B0[17]; c01 = a0 * b1; C01 += c01;
                b0 = B0[2];

                // 3
                a0 = A0[2];  c00 = a0 * b0; C00 += c00;
                b1 = B0[18]; c01 = a0 * b1; C01 += c01;
                b0 = B0[3];

                // 4
                a0 = A0[3];  c00 = a0 * b0; C00 += c00;
                b1 = B0[19]; c01 = a0 * b1; C01 += c01;
                b0 = B0[4];

                // 5
                a0 = A0[4];  c00 = a0 * b0; C00 += c00;
                b1 = B0[20]; c01 = a0 * b1; C01 += c01;
                b0 = B0[5];

                // 6
                a0 = A0[5];  c00 = a0 * b0; C00 += c00;
                b1 = B0[21]; c01 = a0 * b1; C01 += c01;
                b0 = B0[6];

                // 7
                a0 = A0[6];  c00 = a0 * b0; C00 += c00;
                b1 = B0[22]; c01 = a0 * b1; C01 += c01;
                b0 = B0[7];

                // 8
                a0 = A0[7];  c00 = a0 * b0; C00 += c00;
                b1 = B0[23]; c01 = a0 * b1; C01 += c01;
                b0 = B0[8];

                // 9
                a0 = A0[8];  c00 = a0 * b0; C00 += c00;
                b1 = B0[24]; c01 = a0 * b1; C01 += c01;
                b0 = B0[9];

                // 10
                a0 = A0[9];  c00 = a0 * b0; C00 += c00;
                b1 = B0[25]; c01 = a0 * b1; C01 += c01;
                b0 = B0[10];

                // 11
                a0 = A0[10]; c00 = a0 * b0; C00 += c00;
                b1 = B0[26]; c01 = a0 * b1; C01 += c01;
                b0 = B0[11];

                // 12
                a0 = A0[11]; c00 = a0 * b0; C00 += c00;
                b1 = B0[27]; c01 = a0 * b1; C01 += c01;
                b0 = B0[12];

                // 13
                a0 = A0[12]; c00 = a0 * b0; C00 += c00;
                b1 = B0[28]; c01 = a0 * b1; C01 += c01;
                b0 = B0[13];

                // 14
                a0 = A0[13]; c00 = a0 * b0; C00 += c00;
                b1 = B0[29]; c01 = a0 * b1; C01 += c01;
                b0 = B0[14];

                // 15
                a0 = A0[14]; c00 = a0 * b0; C00 += c00;
                b1 = B0[30]; c01 = a0 * b1; C01 += c01;
                b0 = B0[15];

                // 16
                a0 = A0[15]; c00 = a0 * b0; C00 += c00; A0 += 16;
                *C0++ += C00;
                b1 = B0[31]; c01 = a0 * b1; C01 += c01;
                *C1++ += C01;
                b0 = B0[0];
            } while (--i);
            A0 -= Mx16;
            B0 += 32;
            C0 += M;
            C1 += M;
        } while (--j);
    }
    template <typename T>
    static void generic_multmm_16_N_16(
        const ptrdiff_t N, const T* A0, const T* B0, T* C0)
    {
        const ptrdiff_t N_2 = N>>1; // N_2 = N/2
        const ptrdiff_t Nc = N-(N_2<<1); // Nc = N%2

        T a0, a1;
        T b0, b1;
        T C00, C01, C10, C11;
        T c00, c01, c10, c11;

        ptrdiff_t j,i;

        j = N_2; if (j) do {
            b0 = B0[0]; 
            b1 = B0[16]; 
            i = 4; do {
                // 1
                a0 = A0[0];  c00 = a0 * b0; C00 = c00;
                a1 = A0[16]; c10 = a1 * b0; C10 = c10;
                c01 = a0 * b1; C01 = c01;
                b0 = B0[1];  c11 = a1 * b1; C11 = c11;

                // 2
                a0 = A0[1];  c00 = a0 * b0; C00 += c00;
                a1 = A0[17]; c10 = a1 * b0; C10 += c10;
                b1 = B0[17]; c01 = a0 * b1; C01 += c01;
                b0 = B0[2];  c11 = a1 * b1; C11 += c11;

                // 3
                a0 = A0[2];  c00 = a0 * b0; C00 += c00;
                a1 = A0[18]; c10 = a1 * b0; C10 += c10;
                b1 = B0[18]; c01 = a0 * b1; C01 += c01;
                b0 = B0[3];  c11 = a1 * b1; C11 += c11;

                // 4
                a0 = A0[3];  c00 = a0 * b0; C00 += c00;
                a1 = A0[19]; c10 = a1 * b0; C10 += c10;
                b1 = B0[19]; c01 = a0 * b1; C01 += c01;
                b0 = B0[4];  c11 = a1 * b1; C11 += c11;

                // 5
                a0 = A0[4];  c00 = a0 * b0; C00 += c00;
                a1 = A0[20]; c10 = a1 * b0; C10 += c10;
                b1 = B0[20]; c01 = a0 * b1; C01 += c01;
                b0 = B0[5];  c11 = a1 * b1; C11 += c11;

                // 6
                a0 = A0[5];  c00 = a0 * b0; C00 += c00;
                a1 = A0[21]; c10 = a1 * b0; C10 += c10;
                b1 = B0[21]; c01 = a0 * b1; C01 += c01;
                b0 = B0[6];  c11 = a1 * b1; C11 += c11;

                // 7
                a0 = A0[6];  c00 = a0 * b0; C00 += c00;
                a1 = A0[22]; c10 = a1 * b0; C10 += c10;
                b1 = B0[22]; c01 = a0 * b1; C01 += c01;
                b0 = B0[7];  c11 = a1 * b1; C11 += c11;

                // 8
                a0 = A0[7];  c00 = a0 * b0; C00 += c00;
                a1 = A0[23]; c10 = a1 * b0; C10 += c10;
                b1 = B0[23]; c01 = a0 * b1; C01 += c01;
                b0 = B0[8];  c11 = a1 * b1; C11 += c11;

                // 9
                a0 = A0[8];  c00 = a0 * b0; C00 += c00;
                a1 = A0[24]; c10 = a1 * b0; C10 += c10;
                b1 = B0[24]; c01 = a0 * b1; C01 += c01;
                b0 = B0[9];  c11 = a1 * b1; C11 += c11;

                // 10
                a0 = A0[9];  c00 = a0 * b0; C00 += c00;
                a1 = A0[25]; c10 = a1 * b0; C10 += c10;
                b1 = B0[25]; c01 = a0 * b1; C01 += c01;
                b0 = B0[10]; c11 = a1 * b1; C11 += c11;

                // 11
                a0 = A0[10]; c00 = a0 * b0; C00 += c00;
                a1 = A0[26]; c10 = a1 * b0; C10 += c10;
                b1 = B0[26]; c01 = a0 * b1; C01 += c01;
                b0 = B0[11]; c11 = a1 * b1; C11 += c11;

                // 12
                a0 = A0[11]; c00 = a0 * b0; C00 += c00;
                a1 = A0[27]; c10 = a1 * b0; C10 += c10;
                b1 = B0[27]; c01 = a0 * b1; C01 += c01;
                b0 = B0[12]; c11 = a1 * b1; C11 += c11;

                // 13
                a0 = A0[12]; c00 = a0 * b0; C00 += c00;
                a1 = A0[28]; c10 = a1 * b0; C10 += c10;
                b1 = B0[28]; c01 = a0 * b1; C01 += c01;
                b0 = B0[13]; c11 = a1 * b1; C11 += c11;

                // 14
                a0 = A0[13]; c00 = a0 * b0; C00 += c00;
                a1 = A0[29]; c10 = a1 * b0; C10 += c10;
                b1 = B0[29]; c01 = a0 * b1; C01 += c01;
                b0 = B0[14]; c11 = a1 * b1; C11 += c11;

                // 15
                a0 = A0[14]; c00 = a0 * b0; C00 += c00;
                a1 = A0[30]; c10 = a1 * b0; C10 += c10;
                b1 = B0[30]; c01 = a0 * b1; C01 += c01;
                b0 = B0[15]; c11 = a1 * b1; C11 += c11;

                // 16
                a0 = A0[15]; c00 = a0 * b0; C00 += c00;
                C0[0] += C00;
                a1 = A0[31]; c10 = a1 * b0; C10 += c10;
                C0[1] += C10;
                b1 = B0[31]; c01 = a0 * b1; C01 += c01;
                C0[16] += C01;
                c11 = a1 * b1; C11 += c11;
                C0[17] += C11;

                // 16
                a0 = A0[47]; c00 = a0 * b0; C00 = c00;
                a1 = A0[63]; c10 = a1 * b0; C10 = c10;
                c01 = a0 * b1; C01 = c01;
                b0 = B0[14]; c11 = a1 * b1; C11 = c11;

                // 15
                a0 = A0[46]; c00 = a0 * b0; C00 += c00;
                a1 = A0[62]; c10 = a1 * b0; C10 += c10;
                b1 = B0[30]; c01 = a0 * b1; C01 += c01;
                b0 = B0[13]; c11 = a1 * b1; C11 += c11;

                // 14
                a0 = A0[45]; c00 = a0 * b0; C00 += c00;
                a1 = A0[61]; c10 = a1 * b0; C10 += c10;
                b1 = B0[29]; c01 = a0 * b1; C01 += c01;
                b0 = B0[12]; c11 = a1 * b1; C11 += c11;

                // 13
                a0 = A0[44]; c00 = a0 * b0; C00 += c00;
                a1 = A0[60]; c10 = a1 * b0; C10 += c10;
                b1 = B0[28]; c01 = a0 * b1; C01 += c01;
                b0 = B0[11]; c11 = a1 * b1; C11 += c11;

                // 12
                a0 = A0[43]; c00 = a0 * b0; C00 += c00;
                a1 = A0[59]; c10 = a1 * b0; C10 += c10;
                b1 = B0[27]; c01 = a0 * b1; C01 += c01;
                b0 = B0[10]; c11 = a1 * b1; C11 += c11;

                // 11
                a0 = A0[42]; c00 = a0 * b0; C00 += c00;
                a1 = A0[58]; c10 = a1 * b0; C10 += c10;
                b1 = B0[26]; c01 = a0 * b1; C01 += c01;
                b0 = B0[9];  c11 = a1 * b1; C11 += c11;

                // 10
                a0 = A0[41]; c00 = a0 * b0; C00 += c00;
                a1 = A0[57]; c10 = a1 * b0; C10 += c10;
                b1 = B0[25]; c01 = a0 * b1; C01 += c01;
                b0 = B0[8];  c11 = a1 * b1; C11 += c11;

                // 9
                a0 = A0[40]; c00 = a0 * b0; C00 += c00;
                a1 = A0[56]; c10 = a1 * b0; C10 += c10;
                b1 = B0[24]; c01 = a0 * b1; C01 += c01;
                b0 = B0[7];  c11 = a1 * b1; C11 += c11;

                // 8
                a0 = A0[39]; c00 = a0 * b0; C00 += c00;
                a1 = A0[55]; c10 = a1 * b0; C10 += c10;
                b1 = B0[23]; c01 = a0 * b1; C01 += c01;
                b0 = B0[6];  c11 = a1 * b1; C11 += c11;

                // 7
                a0 = A0[38]; c00 = a0 * b0; C00 += c00;
                a1 = A0[54]; c10 = a1 * b0; C10 += c10;
                b1 = B0[22]; c01 = a0 * b1; C01 += c01;
                b0 = B0[5];  c11 = a1 * b1; C11 += c11;

                // 6
                a0 = A0[37]; c00 = a0 * b0; C00 += c00;
                a1 = A0[53]; c10 = a1 * b0; C10 += c10;
                b1 = B0[21]; c01 = a0 * b1; C01 += c01;
                b0 = B0[4];  c11 = a1 * b1; C11 += c11;

                // 5
                a0 = A0[36]; c00 = a0 * b0; C00 += c00;
                a1 = A0[52]; c10 = a1 * b0; C10 += c10;
                b1 = B0[20]; c01 = a0 * b1; C01 += c01;
                b0 = B0[3];  c11 = a1 * b1; C11 += c11;

                // 4
                a0 = A0[35]; c00 = a0 * b0; C00 += c00;
                a1 = A0[51]; c10 = a1 * b0; C10 += c10;
                b1 = B0[19]; c01 = a0 * b1; C01 += c01;
                b0 = B0[2];  c11 = a1 * b1; C11 += c11;

                // 3
                a0 = A0[34]; c00 = a0 * b0; C00 += c00;
                a1 = A0[50]; c10 = a1 * b0; C10 += c10;
                b1 = B0[18]; c01 = a0 * b1; C01 += c01;
                b0 = B0[1];  c11 = a1 * b1; C11 += c11;

                // 2
                a0 = A0[33]; c00 = a0 * b0; C00 += c00;
                a1 = A0[49]; c10 = a1 * b0; C10 += c10;
                b1 = B0[17]; c01 = a0 * b1; C01 += c01;
                b0 = B0[0];  c11 = a1 * b1; C11 += c11;

                // 1
                a0 = A0[32]; c00 = a0 * b0; C00 += c00;
                C0[ 2] += C00;
                a1 = A0[48]; c10 = a1 * b0; C10 += c10;
                C0[ 3] += C10;
                b1 = B0[16]; c01 = a0 * b1; C01 += c01;
                C0[18] += C01;
                A0 += 64;    c11 = a1 * b1; C11 += c11;
                C0[19] += C11;
                C0 += 4;
            } while (--i);
            A0 -= 256;
            B0 += 32;
            C0 += 16;
        } while (--j);
        if (Nc) {
            b0 = B0[0]; 
            i = 4; do {
                // 1
                a0 = A0[0];  c00 = a0 * b0; C00 = c00;
                a1 = A0[16]; c10 = a1 * b0; C10 = c10;
                b0 = B0[1];  

                // 2
                a0 = A0[1];  c00 = a0 * b0; C00 += c00;
                a1 = A0[17]; c10 = a1 * b0; C10 += c10;
                b0 = B0[2];  

                // 3
                a0 = A0[2];  c00 = a0 * b0; C00 += c00;
                a1 = A0[18]; c10 = a1 * b0; C10 += c10;
                b0 = B0[3];  

                // 4
                a0 = A0[3];  c00 = a0 * b0; C00 += c00;
                a1 = A0[19]; c10 = a1 * b0; C10 += c10;
                b0 = B0[4];  

                // 5
                a0 = A0[4];  c00 = a0 * b0; C00 += c00;
                a1 = A0[20]; c10 = a1 * b0; C10 += c10;
                b0 = B0[5];  

                // 6
                a0 = A0[5];  c00 = a0 * b0; C00 += c00;
                a1 = A0[21]; c10 = a1 * b0; C10 += c10;
                b0 = B0[6];  

                // 7
                a0 = A0[6];  c00 = a0 * b0; C00 += c00;
                a1 = A0[22]; c10 = a1 * b0; C10 += c10;
                b0 = B0[7];  

                // 8
                a0 = A0[7];  c00 = a0 * b0; C00 += c00;
                a1 = A0[23]; c10 = a1 * b0; C10 += c10;
                b0 = B0[8];  

                // 9
                a0 = A0[8];  c00 = a0 * b0; C00 += c00;
                a1 = A0[24]; c10 = a1 * b0; C10 += c10;
                b0 = B0[9];  

                // 10
                a0 = A0[9];  c00 = a0 * b0; C00 += c00;
                a1 = A0[25]; c10 = a1 * b0; C10 += c10;
                b0 = B0[10]; 

                // 11
                a0 = A0[10]; c00 = a0 * b0; C00 += c00;
                a1 = A0[26]; c10 = a1 * b0; C10 += c10;
                b0 = B0[11]; 

                // 12
                a0 = A0[11]; c00 = a0 * b0; C00 += c00;
                a1 = A0[27]; c10 = a1 * b0; C10 += c10;
                b0 = B0[12]; 

                // 13
                a0 = A0[12]; c00 = a0 * b0; C00 += c00;
                a1 = A0[28]; c10 = a1 * b0; C10 += c10;
                b0 = B0[13]; 

                // 14
                a0 = A0[13]; c00 = a0 * b0; C00 += c00;
                a1 = A0[29]; c10 = a1 * b0; C10 += c10;
                b0 = B0[14]; 

                // 15
                a0 = A0[14]; c00 = a0 * b0; C00 += c00;
                a1 = A0[30]; c10 = a1 * b0; C10 += c10;
                b0 = B0[15]; 

                // 16
                a0 = A0[15]; c00 = a0 * b0; C00 += c00;
                C0[0] += C00;
                a1 = A0[31]; c10 = a1 * b0; C10 += c10;
                C0[1] += C10;

                // 16
                a0 = A0[47]; c00 = a0 * b0; C00 = c00;
                a1 = A0[63]; c10 = a1 * b0; C10 = c10;
                b0 = B0[14]; 

                // 15
                a0 = A0[46]; c00 = a0 * b0; C00 += c00;
                a1 = A0[62]; c10 = a1 * b0; C10 += c10;
                b0 = B0[13]; 

                // 14
                a0 = A0[45]; c00 = a0 * b0; C00 += c00;
                a1 = A0[61]; c10 = a1 * b0; C10 += c10;
                b0 = B0[12]; 

                // 13
                a0 = A0[44]; c00 = a0 * b0; C00 += c00;
                a1 = A0[60]; c10 = a1 * b0; C10 += c10;
                b0 = B0[11]; 

                // 12
                a0 = A0[43]; c00 = a0 * b0; C00 += c00;
                a1 = A0[59]; c10 = a1 * b0; C10 += c10;
                b0 = B0[10]; 

                // 11
                a0 = A0[42]; c00 = a0 * b0; C00 += c00;
                a1 = A0[58]; c10 = a1 * b0; C10 += c10;
                b0 = B0[9];  

                // 10
                a0 = A0[41]; c00 = a0 * b0; C00 += c00;
                a1 = A0[57]; c10 = a1 * b0; C10 += c10;
                b0 = B0[8];  

                // 9
                a0 = A0[40]; c00 = a0 * b0; C00 += c00;
                a1 = A0[56]; c10 = a1 * b0; C10 += c10;
                b0 = B0[7];  

                // 8
                a0 = A0[39]; c00 = a0 * b0; C00 += c00;
                a1 = A0[55]; c10 = a1 * b0; C10 += c10;
                b0 = B0[6];  

                // 7
                a0 = A0[38]; c00 = a0 * b0; C00 += c00;
                a1 = A0[54]; c10 = a1 * b0; C10 += c10;
                b0 = B0[5];  

                // 6
                a0 = A0[37]; c00 = a0 * b0; C00 += c00;
                a1 = A0[53]; c10 = a1 * b0; C10 += c10;
                b0 = B0[4];  

                // 5
                a0 = A0[36]; c00 = a0 * b0; C00 += c00;
                a1 = A0[52]; c10 = a1 * b0; C10 += c10;
                b0 = B0[3];  

                // 4
                a0 = A0[35]; c00 = a0 * b0; C00 += c00;
                a1 = A0[51]; c10 = a1 * b0; C10 += c10;
                b0 = B0[2];  

                // 3
                a0 = A0[34]; c00 = a0 * b0; C00 += c00;
                a1 = A0[50]; c10 = a1 * b0; C10 += c10;
                b0 = B0[1];  

                // 2
                a0 = A0[33]; c00 = a0 * b0; C00 += c00;
                a1 = A0[49]; c10 = a1 * b0; C10 += c10;
                b0 = B0[0];  

                // 1
                a0 = A0[32]; c00 = a0 * b0; C00 += c00;
                C0[ 2] += C00;
                a1 = A0[48]; c10 = a1 * b0; C10 += c10;
                C0[ 3] += C10;
                A0 += 64;    
                C0 += 4;
            } while (--i);
        }
    }
    template <typename T>
    static void generic_multmm_16_16_K(
        const ptrdiff_t K, const T* A0, const T* B0, T* C0)
    {
        TMVAssert(K > 0);
        const ptrdiff_t K_4 = (K-1)>>2; // K_4 = (K-1)/4 
        // The K-1 here  ^^^ is to make sure Kc > 0, since the regular
        // loop accesses memory past the end, so we use the Kc loop
        // to make sure that doesn't lead to a seg fault.
        const ptrdiff_t Kc = K-(K_4<<2); // = K % 4
        const ptrdiff_t Kcm1 = Kc-1;
        const ptrdiff_t Kx16 = K<<4; // = K * 16
        const ptrdiff_t Kx2 = K<<1; // = K * 2

        T a0, a1;
        T b0, b1;
        T C00, C01, C10, C11;
        T c00, c01, c10, c11;

        const T* A1 = A0 + K;
        const T* B1 = B0 + K;
        T* C1 = C0 + 16;

        ptrdiff_t j,i,k;

        j = 8; do {
            i = 8; do {
                C00 = T(0);
                C01 = T(0);
                C10 = T(0);
                C11 = T(0);
                b0 = B0[0];
                k = K_4; if (k) do {
                    a0 = A0[0]; c00 = a0 * b0; C00 += c00;
                    a1 = A1[0]; c10 = a1 * b0; C10 += c10;
                    b1 = B1[0]; c01 = a0 * b1; C01 += c01;
                    b0 = B0[1]; c11 = a1 * b1; C11 += c11;
                    a0 = A0[1]; c00 = a0 * b0; C00 += c00;
                    a1 = A1[1]; c10 = a1 * b0; C10 += c10;
                    b1 = B1[1]; c01 = a0 * b1; C01 += c01;
                    b0 = B0[2]; c11 = a1 * b1; C11 += c11;
                    a0 = A0[2]; c00 = a0 * b0; C00 += c00;
                    a1 = A1[2]; c10 = a1 * b0; C10 += c10;
                    b1 = B1[2]; c01 = a0 * b1; C01 += c01;
                    b0 = B0[3]; c11 = a1 * b1; C11 += c11;
                    a0 = A0[3]; c00 = a0 * b0; C00 += c00; A0 += 4;
                    a1 = A1[3]; c10 = a1 * b0; C10 += c10; A1 += 4;
                    b1 = B1[3]; c01 = a0 * b1; C01 += c01; B1 += 4;
                    b0 = B0[4]; c11 = a1 * b1; C11 += c11; B0 += 4;
                } while (--k);
                k = Kcm1; if (k) do {
                    a0 = *A0++; c00 = a0 * b0; C00 += c00;
                    a1 = *A1++; c10 = a1 * b0; C10 += c10;
                    b1 = *B1++; c01 = a0 * b1; C01 += c01;
                    b0 = *++B0; c11 = a1 * b1; C11 += c11;
                } while (--k);
                a0 = *A0++; c00 = a0 * b0; C00 += c00;
                a1 = *A1++; c10 = a1 * b0; C10 += c10;
                b1 = *B1++; c01 = a0 * b1; C01 += c01;
                ++B0;       c11 = a1 * b1; C11 += c11;
                C0[0] += C00;
                C0[1] += C10; C0 += 2;
                C1[0] += C01;
                C1[1] += C11; C1 += 2;
                A0 += K; A1 += K;
                B0 -= K; B1 -= K;
            } while (--i);
            A0 -= Kx16; A1 -= Kx16;
            B0 += Kx2; B1 += Kx2;
            C0 += 16; C1 += 16;
        } while (--j);
    }
    template <typename T>
    static void generic_multmm_M_N_K(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K, const ptrdiff_t Kd,
        const T* A0, const T* B0, T* C0)
    {
        const ptrdiff_t M_2 = M>>1; // M_2 = M/2
        const ptrdiff_t Mb = M_2<<1; // Mb = 2*(M/2)
        const ptrdiff_t Mc = M-Mb; // Mc = M%2

        const ptrdiff_t N_2 = N>>1; // N_2 = N/2
        const ptrdiff_t Nc = N-(N_2<<1); // Nc = N%2

        const ptrdiff_t K_4 = (K-1)>>2; // K_4 = (K-1)/4 
        // The K-1 here  ^^^ is to make sure Kc > 0, since the regular
        // loop accesses memory past the end, so we use the Kc loop
        // to make sure that doesn't lead to a seg fault.
        const ptrdiff_t Kc = K-(K_4<<2); // = K % 4
        const ptrdiff_t Kcm1 = Kc-1;
        const ptrdiff_t KdxMb = Kd*Mb; // = Kd * Mb
        const ptrdiff_t Kdx2 = Kd<<1; // = Kd * 2
        const ptrdiff_t Kdx2mK = Kdx2-K; // = Kd * 2 - K


        T a0, a1;
        T b0, b1;
        T C00, C01, C10, C11;
        T c00, c01, c10, c11;

        const T* A1 = A0 + Kd;
        const T* B1 = B0 + Kd;
        T* C1 = C0 + M;

        ptrdiff_t j,i,k;

        j = N_2; if (j) do {
            i = M_2; if (i) do {
                b0 = *B0++;
                C00 = T(0);
                C01 = T(0);
                C10 = T(0);
                C11 = T(0);
                k = K_4; if (k) do {
                    a0 = A0[0]; c00 = a0 * b0; C00 += c00;
                    a1 = A1[0]; c10 = a1 * b0; C10 += c10;
                    b1 = B1[0]; c01 = a0 * b1; C01 += c01;
                    b0 = B0[0]; c11 = a1 * b1; C11 += c11;
                    a0 = A0[1]; c00 = a0 * b0; C00 += c00;
                    a1 = A1[1]; c10 = a1 * b0; C10 += c10;
                    b1 = B1[1]; c01 = a0 * b1; C01 += c01;
                    b0 = B0[1]; c11 = a1 * b1; C11 += c11;
                    a0 = A0[2]; c00 = a0 * b0; C00 += c00;
                    a1 = A1[2]; c10 = a1 * b0; C10 += c10;
                    b1 = B1[2]; c01 = a0 * b1; C01 += c01;
                    b0 = B0[2]; c11 = a1 * b1; C11 += c11;
                    a0 = A0[3]; c00 = a0 * b0; C00 += c00; A0 += 4;
                    a1 = A1[3]; c10 = a1 * b0; C10 += c10; A1 += 4;
                    b1 = B1[3]; c01 = a0 * b1; C01 += c01; B1 += 4;
                    b0 = B0[3]; c11 = a1 * b1; C11 += c11; B0 += 4;
                } while (--k);
                k = Kcm1; if (k) do {
                    a0 = *A0++; c00 = a0 * b0; C00 += c00;
                    a1 = *A1++; c10 = a1 * b0; C10 += c10;
                    b1 = *B1++; c01 = a0 * b1; C01 += c01;
                    b0 = *B0++; c11 = a1 * b1; C11 += c11;
                } while (--k);
                a0 = *A0++; c00 = a0 * b0; C00 += c00;
                a1 = *A1++; c10 = a1 * b0; C10 += c10;
                b1 = *B1++; c01 = a0 * b1; C01 += c01;
                c11 = a1 * b1; C11 += c11;
                C0[0] += C00;
                C0[1] += C10; C0 += 2;
                C1[0] += C01;
                C1[1] += C11; C1 += 2;
                A0 += Kdx2mK; A1 += Kdx2mK;
                B0 -= K; B1 -= K;
            } while (--i);
            if (Mc) {
                b0 = *B0++;
                C00 = T(0);
                C01 = T(0);
                k = K_4; if (k) do {
                    a0 = A0[0]; c00 = a0 * b0; C00 += c00;
                    b1 = B1[0]; c01 = a0 * b1; C01 += c01;
                    b0 = B0[0]; 
                    a0 = A0[1]; c00 = a0 * b0; C00 += c00;
                    b1 = B1[1]; c01 = a0 * b1; C01 += c01;
                    b0 = B0[1]; 
                    a0 = A0[2]; c00 = a0 * b0; C00 += c00;
                    b1 = B1[2]; c01 = a0 * b1; C01 += c01;
                    b0 = B0[2]; 
                    a0 = A0[3]; c00 = a0 * b0; C00 += c00; A0 += 4;
                    b1 = B1[3]; c01 = a0 * b1; C01 += c01; B1 += 4;
                    b0 = B0[3]; B0 += 4;
                } while (--k);
                k = Kcm1; if (k) do {
                    a0 = *A0++; c00 = a0 * b0; C00 += c00;
                    b1 = *B1++; c01 = a0 * b1; C01 += c01;
                    b0 = *B0++;
                } while (--k);
                a0 = *A0++; c00 = a0 * b0; C00 += c00;
                b1 = *B1++; c01 = a0 * b1; C01 += c01;
                *C0++ += C00;
                *C1++ += C01;
                A0 -= K;
                B0 -= K; B1 -= K;
            }
            A0 -= KdxMb; A1 -= KdxMb;
            B0 += Kdx2; B1 += Kdx2;
            C0 += M; C1 += M;
        } while (--j);
        if (Nc) {
            i = M_2; if (i) do {
                b0 = *B0++;
                C00 = T(0);
                C10 = T(0);
                k = K_4; if (k) do {
                    a0 = A0[0]; c00 = a0 * b0; C00 += c00;
                    a1 = A1[0]; c10 = a1 * b0; C10 += c10;
                    b0 = B0[0]; 
                    a0 = A0[1]; c00 = a0 * b0; C00 += c00;
                    a1 = A1[1]; c10 = a1 * b0; C10 += c10;
                    b0 = B0[1]; 
                    a0 = A0[2]; c00 = a0 * b0; C00 += c00;
                    a1 = A1[2]; c10 = a1 * b0; C10 += c10;
                    b0 = B0[2]; 
                    a0 = A0[3]; c00 = a0 * b0; C00 += c00; A0 += 4;
                    a1 = A1[3]; c10 = a1 * b0; C10 += c10; A1 += 4;
                    b0 = B0[3]; B0 += 4;
                } while (--k);
                k = Kcm1; if (k) do {
                    a0 = *A0++; c00 = a0 * b0; C00 += c00;
                    a1 = *A1++; c10 = a1 * b0; C10 += c10;
                    b0 = *B0++; 
                } while (--k);
                a0 = *A0++; c00 = a0 * b0; C00 += c00;
                a1 = *A1++; c10 = a1 * b0; C10 += c10;
                C0[0] += C00;
                C0[1] += C10; C0 += 2;
                A0 += Kdx2mK; A1 += Kdx2mK;
                B0 -= K;
            } while (--i);
            if (Mc) {
                b0 = *B0++;
                C00 = T(0);
                k = K_4; if (k) do {
                    a0 = A0[0]; c00 = a0 * b0; C00 += c00;
                    b0 = B0[0]; 
                    a0 = A0[1]; c00 = a0 * b0; C00 += c00;
                    b0 = B0[1]; 
                    a0 = A0[2]; c00 = a0 * b0; C00 += c00;
                    b0 = B0[2]; 
                    a0 = A0[3]; c00 = a0 * b0; C00 += c00; A0 += 4;
                    b0 = B0[3]; B0 += 4;
                } while (--k);
                k = Kcm1; if (k) do {
                    a0 = *A0++; c00 = a0 * b0; C00 += c00;
                    b0 = *B0++;
                } while (--k);
                a0 = *A0++; c00 = a0 * b0; C00 += c00;
                *C0++ += C00;
            }
        }
    }

    // First the ones that that just call the generic functions.
    // There are always inlined, but get overloaded below
    // as necessary depending on the INST_* definitions.
    template <typename T>
    static void multmm_16_16_16(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t ,
        const T* A0, const T* B0, T* C0)
    { generic_multmm_16_16_16(A0,B0,C0); }

    template <typename T>
    static void multmm_M_16_16(
        const ptrdiff_t M, const ptrdiff_t , const ptrdiff_t ,
        const T* A0, const T* B0, T* C0)
    { generic_multmm_M_16_16(M,A0,B0,C0); }

    template <typename T>
    static void multmm_16_N_16(
        const ptrdiff_t , const ptrdiff_t N, const ptrdiff_t ,
        const T* A0, const T* B0, T* C0)
    { generic_multmm_16_N_16(N,A0,B0,C0); }

    template <typename T>
    static void multmm_16_16_K(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t K,
        const T* A0, const T* B0, T* C0)
    { generic_multmm_16_16_K(K,A0,B0,C0); }

    template <typename T>
    static void multmm_M_N_K(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const T* A0, const T* B0, T* C0)
    { generic_multmm_M_N_K(M,N,K,K,A0,B0,C0); }


#ifndef BLAS
#ifdef INST_FLOAT
    //
    // float
    //

#ifdef __SSE__
    // If we have SSE commands available, we have special versions for float
    static void sse_multmm_16_16_64(
        const float* A, const float* B, float* C0)
    {
        TMVAssert( ((size_t)(A) & 0xf) == 0 );
        TMVAssert( ((size_t)(B) & 0xf) == 0 );
        TMVAssert( ((size_t)(C0) & 0xf) == 0 );

        __m128 C00, C01, C10, C11;
        __m128 c00, c01, c10, c11;
        union { __m128 xm; float xf[4]; } d00, d01, d10, d11;

        const __m128* A0 = (const __m128*) A;
        const __m128* B0 = (const __m128*) B;

        ptrdiff_t i2,j2;

        j2 = 8; do {
            i2 = 4; do {
                C00 = _mm_mul_ps(A0[ 0],B0[ 0]); 
                C10 = _mm_mul_ps(A0[16],B0[ 0]); 
                C01 = _mm_mul_ps(A0[ 0],B0[16]); 
                C11 = _mm_mul_ps(A0[16],B0[16]); 

                c00 = _mm_mul_ps(A0[ 1],B0[ 1]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[17],B0[ 1]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 1],B0[17]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[17],B0[17]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 2],B0[ 2]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[18],B0[ 2]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 2],B0[18]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[18],B0[18]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 3],B0[ 3]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[19],B0[ 3]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 3],B0[19]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[19],B0[19]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 4],B0[ 4]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[20],B0[ 4]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 4],B0[20]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[20],B0[20]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 5],B0[ 5]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[21],B0[ 5]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 5],B0[21]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[21],B0[21]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 6],B0[ 6]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[22],B0[ 6]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 6],B0[22]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[22],B0[22]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 7],B0[ 7]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[23],B0[ 7]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 7],B0[23]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[23],B0[23]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 8],B0[ 8]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[24],B0[ 8]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 8],B0[24]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[24],B0[24]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 9],B0[ 9]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[25],B0[ 9]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 9],B0[25]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[25],B0[25]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[10],B0[10]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[26],B0[10]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[10],B0[26]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[26],B0[26]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[11],B0[11]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[27],B0[11]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[11],B0[27]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[27],B0[27]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[12],B0[12]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[28],B0[12]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[12],B0[28]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[28],B0[28]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[13],B0[13]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[29],B0[13]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[13],B0[29]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[29],B0[29]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[14],B0[14]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[30],B0[14]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[14],B0[30]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[30],B0[30]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[15],B0[15]); d00.xm = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[31],B0[15]); d10.xm = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[15],B0[31]); d01.xm = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[31],B0[31]); d11.xm = _mm_add_ps(C11,c11);

                C0[ 0] += d00.xf[0] + d00.xf[1] + d00.xf[2] + d00.xf[3];
                C00 = _mm_mul_ps(A0[47],B0[15]);
                C10 = _mm_mul_ps(A0[63],B0[15]);
                C01 = _mm_mul_ps(A0[47],B0[31]);
                C11 = _mm_mul_ps(A0[63],B0[31]);

                C0[ 1] += d10.xf[0] + d10.xf[1] + d10.xf[2] + d10.xf[3];
                c00 = _mm_mul_ps(A0[46],B0[14]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[62],B0[14]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[46],B0[30]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[62],B0[30]); C11 = _mm_add_ps(C11,c11);

                C0[16] += d01.xf[0] + d01.xf[1] + d01.xf[2] + d01.xf[3];
                c00 = _mm_mul_ps(A0[45],B0[13]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[61],B0[13]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[45],B0[29]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[61],B0[29]); C11 = _mm_add_ps(C11,c11);

                C0[17] += d11.xf[0] + d11.xf[1] + d11.xf[2] + d11.xf[3];
                c00 = _mm_mul_ps(A0[44],B0[12]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[60],B0[12]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[44],B0[28]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[60],B0[28]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[43],B0[11]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[59],B0[11]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[43],B0[27]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[59],B0[27]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[42],B0[10]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[58],B0[10]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[42],B0[26]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[58],B0[26]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[41],B0[ 9]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[57],B0[ 9]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[41],B0[25]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[57],B0[25]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[40],B0[ 8]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[56],B0[ 8]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[40],B0[24]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[56],B0[24]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[39],B0[ 7]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[55],B0[ 7]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[39],B0[23]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[55],B0[23]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[38],B0[ 6]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[54],B0[ 6]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[38],B0[22]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[54],B0[22]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[37],B0[ 5]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[53],B0[ 5]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[37],B0[21]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[53],B0[21]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[36],B0[ 4]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[52],B0[ 4]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[36],B0[20]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[52],B0[20]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[35],B0[ 3]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[51],B0[ 3]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[35],B0[19]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[51],B0[19]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[34],B0[ 2]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[50],B0[ 2]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[34],B0[18]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[50],B0[18]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[33],B0[ 1]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[49],B0[ 1]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[33],B0[17]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[49],B0[17]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[32],B0[ 0]); d00.xm = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[48],B0[ 0]); d10.xm = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[32],B0[16]); d01.xm = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[48],B0[16]); d11.xm = _mm_add_ps(C11,c11);

                C0[ 2] += d00.xf[0] + d00.xf[1] + d00.xf[2] + d00.xf[3];
                C0[ 3] += d10.xf[0] + d10.xf[1] + d10.xf[2] + d10.xf[3];
                C0[18] += d01.xf[0] + d01.xf[1] + d01.xf[2] + d01.xf[3];
                C0[19] += d11.xf[0] + d11.xf[1] + d11.xf[2] + d11.xf[3];

                A0 += 64;
                C0 += 4;
            } while (--i2);
            A0 -= 256;
            B0 += 32;
            C0 += 16;
        } while (--j2);
    }
    static void sse_multmm_M_16_64(
        const ptrdiff_t M, const float* A, const float* B, float* C0)
    {
        TMVAssert( ((size_t)(A) & 0xf) == 0 );
        TMVAssert( ((size_t)(B) & 0xf) == 0 );
        TMVAssert( ((size_t)(C0) & 0xf) == 0 );
        const ptrdiff_t Mx16 = (M<<4);      // = M*16
        const ptrdiff_t M_4 = M>>2;        // = M/4
        const ptrdiff_t Mc = (M-(M_4<<2)); // = M%4
        const ptrdiff_t Mc1 = Mc>>1;       // = M%4 == 2 or 3
        const ptrdiff_t Mc2 = Mc-(Mc1<<1); // = M%4 == 1 or 3

        __m128 C00, C01, C10, C11;
        __m128 c00, c01, c10, c11;
        union { __m128 xm; float xf[4]; } d00, d01, d10, d11;

        const __m128* A0 = (const __m128*) A;
        const __m128* B0 = (const __m128*) B;
        float* C1 = C0 + M;

        ptrdiff_t i2,j2;

        j2 = 8; do {
            i2 = M_4; if (i2) do {
                C00 = _mm_mul_ps(A0[ 0],B0[ 0]); 
                C10 = _mm_mul_ps(A0[16],B0[ 0]); 
                C01 = _mm_mul_ps(A0[ 0],B0[16]); 
                C11 = _mm_mul_ps(A0[16],B0[16]); 

                c00 = _mm_mul_ps(A0[ 1],B0[ 1]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[17],B0[ 1]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 1],B0[17]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[17],B0[17]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 2],B0[ 2]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[18],B0[ 2]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 2],B0[18]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[18],B0[18]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 3],B0[ 3]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[19],B0[ 3]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 3],B0[19]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[19],B0[19]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 4],B0[ 4]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[20],B0[ 4]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 4],B0[20]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[20],B0[20]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 5],B0[ 5]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[21],B0[ 5]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 5],B0[21]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[21],B0[21]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 6],B0[ 6]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[22],B0[ 6]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 6],B0[22]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[22],B0[22]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 7],B0[ 7]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[23],B0[ 7]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 7],B0[23]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[23],B0[23]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 8],B0[ 8]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[24],B0[ 8]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 8],B0[24]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[24],B0[24]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 9],B0[ 9]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[25],B0[ 9]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 9],B0[25]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[25],B0[25]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[10],B0[10]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[26],B0[10]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[10],B0[26]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[26],B0[26]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[11],B0[11]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[27],B0[11]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[11],B0[27]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[27],B0[27]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[12],B0[12]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[28],B0[12]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[12],B0[28]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[28],B0[28]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[13],B0[13]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[29],B0[13]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[13],B0[29]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[29],B0[29]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[14],B0[14]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[30],B0[14]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[14],B0[30]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[30],B0[30]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[15],B0[15]); d00.xm = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[31],B0[15]); d10.xm = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[15],B0[31]); d01.xm = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[31],B0[31]); d11.xm = _mm_add_ps(C11,c11);

                C0[ 0] += d00.xf[0] + d00.xf[1] + d00.xf[2] + d00.xf[3];
                C00 = _mm_mul_ps(A0[47],B0[15]);
                C10 = _mm_mul_ps(A0[63],B0[15]);
                C01 = _mm_mul_ps(A0[47],B0[31]);
                C11 = _mm_mul_ps(A0[63],B0[31]);

                C0[ 1] += d10.xf[0] + d10.xf[1] + d10.xf[2] + d10.xf[3];
                c00 = _mm_mul_ps(A0[46],B0[14]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[62],B0[14]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[46],B0[30]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[62],B0[30]); C11 = _mm_add_ps(C11,c11);

                C1[ 0] += d01.xf[0] + d01.xf[1] + d01.xf[2] + d01.xf[3];
                c00 = _mm_mul_ps(A0[45],B0[13]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[61],B0[13]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[45],B0[29]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[61],B0[29]); C11 = _mm_add_ps(C11,c11);

                C1[ 1] += d11.xf[0] + d11.xf[1] + d11.xf[2] + d11.xf[3];
                c00 = _mm_mul_ps(A0[44],B0[12]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[60],B0[12]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[44],B0[28]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[60],B0[28]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[43],B0[11]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[59],B0[11]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[43],B0[27]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[59],B0[27]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[42],B0[10]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[58],B0[10]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[42],B0[26]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[58],B0[26]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[41],B0[ 9]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[57],B0[ 9]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[41],B0[25]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[57],B0[25]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[40],B0[ 8]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[56],B0[ 8]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[40],B0[24]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[56],B0[24]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[39],B0[ 7]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[55],B0[ 7]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[39],B0[23]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[55],B0[23]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[38],B0[ 6]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[54],B0[ 6]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[38],B0[22]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[54],B0[22]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[37],B0[ 5]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[53],B0[ 5]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[37],B0[21]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[53],B0[21]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[36],B0[ 4]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[52],B0[ 4]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[36],B0[20]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[52],B0[20]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[35],B0[ 3]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[51],B0[ 3]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[35],B0[19]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[51],B0[19]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[34],B0[ 2]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[50],B0[ 2]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[34],B0[18]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[50],B0[18]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[33],B0[ 1]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[49],B0[ 1]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[33],B0[17]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[49],B0[17]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[32],B0[ 0]); d00.xm = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[48],B0[ 0]); d10.xm = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[32],B0[16]); d01.xm = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[48],B0[16]); d11.xm = _mm_add_ps(C11,c11);

                C0[ 2] += d00.xf[0] + d00.xf[1] + d00.xf[2] + d00.xf[3];
                C0[ 3] += d10.xf[0] + d10.xf[1] + d10.xf[2] + d10.xf[3];
                C1[ 2] += d01.xf[0] + d01.xf[1] + d01.xf[2] + d01.xf[3];
                C1[ 3] += d11.xf[0] + d11.xf[1] + d11.xf[2] + d11.xf[3];

                A0 += 64;
                C0 += 4;
                C1 += 4;
            } while (--i2);
            if (Mc1) {
                C00 = _mm_mul_ps(A0[ 0],B0[ 0]); 
                C10 = _mm_mul_ps(A0[16],B0[ 0]); 
                C01 = _mm_mul_ps(A0[ 0],B0[16]); 
                C11 = _mm_mul_ps(A0[16],B0[16]); 

                c00 = _mm_mul_ps(A0[ 1],B0[ 1]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[17],B0[ 1]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 1],B0[17]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[17],B0[17]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 2],B0[ 2]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[18],B0[ 2]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 2],B0[18]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[18],B0[18]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 3],B0[ 3]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[19],B0[ 3]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 3],B0[19]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[19],B0[19]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 4],B0[ 4]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[20],B0[ 4]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 4],B0[20]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[20],B0[20]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 5],B0[ 5]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[21],B0[ 5]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 5],B0[21]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[21],B0[21]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 6],B0[ 6]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[22],B0[ 6]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 6],B0[22]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[22],B0[22]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 7],B0[ 7]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[23],B0[ 7]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 7],B0[23]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[23],B0[23]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 8],B0[ 8]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[24],B0[ 8]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 8],B0[24]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[24],B0[24]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 9],B0[ 9]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[25],B0[ 9]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 9],B0[25]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[25],B0[25]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[10],B0[10]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[26],B0[10]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[10],B0[26]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[26],B0[26]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[11],B0[11]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[27],B0[11]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[11],B0[27]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[27],B0[27]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[12],B0[12]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[28],B0[12]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[12],B0[28]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[28],B0[28]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[13],B0[13]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[29],B0[13]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[13],B0[29]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[29],B0[29]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[14],B0[14]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[30],B0[14]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[14],B0[30]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[30],B0[30]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[15],B0[15]); d00.xm = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[31],B0[15]); d10.xm = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[15],B0[31]); d01.xm = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[31],B0[31]); d11.xm = _mm_add_ps(C11,c11);

                C0[ 0] += d00.xf[0] + d00.xf[1] + d00.xf[2] + d00.xf[3];
                C0[ 1] += d10.xf[0] + d10.xf[1] + d10.xf[2] + d10.xf[3];
                C1[ 0] += d01.xf[0] + d01.xf[1] + d01.xf[2] + d01.xf[3];
                C1[ 1] += d11.xf[0] + d11.xf[1] + d11.xf[2] + d11.xf[3];
                A0 += 32;
                C0 += 2;
                C1 += 2;
            }
            if (Mc2) {
                C00 = _mm_mul_ps(A0[ 0],B0[ 0]); 
                C01 = _mm_mul_ps(A0[ 0],B0[16]); 

                c00 = _mm_mul_ps(A0[ 1],B0[ 1]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[ 1],B0[17]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[ 2],B0[ 2]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[ 2],B0[18]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[ 3],B0[ 3]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[ 3],B0[19]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[ 4],B0[ 4]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[ 4],B0[20]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[ 5],B0[ 5]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[ 5],B0[21]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[ 6],B0[ 6]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[ 6],B0[22]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[ 7],B0[ 7]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[ 7],B0[23]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[ 8],B0[ 8]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[ 8],B0[24]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[ 9],B0[ 9]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[ 9],B0[25]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[10],B0[10]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[10],B0[26]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[11],B0[11]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[11],B0[27]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[12],B0[12]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[12],B0[28]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[13],B0[13]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[13],B0[29]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[14],B0[14]); C00 = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[14],B0[30]); C01 = _mm_add_ps(C01,c01);

                c00 = _mm_mul_ps(A0[15],B0[15]); d00.xm = _mm_add_ps(C00,c00);
                c01 = _mm_mul_ps(A0[15],B0[31]); d01.xm = _mm_add_ps(C01,c01);

                *C0++ += d00.xf[0] + d00.xf[1] + d00.xf[2] + d00.xf[3];
                *C1++ += d01.xf[0] + d01.xf[1] + d01.xf[2] + d01.xf[3];
                A0 += 16;
            }
            A0 -= Mx16;
            B0 += 32;
            C0 += M;
            C1 += M;
        } while (--j2);
    }
    static void sse_multmm_16_N_64(
        const ptrdiff_t N, const float* A, const float* B, float* C0)
    {
        TMVAssert( ((size_t)(A) & 0xf) == 0 );
        TMVAssert( ((size_t)(B) & 0xf) == 0 );
        TMVAssert( ((size_t)(C0) & 0xf) == 0 );
        const ptrdiff_t N_2 = N>>1;        // = N/2
        const ptrdiff_t Nc = (N-(N_2<<1)); // = N%2

        __m128 C00, C01, C10, C11;
        __m128 c00, c01, c10, c11;
        union { __m128 xm; float xf[4]; } d00, d01, d10, d11;

        const __m128* A0 = (const __m128*) A;
        const __m128* B0 = (const __m128*) B;

        ptrdiff_t i2,j2;

        j2 = N_2; if (j2) do {
            i2 = 4; do {
                C00 = _mm_mul_ps(A0[ 0],B0[ 0]); 
                C10 = _mm_mul_ps(A0[16],B0[ 0]); 
                C01 = _mm_mul_ps(A0[ 0],B0[16]); 
                C11 = _mm_mul_ps(A0[16],B0[16]); 

                c00 = _mm_mul_ps(A0[ 1],B0[ 1]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[17],B0[ 1]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 1],B0[17]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[17],B0[17]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 2],B0[ 2]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[18],B0[ 2]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 2],B0[18]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[18],B0[18]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 3],B0[ 3]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[19],B0[ 3]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 3],B0[19]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[19],B0[19]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 4],B0[ 4]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[20],B0[ 4]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 4],B0[20]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[20],B0[20]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 5],B0[ 5]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[21],B0[ 5]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 5],B0[21]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[21],B0[21]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 6],B0[ 6]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[22],B0[ 6]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 6],B0[22]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[22],B0[22]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 7],B0[ 7]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[23],B0[ 7]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 7],B0[23]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[23],B0[23]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 8],B0[ 8]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[24],B0[ 8]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 8],B0[24]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[24],B0[24]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[ 9],B0[ 9]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[25],B0[ 9]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[ 9],B0[25]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[25],B0[25]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[10],B0[10]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[26],B0[10]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[10],B0[26]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[26],B0[26]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[11],B0[11]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[27],B0[11]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[11],B0[27]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[27],B0[27]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[12],B0[12]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[28],B0[12]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[12],B0[28]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[28],B0[28]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[13],B0[13]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[29],B0[13]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[13],B0[29]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[29],B0[29]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[14],B0[14]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[30],B0[14]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[14],B0[30]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[30],B0[30]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[15],B0[15]); d00.xm = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[31],B0[15]); d10.xm = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[15],B0[31]); d01.xm = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[31],B0[31]); d11.xm = _mm_add_ps(C11,c11);

                C0[ 0] += d00.xf[0] + d00.xf[1] + d00.xf[2] + d00.xf[3];
                C00 = _mm_mul_ps(A0[47],B0[15]);
                C10 = _mm_mul_ps(A0[63],B0[15]);
                C01 = _mm_mul_ps(A0[47],B0[31]);
                C11 = _mm_mul_ps(A0[63],B0[31]);

                C0[ 1] += d10.xf[0] + d10.xf[1] + d10.xf[2] + d10.xf[3];
                c00 = _mm_mul_ps(A0[46],B0[14]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[62],B0[14]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[46],B0[30]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[62],B0[30]); C11 = _mm_add_ps(C11,c11);

                C0[16] += d01.xf[0] + d01.xf[1] + d01.xf[2] + d01.xf[3];
                c00 = _mm_mul_ps(A0[45],B0[13]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[61],B0[13]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[45],B0[29]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[61],B0[29]); C11 = _mm_add_ps(C11,c11);

                C0[17] += d11.xf[0] + d11.xf[1] + d11.xf[2] + d11.xf[3];
                c00 = _mm_mul_ps(A0[44],B0[12]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[60],B0[12]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[44],B0[28]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[60],B0[28]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[43],B0[11]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[59],B0[11]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[43],B0[27]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[59],B0[27]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[42],B0[10]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[58],B0[10]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[42],B0[26]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[58],B0[26]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[41],B0[ 9]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[57],B0[ 9]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[41],B0[25]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[57],B0[25]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[40],B0[ 8]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[56],B0[ 8]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[40],B0[24]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[56],B0[24]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[39],B0[ 7]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[55],B0[ 7]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[39],B0[23]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[55],B0[23]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[38],B0[ 6]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[54],B0[ 6]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[38],B0[22]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[54],B0[22]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[37],B0[ 5]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[53],B0[ 5]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[37],B0[21]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[53],B0[21]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[36],B0[ 4]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[52],B0[ 4]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[36],B0[20]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[52],B0[20]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[35],B0[ 3]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[51],B0[ 3]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[35],B0[19]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[51],B0[19]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[34],B0[ 2]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[50],B0[ 2]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[34],B0[18]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[50],B0[18]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[33],B0[ 1]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[49],B0[ 1]); C10 = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[33],B0[17]); C01 = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[49],B0[17]); C11 = _mm_add_ps(C11,c11);

                c00 = _mm_mul_ps(A0[32],B0[ 0]); d00.xm = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[48],B0[ 0]); d10.xm = _mm_add_ps(C10,c10);
                c01 = _mm_mul_ps(A0[32],B0[16]); d01.xm = _mm_add_ps(C01,c01);
                c11 = _mm_mul_ps(A0[48],B0[16]); d11.xm = _mm_add_ps(C11,c11);

                C0[ 2] += d00.xf[0] + d00.xf[1] + d00.xf[2] + d00.xf[3];
                C0[ 3] += d10.xf[0] + d10.xf[1] + d10.xf[2] + d10.xf[3];
                C0[18] += d01.xf[0] + d01.xf[1] + d01.xf[2] + d01.xf[3];
                C0[19] += d11.xf[0] + d11.xf[1] + d11.xf[2] + d11.xf[3];

                A0 += 64;
                C0 += 4;
            } while (--i2);
            A0 -= 256;
            B0 += 32;
            C0 += 16;
        } while (--j2);
        if (Nc) {
            i2 = 4; do {
                C00 = _mm_mul_ps(A0[ 0],B0[ 0]); 
                C10 = _mm_mul_ps(A0[16],B0[ 0]); 

                c00 = _mm_mul_ps(A0[ 1],B0[ 1]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[17],B0[ 1]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[ 2],B0[ 2]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[18],B0[ 2]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[ 3],B0[ 3]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[19],B0[ 3]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[ 4],B0[ 4]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[20],B0[ 4]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[ 5],B0[ 5]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[21],B0[ 5]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[ 6],B0[ 6]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[22],B0[ 6]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[ 7],B0[ 7]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[23],B0[ 7]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[ 8],B0[ 8]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[24],B0[ 8]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[ 9],B0[ 9]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[25],B0[ 9]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[10],B0[10]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[26],B0[10]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[11],B0[11]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[27],B0[11]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[12],B0[12]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[28],B0[12]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[13],B0[13]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[29],B0[13]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[14],B0[14]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[30],B0[14]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[15],B0[15]); d00.xm = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[31],B0[15]); d10.xm = _mm_add_ps(C10,c10);

                C0[ 0] += d00.xf[0] + d00.xf[1] + d00.xf[2] + d00.xf[3];
                C00 = _mm_mul_ps(A0[47],B0[15]);
                C10 = _mm_mul_ps(A0[63],B0[15]);

                C0[ 1] += d10.xf[0] + d10.xf[1] + d10.xf[2] + d10.xf[3];
                c00 = _mm_mul_ps(A0[46],B0[14]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[62],B0[14]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[45],B0[13]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[61],B0[13]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[44],B0[12]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[60],B0[12]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[43],B0[11]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[59],B0[11]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[42],B0[10]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[58],B0[10]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[41],B0[ 9]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[57],B0[ 9]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[40],B0[ 8]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[56],B0[ 8]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[39],B0[ 7]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[55],B0[ 7]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[38],B0[ 6]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[54],B0[ 6]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[37],B0[ 5]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[53],B0[ 5]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[36],B0[ 4]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[52],B0[ 4]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[35],B0[ 3]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[51],B0[ 3]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[34],B0[ 2]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[50],B0[ 2]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[33],B0[ 1]); C00 = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[49],B0[ 1]); C10 = _mm_add_ps(C10,c10);

                c00 = _mm_mul_ps(A0[32],B0[ 0]); d00.xm = _mm_add_ps(C00,c00);
                c10 = _mm_mul_ps(A0[48],B0[ 0]); d10.xm = _mm_add_ps(C10,c10);

                C0[ 2] += d00.xf[0] + d00.xf[1] + d00.xf[2] + d00.xf[3];
                C0[ 3] += d10.xf[0] + d10.xf[1] + d10.xf[2] + d10.xf[3];

                A0 += 64;
                C0 += 4;
            } while (--i2);
        }
    }
    static void sse_multmm_16_16_K(
        const ptrdiff_t K, const float* A, const float* B, float* C0)
    {
        TMVAssert( ((size_t)(A) & 0xf) == 0 );
        TMVAssert( ((size_t)(B) & 0xf) == 0 );
        TMVAssert( ((size_t)(C0) & 0xf) == 0 );
        TMVAssert(K > 0);
        const ptrdiff_t K_4 = K>>2; // K_4 = K/4 
        const ptrdiff_t K_4m2 = K_4-2; 
        const ptrdiff_t Kc = K-(K_4<<2);
        const ptrdiff_t Kd_4 = ((K-1)>>2)+1;  // = Kd/4
        //const ptrdiff_t Kd = Kd_4<<2; // = K rounded up to next multiple of 4
        const ptrdiff_t Kd_2 = Kd_4<<1;  // = Kd/2
        const ptrdiff_t Kd16 = Kd_4<<4; // = (Kd/4) * 16

        union ff { __m128 xm; float xf[4]; }; // ff = four floats
        __m128 C00, C01, C10, C11;
        __m128 c00, c01, c10, c11;
        ff d00, d01, d10, d11;
        float e00, e01, e10, e11;

        const ff* A0 = (const ff*) A;
        const ff* A1 = A0 + Kd_4;

        const ff* B0 = (const ff*) B;
        const ff* B1 = B0 + Kd_4;

        float* C1 = C0 + 16;

        ptrdiff_t j,i,k;

        if (K_4 >= 2) {
            // this is easier if we know that there are at least two m128's
            // so we do this first and do K_4 = 0,1 below.
            j = 8; do {
                i = 8; do {
                    C00 = _mm_mul_ps((*A0).xm,(*B0).xm); 
                    C10 = _mm_mul_ps((*A1).xm,(*B0++).xm); 
                    C01 = _mm_mul_ps((*A0++).xm,(*B1).xm); 
                    C11 = _mm_mul_ps((*A1++).xm,(*B1++).xm); 

                    k = K_4m2; if (k) do {
                        c00 = _mm_mul_ps((*A0).xm,(*B0).xm); 
                        C00 = _mm_add_ps(C00,c00);
                        c10 = _mm_mul_ps((*A1).xm,(*B0++).xm);
                        C10 = _mm_add_ps(C10,c10);
                        c01 = _mm_mul_ps((*A0++).xm,(*B1).xm);
                        C01 = _mm_add_ps(C01,c01);
                        c11 = _mm_mul_ps((*A1++).xm,(*B1++).xm);
                        C11 = _mm_add_ps(C11,c11);
                    } while (--k);

                    c00 = _mm_mul_ps((*A0).xm,(*B0).xm); 
                    d00.xm = _mm_add_ps(C00,c00);
                    c10 = _mm_mul_ps((*A1).xm,(*B0++).xm);
                    d10.xm = _mm_add_ps(C10,c10);
                    c01 = _mm_mul_ps((*A0++).xm,(*B1).xm);
                    d01.xm = _mm_add_ps(C01,c01);
                    c11 = _mm_mul_ps((*A1++).xm,(*B1++).xm);
                    d11.xm = _mm_add_ps(C11,c11);

                    e00 = d00.xf[0] + d00.xf[1] + d00.xf[2] + d00.xf[3];
                    e10 = d10.xf[0] + d10.xf[1] + d10.xf[2] + d10.xf[3];
                    e01 = d01.xf[0] + d01.xf[1] + d01.xf[2] + d01.xf[3];
                    e11 = d11.xf[0] + d11.xf[1] + d11.xf[2] + d11.xf[3];

                    if (Kc) {
                        d00.xm = _mm_mul_ps((*A0).xm,(*B0).xm); 
                        d10.xm = _mm_mul_ps((*A1).xm,(*B0++).xm);
                        d01.xm = _mm_mul_ps((*A0++).xm,(*B1).xm);
                        d11.xm = _mm_mul_ps((*A1++).xm,(*B1++).xm);

                        for(k=0;k<Kc;++k) {
                            e00 += d00.xf[k];
                            e10 += d10.xf[k];
                            e01 += d01.xf[k];
                            e11 += d11.xf[k];
                        }
                    }

                    C0[0] += e00;
                    C0[1] += e10; C0 += 2;
                    C1[0] += e01;
                    C1[1] += e11; C1 += 2;

                    A0 += Kd_4; A1 += Kd_4;
                    B0 -= Kd_4; B1 -= Kd_4;
                } while (--i);
                A0 -= Kd16; A1 -= Kd16;
                B0 += Kd_2; B1 += Kd_2;
                C0 += 16; C1 += 16;
            } while (--j);
        } else if (K_4 == 1) { // K_4 == 1
            j = 8; do {
                i = 8; do {
                    d00.xm = _mm_mul_ps((*A0).xm,(*B0).xm); 
                    d10.xm = _mm_mul_ps((*A1).xm,(*B0++).xm); 
                    d01.xm = _mm_mul_ps((*A0++).xm,(*B1).xm); 
                    d11.xm = _mm_mul_ps((*A1++).xm,(*B1++).xm); 

                    e00 = d00.xf[0] + d00.xf[1] + d00.xf[2] + d00.xf[3];
                    e10 = d10.xf[0] + d10.xf[1] + d10.xf[2] + d10.xf[3];
                    e01 = d01.xf[0] + d01.xf[1] + d01.xf[2] + d01.xf[3];
                    e11 = d11.xf[0] + d11.xf[1] + d11.xf[2] + d11.xf[3];

                    if (Kc) {
                        d00.xm = _mm_mul_ps((*A0).xm,(*B0).xm); 
                        d10.xm = _mm_mul_ps((*A1).xm,(*B0++).xm);
                        d01.xm = _mm_mul_ps((*A0++).xm,(*B1).xm);
                        d11.xm = _mm_mul_ps((*A1++).xm,(*B1++).xm);

                        for(k=0;k<Kc;++k) {
                            e00 += d00.xf[k];
                            e10 += d10.xf[k];
                            e01 += d01.xf[k];
                            e11 += d11.xf[k];
                        }
                    }

                    C0[0] += e00;
                    C0[1] += e10; C0 += 2;
                    C1[0] += e01;
                    C1[1] += e11; C1 += 2;
                    A0 += Kd_4; A1 += Kd_4;
                    B0 -= Kd_4; B1 -= Kd_4;
                } while (--i);
                A0 -= Kd16; A1 -= Kd16;
                B0 += Kd_2; B1 += Kd_2;
                C0 += 16; C1 += 16;
            } while (--j);
        } else { // K_4 == 0 and Kc > 0
            j = 8; do {
                i = 8; do {
                    d00.xm = _mm_mul_ps((*A0).xm,(*B0).xm); 
                    d10.xm = _mm_mul_ps((*A1).xm,(*B0++).xm);
                    d01.xm = _mm_mul_ps((*A0++).xm,(*B1).xm);
                    d11.xm = _mm_mul_ps((*A1++).xm,(*B1++).xm);

                    e00 = d00.xf[0];
                    e10 = d10.xf[0];
                    e01 = d01.xf[0];
                    e11 = d11.xf[0];
                    for(k=1;k<Kc;++k) {
                        e00 += d00.xf[k];
                        e10 += d10.xf[k];
                        e01 += d01.xf[k];
                        e11 += d11.xf[k];
                    }
                    C0[0] += e00;
                    C0[1] += e10; C0 += 2;
                    C1[0] += e01;
                    C1[1] += e11; C1 += 2;
                    A0 += Kd_4; A1 += Kd_4;
                    B0 -= Kd_4; B1 -= Kd_4;
                } while (--i);
                A0 -= Kd16; A1 -= Kd16;
                B0 += Kd_2; B1 += Kd_2;
                C0 += 16; C1 += 16;
            } while (--j);
        }
    }
#endif

#ifdef __SSE__
    static void multmm_16_16_K(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t K,
        const float* A0, const float* B0, float* C0)
    { sse_multmm_16_16_K(K,A0,B0,C0); }

    static void multmm_M_N_K(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const float* A0, const float* B0, float* C0)
    {
        const ptrdiff_t Kd = (((K-1)>>2)+1)<<2;
        generic_multmm_M_N_K(M,N,K,Kd,A0,B0,C0); 
    }

    static void multmm_16_16_64(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t ,
        const float* A0, const float* B0, float* C0)
    { sse_multmm_16_16_64(A0,B0,C0); }

    static void multmm_M_16_64(
        const ptrdiff_t M, const ptrdiff_t , const ptrdiff_t ,
        const float* A0, const float* B0, float* C0)
    { sse_multmm_M_16_64(M,A0,B0,C0); }

    static void multmm_16_N_64(
        const ptrdiff_t , const ptrdiff_t N, const ptrdiff_t ,
        const float* A0, const float* B0, float* C0)
    { sse_multmm_16_N_64(N,A0,B0,C0); }

#else
    // If no SSE, then repeat the 16 block version to call generic.
    static void multmm_16_16_16(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t ,
        const float* A0, const float* B0, float* C0)
    { generic_multmm_16_16_16(A0,B0,C0); }

    static void multmm_M_16_16(
        const ptrdiff_t M, const ptrdiff_t , const ptrdiff_t ,
        const float* A0, const float* B0, float* C0)
    { generic_multmm_M_16_16(M,A0,B0,C0); }

    static void multmm_16_N_16(
        const ptrdiff_t , const ptrdiff_t N, const ptrdiff_t ,
        const float* A0, const float* B0, float* C0)
    { generic_multmm_16_N_16(N,A0,B0,C0); }

    static void multmm_16_16_K(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t K,
        const float* A0, const float* B0, float* C0)
    { generic_multmm_16_16_K(K,A0,B0,C0); }

    static void multmm_M_N_K(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const float* A0, const float* B0, float* C0)
    { generic_multmm_M_N_K(M,N,K,K,A0,B0,C0); }

#endif

#endif

#ifdef INST_DOUBLE
    //
    // double
    //

#ifdef __SSE2__
    // If we have SSE2 commands available, we have special versions for double
    static void sse2_multmm_16_16_32(
        const double* A, const double* B, double* C0)
    {
        TMVAssert( ((size_t)(A) & 0xf) == 0 );
        TMVAssert( ((size_t)(B) & 0xf) == 0 );
        TMVAssert( ((size_t)(C0) & 0xf) == 0 );

        __m128d C00, C01, C10, C11;
        __m128d c00, c01, c10, c11;
        union { __m128d xm; double xd[2]; } d00, d01, d10, d11;

        const __m128d* A0 = (const __m128d*) A;
        const __m128d* B0 = (const __m128d*) B;

        ptrdiff_t i2,j2;

        j2 = 8; do {
            i2 = 4; do {
                C00 = _mm_mul_pd(A0[ 0],B0[ 0]); 
                C10 = _mm_mul_pd(A0[16],B0[ 0]); 
                C01 = _mm_mul_pd(A0[ 0],B0[16]); 
                C11 = _mm_mul_pd(A0[16],B0[16]); 

                c00 = _mm_mul_pd(A0[ 1],B0[ 1]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[17],B0[ 1]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 1],B0[17]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[17],B0[17]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 2],B0[ 2]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[18],B0[ 2]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 2],B0[18]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[18],B0[18]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 3],B0[ 3]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[19],B0[ 3]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 3],B0[19]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[19],B0[19]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 4],B0[ 4]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[20],B0[ 4]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 4],B0[20]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[20],B0[20]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 5],B0[ 5]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[21],B0[ 5]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 5],B0[21]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[21],B0[21]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 6],B0[ 6]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[22],B0[ 6]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 6],B0[22]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[22],B0[22]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 7],B0[ 7]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[23],B0[ 7]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 7],B0[23]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[23],B0[23]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 8],B0[ 8]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[24],B0[ 8]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 8],B0[24]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[24],B0[24]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 9],B0[ 9]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[25],B0[ 9]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 9],B0[25]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[25],B0[25]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[10],B0[10]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[26],B0[10]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[10],B0[26]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[26],B0[26]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[11],B0[11]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[27],B0[11]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[11],B0[27]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[27],B0[27]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[12],B0[12]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[28],B0[12]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[12],B0[28]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[28],B0[28]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[13],B0[13]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[29],B0[13]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[13],B0[29]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[29],B0[29]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[14],B0[14]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[30],B0[14]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[14],B0[30]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[30],B0[30]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[15],B0[15]); d00.xm = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[31],B0[15]); d10.xm = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[15],B0[31]); d01.xm = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[31],B0[31]); d11.xm = _mm_add_pd(C11,c11);

                C0[ 0] += d00.xd[0] + d00.xd[1];
                C00 = _mm_mul_pd(A0[47],B0[15]);
                C10 = _mm_mul_pd(A0[63],B0[15]);
                C01 = _mm_mul_pd(A0[47],B0[31]);
                C11 = _mm_mul_pd(A0[63],B0[31]);

                C0[ 1] += d10.xd[0] + d10.xd[1];
                c00 = _mm_mul_pd(A0[46],B0[14]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[62],B0[14]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[46],B0[30]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[62],B0[30]); C11 = _mm_add_pd(C11,c11);

                C0[16] += d01.xd[0] + d01.xd[1];
                c00 = _mm_mul_pd(A0[45],B0[13]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[61],B0[13]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[45],B0[29]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[61],B0[29]); C11 = _mm_add_pd(C11,c11);

                C0[17] += d11.xd[0] + d11.xd[1];
                c00 = _mm_mul_pd(A0[44],B0[12]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[60],B0[12]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[44],B0[28]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[60],B0[28]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[43],B0[11]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[59],B0[11]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[43],B0[27]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[59],B0[27]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[42],B0[10]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[58],B0[10]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[42],B0[26]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[58],B0[26]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[41],B0[ 9]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[57],B0[ 9]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[41],B0[25]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[57],B0[25]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[40],B0[ 8]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[56],B0[ 8]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[40],B0[24]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[56],B0[24]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[39],B0[ 7]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[55],B0[ 7]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[39],B0[23]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[55],B0[23]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[38],B0[ 6]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[54],B0[ 6]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[38],B0[22]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[54],B0[22]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[37],B0[ 5]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[53],B0[ 5]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[37],B0[21]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[53],B0[21]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[36],B0[ 4]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[52],B0[ 4]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[36],B0[20]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[52],B0[20]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[35],B0[ 3]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[51],B0[ 3]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[35],B0[19]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[51],B0[19]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[34],B0[ 2]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[50],B0[ 2]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[34],B0[18]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[50],B0[18]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[33],B0[ 1]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[49],B0[ 1]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[33],B0[17]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[49],B0[17]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[32],B0[ 0]); d00.xm = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[48],B0[ 0]); d10.xm = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[32],B0[16]); d01.xm = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[48],B0[16]); d11.xm = _mm_add_pd(C11,c11);

                C0[ 2] += d00.xd[0] + d00.xd[1];
                C0[ 3] += d10.xd[0] + d10.xd[1];
                C0[18] += d01.xd[0] + d01.xd[1];
                C0[19] += d11.xd[0] + d11.xd[1];

                A0 += 64;
                C0 += 4;
            } while (--i2);
            A0 -= 256;
            B0 += 32;
            C0 += 16; 
        } while (--j2);
    }
    static void sse2_multmm_M_16_32(
        const ptrdiff_t M, const double* A, const double* B, double* C0)
    {
        TMVAssert( ((size_t)(A) & 0xf) == 0 );
        TMVAssert( ((size_t)(B) & 0xf) == 0 );
        TMVAssert( ((size_t)(C0) & 0xf) == 0 );
        const ptrdiff_t Mx16 = (M<<4);      // = M*16
        const ptrdiff_t M_4 = M>>2;        // = M/4
        const ptrdiff_t Mc = (M-(M_4<<2)); // = M%4
        const ptrdiff_t Mc1 = Mc>>1;       // = M%4 == 2 or 3
        const ptrdiff_t Mc2 = Mc-(Mc1<<1); // = M%4 == 1 or 3

        __m128d C00, C01, C10, C11;
        __m128d c00, c01, c10, c11;
        union { __m128d xm; double xd[2]; } d00, d01, d10, d11;

        const __m128d* A0 = (const __m128d*) A;
        const __m128d* B0 = (const __m128d*) B;
        double* C1 = C0 + M;

        ptrdiff_t i2,j2;

        j2 = 8; do {
            i2 = M_4; if (i2) do {
                C00 = _mm_mul_pd(A0[ 0],B0[ 0]); 
                C10 = _mm_mul_pd(A0[16],B0[ 0]); 
                C01 = _mm_mul_pd(A0[ 0],B0[16]); 
                C11 = _mm_mul_pd(A0[16],B0[16]); 

                c00 = _mm_mul_pd(A0[ 1],B0[ 1]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[17],B0[ 1]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 1],B0[17]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[17],B0[17]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 2],B0[ 2]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[18],B0[ 2]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 2],B0[18]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[18],B0[18]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 3],B0[ 3]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[19],B0[ 3]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 3],B0[19]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[19],B0[19]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 4],B0[ 4]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[20],B0[ 4]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 4],B0[20]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[20],B0[20]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 5],B0[ 5]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[21],B0[ 5]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 5],B0[21]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[21],B0[21]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 6],B0[ 6]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[22],B0[ 6]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 6],B0[22]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[22],B0[22]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 7],B0[ 7]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[23],B0[ 7]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 7],B0[23]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[23],B0[23]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 8],B0[ 8]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[24],B0[ 8]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 8],B0[24]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[24],B0[24]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 9],B0[ 9]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[25],B0[ 9]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 9],B0[25]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[25],B0[25]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[10],B0[10]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[26],B0[10]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[10],B0[26]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[26],B0[26]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[11],B0[11]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[27],B0[11]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[11],B0[27]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[27],B0[27]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[12],B0[12]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[28],B0[12]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[12],B0[28]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[28],B0[28]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[13],B0[13]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[29],B0[13]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[13],B0[29]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[29],B0[29]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[14],B0[14]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[30],B0[14]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[14],B0[30]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[30],B0[30]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[15],B0[15]); d00.xm = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[31],B0[15]); d10.xm = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[15],B0[31]); d01.xm = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[31],B0[31]); d11.xm = _mm_add_pd(C11,c11);

                C0[ 0] += d00.xd[0] + d00.xd[1];
                C00 = _mm_mul_pd(A0[47],B0[15]);
                C10 = _mm_mul_pd(A0[63],B0[15]);
                C01 = _mm_mul_pd(A0[47],B0[31]);
                C11 = _mm_mul_pd(A0[63],B0[31]);

                C0[ 1] += d10.xd[0] + d10.xd[1];
                c00 = _mm_mul_pd(A0[46],B0[14]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[62],B0[14]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[46],B0[30]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[62],B0[30]); C11 = _mm_add_pd(C11,c11);

                C1[ 0] += d01.xd[0] + d01.xd[1];
                c00 = _mm_mul_pd(A0[45],B0[13]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[61],B0[13]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[45],B0[29]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[61],B0[29]); C11 = _mm_add_pd(C11,c11);

                C1[ 1] += d11.xd[0] + d11.xd[1];
                c00 = _mm_mul_pd(A0[44],B0[12]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[60],B0[12]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[44],B0[28]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[60],B0[28]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[43],B0[11]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[59],B0[11]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[43],B0[27]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[59],B0[27]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[42],B0[10]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[58],B0[10]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[42],B0[26]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[58],B0[26]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[41],B0[ 9]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[57],B0[ 9]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[41],B0[25]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[57],B0[25]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[40],B0[ 8]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[56],B0[ 8]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[40],B0[24]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[56],B0[24]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[39],B0[ 7]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[55],B0[ 7]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[39],B0[23]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[55],B0[23]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[38],B0[ 6]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[54],B0[ 6]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[38],B0[22]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[54],B0[22]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[37],B0[ 5]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[53],B0[ 5]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[37],B0[21]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[53],B0[21]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[36],B0[ 4]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[52],B0[ 4]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[36],B0[20]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[52],B0[20]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[35],B0[ 3]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[51],B0[ 3]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[35],B0[19]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[51],B0[19]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[34],B0[ 2]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[50],B0[ 2]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[34],B0[18]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[50],B0[18]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[33],B0[ 1]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[49],B0[ 1]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[33],B0[17]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[49],B0[17]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[32],B0[ 0]); d00.xm = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[48],B0[ 0]); d10.xm = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[32],B0[16]); d01.xm = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[48],B0[16]); d11.xm = _mm_add_pd(C11,c11);

                C0[ 2] += d00.xd[0] + d00.xd[1];
                C0[ 3] += d10.xd[0] + d10.xd[1];
                C1[ 2] += d01.xd[0] + d01.xd[1];
                C1[ 3] += d11.xd[0] + d11.xd[1];

                A0 += 64;
                C0 += 4;
                C1 += 4;
            } while (--i2);
            if (Mc1) {
                C00 = _mm_mul_pd(A0[ 0],B0[ 0]); 
                C10 = _mm_mul_pd(A0[16],B0[ 0]); 
                C01 = _mm_mul_pd(A0[ 0],B0[16]); 
                C11 = _mm_mul_pd(A0[16],B0[16]); 

                c00 = _mm_mul_pd(A0[ 1],B0[ 1]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[17],B0[ 1]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 1],B0[17]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[17],B0[17]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 2],B0[ 2]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[18],B0[ 2]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 2],B0[18]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[18],B0[18]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 3],B0[ 3]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[19],B0[ 3]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 3],B0[19]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[19],B0[19]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 4],B0[ 4]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[20],B0[ 4]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 4],B0[20]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[20],B0[20]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 5],B0[ 5]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[21],B0[ 5]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 5],B0[21]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[21],B0[21]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 6],B0[ 6]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[22],B0[ 6]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 6],B0[22]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[22],B0[22]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 7],B0[ 7]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[23],B0[ 7]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 7],B0[23]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[23],B0[23]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 8],B0[ 8]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[24],B0[ 8]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 8],B0[24]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[24],B0[24]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 9],B0[ 9]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[25],B0[ 9]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 9],B0[25]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[25],B0[25]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[10],B0[10]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[26],B0[10]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[10],B0[26]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[26],B0[26]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[11],B0[11]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[27],B0[11]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[11],B0[27]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[27],B0[27]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[12],B0[12]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[28],B0[12]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[12],B0[28]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[28],B0[28]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[13],B0[13]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[29],B0[13]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[13],B0[29]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[29],B0[29]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[14],B0[14]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[30],B0[14]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[14],B0[30]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[30],B0[30]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[15],B0[15]); d00.xm = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[31],B0[15]); d10.xm = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[15],B0[31]); d01.xm = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[31],B0[31]); d11.xm = _mm_add_pd(C11,c11);

                C0[ 0] += d00.xd[0] + d00.xd[1];
                C0[ 1] += d10.xd[0] + d10.xd[1];
                C1[ 0] += d01.xd[0] + d01.xd[1];
                C1[ 1] += d11.xd[0] + d11.xd[1];

                A0 += 32;
                C0 += 2;
                C1 += 2;
            }
            if (Mc2) {
                C00 = _mm_mul_pd(A0[ 0],B0[ 0]); 
                C01 = _mm_mul_pd(A0[ 0],B0[16]); 

                c00 = _mm_mul_pd(A0[ 1],B0[ 1]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[ 1],B0[17]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[ 2],B0[ 2]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[ 2],B0[18]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[ 3],B0[ 3]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[ 3],B0[19]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[ 4],B0[ 4]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[ 4],B0[20]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[ 5],B0[ 5]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[ 5],B0[21]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[ 6],B0[ 6]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[ 6],B0[22]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[ 7],B0[ 7]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[ 7],B0[23]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[ 8],B0[ 8]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[ 8],B0[24]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[ 9],B0[ 9]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[ 9],B0[25]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[10],B0[10]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[10],B0[26]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[11],B0[11]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[11],B0[27]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[12],B0[12]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[12],B0[28]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[13],B0[13]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[13],B0[29]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[14],B0[14]); C00 = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[14],B0[30]); C01 = _mm_add_pd(C01,c01);

                c00 = _mm_mul_pd(A0[15],B0[15]); d00.xm = _mm_add_pd(C00,c00);
                c01 = _mm_mul_pd(A0[15],B0[31]); d01.xm = _mm_add_pd(C01,c01);

                *C0++ += d00.xd[0] + d00.xd[1];
                *C1++ += d01.xd[0] + d01.xd[1];

                A0 += 16;
            }
            A0 -= Mx16;
            B0 += 32;
            C0 += M; 
            C1 += M; 
        } while (--j2);
    }
    static void sse2_multmm_16_N_32(
        const ptrdiff_t N, const double* A, const double* B, double* C0)
    {
        TMVAssert( ((size_t)(A) & 0xf) == 0 );
        TMVAssert( ((size_t)(B) & 0xf) == 0 );
        TMVAssert( ((size_t)(C0) & 0xf) == 0 );
        const ptrdiff_t N_2 = N>>1;        // = N/2
        const ptrdiff_t Nc = (N-(N_2<<1)); // = N%2

        __m128d C00, C01, C10, C11;
        __m128d c00, c01, c10, c11;
        union { __m128d xm; double xd[2]; } d00, d01, d10, d11;

        const __m128d* A0 = (const __m128d*) A;
        const __m128d* B0 = (const __m128d*) B;

        ptrdiff_t i2,j2;

        j2 = N_2; if (j2) do {
            i2 = 4; do {
                C00 = _mm_mul_pd(A0[ 0],B0[ 0]); 
                C10 = _mm_mul_pd(A0[16],B0[ 0]); 
                C01 = _mm_mul_pd(A0[ 0],B0[16]); 
                C11 = _mm_mul_pd(A0[16],B0[16]); 

                c00 = _mm_mul_pd(A0[ 1],B0[ 1]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[17],B0[ 1]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 1],B0[17]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[17],B0[17]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 2],B0[ 2]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[18],B0[ 2]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 2],B0[18]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[18],B0[18]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 3],B0[ 3]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[19],B0[ 3]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 3],B0[19]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[19],B0[19]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 4],B0[ 4]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[20],B0[ 4]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 4],B0[20]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[20],B0[20]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 5],B0[ 5]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[21],B0[ 5]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 5],B0[21]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[21],B0[21]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 6],B0[ 6]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[22],B0[ 6]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 6],B0[22]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[22],B0[22]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 7],B0[ 7]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[23],B0[ 7]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 7],B0[23]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[23],B0[23]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 8],B0[ 8]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[24],B0[ 8]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 8],B0[24]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[24],B0[24]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[ 9],B0[ 9]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[25],B0[ 9]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[ 9],B0[25]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[25],B0[25]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[10],B0[10]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[26],B0[10]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[10],B0[26]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[26],B0[26]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[11],B0[11]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[27],B0[11]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[11],B0[27]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[27],B0[27]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[12],B0[12]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[28],B0[12]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[12],B0[28]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[28],B0[28]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[13],B0[13]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[29],B0[13]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[13],B0[29]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[29],B0[29]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[14],B0[14]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[30],B0[14]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[14],B0[30]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[30],B0[30]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[15],B0[15]); d00.xm = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[31],B0[15]); d10.xm = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[15],B0[31]); d01.xm = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[31],B0[31]); d11.xm = _mm_add_pd(C11,c11);

                C0[ 0] += d00.xd[0] + d00.xd[1];
                C00 = _mm_mul_pd(A0[47],B0[15]);
                C10 = _mm_mul_pd(A0[63],B0[15]);
                C01 = _mm_mul_pd(A0[47],B0[31]);
                C11 = _mm_mul_pd(A0[63],B0[31]);

                C0[ 1] += d10.xd[0] + d10.xd[1];
                c00 = _mm_mul_pd(A0[46],B0[14]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[62],B0[14]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[46],B0[30]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[62],B0[30]); C11 = _mm_add_pd(C11,c11);

                C0[16] += d01.xd[0] + d01.xd[1];
                c00 = _mm_mul_pd(A0[45],B0[13]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[61],B0[13]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[45],B0[29]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[61],B0[29]); C11 = _mm_add_pd(C11,c11);

                C0[17] += d11.xd[0] + d11.xd[1];
                c00 = _mm_mul_pd(A0[44],B0[12]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[60],B0[12]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[44],B0[28]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[60],B0[28]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[43],B0[11]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[59],B0[11]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[43],B0[27]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[59],B0[27]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[42],B0[10]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[58],B0[10]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[42],B0[26]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[58],B0[26]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[41],B0[ 9]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[57],B0[ 9]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[41],B0[25]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[57],B0[25]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[40],B0[ 8]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[56],B0[ 8]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[40],B0[24]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[56],B0[24]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[39],B0[ 7]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[55],B0[ 7]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[39],B0[23]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[55],B0[23]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[38],B0[ 6]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[54],B0[ 6]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[38],B0[22]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[54],B0[22]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[37],B0[ 5]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[53],B0[ 5]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[37],B0[21]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[53],B0[21]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[36],B0[ 4]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[52],B0[ 4]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[36],B0[20]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[52],B0[20]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[35],B0[ 3]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[51],B0[ 3]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[35],B0[19]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[51],B0[19]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[34],B0[ 2]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[50],B0[ 2]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[34],B0[18]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[50],B0[18]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[33],B0[ 1]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[49],B0[ 1]); C10 = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[33],B0[17]); C01 = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[49],B0[17]); C11 = _mm_add_pd(C11,c11);

                c00 = _mm_mul_pd(A0[32],B0[ 0]); d00.xm = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[48],B0[ 0]); d10.xm = _mm_add_pd(C10,c10);
                c01 = _mm_mul_pd(A0[32],B0[16]); d01.xm = _mm_add_pd(C01,c01);
                c11 = _mm_mul_pd(A0[48],B0[16]); d11.xm = _mm_add_pd(C11,c11);

                C0[ 2] += d00.xd[0] + d00.xd[1];
                C0[ 3] += d10.xd[0] + d10.xd[1];
                C0[18] += d01.xd[0] + d01.xd[1];
                C0[19] += d11.xd[0] + d11.xd[1];

                A0 += 64;
                C0 += 4;
            } while (--i2);
            A0 -= 256;
            B0 += 32;
            C0 += 16; 
        } while (--j2);
        if (Nc) {
            i2 = 4; do {
                C00 = _mm_mul_pd(A0[ 0],B0[ 0]); 
                C10 = _mm_mul_pd(A0[16],B0[ 0]); 

                c00 = _mm_mul_pd(A0[ 1],B0[ 1]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[17],B0[ 1]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[ 2],B0[ 2]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[18],B0[ 2]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[ 3],B0[ 3]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[19],B0[ 3]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[ 4],B0[ 4]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[20],B0[ 4]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[ 5],B0[ 5]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[21],B0[ 5]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[ 6],B0[ 6]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[22],B0[ 6]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[ 7],B0[ 7]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[23],B0[ 7]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[ 8],B0[ 8]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[24],B0[ 8]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[ 9],B0[ 9]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[25],B0[ 9]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[10],B0[10]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[26],B0[10]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[11],B0[11]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[27],B0[11]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[12],B0[12]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[28],B0[12]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[13],B0[13]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[29],B0[13]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[14],B0[14]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[30],B0[14]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[15],B0[15]); d00.xm = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[31],B0[15]); d10.xm = _mm_add_pd(C10,c10);

                C0[ 0] += d00.xd[0] + d00.xd[1];
                C00 = _mm_mul_pd(A0[47],B0[15]);
                C10 = _mm_mul_pd(A0[63],B0[15]);

                C0[ 1] += d10.xd[0] + d10.xd[1];
                c00 = _mm_mul_pd(A0[46],B0[14]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[62],B0[14]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[45],B0[13]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[61],B0[13]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[44],B0[12]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[60],B0[12]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[43],B0[11]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[59],B0[11]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[42],B0[10]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[58],B0[10]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[41],B0[ 9]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[57],B0[ 9]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[40],B0[ 8]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[56],B0[ 8]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[39],B0[ 7]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[55],B0[ 7]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[38],B0[ 6]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[54],B0[ 6]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[37],B0[ 5]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[53],B0[ 5]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[36],B0[ 4]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[52],B0[ 4]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[35],B0[ 3]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[51],B0[ 3]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[34],B0[ 2]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[50],B0[ 2]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[33],B0[ 1]); C00 = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[49],B0[ 1]); C10 = _mm_add_pd(C10,c10);

                c00 = _mm_mul_pd(A0[32],B0[ 0]); d00.xm = _mm_add_pd(C00,c00);
                c10 = _mm_mul_pd(A0[48],B0[ 0]); d10.xm = _mm_add_pd(C10,c10);

                C0[ 2] += d00.xd[0] + d00.xd[1];
                C0[ 3] += d10.xd[0] + d10.xd[1];

                A0 += 64;
                C0 += 4;
            } while (--i2);
        }
    }
    static void sse2_multmm_16_16_K(
        const ptrdiff_t K, const double* A, const double* B, double* C0)
    {
        TMVAssert( ((size_t)(A) & 0xf) == 0 );
        TMVAssert( ((size_t)(B) & 0xf) == 0 );
        TMVAssert( ((size_t)(C0) & 0xf) == 0 );
        TMVAssert(K > 0);
        const ptrdiff_t K_2 = K>>1; // K_2 = K/2 
        const ptrdiff_t K_2m2 = K_2-2; 
        const ptrdiff_t Kc = K-(K_2<<1);
        const ptrdiff_t Kd = K + Kc; // = K rounded up to next multiple of 2
        const ptrdiff_t Kd_2 = Kd>>1; // = Kd/2
        const ptrdiff_t Kd16 = Kd_2<<4; // = (Kd/2)*16

        union dd { __m128d xm; double xd[2]; }; // dd = two doubles
        __m128d C00, C01, C10, C11;
        __m128d c00, c01, c10, c11;
        dd d00, d01, d10, d11;
        double e00, e01, e10, e11;

        const dd* A0 = (const dd*) A;
        const dd* A1 = A0 + Kd_2;

        const dd* B0 = (const dd*) B;
        const dd* B1 = B0 + Kd_2;

        double* C1 = C0 + 16;

        ptrdiff_t j,i,k;

        if (K_2 >= 2) {
            // this is easier if we know that there are at least two m128's
            // so we do this first and do K_2 = 1 below.
            j = 8; do {
                i = 8; do {
                    C00 = _mm_mul_pd((*A0).xm,(*B0).xm); 
                    C10 = _mm_mul_pd((*A1).xm,(*B0++).xm); 
                    C01 = _mm_mul_pd((*A0++).xm,(*B1).xm); 
                    C11 = _mm_mul_pd((*A1++).xm,(*B1++).xm); 

                    k = K_2m2; if (k) do {
                        c00 = _mm_mul_pd((*A0).xm,(*B0).xm); 
                        C00 = _mm_add_pd(C00,c00);
                        c10 = _mm_mul_pd((*A1).xm,(*B0++).xm);
                        C10 = _mm_add_pd(C10,c10);
                        c01 = _mm_mul_pd((*A0++).xm,(*B1).xm);
                        C01 = _mm_add_pd(C01,c01);
                        c11 = _mm_mul_pd((*A1++).xm,(*B1++).xm);
                        C11 = _mm_add_pd(C11,c11);
                    } while (--k);
                    c00 = _mm_mul_pd((*A0).xm,(*B0).xm); 
                    d00.xm = _mm_add_pd(C00,c00);
                    c10 = _mm_mul_pd((*A1).xm,(*B0++).xm);
                    d10.xm = _mm_add_pd(C10,c10);
                    c01 = _mm_mul_pd((*A0++).xm,(*B1).xm);
                    d01.xm = _mm_add_pd(C01,c01);
                    c11 = _mm_mul_pd((*A1++).xm,(*B1++).xm);
                    d11.xm = _mm_add_pd(C11,c11);

                    e00 = d00.xd[0] + d00.xd[1];
                    e10 = d10.xd[0] + d10.xd[1];
                    e01 = d01.xd[0] + d01.xd[1];
                    e11 = d11.xd[0] + d11.xd[1];
                    if (Kc) {
                        d00.xm = _mm_mul_pd((*A0).xm,(*B0).xm); 
                        d10.xm = _mm_mul_pd((*A1).xm,(*B0++).xm);
                        d01.xm = _mm_mul_pd((*A0++).xm,(*B1).xm);
                        d11.xm = _mm_mul_pd((*A1++).xm,(*B1++).xm);
                        e00 += d00.xd[0];
                        e10 += d10.xd[0];
                        e01 += d01.xd[0];
                        e11 += d11.xd[0];
                    }
                    C0[0] += e00;
                    C0[1] += e10; C0 += 2;
                    C1[0] += e01;
                    C1[1] += e11; C1 += 2;
                    A0 += Kd_2; A1 += Kd_2;
                    B0 -= Kd_2; B1 -= Kd_2;
                } while (--i);
                A0 -= Kd16; A1 -= Kd16;
                B0 += Kd; B1 += Kd;
                C0 += 16; C1 += 16;
            } while (--j);
        } else if (K_2 == 1) {
            j = 8; do {
                i = 8; do {
                    d00.xm = _mm_mul_pd((*A0).xm,(*B0).xm); 
                    d10.xm = _mm_mul_pd((*A1).xm,(*B0++).xm); 
                    d01.xm = _mm_mul_pd((*A0++).xm,(*B1).xm); 
                    d11.xm = _mm_mul_pd((*A1++).xm,(*B1++).xm); 

                    e00 = d00.xd[0] + d00.xd[1];
                    e10 = d10.xd[0] + d10.xd[1];
                    e01 = d01.xd[0] + d01.xd[1];
                    e11 = d11.xd[0] + d11.xd[1];

                    if (Kc) {
                        d00.xm = _mm_mul_pd((*A0).xm,(*B0).xm); 
                        d10.xm = _mm_mul_pd((*A1).xm,(*B0++).xm);
                        d01.xm = _mm_mul_pd((*A0++).xm,(*B1).xm);
                        d11.xm = _mm_mul_pd((*A1++).xm,(*B1++).xm);

                        e00 += d00.xd[0];
                        e10 += d10.xd[0];
                        e01 += d01.xd[0];
                        e11 += d11.xd[0];
                    }

                    C0[0] += e00;
                    C0[1] += e10; C0 += 2;
                    C1[0] += e01;
                    C1[1] += e11; C1 += 2;
                    A0 += Kd_2; A1 += Kd_2;
                    B0 -= Kd_2; B1 -= Kd_2;
                } while (--i);
                A0 -= Kd16; A1 -= Kd16;
                B0 += Kd; B1 += Kd;
                C0 += 16; C1 += 16;
            } while (--j);
        } else { // K_2 == 0 and Kc == 1
            j = 8; do {
                i = 8; do {
                    e00 = (*A0).xd[0] * (*B0).xd[0];
                    e10 = (*A1).xd[0] * (*B0++).xd[0];
                    e01 = (*A0++).xd[0] * (*B1).xd[0];
                    e11 = (*A1++).xd[0] * (*B1++).xd[0];

                    C0[0] += e00;
                    C0[1] += e10; C0 += 2;
                    C1[0] += e01;
                    C1[1] += e11; C1 += 2;
                    A0 += Kd_2; A1 += Kd_2;
                    B0 -= Kd_2; B1 -= Kd_2;
                } while (--i);
                A0 -= Kd16; A1 -= Kd16;
                B0 += Kd; B1 += Kd;
                C0 += 16; C1 += 16;
            } while (--j);
        }
    }
#endif

#ifdef __SSE2__
    static void multmm_16_16_32(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t ,
        const double* A0, const double* B0, double* C0)
    { sse2_multmm_16_16_32(A0,B0,C0); }

    static void multmm_M_16_32(
        const ptrdiff_t M, const ptrdiff_t , const ptrdiff_t ,
        const double* A0, const double* B0, double* C0)
    { sse2_multmm_M_16_32(M,A0,B0,C0); }

    static void multmm_16_N_32(
        const ptrdiff_t , const ptrdiff_t N, const ptrdiff_t ,
        const double* A0, const double* B0, double* C0)
    { sse2_multmm_16_N_32(N,A0,B0,C0); }

    static void multmm_16_16_K(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t K,
        const double* A0, const double* B0, double* C0)
    { sse2_multmm_16_16_K(K,A0,B0,C0); }

    static void multmm_M_N_K(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const double* A0, const double* B0, double* C0)
    {
        const ptrdiff_t Kd = (((K-1)>>1)+1)<<1;
        generic_multmm_M_N_K(M,N,K,Kd,A0,B0,C0); 
    }

#else
    // If no SSE2, then repeat the 16 block version to call generic.
    static void multmm_16_16_16(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t ,
        const double* A0, const double* B0, double* C0)
    { generic_multmm_16_16_16(A0,B0,C0); }

    static void multmm_M_16_16(
        const ptrdiff_t M, const ptrdiff_t , const ptrdiff_t ,
        const double* A0, const double* B0, double* C0)
    { generic_multmm_M_16_16(M,A0,B0,C0); }

    static void multmm_16_N_16(
        const ptrdiff_t , const ptrdiff_t N, const ptrdiff_t ,
        const double* A0, const double* B0, double* C0)
    { generic_multmm_16_N_16(N,A0,B0,C0); }

    static void multmm_16_16_K(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t K,
        const double* A0, const double* B0, double* C0)
    { generic_multmm_16_16_K(K,A0,B0,C0); }

    static void multmm_M_N_K(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const double* A0, const double* B0, double* C0)
    { generic_multmm_M_N_K(M,N,K,K,A0,B0,C0); }

#endif

#endif
#endif

#ifdef INST_INT
    //
    // int
    //

    static void multmm_16_16_16(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t ,
        const int* A0, const int* B0, int* C0)
    { generic_multmm_16_16_16(A0,B0,C0); }

    static void multmm_M_16_16(
        const ptrdiff_t M, const ptrdiff_t , const ptrdiff_t ,
        const int* A0, const int* B0, int* C0)
    { generic_multmm_M_16_16(M,A0,B0,C0); }

    static void multmm_16_N_16(
        const ptrdiff_t , const ptrdiff_t N, const ptrdiff_t ,
        const int* A0, const int* B0, int* C0)
    { generic_multmm_16_N_16(N,A0,B0,C0); }

    static void multmm_16_16_K(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t K,
        const int* A0, const int* B0, int* C0)
    { generic_multmm_16_16_K(K,A0,B0,C0); }

    static void multmm_M_N_K(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const int* A0, const int* B0, int* C0)
    { generic_multmm_M_N_K(M,N,K,K,A0,B0,C0); }

#endif

#ifdef INST_LONGDOUBLE

    //
    // long double
    //

    static void multmm_16_16_16(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t ,
        const long double* A0, const long double* B0, long double* C0)
    { generic_multmm_16_16_16(A0,B0,C0); }

    static void multmm_M_16_16(
        const ptrdiff_t M, const ptrdiff_t , const ptrdiff_t ,
        const long double* A0, const long double* B0, long double* C0)
    { generic_multmm_M_16_16(M,A0,B0,C0); }

    static void multmm_16_N_16(
        const ptrdiff_t , const ptrdiff_t N, const ptrdiff_t ,
        const long double* A0, const long double* B0, long double* C0)
    { generic_multmm_16_N_16(N,A0,B0,C0); }

    static void multmm_16_16_K(
        const ptrdiff_t , const ptrdiff_t , const ptrdiff_t K,
        const long double* A0, const long double* B0, long double* C0)
    { generic_multmm_16_16_K(K,A0,B0,C0); }

    static void multmm_M_N_K(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const long double* A0, const long double* B0, long double* C0)
    { generic_multmm_M_N_K(M,N,K,K,A0,B0,C0); }

#endif

    // Next, we have versions of these that always have ix and T
    // as template parameters, which simply call the correct overloaded
    // funciton.  This is designed to get around the C++ limitation of 
    // not being able to partially specialize template functions, only
    // overload them.  Normally overloading is sufficient, but when we
    // need to take an address of a function and we are writing in terms of
    // a generic T, then the overloaded functions cannot be specified.
    // So we specify this one instead.

    template <typename T>
    static void call_multmm_16_16_16(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const T* A0, const T* B0, T* C0)
    { multmm_16_16_16(M,N,K,A0,B0,C0); }

    template <typename T>
    static void call_multmm_M_16_16(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const T* A0, const T* B0, T* C0)
    { multmm_M_16_16(M,N,K,A0,B0,C0); }

    template <typename T>
    static void call_multmm_16_N_16(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const T* A0, const T* B0, T* C0)
    { multmm_16_N_16(M,N,K,A0,B0,C0); }

    template <typename T>
    static void call_multmm_16_16_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const T* A0, const T* B0, T* C0)
    { multmm_16_16_32(M,N,K,A0,B0,C0); }

    template <typename T>
    static void call_multmm_M_16_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const T* A0, const T* B0, T* C0)
    { multmm_M_16_32(M,N,K,A0,B0,C0); }

    template <typename T>
    static void call_multmm_16_N_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const T* A0, const T* B0, T* C0)
    { multmm_16_N_32(M,N,K,A0,B0,C0); }

    template <typename T>
    static void call_multmm_16_16_64(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const T* A0, const T* B0, T* C0)
    { multmm_16_16_64(M,N,K,A0,B0,C0); }

    template <typename T>
    static void call_multmm_M_16_64(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const T* A0, const T* B0, T* C0)
    { multmm_M_16_64(M,N,K,A0,B0,C0); }

    template <typename T>
    static void call_multmm_16_N_64(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const T* A0, const T* B0, T* C0)
    { multmm_16_N_64(M,N,K,A0,B0,C0); }

    template <typename T>
    static void call_multmm_16_16_K(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const T* A0, const T* B0, T* C0)
    { multmm_16_16_K(M,N,K,A0,B0,C0); }

    template <typename T>
    static void call_multmm_M_N_K(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const T* A0, const T* B0, T* C0)
    { multmm_M_N_K(M,N,K,A0,B0,C0); }


} // namespace tmv

#endif 
