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
#ifndef TMV_MultMM_H
#define TMV_MultMM_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

#ifdef TMV_BLOCKSIZE
    const int MM_BLOCKSIZE = TMV_BLOCKSIZE;
#else
    const int MM_BLOCKSIZE = 64;
#endif

#ifdef TMV_BLOCKSIZE
    const int MMCCC_BLOCKSIZE_M = 4*TMV_BLOCKSIZE;
    const int MMCCC_BLOCKSIZE_N = TMV_BLOCKSIZE;
#else
    const int MMCCC_BLOCKSIZE_M = 4*64;
    const int MMCCC_BLOCKSIZE_N = 64;
#endif

    template <bool ca, bool cb, class T, class Ta, class Tb> 
    void RecursiveCCCMultMM(
        const int M, const int N, const int K,
        const Ta* A, const Tb* B, T* C,
        const int Ask, const int Bsj, const int Csj);

    template <bool add, bool ca, bool cb, class T, class Ta, class Tb> 
    void DoCCCMultMM( 
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, const MatrixView<T>& C);

    template <bool add, class T, class Ta, class Tb> 
    void CCCMultMM(const T alpha, const GenMatrix<Ta>& A,
                   const GenMatrix<Tb>& B, const MatrixView<T>& C);

#ifdef TMV_BLOCKSIZE
    const int MMCRC_BLOCKSIZE_M = 4*TMV_BLOCKSIZE;
    const int MMCRC_BLOCKSIZE_N = TMV_BLOCKSIZE;
#else
    const int MMCRC_BLOCKSIZE_M = 4*64;
    const int MMCRC_BLOCKSIZE_N = 64;
#endif

    template <bool ca, bool cb, class T, class Ta, class Tb> 
    void RecursiveCRCMultMM(
        const int M, const int N, const int K,
        const Ta* A, const Tb* B, T* C,
        const int Ask, const int Bsk, const int Csj);

    template <bool add, bool ca, bool cb, class T, class Ta, class Tb> 
    void DoCRCMultMM(const T alpha, const GenMatrix<Ta>& A,
                     const GenMatrix<Tb>& B, const MatrixView<T>& C);

    template <bool add, class T, class Ta, class Tb> 
    void CRCMultMM(const T alpha, const GenMatrix<Ta>& A,
                   const GenMatrix<Tb>& B, const MatrixView<T>& C);

#ifdef TMV_BLOCKSIZE
    const int MMRCC_BLOCKSIZE_K = 4*TMV_BLOCKSIZE;
    const int MMRCC_BLOCKSIZE_MN = TMV_BLOCKSIZE;
#else
    const int MMRCC_BLOCKSIZE_K = 4*64;
    const int MMRCC_BLOCKSIZE_MN = 64;
#endif

    template <bool ca, bool cb, class T, class Ta, class Tb>
    void RecursiveRCCMultMM(
        const int M, const int N, const int K,
        const Ta* A, const Tb* B, T* C,
        const int Asi, const int Bsj, const int Csj);

    template <bool add, bool ca, bool cb, class T, class Ta, class Tb> 
    void DoRCCMultMM(const T alpha, const GenMatrix<Ta>& A,
                     const GenMatrix<Tb>& B, const MatrixView<T>& C);

    template <bool add, class T, class Ta, class Tb> 
    void RCCMultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C);

}

#endif
