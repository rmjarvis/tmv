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


#ifndef TMV_DiagTriArithFunc_H
#define TMV_DiagTriArithFunc_H

#include "tmv/TMV_BaseDiagMatrix.h"
#include "tmv/TMV_BaseTriMatrix.h"

#define CT std::complex<T>

namespace tmv {

    // C (+)= alpha * A * B
    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, 
        const GenUpperTriMatrix<Tb>& B, const UpperTriMatrixView<T>& C);

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, 
        const GenLowerTriMatrix<Tb>& B, const LowerTriMatrixView<T>& C);

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenDiagMatrix<Tb>& B, const UpperTriMatrixView<T>& C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenDiagMatrix<Tb>& B, const LowerTriMatrixView<T>& C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, 
        const GenUpperTriMatrix<Tb>& B, const UpperTriMatrixView<CT>& C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, 
        const GenLowerTriMatrix<Tb>& B, const LowerTriMatrixView<CT>& C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenDiagMatrix<Tb>& B, const UpperTriMatrixView<CT>& C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenDiagMatrix<Tb>& B, const LowerTriMatrixView<CT>& C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenDiagMatrix<Ta>& , 
        const GenUpperTriMatrix<Tb>& , const UpperTriMatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenDiagMatrix<Ta>& , 
        const GenLowerTriMatrix<Tb>& , const LowerTriMatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenUpperTriMatrix<Ta>& , 
        const GenDiagMatrix<Tb>& , const UpperTriMatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenLowerTriMatrix<Ta>& , 
        const GenDiagMatrix<Tb>& , const LowerTriMatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

} // namespace tmv

#undef CT

#endif
