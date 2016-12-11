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


#ifndef TMV_DiagTriArithFunc_H
#define TMV_DiagTriArithFunc_H

#include "tmv/TMV_BaseDiagMatrix.h"
#include "tmv/TMV_BaseTriMatrix.h"

#define CT std::complex<T>

namespace tmv {

    // C (+)= alpha * A * B
    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<T> C);

    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, LowerTriMatrixView<T> C);

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenDiagMatrix<Tb>& B, UpperTriMatrixView<T> C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenDiagMatrix<Tb>& B, LowerTriMatrixView<T> C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, LowerTriMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenDiagMatrix<Tb>& B, UpperTriMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenDiagMatrix<Tb>& B, LowerTriMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenDiagMatrix<Ta>& ,
        const GenUpperTriMatrix<Tb>& , UpperTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenDiagMatrix<Ta>& ,
        const GenLowerTriMatrix<Tb>& , LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenUpperTriMatrix<Ta>& ,
        const GenDiagMatrix<Tb>& , UpperTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenLowerTriMatrix<Ta>& ,
        const GenDiagMatrix<Tb>& , LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

} // namespace tmv

#undef CT

#endif
