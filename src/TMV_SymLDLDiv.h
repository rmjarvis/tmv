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

#ifndef TMV_SymLDLDiv_H
#define TMV_SymLDLDiv_H

#include "tmv/TMV_BaseSymMatrix.h"

namespace tmv {

    // These are in TMV_SymLDLDiv.cpp:
    template <typename T, typename T1> 
    void LDL_LDivEq(
        const GenSymMatrix<T1>& L, const GenVector<T1>& xD,
        const ptrdiff_t* P, MatrixView<T> m);
    template <typename T, typename T1> 
    void LDL_RDivEq(
        const GenSymMatrix<T1>& L, const GenVector<T1>& xD,
        const ptrdiff_t* P, MatrixView<T> m);

    // This is in TMV_SymLDLInverse.cpp:
    template <typename T, typename T1> 
    void LDL_Inverse(
        const GenSymMatrix<T>& L, const GenVector<T>& xD,
        const ptrdiff_t* P, SymMatrixView<T1> sinv);

    // A quick helper class
    template <bool herm, typename T>
    struct Sym2x2_Helper;

    template <typename T>
    struct Sym2x2_Helper<true,T>
    {
        typedef typename Traits<T>::real_type d_type;
        static d_type calculateD(T a, T b, T c)
        { return TMV_REAL(a)*TMV_REAL(b) - TMV_NORM(c); }
    };

    template <typename T>
    struct Sym2x2_Helper<false,T>
    {
        typedef T d_type;
        static d_type calculateD(T a, T b, T c)
        { return a*b-c*c; }
    };

    // These are in TMV_SymLDLPseudo.cpp:
    template <bool herm, typename T> 
    typename Sym2x2_Helper<herm,T>::d_type SymInvert_2x2(T& a, T& b, T& c);

    template <bool herm, typename T, typename T1> 
    void PseudoDiag_LDivEq(
        const GenVector<T1>& D, const GenVector<T1>& xD, 
        MatrixView<T> m);
    template <bool herm, typename T, typename T1> 
    void PseudoDiag_LMultEq(
        const GenVector<T1>& D, const GenVector<T1>& xD, 
        MatrixView<T> m);

    // Specialize disallowed complex combinations:
#define CT std::complex<T>

    template <typename T>
    inline void LDL_LDivEq(
        const GenSymMatrix<CT>& , const GenVector<CT>& ,
        const ptrdiff_t* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void LDL_RDivEq(
        const GenSymMatrix<CT>& , const GenVector<CT>& ,
        const ptrdiff_t* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void LDL_Inverse(
        const GenSymMatrix<CT>& , const GenVector<CT>& ,
        const ptrdiff_t* , SymMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool herm, typename T>
    inline void PseudoDiag_LDivEq(
        const GenVector<CT>& , const GenVector<CT>& , 
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool herm, typename T>
    inline void PseudoDiag_LMultEq(
        const GenVector<CT>& , const GenVector<CT>& , 
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

#undef CT

}

#endif
