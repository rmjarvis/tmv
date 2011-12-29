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
#ifndef TMV_SymLDLDiv_H
#define TMV_SymLDLDiv_H

#include "tmv/TMV_BaseSymMatrix.h"

namespace tmv {

    // These are in TMV_SymLDLDiv.cpp:
    template <class T, class T1> 
    void LDL_LDivEq(
        const GenSymMatrix<T1>& L, const GenVector<T1>& xD,
        const int* P, const MatrixView<T>& m);
    template <class T, class T1> 
    void LDL_RDivEq(
        const GenSymMatrix<T1>& L, const GenVector<T1>& xD,
        const int* P, const MatrixView<T>& m);

    // This is in TMV_SymLDLInverse.cpp:
    template <class T, class T1> 
    void LDL_Inverse(
        const GenSymMatrix<T>& L, const GenVector<T>& xD,
        const int* P, const SymMatrixView<T1>& sinv);

    // A quick helper class
    template <bool herm, class T>
    struct Sym2x2_Helper;

    template <class T>
    struct Sym2x2_Helper<true,T>
    {
        typedef typename Traits<T>::real_type d_type;
        static d_type calculateD(T a, T b, T c)
        { return TMV_REAL(a)*TMV_REAL(b) - TMV_NORM(c); }
    };

    template <class T>
    struct Sym2x2_Helper<false,T>
    {
        typedef T d_type;
        static d_type calculateD(T a, T b, T c)
        { return a*b-c*c; }
    };

    // These are in TMV_SymLDLPseudo.cpp:
    template <bool herm, class T> 
    typename Sym2x2_Helper<herm,T>::d_type SymInvert_2x2(T& a, T& b, T& c);

    template <bool herm, class T, class T1> 
    void PseudoDiag_LDivEq(
        const GenVector<T1>& D, const GenVector<T1>& xD, 
        const MatrixView<T>& m);
    template <bool herm, class T, class T1> 
    void PseudoDiag_LMultEq(
        const GenVector<T1>& D, const GenVector<T1>& xD, 
        const MatrixView<T>& m);

    // Specialize disallowed complex combinations:
#define CT std::complex<T>

    template <class T>
    inline void LDL_LDivEq(
        const GenSymMatrix<CT>& , const GenVector<CT>& ,
        const int* , const MatrixView<T>& )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void LDL_RDivEq(
        const GenSymMatrix<CT>& , const GenVector<CT>& ,
        const int* , const MatrixView<T>& )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void LDL_Inverse(
        const GenSymMatrix<CT>& , const GenVector<CT>& ,
        const int* , const SymMatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

    template <bool herm, class T>
    inline void PseudoDiag_LDivEq(
        const GenVector<CT>& , const GenVector<CT>& , 
        const MatrixView<T>& )
    { TMVAssert(TMV_FALSE); }
    template <bool herm, class T>
    inline void PseudoDiag_LMultEq(
        const GenVector<CT>& , const GenVector<CT>& , 
        const MatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

#undef CT

}

#endif
