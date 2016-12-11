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


#ifndef TMV_BaseSymBandMatrix_H
#define TMV_BaseSymBandMatrix_H

#include "tmv/TMV_BaseSymMatrix.h"
#include "tmv/TMV_BaseBandMatrix.h"

namespace tmv {

    template <typename T>
    class GenSymBandMatrix;

    template <typename T, int A=0>
    class ConstSymBandMatrixView;

    template <typename T, int A=0>
    class SymBandMatrixView;

    template <typename T, int A=0>
    class SymBandMatrix;

    template <typename T, int A=0>
    class HermBandMatrix;

    template <typename T>
    class HermBandCHDiv;

    template <typename T>
    class HermBandSVDiv;

    template <typename T>
    class SymBandSVDiv;

    template <typename T, typename Tm>
    class QuotXsB;

    template <typename T, typename T1, typename T2>
    class ProdBB;

    template <typename T>
    struct AssignableToSymBandMatrix :
        virtual public AssignableToSymMatrix<T>,
        virtual public AssignableToBandMatrix<T>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        virtual void assignTosB(SymBandMatrixView<RT> m) const = 0;
        virtual void assignTosB(SymBandMatrixView<CT> m) const = 0;
        virtual inline ~AssignableToSymBandMatrix() {}
    };


    template <typename T, int A>
    inline std::string TMV_Text(const SymBandMatrix<T,A>& )
    {
        return std::string("SymBandMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <typename T, int A>
    inline std::string TMV_Text(const HermBandMatrix<T,A>& )
    {
        return std::string("HermBandMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <typename T>
    inline std::string TMV_Text(const GenSymBandMatrix<T>& )
    {
        return std::string("GenSymBandMatrix<") + TMV_Text(T()) + ">";
    }

    template <typename T, int A>
    inline std::string TMV_Text(const ConstSymBandMatrixView<T,A>& )
    {
        return std::string("ConstSymBandMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <typename T, int A>
    inline std::string TMV_Text(const SymBandMatrixView<T,A>& )
    {
        return std::string("SymBandMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

}

#endif
