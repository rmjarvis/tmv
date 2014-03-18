///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2014                                                 //
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

    template <class T> 
    class GenSymBandMatrix;

    template <class T, int A=0>
    class ConstSymBandMatrixView;

    template <class T, int A=0>
    class SymBandMatrixView;

    template <class T, int A=0>
    class SymBandMatrix;

    template <class T, int A=0>
    class HermBandMatrix;

    template <class T> 
    class HermBandCHDiv;

    template <class T> 
    class HermBandSVDiv;

    template <class T> 
    class SymBandSVDiv;

    template <class T, class Tm> 
    class QuotXsB;

    template <class T, class T1, class T2> 
    class ProdBB;

    template <class T> 
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


    template <class T, int A>
    inline std::string TMV_Text(const SymBandMatrix<T,A>& )
    {
        return std::string("SymBandMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <class T, int A>
    inline std::string TMV_Text(const HermBandMatrix<T,A>& )
    {
        return std::string("HermBandMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <class T> 
    inline std::string TMV_Text(const GenSymBandMatrix<T>& )
    {
        return std::string("GenSymBandMatrix<") + TMV_Text(T()) + ">";
    }

    template <class T, int A>
    inline std::string TMV_Text(const ConstSymBandMatrixView<T,A>& )
    {
        return std::string("ConstSymBandMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <class T, int A>
    inline std::string TMV_Text(const SymBandMatrixView<T,A>& )
    {
        return std::string("SymBandMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

}

#endif
