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


#ifndef TMV_BaseSymMatrix_H
#define TMV_BaseSymMatrix_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    template <typename T>
    class GenSymMatrix;

    template <typename T, int A=0>
    class ConstSymMatrixView;

    template <typename T, int A=0>
    class SymMatrixView;

    template <typename T, int A=0>
    class SymMatrix;

    template <typename T, int A=0>
    class HermMatrix;

    template <typename T>
    class SymDivider;

    template <typename T>
    class SymLDLDiv;

    template <typename T>
    class HermCHDiv;

    template <typename T>
    class HermSVDiv;

    template <typename T>
    class SymSVDiv;

    template <typename T, typename Tm>
    class QuotXS;

    template <typename T>
    struct AssignableToSymMatrix :
        virtual public AssignableToMatrix<T>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        virtual ptrdiff_t size() const = 0;
        virtual SymType sym() const = 0;
        inline bool issym() const { return isReal(T()) || (sym() == Sym); }
        inline bool isherm() const
        { return isReal(T()) || (sym() == Herm); }
        virtual void assignToS(SymMatrixView<RT> m) const = 0;
        virtual void assignToS(SymMatrixView<CT> m) const = 0;
        virtual inline ~AssignableToSymMatrix() {}
    };

    template <typename T, int A>
    inline std::string TMV_Text(const SymMatrix<T,A>& )
    {
        return std::string("SymMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <typename T, int A>
    inline std::string TMV_Text(const HermMatrix<T,A>& )
    {
        return std::string("HermMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <typename T>
    inline std::string TMV_Text(const GenSymMatrix<T>& )
    {
        return std::string("GenSymMatrix<") + TMV_Text(T()) + ">";
    }

    template <typename T, int A>
    inline std::string TMV_Text(const ConstSymMatrixView<T,A>& )
    {
        return std::string("ConstSymMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <typename T, int A>
    inline std::string TMV_Text(const SymMatrixView<T,A>& )
    {
        return std::string("SymMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

}

#endif
