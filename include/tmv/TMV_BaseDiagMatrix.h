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


#ifndef TMV_BaseDiagMatrix_H
#define TMV_BaseDiagMatrix_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    template <typename T>
    class GenDiagMatrix;

    template <typename T, int A=0>
    class ConstDiagMatrixView;

    template <typename T, int A=0>
    class DiagMatrixView;

    template <typename T, int A=0>
    class DiagMatrix;

    template <typename T, typename Tm>
    class QuotXD;

    template <typename T>
    struct AssignableToDiagMatrix :
        virtual public AssignableToMatrix<T>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        virtual ptrdiff_t size() const = 0;
        virtual void assignToD(DiagMatrixView<RT> m) const = 0;
        virtual void assignToD(DiagMatrixView<CT> m) const = 0;
        virtual inline ~AssignableToDiagMatrix() {}
    };

    template <typename T, int A>
    inline std::string TMV_Text(const DiagMatrix<T,A>& )
    {
        return std::string("DiagMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::vtext() + ">";
    }

    template <typename T>
    inline std::string TMV_Text(const GenDiagMatrix<T>& )
    {
        return std::string("GenDiagMatrix<") + TMV_Text(T()) + ">";
    }

    template <typename T, int A>
    inline std::string TMV_Text(const ConstDiagMatrixView<T,A>& )
    {
        return std::string("ConstDiagMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::vtext() + ">";
    }

    template <typename T, int A>
    inline std::string TMV_Text(const DiagMatrixView<T,A>& )
    {
        return std::string("DiagMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::vtext() + ">";
    }

} // namespace tmv

#endif
