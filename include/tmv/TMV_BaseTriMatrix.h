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


#ifndef TMV_BaseTriMatrix_H
#define TMV_BaseTriMatrix_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    template <class T> 
    class GenUpperTriMatrix;

    template <class T> 
    class GenLowerTriMatrix;

    template <class T, int A=0>
    class ConstUpperTriMatrixView;

    template <class T, int A=0>
    class ConstLowerTriMatrixView;

    template <class T, int A=0>
    class UpperTriMatrixView;

    template <class T, int A=0>
    class LowerTriMatrixView;

    template <class T, int A=0>
    class UpperTriMatrix;

    template <class T, int A=0>
    class LowerTriMatrix;

    template <class T> 
    class UpperTriDiv;

    template <class T> 
    class LowerTriDiv;

    template <class T, class Tm> 
    class QuotXU;

    template <class T, class Tm> 
    class QuotXL;

    template <class T1, class T2> 
    void Copy(
        const GenUpperTriMatrix<T1>& m1, UpperTriMatrixView<T2> m2);

    template <class T> 
    struct AssignableToUpperTriMatrix :
        virtual public AssignableToMatrix<T>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        virtual ptrdiff_t size() const = 0;
        virtual DiagType dt() const = 0;
        virtual void assignToU(UpperTriMatrixView<RT> m) const = 0;
        virtual void assignToU(UpperTriMatrixView<CT> m) const = 0;
        virtual inline ~AssignableToUpperTriMatrix() {}
    };

    template <class T> 
    struct AssignableToLowerTriMatrix :
        virtual public AssignableToMatrix<T>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        virtual ptrdiff_t size() const = 0;
        virtual DiagType dt() const = 0;
        virtual void assignToL(LowerTriMatrixView<RT> m) const = 0;
        virtual void assignToL(LowerTriMatrixView<CT> m) const = 0;
        virtual inline ~AssignableToLowerTriMatrix() {}
    };

    template <class T, int A>
    inline std::string TMV_Text(const UpperTriMatrix<T,A>& )
    {
        return std::string("UpperTriMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }
    template <class T, int A>
    inline std::string TMV_Text(const LowerTriMatrix<T,A>& )
    {
        return std::string("LowerTriMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <class T> 
    inline std::string TMV_Text(const GenUpperTriMatrix<T>& )
    {
        return std::string("GenUpperTriMatrix<") + TMV_Text(T()) + ">";
    }
    template <class T> 
    inline std::string TMV_Text(const GenLowerTriMatrix<T>& )
    {
        return std::string("GenLowerTriMatrix<") + TMV_Text(T()) + ">";
    }

    template <class T, int A>
    inline std::string TMV_Text(const ConstUpperTriMatrixView<T,A>& )
    {
        return std::string("ConstUpperTriMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }
    template <class T, int A>
    inline std::string TMV_Text(const ConstLowerTriMatrixView<T,A>& )
    {
        return std::string("ConstLowerTriMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <class T, int A>
    inline std::string TMV_Text(const UpperTriMatrixView<T,A>& )
    {
        return std::string("UpperTriMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }
    template <class T, int A>
    inline std::string TMV_Text(const LowerTriMatrixView<T,A>& )
    {
        return std::string("LowerTriMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

} // namespace tmv

#endif
