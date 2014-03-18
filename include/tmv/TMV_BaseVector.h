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


#ifndef TMV_BaseVector_H
#define TMV_BaseVector_H

#include "tmv/TMV_Base.h"

namespace tmv {

    template <class T> 
    class GenVector;

    template <class T, int A=0> 
    class ConstVectorView;

    template <class T, int A=0> 
    class VectorView;

    template <class T, int A=0> 
    class Vector;

    template <class T, ptrdiff_t N, int A=0> 
    class SmallVector;

    template <class T, ptrdiff_t N> 
    class SmallVectorComposite;

    // Things that inherit from this can be assigned to a Vector
    // and explicit constructed to a Vector.
    // But they are not implicitly converted to a Vector.
    template <class T> 
    struct AssignableToVector
    {
        virtual ptrdiff_t size() const = 0;
        virtual void assignToV(VectorView<TMV_RealType(T)> rhs) const = 0;
        virtual void assignToV(VectorView<TMV_ComplexType(T)> rhs) const = 0;
        virtual inline ~AssignableToVector() {}
    };

    // TMV_Text definitions:

    inline std::string TMV_Text(ADType ad)
    { return ad == Ascend ? "Ascend" : "Descend"; }

    inline std::string TMV_Text(CompType comp)
    { 
        return comp == RealComp ? "Real" : comp == AbsComp ? "Abs" :
            comp == ImagComp ? "Imag" : "Arg"; 
    }

    template <class T, int A> 
    inline std::string TMV_Text(const Vector<T,A>& )
    { return std::string("Vector<")+TMV_Text(T())+","+Attrib<A>::vtext()+">"; }

    template <class T> 
    inline std::string TMV_Text(const GenVector<T>& )
    { return std::string("GenVector<")+TMV_Text(T())+">"; }

    template <class T, int A> 
    inline std::string TMV_Text(const ConstVectorView<T,A>& )
    { 
        return std::string("ConstVectorView<")+TMV_Text(T())+","+
            Attrib<A>::vtext()+">"; 
    }

    template <class T, int A> 
    inline std::string TMV_Text(const VectorView<T,A>& )
    { 
        return std::string("VectorView<")+TMV_Text(T())+","+
            Attrib<A>::vtext()+">"; 
    }

    template <class T, ptrdiff_t N, int A> 
    inline std::string TMV_Text(const SmallVector<T,N,A>& )
    { 
        std::ostringstream s;
        s << "SmallVector<"<<TMV_Text(T())<<","<<N<<","<<
            Attrib<A>::vtext()<<">";
        return s.str();
    }

} // namespace tmv

#endif
