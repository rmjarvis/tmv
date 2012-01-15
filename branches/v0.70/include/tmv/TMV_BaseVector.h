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

    template <class T, int N, int A=0> 
    class SmallVector;

    template <class T, int N> 
    class SmallVectorComposite;

    // Things that inherit from this can be assigned to a Vector
    // and explicit constructed to a Vector.
    // But they are not implicitly converted to a Vector.
    template <class T> 
    struct AssignableToVector
    {
        virtual int size() const = 0;
        virtual void assignToV(
            const VectorView<TMV_RealType(T)>& rhs) const = 0;
        virtual void assignToV(
            const VectorView<TMV_ComplexType(T)>& rhs) const = 0;
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

    template <class T, int N, int A> 
    inline std::string TMV_Text(const SmallVector<T,N,A>& )
    { 
        std::ostringstream s;
        s << "SmallVector<"<<TMV_Text(T())<<","<<N<<","<<
            Attrib<A>::vtext()<<">";
        return s.str();
    }

} // namespace tmv

#endif
