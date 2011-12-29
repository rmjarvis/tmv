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

    template <class T, IndexStyle I=CStyle> 
    class ConstVectorView;

    template <class T, IndexStyle I=CStyle> 
    class VectorView;

    template <class T, IndexStyle I=CStyle> 
    class Vector;

    template <class T, int N, IndexStyle I=CStyle> 
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
        virtual void assignToV(const VectorView<TMV_RealType(T)>& rhs) const = 0;
        virtual void assignToV(const VectorView<TMV_ComplexType(T)>& rhs) const = 0;
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

    inline std::string TMV_Text(IndexStyle I)
    { return I == CStyle ? "CStyle" : "FortranStyle"; }

    template <class T, IndexStyle I> 
    inline std::string TMV_Text(const Vector<T,I>& )
    { return std::string("Vector<")+TMV_Text(T())+","+TMV_Text(I)+">"; }

    template <class T> 
    inline std::string TMV_Text(const GenVector<T>& v)
    { return std::string("GenVector<")+TMV_Text(T())+","+TMV_Text(v.ct())+">"; }

    template <class T, IndexStyle I> 
    inline std::string TMV_Text(const ConstVectorView<T,I>& v)
    { 
        return std::string("ConstVectorView<")+TMV_Text(T())+","+TMV_Text(I)+","+
            TMV_Text(v.ct())+">"; 
    }

    template <class T, IndexStyle I> 
    inline std::string TMV_Text(const VectorView<T,I>& v)
    { 
        return std::string("VectorView<")+TMV_Text(T())+","+TMV_Text(I)+","
            +TMV_Text(v.ct())+">"; 
    }

} // namespace tmv

#endif
