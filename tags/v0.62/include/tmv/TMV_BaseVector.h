///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
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
    virtual size_t size() const = 0;
    virtual void AssignToV(const VectorView<RealType(T)>& rhs) const = 0;
    virtual void AssignToV(const VectorView<ComplexType(T)>& rhs) const = 0;
    virtual inline ~AssignableToVector() {}
  };

  inline std::string Text(ADType ad)
  { return ad == ASCEND ? "Ascend" : "Descend"; }

  inline std::string Text(COMPType comp)
  { 
    return comp == REAL_COMP ? "Real" : comp == ABS_COMP ? "Abs" :
      comp == IMAG_COMP ? "Imag" : "Arg"; 
  }

  inline std::string Text(IndexStyle I)
  { return I == CStyle ? "CStyle" : "FortranStyle"; }

  template <class T, IndexStyle I> 
  inline std::string TypeText(const Vector<T,I>& )
  { return std::string("Vector<")+TypeText(T())+","+Text(I)+">"; }

  template <class T> 
  inline std::string TypeText(const GenVector<T>& v)
  { return std::string("GenVector<")+TypeText(T())+","+Text(v.ct())+">"; }

  template <class T, IndexStyle I> 
  inline std::string TypeText(
      const ConstVectorView<T,I>& v)
  { 
    return std::string("ConstVectorView<")+TypeText(T())+","+Text(I)+","+
      Text(v.ct())+">"; 
  }

  template <class T, IndexStyle I> 
  inline std::string TypeText(
      const VectorView<T,I>& v)
  { 
    return std::string("VectorView<")+TypeText(T())+","+Text(I)+","
      +Text(v.ct())+">"; 
  }

} // namespace tmv

#endif
