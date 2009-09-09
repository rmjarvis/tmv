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

//
//
// This file defines the StackArray class used by SmallVector and SmallMatrix.
// It basically emulates a normal C array on the stack:
// T p[N];
// However, it uses the heap when N is large.

#ifndef StackArray_H
#define StackArray_H

const int TMV_MaxStack = 1024;

namespace tmv
{
  template <class T, int N, bool bigN> struct StackArray2;

  template <class T, int N>
  struct StackArray2<T,N,false>
  {
    T p[N];
    inline StackArray2() {}
    inline ~StackArray2() {}
  };

  template <class T, int N>
  struct StackArray2<T,N,true>
  {
    T*const p;
    inline StackArray2() : p(new T[N]) {}
    inline ~StackArray2() { delete [] p; }
  };

  template <class T, int N>
  class StackArray
  {
  private :
    StackArray2<T,N,(N>TMV_MaxStack)> p;

    inline StackArray& operator=(StackArray& p2);
    inline StackArray(const StackArray& p2);

  public :
    inline StackArray() {}
    inline ~StackArray() {}

    inline T& operator*() { return *(p.p); }
    inline T* operator->() { return (p.p); }
    inline operator T*() { return (p.p); }

    inline const T& operator*() const { return *(p.p); }
    inline const T* operator->() const { return (p.p); }
    inline operator const T*() const { return (p.p); }
  };

} // namespace tmv

#endif
