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


#ifndef TMV_Scaling_H
#define TMV_Scaling_H

#include "TMV_Base.h"

namespace tmv {

  // This structure implements a scaling factor with specialization
  // for x = 1 or -1
  // The ix is a flag integer that can be 1, -1 or 0.
  // When ix is 1 or -1 then x is 1 or -1 and isn't stored.
  // When ix is 0, then x is stored and acts normally.
  template <int ix, class T> struct Scaling;

  template <class T>
  struct Scaling<0,T>
  {
#ifdef TMV_DEBUG
    T x;
#else
    const T x;
#endif
    inline Scaling(const T _x) : x(_x) {}
    inline ~Scaling() {
#ifdef TMV_DEBUG
      x = T(999);
#endif
    }
    operator T() const { return x; }
  };
  // We enforce that T is real for 1 and -1 to avoid accidentally 
  // having two different template paths produce functionally 
  // identical code.  
  // Always use real T for 1 and -1.
  template <class T>
  struct Scaling<1,T> // ix = 1
  {
    inline Scaling() { TMVStaticAssert(Traits<T>::isreal); }
    inline Scaling(const T _x) 
    { TMVStaticAssert(Traits<T>::isreal); TMVAssert(_x == T(1));  }
    operator T() const { return T(1); }
  };
  template <class T>
  struct Scaling<-1,T> // ix = -1
  {
    inline Scaling() { TMVStaticAssert(Traits<T>::isreal); }
    inline Scaling(const T _x) 
    { TMVStaticAssert(Traits<T>::isreal); TMVAssert(_x == T(-1));  }
    operator T() const { return T(-1); }
  };

  template <int ix, class T> 
  inline Scaling<ix,T> TMV_CONJ(const Scaling<ix,T>& x)
  { return x; }
  template <class T> 
  inline Scaling<0,T> TMV_CONJ(const Scaling<0,T>& x)
  { return Scaling<0,T>(TMV_CONJ(x.x)); }

#define RT typename Traits<T>::real_type
#define CT typename Traits<T>::complex_type
#define CCT ConjRef<CT>

  // Define how a Scaling object multiplies a regular number:
  template <int ix, class T>
  inline RT& operator*=(RT& y, const Scaling<ix,T>& x)
  { y *= x.x; return y; }
  template <class T>
  inline RT& operator*=(RT& y, const Scaling<1,T>& x)
  { return y; }
  template <class T>
  inline RT& operator*=(RT& y, const Scaling<-1,T>& x)
  { y = -y; return y; }
  template <int ix, class T>
  inline RT& operator/=(RT& y, const Scaling<ix,T>& x)
  { y /= x.x; return y; }
  template <class T>
  inline RT& operator/=(RT& y, const Scaling<1,T>& x)
  { return y; }
  template <class T>
  inline RT& operator/=(RT& y, const Scaling<-1,T>& x)
  { y = -y; return y; }

  template <int ix, class T>
  inline CT& operator*=(CT& y, const Scaling<ix,T>& x)
  { y *= x.x; return y; }
  template <class T>
  inline CT& operator*=(CT& y, const Scaling<1,T>& x)
  { return y; }
  template <class T>
  inline CT& operator*=(CT& y, const Scaling<-1,T>& x)
  { y = -y; return y; }
  template <int ix, class T>
  inline CT& operator/=(CT& y, const Scaling<ix,T>& x)
  { y /= x.x; return y; }
  template <class T>
  inline CT& operator/=(CT& y, const Scaling<1,T>& x)
  { return y; }
  template <class T>
  inline CT& operator/=(CT& y, const Scaling<-1,T>& x)
  { y = -y; return y; }

  template <int ix, class T>
  inline CCT& operator*=(CCT& y, const Scaling<ix,T>& x)
  { y *= x.x; return y; }
  template <class T>
  inline CCT& operator*=(CCT& y, const Scaling<1,T>& x)
  { return y; }
  template <class T>
  inline CCT& operator*=(CCT& y, const Scaling<-1,T>& x)
  { y = -y; return y; }
  template <int ix, class T>
  inline CCT& operator/=(CCT& y, const Scaling<ix,T>& x)
  { y /= x.x; return y; }
  template <class T>
  inline CCT& operator/=(CCT& y, const Scaling<1,T>& x)
  { return y; }
  template <class T>
  inline CCT& operator/=(CCT& y, const Scaling<-1,T>& x)
  { y = -y; return y; }


  template <int ix, class T>
  inline T operator*(const Scaling<ix,T>& x, const RT& y)
  { return x.x * y; }
  template <class T>
  inline RT operator*(const Scaling<1,T>& x, const RT& y)
  { return y; }
  template <class T>
  inline RT operator*(const Scaling<-1,T>& x, const RT& y)
  { return -y; }

  template <int ix, class T>
  inline T operator/(const RT& y, const Scaling<ix,T>& x)
  { return y / x.x; }
  template <class T>
  inline RT operator/(const RT& y, const Scaling<1,T>& x)
  { return y; }
  template <class T>
  inline RT operator/(const RT& y, const Scaling<-1,T>& x)
  { return -y; }

  template <int ix, class T>
  inline CT operator*(const Scaling<ix,T>& x, const CT& y)
  { return ZProd<false,false>::prod(x.x,y); }
  template <class T>
  inline CT operator*(const Scaling<1,T>& x, const CT& y)
  { return y; }
  template <class T>
  inline CT operator*(const Scaling<-1,T>& x, const CT& y)
  { return -y; }

  template <int ix, class T>
  inline CT operator/(const CT& y, const Scaling<ix,T>& x)
  { return y / x.x; }
  template <class T>
  inline CT operator/(const CT& y, const Scaling<1,T>& x)
  { return y; }
  template <class T>
  inline CT operator/(const CT& y, const Scaling<-1,T>& x)
  { return -y; }

  template <int ix, class T>
  inline CT operator*(const Scaling<ix,T>& x, const CCT& y)
  { return x * CT(y); }
  template <class T>
  inline CT operator*(const Scaling<1,T>& x, const CCT& y)
  { return CT(y); }
  template <class T>
  inline CT operator*(const Scaling<-1,T>& x, const CCT& y)
  { return -CT(y); }

  template <int ix, class T>
  inline CT operator/(const CCT& y, const Scaling<ix,T>& x)
  { return CT(y) / x.x; }
  template <class T>
  inline CT operator/(const CCT& y, const Scaling<1,T>& x)
  { return CT(y); }
  template <class T>
  inline CT operator/(const CCT& y, const Scaling<-1,T>& x)
  { return -CT(y); }

  // Other * and / remap into the above operations:
  template <int ix, class T>
  inline T operator*(const RT& y, const Scaling<ix,T>& x)
  { return x * y; }
  template <int ix, class T>
  inline RT operator*(const RT& y, const Scaling<1,T>& x)
  { return x * y; }
  template <int ix, class T>
  inline RT operator*(const RT& y, const Scaling<-1,T>& x)
  { return x * y; }
  template <int ix, class T>
  inline T operator/(const Scaling<ix,T>& x, const RT& y)
  { return T(x) / y; }
  template <int ix, class T>
  inline RT operator/(const Scaling<1,T>& x, const RT& y)
  { return RT(1) / y; }
  template <int ix, class T>
  inline RT operator/(const Scaling<-1,T>& x, const RT& y)
  { return RT(-1) / y; }

  template <int ix, class T>
  inline CT operator*(const CT& y, const Scaling<ix,T>& x)
  { return x * y; }
  template <int ix, class T>
  inline CT operator/(const Scaling<ix,T>& x, const CT& y)
  { return T(x) / y; }

  template <int ix, class T>
  inline CT operator*(const CCT& y, const Scaling<ix,T>& x)
  { return x * CT(y); }
  template <int ix, class T>
  inline CT operator/(const Scaling<ix,T>& x, const CCT& y)
  { return T(x) / CT(y); }

#undef RT
#undef CT
#undef CCT

  template <int ix, class T>
  inline Scaling<-ix,T> operator-(const Scaling<ix,T>& x)
  { return Scaling<-ix,T>(-T(x)); }

  // Scaling * Scaling:
  template <class T1, class T2>
  inline typename Traits2<T1,T2>::type operator*(
      const Scaling<0,T1>& x, const Scaling<0,T2>& y)
  { return ZProd<false,false>::prod(x.x,y.x); }
  template <class T1, class T2>
  inline T2 operator*(const Scaling<1,T1>& x, const Scaling<0,T2>& y)
  { return y.x; }
  template <class T1, class T2>
  inline T2 operator*(const Scaling<-1,T1>& x, const Scaling<0,T2>& y)
  { return -y.x; }

  template <class T1, class T2>
  inline T1 operator*(const Scaling<0,T1>& x, const Scaling<1,T2>& y)
  { return x; }
  template <class T1, class T2>
  inline Scaling<1,T2> operator*(
      const Scaling<1,T1>& x, const Scaling<1,T2>& y)
  { return Scaling<1,T2>(); }
  template <class T1, class T2>
  inline Scaling<-1,T2> operator*(
      const Scaling<-1,T1>& x, const Scaling<1,T2>& y)
  { return Scaling<-1,T2>(); }

  template <class T1, class T2>
  inline T1 operator*(const Scaling<0,T1>& x, const Scaling<-1,T2>& y)
  { return -x.x; }
  template <class T1, class T2>
  inline Scaling<-1,T2> operator*(
      const Scaling<1,T1>& x, const Scaling<-1,T2>& y)
  { return Scaling<-1,T2>(); }
  template <class T1, class T2>
  inline Scaling<1,T2> operator*(
      const Scaling<-1,T1>& x, const Scaling<-1,T2>& y)
  { return Scaling<1,T2>(); }


  // Scaling / Scaling:
  template <class T1, class T2>
  inline typename Traits2<T1,T2>::type operator/(
      const Scaling<0,T1>& x, const Scaling<0,T2>& y)
  { return x.x/y.x; }
  template <class T1, class T2>
  inline T2 operator/(const Scaling<1,T1>& x, const Scaling<0,T2>& y)
  { return T1(1)/y.x; }
  template <class T1, class T2>
  inline T2 operator/(const Scaling<-1,T1>& x, const Scaling<0,T2>& y)
  { return T1(-1)/y.x; }

  template <class T1, class T2>
  inline T1 operator/(const Scaling<0,T1>& x, const Scaling<1,T2>& y)
  { return x; }
  template <class T1, class T2>
  inline Scaling<1,T2> operator/(
      const Scaling<1,T1>& x, const Scaling<1,T2>& y)
  { return Scaling<1,T2>(); }
  template <class T1, class T2>
  inline Scaling<-1,T2> operator/(
      const Scaling<-1,T1>& x, const Scaling<1,T2>& y)
  { return Scaling<-1,T2>(); }

  template <class T1, class T2>
  inline T1 operator/(const Scaling<0,T1>& x, const Scaling<-1,T2>& y)
  { return -x.x; }
  template <class T1, class T2>
  inline Scaling<-1,T2> operator/(
      const Scaling<1,T1>& x, const Scaling<-1,T2>& y)
  { return Scaling<-1,T2>(); }
  template <class T1, class T2>
  inline Scaling<1,T2> operator/(
      const Scaling<-1,T1>& x, const Scaling<-1,T2>& y)
  { return Scaling<1,T2>(); }


};

#endif
