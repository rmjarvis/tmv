///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
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


// Things that need to be #defined on entry:
// (The values for a normal Matrix+Matrix are given)
//
// SUMMM	SumMM
// GENMATRIX1	GenMatrix
// GENMATRIX2	GenMatrix
// PRODXM1	ProdXM
// PRODXM2	ProdXM

#ifndef X3
#define X3
#endif

#ifndef Y
#define Y
#endif

#ifndef CT
#define CT std::complex<T>
#endif
#ifndef CCT
#define CCT ConjRef<std::complex<T> >
#endif
#ifndef VCT
#define VCT VarConjRef<std::complex<T> >
#endif

// -(x*m+x*m)
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator-(
    const SUMMM<T,T1,T2 X3>& smm)
{ 
  return SUMMM<T,T1,T2 X3>(-smm.GetX1(),smm.GetM1(),-smm.GetX2,
      smm.GetM2); 
}

// x*(x*m+x*m)
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator*(
    const T x, const SUMMM<T,T1,T2 X3>& smm)
{
  return SUMMM<T,T1,T2 X3>(smm.GetX1()*x,smm.GetM1(),smm.GetX2()*x,
      smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const T x, const SUMMM<CT,T1,T2 X3>& smm)
{
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()*x,smm.GetM1(),smm.GetX2()*x,
      smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const CT x, const SUMMM<T,T,T X3>& smm)
{
  return SUMMM<CT,T,T X3>(smm.GetX1()*x,smm.GetM1(),smm.GetX2()*x,
      smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const CCT x, const SUMMM<T,T,T X3>& smm)
{ 
  return SUMMM<CT,T,T X3>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const VCT x, const SUMMM<T,T,T X3>& smm)
{ 
  return SUMMM<CT,T,T X3>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const CCT x, const SUMMM<CT,T1,T2 X3>& smm)
{ 
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const VCT x, const SUMMM<CT,T1,T2 X3>& smm)
{ 
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}


// (x*m+x*m)*x
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator*(
    const SUMMM<T,T1,T2 X3>& smm, const T x)
{
  return SUMMM<T,T1,T2 X3>(smm.GetX1()*x,smm.GetM1(),smm.GetX2()*x,
      smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM<CT,T1,T2 X3>& smm, const T x)
{
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()*x,smm.GetM1(),smm.GetX2()*x,
      smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM<T,T,T X3>& smm, const CT x)
{
  return SUMMM<CT,T,T X3>(smm.GetX1()*x,smm.GetM1(),smm.GetX2()*x,
      smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM<T,T,T X3>& smm, const CCT x)
{ 
  return SUMMM<CT,T,T X3>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM<T,T,T X3>& smm, const VCT x)
{ 
  return SUMMM<CT,T,T X3>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM<CT,T1,T2 X3>& smm, const CCT x)
{ 
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM<CT,T1,T2 X3>& smm, const VCT x)
{ 
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

// (x*m+x*m)/x
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator/(
    const SUMMM<T,T1,T2 X3>& smm, const T x)
{
  return SUMMM<T,T1,T2 X3>(smm.GetX1()/x,smm.GetM1(),smm.GetX2()/x,
      smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM<CT,T1,T2 X3>& smm, const T x)
{
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()/x,smm.GetM1(),smm.GetX2()/x,
      smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM<T,T,T X3>& smm, const CT x)
{
  return SUMMM<CT,T,T X3>(smm.GetX1()/x,smm.GetM1(),smm.GetX2()/x,
      smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM<T,T,T X3>& smm, const CCT x)
{ 
  return SUMMM<CT,T,T X3>(smm.GetX1()/CT(x),smm.GetM1(),
      smm.GetX2()/CT(x),smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM<T,T,T X3>& smm, const VCT x)
{ 
  return SUMMM<CT,T,T X3>(smm.GetX1()/CT(x),smm.GetM1(),
      smm.GetX2()/CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM<CT,T1,T2 X3>& smm, const CCT x)
{ 
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()/CT(x),smm.GetM1(),
      smm.GetX2()/CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM<CT,T1,T2 X3>& smm, const VCT x)
{ 
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()/CT(x),smm.GetM1(),
      smm.GetX2()/CT(x),smm.GetM2()); 
}

#ifdef SUMMM_1_1

// -(x*m+x*m)
template <class T, class T1, class T2 Y> inline SUMMM_x_m1<T,T1,T2 X3> operator-(
    const SUMMM_1_1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(T(-1),smm.GetM1(),T(-1),smm.GetM2); }

// x*(x*m+x*m)
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator*(
    const T x, const SUMMM_1_1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(x,smm.GetM1(),x,smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const T x, const SUMMM_1_1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const CT x, const SUMMM_1_1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(x,smm.GetM1(),x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const CCT x, const SUMMM_1_1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const VCT x, const SUMMM_1_1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const CCT x, const SUMMM_1_1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const VCT x, const SUMMM_1_1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),CT(x),smm.GetM2()); }


// (x*m+x*m)*x
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator*(
    const SUMMM_1_1<T,T1,T2 X3>& smm, const T x)
{ return SUMMM<T,T1,T2 X3>(x,smm.GetM1(),x,smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_1<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_1_1<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(x,smm.GetM1(),x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_1_1<T,T,T X3>& smm, const CCT x)
{ return SUMMM<CT,T,T X3>(CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_1_1<T,T,T X3>& smm, const VCT x)
{ return SUMMM<CT,T,T X3>(CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_1<CT,T1,T2 X3>& smm, const CCT x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_1<CT,T1,T2 X3>& smm, const VCT x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

// (x*m+x*m)/x
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator/(
    const SUMMM_1_1<T,T1,T2 X3>& smm, const T x)
{
  return SUMMM<T,T1,T2 X3>(RealType(T)(1)/x,smm.GetM1(),RealType(T)(1)/x,
      smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_1<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(CT(T(1)/x),smm.GetM1(),CT(T(1)/x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_1_1<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(T(1)/x,smm.GetM1(),T(1)/x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_1_1<T,T,T X3>& smm, const CCT x)
{ return SUMMM<CT,T,T X3>(T(1)/CT(x),smm.GetM1(),T(1)/CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_1_1<T,T,T X3>& smm, const VCT x)
{ return SUMMM<CT,T,T X3>(T(1)/CT(x),smm.GetM1(),T(1)/CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_1<CT,T1,T2 X3>& smm, const CCT x)
{ return SUMMM<CT,T1,T2 X3>(T(1)/CT(x),smm.GetM1(),T(1)/CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_1<CT,T1,T2 X3>& smm, const VCT x)
{ return SUMMM<CT,T1,T2 X3>(T(1)/CT(x),smm.GetM1(),T(1)/CT(x),smm.GetM2()); }


// SUMMM_1_m1

// -(x*m+x*m)
template <class T, class T1, class T2 Y> inline SUMMM_x_1<T,T1,T2 X3> operator-(
    const SUMMM_1_m1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(T(-1),smm.GetM1(),T(1),smm.GetM2); }

// x*(x*m+x*m)
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator*(
    const T x, const SUMMM_1_m1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(x,smm.GetM1(),-x,smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const T x, const SUMMM_1_m1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),CT(-x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const CT x, const SUMMM_1_m1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(x,smm.GetM1(),-x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const CCT x, const SUMMM_1_m1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(CT(x),smm.GetM1(),-CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const VCT x, const SUMMM_1_m1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(CT(x),smm.GetM1(),-CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const CCT x, const SUMMM_1_m1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),-CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const VCT x, const SUMMM_1_m1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),-CT(x),smm.GetM2()); }


// (x*m+x*m)*x
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator*(
    const SUMMM_1_m1<T,T1,T2 X3>& smm, const T x)
{ return SUMMM<T,T1,T2 X3>(x,smm.GetM1(),-x,smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_m1<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(x,smm.GetM1(),CT(-x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_1_m1<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(x,smm.GetM1(),-x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_1_m1<T,T,T X3>& smm, const CCT x)
{ return SUMMM<CT,T,T X3>(CT(x),smm.GetM1(),-CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_1_m1<T,T,T X3>& smm, const VCT x)
{ return SUMMM<CT,T,T X3>(CT(x),smm.GetM1(),-CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_m1<CT,T1,T2 X3>& smm, const CCT x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),-CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_m1<CT,T1,T2 X3>& smm, const VCT x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),-CT(x),smm.GetM2()); }

// (x*m+x*m)/x
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator/(
    const SUMMM_1_m1<T,T1,T2 X3>& smm, const T x)
{
  return SUMMM<T,T1,T2 X3>(RealType(T)(1)/x,smm.GetM1(),
      RealType(T)(-1)/x,smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_m1<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(T(1)/x,smm.GetM1(),T(-1)/x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_1_m1<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(T(1)/x,smm.GetM1(),T(-1)/x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_1_m1<T,T,T X3>& smm, const CCT x)
{ return SUMMM<CT,T,T X3>(T(1)/CT(x),smm.GetM1(),T(-1)/CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_1_m1<T,T,T X3>& smm, const VCT x)
{ return SUMMM<CT,T,T X3>(T(1)/CT(x),smm.GetM1(),T(-1)/CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_m1<CT,T1,T2 X3>& smm, const CCT x)
{ return SUMMM<CT,T1,T2 X3>(T(1)/CT(x),smm.GetM1(),T(-1)/CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_m1<CT,T1,T2 X3>& smm, const VCT x)
{ return SUMMM<CT,T1,T2 X3>(T(1)/CT(x),smm.GetM1(),T(-1)/CT(x),smm.GetM2()); }

// SUMMM_1_m1

// -(x*m+x*m)
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator-(
    const SUMMM_1_x<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(T(-1),smm.GetM1(),-smm.GetX2,smm.GetM2); }

// x*(x*m+x*m)
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator*(
    const T x, const SUMMM_1_x<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(x,smm.GetM1(),smm.GetX2()*x,smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const T x, const SUMMM_1_x<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),smm.GetX2()*x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const CT x, const SUMMM_1_x<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(x,smm.GetM1(),smm.GetX2()*x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const CCT x, const SUMMM_1_x<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(CT(x),smm.GetM1(),smm.GetX2()*CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const VCT x, const SUMMM_1_x<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(CT(x),smm.GetM1(),smm.GetX2()*CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const CCT x, const SUMMM_1_x<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),smm.GetX2()*CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const VCT x, const SUMMM_1_x<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),smm.GetX2()*CT(x),smm.GetM2()); }


// (x*m+x*m)*x
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator*(
    const SUMMM_1_x<T,T1,T2 X3>& smm, const T x)
{ return SUMMM<T,T1,T2 X3>(x,smm.GetM1(),smm.GetX2()*x,smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_x<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),smm.GetX2()*x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_1_x<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(x,smm.GetM1(),smm.GetX2()*x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_1_x<T,T,T X3>& smm, const CCT x)
{ return SUMMM<CT,T,T X3>(CT(x),smm.GetM1(),smm.GetX2()*CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_1_x<T,T,T X3>& smm, const VCT x)
{ return SUMMM<CT,T,T X3>(CT(x),smm.GetM1(),smm.GetX2()*CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_x<CT,T1,T2 X3>& smm, const CCT x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),smm.GetX2()*CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_x<CT,T1,T2 X3>& smm, const VCT x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm.GetM1(),smm.GetX2()*CT(x),smm.GetM2()); }

// (x*m+x*m)/x
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator/(
    const SUMMM_1_x<T,T1,T2 X3>& smm, const T x)
{
  return SUMMM<T,T1,T2 X3>(RealType(T)(1)/x,smm.GetM1(),smm.GetX2()/x,
      smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_x<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(CT(T(1)/x),smm.GetM1(),smm.GetX2()/x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_1_x<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(T(1)/x,smm.GetM1(),smm.GetX2()/x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_1_x<T,T,T X3>& smm, const CCT x)
{ 
  return SUMMM<CT,T,T X3>(T(1)/CT(x),smm.GetM1(),
      smm.GetX2()/CT(x),smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_1_x<T,T,T X3>& smm, const VCT x)
{ 
  return SUMMM<CT,T,T X3>(T(1)/CT(x),smm.GetM1(),
      smm.GetX2()/CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_x<CT,T1,T2 X3>& smm, const CCT x)
{ 
  return SUMMM<CT,T1,T2 X3>(T(1)/CT(x),smm.GetM1(),
      smm.GetX2()/CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_x<CT,T1,T2 X3>& smm, const VCT x)
{ 
  return SUMMM<CT,T1,T2 X3>(T(1)/CT(x),smm.GetM1(),
      smm.GetX2()/CT(x),smm.GetM2()); 
}

// SUMMM_1_m1

// -(x*m+x*m)
template <class T, class T1, class T2 Y> inline SUMMM_x_m1<T,T1,T2 X3> operator-(
    const SUMMM_x_1<T,T1,T2 X3>& smm)
{ return SUMMM_x_m1<T,T1,T2 X3>(-smm.GetX1(),smm.GetM1(),T(-1),smm.GetM2); }

// x*(x*m+x*m)
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator*(
    const T x, const SUMMM_x_1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(smm.GetX1()*x,smm.GetM1(),x,smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const T x, const SUMMM_x_1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(smm.GetX1()*x,smm.GetM1(),CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const CT x, const SUMMM_x_1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(smm.GetX1()*x,smm.GetM1(),x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const CCT x, const SUMMM_x_1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const VCT x, const SUMMM_x_1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const CCT x, const SUMMM_x_1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const VCT x, const SUMMM_x_1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }


// (x*m+x*m)*x
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator*(
    const SUMMM_x_1<T,T1,T2 X3>& smm, const T x)
{ return SUMMM<T,T1,T2 X3>(smm.GetX1()*x,smm.GetM1(),x,smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_x_1<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(smm.GetX1()*x,smm.GetM1(),CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_x_1<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(smm.GetX1()*x,smm.GetM1(),x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_x_1<T,T,T X3>& smm, const CCT x)
{ return SUMMM<CT,T,T X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_x_1<T,T,T X3>& smm, const VCT x)
{ return SUMMM<CT,T,T X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_x_1<CT,T1,T2 X3>& smm, const CCT x)
{ return SUMMM<CT,T1,T2 X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_x_1<CT,T1,T2 X3>& smm, const VCT x)
{ return SUMMM<CT,T1,T2 X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

// (x*m+x*m)/x
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator/(
    const SUMMM_x_1<T,T1,T2 X3>& smm, const T x)
{
  return SUMMM<T,T1,T2 X3>(smm.GetX1()/x,smm.GetM1(),RealType(T)(1)/x,
    smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_x_1<CT,T1,T2 X3>& smm, const T x)
{
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()/x,smm.GetM1(),CT(T(1)/x),
      smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_x_1<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(smm.GetX1()/x,smm.GetM1(),T(1)/x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_x_1<T,T,T X3>& smm, const CCT x)
{ 
  return SUMMM<CT,T,T X3>(smm.GetX1()/CT(x),smm.GetM1(),
      T(1)/CT(x),smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_x_1<T,T,T X3>& smm, const VCT x)
{ 
  return SUMMM<CT,T,T X3>(smm.GetX1()/CT(x),smm.GetM1(),
      T(1)/CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_x_1<CT,T1,T2 X3>& smm, const CCT x)
{ 
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()/CT(x),smm.GetM1(),
      T(1)/CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_x_1<CT,T1,T2 X3>& smm, const VCT x)
{ 
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()/CT(x),smm.GetM1(),
      T(1)/CT(x),smm.GetM2()); 
}

// SUMMM_1_m1

// -(x*m+x*m)
template <class T, class T1, class T2 Y> inline SUMMM_x_1<T,T1,T2 X3> operator-(
    const SUMMM_x_m1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(-smm.GetX1(),smm.GetM1(),T(1),smm.GetM2); }

// x*(x*m+x*m)
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator*(
    const T x, const SUMMM_x_m1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(smm.GetX1()*x,smm.GetM1(),x,smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const T x, const SUMMM_x_m1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(smm.GetX1()*x,smm.GetM1(),CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const CT x, const SUMMM_x_m1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(smm.GetX1()*x,smm.GetM1(),x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const CCT x, const SUMMM_x_m1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const VCT x, const SUMMM_x_m1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const CCT x, const SUMMM_x_m1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const VCT x, const SUMMM_x_m1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }


// (x*m+x*m)*x
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator*(
    const SUMMM_x_m1<T,T1,T2 X3>& smm, const T x)
{ return SUMMM<T,T1,T2 X3>(smm.GetX1()*x,smm.GetM1(),x,smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_x_m1<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(smm.GetX1()*x,smm.GetM1(),CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_x_m1<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(smm.GetX1()*x,smm.GetM1(),x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_x_m1<T,T,T X3>& smm, const CCT x)
{ return SUMMM<CT,T,T X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator*(
    const SUMMM_x_m1<T,T,T X3>& smm, const VCT x)
{ return SUMMM<CT,T,T X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_x_m1<CT,T1,T2 X3>& smm, const CCT x)
{ return SUMMM<CT,T1,T2 X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_x_m1<CT,T1,T2 X3>& smm, const VCT x)
{ return SUMMM<CT,T1,T2 X3>(smm.GetX1()*CT(x),smm.GetM1(),CT(x),smm.GetM2()); }

// (x*m+x*m)/x
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator/(
    const SUMMM_x_m1<T,T1,T2 X3>& smm, const T x)
{ 
  return SUMMM<T,T1,T2 X3>(smm.GetX1()/x,smm.GetM1(),
      RealType(T)(-1)/x,smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_x_m1<CT,T1,T2 X3>& smm, const T x)
{
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()/x,smm.GetM1(),
    CT(T(-1)/x),smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_x_m1<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(smm.GetX1()/x,smm.GetM1(),T(-1)/x,smm.GetM2()); }

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_x_m1<T,T,T X3>& smm, const CCT x)
{ 
  return SUMMM<CT,T,T X3>(smm.GetX1()/CT(x),smm.GetM1(),
      T(-1)/CT(x),smm.GetM2()); 
}

template <class T Y> inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_x_m1<T,T,T X3>& smm, const VCT x)
{ 
  return SUMMM<CT,T,T X3>(smm.GetX1()/CT(x),smm.GetM1(),
      T(-1)/CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_x_m1<CT,T1,T2 X3>& smm, const CCT x)
{ 
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()/CT(x),smm.GetM1(),
      T(-1)/CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2 Y> inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_x_m1<CT,T1,T2 X3>& smm, const VCT x)
{ 
  return SUMMM<CT,T1,T2 X3>(smm.GetX1()/CT(x),smm.GetM1(),
      T(-1)/CT(x),smm.GetM2()); 
}
#undef SUMMM_1_1
#undef SUMMM_1_m1
#undef SUMMM_1_x
#undef SUMMM_x_m1
#undef SUMMM_x_1
#endif

#undef X3
#undef Y

