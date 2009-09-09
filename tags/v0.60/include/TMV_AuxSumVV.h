///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// Things that need to be #defined on entry:
// (The values for a normal Vector+Vector are given)
//
// SUMVV	SumVV
// GENVECTOR1	GenVector
// GENVECTOR2	GenVector
// PRODXV1	ProdXV
// PRODXV2	ProdXV
//
// The following are null for regular Vectors.
// See TMV_SmallVectorArith.h for where these are important
// X1
// X2
// X3
// Y

#ifndef X1
#define X1
#endif

#ifndef X2
#define X2
#endif

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

// v+v
template <class T Y> inline SUMVV<T,T,T X3> operator+(
    const GENVECTOR1<T X1>& v1, const GENVECTOR2<T X2>& v2)
{ return SUMVV<T,T,T X3>(T(1),v1,T(1),v2); }

template <class T Y> inline SUMVV<CT,CT,T X3> operator+(
    const GENVECTOR1<CT X1>& v1,const GENVECTOR2<T X2>& v2)
{ return SUMVV<CT,CT,T X3>(CT(1),v1,CT(1),v2); }

template <class T Y> inline SUMVV<CT,T,CT X3> operator+(
    const GENVECTOR1<T X1>& v1,const GENVECTOR2<CT X2>& v2)
{ return SUMVV<CT,T,CT X3>(CT(1),v1,CT(1),v2); }

// v-v
template <class T Y> inline SUMVV<T,T,T X3> operator-(
    const GENVECTOR1<T X1>& v1, const GENVECTOR2<T X2>& v2)
{ return SUMVV<T,T,T X3>(T(1),v1,T(-1),v2); }

template <class T Y> inline SUMVV<CT,CT,T X3> operator-(
    const GENVECTOR1<CT X1>& v1,const GENVECTOR2<T X2>& v2)
{ return SUMVV<CT,CT,T X3>(CT(1),v1,CT(-1),v2); }

template <class T Y> inline SUMVV<CT,T,CT X3> operator-(
    const GENVECTOR1<T X1>& v1,const GENVECTOR2<CT X2>& v2)
{ return SUMVV<CT,T,CT X3>(CT(1),v1,CT(-1),v2); }

// -(x*v+x*v)
template <class T, class T1, class T2 Y> 
inline SUMVV<T,T1,T2 X3> operator-(
    const SUMVV<T,T1,T2 X3>& svv)
{ return SUMVV<T,T1,T2 X3>(-svv.GetX1(),svv.GetV1(),-svv.GetX2,svv.GetV2); }

// x*(x*v+x*v)
template <class T, class T1, class T2 Y> 
inline SUMVV<T,T1,T2 X3> operator*(
    const T x, const SUMVV<T,T1,T2 X3>& svv)
{
  return SUMVV<T,T1,T2 X3>(svv.GetX1()*x,svv.GetV1(),
      svv.GetX2()*x,svv.GetV2()); 
}

template <class T, class T1, class T2 Y> 
inline SUMVV<CT,T1,T2 X3> operator*(
    const T x, const SUMVV<CT,T1,T2 X3>& svv)
{
  return SUMVV<CT,T1,T2 X3>(svv.GetX1()*x,svv.GetV1(),
      svv.GetX2()*x,svv.GetV2()); 
}

template <class T Y> 
inline SUMVV<CT,T,T X3> operator*(
    const CT x, const SUMVV<T,T,T X3>& svv)
{ 
  return SUMVV<CT,T,T X3>(svv.GetX1()*x,svv.GetV1(),
      svv.GetX2()*x,svv.GetV2()); 
}

template <class T Y> inline SUMVV<CT,T,T X3> operator*(
    const CCT x, const SUMVV<T,T,T X3>& svv)
{ 
  return SUMVV<CT,T,T X3>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

template <class T Y> inline SUMVV<CT,T,T X3> operator*(
    const VCT x, const SUMVV<T,T,T X3>& svv)
{ 
  return SUMVV<CT,T,T X3>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

template <class T, class T1, class T2 Y> 
inline SUMVV<CT,T1,T2 X3> operator*(
    const CCT x, const SUMVV<CT,T1,T2 X3>& svv)
{ 
  return SUMVV<CT,T1,T2 X3>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

template <class T, class T1, class T2 Y> 
inline SUMVV<CT,T1,T2 X3> operator*(
    const VCT x, const SUMVV<CT,T1,T2 X3>& svv)
{ 
  return SUMVV<CT,T1,T2 X3>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

// (x*v+x*v)*x
template <class T, class T1, class T2 Y> 
inline SUMVV<T,T1,T2 X3> operator*(
    const SUMVV<T,T1,T2 X3>& svv, const T x)
{
  return SUMVV<T,T1,T2 X3>(svv.GetX1()*x,svv.GetV1(),
      svv.GetX2()*x,svv.GetV2()); 
}

template <class T, class T1, class T2 Y> 
inline SUMVV<CT,T1,T2 X3> operator*(
    const SUMVV<CT,T1,T2 X3>& svv, const T x)
{
  return SUMVV<CT,T1,T2 X3>(svv.GetX1()*x,svv.GetV1(),
      svv.GetX2()*x,svv.GetV2()); 
}

template <class T Y> inline SUMVV<CT,T,T X3> operator*(
    const SUMVV<T,T,T X3>& svv, const CT x)
{
  return SUMVV<CT,T,T X3>(svv.GetX1()*x,svv.GetV1(),
      svv.GetX2()*x,svv.GetV2()); 
}

template <class T Y> inline SUMVV<CT,T,T X3> operator*(
    const SUMVV<T,T,T X3>& svv, const CCT x)
{ 
  return SUMVV<CT,T,T X3>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

template <class T Y> inline SUMVV<CT,T,T X3> operator*(
    const SUMVV<T,T,T X3>& svv, const VCT x)
{ 
  return SUMVV<CT,T,T X3>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

template <class T, class T1, class T2 Y> 
inline SUMVV<CT,T1,T2 X3> operator*(
    const SUMVV<CT,T1,T2 X3>& svv, const CCT x)
{ 
  return SUMVV<CT,T1,T2 X3>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

template <class T, class T1, class T2 Y> 
inline SUMVV<CT,T1,T2 X3> operator*(
    const SUMVV<CT,T1,T2 X3>& svv, const VCT x)
{ 
  return SUMVV<CT,T1,T2 X3>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

// (x*v+x*v)/x
template <class T, class T1, class T2 Y> 
inline SUMVV<T,T1,T2 X3> operator/(
    const SUMVV<T,T1,T2 X3>& svv, const T x)
{ 
  return SUMVV<T,T1,T2 X3>(svv.GetX1()/x,svv.GetV1(),
      svv.GetX2()/x,svv.GetV2()); 
}

template <class T, class T1, class T2 Y> 
inline SUMVV<CT,T1,T2 X3> operator/(
    const SUMVV<CT,T1,T2 X3>& svv, const T x)
{
  return SUMVV<CT,T1,T2 X3>(svv.GetX1()/x,svv.GetV1(),
      svv.GetX2()/x,svv.GetV2()); 
}

template <class T Y> inline SUMVV<CT,T,T X3> operator/(
    const SUMVV<T,T,T X3>& svv, const CT x)
{ 
  return SUMVV<CT,T,T X3>(svv.GetX1()/x,svv.GetV1(),
      svv.GetX2()/x,svv.GetV2()); 
}

template <class T Y> inline SUMVV<CT,T,T X3> operator/(
    const SUMVV<T,T,T X3>& svv, const CCT x)
{ 
  return SUMVV<CT,T,T X3>(svv.GetX1()/CT(x),svv.GetV1(),
      svv.GetX2()/CT(x),svv.GetV2()); 
}

template <class T Y> inline SUMVV<CT,T,T X3> operator/(
    const SUMVV<T,T,T X3>& svv, const VCT x)
{ 
  return SUMVV<CT,T,T X3>(svv.GetX1()/CT(x),svv.GetV1(),
      svv.GetX2()/CT(x),svv.GetV2()); 
}

template <class T, class T1, class T2 Y> 
inline SUMVV<CT,T1,T2 X3> operator/(
    const SUMVV<CT,T1,T2 X3>& svv, const CCT x)
{ 
  return SUMVV<CT,T1,T2 X3>(svv.GetX1()/CT(x),svv.GetV1(),
      svv.GetX2()/CT(x),svv.GetV2()); 
}

template <class T, class T1, class T2 Y> 
inline SUMVV<CT,T1,T2 X3> operator/(
    const SUMVV<CT,T1,T2 X3>& svv, const VCT x)
{ 
  return SUMVV<CT,T1,T2 X3>(svv.GetX1()/CT(x),svv.GetV1(),
      svv.GetX2()/CT(x),svv.GetV2()); 
}

// (x*v)+v
template <class T Y> inline SUMVV<T,T,T X3> operator+(
    const PRODXV1<T,T X1>& pxv, const GENVECTOR2<T X2>& v)
{ return SUMVV<T,T,T X3>(pxv.GetX(),pxv.GetV(),T(1),v); }

template <class T Y> inline SUMVV<CT,T,CT X3> operator+(
    const PRODXV1<T,T X1>& pxv, const GENVECTOR2<CT X2>& v)
{ return SUMVV<CT,T,CT X3>(CT(pxv.GetX()),pxv.GetV(),CT(1),v); }

template <class T, class T1 Y> inline SUMVV<CT,T1,T X3> operator+(
    const PRODXV1<CT,T1 X1>& pxv, const GENVECTOR2<T X2>& v)
{ return SUMVV<CT,T1,T X3>(pxv.GetX(),pxv.GetV(),CT(1),v); }

// v+(x*v)
template <class T, class T2 Y> inline SUMVV<T,T,T2 X3> operator+(
    const GENVECTOR1<T X1>& v, const PRODXV2<T,T2 X2>& pxv)
{ return SUMVV<T,T,T2 X3>(T(1),v,pxv.GetX(),pxv.GetV()); }

template <class T Y> inline SUMVV<CT,CT,T X3> operator+(
    const GENVECTOR1<CT X1>& v, const PRODXV2<T,T X2>& pxv)
{ return SUMVV<CT,CT,T X3>(CT(1),v,CT(pxv.GetX()),pxv.GetV()); }

template <class T, class T2 Y> inline SUMVV<CT,T,T2 X3> operator+(
    const GENVECTOR1<T X1>& v, const PRODXV2<CT,T2 X2>& pxv)
{ return SUMVV<CT,T,T2 X3>(CT(1),v,pxv.GetX(),pxv.GetV()); }

// (x*v)-v
template <class T, class T1 Y> inline SUMVV<T,T1,T X3> operator-(
    const PRODXV1<T,T1 X1>& pxv, const GENVECTOR2<T X2>& v)
{ return SUMVV<T,T1,T X3>(pxv.GetX(),pxv.GetV(),T(-1),v); }

template <class T Y> inline SUMVV<CT,T,CT X3> operator-(
    const PRODXV1<T,T X1>& pxv, const GENVECTOR2<CT X2>& v)
{ return SUMVV<CT,T,CT X3>(CT(pxv.GetX()),pxv.GetV(),CT(-1),v); }

template <class T, class T1 Y> inline SUMVV<CT,T1,T X3> operator-(
    const PRODXV1<CT,T1 X1>& pxv, const GENVECTOR2<T X2>& v)
{ return SUMVV<CT,T1,T X3>(pxv.GetX(),pxv.GetV(),CT(-1),v); }

// v-(x*v)
template <class T, class T2 Y> inline SUMVV<T,T,T2 X3> operator-(
    const GENVECTOR1<T X1>& v, const PRODXV2<T,T2 X2>& pxv)
{ return SUMVV<T,T,T2 X3>(T(1),v,-pxv.GetX(),pxv.GetV()); }

template <class T Y> inline SUMVV<CT,CT,T X3> operator-(
    const GENVECTOR1<CT X1>& v, const PRODXV2<T,T X2>& pxv)
{ return SUMVV<CT,CT,T X3>(CT(1),v,CT(-pxv.GetX()),pxv.GetV()); }

template <class T, class T2 Y> inline SUMVV<CT,T,T2 X3> operator-(
    const GENVECTOR1<T X1>& v, const PRODXV2<CT,T2 X2>& pxv)
{ return SUMVV<CT,T,T2 X3>(CT(1),v,-pxv.GetX(),pxv.GetV()); }

// (x*v)+(x*v)
template <class T, class T1, class T2 Y> 
inline SUMVV<T,T1,T2 X3> operator+(
    const PRODXV1<T,T1 X1>& pxv1, const PRODXV2<T,T2 X2>& pxv2)
{ return SUMVV<T,T1,T2 X3>(pxv1.GetX(),pxv1.GetV(),pxv2.GetX(),pxv2.GetV()); }

template <class T, class T2 Y> inline SUMVV<CT,T,T2 X3> operator+(
    const PRODXV1<T,T X1>& pxv1, const PRODXV2<CT,T2 X2>& pxv2)
{ return SUMVV<CT,T,T2 X3>(CT(pxv1.GetX()),pxv1.GetV(),pxv2.GetX(),pxv2.GetV()); }

template <class T, class T1 Y> inline SUMVV<CT,T1,T X3> operator+(
    const PRODXV1<CT,T1 X1>& pxv1, const PRODXV2<T,T X2>& pxv2)
{ return SUMVV<CT,T1,T X3>(pxv1.GetX(),pxv1.GetV(),CT(pxv2.GetX()),pxv2.GetV()); }

// (x*v)-(x*v)
template <class T, class T1, class T2 Y> 
inline SUMVV<T,T1,T2 X3> operator-(
    const PRODXV1<T,T1 X1>& pxv1, const PRODXV2<T,T2 X2>& pxv2)
{
  return SUMVV<T,T1,T2 X3>(pxv1.GetX(),pxv1.GetV(),
      -pxv2.GetX(),pxv2.GetV()); 
}

template <class T, class T2 Y> inline SUMVV<CT,T,T2 X3> operator-(
    const PRODXV1<T,T X1>& pxv1, const PRODXV2<CT,T2 X2>& pxv2)
{
  return SUMVV<CT,T,T2 X3>(CT(pxv1.GetX()),pxv1.GetV(),
      -pxv2.GetX(),pxv2.GetV()); 
}

template <class T, class T1 Y> inline SUMVV<CT,T1,T X3> operator-(
    const PRODXV1<CT,T1 X1>& pxv1, const PRODXV2<T,T X2>& pxv2)
{
  return SUMVV<CT,T1,T X3>(pxv1.GetX(),pxv1.GetV(),
      CT(-pxv2.GetX()),pxv2.GetV()); 
}

#undef X1
#undef X2
#undef X3
#undef Y
