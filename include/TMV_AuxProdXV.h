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
// (The values for a normal Vector are given)
//
// PRODXV	ProdXV
// GENVECTOR	GenVector

#ifndef X
#define X
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

// -v
template <class T Y> inline PRODXV<T,T X> operator-(
    const GENVECTOR<T X>& v)
{ return PRODXV<T,T X>(T(-1),v); }

// x*v
template <class T Y> inline PRODXV<T,T X> operator*(
    T x, const GENVECTOR<T X>& v) 
{ return PRODXV<T,T X>(x,v); }

template <class T Y> inline PRODXV<CT,T X> operator*(
    CT x, const GENVECTOR<T X>& v)
{ return PRODXV<CT,T X>(x,v); }

template <class T Y> inline PRODXV<CT,T X> operator*(
    CCT x, const GENVECTOR<T X>& v)
{ return PRODXV<CT,T X>(CT(x),v); }

template <class T Y> inline PRODXV<CT,T X> operator*(
    VCT x, const GENVECTOR<T X>& v)
{ return PRODXV<CT,T X>(CT(x),v); }

template <class T Y> inline PRODXV<CT,CT X> operator*(
    T x, const GENVECTOR<CT X>& v) 
{ return PRODXV<CT,CT X>(x,v); }

template <class T Y> inline PRODXV<CT,CT X> operator*(
    CCT x, const GENVECTOR<CT X>& v)
{ return PRODXV<CT,CT X>(CT(x),v); }

template <class T Y> inline PRODXV<CT,CT X> operator*(
    VCT x, const GENVECTOR<CT X>& v)
{ return PRODXV<CT,CT X>(CT(x),v); }

// v*x
template <class T Y> inline PRODXV<T,T X> operator*(
    const GENVECTOR<T X>& v, T x) 
{ return PRODXV<T,T X>(x,v); }

template <class T Y> inline PRODXV<CT,T X> operator*(
    const GENVECTOR<T X>& v, CT x)
{ return PRODXV<CT,T X>(x,v); }

template <class T Y> inline PRODXV<CT,T X> operator*(
    const GENVECTOR<T X>& v, CCT x)
{ return PRODXV<CT,T X>(CT(x),v); }

template <class T Y> inline PRODXV<CT,T X> operator*(
    const GENVECTOR<T X>& v, VCT x)
{ return PRODXV<CT,T X>(CT(x),v); }

template <class T Y> inline PRODXV<CT,CT X> operator*(
    const GENVECTOR<CT X>& v, T x) 
{ return PRODXV<CT,CT X>(x,v); }

template <class T Y> inline PRODXV<CT,CT X> operator*(
    const GENVECTOR<CT X>& v, CCT x)
{ return PRODXV<CT,CT X>(CT(x),v); }

template <class T Y> inline PRODXV<CT,CT X> operator*(
    const GENVECTOR<CT X>& v, VCT x)
{ return PRODXV<CT,CT X>(CT(x),v); }

// v/x
template <class T Y> inline PRODXV<T,T X> operator/(
    const GENVECTOR<T X>& v, T x) 
{ return PRODXV<T,T X>(RealType(T)(1)/x,v); }

template <class T Y> inline PRODXV<CT,T X> operator/(
    const GENVECTOR<T X>& v, CT x)
{ return PRODXV<CT,T X>(T(1)/x,v); }

template <class T Y> inline PRODXV<CT,T X> operator/(
    const GENVECTOR<T X>& v, CCT x)
{ return PRODXV<CT,T X>(T(1)/CT(x),v); }

template <class T Y> inline PRODXV<CT,T X> operator/(
    const GENVECTOR<T X>& v, VCT x)
{ return PRODXV<CT,T X>(T(1)/CT(x),v); }

template <class T Y> inline PRODXV<CT,CT X> operator/(
    const GENVECTOR<CT X>& v, T x)
{ return PRODXV<CT,CT X>(T(1)/x,v); }

template <class T Y> inline PRODXV<CT,CT X> operator/(
    const GENVECTOR<CT X>& v, CCT x)
{ return PRODXV<CT,CT X>(T(1)/CT(x),v); }

template <class T Y> inline PRODXV<CT,CT X> operator/(
    const GENVECTOR<CT X>& v, VCT x)
{ return PRODXV<CT,CT X>(T(1)/CT(x),v); }

// -(x*v)
template <class T, class T2 Y> inline PRODXV<T,T2 X> operator-(
    const PRODXV<T,T2 X>& pxv)
{ return PRODXV<T,T2 X>(-pxv.GetX(),pxv.GetV()); }

// x*(x*v)
template <class T, class T2 Y> inline PRODXV<T,T2 X> operator*(
    const T x, const PRODXV<T,T2 X>& pxv)
{ return PRODXV<T,T2 X>(x*pxv.GetX(),pxv.GetV()); }

template <class T, class T2 Y> inline PRODXV<CT,T2 X> operator*(
    const T x, const PRODXV<CT,T2 X>& pxv)
{ return PRODXV<CT,T2 X>(x*pxv.GetX(),pxv.GetV()); }

template <class T Y> inline PRODXV<CT,T X> operator*(
    const CT x, const PRODXV<T,T X>& pxv)
{ return PRODXV<CT,T X>(x*pxv.GetX(),pxv.GetV()); }

template <class T Y> inline PRODXV<CT,T X> operator*(
    const CCT x, const PRODXV<T,T X>& pxv)
{ return PRODXV<CT,T X>(CT(x)*pxv.GetX(),pxv.GetV()); }

template <class T Y> inline PRODXV<CT,T X> operator*(
    const VCT x, const PRODXV<T,T X>& pxv)
{ return PRODXV<CT,T X>(CT(x)*pxv.GetX(),pxv.GetV()); }

template <class T, class T2 Y> inline PRODXV<CT,T2 X> operator*(
    const CCT x, const PRODXV<CT,T2 X>& pxv)
{ return PRODXV<CT,T2 X>(CT(x)*pxv.GetX(),pxv.GetV()); }

template <class T, class T2 Y> inline PRODXV<CT,T2 X> operator*(
    const VCT x, const PRODXV<CT,T2 X>& pxv)
{ return PRODXV<CT,T2 X>(CT(x)*pxv.GetX(),pxv.GetV()); }

// (x*v)*x
template <class T, class T2 Y> inline PRODXV<T,T2 X> operator*(
    const PRODXV<T,T2 X>& pxv, const T x)
{ return PRODXV<T,T2 X>(x*pxv.GetX(),pxv.GetV()); }

template <class T, class T2 Y> inline PRODXV<CT,T2 X> operator*(
    const PRODXV<CT,T2 X>& pxv, const T x)
{ return PRODXV<CT,T2 X>(x*pxv.GetX(),pxv.GetV()); }

template <class T Y> inline PRODXV<CT,T X> operator*(
    const PRODXV<T,T X>& pxv, const CT x)
{ return PRODXV<CT,T X>(x*pxv.GetX(),pxv.GetV()); }

template <class T Y> inline PRODXV<CT,T X> operator*(
    const PRODXV<T,T X>& pxv, const CCT x)
{ return PRODXV<CT,T X>(CT(x)*pxv.GetX(),pxv.GetV()); }

template <class T Y> inline PRODXV<CT,T X> operator*(
    const PRODXV<T,T X>& pxv, const VCT x)
{ return PRODXV<CT,T X>(CT(x)*pxv.GetX(),pxv.GetV()); }

template <class T, class T2 Y> inline PRODXV<CT,T2 X> operator*(
    const PRODXV<CT,T2 X>& pxv, const CCT x)
{ return PRODXV<CT,T2 X>(CT(x)*pxv.GetX(),pxv.GetV()); }

template <class T, class T2 Y> inline PRODXV<CT,T2 X> operator*(
    const PRODXV<CT,T2 X>& pxv, const VCT x)
{ return PRODXV<CT,T2 X>(CT(x)*pxv.GetX(),pxv.GetV()); }

// (x*v)/x
template <class T, class T2 Y> inline PRODXV<T,T X> operator/(
    const PRODXV<T,T2 X>& pxv, const T x)
{ return PRODXV<T,T2 X>(pxv.GetX()/x,pxv.GetV()); }

template <class T, class T2 Y> inline PRODXV<CT,T2 X> operator/(
    const PRODXV<CT,T2 X>& pxv, const T x)
{ return PRODXV<CT,T2 X>(pxv.GetX()/x,pxv.GetV()); }

template <class T Y> inline PRODXV<CT,T X> operator/(
    const PRODXV<T,T X>& pxv, const CT x)
{ return PRODXV<CT,T X>(pxv.GetX()/x,pxv.GetV()); }

template <class T Y> inline PRODXV<CT,T X> operator/(
    const PRODXV<T,T X>& pxv, const CCT x)
{ return PRODXV<CT,T X>(pxv.GetX()/CT(x),pxv.GetV()); }

template <class T Y> inline PRODXV<CT,T X> operator/(
    const PRODXV<T,T X>& pxv, const VCT x)
{ return PRODXV<CT,T X>(pxv.GetX()/CT(x),pxv.GetV()); }

template <class T, class T2 Y> inline PRODXV<CT,T2 X> operator/(
    const PRODXV<CT,T2 X>& pxv, const CCT x)
{ return PRODXV<CT,T2 X>(pxv.GetX()/CT(x),pxv.GetV()); }

template <class T, class T2 Y> inline PRODXV<CT,T2 X> operator/(
    const PRODXV<CT,T2 X>& pxv, const VCT x)
{ return PRODXV<CT,T2 X>(pxv.GetX()/CT(x),pxv.GetV()); }

#undef X
#undef Y
