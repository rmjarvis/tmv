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


// Things that need to be #defined on entry:
// (The values for a normal Matrix are given)
//
// PRODXM	ProdXM
// GENMATRIX	GenMatrix

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

#ifndef GETM
#define GETM .GetM()
#endif

// -m
template <class T Y> 
inline PRODXM<T,T X> operator-(const GENMATRIX<T X>& m)
{ return PRODXM<T,T X>(T(-1),m); }

// x*m
template <class T Y> 
inline PRODXM<T,T X> operator*(T x, const GENMATRIX<T X>& m) 
{ return PRODXM<T,T X>(x,m); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(CT x, const GENMATRIX<T X>& m)
{ return PRODXM<CT,T X>(x,m); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(CCT x, const GENMATRIX<T X>& m)
{ return PRODXM<CT,T X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(VCT x, const GENMATRIX<T X>& m)
{ return PRODXM<CT,T X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,CT X> operator*(T x, const GENMATRIX<CT X>& m) 
{ return PRODXM<CT,CT X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,CT X> operator*(CCT x, const GENMATRIX<CT X>& m)
{ return PRODXM<CT,CT X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,CT X> operator*(VCT x, const GENMATRIX<CT X>& m)
{ return PRODXM<CT,CT X>(CT(x),m); }

// m*x
template <class T Y> 
inline PRODXM<T,T X> operator*(const GENMATRIX<T X>& m, T x) 
{ return PRODXM<T,T X>(x,m); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const GENMATRIX<T X>& m,CT x)
{ return PRODXM<CT,T X>(x,m); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const GENMATRIX<T X>& m, CCT x)
{ return PRODXM<CT,T X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const GENMATRIX<T X>& m, VCT x)
{ return PRODXM<CT,T X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,CT X> operator*(const GENMATRIX<CT X>& m, CCT x)
{ return PRODXM<CT,CT X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,CT X> operator*(const GENMATRIX<CT X>& m, VCT x)
{ return PRODXM<CT,CT X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,CT X> operator*(const GENMATRIX<CT X>& m, T x) 
{ return PRODXM<CT,CT X>(CT(x),m); }

// m/x
template <class T Y> 
inline PRODXM<T,T X> operator/(const GENMATRIX<T X>& m, T x) 
{ return PRODXM<T,T X>(RealType(T)(1)/x,m); }

template <class T Y> 
inline PRODXM<CT,T X> operator/(const GENMATRIX<T X>& m, CT x)
{ return PRODXM<CT,T X>(T(1)/x,m); }

template <class T Y> 
inline PRODXM<CT,T X> operator/(const GENMATRIX<T X>& m, CCT x)
{ return PRODXM<CT,T X>(T(1)/CT(x),m); }

template <class T Y> 
inline PRODXM<CT,T X> operator/(const GENMATRIX<T X>& m, VCT x)
{ return PRODXM<CT,T X>(T(1)/CT(x),m); }

template <class T Y> 
inline PRODXM<CT,CT X> operator/(const GENMATRIX<CT X>& m, CCT x)
{ return PRODXM<CT,CT X>(T(1)/CT(x),m); }

template <class T Y> 
inline PRODXM<CT,CT X> operator/(const GENMATRIX<CT X>& m, VCT x)
{ return PRODXM<CT,CT X>(T(1)/CT(x),m); }

template <class T Y> 
inline PRODXM<CT,CT X> operator/(const GENMATRIX<CT X>& m, T x)
{ return PRODXM<CT,CT X>(CT(T(1)/x),m); }

// -(x*m)
template <class T, class T2 Y> 
inline PRODXM<T,T2 X> operator-(const PRODXM<T,T2 X>& pxm)
{ return PRODXM<T,T2 X>(-pxm.GetX(),pxm GETM); }

// x*(x*m)
template <class T, class T2 Y> 
inline PRODXM<T,T2 X> operator*(const T x, const PRODXM<T,T2 X>& pxm)
{ return PRODXM<T,T2 X>(x*pxm.GetX(),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator*(const T x, const PRODXM<CT,T2 X>& pxm)
{ return PRODXM<CT,T2 X>(x*pxm.GetX(),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const CT x, const PRODXM<T,T X>& pxm)
{ return PRODXM<CT,T X>(x*pxm.GetX(),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const CCT x, const PRODXM<T,T X>& pxm)
{ return PRODXM<CT,T X>(CT(x)*pxm.GetX(),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const VCT x, const PRODXM<T,T X>& pxm)
{ return PRODXM<CT,T X>(CT(x)*pxm.GetX(),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator*(const CCT x, const PRODXM<CT,T2 X>& pxm)
{ return PRODXM<CT,T2 X>(CT(x)*pxm.GetX(),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator*(const VCT x, const PRODXM<CT,T2 X>& pxm)
{ return PRODXM<CT,T2 X>(CT(x)*pxm.GetX(),pxm GETM); }

// (x*m)*x
template <class T, class T2 Y> 
inline PRODXM<T,T2 X> operator*(const PRODXM<T,T2 X>& pxm, const T x)
{ return PRODXM<T,T2 X>(x*pxm.GetX(),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator*(const PRODXM<CT,T2 X>& pxm, const T x)
{ return PRODXM<CT,T2 X>(x*pxm.GetX(),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const PRODXM<T,T X>& pxm, const CT x)
{ return PRODXM<CT,T X>(x*pxm.GetX(),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const PRODXM<T,T X>& pxm, const CCT x)
{ return PRODXM<CT,T X>(CT(x)*pxm.GetX(),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const PRODXM<T,T X>& pxm, const VCT x)
{ return PRODXM<CT,T X>(CT(x)*pxm.GetX(),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator*(const PRODXM<CT,T2 X>& pxm, const CCT x)
{ return PRODXM<CT,T2 X>(CT(x)*pxm.GetX(),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator*(const PRODXM<CT,T2 X>& pxm, const VCT x)
{ return PRODXM<CT,T2 X>(CT(x)*pxm.GetX(),pxm GETM); }

// (x*m)/x
template <class T, class T2 Y> 
inline PRODXM<T,T2 X> operator/(const PRODXM<T,T2 X>& pxm, const T x)
{ return PRODXM<T,T2 X>(pxm.GetX()/x,pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator/(const PRODXM<CT,T2 X>& pxm, const T x)
{ return PRODXM<CT,T2 X>(pxm.GetX()/x,pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator/(const PRODXM<T,T X>& pxm, const CT x)
{ return PRODXM<CT,T X>(pxm.GetX()/x,pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator/(const PRODXM<T,T X>& pxm, const CCT x)
{ return PRODXM<CT,T X>(pxm.GetX()/CT(x),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator/(const PRODXM<T,T X>& pxm, const VCT x)
{ return PRODXM<CT,T X>(pxm.GetX()/CT(x),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator/(const PRODXM<CT,T2 X>& pxm, const CCT x)
{ return PRODXM<CT,T2 X>(pxm.GetX()/CT(x),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator/(const PRODXM<CT,T2 X>& pxm, const VCT x)
{ return PRODXM<CT,T2 X>(pxm.GetX()/CT(x),pxm GETM); }

#undef X
#undef Y
#undef GETM
