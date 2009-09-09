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


// Need to define the following with #define statements.
// (The given definition is for a regular Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX GenMatrix
// #define QUOTXM QuotXM

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

// -(x/m)
template <class T, class T2 Y> 
inline QUOTXM<T,T2 X> operator-(const QUOTXM<T,T2 X>& qxm)
{ return QUOTXM<T,T2 X>(-qxm.GetX(),qxm GETM); }

// x*(x/m)
template <class T, class T2 Y> 
inline QUOTXM<T,T2 X> operator*(const T x, const QUOTXM<T,T2 X>& qxm)
{ return QUOTXM<T,T2 X>(x*qxm.GetX(),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator*(const T x, const QUOTXM<CT,T2 X>& qxm)
{ return QUOTXM<CT,T2 X>(x*qxm.GetX(),qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator*(const CT x, const QUOTXM<T,T X>& qxm)
{ return QUOTXM<CT,T X>(x*qxm.GetX(),qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator*(const CCT x, const QUOTXM<T,T X>& qxm)
{ return QUOTXM<CT,T X>(CT(x)*qxm.GetX(),qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator*(const VCT x, const QUOTXM<T,T X>& qxm)
{ return QUOTXM<CT,T X>(CT(x)*qxm.GetX(),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator*(const CCT x, const QUOTXM<CT,T2 X>& qxm)
{ return QUOTXM<CT,T2 X>(CT(x)*qxm.GetX(),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator*(const VCT x, const QUOTXM<CT,T2 X>& qxm)
{ return QUOTXM<CT,T2 X>(CT(x)*qxm.GetX(),qxm GETM); }

// (x/m)*x
template <class T, class T2 Y> 
inline QUOTXM<T,T2 X> operator*(const QUOTXM<T,T2 X>& qxm, const T x)
{ return QUOTXM<T,T2 X>(x*qxm.GetX(),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator*(const QUOTXM<CT,T2 X>& qxm, const T x)
{ return QUOTXM<CT,T2 X>(x*qxm.GetX(),qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator*(const QUOTXM<T,T X>& qxm, const CT x)
{ return QUOTXM<CT,T X>(x*qxm.GetX(),qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator*(const QUOTXM<T,T X>& qxm, const CCT x)
{ return QUOTXM<CT,T X>(CT(x)*qxm.GetX(),qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator*(const QUOTXM<T,T X>& qxm, const VCT x)
{ return QUOTXM<CT,T X>(CT(x)*qxm.GetX(),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator*(const QUOTXM<CT,T2 X>& qxm, const CCT x)
{ return QUOTXM<CT,T2 X>(CT(x)*qxm.GetX(),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator*(const QUOTXM<CT,T2 X>& qxm, const VCT x)
{ return QUOTXM<CT,T2 X>(CT(x)*qxm.GetX(),qxm GETM); }

// (x/m)/x
template <class T, class T2 Y> 
inline QUOTXM<T,T2 X> operator/(const QUOTXM<T,T2 X>& qxm, const T x)
{ return QUOTXM<T,T2 X>(qxm.GetX()/x,qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator/(const QUOTXM<CT,T2 X>& qxm, const T x)
{ return QUOTXM<CT,T2 X>(qxm.GetX()/x,qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator/(const QUOTXM<T,T X>& qxm, const CT x)
{ return QUOTXM<CT,T X>(qxm.GetX()/x,qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator/(const QUOTXM<T,T X>& qxm, const CCT x)
{ return QUOTXM<CT,T X>(qxm.GetX()/CT(x),qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator/(const QUOTXM<T,T X>& qxm, const VCT x)
{ return QUOTXM<CT,T X>(qxm.GetX()/CT(x),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator/(const QUOTXM<CT,T2 X>& qxm, const CCT x)
{ return QUOTXM<CT,T2 X>(qxm.GetX()/CT(x),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator/(const QUOTXM<CT,T2 X>& qxm, const VCT x)
{ return QUOTXM<CT,T2 X>(qxm.GetX()/CT(x),qxm GETM); }

#ifdef QUOTXM_1

// -(1/m)
template <class T, class T2 Y> 
inline QUOTXM<T,T2 X> operator-(const QUOTXM_1<T,T2 X>& qxm)
{ return QUOTXM<T,T2 X>(T(-1),qxm GETM); }

// x*(1/m)
template <class T, class T2 Y> 
inline QUOTXM<T,T2 X> operator*(const T x, const QUOTXM_1<T,T2 X>& qxm)
{ return QUOTXM<T,T2 X>(x,qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator*(const T x, const QUOTXM_1<CT,T2 X>& qxm)
{ return QUOTXM<CT,T2 X>(x,qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator*(const CT x, const QUOTXM_1<T,T X>& qxm)
{ return QUOTXM<CT,T X>(x,qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator*(const CCT x, const QUOTXM_1<T,T X>& qxm)
{ return QUOTXM<CT,T X>(CT(x),qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator*(const VCT x, const QUOTXM_1<T,T X>& qxm)
{ return QUOTXM<CT,T X>(CT(x),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator*(const CCT x, const QUOTXM_1<CT,T2 X>& qxm)
{ return QUOTXM<CT,T2 X>(CT(x),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator*(const VCT x, const QUOTXM_1<CT,T2 X>& qxm)
{ return QUOTXM<CT,T2 X>(CT(x),qxm GETM); }

// (1/m)*x
template <class T, class T2 Y> 
inline QUOTXM<T,T2 X> operator*(const QUOTXM_1<T,T2 X>& qxm, const T x)
{ return QUOTXM<T,T2 X>(x,qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator*(const QUOTXM_1<CT,T2 X>& qxm, const T x)
{ return QUOTXM<CT,T2 X>(x,qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator*(const QUOTXM_1<T,T X>& qxm, const CT x)
{ return QUOTXM<CT,T X>(x,qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator*(const QUOTXM_1<T,T X>& qxm, const CCT x)
{ return QUOTXM<CT,T X>(CT(x),qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator*(const QUOTXM_1<T,T X>& qxm, const VCT x)
{ return QUOTXM<CT,T X>(CT(x),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator*(const QUOTXM_1<CT,T2 X>& qxm, const CCT x)
{ return QUOTXM<CT,T2 X>(CT(x),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator*(const QUOTXM_1<CT,T2 X>& qxm, const VCT x)
{ return QUOTXM<CT,T2 X>(CT(x),qxm GETM); }

// (1/m)/x
template <class T, class T2 Y> 
inline QUOTXM<T,T2 X> operator/(const QUOTXM_1<T,T2 X>& qxm, const T x)
{ return QUOTXM<T,T2 X>(RealType(T)(1)/x,qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator/(const QUOTXM_1<CT,T2 X>& qxm, const T x)
{ return QUOTXM<CT,T2 X>(T(1)/x,qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator/(const QUOTXM_1<T,T X>& qxm, const CT x)
{ return QUOTXM<CT,T X>(T(1)/x,qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator/(const QUOTXM_1<T,T X>& qxm, const CCT x)
{ return QUOTXM<CT,T X>(T(1)/CT(x),qxm GETM); }

template <class T Y> 
inline QUOTXM<CT,T X> operator/(const QUOTXM_1<T,T X>& qxm, const VCT x)
{ return QUOTXM<CT,T X>(T(1)/CT(x),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator/(const QUOTXM_1<CT,T2 X>& qxm, const CCT x)
{ return QUOTXM<CT,T2 X>(T(1)/CT(x),qxm GETM); }

template <class T, class T2 Y> 
inline QUOTXM<CT,T2 X> operator/(const QUOTXM_1<CT,T2 X>& qxm, const VCT x)
{ return QUOTXM<CT,T2 X>(T(1)/CT(x),qxm GETM); }

#undef QUOTXM_1
#endif

#undef X
#undef Y
#undef GETM
