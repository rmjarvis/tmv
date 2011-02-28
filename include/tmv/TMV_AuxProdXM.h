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
#define GETM .getM()
#endif

// -m

template <class T Y> 
inline PRODXM<T,T X> operator-(const GENMATRIX<T X>& m)
{ return PRODXM<T,T X>(T(-1),m); }

// x*m
#ifdef INTT
inline PRODXM<int,int> operator*(int x, const GENMATRIX<int>& m) 
{ return PRODXM<int,int>(x,m); }

inline PRODXM<float,int> operator*(float x, const GENMATRIX<int>& m) 
{ return PRODXM<float,int>(x,m); }

inline PRODXM<double,int> operator*(double x, const GENMATRIX<int>& m) 
{ return PRODXM<double,int>(x,m); }

inline PRODXM<long double,int> operator*(
    long double x, const GENMATRIX<int>& m) 
{ return PRODXM<long double,int>(x,m); }

template <class T Y> 
inline PRODXM<CT,int X> operator*(CT x, const GENMATRIX<int X>& m)
{ return PRODXM<CT,int X>(x,m); }

template <class T Y> 
inline PRODXM<CT,int X> operator*(CCT x, const GENMATRIX<int X>& m)
{ return PRODXM<CT,int X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,int X> operator*(VCT x, const GENMATRIX<int X>& m)
{ return PRODXM<CT,int X>(CT(x),m); }
#else
template <class T Y> 
inline PRODXM<T,T X> operator*(T x, const GENMATRIX<T X>& m) 
{ return PRODXM<T,T X>(x,m); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(CT x, const GENMATRIX<T X>& m)
{ return PRODXM<CT,T X>(x,m); }

template <class T Y> 
inline PRODXM<CT,CT X> operator*(T x, const GENMATRIX<CT X>& m) 
{ return PRODXM<CT,CT X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(CCT x, const GENMATRIX<T X>& m)
{ return PRODXM<CT,T X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(VCT x, const GENMATRIX<T X>& m)
{ return PRODXM<CT,T X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,CT X> operator*(CCT x, const GENMATRIX<CT X>& m)
{ return PRODXM<CT,CT X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,CT X> operator*(VCT x, const GENMATRIX<CT X>& m)
{ return PRODXM<CT,CT X>(CT(x),m); }
#endif

// m*x
#ifdef INTT
inline PRODXM<int,int> operator*(const GENMATRIX<int>& m, int x) 
{ return PRODXM<int,int>(x,m); }

inline PRODXM<float,int> operator*(const GENMATRIX<int>& m, float x) 
{ return PRODXM<float,int>(x,m); }

inline PRODXM<double,int> operator*(const GENMATRIX<int>& m, double x) 
{ return PRODXM<double,int>(x,m); }

inline PRODXM<long double,int> operator*(
    const GENMATRIX<int>& m, long double x) 
{ return PRODXM<long double,int>(x,m); }

template <class T Y> 
inline PRODXM<CT,int X> operator*(const GENMATRIX<int X>& m, CT x)
{ return PRODXM<CT,int X>(x,m); }

template <class T Y> 
inline PRODXM<CT,int X> operator*(const GENMATRIX<int X>& m, CCT x)
{ return PRODXM<CT,int X>(CT(x),m); }

template <class T Y> 
inline PRODXM<CT,int X> operator*(const GENMATRIX<int X>& m, VCT x)
{ return PRODXM<CT,int X>(CT(x),m); }
#else
template <class T Y> 
inline PRODXM<T,T X> operator*(const GENMATRIX<T X>& m, T x) 
{ return PRODXM<T,T X>(x,m); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const GENMATRIX<T X>& m, CT x)
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
#endif

// m/x
#ifdef INTT
inline PRODXM<int,int> operator/(const GENMATRIX<int>& m, int x) 
{ return PRODXM<int,int>(int(1)/int(x),m); }

inline PRODXM<float,int> operator/(const GENMATRIX<int>& m, float x) 
{ return PRODXM<float,int>(float(1)/float(x),m); }

inline PRODXM<double,int> operator/(const GENMATRIX<int>& m, double x) 
{ return PRODXM<double,int>(double(1)/double(x),m); }

inline PRODXM<long double,int> operator/(
    const GENMATRIX<int>& m, long double x) 
{ return PRODXM<long double,int>((long double)(1)/(long double)(x),m); }

template <class T Y> 
inline PRODXM<CT,int X> operator/(const GENMATRIX<int X>& m, CT x)
{ return PRODXM<CT,int X>(T(1)/x,m); }

template <class T Y> 
inline PRODXM<CT,T X> operator/(const GENMATRIX<int X>& m, CCT x)
{ return PRODXM<CT,T X>(T(1)/CT(x),m); }

template <class T Y> 
inline PRODXM<CT,int X> operator/(const GENMATRIX<int X>& m, VCT x)
{ return PRODXM<CT,int X>(T(1)/CT(x),m); }
#undef INTT
#else
template <class T Y> 
inline PRODXM<T,T X> operator/(const GENMATRIX<T X>& m, T x) 
{ return PRODXM<T,T X>(TMV_RealType(T)(1)/x,m); }

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
#endif

// -(x*m)

template <class T, class T2 Y> 
inline PRODXM<T,T2 X> operator-(const PRODXM<T,T2 X>& pxm)
{ return PRODXM<T,T2 X>(-pxm.getX(),pxm GETM); }

// x*(x*m)

template <class T, class T2 Y> 
inline PRODXM<T,T2 X> operator*(const T x, const PRODXM<T,T2 X>& pxm)
{ return PRODXM<T,T2 X>(x*pxm.getX(),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator*(const T x, const PRODXM<CT,T2 X>& pxm)
{ return PRODXM<CT,T2 X>(x*pxm.getX(),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const CT x, const PRODXM<T,T X>& pxm)
{ return PRODXM<CT,T X>(x*pxm.getX(),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const CCT x, const PRODXM<T,T X>& pxm)
{ return PRODXM<CT,T X>(CT(x)*pxm.getX(),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const VCT x, const PRODXM<T,T X>& pxm)
{ return PRODXM<CT,T X>(CT(x)*pxm.getX(),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator*(const CCT x, const PRODXM<CT,T2 X>& pxm)
{ return PRODXM<CT,T2 X>(CT(x)*pxm.getX(),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator*(const VCT x, const PRODXM<CT,T2 X>& pxm)
{ return PRODXM<CT,T2 X>(CT(x)*pxm.getX(),pxm GETM); }

// (x*m)*x

template <class T, class T2 Y> 
inline PRODXM<T,T2 X> operator*(const PRODXM<T,T2 X>& pxm, const T x)
{ return PRODXM<T,T2 X>(x*pxm.getX(),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator*(const PRODXM<CT,T2 X>& pxm, const T x)
{ return PRODXM<CT,T2 X>(x*pxm.getX(),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const PRODXM<T,T X>& pxm, const CT x)
{ return PRODXM<CT,T X>(x*pxm.getX(),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const PRODXM<T,T X>& pxm, const CCT x)
{ return PRODXM<CT,T X>(CT(x)*pxm.getX(),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator*(const PRODXM<T,T X>& pxm, const VCT x)
{ return PRODXM<CT,T X>(CT(x)*pxm.getX(),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator*(const PRODXM<CT,T2 X>& pxm, const CCT x)
{ return PRODXM<CT,T2 X>(CT(x)*pxm.getX(),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator*(const PRODXM<CT,T2 X>& pxm, const VCT x)
{ return PRODXM<CT,T2 X>(CT(x)*pxm.getX(),pxm GETM); }

// (x*m)/x

template <class T, class T2 Y> 
inline PRODXM<T,T2 X> operator/(const PRODXM<T,T2 X>& pxm, const T x)
{ return PRODXM<T,T2 X>(pxm.getX()/x,pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator/(const PRODXM<CT,T2 X>& pxm, const T x)
{ return PRODXM<CT,T2 X>(pxm.getX()/x,pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator/(const PRODXM<T,T X>& pxm, const CT x)
{ return PRODXM<CT,T X>(pxm.getX()/x,pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator/(const PRODXM<T,T X>& pxm, const CCT x)
{ return PRODXM<CT,T X>(pxm.getX()/CT(x),pxm GETM); }

template <class T Y> 
inline PRODXM<CT,T X> operator/(const PRODXM<T,T X>& pxm, const VCT x)
{ return PRODXM<CT,T X>(pxm.getX()/CT(x),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator/(const PRODXM<CT,T2 X>& pxm, const CCT x)
{ return PRODXM<CT,T2 X>(pxm.getX()/CT(x),pxm GETM); }

template <class T, class T2 Y> 
inline PRODXM<CT,T2 X> operator/(const PRODXM<CT,T2 X>& pxm, const VCT x)
{ return PRODXM<CT,T2 X>(pxm.getX()/CT(x),pxm GETM); }

#undef X
#undef Y
#undef GETM
