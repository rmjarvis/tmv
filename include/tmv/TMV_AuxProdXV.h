///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
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

template <typename T Y>
inline PRODXV<T,T X> operator-(const GENVECTOR<T X>& v)
{ return PRODXV<T,T X>(T(-1),v); }

// x*v

template <typename T Y>
inline PRODXV<T,T X> operator*(T x, const GENVECTOR<T X>& v)
{ return PRODXV<T,T X>(x,v); }

template <typename T Y>
inline PRODXV<CT,T X> operator*(CT x, const GENVECTOR<T X>& v)
{ return PRODXV<CT,T X>(x,v); }

template <typename T Y>
inline PRODXV<CT,T X> operator*(CCT x, const GENVECTOR<T X>& v)
{ return PRODXV<CT,T X>(CT(x),v); }

template <typename T Y>
inline PRODXV<CT,T X> operator*(VCT x, const GENVECTOR<T X>& v)
{ return PRODXV<CT,T X>(CT(x),v); }

template <typename T Y>
inline PRODXV<CT,CT X> operator*(T x, const GENVECTOR<CT X>& v)
{ return PRODXV<CT,CT X>(x,v); }

template <typename T Y>
inline PRODXV<CT,CT X> operator*(CCT x, const GENVECTOR<CT X>& v)
{ return PRODXV<CT,CT X>(CT(x),v); }

template <typename T Y>
inline PRODXV<CT,CT X> operator*(VCT x, const GENVECTOR<CT X>& v)
{ return PRODXV<CT,CT X>(CT(x),v); }

// v*x

template <typename T Y>
inline PRODXV<T,T X> operator*(const GENVECTOR<T X>& v, T x)
{ return PRODXV<T,T X>(x,v); }

template <typename T Y>
inline PRODXV<CT,T X> operator*(const GENVECTOR<T X>& v, CT x)
{ return PRODXV<CT,T X>(x,v); }

template <typename T Y>
inline PRODXV<CT,T X> operator*(const GENVECTOR<T X>& v, CCT x)
{ return PRODXV<CT,T X>(CT(x),v); }

template <typename T Y>
inline PRODXV<CT,T X> operator*(const GENVECTOR<T X>& v, VCT x)
{ return PRODXV<CT,T X>(CT(x),v); }

template <typename T Y>
inline PRODXV<CT,CT X> operator*(const GENVECTOR<CT X>& v, T x)
{ return PRODXV<CT,CT X>(x,v); }

template <typename T Y>
inline PRODXV<CT,CT X> operator*(const GENVECTOR<CT X>& v, CCT x)
{ return PRODXV<CT,CT X>(CT(x),v); }

template <typename T Y>
inline PRODXV<CT,CT X> operator*(const GENVECTOR<CT X>& v, VCT x)
{ return PRODXV<CT,CT X>(CT(x),v); }

// v/x

template <typename T Y>
inline PRODXV<T,T X> operator/(const GENVECTOR<T X>& v, T x)
{ return PRODXV<T,T X>(TMV_InverseOf(x),v); }

template <typename T Y>
inline PRODXV<CT,T X> operator/(const GENVECTOR<T X>& v, CT x)
{ return PRODXV<CT,T X>(TMV_InverseOf(x),v); }

template <typename T Y>
inline PRODXV<CT,T X> operator/(const GENVECTOR<T X>& v, CCT x)
{ return PRODXV<CT,T X>(TMV_InverseOf(CT(x)),v); }

template <typename T Y>
inline PRODXV<CT,T X> operator/(const GENVECTOR<T X>& v, VCT x)
{ return PRODXV<CT,T X>(TMV_InverseOf(CT(x)),v); }

template <typename T Y>
inline PRODXV<CT,CT X> operator/(const GENVECTOR<CT X>& v, T x)
{ return PRODXV<CT,CT X>(TMV_InverseOf(x),v); }

template <typename T Y>
inline PRODXV<CT,CT X> operator/(const GENVECTOR<CT X>& v, CCT x)
{ return PRODXV<CT,CT X>(TMV_InverseOf(CT(x)),v); }

template <typename T Y>
inline PRODXV<CT,CT X> operator/(const GENVECTOR<CT X>& v, VCT x)
{ return PRODXV<CT,CT X>(TMV_InverseOf(CT(x)),v); }

// -(x*v)

template <typename T, typename T2 Y>
inline PRODXV<T,T2 X> operator-(const PRODXV<T,T2 X>& pxv)
{ return PRODXV<T,T2 X>(-pxv.getX(),pxv.getV()); }

// x*(x*v)

template <typename T, typename T2 Y>
inline PRODXV<T,T2 X> operator*(const T x, const PRODXV<T,T2 X>& pxv)
{ return PRODXV<T,T2 X>(x*pxv.getX(),pxv.getV()); }

template <typename T, typename T2 Y>
inline PRODXV<CT,T2 X> operator*(const T x, const PRODXV<CT,T2 X>& pxv)
{ return PRODXV<CT,T2 X>(x*pxv.getX(),pxv.getV()); }

template <typename T Y>
inline PRODXV<CT,T X> operator*(const CT x, const PRODXV<T,T X>& pxv)
{ return PRODXV<CT,T X>(x*pxv.getX(),pxv.getV()); }

template <typename T Y>
inline PRODXV<CT,T X> operator*(const CCT x, const PRODXV<T,T X>& pxv)
{ return PRODXV<CT,T X>(CT(x)*pxv.getX(),pxv.getV()); }

template <typename T Y>
inline PRODXV<CT,T X> operator*(const VCT x, const PRODXV<T,T X>& pxv)
{ return PRODXV<CT,T X>(CT(x)*pxv.getX(),pxv.getV()); }

template <typename T, typename T2 Y>
inline PRODXV<CT,T2 X> operator*(const CCT x, const PRODXV<CT,T2 X>& pxv)
{ return PRODXV<CT,T2 X>(CT(x)*pxv.getX(),pxv.getV()); }

template <typename T, typename T2 Y>
inline PRODXV<CT,T2 X> operator*(const VCT x, const PRODXV<CT,T2 X>& pxv)
{ return PRODXV<CT,T2 X>(CT(x)*pxv.getX(),pxv.getV()); }

// (x*v)*x

template <typename T, typename T2 Y>
inline PRODXV<T,T2 X> operator*(const PRODXV<T,T2 X>& pxv, const T x)
{ return PRODXV<T,T2 X>(x*pxv.getX(),pxv.getV()); }

template <typename T, typename T2 Y>
inline PRODXV<CT,T2 X> operator*(const PRODXV<CT,T2 X>& pxv, const T x)
{ return PRODXV<CT,T2 X>(x*pxv.getX(),pxv.getV()); }

template <typename T Y>
inline PRODXV<CT,T X> operator*(const PRODXV<T,T X>& pxv, const CT x)
{ return PRODXV<CT,T X>(x*pxv.getX(),pxv.getV()); }

template <typename T Y>
inline PRODXV<CT,T X> operator*(const PRODXV<T,T X>& pxv, const CCT x)
{ return PRODXV<CT,T X>(CT(x)*pxv.getX(),pxv.getV()); }

template <typename T Y>
inline PRODXV<CT,T X> operator*(const PRODXV<T,T X>& pxv, const VCT x)
{ return PRODXV<CT,T X>(CT(x)*pxv.getX(),pxv.getV()); }

template <typename T, typename T2 Y>
inline PRODXV<CT,T2 X> operator*(const PRODXV<CT,T2 X>& pxv, const CCT x)
{ return PRODXV<CT,T2 X>(CT(x)*pxv.getX(),pxv.getV()); }

template <typename T, typename T2 Y>
inline PRODXV<CT,T2 X> operator*(const PRODXV<CT,T2 X>& pxv, const VCT x)
{ return PRODXV<CT,T2 X>(CT(x)*pxv.getX(),pxv.getV()); }

// (x*v)/x

template <typename T, typename T2 Y>
inline PRODXV<T,T X> operator/(const PRODXV<T,T2 X>& pxv, const T x)
{ return PRODXV<T,T2 X>(TMV_Divide(pxv.getX(),x),pxv.getV()); }

template <typename T, typename T2 Y>
inline PRODXV<CT,T2 X> operator/(const PRODXV<CT,T2 X>& pxv, const T x)
{ return PRODXV<CT,T2 X>(TMV_Divide(pxv.getX(),x),pxv.getV()); }

template <typename T Y>
inline PRODXV<CT,T X> operator/(const PRODXV<T,T X>& pxv, const CT x)
{ return PRODXV<CT,T X>(TMV_Divide(pxv.getX(),x),pxv.getV()); }

template <typename T Y>
inline PRODXV<CT,T X> operator/(const PRODXV<T,T X>& pxv, const CCT x)
{ return PRODXV<CT,T X>(TMV_Divide(pxv.getX(),CT(x)),pxv.getV()); }

template <typename T Y>
inline PRODXV<CT,T X> operator/(const PRODXV<T,T X>& pxv, const VCT x)
{ return PRODXV<CT,T X>(TMV_Divide(pxv.getX(),CT(x)),pxv.getV()); }

template <typename T, typename T2 Y>
inline PRODXV<CT,T2 X> operator/(const PRODXV<CT,T2 X>& pxv, const CCT x)
{ return PRODXV<CT,T2 X>(TMV_Divide(pxv.getX(),CT(x)),pxv.getV()); }

template <typename T, typename T2 Y>
inline PRODXV<CT,T2 X> operator/(const PRODXV<CT,T2 X>& pxv, const VCT x)
{ return PRODXV<CT,T2 X>(TMV_Divide(pxv.getX(),CT(x)),pxv.getV()); }

#undef X
#undef Y
