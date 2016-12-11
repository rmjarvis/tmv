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

template <typename T Y>
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

template <typename T Y>
inline PRODXM<CT,int X> operator*(CT x, const GENMATRIX<int X>& m)
{ return PRODXM<CT,int X>(x,m); }

template <typename T Y>
inline PRODXM<CT,int X> operator*(CCT x, const GENMATRIX<int X>& m)
{ return PRODXM<CT,int X>(CT(x),m); }

template <typename T Y>
inline PRODXM<CT,int X> operator*(VCT x, const GENMATRIX<int X>& m)
{ return PRODXM<CT,int X>(CT(x),m); }
#else
template <typename T Y>
inline PRODXM<T,T X> operator*(T x, const GENMATRIX<T X>& m)
{ return PRODXM<T,T X>(x,m); }

template <typename T Y>
inline PRODXM<CT,T X> operator*(CT x, const GENMATRIX<T X>& m)
{ return PRODXM<CT,T X>(x,m); }

template <typename T Y>
inline PRODXM<CT,CT X> operator*(T x, const GENMATRIX<CT X>& m)
{ return PRODXM<CT,CT X>(CT(x),m); }

template <typename T Y>
inline PRODXM<CT,T X> operator*(CCT x, const GENMATRIX<T X>& m)
{ return PRODXM<CT,T X>(CT(x),m); }

template <typename T Y>
inline PRODXM<CT,T X> operator*(VCT x, const GENMATRIX<T X>& m)
{ return PRODXM<CT,T X>(CT(x),m); }

template <typename T Y>
inline PRODXM<CT,CT X> operator*(CCT x, const GENMATRIX<CT X>& m)
{ return PRODXM<CT,CT X>(CT(x),m); }

template <typename T Y>
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

template <typename T Y>
inline PRODXM<CT,int X> operator*(const GENMATRIX<int X>& m, CT x)
{ return PRODXM<CT,int X>(x,m); }

template <typename T Y>
inline PRODXM<CT,int X> operator*(const GENMATRIX<int X>& m, CCT x)
{ return PRODXM<CT,int X>(CT(x),m); }

template <typename T Y>
inline PRODXM<CT,int X> operator*(const GENMATRIX<int X>& m, VCT x)
{ return PRODXM<CT,int X>(CT(x),m); }
#else
template <typename T Y>
inline PRODXM<T,T X> operator*(const GENMATRIX<T X>& m, T x)
{ return PRODXM<T,T X>(x,m); }

template <typename T Y>
inline PRODXM<CT,T X> operator*(const GENMATRIX<T X>& m, CT x)
{ return PRODXM<CT,T X>(x,m); }

template <typename T Y>
inline PRODXM<CT,T X> operator*(const GENMATRIX<T X>& m, CCT x)
{ return PRODXM<CT,T X>(CT(x),m); }

template <typename T Y>
inline PRODXM<CT,T X> operator*(const GENMATRIX<T X>& m, VCT x)
{ return PRODXM<CT,T X>(CT(x),m); }

template <typename T Y>
inline PRODXM<CT,CT X> operator*(const GENMATRIX<CT X>& m, CCT x)
{ return PRODXM<CT,CT X>(CT(x),m); }

template <typename T Y>
inline PRODXM<CT,CT X> operator*(const GENMATRIX<CT X>& m, VCT x)
{ return PRODXM<CT,CT X>(CT(x),m); }

template <typename T Y>
inline PRODXM<CT,CT X> operator*(const GENMATRIX<CT X>& m, T x)
{ return PRODXM<CT,CT X>(CT(x),m); }
#endif

// m/x
#ifdef INTT
inline PRODXM<double,int> operator/(const GENMATRIX<int>& m, int x)
{ return PRODXM<double,int>(TMV_InverseOf(double(x)),v); }

inline PRODXM<float,int> operator/(const GENMATRIX<int>& m, float x)
{ return PRODXM<float,int>(TMV_InverseOf(float(x)),v); }

inline PRODXM<double,int> operator/(const GENMATRIX<int>& m, double x)
{ return PRODXM<double,int>(TMV_InverseOf(x),v); }

inline PRODXM<long double,int> operator/(
    const GENMATRIX<int>& m, long double x)
{ return PRODXM<long double,int>(TMV_InverseOf(x),v); }

template <typename T Y>
inline PRODXM<CT,int X> operator/(const GENMATRIX<int X>& m, CT x)
{ return PRODXM<CT,int X>(TMV_InverseOf(x),m); }

template <typename T Y>
inline PRODXM<CT,int X> operator/(const GENMATRIX<int X>& m, CCT x)
{ return PRODXM<CT,int X>(TMV_InverseOf(CT(x)),m); }

template <typename T Y>
inline PRODXM<CT,int X> operator/(const GENMATRIX<int X>& m, VCT x)
{ return PRODXM<CT,int X>(TMV_InverseOf(CT(x)),m); }
#undef INTT
#else
template <typename T Y>
inline PRODXM<T,T X> operator/(const GENMATRIX<T X>& m, T x)
{ return PRODXM<T,T X>(TMV_InverseOf(x),m); }

template <typename T Y>
inline PRODXM<CT,T X> operator/(const GENMATRIX<T X>& m, CT x)
{ return PRODXM<CT,T X>(TMV_InverseOf(x),m); }

template <typename T Y>
inline PRODXM<CT,T X> operator/(const GENMATRIX<T X>& m, CCT x)
{ return PRODXM<CT,T X>(TMV_InverseOf(CT(x)),m); }

template <typename T Y>
inline PRODXM<CT,T X> operator/(const GENMATRIX<T X>& m, VCT x)
{ return PRODXM<CT,T X>(TMV_InverseOf(CT(x)),m); }

template <typename T Y>
inline PRODXM<CT,CT X> operator/(const GENMATRIX<CT X>& m, CCT x)
{ return PRODXM<CT,CT X>(TMV_InverseOf(CT(x)),m); }

template <typename T Y>
inline PRODXM<CT,CT X> operator/(const GENMATRIX<CT X>& m, VCT x)
{ return PRODXM<CT,CT X>(TMV_InverseOf(CT(x)),m); }

template <typename T Y>
inline PRODXM<CT,CT X> operator/(const GENMATRIX<CT X>& m, T x)
{ return PRODXM<CT,CT X>(TMV_InverseOf(x),m); }

#endif

// -(x*m)

template <typename T, typename T2 Y>
inline PRODXM<T,T2 X> operator-(const PRODXM<T,T2 X>& pxm)
{ return PRODXM<T,T2 X>(-pxm.getX(),pxm GETM); }

// x*(x*m)

template <typename T, typename T2 Y>
inline PRODXM<T,T2 X> operator*(const T x, const PRODXM<T,T2 X>& pxm)
{ return PRODXM<T,T2 X>(x*pxm.getX(),pxm GETM); }

template <typename T, typename T2 Y>
inline PRODXM<CT,T2 X> operator*(const T x, const PRODXM<CT,T2 X>& pxm)
{ return PRODXM<CT,T2 X>(x*pxm.getX(),pxm GETM); }

template <typename T Y>
inline PRODXM<CT,T X> operator*(const CT x, const PRODXM<T,T X>& pxm)
{ return PRODXM<CT,T X>(x*pxm.getX(),pxm GETM); }

template <typename T Y>
inline PRODXM<CT,T X> operator*(const CCT x, const PRODXM<T,T X>& pxm)
{ return PRODXM<CT,T X>(CT(x)*pxm.getX(),pxm GETM); }

template <typename T Y>
inline PRODXM<CT,T X> operator*(const VCT x, const PRODXM<T,T X>& pxm)
{ return PRODXM<CT,T X>(CT(x)*pxm.getX(),pxm GETM); }

template <typename T, typename T2 Y>
inline PRODXM<CT,T2 X> operator*(const CCT x, const PRODXM<CT,T2 X>& pxm)
{ return PRODXM<CT,T2 X>(CT(x)*pxm.getX(),pxm GETM); }

template <typename T, typename T2 Y>
inline PRODXM<CT,T2 X> operator*(const VCT x, const PRODXM<CT,T2 X>& pxm)
{ return PRODXM<CT,T2 X>(CT(x)*pxm.getX(),pxm GETM); }

// (x*m)*x

template <typename T, typename T2 Y>
inline PRODXM<T,T2 X> operator*(const PRODXM<T,T2 X>& pxm, const T x)
{ return PRODXM<T,T2 X>(x*pxm.getX(),pxm GETM); }

template <typename T, typename T2 Y>
inline PRODXM<CT,T2 X> operator*(const PRODXM<CT,T2 X>& pxm, const T x)
{ return PRODXM<CT,T2 X>(x*pxm.getX(),pxm GETM); }

template <typename T Y>
inline PRODXM<CT,T X> operator*(const PRODXM<T,T X>& pxm, const CT x)
{ return PRODXM<CT,T X>(x*pxm.getX(),pxm GETM); }

template <typename T Y>
inline PRODXM<CT,T X> operator*(const PRODXM<T,T X>& pxm, const CCT x)
{ return PRODXM<CT,T X>(CT(x)*pxm.getX(),pxm GETM); }

template <typename T Y>
inline PRODXM<CT,T X> operator*(const PRODXM<T,T X>& pxm, const VCT x)
{ return PRODXM<CT,T X>(CT(x)*pxm.getX(),pxm GETM); }

template <typename T, typename T2 Y>
inline PRODXM<CT,T2 X> operator*(const PRODXM<CT,T2 X>& pxm, const CCT x)
{ return PRODXM<CT,T2 X>(CT(x)*pxm.getX(),pxm GETM); }

template <typename T, typename T2 Y>
inline PRODXM<CT,T2 X> operator*(const PRODXM<CT,T2 X>& pxm, const VCT x)
{ return PRODXM<CT,T2 X>(CT(x)*pxm.getX(),pxm GETM); }

// (x*m)/x

template <typename T, typename T2 Y>
inline PRODXM<T,T2 X> operator/(const PRODXM<T,T2 X>& pxm, const T x)
{ return PRODXM<T,T2 X>(TMV_Divide(pxm.getX(),x),pxm GETM); }

template <typename T, typename T2 Y>
inline PRODXM<CT,T2 X> operator/(const PRODXM<CT,T2 X>& pxm, const T x)
{ return PRODXM<CT,T2 X>(TMV_Divide(pxm.getX(),x),pxm GETM); }

template <typename T Y>
inline PRODXM<CT,T X> operator/(const PRODXM<T,T X>& pxm, const CT x)
{ return PRODXM<CT,T X>(TMV_Divide(pxm.getX(),x),pxm GETM); }

template <typename T Y>
inline PRODXM<CT,T X> operator/(const PRODXM<T,T X>& pxm, const CCT x)
{ return PRODXM<CT,T X>(TMV_Divide(pxm.getX(),CT(x)),pxm GETM); }

template <typename T Y>
inline PRODXM<CT,T X> operator/(const PRODXM<T,T X>& pxm, const VCT x)
{ return PRODXM<CT,T X>(TMV_Divide(pxm.getX(),CT(x)),pxm GETM); }

template <typename T, typename T2 Y>
inline PRODXM<CT,T2 X> operator/(const PRODXM<CT,T2 X>& pxm, const CCT x)
{ return PRODXM<CT,T2 X>(TMV_Divide(pxm.getX(),CT(x)),pxm GETM); }

template <typename T, typename T2 Y>
inline PRODXM<CT,T2 X> operator/(const PRODXM<CT,T2 X>& pxm, const VCT x)
{ return PRODXM<CT,T2 X>(TMV_Divide(pxm.getX(),CT(x)),pxm GETM); }

#undef X
#undef Y
#undef GETM
