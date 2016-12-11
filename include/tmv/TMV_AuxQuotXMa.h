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
#define GETM .getM()
#endif

// -(x/m)

template <typename T, typename T2 Y>
inline QUOTXM<T,T2 X> operator-(const QUOTXM<T,T2 X>& qxm)
{ return QUOTXM<T,T2 X>(-qxm.getX(),qxm GETM); }

// x*(x/m)

template <typename T, typename T2 Y>
inline QUOTXM<T,T2 X> operator*(const T x, const QUOTXM<T,T2 X>& qxm)
{ return QUOTXM<T,T2 X>(x*qxm.getX(),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator*(const T x, const QUOTXM<CT,T2 X>& qxm)
{ return QUOTXM<CT,T2 X>(x*qxm.getX(),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator*(const CT x, const QUOTXM<T,T X>& qxm)
{ return QUOTXM<CT,T X>(x*qxm.getX(),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator*(const CCT x, const QUOTXM<T,T X>& qxm)
{ return QUOTXM<CT,T X>(CT(x)*qxm.getX(),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator*(const VCT x, const QUOTXM<T,T X>& qxm)
{ return QUOTXM<CT,T X>(CT(x)*qxm.getX(),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator*(const CCT x, const QUOTXM<CT,T2 X>& qxm)
{ return QUOTXM<CT,T2 X>(CT(x)*qxm.getX(),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator*(const VCT x, const QUOTXM<CT,T2 X>& qxm)
{ return QUOTXM<CT,T2 X>(CT(x)*qxm.getX(),qxm GETM); }

// (x/m)*x

template <typename T, typename T2 Y>
inline QUOTXM<T,T2 X> operator*(const QUOTXM<T,T2 X>& qxm, const T x)
{ return QUOTXM<T,T2 X>(x*qxm.getX(),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator*(const QUOTXM<CT,T2 X>& qxm, const T x)
{ return QUOTXM<CT,T2 X>(x*qxm.getX(),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator*(const QUOTXM<T,T X>& qxm, const CT x)
{ return QUOTXM<CT,T X>(x*qxm.getX(),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator*(const QUOTXM<T,T X>& qxm, const CCT x)
{ return QUOTXM<CT,T X>(CT(x)*qxm.getX(),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator*(const QUOTXM<T,T X>& qxm, const VCT x)
{ return QUOTXM<CT,T X>(CT(x)*qxm.getX(),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator*(const QUOTXM<CT,T2 X>& qxm, const CCT x)
{ return QUOTXM<CT,T2 X>(CT(x)*qxm.getX(),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator*(const QUOTXM<CT,T2 X>& qxm, const VCT x)
{ return QUOTXM<CT,T2 X>(CT(x)*qxm.getX(),qxm GETM); }

// (x/m)/x

template <typename T, typename T2 Y>
inline QUOTXM<T,T2 X> operator/(const QUOTXM<T,T2 X>& qxm, const T x)
{ return QUOTXM<T,T2 X>(TMV_Divide(qxm.getX(),x),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator/(const QUOTXM<CT,T2 X>& qxm, const T x)
{ return QUOTXM<CT,T2 X>(TMV_Divide(qxm.getX(),x),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator/(const QUOTXM<T,T X>& qxm, const CT x)
{ return QUOTXM<CT,T X>(TMV_Divide(qxm.getX(),x),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator/(const QUOTXM<T,T X>& qxm, const CCT x)
{ return QUOTXM<CT,T X>(TMV_Divide(qxm.getX(),CT(x)),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator/(const QUOTXM<T,T X>& qxm, const VCT x)
{ return QUOTXM<CT,T X>(TMV_Divide(qxm.getX(),CT(x)),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator/(const QUOTXM<CT,T2 X>& qxm, const CCT x)
{ return QUOTXM<CT,T2 X>(TMV_Divide(qxm.getX(),CT(x)),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator/(const QUOTXM<CT,T2 X>& qxm, const VCT x)
{ return QUOTXM<CT,T2 X>(TMV_Divide(qxm.getX(),CT(x)),qxm GETM); }

#ifdef QUOTXM_1

// -(1/m)

template <typename T, typename T2 Y>
inline QUOTXM<T,T2 X> operator-(const QUOTXM_1<T,T2 X>& qxm)
{ return QUOTXM<T,T2 X>(T(-1),qxm GETM); }

// x*(1/m)

template <typename T, typename T2 Y>
inline QUOTXM<T,T2 X> operator*(const T x, const QUOTXM_1<T,T2 X>& qxm)
{ return QUOTXM<T,T2 X>(x,qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator*(const T x, const QUOTXM_1<CT,T2 X>& qxm)
{ return QUOTXM<CT,T2 X>(x,qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator*(const CT x, const QUOTXM_1<T,T X>& qxm)
{ return QUOTXM<CT,T X>(x,qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator*(const CCT x, const QUOTXM_1<T,T X>& qxm)
{ return QUOTXM<CT,T X>(CT(x),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator*(const VCT x, const QUOTXM_1<T,T X>& qxm)
{ return QUOTXM<CT,T X>(CT(x),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator*(const CCT x, const QUOTXM_1<CT,T2 X>& qxm)
{ return QUOTXM<CT,T2 X>(CT(x),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator*(const VCT x, const QUOTXM_1<CT,T2 X>& qxm)
{ return QUOTXM<CT,T2 X>(CT(x),qxm GETM); }

// (1/m)*x

template <typename T, typename T2 Y>
inline QUOTXM<T,T2 X> operator*(const QUOTXM_1<T,T2 X>& qxm, const T x)
{ return QUOTXM<T,T2 X>(x,qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator*(const QUOTXM_1<CT,T2 X>& qxm, const T x)
{ return QUOTXM<CT,T2 X>(x,qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator*(const QUOTXM_1<T,T X>& qxm, const CT x)
{ return QUOTXM<CT,T X>(x,qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator*(const QUOTXM_1<T,T X>& qxm, const CCT x)
{ return QUOTXM<CT,T X>(CT(x),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator*(const QUOTXM_1<T,T X>& qxm, const VCT x)
{ return QUOTXM<CT,T X>(CT(x),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator*(const QUOTXM_1<CT,T2 X>& qxm, const CCT x)
{ return QUOTXM<CT,T2 X>(CT(x),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator*(const QUOTXM_1<CT,T2 X>& qxm, const VCT x)
{ return QUOTXM<CT,T2 X>(CT(x),qxm GETM); }

// (1/m)/x

template <typename T, typename T2 Y>
inline QUOTXM<T,T2 X> operator/(const QUOTXM_1<T,T2 X>& qxm, const T x)
{ return QUOTXM<T,T2 X>(TMV_InverseOf(x),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator/(const QUOTXM_1<CT,T2 X>& qxm, const T x)
{ return QUOTXM<CT,T2 X>(TMV_InverseOf(x),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator/(const QUOTXM_1<T,T X>& qxm, const CT x)
{ return QUOTXM<CT,T X>(TMV_InverseOf(x),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator/(const QUOTXM_1<T,T X>& qxm, const CCT x)
{ return QUOTXM<CT,T X>(TMV_InverseOf(CT(x)),qxm GETM); }

template <typename T Y>
inline QUOTXM<CT,T X> operator/(const QUOTXM_1<T,T X>& qxm, const VCT x)
{ return QUOTXM<CT,T X>(TMV_InverseOf(CT(x)),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator/(const QUOTXM_1<CT,T2 X>& qxm, const CCT x)
{ return QUOTXM<CT,T2 X>(TMV_InverseOf(CT(x)),qxm GETM); }

template <typename T, typename T2 Y>
inline QUOTXM<CT,T2 X> operator/(const QUOTXM_1<CT,T2 X>& qxm, const VCT x)
{ return QUOTXM<CT,T2 X>(TMV_InverseOf(CT(x)),qxm GETM); }

#undef QUOTXM_1
#endif

#undef X
#undef Y
#undef GETM
