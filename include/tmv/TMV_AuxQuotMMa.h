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

#ifndef GETM1
#define GETM1 .getM1()
#endif

#ifndef GETM2
#define GETM2 .getM2()
#endif

// x*(x*m/m)

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<T,T1,T2 X3> operator*(T x, const QUOTMM<T,T1,T2 X3>& qmm)
{ return QUOTMM<T,T1,T2 X3>(qmm.getX()*x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator*(CT x, const QUOTMM<CT,T,T X3>& qmm)
{ return QUOTMM<CT,T,T X3>(qmm.getX()*x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator*(CCT x, const QUOTMM<CT,T,T X3>& qmm)
{ return QUOTMM<CT,T,T X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator*(VCT x, const QUOTMM<CT,T,T X3>& qmm)
{ return QUOTMM<CT,T,T X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator*(T x, const QUOTMM<CT,T1,T2 X3>& qmm)
{ return QUOTMM<CT,T1,T2 X3>(qmm.getX()*x,qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator*(CCT x, const QUOTMM<CT,T1,T2 X3>& qmm)
{ return QUOTMM<CT,T1,T2 X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator*(VCT x, const QUOTMM<CT,T1,T2 X3>& qmm)
{ return QUOTMM<CT,T1,T2 X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

// (x*m/m)*x

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<T,T1,T2 X3> operator*(const QUOTMM<T,T1,T2 X3>& qmm, T x)
{ return QUOTMM<T,T1,T2 X3>(qmm.getX()*x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator*(const QUOTMM<CT,T,T X3>& qmm, CT x)
{ return QUOTMM<CT,T,T X3>(qmm.getX()*x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator*(const QUOTMM<CT,T,T X3>& qmm, CCT x)
{ return QUOTMM<CT,T,T X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator*(const QUOTMM<CT,T,T X3>& qmm, VCT x)
{ return QUOTMM<CT,T,T X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator*(const QUOTMM<CT,T1,T2 X3>& qmm, T x)
{ return QUOTMM<CT,T1,T2 X3>(qmm.getX()*x,qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator*(const QUOTMM<CT,T1,T2 X3>& qmm, CCT x)
{ return QUOTMM<CT,T1,T2 X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator*(const QUOTMM<CT,T1,T2 X3>& qmm, VCT x)
{ return QUOTMM<CT,T1,T2 X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

// (x*m/m)/x

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<T,T1,T2 X3> operator/(const QUOTMM<T,T1,T2 X3>& qmm, T x)
{ return QUOTMM<T,T1,T2 X3>(TMV_Divide(qmm.getX(),x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator/(const QUOTMM<CT,T,T X3>& qmm, CT x)
{ return QUOTMM<CT,T,T X3>(TMV_Divide(qmm.getX(),x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator/(const QUOTMM<CT,T,T X3>& qmm, CCT x)
{ return QUOTMM<CT,T,T X3>(TMV_Divide(qmm.getX(),CT(x)),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator/(const QUOTMM<CT,T,T X3>& qmm, VCT x)
{ return QUOTMM<CT,T,T X3>(TMV_Divide(qmm.getX(),CT(x)),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator/(const QUOTMM<CT,T1,T2 X3>& qmm, T x)
{ return QUOTMM<CT,T1,T2 X3>(TMV_Divide(qmm.getX(),x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator/(const QUOTMM<CT,T1,T2 X3>& qmm, CCT x)
{ return QUOTMM<CT,T1,T2 X3>(TMV_Divide(qmm.getX(),CT(x)),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator/(const QUOTMM<CT,T1,T2 X3>& qmm, VCT x)
{ return QUOTMM<CT,T1,T2 X3>(TMV_Divide(qmm.getX(),CT(x)),qmm GETM1,qmm GETM2); }

// x*(x*m%m)

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<T,T1,T2 X3> operator*(T x, const RQUOTMM<T,T1,T2 X3>& qmm)
{ return RQUOTMM<T,T1,T2 X3>(qmm.getX()*x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator*(CT x, const RQUOTMM<CT,T,T X3>& qmm)
{ return RQUOTMM<CT,T,T X3>(qmm.getX()*x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator*(CCT x, const RQUOTMM<CT,T,T X3>& qmm)
{ return RQUOTMM<CT,T,T X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator*(VCT x, const RQUOTMM<CT,T,T X3>& qmm)
{ return RQUOTMM<CT,T,T X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator*(T x, const RQUOTMM<CT,T1,T2 X3>& qmm)
{ return RQUOTMM<CT,T1,T2 X3>(qmm.getX()*x,qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator*(CCT x, const RQUOTMM<CT,T1,T2 X3>& qmm)
{ return RQUOTMM<CT,T1,T2 X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator*(VCT x, const RQUOTMM<CT,T1,T2 X3>& qmm)
{ return RQUOTMM<CT,T1,T2 X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

// (x*m%m)*x

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<T,T1,T2 X3> operator*(const RQUOTMM<T,T1,T2 X3>& qmm, T x)
{ return RQUOTMM<T,T1,T2 X3>(qmm.getX()*x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator*(const RQUOTMM<CT,T,T X3>& qmm, CT x)
{ return RQUOTMM<CT,T,T X3>(qmm.getX()*x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator*(const RQUOTMM<CT,T,T X3>& qmm, CCT x)
{ return RQUOTMM<CT,T,T X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator*(const RQUOTMM<CT,T,T X3>& qmm, VCT x)
{ return RQUOTMM<CT,T,T X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator*(const RQUOTMM<CT,T1,T2 X3>& qmm, T x)
{ return RQUOTMM<CT,T1,T2 X3>(qmm.getX()*x,qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator*(const RQUOTMM<CT,T1,T2 X3>& qmm, CCT x)
{ return RQUOTMM<CT,T1,T2 X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator*(const RQUOTMM<CT,T1,T2 X3>& qmm, VCT x)
{ return RQUOTMM<CT,T1,T2 X3>(qmm.getX()*CT(x),qmm GETM1,qmm GETM2); }

// (x*m%m)/x

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<T,T1,T2 X3> operator/(const RQUOTMM<T,T1,T2 X3>& qmm, T x)
{ return RQUOTMM<T,T1,T2 X3>(TMV_Divide(qmm.getX(),x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator/(const RQUOTMM<CT,T,T X3>& qmm, CT x)
{ return RQUOTMM<CT,T,T X3>(TMV_Divide(qmm.getX(),x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator/(const RQUOTMM<CT,T,T X3>& qmm, CCT x)
{ return RQUOTMM<CT,T,T X3>(TMV_Divide(qmm.getX(),CT(x)),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator/(const RQUOTMM<CT,T,T X3>& qmm, VCT x)
{ return RQUOTMM<CT,T,T X3>(TMV_Divide(qmm.getX(),CT(x)),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator/(const RQUOTMM<CT,T1,T2 X3>& qmm, T x)
{ return RQUOTMM<CT,T1,T2 X3>(TMV_Divide(qmm.getX(),x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator/(const RQUOTMM<CT,T1,T2 X3>& qmm, CCT x)
{ return RQUOTMM<CT,T1,T2 X3>(TMV_Divide(qmm.getX(),CT(x)),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator/(const RQUOTMM<CT,T1,T2 X3>& qmm, VCT x)
{ return RQUOTMM<CT,T1,T2 X3>(TMV_Divide(qmm.getX(),CT(x)),qmm GETM1,qmm GETM2); }

#ifdef QUOTMM_1

// x*(m/m)

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<T,T1,T2 X3> operator*(T x, const QUOTMM_1<T,T1,T2 X3>& qmm)
{ return QUOTMM<T,T1,T2 X3>(x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator*(CT x, const QUOTMM_1<CT,T,T X3>& qmm)
{ return QUOTMM<CT,T,T X3>(x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator*(CCT x, const QUOTMM_1<CT,T,T X3>& qmm)
{ return QUOTMM<CT,T,T X3>(CT(x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator*(VCT x, const QUOTMM_1<CT,T,T X3>& qmm)
{ return QUOTMM<CT,T,T X3>(CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator*(T x, const QUOTMM_1<CT,T1,T2 X3>& qmm)
{ return QUOTMM<CT,T1,T2 X3>(x,qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator*(CCT x, const QUOTMM_1<CT,T1,T2 X3>& qmm)
{ return QUOTMM<CT,T1,T2 X3>(CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator*(VCT x, const QUOTMM_1<CT,T1,T2 X3>& qmm)
{ return QUOTMM<CT,T1,T2 X3>(CT(x),qmm GETM1,qmm GETM2); }

// (m/m)*x

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<T,T1,T2 X3> operator*(const QUOTMM_1<T,T1,T2 X3>& qmm, T x)
{ return QUOTMM<T,T1,T2 X3>(x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator*(const QUOTMM_1<CT,T,T X3>& qmm, CT x)
{ return QUOTMM<CT,T,T X3>(x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator*(const QUOTMM_1<CT,T,T X3>& qmm, CCT x)
{ return QUOTMM<CT,T,T X3>(CT(x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator*(const QUOTMM_1<CT,T,T X3>& qmm, VCT x)
{ return QUOTMM<CT,T,T X3>(CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator*(const QUOTMM_1<CT,T1,T2 X3>& qmm, T x)
{ return QUOTMM<CT,T1,T2 X3>(x,qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator*(const QUOTMM_1<CT,T1,T2 X3>& qmm, CCT x)
{ return QUOTMM<CT,T1,T2 X3>(CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator*(const QUOTMM_1<CT,T1,T2 X3>& qmm, VCT x)
{ return QUOTMM<CT,T1,T2 X3>(CT(x),qmm GETM1,qmm GETM2); }

// (m/m)/x

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<T,T1,T2 X3> operator/(const QUOTMM_1<T,T1,T2 X3>& qmm, T x)
{ return QUOTMM<T,T1,T2 X3>(TMV_InverseOf(x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator/(const QUOTMM_1<CT,T,T X3>& qmm, CT x)
{ return QUOTMM<CT,T,T X3>(TMV_InverseOf(x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator/(const QUOTMM_1<CT,T,T X3>& qmm, CCT x)
{ return QUOTMM<CT,T,T X3>(TMV_InverseOf(CT(x)),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline QUOTMM<CT,T,T X3> operator/(const QUOTMM_1<CT,T,T X3>& qmm, VCT x)
{ return QUOTMM<CT,T,T X3>(TMV_InverseOf(CT(x)),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator/(const QUOTMM_1<CT,T1,T2 X3>& qmm, T x)
{ return QUOTMM<CT,T1,T2 X3>(TMV_InverseOf(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator/(const QUOTMM_1<CT,T1,T2 X3>& qmm, CCT x)
{ return QUOTMM<CT,T1,T2 X3>(TMV_InverseOf(CT(x)),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline QUOTMM<CT,T1,T2 X3> operator/(const QUOTMM_1<CT,T1,T2 X3>& qmm, VCT x)
{ return QUOTMM<CT,T1,T2 X3>(TMV_InverseOf(CT(x)),qmm GETM1,qmm GETM2); }

// x*(m%m)

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<T,T1,T2 X3> operator*(T x, const RQUOTMM_1<T,T1,T2 X3>& qmm)
{ return RQUOTMM<T,T1,T2 X3>(x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator*(CT x, const RQUOTMM_1<CT,T,T X3>& qmm)
{ return RQUOTMM<CT,T,T X3>(x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator*(CCT x, const RQUOTMM_1<CT,T,T X3>& qmm)
{ return RQUOTMM<CT,T,T X3>(CT(x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator*(VCT x, const RQUOTMM_1<CT,T,T X3>& qmm)
{ return RQUOTMM<CT,T,T X3>(CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator*(T x, const RQUOTMM_1<CT,T1,T2 X3>& qmm)
{ return RQUOTMM<CT,T1,T2 X3>(x,qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator*(CCT x, const RQUOTMM_1<CT,T1,T2 X3>& qmm)
{ return RQUOTMM<CT,T1,T2 X3>(CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator*(VCT x, const RQUOTMM_1<CT,T1,T2 X3>& qmm)
{ return RQUOTMM<CT,T1,T2 X3>(CT(x),qmm GETM1,qmm GETM2); }

// (m%m)*x

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<T,T1,T2 X3> operator*(const RQUOTMM_1<T,T1,T2 X3>& qmm, T x)
{ return RQUOTMM<T,T1,T2 X3>(x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator*(const RQUOTMM_1<CT,T,T X3>& qmm, CT x)
{ return RQUOTMM<CT,T,T X3>(x,qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator*(const RQUOTMM_1<CT,T,T X3>& qmm, CCT x)
{ return RQUOTMM<CT,T,T X3>(CT(x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator*(const RQUOTMM_1<CT,T,T X3>& qmm, VCT x)
{ return RQUOTMM<CT,T,T X3>(CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator*(const RQUOTMM_1<CT,T1,T2 X3>& qmm, T x)
{ return RQUOTMM<CT,T1,T2 X3>(x,qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator*(const RQUOTMM_1<CT,T1,T2 X3>& qmm, CCT x)
{ return RQUOTMM<CT,T1,T2 X3>(CT(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator*(const RQUOTMM_1<CT,T1,T2 X3>& qmm, VCT x)
{ return RQUOTMM<CT,T1,T2 X3>(CT(x),qmm GETM1,qmm GETM2); }

// (m%m)/x

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<T,T1,T2 X3> operator/(const RQUOTMM_1<T,T1,T2 X3>& qmm, T x)
{ return RQUOTMM<T,T1,T2 X3>(TMV_InverseOf(x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator/(const RQUOTMM_1<CT,T,T X3>& qmm, CT x)
{ return RQUOTMM<CT,T,T X3>(TMV_InverseOf(x),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator/(const RQUOTMM_1<CT,T,T X3>& qmm, CCT x)
{ return RQUOTMM<CT,T,T X3>(TMV_InverseOf(CT(x)),qmm GETM1,qmm GETM2); }

template <typename T Y>
inline RQUOTMM<CT,T,T X3> operator/(const RQUOTMM_1<CT,T,T X3>& qmm, VCT x)
{ return RQUOTMM<CT,T,T X3>(TMV_InverseOf(CT(x)),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator/(const RQUOTMM_1<CT,T1,T2 X3>& qmm, T x)
{ return RQUOTMM<CT,T1,T2 X3>(TMV_InverseOf(x),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator/(const RQUOTMM_1<CT,T1,T2 X3>& qmm, CCT x)
{ return RQUOTMM<CT,T1,T2 X3>(TMV_InverseOf(CT(x)),qmm GETM1,qmm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline RQUOTMM<CT,T1,T2 X3> operator/(const RQUOTMM_1<CT,T1,T2 X3>& qmm, VCT x)
{ return RQUOTMM<CT,T1,T2 X3>(TMV_InverseOf(CT(x)),qmm GETM1,qmm GETM2); }

#undef QUOTMM_1
#endif

#undef X3
#undef Y
#undef GETM1
#undef GETM2


