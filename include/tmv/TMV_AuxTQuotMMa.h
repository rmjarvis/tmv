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


#ifndef CT
#define CT std::complex<T>
#endif
#ifndef CCT
#define CCT ConjRef<std::complex<T> >
#endif
#ifndef VCT
#define VCT VarConjRef<std::complex<T> >
#endif

// x*(x*m/m)

template <typename T, typename T1, typename T2>
inline TQUOTMM<T,T1,T2> operator*(T x, const TQUOTMM<T,T1,T2>& qmm)
{ return TQUOTMM<T,T1,T2>(qmm.getX()*x,qmm.getP(),qmm.getM2()); }

template <typename T>
inline TQUOTMM<CT,T,T> operator*(CT x, const TQUOTMM<CT,T,T>& qmm)
{ return TQUOTMM<CT,T,T>(qmm.getX()*x,qmm.getP(),qmm.getM2()); }

template <typename T>
inline TQUOTMM<CT,T,T> operator*(CCT x, const TQUOTMM<CT,T,T>& qmm)
{ return TQUOTMM<CT,T,T>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

template <typename T>
inline TQUOTMM<CT,T,T> operator*(VCT x, const TQUOTMM<CT,T,T>& qmm)
{ return TQUOTMM<CT,T,T>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TQUOTMM<CT,T1,T2> operator*(T x, const TQUOTMM<CT,T1,T2>& qmm)
{ return TQUOTMM<CT,T1,T2>(qmm.getX()*x,qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TQUOTMM<CT,T1,T2> operator*(CCT x, const TQUOTMM<CT,T1,T2>& qmm)
{ return TQUOTMM<CT,T1,T2>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TQUOTMM<CT,T1,T2> operator*(VCT x, const TQUOTMM<CT,T1,T2>& qmm)
{ return TQUOTMM<CT,T1,T2>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

// (x*m/m)*x

template <typename T, typename T1, typename T2>
inline TQUOTMM<T,T1,T2> operator*(const TQUOTMM<T,T1,T2>& qmm, T x)
{ return TQUOTMM<T,T1,T2>(qmm.getX()*x,qmm.getP(),qmm.getM2()); }

template <typename T>
inline TQUOTMM<CT,T,T> operator*(const TQUOTMM<CT,T,T>& qmm, CT x)
{ return TQUOTMM<CT,T,T>(qmm.getX()*x,qmm.getP(),qmm.getM2()); }

template <typename T>
inline TQUOTMM<CT,T,T> operator*(const TQUOTMM<CT,T,T>& qmm, CCT x)
{ return TQUOTMM<CT,T,T>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

template <typename T>
inline TQUOTMM<CT,T,T> operator*(const TQUOTMM<CT,T,T>& qmm, VCT x)
{ return TQUOTMM<CT,T,T>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TQUOTMM<CT,T1,T2> operator*(const TQUOTMM<CT,T1,T2>& qmm, T x)
{ return TQUOTMM<CT,T1,T2>(qmm.getX()*x,qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TQUOTMM<CT,T1,T2> operator*(const TQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return TQUOTMM<CT,T1,T2>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TQUOTMM<CT,T1,T2> operator*(const TQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return TQUOTMM<CT,T1,T2>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

// (x*m/m)/x

template <typename T, typename T1, typename T2>
inline TQUOTMM<T,T1,T2> operator/(const TQUOTMM<T,T1,T2>& qmm, T x)
{ return TQUOTMM<T,T1,T2>(TMV_Divide(qmm.getX(),x),qmm.getP(),qmm.getM2()); }

template <typename T>
inline TQUOTMM<CT,T,T> operator/(const TQUOTMM<CT,T,T>& qmm, CT x)
{ return TQUOTMM<CT,T,T>(TMV_Divide(qmm.getX(),x),qmm.getP(),qmm.getM2()); }

template <typename T>
inline TQUOTMM<CT,T,T> operator/(const TQUOTMM<CT,T,T>& qmm, CCT x)
{ return TQUOTMM<CT,T,T>(TMV_Divide(qmm.getX(),CT(x)),qmm.getP(),qmm.getM2()); }

template <typename T>
inline TQUOTMM<CT,T,T> operator/(const TQUOTMM<CT,T,T>& qmm, VCT x)
{ return TQUOTMM<CT,T,T>(TMV_Divide(qmm.getX(),CT(x)),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TQUOTMM<CT,T1,T2> operator/(const TQUOTMM<CT,T1,T2>& qmm, T x)
{ return TQUOTMM<CT,T1,T2>(TMV_Divide(qmm.getX(),x),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TQUOTMM<CT,T1,T2> operator/(const TQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return TQUOTMM<CT,T1,T2>(TMV_Divide(qmm.getX(),CT(x)),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TQUOTMM<CT,T1,T2> operator/(const TQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return TQUOTMM<CT,T1,T2>(TMV_Divide(qmm.getX(),CT(x)),qmm.getP(),qmm.getM2()); }

// x*(x*m%m)

template <typename T, typename T1, typename T2>
inline TRQUOTMM<T,T1,T2> operator*(T x, const TRQUOTMM<T,T1,T2>& qmm)
{ return TRQUOTMM<T,T1,T2>(qmm.getX()*x,qmm.getP(),qmm.getM2()); }

template <typename T>
inline TRQUOTMM<CT,T,T> operator*(CT x, const TRQUOTMM<CT,T,T>& qmm)
{ return TRQUOTMM<CT,T,T>(qmm.getX()*x,qmm.getP(),qmm.getM2()); }

template <typename T>
inline TRQUOTMM<CT,T,T> operator*(CCT x, const TRQUOTMM<CT,T,T>& qmm)
{ return TRQUOTMM<CT,T,T>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

template <typename T>
inline TRQUOTMM<CT,T,T> operator*(VCT x, const TRQUOTMM<CT,T,T>& qmm)
{ return TRQUOTMM<CT,T,T>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TRQUOTMM<CT,T1,T2> operator*(T x, const TRQUOTMM<CT,T1,T2>& qmm)
{ return TRQUOTMM<CT,T1,T2>(qmm.getX()*x,qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TRQUOTMM<CT,T1,T2> operator*(CCT x, const TRQUOTMM<CT,T1,T2>& qmm)
{ return TRQUOTMM<CT,T1,T2>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TRQUOTMM<CT,T1,T2> operator*(VCT x, const TRQUOTMM<CT,T1,T2>& qmm)
{ return TRQUOTMM<CT,T1,T2>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

// (x*m%m)*x

template <typename T, typename T1, typename T2>
inline TRQUOTMM<T,T1,T2> operator*(const TRQUOTMM<T,T1,T2>& qmm, T x)
{ return TRQUOTMM<T,T1,T2>(qmm.getX()*x,qmm.getP(),qmm.getM2()); }

template <typename T>
inline TRQUOTMM<CT,T,T> operator*(const TRQUOTMM<CT,T,T>& qmm, CT x)
{ return TRQUOTMM<CT,T,T>(qmm.getX()*x,qmm.getP(),qmm.getM2()); }

template <typename T>
inline TRQUOTMM<CT,T,T> operator*(const TRQUOTMM<CT,T,T>& qmm, CCT x)
{ return TRQUOTMM<CT,T,T>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

template <typename T>
inline TRQUOTMM<CT,T,T> operator*(const TRQUOTMM<CT,T,T>& qmm, VCT x)
{ return TRQUOTMM<CT,T,T>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TRQUOTMM<CT,T1,T2> operator*(const TRQUOTMM<CT,T1,T2>& qmm, T x)
{ return TRQUOTMM<CT,T1,T2>(qmm.getX()*x,qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TRQUOTMM<CT,T1,T2> operator*(const TRQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return TRQUOTMM<CT,T1,T2>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TRQUOTMM<CT,T1,T2> operator*(const TRQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return TRQUOTMM<CT,T1,T2>(qmm.getX()*CT(x),qmm.getP(),qmm.getM2()); }

// (x*m%m)/x

template <typename T, typename T1, typename T2>
inline TRQUOTMM<T,T1,T2> operator/(const TRQUOTMM<T,T1,T2>& qmm, T x)
{ return TRQUOTMM<T,T1,T2>(TMV_Divide(qmm.getX(),x),qmm.getP(),qmm.getM2()); }

template <typename T>
inline TRQUOTMM<CT,T,T> operator/(const TRQUOTMM<CT,T,T>& qmm, CT x)
{ return TRQUOTMM<CT,T,T>(TMV_Divide(qmm.getX(),x),qmm.getP(),qmm.getM2()); }

template <typename T>
inline TRQUOTMM<CT,T,T> operator/(const TRQUOTMM<CT,T,T>& qmm, CCT x)
{ return TRQUOTMM<CT,T,T>(TMV_Divide(qmm.getX(),CT(x)),qmm.getP(),qmm.getM2()); }

template <typename T>
inline TRQUOTMM<CT,T,T> operator/(const TRQUOTMM<CT,T,T>& qmm, VCT x)
{ return TRQUOTMM<CT,T,T>(TMV_Divide(qmm.getX(),CT(x)),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TRQUOTMM<CT,T1,T2> operator/(const TRQUOTMM<CT,T1,T2>& qmm, T x)
{ return TRQUOTMM<CT,T1,T2>(TMV_Divide(qmm.getX(),x),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TRQUOTMM<CT,T1,T2> operator/(const TRQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return TRQUOTMM<CT,T1,T2>(TMV_Divide(qmm.getX(),CT(x)),qmm.getP(),qmm.getM2()); }

template <typename T, typename T1, typename T2>
inline TRQUOTMM<CT,T1,T2> operator/(const TRQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return TRQUOTMM<CT,T1,T2>(TMV_Divide(qmm.getX(),CT(x)),qmm.getP(),qmm.getM2()); }

