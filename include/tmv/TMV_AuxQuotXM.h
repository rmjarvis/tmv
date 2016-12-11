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

// x/m

template <typename T Y>
inline QUOTXM<T,T X> operator/(T x, const GENMATRIX<T X>& m)
{ return QUOTXM<T,T X>(x,m); }

template <typename T Y>
inline QUOTXM<CT,T X> operator/(CT x, const GENMATRIX<T X>& m)
{ return QUOTXM<CT,T X>(x,m); }

template <typename T Y>
inline QUOTXM<CT,T X> operator/(CCT x, const GENMATRIX<T X>& m)
{ return QUOTXM<CT,T X>(CT(x),m); }

template <typename T Y>
inline QUOTXM<CT,T X> operator/(VCT x, const GENMATRIX<T X>& m)
{ return QUOTXM<CT,T X>(CT(x),m); }

template <typename T Y>
inline QUOTXM<CT,CT X> operator/(T x, const GENMATRIX<CT X>& m)
{ return QUOTXM<CT,CT X>(CT(x),m); }

template <typename T Y>
inline QUOTXM<CT,CT X> operator/(CCT x, const GENMATRIX<CT X>& m)
{ return QUOTXM<CT,CT X>(CT(x),m); }

template <typename T Y>
inline QUOTXM<CT,CT X> operator/(VCT x, const GENMATRIX<CT X>& m)
{ return QUOTXM<CT,CT X>(CT(x),m); }

// x%m

template <typename T Y>
inline QUOTXM<T,T X> operator%(T x, const GENMATRIX<T X>& m)
{ return QUOTXM<T,T X>(x,m); }

template <typename T Y>
inline QUOTXM<CT,T X> operator%(CT x, const GENMATRIX<T X>& m)
{ return QUOTXM<CT,T X>(x,m); }

template <typename T Y>
inline QUOTXM<CT,T X> operator%(CCT x, const GENMATRIX<T X>& m)
{ return QUOTXM<CT,T X>(CT(x),m); }

template <typename T Y>
inline QUOTXM<CT,T X> operator%(VCT x, const GENMATRIX<T X>& m)
{ return QUOTXM<CT,T X>(CT(x),m); }

template <typename T Y>
inline QUOTXM<CT,CT X> operator%(T x, const GENMATRIX<CT X>& m)
{ return QUOTXM<CT,CT X>(CT(x),m); }

template <typename T Y>
inline QUOTXM<CT,CT X> operator%(CCT x, const GENMATRIX<CT X>& m)
{ return QUOTXM<CT,CT X>(CT(x),m); }

template <typename T Y>
inline QUOTXM<CT,CT X> operator%(VCT x, const GENMATRIX<CT X>& m)
{ return QUOTXM<CT,CT X>(CT(x),m); }


#undef X
#undef Y
#undef GETM
