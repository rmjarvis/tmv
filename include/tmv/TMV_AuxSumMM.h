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
// (The values for a normal Matrix+Matrix are given)
//
// SUMMM	SumMM
// GENMATRIX1	GenMatrix
// GENMATRIX2	GenMatrix
// PRODXM1	ProdXM
// PRODXM2	ProdXM

#ifndef X1
#define X1
#endif

#ifndef X2
#define X2
#endif

#ifndef X3
#define X3
#endif

#ifndef Y
#define Y
#endif

#ifndef SUMMM_1_1
#define SUMMM_1_1 SUMMM
#define SUMMM_1_m1 SUMMM
#define SUMMM_1_x SUMMM
#define SUMMM_x_1 SUMMM
#define SUMMM_x_m1 SUMMM
#endif

#ifndef GETM1
#define GETM1 .getM()
#endif

#ifndef GETM2
#define GETM2 .getM()
#endif

// m+m

#ifdef INTT1
template <typename T2 Y>
inline SUMMM_1_1<T2,int,T2 X3> operator+(
    const GENMATRIX1<int X1>& m1, const GENMATRIX2<T2 X2>& m2)
{ return SUMMM_1_1<T2,int,T2 X3>(T2(1),m1,T2(1),m2); }
#elif defined(INTT2)
template <typename T1 Y>
inline SUMMM_1_1<T1,T1,int X3> operator+(
    const GENMATRIX1<T1 X1>& m1, const GENMATRIX2<int X2>& m2)
{ return SUMMM_1_1<T1,T1,int X3>(T1(1),m1,T1(1),m2); }
#else
template <typename T Y>
inline SUMMM_1_1<T,T,T X3> operator+(
    const GENMATRIX1<T X1>& m1, const GENMATRIX2<T X2>& m2)
{ return SUMMM_1_1<T,T,T X3>(T(1),m1,T(1),m2); }

template <typename T Y>
inline SUMMM_1_1<CT,CT,T X3> operator+(
    const GENMATRIX1<CT X1>& m1,const GENMATRIX2<T X2>& m2)
{ return SUMMM_1_1<CT,CT,T X3>(CT(1),m1,CT(1),m2); }

template <typename T Y>
inline SUMMM_1_1<CT,T,CT X3> operator+(
    const GENMATRIX1<T X1>& m1,const GENMATRIX2<CT X2>& m2)
{ return SUMMM_1_1<CT,T,CT X3>(CT(1),m1,CT(1),m2); }
#endif

// m-m

#ifdef INTT1
template <typename T2 Y>
inline SUMMM_1_m1<T2,int,T2 X3> operator-(
    const GENMATRIX1<int X1>& m1, const GENMATRIX2<T2 X2>& m2)
{ return SUMMM_1_m1<T2,int,T2 X3>(T2(1),m1,T2(-1),m2); }
#undef INTT1
#elif defined(INTT2)
template <typename T1 Y>
inline SUMMM_1_m1<T1,T1,int X3> operator-(
    const GENMATRIX1<T1 X1>& m1, const GENMATRIX2<int X2>& m2)
{ return SUMMM_1_m1<T1,T1,int X3>(T1(1),m1,T1(-1),m2); }
#undef INTT2
#else
template <typename T Y>
inline SUMMM_1_m1<T,T,T X3> operator-(
    const GENMATRIX1<T X1>& m1, const GENMATRIX2<T X2>& m2)
{ return SUMMM_1_m1<T,T,T X3>(T(1),m1,T(-1),m2); }

template <typename T Y>
inline SUMMM_1_m1<CT,CT,T X3> operator-(
    const GENMATRIX1<CT X1>& m1,const GENMATRIX2<T X2>& m2)
{ return SUMMM_1_m1<CT,CT,T X3>(CT(1),m1,CT(-1),m2); }

template <typename T Y>
inline SUMMM_1_m1<CT,T,CT X3> operator-(
    const GENMATRIX1<T X1>& m1,const GENMATRIX2<CT X2>& m2)
{ return SUMMM_1_m1<CT,T,CT X3>(CT(1),m1,CT(-1),m2); }
#endif

// (x*m)+m

template <typename T, typename T1 Y>
inline SUMMM_x_1<T,T1,T X3> operator+(
    const PRODXM1<T,T1 X1>& pxm, const GENMATRIX2<T X2>& m)
{ return SUMMM_x_1<T,T1,T X3>(pxm.getX(),pxm GETM1,T(1),m); }

template <typename T Y>
inline SUMMM_x_1<CT,T,CT X3> operator+(
    const PRODXM1<T,T X1>& pxm, const GENMATRIX2<CT X2>& m)
{ return SUMMM_x_1<CT,T,CT X3>(CT(pxm.getX()),pxm GETM1,CT(1),m); }

template <typename T, typename T1 Y>
inline SUMMM_x_1<CT,T1,T X3> operator+(
    const PRODXM1<CT,T1 X1>& pxm, const GENMATRIX2<T X2>& m)
{ return SUMMM_x_1<CT,T1,T X3>(pxm.getX(),pxm GETM1,CT(1),m); }

// m+(x*m)

template <typename T, typename T2 Y>
inline SUMMM_1_x<T,T,T2 X3> operator+(
    const GENMATRIX1<T X1>& m, const PRODXM2<T,T2 X2>& pxm)
{ return SUMMM_1_x<T,T,T2 X3>(T(1),m,pxm.getX(),pxm GETM2); }

template <typename T Y>
inline SUMMM_1_x<CT,CT,T X3> operator+(
    const GENMATRIX1<CT X1>& m, const PRODXM2<T,T X2>& pxm)
{ return SUMMM_1_x<CT,CT,T X3>(CT(1),m,CT(pxm.getX()),pxm GETM2); }

template <typename T, typename T2 Y>
inline SUMMM_1_x<CT,T,T2 X3> operator+(
    const GENMATRIX1<T X1>& m, const PRODXM2<CT,T2 X2>& pxm)
{ return SUMMM_1_x<CT,T,T2 X3>(CT(1),m,pxm.getX(),pxm GETM2); }

// (x*m)-m

template <typename T, typename T1 Y>
inline SUMMM_x_m1<T,T1,T X3> operator-(
    const PRODXM1<T,T1 X1>& pxm, const GENMATRIX2<T X2>& m)
{ return SUMMM_x_m1<T,T1,T X3>(pxm.getX(),pxm GETM1,T(-1),m); }

template <typename T Y>
inline SUMMM_x_m1<CT,T,CT X3> operator-(
    const PRODXM1<T,T X1>& pxm, const GENMATRIX2<CT X2>& m)
{ return SUMMM_x_m1<CT,T,CT X3>(CT(pxm.getX()),pxm GETM1,CT(-1),m); }

template <typename T, typename T1 Y>
inline SUMMM_x_m1<CT,T1,T X3> operator-(
    const PRODXM1<CT,T1 X1>& pxm, const GENMATRIX2<T X2>& m)
{ return SUMMM_x_m1<CT,T1,T X3>(pxm.getX(),pxm GETM1,CT(-1),m); }

// m-(x*m)

template <typename T, typename T2 Y>
inline SUMMM_1_x<T,T,T2 X3> operator-(
    const GENMATRIX1<T X1>& m, const PRODXM2<T,T2 X2>& pxm)
{ return SUMMM_1_x<T,T,T2 X3>(T(1),m,-pxm.getX(),pxm GETM2); }

template <typename T Y>
inline SUMMM_1_x<CT,CT,T X3> operator-(
    const GENMATRIX1<CT X1>& m, const PRODXM2<T,T X2>& pxm)
{ return SUMMM_1_x<CT,CT,T X3>(CT(1),m,CT(-pxm.getX()),pxm GETM2); }

template <typename T, typename T2 Y>
inline SUMMM_1_x<CT,T,T2 X3> operator-(
    const GENMATRIX1<T X1>& m, const PRODXM2<CT,T2 X2>& pxm)
{ return SUMMM_1_x<CT,T,T2 X3>(CT(1),m,-pxm.getX(),pxm GETM2); }

// (x*m)+(x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator+(
    const PRODXM1<T,T1 X1>& pxm1, const PRODXM2<T,T2 X2>& pxm2)
{ return SUMMM<T,T1,T2 X3>(pxm1.getX(),pxm1 GETM1,pxm2.getX(),pxm2 GETM2); }

template <typename T, typename T2 Y>
inline SUMMM<CT,T,T2 X3> operator+(
    const PRODXM1<T,T X1>& pxm1, const PRODXM2<CT,T2 X2>& pxm2)
{ return SUMMM<CT,T,T2 X3>(CT(pxm1.getX()),pxm1 GETM1,pxm2.getX(),pxm2 GETM2); }

template <typename T, typename T1 Y>
inline SUMMM<CT,T1,T X3> operator+(
    const PRODXM1<CT,T1 X1>& pxm1, const PRODXM2<T,T X2>& pxm2)
{ return SUMMM<CT,T1,T X3>(pxm1.getX(),pxm1 GETM1,CT(pxm2.getX()),pxm2 GETM2); }

// (x*m)-(x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator-(
    const PRODXM1<T,T1 X1>& pxm1, const PRODXM2<T,T2 X2>& pxm2)
{ return SUMMM<T,T1,T2 X3>(pxm1.getX(),pxm1 GETM1,-pxm2.getX(),pxm2 GETM2); }

template <typename T, typename T2 Y>
inline SUMMM<CT,T,T2 X3> operator-(
    const PRODXM1<T,T X1>& pxm1, const PRODXM2<CT,T2 X2>& pxm2)
{ return SUMMM<CT,T,T2 X3>(CT(pxm1.getX()),pxm1 GETM1,-pxm2.getX(),pxm2 GETM2); }

template <typename T, typename T1 Y>
inline SUMMM<CT,T1,T X3> operator-(
    const PRODXM1<CT,T1 X1>& pxm1, const PRODXM2<T,T X2>& pxm2)
{ return SUMMM<CT,T1,T X3>(pxm1.getX(),pxm1 GETM1,CT(-pxm2.getX()),pxm2 GETM2); }

#undef X1
#undef X2
#undef X3
#undef Y

#undef SUMMM_1_1
#undef SUMMM_1_m1
#undef SUMMM_1_x
#undef SUMMM_x_1
#undef SUMMM_x_m1

#undef GETM1
#undef GETM2
