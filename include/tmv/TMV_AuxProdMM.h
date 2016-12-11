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
// (The given definition is for a regular Matrix*Matrix.  Modify as
// appropriate for the various other matrices.)
//
// #define GENMATRIX1 GenMatrix
// #define GENMATRIX2 GenMatrix
// #define PRODMM ProdMM
// #define PRODXM1 ProdXM
// #define PRODXM2 ProdXM

#ifndef OP
#define OP operator*
#endif

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

#ifndef PRODMM_1
#define PRODMM_1 PRODMM
#endif

#ifndef GETM1
#define GETM1 .getM()
#endif

#ifndef GETM2
#define GETM2 .getM()
#endif

// m*m
#ifdef INTT1
template <typename T2 Y>
inline PRODMM_1<T2,int,T2 X3> OP(
    const GENMATRIX1<int X1>& m1, const GENMATRIX2<T2 X2>& m2)
{ return PRODMM_1<T2,int,T2 X3>(T2(1),m1,m2); }
#undef INTT1
#elif defined(INTT2)
template <typename T1 Y>
inline PRODMM_1<T1,T1,int X3> OP(
    const GENMATRIX1<T1 X1>& m1, const GENMATRIX2<int X2>& m2)
{ return PRODMM_1<T1,T1,int X3>(T1(1),m1,m2); }
#undef INTT2
#else
template <typename T Y>
inline PRODMM_1<T,T,T X3> OP(
    const GENMATRIX1<T X1>& m1, const GENMATRIX2<T X2>& m2)
{ return PRODMM_1<T,T,T X3>(T(1),m1,m2); }

template <typename T Y>
inline PRODMM_1<CT,CT,T X3> OP(
    const GENMATRIX1<CT X1>& m1, const GENMATRIX2<T X2>& m2)
{ return PRODMM_1<CT,CT,T X3>(CT(1),m1,m2); }

template <typename T Y>
inline PRODMM_1<CT,T,CT X3> OP(
    const GENMATRIX1<T X1>& m1, const GENMATRIX2<CT X2>& m2)
{ return PRODMM_1<CT,T,CT X3>(CT(1),m1,m2); }
#endif

// m*(x*m)
template <typename T, typename T2 Y>
inline PRODMM<T,T,T2 X3> OP(
    const GENMATRIX1<T X1>& m, const PRODXM2<T,T2 X2>& pxm)
{ return PRODMM<T,T,T2 X3>(pxm.getX(),m,pxm GETM2); }

template <typename T Y>
inline PRODMM<CT,CT,T X3> OP(
    const GENMATRIX1<CT X1>& m, const PRODXM2<T,T X2>& pxm)
{ return PRODMM<CT,CT,T X3>(CT(pxm.getX()),m,pxm GETM2); }

template <typename T, typename T2 Y>
inline PRODMM<CT,T,T2 X3> OP(
    const GENMATRIX1<T X1>& m, const PRODXM2<CT,T2 X2>& pxm)
{ return PRODMM<CT,T,T2 X3>(pxm.getX(),m,pxm GETM2); }

// (x*m)*m
template <typename T, typename T1 Y>
inline PRODMM<T,T1,T X3> OP(
    const PRODXM1<T,T1 X1>& pxm, const GENMATRIX2<T X2>& m)
{ return PRODMM<T,T1,T X3>(pxm.getX(),pxm GETM1,m); }

template <typename T Y>
inline PRODMM<CT,T,CT X3> OP(
    const PRODXM1<T,T X1>& pxm, const GENMATRIX2<CT X2>& m)
{ return PRODMM<CT,T,CT X3>(CT(pxm.getX()),pxm GETM1,m); }

template <typename T, typename T1 Y>
inline PRODMM<CT,T1,T X3> OP(
    const PRODXM1<CT,T1 X1>& pxm, const GENMATRIX2<T X2>& m)
{ return PRODMM<CT,T1,T X3>(pxm.getX(),pxm GETM1,m); }

// (x*m)*(x*m)
template <typename T, typename T1, typename T2 Y>
inline PRODMM<T,T1,T2 X3> OP(
    const PRODXM1<T,T1 X1>& pxm1, const PRODXM2<T,T2 X2>& pxm2)
{
    return PRODMM<T,T1,T2 X3>(
        pxm1.getX()*pxm2.getX(),pxm1 GETM1,pxm2 GETM2);
}

template <typename T, typename T1 Y>
inline PRODMM<CT,T1,T X3> OP(
    const PRODXM1<CT,T1 X1>& pxm1, const PRODXM2<T,T X2>& pxm2)
{
    return PRODMM<CT,T1,T X3>(
        pxm1.getX()*pxm2.getX(),pxm1 GETM1,pxm2 GETM2);
}

template <typename T, typename T2 Y>
inline PRODMM<CT,T,T2 X3> OP(
    const PRODXM1<T,T X1>& pxm1, const PRODXM2<CT,T2 X2>& pxm2)
{
    return PRODMM<CT,T,T2 X3>(
        pxm1.getX()*pxm2.getX(),pxm1 GETM1,pxm2 GETM2);
}

#undef OP
#undef X1
#undef X2
#undef X3
#undef Y
#undef PRODMM_1
#undef GETM1
#undef GETM2
