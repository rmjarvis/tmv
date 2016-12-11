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

// m+m

template <typename T Y>
inline SUMMM<T,T,T X3> operator+(
    const GENMATRIX2<T X2>& m2, const GENMATRIX1<T X1>& m1)
{ return SUMMM<T,T,T X3>(T(1),m1,T(1),m2); }

template <typename T Y>
inline SUMMM<CT,T,CT X3> operator+(
    const GENMATRIX2<CT X2>& m2,const GENMATRIX1<T X1>& m1)
{ return SUMMM<CT,T,CT X3>(CT(1),m1,CT(1),m2); }

template <typename T Y>
inline SUMMM<CT,CT,T X3> operator+(
    const GENMATRIX2<T X2>& m2,const GENMATRIX1<CT X1>& m1)
{ return SUMMM<CT,CT,T X3>(CT(1),m1,CT(1),m2); }

// m-m

template <typename T Y>
inline SUMMM<T,T,T X3> operator-(
    const GENMATRIX2<T X2>& m2, const GENMATRIX1<T X1>& m1)
{ return SUMMM<T,T,T X3>(T(-1),m1,T(1),m2); }

template <typename T Y>
inline SUMMM<CT,T,CT X3> operator-(
    const GENMATRIX2<CT X2>& m2,const GENMATRIX1<T X1>& m1)
{ return SUMMM<CT,T,CT X3>(CT(-1),m1,CT(1),m2); }

template <typename T Y>
inline SUMMM<CT,CT,T X3> operator-(
    const GENMATRIX2<T X2>& m2,const GENMATRIX1<CT X1>& m1)
{ return SUMMM<CT,CT,T X3>(CT(-1),m1,CT(1),m2); }

// (x*m)+m

template <typename T, typename T1 Y>
inline SUMMM<T,T1,T X3> operator+(
    const GENMATRIX2<T X2>& m, const PRODXM1<T,T1 X1>& pxm)
{ return SUMMM<T,T1,T X3>(pxm.getX(),pxm.getM(),T(1),m); }

template <typename T Y>
inline SUMMM<CT,T,CT X3> operator+(
    const GENMATRIX2<CT X2>& m, const PRODXM1<T,T X1>& pxm)
{ return SUMMM<CT,T,CT X3>(CT(pxm.getX()),pxm.getM(),CT(1),m); }

template <typename T, typename T1 Y>
inline SUMMM<CT,T1,T X3> operator+(
    const GENMATRIX2<T X2>& m, const PRODXM1<CT,T1 X1>& pxm)
{ return SUMMM<CT,T1,T X3>(pxm.getX(),pxm.getM(),CT(1),m); }

// m+(x*m)

template <typename T, typename T2 Y>
inline SUMMM<T,T,T2 X3> operator+(
    const PRODXM2<T,T2 X2>& pxm, const GENMATRIX1<T X1>& m)
{ return SUMMM<T,T,T2 X3>(T(1),m,pxm.getX(),pxm.getM()); }

template <typename T Y>
inline SUMMM<CT,CT,T X3> operator+(
    const PRODXM2<T,T X2>& pxm, const GENMATRIX1<CT X1>& m)
{ return SUMMM<CT,CT,T X3>(CT(1),m,CT(pxm.getX()),pxm.getM()); }

template <typename T, typename T2 Y>
inline SUMMM<CT,T,T2 X3> operator+(
    const PRODXM2<CT,T2 X2>& pxm, const GENMATRIX1<T X1>& m)
{ return SUMMM<CT,T,T2 X3>(CT(1),m,pxm.getX(),pxm.getM()); }

// (x*m)-m

template <typename T, typename T1 Y>
inline SUMMM<T,T1,T X3> operator-(
    const GENMATRIX2<T X2>& m, const PRODXM1<T,T1 X1>& pxm)
{ return SUMMM<T,T1,T X3>(-pxm.getX(),pxm.getM(),T(1),m); }

template <typename T Y>
inline SUMMM<CT,T,CT X3> operator-(
    const GENMATRIX2<CT X2>& m, const PRODXM1<T,T X1>& pxm)
{ return SUMMM<CT,T,CT X3>(CT(-pxm.getX()),pxm.getM(),CT(1),m); }

template <typename T, typename T1 Y>
inline SUMMM<CT,T1,T X3> operator-(
    const GENMATRIX2<T X2>& m, const PRODXM1<CT,T1 X1>& pxm)
{ return SUMMM<CT,T1,T X3>(-pxm.getX(),pxm.getM(),CT(1),m); }

// m-(x*m)

template <typename T, typename T2 Y>
inline SUMMM<T,T,T2 X3> operator-(
    const PRODXM2<T,T2 X2>& pxm, const GENMATRIX1<T X1>& m)
{ return SUMMM<T,T,T2 X3>(T(-1),m,pxm.getX(),pxm.getM()); }

template <typename T Y>
inline SUMMM<CT,CT,T X3> operator-(
    const PRODXM2<T,T X2>& pxm, const GENMATRIX1<CT X1>& m)
{ return SUMMM<CT,CT,T X3>(CT(-1),m,CT(pxm.getX()),pxm.getM()); }

template <typename T, typename T2 Y>
inline SUMMM<CT,T,T2 X3> operator-(
    const PRODXM2<CT,T2 X2>& pxm, const GENMATRIX1<T X1>& m)
{ return SUMMM<CT,T,T2 X3>(CT(-1),m,pxm.getX(),pxm.getM()); }

// (x*m)+(x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator+(
    const PRODXM2<T,T2 X2>& pxm2, const PRODXM1<T,T1 X1>& pxm1)
{ return SUMMM<T,T1,T2 X3>(pxm1.getX(),pxm1.getM(),pxm2.getX(),pxm2.getM()); }

template <typename T, typename T1 Y>
inline SUMMM<CT,T1,T X3> operator+(
    const PRODXM2<T,T X2>& pxm2, const PRODXM1<CT,T1 X1>& pxm1)
{
    return SUMMM<CT,T1,T X3>(
        pxm1.getX(),pxm1.getM(),CT(pxm2.getX()),pxm2.getM());
}

template <typename T, typename T2 Y>
inline SUMMM<CT,T,T2 X3> operator+(
    const PRODXM2<CT,T2 X2>& pxm2, const PRODXM1<T,T X1>& pxm1)
{
    return SUMMM<CT,T,T2 X3>(
        CT(pxm1.getX()),pxm1.getM(),pxm2.getX(),pxm2.getM());
}


// (x*m)-(x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator-(
    const PRODXM2<T,T2 X2>& pxm2, const PRODXM1<T,T1 X1>& pxm1)
{ return SUMMM<T,T1,T2 X3>(-pxm1.getX(),pxm1.getM(),pxm2.getX(),pxm2.getM()); }

template <typename T, typename T1 Y>
inline SUMMM<CT,T1,T X3> operator-(
    const PRODXM2<T,T X2>& pxm2, const PRODXM1<CT,T1 X1>& pxm1)
{
    return SUMMM<CT,T1,T X3>(
        -pxm1.getX(),pxm1.getM(),CT(pxm2.getX()),pxm2.getM());
}

template <typename T, typename T2 Y>
inline SUMMM<CT,T,T2 X3> operator-(
    const PRODXM2<CT,T2 X2>& pxm2, const PRODXM1<T,T X1>& pxm1)
{
    return SUMMM<CT,T,T2 X3>(
        CT(-pxm1.getX()),pxm1.getM(),pxm2.getX(),pxm2.getM());
}

#undef X1
#undef X2
#undef X3
#undef Y
