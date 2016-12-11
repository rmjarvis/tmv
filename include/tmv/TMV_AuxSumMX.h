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
// SUMMX	SumMX
// PRODXM	ProdXM
// GENMATRIX	GenMatrix

#ifndef X1
#define X1
#endif

#ifndef X2
#define X2
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

#ifndef SUMMX_1
#define SUMMX_1 SUMMX
#define NO_1
#endif

// m+x
#ifdef INTT
inline SUMMX_1<int,int> operator+(const GENMATRIX<int>& m1, int x2)
{ return SUMMX_1<int,int>(int(1),m1,x2); }

inline SUMMX_1<float,int> operator+(const GENMATRIX<int>& m1, float x2)
{ return SUMMX_1<float,int>(float(1),m1,x2); }

inline SUMMX_1<double,int> operator+(const GENMATRIX<int>& m1, double x2)
{ return SUMMX_1<double,int>(double(1),m1,x2); }

inline SUMMX_1<long double,int> operator+(
    const GENMATRIX<int>& m1, long double x2)
{ return SUMMX_1<long double,int>((long double)(1),m1,x2); }

template <typename T Y>
inline SUMMX_1<CT,int X2> operator+(const GENMATRIX<int X1>& m1, CCT x2)
{ return SUMMX_1<CT,int X2>(CT(1),m1,CT(x2)); }

template <typename T Y>
inline SUMMX_1<CT,int X2> operator+(const GENMATRIX<int X1>& m1, VCT x2)
{ return SUMMX_1<CT,int X2>(CT(1),m1,CT(x2)); }
#else
template <typename T Y>
inline SUMMX_1<T,T X2> operator+(const GENMATRIX<T X1>& m1, T x2)
{ return SUMMX_1<T,T X2>(T(1),m1,x2); }

template <typename T Y>
inline SUMMX_1<CT,T X2> operator+(const GENMATRIX<T X1>& m1, CT x2)
{ return SUMMX_1<CT,T X2>(CT(1),m1,x2); }

template <typename T Y>
inline SUMMX_1<CT,T X2> operator+(const GENMATRIX<T X1>& m1, CCT x2)
{ return SUMMX_1<CT,T X2>(CT(1),m1,CT(x2)); }

template <typename T Y>
inline SUMMX_1<CT,T X2> operator+(const GENMATRIX<T X1>& m1, VCT x2)
{ return SUMMX_1<CT,T X2>(CT(1),m1,CT(x2)); }

template <typename T Y>
inline SUMMX_1<CT,CT X2> operator+(const GENMATRIX<CT X1>& m1, T x2)
{ return SUMMX_1<CT,CT X2>(CT(1),m1,x2); }

template <typename T Y>
inline SUMMX_1<CT,CT X2> operator+(const GENMATRIX<CT X1>& m1, CCT x2)
{ return SUMMX_1<CT,CT X2>(CT(1),m1,CT(x2)); }

template <typename T Y>
inline SUMMX_1<CT,CT X2> operator+(const GENMATRIX<CT X1>& m1, VCT x2)
{ return SUMMX_1<CT,CT X2>(CT(1),m1,CT(x2)); }
#endif

// x+m
#ifdef INTT
inline SUMMX_1<int,int> operator+(int x1, const GENMATRIX<int>& m2)
{ return SUMMX_1<int,int>(int(1),m2,x1); }

inline SUMMX_1<float,int> operator+(float x1, const GENMATRIX<int>& m2)
{ return SUMMX_1<float,int>(float(1),m2,x1); }

inline SUMMX_1<double,int> operator+(double x1, const GENMATRIX<int>& m2)
{ return SUMMX_1<double,int>(double(1),m2,x1); }

inline SUMMX_1<long double,int> operator+(
    long double x1, const GENMATRIX<int>& m2)
{ return SUMMX_1<long double,int>((long double)(1),m2,x1); }

template <typename T Y>
inline SUMMX_1<CT,int X2> operator+(CCT x1, const GENMATRIX<int X1>& m2)
{ return SUMMX_1<CT,int X2>(CT(1),m2,CT(x1)); }

template <typename T Y>
inline SUMMX_1<CT,int X2> operator+(VCT x1, const GENMATRIX<int X1>& m2)
{ return SUMMX_1<CT,int X2>(CT(1),m2,CT(x1)); }
#else
template <typename T Y>
inline SUMMX_1<T,T X2> operator+(T x1, const GENMATRIX<T X1>& m2)
{ return SUMMX_1<T,T X2>(T(1),m2,x1); }

template <typename T Y>
inline SUMMX_1<CT,T X2> operator+(CT x1, const GENMATRIX<T X1>& m2)
{ return SUMMX_1<CT,T X2>(CT(1),m2,x1); }

template <typename T Y>
inline SUMMX_1<CT,T X2> operator+(CCT x1, const GENMATRIX<T X1>& m2)
{ return SUMMX_1<CT,T X2>(CT(1),m2,CT(x1)); }

template <typename T Y>
inline SUMMX_1<CT,T X2> operator+(VCT x1, const GENMATRIX<T X1>& m2)
{ return SUMMX_1<CT,T X2>(CT(1),m2,CT(x1)); }

template <typename T Y>
inline SUMMX_1<CT,CT X2> operator+(T x1, const GENMATRIX<CT X1>& m2)
{ return SUMMX_1<CT,CT X2>(CT(1),m2,x1); }

template <typename T Y>
inline SUMMX_1<CT,CT X2> operator+(CCT x1, const GENMATRIX<CT X1>& m2)
{ return SUMMX_1<CT,CT X2>(CT(1),m2,CT(x1)); }

template <typename T Y>
inline SUMMX_1<CT,CT X2> operator+(VCT x1, const GENMATRIX<CT X1>& m2)
{ return SUMMX_1<CT,CT X2>(CT(1),m2,CT(x1)); }
#endif

// m-x
#ifdef INTT
inline SUMMX_1<int,int> operator-(const GENMATRIX<int>& m1, int x2)
{ return SUMMX_1<int,int>(int(1),m1,-x2); }

inline SUMMX_1<float,int> operator-(const GENMATRIX<int>& m1, float x2)
{ return SUMMX_1<float,int>(float(1),m1,-x2); }

inline SUMMX_1<double,int> operator-(const GENMATRIX<int>& m1, double x2)
{ return SUMMX_1<double,int>(double(1),m1,-x2); }

inline SUMMX_1<long double,int> operator-(
    const GENMATRIX<int>& m1, long double x2)
{ return SUMMX_1<long double,int>((long double)(1),m1,-x2); }

template <typename T Y>
inline SUMMX_1<CT,int X2> operator-(const GENMATRIX<int X1>& m1, CCT x2)
{ return SUMMX_1<CT,int X2>(CT(1),m1,-CT(x2)); }

template <typename T Y>
inline SUMMX_1<CT,int X2> operator-(const GENMATRIX<int X1>& m1, VCT x2)
{ return SUMMX_1<CT,int X2>(CT(1),m1,-CT(x2)); }
#else
template <typename T Y>
inline SUMMX_1<T,T X2> operator-(const GENMATRIX<T X1>& m1, T x2)
{ return SUMMX_1<T,T X2>(T(1),m1,-x2); }

template <typename T Y>
inline SUMMX_1<CT,T X2> operator-(const GENMATRIX<T X1>& m1, CT x2)
{ return SUMMX_1<CT,T X2>(CT(1),m1,-x2); }

template <typename T Y>
inline SUMMX_1<CT,T X2> operator-(const GENMATRIX<T X1>& m1, CCT x2)
{ return SUMMX_1<CT,T X2>(CT(1),m1,-CT(x2)); }

template <typename T Y>
inline SUMMX_1<CT,T X2> operator-(const GENMATRIX<T X1>& m1, VCT x2)
{ return SUMMX_1<CT,T X2>(CT(1),m1,-CT(x2)); }

template <typename T Y>
inline SUMMX_1<CT,CT X2> operator-(const GENMATRIX<CT X1>& m1, T x2)
{ return SUMMX_1<CT,CT X2>(CT(1),m1,-x2); }

template <typename T Y>
inline SUMMX_1<CT,CT X2> operator-(const GENMATRIX<CT X1>& m1, CCT x2)
{ return SUMMX_1<CT,CT X2>(CT(1),m1,-CT(x2)); }

template <typename T Y>
inline SUMMX_1<CT,CT X2> operator-(const GENMATRIX<CT X1>& m1, VCT x2)
{ return SUMMX_1<CT,CT X2>(CT(1),m1,-CT(x2)); }
#endif

// x-m
#ifdef INTT
inline SUMMX<int,int> operator-(int x1, const GENMATRIX<int>& m2)
{ return SUMMX<int,int>(int(-1),m2,x1); }

inline SUMMX<float,int> operator-(float x1, const GENMATRIX<int>& m2)
{ return SUMMX<float,int>(float(-1),m2,x1); }

inline SUMMX<double,int> operator-(double x1, const GENMATRIX<int>& m2)
{ return SUMMX<double,int>(double(-1),m2,x1); }

inline SUMMX<long double,int> operator-(
    long double x1, const GENMATRIX<int>& m2)
{ return SUMMX<long double,int>((long double)(-1),m2,x1); }

template <typename T Y>
inline SUMMX<CT,int X2> operator-(CCT x1, const GENMATRIX<int X1>& m2)
{ return SUMMX<CT,int X2>(CT(-1),m2,CT(x1)); }

template <typename T Y>
inline SUMMX<CT,int X2> operator-(VCT x1, const GENMATRIX<int X1>& m2)
{ return SUMMX<CT,int X2>(CT(-1),m2,CT(x1)); }
#undef INTT
#else
template <typename T Y>
inline SUMMX<T,T X2> operator-(T x1, const GENMATRIX<T X1>& m2)
{ return SUMMX<T,T X2>(T(-1),m2,x1); }

template <typename T Y>
inline SUMMX<CT,T X2> operator-(CT x1, const GENMATRIX<T X1>& m2)
{ return SUMMX<CT,T X2>(CT(-1),m2,x1); }

template <typename T Y>
inline SUMMX<CT,T X2> operator-(CCT x1, const GENMATRIX<T X1>& m2)
{ return SUMMX<CT,T X2>(CT(-1),m2,CT(x1)); }

template <typename T Y>
inline SUMMX<CT,T X2> operator-(VCT x1, const GENMATRIX<T X1>& m2)
{ return SUMMX<CT,T X2>(CT(-1),m2,CT(x1)); }

template <typename T Y>
inline SUMMX<CT,CT X2> operator-(T x1, const GENMATRIX<CT X1>& m2)
{ return SUMMX<CT,CT X2>(CT(-1),m2,x1); }

template <typename T Y>
inline SUMMX<CT,CT X2> operator-(CCT x1, const GENMATRIX<CT X1>& m2)
{ return SUMMX<CT,CT X2>(CT(-1),m2,CT(x1)); }

template <typename T Y>
inline SUMMX<CT,CT X2> operator-(VCT x1, const GENMATRIX<CT X1>& m2)
{ return SUMMX<CT,CT X2>(CT(-1),m2,CT(x1)); }
#endif

// -(x*m+x)

template <typename T, typename T1 Y>
inline SUMMX<T,T1 X2> operator-(const SUMMX<T,T1 X2>& smx)
{ return SUMMX<T,T1 X2>(-smx.getX1(),smx.getM(),-smx.getX2); }

// x*(x*m+x)

template <typename T, typename T1 Y>
inline SUMMX<T,T1 X2> operator*(const T x, const SUMMX<T,T1 X2>& smx)
{ return SUMMX<T,T1 X2>(smx.getX1()*x,smx.getM(),smx.getX2()*x); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator*(const T x, const SUMMX<CT,T1 X2>& smx)
{ return SUMMX<CT,T1 X2>(smx.getX1()*x,smx.getM(),smx.getX2()*x); }

template <typename T Y>
inline SUMMX<CT,T X2> operator*(const CT x, const SUMMX<T,T X2>& smx)
{ return SUMMX<CT,T X2>(smx.getX1()*x,smx.getM(),smx.getX2()*x); }

template <typename T Y>
inline SUMMX<CT,T X2> operator*(const CCT x, const SUMMX<T,T X2>& smx)
{ return SUMMX<CT,T X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <typename T Y>
inline SUMMX<CT,T X2> operator*(const VCT x, const SUMMX<T,T X2>& smx)
{ return SUMMX<CT,T X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator*(const CCT x, const SUMMX<CT,T1 X2>& smx)
{ return SUMMX<CT,T1 X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator*(const VCT x, const SUMMX<CT,T1 X2>& smx)
{ return SUMMX<CT,T1 X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

// (x*m+x)*x

template <typename T, typename T1 Y>
inline SUMMX<T,T1 X2> operator*(const SUMMX<T,T1 X2>& smx, const T x)
{ return SUMMX<T,T1 X2>(smx.getX1()*x,smx.getM(),smx.getX2()*x); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator*(const SUMMX<CT,T1 X2>& smx, const T x)
{ return SUMMX<CT,T1 X2>(smx.getX1()*x,smx.getM(),smx.getX2()*x); }

template <typename T Y>
inline SUMMX<CT,T X2> operator*(const SUMMX<T,T X2>& smx, const CT x)
{ return SUMMX<CT,T X2>(smx.getX1()*x,smx.getM(),smx.getX2()*x); }

template <typename T Y>
inline SUMMX<CT,T X2> operator*(const SUMMX<T,T X2>& smx, const CCT x)
{ return SUMMX<CT,T X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <typename T Y>
inline SUMMX<CT,T X2> operator*(const SUMMX<T,T X2>& smx, const VCT x)
{ return SUMMX<CT,T X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator*(const SUMMX<CT,T1 X2>& smx, const CCT x)
{ return SUMMX<CT,T1 X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator*(const SUMMX<CT,T1 X2>& smx, const VCT x)
{ return SUMMX<CT,T1 X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

// (x*m+x)/x

template <typename T, typename T1 Y>
inline SUMMX<T,T1 X2> operator/(const SUMMX<T,T1 X2>& smx, const T x)
{ return SUMMX<T,T1 X2>(TMV_Divide(smx.getX1(),x),smx.getM(),TMV_Divide(smx.getX2(),x)); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator/(const SUMMX<CT,T1 X2>& smx, const T x)
{ return SUMMX<CT,T1 X2>(TMV_Divide(smx.getX1(),x),smx.getM(),TMV_Divide(smx.getX2(),x)); }

template <typename T Y>
inline SUMMX<CT,T X2> operator/(const SUMMX<T,T X2>& smx, const CT x)
{ return SUMMX<CT,T X2>(TMV_Divide(smx.getX1(),x),smx.getM(),TMV_Divide(smx.getX2(),x)); }

template <typename T Y>
inline SUMMX<CT,T X2> operator/(const SUMMX<T,T X2>& smx, const CCT x)
{ return SUMMX<CT,T X2>(TMV_Divide(smx.getX1(),CT(x)),smx.getM(),TMV_Divide(smx.getX2(),CT(x))); }

template <typename T Y>
inline SUMMX<CT,T X2> operator/(const SUMMX<T,T X2>& smx, const VCT x)
{ return SUMMX<CT,T X2>(TMV_Divide(smx.getX1(),CT(x)),smx.getM(),TMV_Divide(smx.getX2(),CT(x))); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator/(const SUMMX<CT,T1 X2>& smx, const CCT x)
{ return SUMMX<CT,T1 X2>(TMV_Divide(smx.getX1(),CT(x)),smx.getM(),TMV_Divide(smx.getX2(),CT(x))); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator/(const SUMMX<CT,T1 X2>& smx, const VCT x)
{ return SUMMX<CT,T1 X2>(TMV_Divide(smx.getX1(),CT(x)),smx.getM(),TMV_Divide(smx.getX2(),CT(x))); }

// x+(x*m)

template <typename T, typename T2 Y>
inline SUMMX<T,T2 X2> operator+(const T x3, const PRODXM<T,T2 X1>& pxm)
{ return SUMMX<T,T2 X2>(pxm.getX1(),pxm.getM2(),x3); }

template <typename T, typename T2 Y>
inline SUMMX<CT,T2 X2> operator+(const T x3, const PRODXM<CT,T2 X1>& pxm)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),x3); }

template <typename T Y>
inline SUMMX<CT,T X2> operator+(const CT x3, const PRODXM<T,T X1>& pxm)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),x3); }

template <typename T Y>
inline SUMMX<CT,T X2> operator+(const CCT x3, const PRODXM<T,T X1>& pxm)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

template <typename T Y>
inline SUMMX<CT,T X2> operator+(const VCT x3, const PRODXM<T,T X1>& pxm)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

template <typename T, typename T2 Y>
inline SUMMX<CT,T2 X2> operator+(const CCT x3, const PRODXM<CT,T2 X1>& pxm)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

template <typename T, typename T2 Y>
inline SUMMX<CT,T2 X2> operator+(const VCT x3, const PRODXM<CT,T2 X1>& pxm)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

// (x*m)+x

template <typename T, typename T2 Y>
inline SUMMX<T,T2 X2> operator+(const PRODXM<T,T2 X1>& pxm, const T x3)
{ return SUMMX<T,T2 X2>(pxm.getX1(),pxm.getM2(),x3); }

template <typename T, typename T2 Y>
inline SUMMX<CT,T2 X2> operator+(const PRODXM<CT,T2 X1>& pxm, const T x3)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),x3); }

template <typename T Y>
inline SUMMX<CT,T X2> operator+(const PRODXM<T,T X1>& pxm, const CT x3)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),x3); }

template <typename T Y>
inline SUMMX<CT,T X2> operator+(const PRODXM<T,T X1>& pxm, const CCT x3)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

template <typename T Y>
inline SUMMX<CT,T X2> operator+(const PRODXM<T,T X1>& pxm, const VCT x3)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

template <typename T, typename T2 Y>
inline SUMMX<CT,T2 X2> operator+(const PRODXM<CT,T2 X1>& pxm, const CCT x3)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

template <typename T, typename T2 Y>
inline SUMMX<CT,T2 X2> operator+(const PRODXM<CT,T2 X1>& pxm, const VCT x3)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

// x-(x*m)

template <typename T, typename T2 Y>
inline SUMMX<T,T2 X2> operator-(const T x3, const PRODXM<T,T2 X1>& pxm)
{ return SUMMX<T,T2 X2>(-pxm.getX1(),pxm.getM2(),x3); }

template <typename T, typename T2 Y>
inline SUMMX<CT,T2 X2> operator-(const T x3, const PRODXM<CT,T2 X1>& pxm)
{ return SUMMX<CT,T2 X2>(-pxm.getX1(),pxm.getM2(),x3); }

template <typename T Y>
inline SUMMX<CT,T X2> operator-(const CT x3, const PRODXM<T,T X1>& pxm)
{ return SUMMX<CT,T X2>(-pxm.getX1(),pxm.getM2(),x3); }

template <typename T Y>
inline SUMMX<CT,T X2> operator-(const CCT x3, const PRODXM<T,T X1>& pxm)
{ return SUMMX<CT,T X2>(-pxm.getX1(),pxm.getM2(),CT(x3)); }

template <typename T Y>
inline SUMMX<CT,T X2> operator-(const VCT x3, const PRODXM<T,T X1>& pxm)
{ return SUMMX<CT,T X2>(-pxm.getX1(),pxm.getM2(),CT(x3)); }

template <typename T, typename T2 Y>
inline SUMMX<CT,T2 X2> operator-(const CCT x3, const PRODXM<CT,T2 X1>& pxm)
{ return SUMMX<CT,T2 X2>(-pxm.getX1(),pxm.getM2(),CT(x3)); }

template <typename T, typename T2 Y>
inline SUMMX<CT,T2 X2> operator-(const VCT x3, const PRODXM<CT,T2 X1>& pxm)
{ return SUMMX<CT,T2 X2>(-pxm.getX1(),pxm.getM2(),CT(x3)); }

// (x*m)-x

template <typename T, typename T2 Y>
inline SUMMX<T,T2 X2> operator-(const PRODXM<T,T2 X1>& pxm, const T x3)
{ return SUMMX<T,T2 X2>(pxm.getX1(),pxm.getM2(),-x3); }

template <typename T, typename T2 Y>
inline SUMMX<CT,T2 X2> operator-(const PRODXM<CT,T2 X1>& pxm, const T x3)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),-x3); }

template <typename T Y>
inline SUMMX<CT,T X2> operator-(const PRODXM<T,T X1>& pxm, const CT x3)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),-x3); }

template <typename T Y>
inline SUMMX<CT,T X2> operator-(const PRODXM<T,T X1>& pxm, const CCT x3)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),-CT(x3)); }

template <typename T Y>
inline SUMMX<CT,T X2> operator-(const PRODXM<T,T X1>& pxm, const VCT x3)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),-CT(x3)); }

template <typename T, typename T2 Y>
inline SUMMX<CT,T2 X2> operator-(const PRODXM<CT,T2 X1>& pxm, const CCT x3)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),-CT(x3)); }

template <typename T, typename T2 Y>
inline SUMMX<CT,T2 X2> operator-(const PRODXM<CT,T2 X1>& pxm, const VCT x3)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),-CT(x3)); }

#ifndef NO_1
// -(x*m+x)

template <typename T, typename T1 Y>
inline SUMMX<T,T1 X2> operator-(const SUMMX_1<T,T1 X2>& smx)
{ return SUMMX<T,T1 X2>(T(-1),smx.getM(),-smx.getX2); }

// x*(x*m+x)

template <typename T, typename T1 Y>
inline SUMMX<T,T1 X2> operator*(const T x, const SUMMX_1<T,T1 X2>& smx)
{ return SUMMX<T,T1 X2>(x,smx.getM(),smx.getX2()*x); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator*(const T x, const SUMMX_1<CT,T1 X2>& smx)
{ return SUMMX<CT,T1 X2>(x,smx.getM(),smx.getX2()*x); }

template <typename T Y>
inline SUMMX<CT,T X2> operator*(const CT x, const SUMMX_1<T,T X2>& smx)
{ return SUMMX<CT,T X2>(x,smx.getM(),smx.getX2()*x); }

template <typename T Y>
inline SUMMX<CT,T X2> operator*(const CCT x, const SUMMX_1<T,T X2>& smx)
{ return SUMMX<CT,T X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <typename T Y>
inline SUMMX<CT,T X2> operator*(const VCT x, const SUMMX_1<T,T X2>& smx)
{ return SUMMX<CT,T X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator*(const CCT x, const SUMMX_1<CT,T1 X2>& smx)
{ return SUMMX<CT,T1 X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator*(const VCT x, const SUMMX_1<CT,T1 X2>& smx)
{ return SUMMX<CT,T1 X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

// (x*m+x)*x

template <typename T, typename T1 Y>
inline SUMMX<T,T1 X2> operator*(const SUMMX_1<T,T1 X2>& smx, const T x)
{ return SUMMX<T,T1 X2>(x,smx.getM(),smx.getX2()*x); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator*(const SUMMX_1<CT,T1 X2>& smx, const T x)
{ return SUMMX<CT,T1 X2>(x,smx.getM(),smx.getX2()*x); }

template <typename T Y>
inline SUMMX<CT,T X2> operator*(const SUMMX_1<T,T X2>& smx, const CT x)
{ return SUMMX<CT,T X2>(x,smx.getM(),smx.getX2()*x); }

template <typename T Y>
inline SUMMX<CT,T X2> operator*(const SUMMX_1<T,T X2>& smx, const CCT x)
{ return SUMMX<CT,T X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <typename T Y>
inline SUMMX<CT,T X2> operator*(const SUMMX_1<T,T X2>& smx, const VCT x)
{ return SUMMX<CT,T X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator*(const SUMMX_1<CT,T1 X2>& smx, const CCT x)
{ return SUMMX<CT,T1 X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator*(const SUMMX_1<CT,T1 X2>& smx, const VCT x)
{ return SUMMX<CT,T1 X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

// (x*m+x)/x

template <typename T, typename T1 Y>
inline SUMMX<T,T1 X2> operator/(const SUMMX_1<T,T1 X2>& smx, const T x)
{ return SUMMX<T,T1 X2>(TMV_InverseOf(x),smx.getM(),TMV_Divide(smx.getX2(),x)); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator/(const SUMMX_1<CT,T1 X2>& smx, const T x)
{ return SUMMX<CT,T1 X2>(TMV_InverseOf(x),smx.getM(),TMV_Divide(smx.getX2(),x)); }

template <typename T Y>
inline SUMMX<CT,T X2> operator/(const SUMMX_1<T,T X2>& smx, const CT x)
{ return SUMMX<CT,T X2>(TMV_InverseOf(x),smx.getM(),TMV_Divide(smx.getX2(),x)); }

template <typename T Y>
inline SUMMX<CT,T X2> operator/(const SUMMX_1<T,T X2>& smx, const CCT x)
{ return SUMMX<CT,T X2>(TMV_InverseOf(CT(x)),smx.getM(),TMV_Divide(smx.getX2(),CT(x))); }

template <typename T Y>
inline SUMMX<CT,T X2> operator/(const SUMMX_1<T,T X2>& smx, const VCT x)
{ return SUMMX<CT,T X2>(TMV_InverseOf(CT(x)),smx.getM(),TMV_Divide(smx.getX2(),CT(x))); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator/(const SUMMX_1<CT,T1 X2>& smx, const CCT x)
{ return SUMMX<CT,T1 X2>(TMV_InverseOf(CT(x)),smx.getM(),TMV_Divide(smx.getX2(),CT(x))); }

template <typename T, typename T1 Y>
inline SUMMX<CT,T1 X2> operator/(const SUMMX_1<CT,T1 X2>& smx, const VCT x)
{ return SUMMX<CT,T1 X2>(TMV_InverseOf(CT(x)),smx.getM(),TMV_Divide(smx.getX2(),CT(x))); }
#endif

#undef X1
#undef X2
#undef Y
#undef SUMMX_1
#ifdef NO_1
#undef NO_1
#endif
