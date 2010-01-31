///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
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

template <class T Y> 
inline SUMMX_1<T,T X2> operator+(const GENMATRIX<T X1>& m1, T x2)
{ return SUMMX_1<T,T X2>(T(1),m1,x2); }

template <class T Y> 
inline SUMMX_1<CT,T X2> operator+(const GENMATRIX<T X1>& m1, CT x2)
{ return SUMMX_1<CT,T X2>(CT(1),m1,x2); }

template <class T Y> 
inline SUMMX_1<CT,T X2> operator+(const GENMATRIX<T X1>& m1, CCT x2)
{ return SUMMX_1<CT,T X2>(CT(1),m1,CT(x2)); }

template <class T Y> 
inline SUMMX_1<CT,T X2> operator+(const GENMATRIX<T X1>& m1, VCT x2)
{ return SUMMX_1<CT,T X2>(CT(1),m1,CT(x2)); }

template <class T Y> 
inline SUMMX_1<CT,CT X2> operator+(const GENMATRIX<CT X1>& m1, T x2)
{ return SUMMX_1<CT,CT X2>(CT(1),m1,x2); }

template <class T Y> 
inline SUMMX_1<CT,CT X2> operator+(const GENMATRIX<CT X1>& m1, CCT x2)
{ return SUMMX_1<CT,CT X2>(CT(1),m1,CT(x2)); }

template <class T Y> 
inline SUMMX_1<CT,CT X2> operator+(const GENMATRIX<CT X1>& m1, VCT x2)
{ return SUMMX_1<CT,CT X2>(CT(1),m1,CT(x2)); }

// x+m

template <class T Y> 
inline SUMMX_1<T,T X2> operator+(T x1, const GENMATRIX<T X1>& m2)
{ return SUMMX_1<T,T X2>(T(1),m2,x1); }

template <class T Y> 
inline SUMMX_1<CT,T X2> operator+(CT x1, const GENMATRIX<T X1>& m2)
{ return SUMMX_1<CT,T X2>(CT(1),m2,x1); }

template <class T Y> 
inline SUMMX_1<CT,T X2> operator+(CCT x1, const GENMATRIX<T X1>& m2)
{ return SUMMX_1<CT,T X2>(CT(1),m2,CT(x1)); }

template <class T Y> 
inline SUMMX_1<CT,T X2> operator+(VCT x1, const GENMATRIX<T X1>& m2)
{ return SUMMX_1<CT,T X2>(CT(1),m2,CT(x1)); }

template <class T Y> 
inline SUMMX_1<CT,CT X2> operator+(T x1, const GENMATRIX<CT X1>& m2)
{ return SUMMX_1<CT,CT X2>(CT(1),m2,x1); }

template <class T Y> 
inline SUMMX_1<CT,CT X2> operator+(CCT x1, const GENMATRIX<CT X1>& m2)
{ return SUMMX_1<CT,CT X2>(CT(1),m2,CT(x1)); }

template <class T Y> 
inline SUMMX_1<CT,CT X2> operator+(VCT x1, const GENMATRIX<CT X1>& m2)
{ return SUMMX_1<CT,CT X2>(CT(1),m2,CT(x1)); }

// m-x

template <class T Y> 
inline SUMMX_1<T,T X2> operator-(const GENMATRIX<T X1>& m1, T x2)
{ return SUMMX_1<T,T X2>(T(1),m1,-x2); }

template <class T Y> 
inline SUMMX_1<CT,T X2> operator-(const GENMATRIX<T X1>& m1, CT x2)
{ return SUMMX_1<CT,T X2>(CT(1),m1,-x2); }

template <class T Y> 
inline SUMMX_1<CT,T X2> operator-(const GENMATRIX<T X1>& m1, CCT x2)
{ return SUMMX_1<CT,T X2>(CT(1),m1,-CT(x2)); }

template <class T Y> 
inline SUMMX_1<CT,T X2> operator-(const GENMATRIX<T X1>& m1, VCT x2)
{ return SUMMX_1<CT,T X2>(CT(1),m1,-CT(x2)); }

template <class T Y> 
inline SUMMX_1<CT,CT X2> operator-(const GENMATRIX<CT X1>& m1, T x2)
{ return SUMMX_1<CT,CT X2>(CT(1),m1,-x2); }

template <class T Y> 
inline SUMMX_1<CT,CT X2> operator-(const GENMATRIX<CT X1>& m1, CCT x2)
{ return SUMMX_1<CT,CT X2>(CT(1),m1,-CT(x2)); }

template <class T Y> 
inline SUMMX_1<CT,CT X2> operator-(const GENMATRIX<CT X1>& m1, VCT x2)
{ return SUMMX_1<CT,CT X2>(CT(1),m1,-CT(x2)); }

// x-m

template <class T Y> 
inline SUMMX<T,T X2> operator-(T x1, const GENMATRIX<T X1>& m2)
{ return SUMMX<T,T X2>(T(-1),m2,x1); }

template <class T Y> 
inline SUMMX<CT,T X2> operator-(CT x1, const GENMATRIX<T X1>& m2)
{ return SUMMX<CT,T X2>(CT(-1),m2,x1); }

template <class T Y> 
inline SUMMX<CT,T X2> operator-(CCT x1, const GENMATRIX<T X1>& m2)
{ return SUMMX<CT,T X2>(CT(-1),m2,CT(x1)); }

template <class T Y> 
inline SUMMX<CT,T X2> operator-(VCT x1, const GENMATRIX<T X1>& m2)
{ return SUMMX<CT,T X2>(CT(-1),m2,CT(x1)); }

template <class T Y> 
inline SUMMX<CT,CT X2> operator-(T x1, const GENMATRIX<CT X1>& m2)
{ return SUMMX<CT,CT X2>(CT(-1),m2,x1); }

template <class T Y> 
inline SUMMX<CT,CT X2> operator-(CCT x1, const GENMATRIX<CT X1>& m2)
{ return SUMMX<CT,CT X2>(CT(-1),m2,CT(x1)); }

template <class T Y> 
inline SUMMX<CT,CT X2> operator-(VCT x1, const GENMATRIX<CT X1>& m2)
{ return SUMMX<CT,CT X2>(CT(-1),m2,CT(x1)); }

// -(x*m+x)

template <class T, class T1 Y> 
inline SUMMX<T,T1 X2> operator-(const SUMMX<T,T1 X2>& smx)
{ return SUMMX<T,T1 X2>(-smx.getX1(),smx.getM(),-smx.getX2); }

// x*(x*m+x)

template <class T, class T1 Y> 
inline SUMMX<T,T1 X2> operator*(const T x, const SUMMX<T,T1 X2>& smx)
{ return SUMMX<T,T1 X2>(smx.getX1()*x,smx.getM(),smx.getX2()*x); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator*(const T x, const SUMMX<CT,T1 X2>& smx)
{ return SUMMX<CT,T1 X2>(smx.getX1()*x,smx.getM(),smx.getX2()*x); }

template <class T Y> 
inline SUMMX<CT,T X2> operator*(const CT x, const SUMMX<T,T X2>& smx)
{ return SUMMX<CT,T X2>(smx.getX1()*x,smx.getM(),smx.getX2()*x); }

template <class T Y> 
inline SUMMX<CT,T X2> operator*(const CCT x, const SUMMX<T,T X2>& smx)
{ return SUMMX<CT,T X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <class T Y> 
inline SUMMX<CT,T X2> operator*(const VCT x, const SUMMX<T,T X2>& smx)
{ return SUMMX<CT,T X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator*(const CCT x, const SUMMX<CT,T1 X2>& smx)
{ return SUMMX<CT,T1 X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator*(const VCT x, const SUMMX<CT,T1 X2>& smx)
{ return SUMMX<CT,T1 X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

// (x*m+x)*x

template <class T, class T1 Y> 
inline SUMMX<T,T1 X2> operator*(const SUMMX<T,T1 X2>& smx, const T x)
{ return SUMMX<T,T1 X2>(smx.getX1()*x,smx.getM(),smx.getX2()*x); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator*(const SUMMX<CT,T1 X2>& smx, const T x)
{ return SUMMX<CT,T1 X2>(smx.getX1()*x,smx.getM(),smx.getX2()*x); }

template <class T Y> 
inline SUMMX<CT,T X2> operator*(const SUMMX<T,T X2>& smx, const CT x)
{ return SUMMX<CT,T X2>(smx.getX1()*x,smx.getM(),smx.getX2()*x); }

template <class T Y> 
inline SUMMX<CT,T X2> operator*(const SUMMX<T,T X2>& smx, const CCT x)
{ return SUMMX<CT,T X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <class T Y> 
inline SUMMX<CT,T X2> operator*(const SUMMX<T,T X2>& smx, const VCT x)
{ return SUMMX<CT,T X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator*(const SUMMX<CT,T1 X2>& smx, const CCT x)
{ return SUMMX<CT,T1 X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator*(const SUMMX<CT,T1 X2>& smx, const VCT x)
{ return SUMMX<CT,T1 X2>(smx.getX1()*CT(x),smx.getM(),smx.getX2()*CT(x)); }

// (x*m+x)/x

template <class T, class T1 Y> 
inline SUMMX<T,T1 X2> operator/(const SUMMX<T,T1 X2>& smx, const T x)
{ return SUMMX<T,T1 X2>(smx.getX1()/x,smx.getM(),smx.getX2()/x); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator/(const SUMMX<CT,T1 X2>& smx, const T x)
{ return SUMMX<CT,T1 X2>(smx.getX1()/x,smx.getM(),smx.getX2()/x); }

template <class T Y> 
inline SUMMX<CT,T X2> operator/(const SUMMX<T,T X2>& smx, const CT x)
{ return SUMMX<CT,T X2>(smx.getX1()/x,smx.getM(),smx.getX2()/x); }

template <class T Y> 
inline SUMMX<CT,T X2> operator/(const SUMMX<T,T X2>& smx, const CCT x)
{ return SUMMX<CT,T X2>(smx.getX1()/CT(x),smx.getM(),smx.getX2()/CT(x)); }

template <class T Y> 
inline SUMMX<CT,T X2> operator/(const SUMMX<T,T X2>& smx, const VCT x)
{ return SUMMX<CT,T X2>(smx.getX1()/CT(x),smx.getM(),smx.getX2()/CT(x)); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator/(const SUMMX<CT,T1 X2>& smx, const CCT x)
{ return SUMMX<CT,T1 X2>(smx.getX1()/CT(x),smx.getM(),smx.getX2()/CT(x)); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator/(const SUMMX<CT,T1 X2>& smx, const VCT x)
{ return SUMMX<CT,T1 X2>(smx.getX1()/CT(x),smx.getM(),smx.getX2()/CT(x)); }

// x+(x*m)

template <class T, class T2 Y> 
inline SUMMX<T,T2 X2> operator+(const T x3, const PRODXM<T,T2 X1>& pxm)
{ return SUMMX<T,T2 X2>(pxm.getX1(),pxm.getM2(),x3); }

template <class T, class T2 Y> 
inline SUMMX<CT,T2 X2> operator+(const T x3, const PRODXM<CT,T2 X1>& pxm)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),x3); }

template <class T Y> 
inline SUMMX<CT,T X2> operator+(const CT x3, const PRODXM<T,T X1>& pxm)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),x3); }

template <class T Y> 
inline SUMMX<CT,T X2> operator+(const CCT x3, const PRODXM<T,T X1>& pxm)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

template <class T Y> 
inline SUMMX<CT,T X2> operator+(const VCT x3, const PRODXM<T,T X1>& pxm)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

template <class T, class T2 Y> 
inline SUMMX<CT,T2 X2> operator+(const CCT x3, const PRODXM<CT,T2 X1>& pxm)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

template <class T, class T2 Y> 
inline SUMMX<CT,T2 X2> operator+(const VCT x3, const PRODXM<CT,T2 X1>& pxm)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

// (x*m)+x

template <class T, class T2 Y> 
inline SUMMX<T,T2 X2> operator+(const PRODXM<T,T2 X1>& pxm, const T x3)
{ return SUMMX<T,T2 X2>(pxm.getX1(),pxm.getM2(),x3); }

template <class T, class T2 Y> 
inline SUMMX<CT,T2 X2> operator+(const PRODXM<CT,T2 X1>& pxm, const T x3)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),x3); }

template <class T Y> 
inline SUMMX<CT,T X2> operator+(const PRODXM<T,T X1>& pxm, const CT x3)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),x3); }

template <class T Y> 
inline SUMMX<CT,T X2> operator+(const PRODXM<T,T X1>& pxm, const CCT x3)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

template <class T Y> 
inline SUMMX<CT,T X2> operator+(const PRODXM<T,T X1>& pxm, const VCT x3)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

template <class T, class T2 Y> 
inline SUMMX<CT,T2 X2> operator+(const PRODXM<CT,T2 X1>& pxm, const CCT x3)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

template <class T, class T2 Y> 
inline SUMMX<CT,T2 X2> operator+(const PRODXM<CT,T2 X1>& pxm, const VCT x3)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),CT(x3)); }

// x-(x*m)

template <class T, class T2 Y> 
inline SUMMX<T,T2 X2> operator-(const T x3, const PRODXM<T,T2 X1>& pxm)
{ return SUMMX<T,T2 X2>(-pxm.getX1(),pxm.getM2(),x3); }

template <class T, class T2 Y> 
inline SUMMX<CT,T2 X2> operator-(const T x3, const PRODXM<CT,T2 X1>& pxm)
{ return SUMMX<CT,T2 X2>(-pxm.getX1(),pxm.getM2(),x3); }

template <class T Y> 
inline SUMMX<CT,T X2> operator-(const CT x3, const PRODXM<T,T X1>& pxm)
{ return SUMMX<CT,T X2>(-pxm.getX1(),pxm.getM2(),x3); }

template <class T Y> 
inline SUMMX<CT,T X2> operator-(const CCT x3, const PRODXM<T,T X1>& pxm)
{ return SUMMX<CT,T X2>(-pxm.getX1(),pxm.getM2(),CT(x3)); }

template <class T Y> 
inline SUMMX<CT,T X2> operator-(const VCT x3, const PRODXM<T,T X1>& pxm)
{ return SUMMX<CT,T X2>(-pxm.getX1(),pxm.getM2(),CT(x3)); }

template <class T, class T2 Y> 
inline SUMMX<CT,T2 X2> operator-(const CCT x3, const PRODXM<CT,T2 X1>& pxm)
{ return SUMMX<CT,T2 X2>(-pxm.getX1(),pxm.getM2(),CT(x3)); }

template <class T, class T2 Y> 
inline SUMMX<CT,T2 X2> operator-(const VCT x3, const PRODXM<CT,T2 X1>& pxm)
{ return SUMMX<CT,T2 X2>(-pxm.getX1(),pxm.getM2(),CT(x3)); }

// (x*m)-x

template <class T, class T2 Y> 
inline SUMMX<T,T2 X2> operator-(const PRODXM<T,T2 X1>& pxm, const T x3)
{ return SUMMX<T,T2 X2>(pxm.getX1(),pxm.getM2(),-x3); }

template <class T, class T2 Y> 
inline SUMMX<CT,T2 X2> operator-(const PRODXM<CT,T2 X1>& pxm, const T x3)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),-x3); }

template <class T Y> 
inline SUMMX<CT,T X2> operator-(const PRODXM<T,T X1>& pxm, const CT x3)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),-x3); }

template <class T Y> 
inline SUMMX<CT,T X2> operator-(const PRODXM<T,T X1>& pxm, const CCT x3)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),-CT(x3)); }

template <class T Y> 
inline SUMMX<CT,T X2> operator-(const PRODXM<T,T X1>& pxm, const VCT x3)
{ return SUMMX<CT,T X2>(pxm.getX1(),pxm.getM2(),-CT(x3)); }

template <class T, class T2 Y> 
inline SUMMX<CT,T2 X2> operator-(const PRODXM<CT,T2 X1>& pxm, const CCT x3)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),-CT(x3)); }

template <class T, class T2 Y> 
inline SUMMX<CT,T2 X2> operator-(const PRODXM<CT,T2 X1>& pxm, const VCT x3)
{ return SUMMX<CT,T2 X2>(pxm.getX1(),pxm.getM2(),-CT(x3)); }

#ifndef NO_1
// -(x*m+x)

template <class T, class T1 Y> 
inline SUMMX<T,T1 X2> operator-(const SUMMX_1<T,T1 X2>& smx)
{ return SUMMX<T,T1 X2>(T(-1),smx.getM(),-smx.getX2); }

// x*(x*m+x)

template <class T, class T1 Y> 
inline SUMMX<T,T1 X2> operator*(const T x, const SUMMX_1<T,T1 X2>& smx)
{ return SUMMX<T,T1 X2>(x,smx.getM(),smx.getX2()*x); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator*(const T x, const SUMMX_1<CT,T1 X2>& smx)
{ return SUMMX<CT,T1 X2>(x,smx.getM(),smx.getX2()*x); }

template <class T Y> 
inline SUMMX<CT,T X2> operator*(const CT x, const SUMMX_1<T,T X2>& smx)
{ return SUMMX<CT,T X2>(x,smx.getM(),smx.getX2()*x); }

template <class T Y> 
inline SUMMX<CT,T X2> operator*(const CCT x, const SUMMX_1<T,T X2>& smx)
{ return SUMMX<CT,T X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <class T Y> 
inline SUMMX<CT,T X2> operator*(const VCT x, const SUMMX_1<T,T X2>& smx)
{ return SUMMX<CT,T X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator*(const CCT x, const SUMMX_1<CT,T1 X2>& smx)
{ return SUMMX<CT,T1 X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator*(const VCT x, const SUMMX_1<CT,T1 X2>& smx)
{ return SUMMX<CT,T1 X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

// (x*m+x)*x

template <class T, class T1 Y> 
inline SUMMX<T,T1 X2> operator*(const SUMMX_1<T,T1 X2>& smx, const T x)
{ return SUMMX<T,T1 X2>(x,smx.getM(),smx.getX2()*x); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator*(const SUMMX_1<CT,T1 X2>& smx, const T x)
{ return SUMMX<CT,T1 X2>(x,smx.getM(),smx.getX2()*x); }

template <class T Y> 
inline SUMMX<CT,T X2> operator*(const SUMMX_1<T,T X2>& smx, const CT x)
{ return SUMMX<CT,T X2>(x,smx.getM(),smx.getX2()*x); }

template <class T Y> 
inline SUMMX<CT,T X2> operator*(const SUMMX_1<T,T X2>& smx, const CCT x)
{ return SUMMX<CT,T X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <class T Y> 
inline SUMMX<CT,T X2> operator*(const SUMMX_1<T,T X2>& smx, const VCT x)
{ return SUMMX<CT,T X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator*(const SUMMX_1<CT,T1 X2>& smx, const CCT x)
{ return SUMMX<CT,T1 X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator*(const SUMMX_1<CT,T1 X2>& smx, const VCT x)
{ return SUMMX<CT,T1 X2>(CT(x),smx.getM(),smx.getX2()*CT(x)); }

// (x*m+x)/x

template <class T, class T1 Y> 
inline SUMMX<T,T1 X2> operator/(const SUMMX_1<T,T1 X2>& smx, const T x)
{ return SUMMX<T,T1 X2>(TMV_RealType(T)(1)/x,smx.getM(),smx.getX2()/x); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator/(const SUMMX_1<CT,T1 X2>& smx, const T x)
{ return SUMMX<CT,T1 X2>(T(1)/x,smx.getM(),smx.getX2()/x); }

template <class T Y> 
inline SUMMX<CT,T X2> operator/(const SUMMX_1<T,T X2>& smx, const CT x)
{ return SUMMX<CT,T X2>(T(1)/x,smx.getM(),smx.getX2()/x); }

template <class T Y> 
inline SUMMX<CT,T X2> operator/(const SUMMX_1<T,T X2>& smx, const CCT x)
{ return SUMMX<CT,T X2>(T(1)/CT(x),smx.getM(),smx.getX2()/CT(x)); }

template <class T Y> 
inline SUMMX<CT,T X2> operator/(const SUMMX_1<T,T X2>& smx, const VCT x)
{ return SUMMX<CT,T X2>(T(1)/CT(x),smx.getM(),smx.getX2()/CT(x)); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator/(const SUMMX_1<CT,T1 X2>& smx, const CCT x)
{ return SUMMX<CT,T1 X2>(T(1)/CT(x),smx.getM(),smx.getX2()/CT(x)); }

template <class T, class T1 Y> 
inline SUMMX<CT,T1 X2> operator/(const SUMMX_1<CT,T1 X2>& smx, const VCT x)
{ return SUMMX<CT,T1 X2>(T(1)/CT(x),smx.getM(),smx.getX2()/CT(x)); }
#endif

#undef X1
#undef X2
#undef Y
#undef SUMMX_1
#ifdef NO_1
#undef NO_1
#endif
