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


// Need to define the following with #define statements.
// (The given definition is for a regular Matrix/Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX1 GenMatrix
// #define GENMATRIX2 GenMatrix
// #define QUOTMM QuotMM
// #define RQUOTMM RQuotMM
// #define PRODXM1 ProdXM
// #define PRODXM2 ProdXM
// #define QUOTXM QuotXM

#ifndef X1
#define X1
#endif

#ifndef X1b
#define X1b X1
#endif

#ifndef X2
#define X2
#endif

#ifndef X2b
#define X2b X2
#endif

#ifndef X3
#define X3
#endif

#ifndef X3b
#define X3b X3
#endif

#ifndef Y
#define Y
#endif

#ifndef Yb
#define Yb Y
#endif

#ifndef GETM1
#define GETM1 .getM()
#endif

#ifndef GETM2
#define GETM2 .getM()
#endif

#ifndef QUOTMM_1
#define QUOTMM_1 QUOTMM
#endif

#ifndef RQUOTMM_1
#define RQUOTMM_1 RQUOTMM
#endif

// m/m

template <class T Y> 
inline QUOTMM_1<T,T,T X3> operator/(
    const GENMATRIX1<T X1>& m1, const GENMATRIX2<T X2>& m2)
{ return QUOTMM_1<T,T,T X3>(T(1),m1,m2); }

template <class T Y> 
inline QUOTMM_1<CT,CT,T X3> operator/(
    const GENMATRIX1<CT X1>& m1, const GENMATRIX2<T X2>& m2)
{ return QUOTMM_1<CT,CT,T X3>(CT(1),m1,m2); }

template <class T Y> 
inline QUOTMM_1<CT,T,CT X3> operator/(
    const GENMATRIX1<T X1>& m1, const GENMATRIX2<CT X2>& m2)
{ return QUOTMM_1<CT,T,CT X3>(CT(1),m1,m2); }

// m%m

template <class T Yb> 
inline RQUOTMM_1<T,T,T X3b> operator%(
    const GENMATRIX1<T X1b>& m1, const GENMATRIX2<T X2b>& m2)
{ return RQUOTMM_1<T,T,T X3b>(T(1),m1,m2); }

template <class T Yb> 
inline RQUOTMM_1<CT,CT,T X3b> operator%(
    const GENMATRIX1<CT X1b>& m1, const GENMATRIX2<T X2b>& m2)
{ return RQUOTMM_1<CT,CT,T X3b>(CT(1),m1,m2); }

template <class T Yb> 
inline RQUOTMM_1<CT,T,CT X3b> operator%(
    const GENMATRIX1<T X1b>& m1, const GENMATRIX2<CT X2b>& m2)
{ return RQUOTMM_1<CT,T,CT X3b>(CT(1),m1,m2); }

// (x*m)/m

template <class T, class T1 Y> 
inline QUOTMM<T,T1,T X3> operator/(
    const PRODXM1<T,T1 X1>& pxm, const GENMATRIX2<T X2>& m)
{ return QUOTMM<T,T1,T X3>(pxm.getX(),pxm GETM1,m); }

template <class T Y> 
inline QUOTMM<CT,T,CT X3> operator/(
    const PRODXM1<T,T X1>& pxm, const GENMATRIX2<CT X2>& m)
{ return QUOTMM<CT,T,CT X3>(CT(pxm.getX()),pxm GETM1,m); }

template <class T, class T1 Y> 
inline QUOTMM<CT,T1,T X3> operator/(
    const PRODXM1<CT,T1 X1>& pxm, const GENMATRIX2<T X2>& m)
{ return QUOTMM<CT,T1,T X3>(pxm.getX(),pxm GETM1,m); }

// (x*m)%m

template <class T, class T1 Yb> 
inline RQUOTMM<T,T1,T X3b> operator%(
    const PRODXM1<T,T1 X1b>& pxm, const GENMATRIX2<T X2b>& m)
{ return RQUOTMM<T,T1,T X3b>(pxm.getX(),pxm GETM1,m); }

template <class T Yb> 
inline RQUOTMM<CT,T,CT X3b> operator%(
    const PRODXM1<T,T X1b>& pxm, const GENMATRIX2<CT X2b>& m)
{ return RQUOTMM<CT,T,CT X3b>(CT(pxm.getX()),pxm GETM1,m); }

template <class T, class T1 Yb> 
inline RQUOTMM<CT,T1,T X3b> operator%(
    const PRODXM1<CT,T1 X1b>& pxm, const GENMATRIX2<T X2b>& m)
{ return RQUOTMM<CT,T1,T X3b>(pxm.getX(),pxm GETM1,m); }

// m/(x*m)

template <class T, class T2 Y> 
inline QUOTMM<T,T,T2 X3> operator/(
    const GENMATRIX1<T X1>& m, const PRODXM2<T,T2 X2>& pxm)
{ return QUOTMM<T,T,T2 X3>(TMV_RealType(T)(1)/pxm.getX(),m,pxm GETM2); }

template <class T Y> 
inline QUOTMM<CT,CT,T X3> operator/(
    const GENMATRIX1<CT X1>& m, const PRODXM2<T,T X2>& pxm)
{ return QUOTMM<CT,CT,T X3>(CT(T(1)/pxm.getX()),m,pxm GETM2); }

template <class T, class T2 Y> 
inline QUOTMM<CT,T,T2 X3> operator/(
    const GENMATRIX1<T X1>& m, const PRODXM2<CT,T2 X2>& pxm)
{ return QUOTMM<CT,T,T2 X3>(T(1)/pxm.getX(),m,pxm GETM2); }

// m%(x*m)

template <class T, class T2 Yb> 
inline RQUOTMM<T,T,T2 X3b> operator%(
    const GENMATRIX1<T X1b>& m, const PRODXM2<T,T2 X2b>& pxm)
{ return RQUOTMM<T,T,T2 X3b>(TMV_RealType(T)(1)/pxm.getX(),m,pxm GETM2); }

template <class T Yb> 
inline RQUOTMM<CT,CT,T X3b> operator%(
    const GENMATRIX1<CT X1b>& m, const PRODXM2<T,T X2b>& pxm)
{ return RQUOTMM<CT,CT,T X3b>(CT(T(1)/pxm.getX()),m,pxm GETM2); }

template <class T, class T2 Yb> 
inline RQUOTMM<CT,T,T2 X3b> operator%(
    const GENMATRIX1<T X1b>& m, const PRODXM2<CT,T2 X2b>& pxm)
{ return RQUOTMM<CT,T,T2 X3b>(T(1)/pxm.getX(),m,pxm GETM2); }

// (x*m)/(x*m)

template <class T, class T1, class T2 Y> 
inline QUOTMM<T,T1,T2 X3> operator/(
    const PRODXM1<T,T1 X1>& pxm1, const PRODXM2<T,T2 X2>& pxm2)
{ return QUOTMM<T,T1,T2 X3>(pxm1.getX()/pxm2.getX(),pxm1 GETM1,pxm2 GETM2); }

template <class T, class T1 Y> 
inline QUOTMM<CT,T1,T X3> operator/(
    const PRODXM1<CT,T1 X1>& pxm1, const PRODXM2<T,T X2>& pxm2)
{ return QUOTMM<CT,T1,T X3>(pxm1.getX()/pxm2.getX(),pxm1 GETM1,pxm2 GETM2); }

template <class T, class T2 Y> 
inline QUOTMM<CT,T,T2 X3> operator/(
    const PRODXM1<T,T X1>& pxm1, const PRODXM2<CT,T2 X2>& pxm2)
{ return QUOTMM<CT,T,T2 X3>(pxm1.getX()/pxm2.getX(),pxm1 GETM1,pxm2 GETM2); }

// (x*m)%(x*m)

template <class T, class T1, class T2 Yb> 
inline RQUOTMM<T,T1,T2 X3b> operator%(
    const PRODXM1<T,T1 X1b>& pxm1, const PRODXM2<T,T2 X2b>& pxm2)
{ return RQUOTMM<T,T1,T2 X3b>(pxm1.getX()/pxm2.getX(),pxm1 GETM1,pxm2 GETM2); }

template <class T, class T1 Yb> 
inline RQUOTMM<CT,T1,T X3b> operator%(
    const PRODXM1<CT,T1 X1b>& pxm1, const PRODXM2<T,T X2b>& pxm2)
{ return RQUOTMM<CT,T1,T X3b>(pxm1.getX()/pxm2.getX(),pxm1 GETM1,pxm2 GETM2); }

template <class T, class T2 Yb> 
inline RQUOTMM<CT,T,T2 X3b> operator%(
    const PRODXM1<T,T X1b>& pxm1, const PRODXM2<CT,T2 X2b>& pxm2)
{ return RQUOTMM<CT,T,T2 X3b>(pxm1.getX()/pxm2.getX(),pxm1 GETM1,pxm2 GETM2); }

#ifdef QUOTXM
// (x/m)*m

template <class T, class T2 Y> 
inline QUOTMM<T,T,T2 X3> operator*(
    const QUOTXM<T,T2 X2>& qxm, const GENMATRIX1<T X1>& m)
{ return QUOTMM<T,T,T2 X3>(qxm.getX(),m,qxm.getM()); }

template <class T Y> 
inline QUOTMM<CT,CT,T X3> operator*(
    const QUOTXM<T,T X2>& qxm, const GENMATRIX1<CT X1>& m)
{ return QUOTMM<CT,CT,T X3>(CT(qxm.getX()),m,qxm.getM()); }

template <class T, class T2 Y> 
inline QUOTMM<CT,T,T2 X3> operator*(
    const QUOTXM<CT,T2 X2>& qxm, const GENMATRIX1<T X1>& m)
{ return QUOTMM<CT,T,T2 X3>(qxm.getX(),m,qxm.getM()); }

// m*(x/m)

template <class T, class T2 Yb> 
inline RQUOTMM<T,T,T2 X3b> operator*(
    const GENMATRIX1<T X1b>& m, const QUOTXM<T,T2 X2b>& qxm)
{ return RQUOTMM<T,T,T2 X3b>(qxm.getX(),m,qxm.getM()); }

template <class T Yb> 
inline RQUOTMM<CT,CT,T X3b> operator*(
    const GENMATRIX1<CT X1b>& m, const QUOTXM<T,T X2b>& qxm)
{ return RQUOTMM<CT,CT,T X3b>(CT(qxm.getX()),m,qxm.getM()); }

template <class T, class T2 Yb> 
inline RQUOTMM<CT,T,T2 X3b> operator*(
    const GENMATRIX1<T X1b>& m, const QUOTXM<CT,T2 X2b>& qxm)
{ return RQUOTMM<CT,T,T2 X3b>(qxm.getX(),m,qxm.getM()); }

// (x/m)*(x*m)

template <class T, class T1, class T2 Y> 
inline QUOTMM<T,T,T2 X3> operator*(
    const QUOTXM<T,T2 X2>& qxm, const PRODXM1<T,T1 X1>& pxm)
{ return QUOTMM<T,T,T2 X3>(pxm.getX()*qxm.getX(),pxm GETM2,qxm.getM()); }

template <class T, class T1 Y> 
inline QUOTMM<CT,CT,T X3> operator*(
    const QUOTXM<T,T X2>& qxm, const PRODXM1<CT,T1 X1>& pxm)
{ return QUOTMM<CT,CT,T X3>(pxm.getX()*qxm.getX(),pxm GETM2,qxm.getM()); }

template <class T, class T2 Y> 
inline QUOTMM<CT,T,T2 X3> operator*(
    const QUOTXM<CT,T2 X2>& qxm, const PRODXM1<T,T X1>& pxm)
{ return QUOTMM<CT,T,T2 X3>(pxm.getX()*qxm.getX(),pxm GETM2,qxm.getM()); }

// (x*m)*(x/m)

template <class T, class T1, class T2 Yb> 
inline RQUOTMM<T,T1,T2 X3b> operator*(
    const PRODXM1<T,T1 X1b>& pxm, const QUOTXM<T,T2 X2b>& qxm)
{ return RQUOTMM<T,T1,T2 X3b>(pxm.getX()*qxm.getX(),pxm GETM1,qxm.getM()); }

template <class T, class T1 Yb> 
inline RQUOTMM<CT,T1,T X3b> operator*(
    const PRODXM1<CT,T1 X1b>& pxm, const QUOTXM<T,T X2b>& qxm)
{ return RQUOTMM<CT,T1,T X3b>(pxm.getX()*qxm.getX(),pxm GETM1,qxm.getM()); }

template <class T, class T2 Yb> 
inline RQUOTMM<CT,T,T2 X3b> operator*(
    const PRODXM1<T,T X1b>& pxm, const QUOTXM<CT,T2 X2b>& qxm)
{ return RQUOTMM<CT,T,T2 X3b>(pxm.getX()*qxm.getX(),pxm GETM1,qxm.getM()); }
#endif

#ifdef QUOTXM_1
// (1/m)*m

template <class T, class T2 Y> 
inline QUOTMM_1<T,T,T2 X3> operator*(
    const QUOTXM_1<T,T2 X2>& qxm, const GENMATRIX1<T X1>& m)
{ return QUOTMM_1<T,T,T2 X3>(T(1),m,qxm.getM()); }

template <class T Y> 
inline QUOTMM_1<CT,CT,T X3> operator*(
    const QUOTXM_1<T,T X2>& qxm, const GENMATRIX1<CT X1>& m)
{ return QUOTMM_1<CT,CT,T X3>(CT(1),m,qxm.getM()); }

template <class T, class T2 Y> 
inline QUOTMM_1<CT,T,T2 X3> operator*(
    const QUOTXM_1<CT,T2 X2>& qxm, const GENMATRIX1<T X1>& m)
{ return QUOTMM_1<CT,T,T2 X3>(CT(1),m,qxm.getM()); }

// m*(1/m)

template <class T, class T2 Yb> 
inline RQUOTMM_1<T,T,T2 X3b> operator*(
    const GENMATRIX1<T X1b>& m, const QUOTXM_1<T,T2 X2b>& qxm)
{ return RQUOTMM_1<T,T,T2 X3b>(T(1),m,qxm.getM()); }

template <class T Yb> 
inline RQUOTMM_1<CT,CT,T X3b> operator*(
    const GENMATRIX1<CT X1b>& m, const QUOTXM_1<T,T X2b>& qxm)
{ return RQUOTMM_1<CT,CT,T X3b>(CT(1),m,qxm.getM()); }

template <class T, class T2 Yb> 
inline RQUOTMM_1<CT,T,T2 X3b> operator*(
    const GENMATRIX1<T X1b>& m, const QUOTXM_1<CT,T2 X2b>& qxm)
{ return RQUOTMM_1<CT,T,T2 X3b>(CT(1),m,qxm.getM()); }

// (1/m)*(x*m)

template <class T, class T1, class T2 Y> 
inline QUOTMM<T,T,T2 X3> operator*(
    const QUOTXM_1<T,T2 X2>& qxm, const PRODXM1<T,T1 X1>& pxm)
{ return QUOTMM<T,T,T2 X3>(pxm.getX(),pxm GETM2,qxm.getM()); }

template <class T, class T1 Y> 
inline QUOTMM<CT,CT,T X3> operator*(
    const QUOTXM_1<T,T X2>& qxm, const PRODXM1<CT,T1 X1>& pxm)
{ return QUOTMM<CT,CT,T X3>(pxm.getX(),pxm GETM2,qxm.getM()); }

template <class T, class T2 Y> 
inline QUOTMM<CT,T,T2 X3> operator*(
    const QUOTXM_1<CT,T2 X2>& qxm, const PRODXM1<T,T X1>& pxm)
{ return QUOTMM<CT,T,T2 X3>(pxm.getX(),pxm GETM2,qxm.getM()); }

// (x*m)*(1/m)

template <class T, class T1, class T2 Yb> 
inline RQUOTMM<T,T1,T2 X3b> operator*(
    const PRODXM1<T,T1 X1b>& pxm, const QUOTXM_1<T,T2 X2b>& qxm)
{ return RQUOTMM<T,T1,T2 X3b>(pxm.getX(),pxm GETM1,qxm.getM()); }

template <class T, class T1 Yb> 
inline RQUOTMM<CT,T1,T X3b> operator*(
    const PRODXM1<CT,T1 X1b>& pxm, const QUOTXM_1<T,T X2b>& qxm)
{ return RQUOTMM<CT,T1,T X3b>(pxm.getX(),pxm GETM1,qxm.getM()); }

template <class T, class T2 Yb> 
inline RQUOTMM<CT,T,T2 X3b> operator*(
    const PRODXM1<T,T X1b>& pxm, const QUOTXM_1<CT,T2 X2b>& qxm)
{ return RQUOTMM<CT,T,T2 X3b>(pxm.getX(),pxm GETM1,qxm.getM()); }

#undef QUOTXM_1
#endif

#undef X1
#undef X1b
#undef X2
#undef X2b
#undef X3
#undef X3b
#undef Y
#undef Yb
#undef GETM1
#undef GETM2
#undef QUOTMM_1
#undef RQUOTMM_1
