///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
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

// m+m
template <class T Y> inline SUMMM_1_1<T,T,T X3> operator+(
    const GENMATRIX1<T X1>& m1, const GENMATRIX2<T X2>& m2)
{ return SUMMM_1_1<T,T,T X3>(T(1),m1,T(1),m2); }

template <class T Y> inline SUMMM_1_1<CT,CT,T X3> operator+(
    const GENMATRIX1<CT X1>& m1,const GENMATRIX2<T X2>& m2)
{ return SUMMM_1_1<CT,CT,T X3>(CT(1),m1,CT(1),m2); }

template <class T Y> inline SUMMM_1_1<CT,T,CT X3> operator+(
    const GENMATRIX1<T X1>& m1,const GENMATRIX2<CT X2>& m2)
{ return SUMMM_1_1<CT,T,CT X3>(CT(1),m1,CT(1),m2); }

// m-m
template <class T Y> inline SUMMM_1_m1<T,T,T X3> operator-(
    const GENMATRIX1<T X1>& m1, const GENMATRIX2<T X2>& m2)
{ return SUMMM_1_m1<T,T,T X3>(T(1),m1,T(-1),m2); }

template <class T Y> inline SUMMM_1_m1<CT,CT,T X3> operator-(
    const GENMATRIX1<CT X1>& m1,const GENMATRIX2<T X2>& m2)
{ return SUMMM_1_m1<CT,CT,T X3>(CT(1),m1,CT(-1),m2); }

template <class T Y> inline SUMMM_1_m1<CT,T,CT X3> operator-(
    const GENMATRIX1<T X1>& m1,const GENMATRIX2<CT X2>& m2)
{ return SUMMM_1_m1<CT,T,CT X3>(CT(1),m1,CT(-1),m2); }

// (x*m)+m
template <class T, class T1 Y> inline SUMMM_x_1<T,T1,T X3> operator+(
    const PRODXM1<T,T1 X1>& pxm, const GENMATRIX2<T X2>& m)
{ return SUMMM_x_1<T,T1,T X3>(pxm.GetX(),pxm.GetM(),T(1),m); }

template <class T Y> inline SUMMM_x_1<CT,T,CT X3> operator+(
    const PRODXM1<T,T X1>& pxm, const GENMATRIX2<CT X2>& m)
{ return SUMMM_x_1<CT,T,CT X3>(CT(pxm.GetX()),pxm.GetM(),CT(1),m); }

template <class T, class T1 Y> inline SUMMM_x_1<CT,T1,T X3> operator+(
    const PRODXM1<CT,T1 X1>& pxm, const GENMATRIX2<T X2>& m)
{ return SUMMM_x_1<CT,T1,T X3>(pxm.GetX(),pxm.GetM(),CT(1),m); }

// m+(x*m)
template <class T, class T2 Y> inline SUMMM_1_x<T,T,T2 X3> operator+(
    const GENMATRIX1<T X1>& m, const PRODXM2<T,T2 X2>& pxm)
{ return SUMMM_1_x<T,T,T2 X3>(T(1),m,pxm.GetX(),pxm.GetM()); }

template <class T Y> inline SUMMM_1_x<CT,CT,T X3> operator+(
    const GENMATRIX1<CT X1>& m, const PRODXM2<T,T X2>& pxm)
{ return SUMMM_1_x<CT,CT,T X3>(CT(1),m,CT(pxm.GetX()),pxm.GetM()); }

template <class T, class T2 Y> inline SUMMM_1_x<CT,T,T2 X3> operator+(
    const GENMATRIX1<T X1>& m, const PRODXM2<CT,T2 X2>& pxm)
{ return SUMMM_1_x<CT,T,T2 X3>(CT(1),m,pxm.GetX(),pxm.GetM()); }

// (x*m)-m
template <class T, class T1 Y> inline SUMMM_x_m1<T,T1,T X3> operator-(
    const PRODXM1<T,T1 X1>& pxm, const GENMATRIX2<T X2>& m)
{ return SUMMM_x_m1<T,T1,T X3>(pxm.GetX(),pxm.GetM(),T(-1),m); }

template <class T Y> inline SUMMM_x_m1<CT,T,CT X3> operator-(
    const PRODXM1<T,T X1>& pxm, const GENMATRIX2<CT X2>& m)
{ return SUMMM_x_m1<CT,T,CT X3>(CT(pxm.GetX()),pxm.GetM(),CT(-1),m); }

template <class T, class T1 Y> inline SUMMM_x_m1<CT,T1,T X3> operator-(
    const PRODXM1<CT,T1 X1>& pxm, const GENMATRIX2<T X2>& m)
{ return SUMMM_x_m1<CT,T1,T X3>(pxm.GetX(),pxm.GetM(),CT(-1),m); }

// m-(x*m)
template <class T, class T2 Y> inline SUMMM_1_x<T,T,T2 X3> operator-(
    const GENMATRIX1<T X1>& m, const PRODXM2<T,T2 X2>& pxm)
{ return SUMMM_1_x<T,T,T2 X3>(T(1),m,-pxm.GetX(),pxm.GetM()); }

template <class T Y> inline SUMMM_1_x<CT,CT,T X3> operator-(
    const GENMATRIX1<CT X1>& m, const PRODXM2<T,T X2>& pxm)
{ return SUMMM_1_x<CT,CT,T X3>(CT(1),m,CT(-pxm.GetX()),pxm.GetM()); }

template <class T, class T2 Y> inline SUMMM_1_x<CT,T,T2 X3> operator-(
    const GENMATRIX1<T X1>& m, const PRODXM2<CT,T2 X2>& pxm)
{ return SUMMM_1_x<CT,T,T2 X3>(CT(1),m,-pxm.GetX(),pxm.GetM()); }

// (x*m)+(x*m)
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator+(
    const PRODXM1<T,T1 X1>& pxm1, const PRODXM2<T,T2 X2>& pxm2)
{ return SUMMM<T,T1,T2 X3>(pxm1.GetX(),pxm1.GetM(),pxm2.GetX(),pxm2.GetM()); }

template <class T, class T2 Y> inline SUMMM<CT,T,T2 X3> operator+(
    const PRODXM1<T,T X1>& pxm1, const PRODXM2<CT,T2 X2>& pxm2)
{ return SUMMM<CT,T,T2 X3>(CT(pxm1.GetX()),pxm1.GetM(),pxm2.GetX(),pxm2.GetM()); }

template <class T, class T1 Y> inline SUMMM<CT,T1,T X3> operator+(
    const PRODXM1<CT,T1 X1>& pxm1, const PRODXM2<T,T X2>& pxm2)
{ return SUMMM<CT,T1,T X3>(pxm1.GetX(),pxm1.GetM(),CT(pxm2.GetX()),pxm2.GetM()); }

// (x*m)-(x*m)
template <class T, class T1, class T2 Y> inline SUMMM<T,T1,T2 X3> operator-(
    const PRODXM1<T,T1 X1>& pxm1, const PRODXM2<T,T2 X2>& pxm2)
{ return SUMMM<T,T1,T2 X3>(pxm1.GetX(),pxm1.GetM(),-pxm2.GetX(),pxm2.GetM()); }

template <class T, class T2 Y> inline SUMMM<CT,T,T2 X3> operator-(
    const PRODXM1<T,T X1>& pxm1, const PRODXM2<CT,T2 X2>& pxm2)
{ return SUMMM<CT,T,T2 X3>(CT(pxm1.GetX()),pxm1.GetM(),-pxm2.GetX(),pxm2.GetM()); }

template <class T, class T1 Y> inline SUMMM<CT,T1,T X3> operator-(
    const PRODXM1<CT,T1 X1>& pxm1, const PRODXM2<T,T X2>& pxm2)
{ return SUMMM<CT,T1,T X3>(pxm1.GetX(),pxm1.GetM(),CT(-pxm2.GetX()),pxm2.GetM()); }

#undef X1
#undef X2
#undef X3
#undef Y

#undef SUMMM_1_1
#undef SUMMM_1_m1
#undef SUMMM_1_x
#undef SUMMM_x_1
#undef SUMMM_x_m1
