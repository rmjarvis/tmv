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

template <class T Y> 
inline SUMMM<T,T,T X3> operator+(
    const GENMATRIX2<T X2>& m2, const GENMATRIX1<T X1>& m1)
{ return SUMMM<T,T,T X3>(T(1),m1,T(1),m2); }

template <class T Y> 
inline SUMMM<CT,T,CT X3> operator+(
    const GENMATRIX2<CT X2>& m2,const GENMATRIX1<T X1>& m1)
{ return SUMMM<CT,T,CT X3>(CT(1),m1,CT(1),m2); }

template <class T Y> 
inline SUMMM<CT,CT,T X3> operator+(
    const GENMATRIX2<T X2>& m2,const GENMATRIX1<CT X1>& m1)
{ return SUMMM<CT,CT,T X3>(CT(1),m1,CT(1),m2); }

// m-m

template <class T Y> 
inline SUMMM<T,T,T X3> operator-(
    const GENMATRIX2<T X2>& m2, const GENMATRIX1<T X1>& m1)
{ return SUMMM<T,T,T X3>(T(-1),m1,T(1),m2); }

template <class T Y> 
inline SUMMM<CT,T,CT X3> operator-(
    const GENMATRIX2<CT X2>& m2,const GENMATRIX1<T X1>& m1)
{ return SUMMM<CT,T,CT X3>(CT(-1),m1,CT(1),m2); }

template <class T Y> 
inline SUMMM<CT,CT,T X3> operator-(
    const GENMATRIX2<T X2>& m2,const GENMATRIX1<CT X1>& m1)
{ return SUMMM<CT,CT,T X3>(CT(-1),m1,CT(1),m2); }

// (x*m)+m

template <class T, class T1 Y> 
inline SUMMM<T,T1,T X3> operator+(
    const GENMATRIX2<T X2>& m, const PRODXM1<T,T1 X1>& pxm) 
{ return SUMMM<T,T1,T X3>(pxm.getX(),pxm.getM(),T(1),m); }

template <class T Y> 
inline SUMMM<CT,T,CT X3> operator+(
    const GENMATRIX2<CT X2>& m, const PRODXM1<T,T X1>& pxm) 
{ return SUMMM<CT,T,CT X3>(CT(pxm.getX()),pxm.getM(),CT(1),m); }

template <class T, class T1 Y> 
inline SUMMM<CT,T1,T X3> operator+(
    const GENMATRIX2<T X2>& m, const PRODXM1<CT,T1 X1>& pxm) 
{ return SUMMM<CT,T1,T X3>(pxm.getX(),pxm.getM(),CT(1),m); }

// m+(x*m)

template <class T, class T2 Y> 
inline SUMMM<T,T,T2 X3> operator+(
    const PRODXM2<T,T2 X2>& pxm, const GENMATRIX1<T X1>& m)
{ return SUMMM<T,T,T2 X3>(T(1),m,pxm.getX(),pxm.getM()); }

template <class T Y> 
inline SUMMM<CT,CT,T X3> operator+(
    const PRODXM2<T,T X2>& pxm, const GENMATRIX1<CT X1>& m) 
{ return SUMMM<CT,CT,T X3>(CT(1),m,CT(pxm.getX()),pxm.getM()); }

template <class T, class T2 Y> 
inline SUMMM<CT,T,T2 X3> operator+(
    const PRODXM2<CT,T2 X2>& pxm, const GENMATRIX1<T X1>& m) 
{ return SUMMM<CT,T,T2 X3>(CT(1),m,pxm.getX(),pxm.getM()); }

// (x*m)-m

template <class T, class T1 Y> 
inline SUMMM<T,T1,T X3> operator-(
    const GENMATRIX2<T X2>& m, const PRODXM1<T,T1 X1>& pxm)
{ return SUMMM<T,T1,T X3>(-pxm.getX(),pxm.getM(),T(1),m); }

template <class T Y> 
inline SUMMM<CT,T,CT X3> operator-(
    const GENMATRIX2<CT X2>& m, const PRODXM1<T,T X1>& pxm) 
{ return SUMMM<CT,T,CT X3>(CT(-pxm.getX()),pxm.getM(),CT(1),m); }

template <class T, class T1 Y> 
inline SUMMM<CT,T1,T X3> operator-(
    const GENMATRIX2<T X2>& m, const PRODXM1<CT,T1 X1>& pxm)
{ return SUMMM<CT,T1,T X3>(-pxm.getX(),pxm.getM(),CT(1),m); }

// m-(x*m)

template <class T, class T2 Y> 
inline SUMMM<T,T,T2 X3> operator-(
    const PRODXM2<T,T2 X2>& pxm, const GENMATRIX1<T X1>& m)
{ return SUMMM<T,T,T2 X3>(T(-1),m,pxm.getX(),pxm.getM()); }

template <class T Y> 
inline SUMMM<CT,CT,T X3> operator-(
    const PRODXM2<T,T X2>& pxm, const GENMATRIX1<CT X1>& m)
{ return SUMMM<CT,CT,T X3>(CT(-1),m,CT(pxm.getX()),pxm.getM()); }

template <class T, class T2 Y> 
inline SUMMM<CT,T,T2 X3> operator-(
    const PRODXM2<CT,T2 X2>& pxm, const GENMATRIX1<T X1>& m)
{ return SUMMM<CT,T,T2 X3>(CT(-1),m,pxm.getX(),pxm.getM()); }

// (x*m)+(x*m)

template <class T, class T1, class T2 Y> 
inline SUMMM<T,T1,T2 X3> operator+(
    const PRODXM2<T,T2 X2>& pxm2, const PRODXM1<T,T1 X1>& pxm1)
{ return SUMMM<T,T1,T2 X3>(pxm1.getX(),pxm1.getM(),pxm2.getX(),pxm2.getM()); }

template <class T, class T1 Y> 
inline SUMMM<CT,T1,T X3> operator+(
    const PRODXM2<T,T X2>& pxm2, const PRODXM1<CT,T1 X1>& pxm1)
{
    return SUMMM<CT,T1,T X3>(
        pxm1.getX(),pxm1.getM(),CT(pxm2.getX()),pxm2.getM()); 
}

template <class T, class T2 Y> 
inline SUMMM<CT,T,T2 X3> operator+(
    const PRODXM2<CT,T2 X2>& pxm2, const PRODXM1<T,T X1>& pxm1)
{
    return SUMMM<CT,T,T2 X3>(
        CT(pxm1.getX()),pxm1.getM(),pxm2.getX(),pxm2.getM()); 
}


// (x*m)-(x*m)

template <class T, class T1, class T2 Y> 
inline SUMMM<T,T1,T2 X3> operator-(
    const PRODXM2<T,T2 X2>& pxm2, const PRODXM1<T,T1 X1>& pxm1)
{ return SUMMM<T,T1,T2 X3>(-pxm1.getX(),pxm1.getM(),pxm2.getX(),pxm2.getM()); }

template <class T, class T1 Y> 
inline SUMMM<CT,T1,T X3> operator-(
    const PRODXM2<T,T X2>& pxm2, const PRODXM1<CT,T1 X1>& pxm1)
{
    return SUMMM<CT,T1,T X3>(
        -pxm1.getX(),pxm1.getM(),CT(pxm2.getX()),pxm2.getM()); 
}

template <class T, class T2 Y> 
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
