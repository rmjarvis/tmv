///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



// Need to define the following with #define statements.
// (The given definition is for a regular Matrix*Vector.  Modify as 
// appropriate for the various other matrices.)
//
// #define PRODMV ProdMV

#ifndef X3
#define X3
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

// x*(x*m*v)
template <class T, class T1, class T2 Y> 
inline PRODMV<T,T1,T2 X3> operator*(T x, const PRODMV<T,T1,T2 X3>& pmv)
{ return PRODMV<T,T1,T2 X3>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator*(
    CT x, const PRODMV<T,T,T X3>& pmv)
{ return PRODMV<CT,T,T X3>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator*(
    CCT x, const PRODMV<T,T,T X3>& pmv)
{ return PRODMV<CT,T,T X3>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator*(
    VCT x, const PRODMV<T,T,T X3>& pmv)
{ return PRODMV<CT,T,T X3>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T2,T1 X3> operator*(T x, const PRODMV<CT,T1,T2 X3>& pmv)
{ return PRODMV<CT,T1,T2 X3>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator*(CCT x, 
    const PRODMV<CT,T1,T2 X3>& pmv)
{ return PRODMV<CT,T1,T2 X3>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator*(VCT x, 
    const PRODMV<CT,T1,T2 X3>& pmv)
{ return PRODMV<CT,T1,T2 X3>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

// (x*m*v)*x
template <class T, class T1, class T2 Y> 
inline PRODMV<T,T1,T2 X3> operator*(const PRODMV<T,T1,T2 X3>& pmv, T x)
{ return PRODMV<T,T1,T2 X3>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator*(
    const PRODMV<T,T,T X3>& pmv, CT x)
{ return PRODMV<CT,T,T X3>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator*(
    const PRODMV<T,T,T X3>& pmv, CCT x)
{ return PRODMV<CT,T,T X3>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator*(
    const PRODMV<T,T,T X3>& pmv, VCT x)
{ return PRODMV<CT,T,T X3>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator*(const PRODMV<CT,T1,T2 X3>& pmv, 
    T x)
{ return PRODMV<CT,T1,T2 X3>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator*(const PRODMV<CT,T1,T2 X3>& pmv, 
    CCT x)
{ return PRODMV<CT,T1,T2 X3>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator*(const PRODMV<CT,T1,T2 X3>& pmv, 
    VCT x)
{ return PRODMV<CT,T1,T2 X3>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

// (x*m*v)/x
template <class T, class T1, class T2 Y> 
inline PRODMV<T,T1,T2 X3> operator/(const PRODMV<T,T1,T2 X3>& pmv, T x)
{ return PRODMV<T,T1,T2 X3>(pmv.GetX()/x,pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator/(
    const PRODMV<T,T,T X3>& pmv, CT x)
{ return PRODMV<CT,T,T X3>(pmv.GetX()/x,pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator/(
    const PRODMV<T,T,T X3>& pmv, CCT x)
{ return PRODMV<CT,T,T X3>(pmv.GetX()/CT(x),pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator/(
    const PRODMV<T,T,T X3>& pmv, VCT x)
{ return PRODMV<CT,T,T X3>(pmv.GetX()/CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator/(const PRODMV<CT,T1,T2 X3>& pmv, 
    T x)
{ return PRODMV<CT,T1,T2 X3>(pmv.GetX()/x,pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator/(const PRODMV<CT,T1,T2 X3>& pmv, 
    CCT x)
{ return PRODMV<CT,T1,T2 X3>(pmv.GetX()/CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator/(const PRODMV<CT,T1,T2 X3>& pmv, 
    VCT x)
{ return PRODMV<CT,T1,T2 X3>(pmv.GetX()/CT(x),pmv.GetM(),pmv.GetV()); }

#ifdef PRODMV_1

// x*(x*m*v)
template <class T, class T1, class T2 Y> 
inline PRODMV<T,T1,T2 X3> operator*(T x, const PRODMV_1<T,T1,T2 X3>& pmv)
{ return PRODMV<T,T1,T2 X3>(x,pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator*(
    CT x, const PRODMV_1<T,T,T X3>& pmv)
{ return PRODMV<CT,T,T X3>(x,pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator*(
    CCT x, const PRODMV_1<T,T,T X3>& pmv)
{ return PRODMV<CT,T,T X3>(CT(x),pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator*(
    VCT x, const PRODMV_1<T,T,T X3>& pmv)
{ return PRODMV<CT,T,T X3>(CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T2,T1 X3> operator*(T x, const PRODMV_1<CT,T1,T2 X3>& pmv)
{ return PRODMV<CT,T1,T2 X3>(CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator*(CCT x, const PRODMV_1<CT,T1,T2 X3>& pmv)
{ return PRODMV<CT,T1,T2 X3>(CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator*(VCT x, const PRODMV_1<CT,T1,T2 X3>& pmv)
{ return PRODMV<CT,T1,T2 X3>(CT(x),pmv.GetM(),pmv.GetV()); }

// (x*m*v)*x
template <class T, class T1, class T2 Y> 
inline PRODMV<T,T1,T2 X3> operator*(const PRODMV_1<T,T1,T2 X3>& pmv, T x)
{ return PRODMV<T,T1,T2 X3>(x,pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator*(
    const PRODMV_1<T,T,T X3>& pmv, CT x)
{ return PRODMV<CT,T,T X3>(x,pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator*(
    const PRODMV_1<T,T,T X3>& pmv, CCT x)
{ return PRODMV<CT,T,T X3>(CT(x),pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator*(
    const PRODMV_1<T,T,T X3>& pmv, VCT x)
{ return PRODMV<CT,T,T X3>(CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator*(const PRODMV_1<CT,T1,T2 X3>& pmv, 
    T x)
{ return PRODMV<CT,T1,T2 X3>(CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator*(const PRODMV_1<CT,T1,T2 X3>& pmv, 
    CCT x)
{ return PRODMV<CT,T1,T2 X3>(CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator*(const PRODMV_1<CT,T1,T2 X3>& pmv, 
    VCT x)
{ return PRODMV<CT,T1,T2 X3>(CT(x),pmv.GetM(),pmv.GetV()); }

// (x*m*v)/x
template <class T, class T1, class T2 Y> 
inline PRODMV<T,T1,T2 X3> operator/(const PRODMV_1<T,T1,T2 X3>& pmv, T x)
{ return PRODMV<T,T1,T2 X3>(RealType(T)(1)/x,pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator/(
    const PRODMV_1<T,T,T X3>& pmv, CT x)
{ return PRODMV<CT,T,T X3>(T(1)/x,pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator/(
    const PRODMV_1<T,T,T X3>& pmv, CCT x)
{ return PRODMV<CT,T,T X3>(T(1)/CT(x),pmv.GetM(),pmv.GetV()); }

template <class T Y> inline PRODMV<CT,T,T X3> operator/(
    const PRODMV_1<T,T,T X3>& pmv, VCT x)
{ return PRODMV<CT,T,T X3>(T(1)/CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator/(const PRODMV_1<CT,T1,T2 X3>& pmv, 
    T x)
{ return PRODMV<CT,T1,T2 X3>(CT(T(1)/x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator/(const PRODMV_1<CT,T1,T2 X3>& pmv, 
    CCT x)
{ return PRODMV<CT,T1,T2 X3>(T(1)/CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, class T2 Y> 
inline PRODMV<CT,T1,T2 X3> operator/(const PRODMV_1<CT,T1,T2 X3>& pmv, 
    VCT x)
{ return PRODMV<CT,T1,T2 X3>(T(1)/CT(x),pmv.GetM(),pmv.GetV()); }

#undef PRODMV_1
#endif

#undef X3
#undef Y
