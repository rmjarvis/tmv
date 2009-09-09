///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
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
#define GETM .GetM()
#endif

// x/m
template <class T Y> 
inline QUOTXM<T,T X> operator/(T x, const GENMATRIX<T X>& m)
{ return QUOTXM<T,T X>(x,m); }

template <class T Y> 
inline QUOTXM<CT,T X> operator/(CT x, const GENMATRIX<T X>& m)
{ return QUOTXM<CT,T X>(x,m); }

template <class T Y> 
inline QUOTXM<CT,T X> operator/(CCT x, const GENMATRIX<T X>& m)
{ return QUOTXM<CT,T X>(CT(x),m); }

template <class T Y> 
inline QUOTXM<CT,T X> operator/(VCT x, const GENMATRIX<T X>& m)
{ return QUOTXM<CT,T X>(CT(x),m); }

template <class T Y> 
inline QUOTXM<CT,CT X> operator/(T x, const GENMATRIX<CT X>& m)
{ return QUOTXM<CT,CT X>(CT(x),m); }

template <class T Y> 
inline QUOTXM<CT,CT X> operator/(CCT x, const GENMATRIX<CT X>& m)
{ return QUOTXM<CT,CT X>(CT(x),m); }

template <class T Y> 
inline QUOTXM<CT,CT X> operator/(VCT x, const GENMATRIX<CT X>& m)
{ return QUOTXM<CT,CT X>(CT(x),m); }

// x%m
template <class T Y> 
inline QUOTXM<T,T X> operator%(T x, const GENMATRIX<T X>& m)
{ return QUOTXM<T,T X>(x,m); }

template <class T Y> 
inline QUOTXM<CT,T X> operator%(CT x, const GENMATRIX<T X>& m)
{ return QUOTXM<CT,T X>(x,m); }

template <class T Y> 
inline QUOTXM<CT,T X> operator%(CCT x, const GENMATRIX<T X>& m)
{ return QUOTXM<CT,T X>(CT(x),m); }

template <class T Y> 
inline QUOTXM<CT,T X> operator%(VCT x, const GENMATRIX<T X>& m)
{ return QUOTXM<CT,T X>(CT(x),m); }

template <class T Y> 
inline QUOTXM<CT,CT X> operator%(T x, const GENMATRIX<CT X>& m)
{ return QUOTXM<CT,CT X>(CT(x),m); }

template <class T Y> 
inline QUOTXM<CT,CT X> operator%(CCT x, const GENMATRIX<CT X>& m)
{ return QUOTXM<CT,CT X>(CT(x),m); }

template <class T Y> 
inline QUOTXM<CT,CT X> operator%(VCT x, const GENMATRIX<CT X>& m)
{ return QUOTXM<CT,CT X>(CT(x),m); }


#undef X
#undef Y
#undef GETM
