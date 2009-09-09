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


#ifndef CT
#define CT std::complex<T>
#endif
#ifndef CCT
#define CCT ConjRef<std::complex<T> >
#endif
#ifndef VCT
#define VCT VarConjRef<std::complex<T> >
#endif

// x*(x*m/m)
template <class T, class T1, class T2> 
inline TQUOTMM<T,T1,T2> operator*(T x, const TQUOTMM<T,T1,T2>& qmm)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TQUOTMM<CT,T,T> operator*(CT x, const TQUOTMM<CT,T,T>& qmm)
{ return TQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TQUOTMM<CT,T,T> operator*(CCT x, const TQUOTMM<CT,T,T>& qmm)
{ return TQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TQUOTMM<CT,T,T> operator*(VCT x, const TQUOTMM<CT,T,T>& qmm)
{ return TQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

template <class T, class T1, class T2> 
inline TQUOTMM<CT,T1,T2> operator*(T x, const TQUOTMM<CT,T1,T2>& qmm)
{ return TQUOTMM<CT,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> 
inline TQUOTMM<CT,T1,T2> operator*(CCT x, const TQUOTMM<CT,T1,T2>& qmm)
{ return TQUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> 
inline TQUOTMM<CT,T1,T2> operator*(VCT x, const TQUOTMM<CT,T1,T2>& qmm)
{ return TQUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

// (x*m/m)*x
template <class T, class T1, class T2> 
inline TQUOTMM<T,T1,T2> operator*(const TQUOTMM<T,T1,T2>& qmm, T x)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TQUOTMM<CT,T,T> operator*(const TQUOTMM<CT,T,T>& qmm, CT x)
{ return TQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TQUOTMM<CT,T,T> operator*(const TQUOTMM<CT,T,T>& qmm, CCT x)
{ return TQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TQUOTMM<CT,T,T> operator*(const TQUOTMM<CT,T,T>& qmm, VCT x)
{ return TQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

template <class T, class T1, class T2> 
inline TQUOTMM<CT,T1,T2> operator*(const TQUOTMM<CT,T1,T2>& qmm, T x)
{ return TQUOTMM<CT,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> 
inline TQUOTMM<CT,T1,T2> operator*(const TQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return TQUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> 
inline TQUOTMM<CT,T1,T2> operator*(const TQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return TQUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

// (x*m/m)/x
template <class T, class T1, class T2> 
inline TQUOTMM<T,T1,T2> operator/(const TQUOTMM<T,T1,T2>& qmm, T x)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TQUOTMM<CT,T,T> operator/(const TQUOTMM<CT,T,T>& qmm, CT x)
{ return TQUOTMM<CT,T,T>(qmm.GetX()/x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TQUOTMM<CT,T,T> operator/(const TQUOTMM<CT,T,T>& qmm, CCT x)
{ return TQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TQUOTMM<CT,T,T> operator/(const TQUOTMM<CT,T,T>& qmm, VCT x)
{ return TQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }

template <class T, class T1, class T2> 
inline TQUOTMM<CT,T1,T2> operator/(const TQUOTMM<CT,T1,T2>& qmm, T x)
{ return TQUOTMM<CT,T1,T2>(qmm.GetX()/x,qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> 
inline TQUOTMM<CT,T1,T2> operator/(const TQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return TQUOTMM<CT,T1,T2>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> 
inline TQUOTMM<CT,T1,T2> operator/(const TQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return TQUOTMM<CT,T1,T2>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }

// x*(x*m%m)
template <class T, class T1, class T2> 
inline TRQUOTMM<T,T1,T2> operator*(T x, const TRQUOTMM<T,T1,T2>& qmm)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TRQUOTMM<CT,T,T> operator*(CT x, const TRQUOTMM<CT,T,T>& qmm)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TRQUOTMM<CT,T,T> operator*(CCT x, const TRQUOTMM<CT,T,T>& qmm)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TRQUOTMM<CT,T,T> operator*(VCT x, const TRQUOTMM<CT,T,T>& qmm)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

template <class T, class T1, class T2> 
inline TRQUOTMM<CT,T1,T2> operator*(T x, const TRQUOTMM<CT,T1,T2>& qmm)
{ return TRQUOTMM<CT,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> 
inline TRQUOTMM<CT,T1,T2> operator*(CCT x, const TRQUOTMM<CT,T1,T2>& qmm)
{ return TRQUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> 
inline TRQUOTMM<CT,T1,T2> operator*(VCT x, const TRQUOTMM<CT,T1,T2>& qmm)
{ return TRQUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

// (x*m%m)*x
template <class T, class T1, class T2> 
inline TRQUOTMM<T,T1,T2> operator*(const TRQUOTMM<T,T1,T2>& qmm, T x)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TRQUOTMM<CT,T,T> operator*(const TRQUOTMM<CT,T,T>& qmm, CT x)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TRQUOTMM<CT,T,T> operator*(const TRQUOTMM<CT,T,T>& qmm, CCT x)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TRQUOTMM<CT,T,T> operator*(const TRQUOTMM<CT,T,T>& qmm, VCT x)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

template <class T, class T1, class T2> 
inline TRQUOTMM<CT,T1,T2> operator*(const TRQUOTMM<CT,T1,T2>& qmm, T x)
{ return TRQUOTMM<CT,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> 
inline TRQUOTMM<CT,T1,T2> operator*(const TRQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return TRQUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> 
inline TRQUOTMM<CT,T1,T2> operator*(const TRQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return TRQUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

// (x*m%m)/x
template <class T, class T1, class T2> 
inline TRQUOTMM<T,T1,T2> operator/(const TRQUOTMM<T,T1,T2>& qmm, T x)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TRQUOTMM<CT,T,T> operator/(const TRQUOTMM<CT,T,T>& qmm, CT x)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()/x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TRQUOTMM<CT,T,T> operator/(const TRQUOTMM<CT,T,T>& qmm, CCT x)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T> 
inline TRQUOTMM<CT,T,T> operator/(const TRQUOTMM<CT,T,T>& qmm, VCT x)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }

template <class T, class T1, class T2> 
inline TRQUOTMM<CT,T1,T2> operator/(const TRQUOTMM<CT,T1,T2>& qmm, T x)
{ return TRQUOTMM<CT,T1,T2>(qmm.GetX()/x,qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> 
inline TRQUOTMM<CT,T1,T2> operator/(const TRQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return TRQUOTMM<CT,T1,T2>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> 
inline TRQUOTMM<CT,T1,T2> operator/(const TRQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return TRQUOTMM<CT,T1,T2>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }

