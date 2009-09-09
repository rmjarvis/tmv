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
template <class T, class T1, class T2> inline QUOTMM<T,T1,T2> operator*(
    T x, const QUOTMM<T,T1,T2>& qmm)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    CT x, const QUOTMM<CT,T,T>& qmm)
{ return QUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    CCT x, const QUOTMM<CT,T,T>& qmm)
{ return QUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    VCT x, const QUOTMM<CT,T,T>& qmm)
{ return QUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    T x, const QUOTMM<CT,T1,T2>& qmm)
{ return QUOTMM<CT,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    CCT x, const QUOTMM<CT,T1,T2>& qmm)
{ return QUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    VCT x, const QUOTMM<CT,T1,T2>& qmm)
{ return QUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

// (x*m/m)*x
template <class T, class T1, class T2> inline QUOTMM<T,T1,T2> operator*(
    const QUOTMM<T,T1,T2>& qmm, T x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    const QUOTMM<CT,T,T>& qmm, CT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    const QUOTMM<CT,T,T>& qmm, CCT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    const QUOTMM<CT,T,T>& qmm, VCT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    const QUOTMM<CT,T1,T2>& qmm, T x)
{ return QUOTMM<CT,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    const QUOTMM<CT,T1,T2>& qmm, CCT x)
{ return QUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    const QUOTMM<CT,T1,T2>& qmm, VCT x)
{ return QUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

// (x*m/m)/x
template <class T, class T1, class T2> inline QUOTMM<T,T1,T2> operator/(
    const QUOTMM<T,T1,T2>& qmm, T x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator/(
    const QUOTMM<CT,T,T>& qmm, CT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator/(
    const QUOTMM<CT,T,T>& qmm, CCT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator/(
    const QUOTMM<CT,T,T>& qmm, VCT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator/(
    const QUOTMM<CT,T1,T2>& qmm, T x)
{ return QUOTMM<CT,T1,T2>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator/(
    const QUOTMM<CT,T1,T2>& qmm, CCT x)
{ return QUOTMM<CT,T1,T2>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator/(
    const QUOTMM<CT,T1,T2>& qmm, VCT x)
{ return QUOTMM<CT,T1,T2>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }

// x*(x*m%m)
template <class T, class T1, class T2> inline RQUOTMM<T,T1,T2> operator*(
    T x, const RQUOTMM<T,T1,T2>& qmm)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    CT x, const RQUOTMM<CT,T,T>& qmm)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    CCT x, const RQUOTMM<CT,T,T>& qmm)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    VCT x, const RQUOTMM<CT,T,T>& qmm)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    T x, const RQUOTMM<CT,T1,T2>& qmm)
{ return RQUOTMM<CT,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    CCT x, const RQUOTMM<CT,T1,T2>& qmm)
{ return RQUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    VCT x, const RQUOTMM<CT,T1,T2>& qmm)
{ return RQUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

// (x*m%m)*x
template <class T, class T1, class T2> inline RQUOTMM<T,T1,T2> operator*(
    const RQUOTMM<T,T1,T2>& qmm, T x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    const RQUOTMM<CT,T,T>& qmm, CT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    const RQUOTMM<CT,T,T>& qmm, CCT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    const RQUOTMM<CT,T,T>& qmm, VCT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    const RQUOTMM<CT,T1,T2>& qmm, T x)
{ return RQUOTMM<CT,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    const RQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return RQUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    const RQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return RQUOTMM<CT,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

// (x*m%m)/x
template <class T, class T1, class T2> inline RQUOTMM<T,T1,T2> operator/(
    const RQUOTMM<T,T1,T2>& qmm, T x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator/(
    const RQUOTMM<CT,T,T>& qmm, CT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator/(
    const RQUOTMM<CT,T,T>& qmm, CCT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator/(
    const RQUOTMM<CT,T,T>& qmm, VCT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator/(
    const RQUOTMM<CT,T1,T2>& qmm, T x)
{ return RQUOTMM<CT,T1,T2>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator/(
    const RQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return RQUOTMM<CT,T1,T2>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator/(
    const RQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return RQUOTMM<CT,T1,T2>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }

