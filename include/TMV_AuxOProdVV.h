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


// Things that need to be #defined on entry:
// (The values for a normal Vector^Vector are given)
//
// OPRODVV	OProdVV
// PRODXV1	ProdXV
// PRODXV2	ProdXV
// GENVECTOR1	GenVector
// GENVECTOR2	GenVector

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

#ifndef CT
#define CT std::complex<T>
#endif
#ifndef CCT
#define CCT ConjRef<std::complex<T> >
#endif
#ifndef VCT
#define VCT VarConjRef<std::complex<T> >
#endif

#ifndef OPRODVV_1
#define OPRODVV_1 OPRODVV
#define NO_1
#endif

// v^v
template <class T Y> inline OPRODVV_1<T,T,T X3> operator^(
    const GENVECTOR1<T X1>& v1, const GENVECTOR2<T X2>& v2)
{ return OPRODVV_1<T,T,T X3>(T(1),v1,v2); }

template <class T Y> inline OPRODVV_1<CT,CT,T X3> operator^(
    const GENVECTOR1<CT X1>& v1, const GENVECTOR2<T X2>& v2)
{ return OPRODVV_1<CT,CT,T X3>(CT(1),v1,v2); }

template <class T Y> inline OPRODVV_1<CT,T,CT X3> operator^(
    const GENVECTOR1<T X1>& v1, const GENVECTOR2<CT X2>& v2)
{ return OPRODVV_1<CT,T,CT X3>(CT(1),v1,v2); }

// x*(x*v^v)
template <class T, class T1, class T2 Y> inline OPRODVV<T,T1,T2 X3> operator*(
    T x, const OPRODVV<T,T1,T2 X3>& opvv)
{ return OPRODVV<T,T1,T2 X3>(x*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator*(
    CT x, const OPRODVV<T,T,T X3>& opvv)
{ return OPRODVV<CT,T,T X3>(x*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator*(CCT x,
    const OPRODVV<T,T,T X3>& opvv)
{ return OPRODVV<CT,T,T X3>(CT(x)*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator*(VCT x,
    const OPRODVV<T,T,T X3>& opvv)
{ return OPRODVV<CT,T,T X3>(CT(x)*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator*(
    T x, const OPRODVV<CT,T1,T2 X3>& opvv)
{ return OPRODVV<CT,T1,T2 X3>(x*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator*(
    CCT x, const OPRODVV<CT,T1,T2 X3>& opvv)
{ return OPRODVV<CT,T1,T2 X3>(CT(x)*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator*(
    VCT x, const OPRODVV<CT,T1,T2 X3>& opvv)
{ return OPRODVV<CT,T1,T2 X3>(CT(x)*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

// (x*v^v)*x
template <class T, class T1, class T2 Y> inline OPRODVV<T,T1,T2 X3> operator*(
    const OPRODVV<T,T1,T2 X3>& opvv, T x)
{ return OPRODVV<T,T1,T2 X3>(x*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator*(
    const OPRODVV<T,T,T X3>& opvv, CT x)
{ return OPRODVV<CT,T,T X3>(x*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator*(
    const OPRODVV<T,T,T X3>& opvv, CCT x)
{ return OPRODVV<CT,T,T X3>(CT(x)*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator*(
    const OPRODVV<T,T,T X3>& opvv, VCT x)
{ return OPRODVV<CT,T,T X3>(CT(x)*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator*(
    const OPRODVV<CT,T1,T2 X3>& opvv, T x)
{ return OPRODVV<CT,T1,T2 X3>(x*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator*(
    const OPRODVV<CT,T1,T2 X3>& opvv, CCT x)
{ return OPRODVV<CT,T1,T2 X3>(CT(x)*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator*(
    const OPRODVV<CT,T1,T2 X3>& opvv, VCT x)
{ return OPRODVV<CT,T1,T2 X3>(CT(x)*opvv.GetX(),opvv.GetV1(),opvv.GetV2()); }

// (x*v^v)/x
template <class T, class T1, class T2 Y> inline OPRODVV<T,T1,T2 X3> operator/(
    const OPRODVV<T,T1,T2 X3>& opvv, T x)
{ return OPRODVV<T,T1,T2 X3>(opvv.GetX()/x,opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator/(
    const OPRODVV<T,T,T X3>& opvv, CT x)
{ return OPRODVV<CT,T,T X3>(opvv.GetX()/x,opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator/(
    const OPRODVV<T,T,T X3>& opvv, CCT x)
{ return OPRODVV<CT,T,T X3>(opvv.GetX()/CT(x),opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator/(
    const OPRODVV<T,T,T X3>& opvv, VCT x)
{ return OPRODVV<CT,T,T X3>(opvv.GetX()/CT(x),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator/(
    const OPRODVV<CT,T1,T2 X3>& opvv, T x)
{ return OPRODVV<CT,T1,T2 X3>(opvv.GetX()/x,opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator/(
    const OPRODVV<CT,T1,T2 X3>& opvv, CCT x)
{ return OPRODVV<CT,T1,T2 X3>(opvv.GetX()/CT(x),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator/(
    const OPRODVV<CT,T1,T2 X3>& opvv, VCT x)
{ return OPRODVV<CT,T1,T2 X3>(opvv.GetX()/CT(x),opvv.GetV1(),opvv.GetV2()); }

// (x*v)^v
template <class T, class T1 Y> inline OPRODVV<T,T1,T X3> operator^(
    const PRODXV1<T,T1 X1>& pxv, const GENVECTOR2<T X2>& v)
{ return OPRODVV<T,T1,T X3>(pxv.GetX(),pxv.GetV(),v); }

template <class T, class T1 Y> inline OPRODVV<CT,T1,T X3> operator^(
    const PRODXV1<CT,T1 X1>& pxv, const GENVECTOR2<T X2>& v)
{ return OPRODVV<CT,T1,T X3>(pxv.GetX(),pxv.GetV(),v); }

template <class T Y> inline OPRODVV<CT,T,CT X3> operator^(
    const PRODXV1<T,T X1>& pxv, const GENVECTOR2<CT X2>& v)
{ return OPRODVV<CT,T,CT X3>(CT(pxv.GetX()),pxv.GetV(),v); }

// v^(x*v)
template <class T, class T2 Y> inline OPRODVV<T,T,T2 X3> operator^(
    const GENVECTOR1<T X1>& v, const PRODXV2<T,T2 X2>& pxv)
{ return OPRODVV<T,T,T2 X3>(pxv.GetX(),v,pxv.GetV()); }

template <class T Y> inline OPRODVV<CT,CT,T X3> operator^(
    const GENVECTOR1<CT X1>& v, const PRODXV2<T,T X2>& pxv)
{ return OPRODVV<CT,CT,T X3>(CT(pxv.GetX()),v,pxv.GetV()); }

template <class T, class T2 Y> inline OPRODVV<CT,T,T2 X3> operator^(
    const GENVECTOR1<T X1>& v, const PRODXV2<CT,T2 X2>& pxv)
{ return OPRODVV<CT,T,T2 X3>(pxv.GetX(),v,pxv.GetV()); }

// (x*v)^(x*v)
template <class T, class T1, class T2 Y> inline OPRODVV<T,T1,T2 X3> operator^(
    const PRODXV1<T,T1 X1>& pxv1, const PRODXV2<T,T2 X2>& pxv2)
{
  return OPRODVV<T,T1,T2 X3>(pxv1.GetX()*pxv2.GetX(),pxv1.GetV(),pxv2.GetV()); 
}

template <class T, class T1 Y> inline OPRODVV<CT,T1,T X3> operator^(
    const PRODXV1<CT,T1 X1>& pxv1, const PRODXV2<T,T X2>& pxv2)
{ 
  return OPRODVV<CT,T1,T X3>(pxv1.GetX()*pxv2.GetX(),pxv1.GetV(),pxv2.GetV()); 
}

template <class T, class T2 Y> inline OPRODVV<CT,T,T2 X3> operator^(
    const PRODXV1<T,T X1>& pxv1, const PRODXV2<CT,T2 X2>& pxv2)
{ 
  return OPRODVV<CT,T,T2 X3>(pxv1.GetX()*pxv2.GetX(),pxv1.GetV(),pxv2.GetV()); 
}

#ifndef NO_1
// x*(x*v^v)
template <class T, class T1, class T2 Y> inline OPRODVV<T,T1,T2 X3> operator*(
    T x, const OPRODVV_1<T,T1,T2 X3>& opvv)
{ return OPRODVV<T,T1,T2 X3>(x,opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator*(
    CT x, const OPRODVV_1<T,T,T X3>& opvv)
{ return OPRODVV<CT,T,T X3>(x,opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator*(CCT x,
    const OPRODVV_1<T,T,T X3>& opvv)
{ return OPRODVV<CT,T,T X3>(CT(x),opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator*(VCT x,
    const OPRODVV_1<T,T,T X3>& opvv)
{ return OPRODVV<CT,T,T X3>(CT(x),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator*(
    T x, const OPRODVV_1<CT,T1,T2 X3>& opvv)
{ return OPRODVV<CT,T1,T2 X3>(x,opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator*(
    CCT x, const OPRODVV_1<CT,T1,T2 X3>& opvv)
{ return OPRODVV<CT,T1,T2 X3>(CT(x),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator*(
    VCT x, const OPRODVV_1<CT,T1,T2 X3>& opvv)
{ return OPRODVV<CT,T1,T2 X3>(CT(x),opvv.GetV1(),opvv.GetV2()); }

// (x*v^v)*x
template <class T, class T1, class T2 Y> inline OPRODVV<T,T1,T2 X3> operator*(
    const OPRODVV_1<T,T1,T2 X3>& opvv, T x)
{ return OPRODVV<T,T1,T2 X3>(x,opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator*(
    const OPRODVV_1<T,T,T X3>& opvv, CT x)
{ return OPRODVV<CT,T,T X3>(x,opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator*(
    const OPRODVV_1<T,T,T X3>& opvv, CCT x)
{ return OPRODVV<CT,T,T X3>(CT(x),opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator*(
    const OPRODVV_1<T,T,T X3>& opvv, VCT x)
{ return OPRODVV<CT,T,T X3>(CT(x),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator*(
    const OPRODVV_1<CT,T1,T2 X3>& opvv, T x)
{ return OPRODVV<CT,T1,T2 X3>(x,opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator*(
    const OPRODVV_1<CT,T1,T2 X3>& opvv, CCT x)
{ return OPRODVV<CT,T1,T2 X3>(CT(x),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator*(
    const OPRODVV_1<CT,T1,T2 X3>& opvv, VCT x)
{ return OPRODVV<CT,T1,T2 X3>(CT(x),opvv.GetV1(),opvv.GetV2()); }

// (x*v^v)/x
template <class T, class T1, class T2 Y> inline OPRODVV<T,T1,T2 X3> operator/(
    const OPRODVV_1<T,T1,T2 X3>& opvv, T x)
{ return OPRODVV<T,T1,T2 X3>(RealType(T)(1)/x,opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator/(
    const OPRODVV_1<T,T,T X3>& opvv, CT x)
{ return OPRODVV<CT,T,T X3>(T(1)/x,opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator/(
    const OPRODVV_1<T,T,T X3>& opvv, CCT x)
{ return OPRODVV<CT,T,T X3>(T(1)/CT(x),opvv.GetV1(),opvv.GetV2()); }

template <class T Y> inline OPRODVV<CT,T,T X3> operator/(
    const OPRODVV_1<T,T,T X3>& opvv, VCT x)
{ return OPRODVV<CT,T,T X3>(T(1)/CT(x),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator/(
    const OPRODVV_1<CT,T1,T2 X3>& opvv, T x)
{ return OPRODVV<CT,T1,T2 X3>(T(1)/x,opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator/(
    const OPRODVV_1<CT,T1,T2 X3>& opvv, CCT x)
{ return OPRODVV<CT,T1,T2 X3>(T(1)/CT(x),opvv.GetV1(),opvv.GetV2()); }

template <class T, class T1, class T2 Y> inline OPRODVV<CT,T1,T2 X3> operator/(
    const OPRODVV_1<CT,T1,T2 X3>& opvv, VCT x)
{ return OPRODVV<CT,T1,T2 X3>(T(1)/CT(x),opvv.GetV1(),opvv.GetV2()); }
#endif

#undef X1
#undef X2
#undef X3
#undef Y
#undef OPRODVV_1
#ifdef NO_1
#undef NO_1
#endif
