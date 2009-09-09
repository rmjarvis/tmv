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



// Need to define the following with #define statements.
// (The given definition is for a regular Matrix*Vector.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENVECTOR GenVector
// #define GENMATRIX GenMatrix
// #define PRODMV ProdMV
// #define PRODXV ProdXV
// #define PRODXM ProdXM

#ifndef X1a
#define X1a
#endif

#ifndef X1b
#define X1b
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

#ifndef PRODMV_1
#define PRODMV_1 PRODMV
#endif

#ifndef PRODMTV
#define PRODMTV PRODMV
#ifndef PRODMTV_1
#define PRODMTV_1 PRODMV
#endif
#else
template <class T, class T1, class T2 Y> class PRODMTV;
#ifndef PRODMTV_1
#define PRODMTV_1 PRODMTV
#endif
#endif

// m*v
template <class T Y> inline PRODMV_1<T,T,T X3> operator*(
    const GENMATRIX<T X1a>& m, const GENVECTOR<T X2>& v)
{ return PRODMV_1<T,T,T X3>(T(1),m,v); }

template <class T Y> inline PRODMV_1<CT,CT,T X3> operator*(
    const GENMATRIX<CT X1a>& m, const GENVECTOR<T X2>& v)
{ return PRODMV_1<CT,CT,T X3>(CT(1),m,v); }

template <class T Y> inline PRODMV_1<CT,T,CT X3> operator*(
    const GENMATRIX<T X1a>& m, const GENVECTOR<CT X2>& v)
{ return PRODMV_1<CT,T,CT X3>(CT(1),m,v); }

// v*m
template <class T Y> inline PRODMTV_1<T,T,T X3> operator*(
    const GENVECTOR<T X2>& v, const GENMATRIX<T X1b>& m)
{ return PRODMTV_1<T,T,T X3>(T(1),m.Transpose(),v); }

template <class T Y> inline PRODMTV_1<CT,T,CT X3> operator*(
    const GENVECTOR<CT X2>& v, const GENMATRIX<T X1b>& m)
{ return PRODMTV_1<CT,T,CT X3>(CT(1),m.Transpose(),v); }

template <class T Y> inline PRODMTV_1<CT,CT,T X3> operator*(
    const GENVECTOR<T X2>& v, const GENMATRIX<CT X1b>& m)
{ return PRODMTV_1<CT,CT,T X3>(CT(1),m.Transpose(),v); }

// m*(x*v)
template <class T, class T2 Y> inline PRODMV<T,T,T2 X3> operator*(
    const GENMATRIX<T X1a>& m, const PRODXV<T,T2 X2>& pxv)
{ return PRODMV<T,T,T2 X3>(pxv.GetX(),m,pxv.GetV()); }

template <class T Y> inline PRODMV<CT,CT,T X3> operator*(
    const GENMATRIX<CT X1a>& m, const PRODXV<T,T X2>& pxv)
{ return PRODMV<CT,CT,T X3>(CT(pxv.GetX()),m,pxv.GetV()); }

template <class T, class T2 Y> inline PRODMV<CT,T,T2 X3> operator*(
    const GENMATRIX<T X1a>& m, const PRODXV<CT,T2 X2>& pxv)
{ return PRODMV<CT,T,T2 X3>(pxv.GetX(),m,pxv.GetV()); }

// (x*v)*m
template <class T, class T2 Y> inline PRODMTV<T,T,T2 X3> operator*(
    const PRODXV<T,T2 X2>& pxv, const GENMATRIX<T X1b>& m)
{ return PRODMTV<T,T,T2 X3>(pxv.GetX(),m.Transpose(),pxv.GetV()); }

template <class T Y> inline PRODMTV<CT,CT,T X3> operator*(
    const PRODXV<T,T X2>& pxv, const GENMATRIX<CT X1b>& m)
{ return PRODMTV<CT,CT,T X3>(CT(pxv.GetX()),m.Transpose(),pxv.GetV()); }

template <class T, class T2 Y> inline PRODMTV<CT,T,T2 X3> operator*(
    const PRODXV<CT,T2 X2>& pxv, const GENMATRIX<T X1b>& m)
{ return PRODMTV<CT,T,T2 X3>(pxv.GetX(),m.Transpose(),pxv.GetV()); }

// (x*m)*v
template <class T, class T1 Y> inline PRODMV<T,T1,T X3> operator*(
    const PRODXM<T,T1 X1a>& pxm, const GENVECTOR<T X2>& v)
{ return PRODMV<T,T1,T X3>(pxm.GetX(),pxm.GetM(),v); }

template <class T Y> inline PRODMV<CT,T,CT X3> operator*(
    const PRODXM<T,T X1a>& pxm, const GENVECTOR<CT X2>& v)
{ return PRODMV<CT,T,CT X3>(CT(pxm.GetX()),pxm.GetM(),v); }

template <class T, class T1 Y> inline PRODMV<CT,T1,T X3> operator*(
    const PRODXM<CT,T1 X1a>& pxm, const GENVECTOR<T X2>& v)
{ return PRODMV<CT,T1,T X3>(pxm.GetX(),pxm.GetM(),v); }

// v*(x*m)
template <class T, class T1 Y> inline PRODMTV<T,T1,T X3> operator*(
    const GENVECTOR<T X2>& v, const PRODXM<T,T1 X1b>& pxm)
{ return PRODMTV<T,T1,T X3>(pxm.GetX(),pxm.GetM().Transpose(),v); }

template <class T Y> inline PRODMTV<CT,T,CT X3> operator*(
    const GENVECTOR<CT X2>& v, const PRODXM<T,T X1b>& pxm)
{ return PRODMTV<CT,T,CT X3>(CT(pxm.GetX()),pxm.GetM().Transpose(),v); }

template <class T, class T1 Y> inline PRODMTV<CT,T1,T X3> operator*(
    const GENVECTOR<T X2>& v, const PRODXM<CT,T1 X1b>& pxm)
{ return PRODMTV<CT,T1,T X3>(pxm.GetX(),pxm.GetM().Transpose(),v); }

// (x*m)*(x*v)
template <class T, class T1, class T2 Y> 
inline PRODMV<T,T1,T2 X3> operator*(
    const PRODXM<T,T1 X1a>& pxm, const PRODXV<T,T2 X2>& pxv) 
{ 
  return PRODMV<T,T1,T2 X3>(pxm.GetX()*pxv.GetX(),pxm.GetM(),pxv.GetV()); 
}

template <class T, class T1 Y> inline PRODMV<CT,T1,T X3> operator*(
    const PRODXM<CT,T1 X1a>& pxm, const PRODXV<T,T X2>& pxv) 
{ 
  return PRODMV<CT,T1,T X3>(pxm.GetX()*pxv.GetX(),pxm.GetM(),pxv.GetV()); 
}

template <class T, class T2 Y> inline PRODMV<CT,T,T2 X3> operator*(
    const PRODXM<T,T X1a>& pxm, const PRODXV<CT,T2 X2>& pxv) 
{ 
  return PRODMV<CT,T,T2 X3>(pxm.GetX()*pxv.GetX(),pxm.GetM(),pxv.GetV()); 
}

// (x*v)*(x*m)
template <class T, class T1, class T2 Y> inline PRODMTV<T,T1,T2 X3> operator*(
    const PRODXV<T,T2 X2>& pxv, const PRODXM<T,T1 X1b>& pxm)
{ 
  return PRODMTV<T,T1,T2 X3>(pxv.GetX()*pxm.GetX(),pxm.GetM().Transpose(),
      pxv.GetV()); 
}

template <class T, class T2 Y> inline PRODMTV<CT,T,T2 X3> operator*(
    const PRODXV<CT,T2 X2>& pxv, const PRODXM<T,T X1b>& pxm)
{ 
  return PRODMTV<CT,T,T2 X3>(pxv.GetX()*pxm.GetX(),pxm.GetM().Transpose(),
      pxv.GetV()); 
}

template <class T, class T1 Y> inline PRODMTV<CT,T1,T X3> operator*(
    const PRODXV<T,T X2>& pxv, const PRODXM<CT,T1 X1b>& pxm)
{ 
  return PRODMTV<CT,T1,T X3>(pxv.GetX()*pxm.GetX(),pxm.GetM().Transpose(),
      pxv.GetV()); 
}

#undef X1a
#undef X1b
#undef X2
#undef X3
#undef Y
#undef PRODMTV
#undef PRODMV_1
#undef PRODMTV_1
