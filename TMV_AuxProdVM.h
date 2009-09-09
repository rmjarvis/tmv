
// Need to define the following with #define statements.
// (The given definition is for a regular Matrix*Vector.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENVECTOR1 GenVector
// #define GENMATRIX2 GenMatrix
// #define PRODMV ProdMV
// #define PRODXV1 ProdXV
// #define PRODXM2 ProdXM

// m*v
template <class T> inline PRODMV<T,T,false,T> operator*(
    const GENMATRIX2<T>& m, const GENVECTOR1<T>& v)
{ return PRODMV<T,T,false,T>(T(1),m,v); }

template <class T> inline PRODMV<CT,CT,false,T> operator*(
    const GENMATRIX2<CT>& m, const GENVECTOR1<T>& v)
{ return PRODMV<CT,CT,false,T>(CT(1),m,v); }

template <class T> inline PRODMV<CT,T,false,CT> operator*(
    const GENMATRIX2<T>& m, const GENVECTOR1<CT>& v)
{ return PRODMV<CT,T,false,CT>(CT(1),m,v); }

// v*m
template <class T> inline PRODMV<T,T,true,T> operator*(
    const GENVECTOR1<T>& v, const GENMATRIX2<T>& m)
{ return PRODMV<T,T,true,T>(T(1),m,v); }

template <class T> inline PRODMV<CT,T,true,CT> operator*(
    const GENVECTOR1<CT>& v, const GENMATRIX2<T>& m)
{ return PRODMV<CT,T,true,CT>(CT(1),m,v); }

template <class T> inline PRODMV<CT,CT,true,T> operator*(
    const GENVECTOR1<T>& v, const GENMATRIX2<CT>& m)
{ return PRODMV<CT,CT,true,T>(CT(1),m,v); }

// x*(x*m*v)
template <class T, class T1, bool tr, class T2> 
inline PRODMV<T,T1,tr,T2> operator*(T x, const PRODMV<T,T1,tr,T2>& pmv)
{ return PRODMV<T,T1,tr,T2>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    CT x, const PRODMV<T,T,tr,T>& pmv)
{ return PRODMV<CT,T,tr,T>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    CCT x, const PRODMV<T,T,tr,T>& pmv)
{ return PRODMV<CT,T,tr,T>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    VCT x, const PRODMV<T,T,tr,T>& pmv)
{ return PRODMV<CT,T,tr,T>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T2,tr,T1> operator*(T x, const PRODMV<CT,T1,tr,T2>& pmv)
{ return PRODMV<CT,T1,tr,T2>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator*(CCT x, const PRODMV<CT,T1,tr,T2>& pmv)
{ return PRODMV<CT,T1,tr,T2>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator*(VCT x, const PRODMV<CT,T1,tr,T2>& pmv)
{ return PRODMV<CT,T1,tr,T2>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

// (x*m*v)*x
template <class T, class T1, bool tr, class T2> 
inline PRODMV<T,T1,tr,T2> operator*(const PRODMV<T,T1,tr,T2>& pmv, T x)
{ return PRODMV<T,T1,tr,T2>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    const PRODMV<T,T,tr,T>& pmv, CT x)
{ return PRODMV<CT,T,tr,T>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    const PRODMV<T,T,tr,T>& pmv, CCT x)
{ return PRODMV<CT,T,tr,T>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    const PRODMV<T,T,tr,T>& pmv, VCT x)
{ return PRODMV<CT,T,tr,T>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator*(const PRODMV<CT,T1,tr,T2>& pmv, T x)
{ return PRODMV<CT,T1,tr,T2>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator*(const PRODMV<CT,T1,tr,T2>& pmv, CCT x)
{ return PRODMV<CT,T1,tr,T2>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator*(const PRODMV<CT,T1,tr,T2>& pmv, VCT x)
{ return PRODMV<CT,T1,tr,T2>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

// (x*m*v)/x
template <class T, class T1, bool tr, class T2> 
inline PRODMV<T,T1,tr,T2> operator/(const PRODMV<T,T1,tr,T2>& pmv, T x)
{ return PRODMV<T,T1,tr,T2>(pmv.GetX()/x,pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator/(
    const PRODMV<T,T,tr,T>& pmv, CT x)
{ return PRODMV<CT,T,tr,T>(pmv.GetX()/x,pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator/(
    const PRODMV<T,T,tr,T>& pmv, CCT x)
{ return PRODMV<CT,T,tr,T>(pmv.GetX()/CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator/(
    const PRODMV<T,T,tr,T>& pmv, VCT x)
{ return PRODMV<CT,T,tr,T>(pmv.GetX()/CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator/(const PRODMV<CT,T1,tr,T2>& pmv, T x)
{ return PRODMV<CT,T1,tr,T2>(pmv.GetX()/x,pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator/(const PRODMV<CT,T1,tr,T2>& pmv, CCT x)
{ return PRODMV<CT,T1,tr,T2>(pmv.GetX()/CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator/(const PRODMV<CT,T1,tr,T2>& pmv, VCT x)
{ return PRODMV<CT,T1,tr,T2>(pmv.GetX()/CT(x),pmv.GetM(),pmv.GetV()); }

// m*(x*v)
template <class T, class T2> inline PRODMV<T,T,false,T2> operator*(
    const GENMATRIX2<T>& m, const PRODXV1<T,T2>& pxv)
{ return PRODMV<T,T,false,T2>(pxv.GetX(),m,pxv.GetV()); }

template <class T> inline PRODMV<CT,CT,false,T> operator*(
    const GENMATRIX2<CT>& m, const PRODXV1<T,T>& pxv)
{ return PRODMV<CT,CT,false,T>(pxv.GetX(),m,pxv.GetV()); }

template <class T, class T2> inline PRODMV<CT,T,false,T2> operator*(
    const GENMATRIX2<T>& m, const PRODXV1<CT,T2>& pxv)
{ return PRODMV<CT,T,false,T2>(pxv.GetX(),m,pxv.GetV()); }

// (x*v)*m
template <class T, class T2> inline PRODMV<T,T,true,T2> operator*(
    const PRODXV1<T,T2>& pxv, const GENMATRIX2<T>& m)
{ return PRODMV<T,T,true,T2>(pxv.GetX(),m,pxv.GetV()); }

template <class T> inline PRODMV<CT,CT,true,T> operator*(
    const PRODXV1<T,T>& pxv, const GENMATRIX2<CT>& m)
{ return PRODMV<CT,CT,true,T>(pxv.GetX(),m,pxv.GetV()); }

template <class T, class T2> inline PRODMV<CT,T,true,T2> operator*(
    const PRODXV1<CT,T2>& pxv, const GENMATRIX2<T>& m)
{ return PRODMV<CT,T,true,T2>(pxv.GetX(),m,pxv.GetV()); }

// (x*m)*v
template <class T, class T1> inline PRODMV<T,T1,false,T> operator*(
    const PRODXM2<T,T1>& pxm, const GENVECTOR1<T>& v)
{ return PRODMV<T,T1,false,T>(pxm.GetX(),pxm.GetM(),v); }

template <class T> inline PRODMV<CT,T,false,CT> operator*(
    const PRODXM2<T,T>& pxm, const GENVECTOR1<CT>& v)
{ return PRODMV<CT,T,false,CT>(pxm.GetX(),pxm.GetM(),v); }

template <class T, class T1> inline PRODMV<CT,T1,false,T> operator*(
    const PRODXM2<CT,T1>& pxm, const GENVECTOR1<T>& v)
{ return PRODMV<CT,T1,false,T>(pxm.GetX(),pxm.GetM(),v); }

// v*(x*m)
template <class T, class T1> inline PRODMV<T,T1,true,T> operator*(
    const GENVECTOR1<T>& v, const PRODXM2<T,T1>& pxm)
{ return PRODMV<T,T1,true,T>(pxm.GetX(),pxm.GetM(),v); }

template <class T> inline PRODMV<CT,T,true,CT> operator*(
    const GENVECTOR1<CT>& v, const PRODXM2<T,T>& pxm)
{ return PRODMV<CT,T,true,CT>(pxm.GetX(),pxm.GetM(),v); }

template <class T, class T1> inline PRODMV<CT,T1,true,T> operator*(
    const GENVECTOR1<T>& v, const PRODXM2<CT,T1>& pxm)
{ return PRODMV<CT,T1,true,T>(pxm.GetX(),pxm.GetM(),v); }

// (x*m)*(x*v)
template <class T, class T1, class T2> inline PRODMV<T,T1,false,T2> operator*(
    const PRODXM2<T,T1>& pxm, const PRODXV1<T,T2>& pxv) 
{ 
  return PRODMV<T,T1,false,T2>(pxm.GetX()*pxv.GetX(),pxm.GetM(),pxv.GetV()); 
}

template <class T, class T1> inline PRODMV<CT,T1,false,T> operator*(
    const PRODXM2<CT,T1>& pxm, const PRODXV1<T,T>& pxv) 
{ 
  return PRODMV<CT,T1,false,T>(pxm.GetX()*pxv.GetX(),pxm.GetM(),pxv.GetV()); 
}

template <class T, class T2> inline PRODMV<CT,T,false,T2> operator*(
    const PRODXM2<T,T>& pxm, const PRODXV1<CT,T2>& pxv) 
{ 
  return PRODMV<CT,T,false,T2>(pxm.GetX()*pxv.GetX(),pxm.GetM(),pxv.GetV()); 
}

// (x*v)*(x*m)
template <class T, class T1, class T2> inline PRODMV<T,T1,true,T2> operator*(
    const PRODXV1<T,T2>& pxv, const PRODXM2<T,T1>& pxm)
{ 
  return PRODMV<T,T1,true,T2>(pxv.GetX()*pxm.GetX(),pxm.GetM(),pxv.GetV()); 
}

template <class T, class T2> inline PRODMV<CT,T,true,T2> operator*(
    const PRODXV1<CT,T2>& pxv, const PRODXM2<T,T>& pxm)
{ 
  return PRODMV<CT,T,true,T2>(pxv.GetX()*pxm.GetX(),pxm.GetM(),pxv.GetV()); 
}

template <class T, class T1> inline PRODMV<CT,T1,true,T> operator*(
    const PRODXV1<T,T>& pxv, const PRODXM2<CT,T1>& pxm)
{ 
  return PRODMV<CT,T1,true,T>(pxv.GetX()*pxm.GetX(),pxm.GetM(),pxv.GetV()); 
}
