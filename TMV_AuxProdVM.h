
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
    const GENMATRIX2<T>& m1, const GENVECTOR1<T>& v2)
{ return PRODMV<T,T,false,T>(T(1),m1,v2); }

template <class T> inline PRODMV<CT,CT,false,T> operator*(
    const GENMATRIX2<CT>& m1, const GENVECTOR1<T>& v2)
{ return PRODMV<CT,CT,false,T>(CT(1),m1,v2); }

template <class T> inline PRODMV<CT,T,false,CT> operator*(
    const GENMATRIX2<T>& m1, const GENVECTOR1<CT>& v2)
{ return PRODMV<CT,T,false,CT>(CT(1),m1,v2); }

// v*m
template <class T> inline PRODMV<T,T,true,T> operator*(
    const GENVECTOR1<T>& v1, const GENMATRIX2<T>& m2)
{ return PRODMV<T,T,true,T>(T(1),m2,v1); }

template <class T> inline PRODMV<CT,T,true,CT> operator*(
    const GENVECTOR1<CT>& v1, const GENMATRIX2<T>& m2)
{ return PRODMV<CT,T,true,CT>(CT(1),m2,v1); }

template <class T> inline PRODMV<CT,CT,true,T> operator*(
    const GENVECTOR1<T>& v1, const GENMATRIX2<CT>& m2)
{ return PRODMV<CT,CT,true,T>(CT(1),m2,v1); }

// x*(x*m*v)
template <class T, class T2, bool tr, class T1> 
inline PRODMV<T,T2,tr,T1> operator*(
    T x1, const PRODMV<T,T2,tr,T1>& pmv2)
{ return PRODMV<T,T2,tr,T1>(x1*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetV3()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    CT x1, const PRODMV<T,T,tr,T>& pmv2)
{ return PRODMV<CT,T,tr,T>(x1*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetV3()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    CCT x1, const PRODMV<T,T,tr,T>& pmv2)
{ return PRODMV<CT,T,tr,T>(CT(x1)*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetV3()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    VCT x1, const PRODMV<T,T,tr,T>& pmv2)
{ return PRODMV<CT,T,tr,T>(CT(x1)*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetV3()); }

template <class T, class T2, bool tr, class T1> 
inline PRODMV<CT,T2,tr,T1> operator*(
    T x1, const PRODMV<CT,T2,tr,T1>& pmv2)
{ return PRODMV<CT,T2,tr,T1>(x1*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetV3()); }

template <class T, class T2, bool tr, class T1> 
inline PRODMV<CT,T2,tr,T1> operator*(
    CCT x1, const PRODMV<CT,T2,tr,T1>& pmv2)
{ return PRODMV<CT,T2,tr,T1>(CT(x1)*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetV3()); }

template <class T, class T2, bool tr, class T1> 
inline PRODMV<CT,T2,tr,T1> operator*(
    VCT x1, const PRODMV<CT,T2,tr,T1>& pmv2)
{ return PRODMV<CT,T2,tr,T1>(CT(x1)*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetV3()); }

// (x*m*v)*x
template <class T, class T2, bool tr, class T1> 
inline PRODMV<T,T2,tr,T1> operator*(
    const PRODMV<T,T2,tr,T1>& pmv1, T x2)
{ return PRODMV<T,T2,tr,T1>(x2*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetV3()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    const PRODMV<T,T,tr,T>& pmv1, CT x2)
{ return PRODMV<CT,T,tr,T>(x2*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetV3()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    const PRODMV<T,T,tr,T>& pmv1, CCT x2)
{ return PRODMV<CT,T,tr,T>(CT(x2)*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetV3()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    const PRODMV<T,T,tr,T>& pmv1, VCT x2)
{ return PRODMV<CT,T,tr,T>(CT(x2)*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetV3()); }

template <class T, class T2, bool tr, class T1> 
inline PRODMV<CT,T2,tr,T1> operator*(
    const PRODMV<CT,T2,tr,T1>& pmv1, T x2)
{ return PRODMV<CT,T2,tr,T1>(x2*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetV3()); }

template <class T, class T2, bool tr, class T1> 
inline PRODMV<CT,T2,tr,T1> operator*(
    const PRODMV<CT,T2,tr,T1>& pmv1, CCT x2)
{ return PRODMV<CT,T2,tr,T1>(CT(x2)*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetV3()); }

template <class T, class T2, bool tr, class T1> 
inline PRODMV<CT,T2,tr,T1> operator*(
    const PRODMV<CT,T2,tr,T1>& pmv1, VCT x2)
{ return PRODMV<CT,T2,tr,T1>(CT(x2)*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetV3()); }

// m*(x*v)
template <class T, class T1> inline PRODMV<T,T,false,T1> operator*(
    const GENMATRIX2<T>& m1, const PRODXV1<T,T1>& pxv2)
{ return PRODMV<T,T,false,T1>(pxv2.GetX1(),m1,pxv2.GetV2()); }

template <class T> inline PRODMV<CT,CT,false,T> operator*(
    const GENMATRIX2<CT>& m1, const PRODXV1<T,T>& pxv2)
{ return PRODMV<CT,CT,false,T>(pxv2.GetX1(),m1,pxv2.GetV2()); }

template <class T, class T1> inline PRODMV<CT,T,false,T1> operator*(
    const GENMATRIX2<T>& m1, const PRODXV1<CT,T1>& pxv2)
{ return PRODMV<CT,T,false,T1>(pxv2.GetX1(),m1,pxv2.GetV2()); }

// (x*v)*m
template <class T, class T1> inline PRODMV<T,T,true,T1> operator*(
    const PRODXV1<T,T1>& pxv1, const GENMATRIX2<T>& m2)
{ return PRODMV<T,T,true,T1>(pxv1.GetX1(),m2,pxv1.GetV2()); }

template <class T> inline PRODMV<CT,CT,true,T> operator*(
    const PRODXV1<T,T>& pxv1, const GENMATRIX2<CT>& m2)
{ return PRODMV<CT,CT,true,T>(pxv1.GetX1(),m2,pxv1.GetV2()); }

template <class T, class T1> inline PRODMV<CT,T,true,T1> operator*(
    const PRODXV1<CT,T1>& pxv1, const GENMATRIX2<T>& m2)
{ return PRODMV<CT,T,true,T1>(pxv1.GetX1(),m2,pxv1.GetV2()); }

// (x*m)*v
template <class T, class T2> inline PRODMV<T,T2,false,T> operator*(
    const PRODXM2<T,T2>& pxm1, const GENVECTOR1<T>& v2)
{ return PRODMV<T,T2,false,T>(pxm1.GetX1(),pxm1.GetM2(),v2); }

template <class T> inline PRODMV<CT,T,false,CT> operator*(
    const PRODXM2<T,T>& pxm1, const GENVECTOR1<CT>& v2)
{ return PRODMV<CT,T,false,CT>(pxm1.GetX1(),pxm1.GetM2(),v2); }

template <class T, class T2> inline PRODMV<CT,T2,false,T> operator*(
    const PRODXM2<CT,T2>& pxm1, const GENVECTOR1<T>& v2)
{ return PRODMV<CT,T2,false,T>(pxm1.GetX1(),pxm1.GetM2(),v2); }

// v*(x*m)
template <class T, class T2> inline PRODMV<T,T2,true,T> operator*(
    const GENVECTOR1<T>& v1, const PRODXM2<T,T2>& pxm2)
{ return PRODMV<T,T2,true,T>(pxm2.GetX1(),pxm2.GetM2(),v1); }

template <class T> inline PRODMV<CT,T,true,CT> operator*(
    const GENVECTOR1<CT>& v1, const PRODXM2<T,T>& pxm2)
{ return PRODMV<CT,T,true,CT>(pxm2.GetX1(),pxm2.GetM2(),v1); }

template <class T, class T2> inline PRODMV<CT,T2,true,T> operator*(
    const GENVECTOR1<T>& v1, const PRODXM2<CT,T2>& pxm2)
{ return PRODMV<CT,T2,true,T>(pxm2.GetX1(),pxm2.GetM2(),v1); }

// (x*m)*(x*v)
template <class T, class T2, class T1> inline PRODMV<T,T2,false,T1> operator*(
    const PRODXM2<T,T2>& pxm1, const PRODXV1<T,T1>& pxv2) 
{ 
  return PRODMV<T,T2,false,T1>(pxm1.GetX1()*pxv2.GetX1(),
      pxm1.GetM2(),pxv2.GetV2()); 
}

template <class T, class T2> inline PRODMV<CT,T2,false,T> operator*(
    const PRODXM2<CT,T2>& pxm1, const PRODXV1<T,T>& pxv2) 
{ 
  return PRODMV<CT,T2,false,T>(pxm1.GetX1()*pxv2.GetX1(),
      pxm1.GetM2(),pxv2.GetV2()); 
}

template <class T, class T1> inline PRODMV<CT,T,false,T1> operator*(
    const PRODXM2<T,T>& pxm1, const PRODXV1<CT,T1>& pxv2) 
{ 
  return PRODMV<CT,T,false,T1>(pxm1.GetX1()*pxv2.GetX1(),
      pxm1.GetM2(),pxv2.GetV2()); 
}

// (x*v)*(x*m)
template <class T, class T1, class T2> inline PRODMV<T,T,true,T> operator*(
    const PRODXV1<T,T1>& pxv1, const PRODXM2<T,T2>& pxm2)
{ 
  return PRODMV<T,T2,true,T1>(pxv1.GetX1()*pxm2.GetX1(),
      pxm2.GetM2(),pxv1.GetV2()); 
}

template <class T, class T1> inline PRODMV<CT,T,true,T1> operator*(
    const PRODXV1<CT,T1>& pxv1, const PRODXM2<T,T>& pxm2)
{ 
  return PRODMV<CT,T,true,T1>(pxv1.GetX1()*pxm2.GetX1(),
      pxm2.GetM2(),pxv1.GetV2()); 
}

template <class T, class T2> inline PRODMV<CT,T2,true,T> operator*(
    const PRODXV1<T,T>& pxv1, const PRODXM2<CT,T2>& pxm2)
{ 
  return PRODMV<CT,T2,true,T>(pxv1.GetX1()*pxm2.GetX1(),
      pxm2.GetM2(),pxv1.GetV2()); 
}
