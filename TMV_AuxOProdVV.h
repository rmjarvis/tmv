// Things that need to be #defined on entry:
// (The values for a normal Vector^Vector are given)
//
// OPRODVV	OProdVV
// PRODXV1	ProdXV
// PRODXV2	ProdXV
// GENVECTOR1	GenVector
// GENVECTOR2	GenVector

// v^v
template <class T> inline OPRODVV<T,T,T> operator^(
    const GENVECTOR1<T>& v1, const GENVECTOR2<T>& v2)
{ return OPRODVV<T,T,T>(T(1),v1,v2); }

template <class T> inline OPRODVV<CT,CT,T> operator^(
    const GENVECTOR1<CT>& v1, const GENVECTOR2<T>& v2)
{ return OPRODVV<CT,CT,T>(CT(1),v1,v2); }

template <class T> inline OPRODVV<CT,T,CT> operator^(
    const GENVECTOR1<T>& v1, const GENVECTOR2<CT>& v2)
{ return OPRODVV<CT,T,CT>(CT(1),v1,v2); }

// x*(x*v^v)
template <class T, class T1, class T2> 
inline OPRODVV<T,T1,T2> operator*(
    T x, const OPRODVV<T,T1,T2>& opvv)
{ return OPRODVV<T,T1,T2>(x*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

template <class T> inline OPRODVV<CT,T,T> operator*(
    CT x, const OPRODVV<T,T,T>& opvv)
{ return OPRODVV<CT,T,T>(x*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

template <class T> inline OPRODVV<CT,T,T> operator*(CCT x,
    const OPRODVV<T,T,T>& opvv)
{ return OPRODVV<CT,T,T>(CT(x)*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

template <class T> inline OPRODVV<CT,T,T> operator*(VCT x,
    const OPRODVV<T,T,T>& opvv)
{ return OPRODVV<CT,T,T>(CT(x)*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

template <class T, class T1, class T2> inline OPRODVV<CT,T1,T2> operator*(
    T x, const OPRODVV<CT,T1,T2>& opvv)
{ return OPRODVV<CT,T1,T2>(x*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

template <class T, class T1, class T2> inline OPRODVV<CT,T1,T2> operator*(
    CCT x, const OPRODVV<CT,T1,T2>& opvv)
{ return OPRODVV<CT,T1,T2>(CT(x)*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

template <class T, class T1, class T2> inline OPRODVV<CT,T1,T2> operator*(
    VCT x, const OPRODVV<CT,T1,T2>& opvv)
{ return OPRODVV<CT,T1,T2>(CT(x)*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

// (x*v^v)*x
template <class T, class T1, class T2> inline OPRODVV<T,T1,T2> operator*(
    const OPRODVV<T,T1,T2>& opvv, T x)
{ return OPRODVV<T,T1,T2>(x*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

template <class T> inline OPRODVV<CT,T,T> operator*(
    const OPRODVV<T,T,T>& opvv, CT x)
{ return OPRODVV<CT,T,T>(x*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

template <class T> inline OPRODVV<CT,T,T> operator*(
    const OPRODVV<T,T,T>& opvv, CCT x)
{ return OPRODVV<CT,T,T>(CT(x)*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

template <class T> inline OPRODVV<CT,T,T> operator*(
    const OPRODVV<T,T,T>& opvv, VCT x)
{ return OPRODVV<CT,T,T>(CT(x)*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

template <class T, class T1, class T2> inline OPRODVV<CT,T1,T2> operator*(
    const OPRODVV<CT,T1,T2>& opvv, T x)
{ return OPRODVV<CT,T1,T2>(x*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

template <class T, class T1, class T2> inline OPRODVV<CT,T1,T2> operator*(
    const OPRODVV<CT,T1,T2>& opvv, CCT x)
{ return OPRODVV<CT,T1,T2>(CT(x)*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

template <class T, class T1, class T2> inline OPRODVV<CT,T1,T2> operator*(
    const OPRODVV<CT,T1,T2>& opvv, VCT x)
{ return OPRODVV<CT,T1,T2>(CT(x)*opvv.GetX1(),opvv.GetV2(),opvv.GetV3()); }

// (x*v^v)/x
template <class T, class T1, class T2> inline OPRODVV<T,T1,T2> operator/(
    const OPRODVV<T,T1,T2>& opvv, T x)
{ return OPRODVV<T,T1,T2>(opvv.GetX1()/x,opvv.GetV2(),opvv.GetV3()); }

template <class T> inline OPRODVV<CT,T,T> operator/(
    const OPRODVV<T,T,T>& opvv, CT x)
{ return OPRODVV<CT,T,T>(opvv.GetX1()/x,opvv.GetV2(),opvv.GetV3()); }

template <class T> inline OPRODVV<CT,T,T> operator/(
    const OPRODVV<T,T,T>& opvv, CCT x)
{ return OPRODVV<CT,T,T>(opvv.GetX1()/CT(x),opvv.GetV2(),opvv.GetV3()); }

template <class T> inline OPRODVV<CT,T,T> operator/(
    const OPRODVV<T,T,T>& opvv, VCT x)
{ return OPRODVV<CT,T,T>(opvv.GetX1()/CT(x),opvv.GetV2(),opvv.GetV3()); }

template <class T, class T1, class T2> inline OPRODVV<CT,T1,T2> operator/(
    const OPRODVV<CT,T1,T2>& opvv, T x)
{ return OPRODVV<CT,T1,T2>(opvv.GetX1()/x,opvv.GetV2(),opvv.GetV3()); }

template <class T, class T1, class T2> inline OPRODVV<CT,T1,T2> operator/(
    const OPRODVV<CT,T1,T2>& opvv, CCT x)
{ return OPRODVV<CT,T1,T2>(opvv.GetX1()/CT(x),opvv.GetV2(),opvv.GetV3()); }

template <class T, class T1, class T2> inline OPRODVV<CT,T1,T2> operator/(
    const OPRODVV<CT,T1,T2>& opvv, VCT x)
{ return OPRODVV<CT,T1,T2>(opvv.GetX1()/CT(x),opvv.GetV2(),opvv.GetV3()); }

// (x*v)^v
template <class T, class T1> inline OPRODVV<T,T1,T> operator^(
    const PRODXV1<T,T1>& pxv, const GENVECTOR2<T>& v)
{ return OPRODVV<T,T1,T>(pxv.GetX1(),pxv.GetV2(),v); }

template <class T, class T1> inline OPRODVV<CT,T1,T> operator^(
    const PRODXV1<CT,T1>& pxv, const GENVECTOR2<T>& v)
{ return OPRODVV<CT,T1,T>(pxv.GetX1(),pxv.GetV2(),v); }

template <class T> inline OPRODVV<CT,T,CT> operator^(
    const PRODXV1<T,T>& pxv, const GENVECTOR2<CT>& v)
{ return OPRODVV<CT,T,CT>(pxv.GetX1(),pxv.GetV2(),v); }

// v^(x*v)
template <class T, class T2> inline OPRODVV<T,T2,T> operator^(
    const GENVECTOR1<T>& v, const PRODXV2<T,T2>& pxv)
{ return OPRODVV<T,T,T2>(pxv.GetX1(),v,pxv.GetV2()); }

template <class T> inline OPRODVV<CT,CT,T> operator^(
    const GENVECTOR1<CT>& v, const PRODXV2<T,T>& pxv)
{ return OPRODVV<CT,CT,T>(pxv.GetX1(),v,pxv.GetV2()); }

template <class T, class T2> inline OPRODVV<CT,T,T2> operator^(
    const GENVECTOR1<T>& v, const PRODXV2<CT,T2>& pxv)
{ return OPRODVV<CT,T,T2>(pxv.GetX1(),v,pxv.GetV2()); }

// (x*v)^(x*v)
template <class T, class T1, class T2> inline OPRODVV<T,T1,T2> operator^(
    const PRODXV1<T,T1>& pxv1, const PRODXV2<T,T2>& pxv2)
{
  return OPRODVV<T,T1,T2>(pxv1.GetX1()*pxv2.GetX1(),pxv1.GetV2(),pxv2.GetV2()); 
}

template <class T, class T1> inline OPRODVV<CT,T1,T> operator^(
    const PRODXV1<CT,T1>& pxv1, const PRODXV2<T,T>& pxv2)
{ 
  return OPRODVV<CT,T1,T>(pxv1.GetX1()*pxv2.GetX1(),pxv1.GetV2(),pxv2.GetV2()); 
}

template <class T, class T2> inline OPRODVV<CT,T,T2> operator^(
    const PRODXV1<T,T>& pxv1, const PRODXV2<CT,T2>& pxv2)
{ 
  return OPRODVV<CT,T,T2>(pxv1.GetX1()*pxv2.GetX1(),pxv1.GetV2(),pxv2.GetV2()); 
}

