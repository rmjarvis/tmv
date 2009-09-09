// Things that need to be #defined on entry:
// (The values for a normal Vector are given)
//
// PRODXV	ProdXV
// GENVECTOR	GenVector

// -v
template <class T> inline PRODXV<T,T> operator-(
    const GENVECTOR<T>& v)
{ return PRODXV<T,T>(T(-1),v); }

// x*v
template <class T> inline PRODXV<T,T> operator*(
    T x1, const GENVECTOR<T>& v2) 
{ return PRODXV<T,T>(x1,v2); }

template <class T> inline PRODXV<CT,T> operator*(
    CT x1, const GENVECTOR<T>& v2)
{ return PRODXV<CT,T>(x1,v2); }

template <class T> inline PRODXV<CT,T> operator*(
    CCT x1, const GENVECTOR<T>& v2)
{ return PRODXV<CT,T>(CT(x1),v2); }

template <class T> inline PRODXV<CT,T> operator*(
    VCT x1, const GENVECTOR<T>& v2)
{ return PRODXV<CT,T>(CT(x1),v2); }

template <class T> inline PRODXV<CT,CT> operator*(
    T x1, const GENVECTOR<CT>& v2) 
{ return PRODXV<CT,CT>(x1,v2); }

template <class T> inline PRODXV<CT,CT> operator*(
    CCT x1, const GENVECTOR<CT>& v2)
{ return PRODXV<CT,CT>(CT(x1),v2); }

template <class T> inline PRODXV<CT,CT> operator*(
    VCT x1, const GENVECTOR<CT>& v2)
{ return PRODXV<CT,CT>(CT(x1),v2); }

// v*x
template <class T> inline PRODXV<T,T> operator*(
    const GENVECTOR<T>& v1, T x2) 
{ return PRODXV<T,T>(x2,v1); }

template <class T> inline PRODXV<CT,T> operator*(
    const GENVECTOR<T>& v1, CT x2)
{ return PRODXV<CT,T>(x2,v1); }

template <class T> inline PRODXV<CT,T> operator*(
    const GENVECTOR<T>& v1, CCT x2)
{ return PRODXV<CT,T>(CT(x2),v1); }

template <class T> inline PRODXV<CT,T> operator*(
    const GENVECTOR<T>& v1, VCT x2)
{ return PRODXV<CT,T>(CT(x2),v1); }

template <class T> inline PRODXV<CT,CT> operator*(
    const GENVECTOR<CT>& v1, T x2) 
{ return PRODXV<CT,CT>(x2,v1); }

template <class T> inline PRODXV<CT,CT> operator*(
    const GENVECTOR<CT>& v1, CCT x2)
{ return PRODXV<CT,CT>(CT(x2),v1); }

template <class T> inline PRODXV<CT,CT> operator*(
    const GENVECTOR<CT>& v1, VCT x2)
{ return PRODXV<CT,CT>(CT(x2),v1); }

// v/x
template <class T> inline PRODXV<T,T> operator/(
    const GENVECTOR<T>& v1, T x2) 
{ return PRODXV<T,T>(RealType(T)(1)/x2,v1); }

template <class T> inline PRODXV<CT,T> operator/(
    const GENVECTOR<T>& v1, CT x2)
{ return PRODXV<CT,T>(T(1)/x2,v1); }

template <class T> inline PRODXV<CT,T> operator/(
    const GENVECTOR<T>& v1, CCT x2)
{ return PRODXV<CT,T>(T(1)/CT(x2),v1); }

template <class T> inline PRODXV<CT,T> operator/(
    const GENVECTOR<T>& v1, VCT x2)
{ return PRODXV<CT,T>(T(1)/CT(x2),v1); }

template <class T> inline PRODXV<CT,CT> operator/(
    const GENVECTOR<CT>& v1, T x2)
{ return PRODXV<CT,CT>(T(1)/x2,v1); }

template <class T> inline PRODXV<CT,CT> operator/(
    const GENVECTOR<CT>& v1, CCT x2)
{ return PRODXV<CT,CT>(T(1)/CT(x2),v1); }

template <class T> inline PRODXV<CT,CT> operator/(
    const GENVECTOR<CT>& v1, VCT x2)
{ return PRODXV<CT,CT>(T(1)/CT(x2),v1); }

// -(x*v)
template <class T, class T2> inline PRODXV<T,T2> operator-(
    const PRODXV<T,T2>& pxv)
{ return PRODXV<T,T2>(-pxv.GetX1(),pxv.GetV2()); }

// x*(x*v)
template <class T, class T2> inline PRODXV<T,T2> operator*(
    const T x3, const PRODXV<T,T2>& pxv)
{ return PRODXV<T,T2>(x3*pxv.GetX1(),pxv.GetV2()); }

template <class T, class T2> inline PRODXV<CT,T2> operator*(
    const T x3, const PRODXV<CT,T2>& pxv)
{ return PRODXV<CT,T2>(x3*pxv.GetX1(),pxv.GetV2()); }

template <class T> inline PRODXV<CT,T> operator*(
    const CT x3, const PRODXV<T,T>& pxv)
{ return PRODXV<CT,T>(x3*pxv.GetX1(),pxv.GetV2()); }

template <class T> inline PRODXV<CT,T> operator*(
    const CCT x3, const PRODXV<T,T>& pxv)
{ return PRODXV<CT,T>(CT(x3)*pxv.GetX1(),pxv.GetV2()); }

template <class T> inline PRODXV<CT,T> operator*(
    const VCT x3, const PRODXV<T,T>& pxv)
{ return PRODXV<CT,T>(CT(x3)*pxv.GetX1(),pxv.GetV2()); }

template <class T, class T2> inline PRODXV<CT,T2> operator*(
    const CCT x3, const PRODXV<CT,T2>& pxv)
{ return PRODXV<CT,T2>(CT(x3)*pxv.GetX1(),pxv.GetV2()); }

template <class T, class T2> inline PRODXV<CT,T2> operator*(
    const VCT x3, const PRODXV<CT,T2>& pxv)
{ return PRODXV<CT,T2>(CT(x3)*pxv.GetX1(),pxv.GetV2()); }

// (x*v)*x
template <class T, class T2> inline PRODXV<T,T2> operator*(
    const PRODXV<T,T2>& pxv, const T x3)
{ return PRODXV<T,T2>(x3*pxv.GetX1(),pxv.GetV2()); }

template <class T, class T2> inline PRODXV<CT,T2> operator*(
    const PRODXV<CT,T2>& pxv, const T x3)
{ return PRODXV<CT,T2>(x3*pxv.GetX1(),pxv.GetV2()); }

template <class T> inline PRODXV<CT,T> operator*(
    const PRODXV<T,T>& pxv, const CT x3)
{ return PRODXV<CT,T>(x3*pxv.GetX1(),pxv.GetV2()); }

template <class T> inline PRODXV<CT,T> operator*(
    const PRODXV<T,T>& pxv, const CCT x3)
{ return PRODXV<CT,T>(CT(x3)*pxv.GetX1(),pxv.GetV2()); }

template <class T> inline PRODXV<CT,T> operator*(
    const PRODXV<T,T>& pxv, const VCT x3)
{ return PRODXV<CT,T>(CT(x3)*pxv.GetX1(),pxv.GetV2()); }

template <class T, class T2> inline PRODXV<CT,T2> operator*(
    const PRODXV<CT,T2>& pxv, const CCT x3)
{ return PRODXV<CT,T2>(CT(x3)*pxv.GetX1(),pxv.GetV2()); }

template <class T, class T2> inline PRODXV<CT,T2> operator*(
    const PRODXV<CT,T2>& pxv, const VCT x3)
{ return PRODXV<CT,T2>(CT(x3)*pxv.GetX1(),pxv.GetV2()); }

// (x*v)/x
template <class T, class T2> inline PRODXV<T,T> operator/(
    const PRODXV<T,T2>& pxv, const T x3)
{ return PRODXV<T,T2>(pxv.GetX1()/x3,pxv.GetV2()); }

template <class T, class T2> inline PRODXV<CT,T2> operator/(
    const PRODXV<CT,T2>& pxv, const T x3)
{ return PRODXV<CT,T2>(pxv.GetX1()/x3,pxv.GetV2()); }

template <class T> inline PRODXV<CT,T> operator/(
    const PRODXV<T,T>& pxv, const CT x3)
{ return PRODXV<CT,T>(pxv.GetX1()/x3,pxv.GetV2()); }

template <class T> inline PRODXV<CT,T> operator/(
    const PRODXV<T,T>& pxv, const CCT x3)
{ return PRODXV<CT,T>(pxv.GetX1()/CT(x3),pxv.GetV2()); }

template <class T> inline PRODXV<CT,T> operator/(
    const PRODXV<T,T>& pxv, const VCT x3)
{ return PRODXV<CT,T>(pxv.GetX1()/CT(x3),pxv.GetV2()); }

template <class T, class T2> inline PRODXV<CT,T2> operator/(
    const PRODXV<CT,T2>& pxv, const CCT x3)
{ return PRODXV<CT,T2>(pxv.GetX1()/CT(x3),pxv.GetV2()); }

template <class T, class T2> inline PRODXV<CT,T2> operator/(
    const PRODXV<CT,T2>& pxv, const VCT x3)
{ return PRODXV<CT,T2>(pxv.GetX1()/CT(x3),pxv.GetV2()); }

