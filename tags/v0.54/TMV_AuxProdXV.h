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
    T x, const GENVECTOR<T>& v) 
{ return PRODXV<T,T>(x,v); }

template <class T> inline PRODXV<CT,T> operator*(
    CT x, const GENVECTOR<T>& v)
{ return PRODXV<CT,T>(x,v); }

template <class T> inline PRODXV<CT,T> operator*(
    CCT x, const GENVECTOR<T>& v)
{ return PRODXV<CT,T>(CT(x),v); }

template <class T> inline PRODXV<CT,T> operator*(
    VCT x, const GENVECTOR<T>& v)
{ return PRODXV<CT,T>(CT(x),v); }

template <class T> inline PRODXV<CT,CT> operator*(
    T x, const GENVECTOR<CT>& v) 
{ return PRODXV<CT,CT>(x,v); }

template <class T> inline PRODXV<CT,CT> operator*(
    CCT x, const GENVECTOR<CT>& v)
{ return PRODXV<CT,CT>(CT(x),v); }

template <class T> inline PRODXV<CT,CT> operator*(
    VCT x, const GENVECTOR<CT>& v)
{ return PRODXV<CT,CT>(CT(x),v); }

// v*x
template <class T> inline PRODXV<T,T> operator*(
    const GENVECTOR<T>& v, T x) 
{ return PRODXV<T,T>(x,v); }

template <class T> inline PRODXV<CT,T> operator*(
    const GENVECTOR<T>& v, CT x)
{ return PRODXV<CT,T>(x,v); }

template <class T> inline PRODXV<CT,T> operator*(
    const GENVECTOR<T>& v, CCT x)
{ return PRODXV<CT,T>(CT(x),v); }

template <class T> inline PRODXV<CT,T> operator*(
    const GENVECTOR<T>& v, VCT x)
{ return PRODXV<CT,T>(CT(x),v); }

template <class T> inline PRODXV<CT,CT> operator*(
    const GENVECTOR<CT>& v, T x) 
{ return PRODXV<CT,CT>(x,v); }

template <class T> inline PRODXV<CT,CT> operator*(
    const GENVECTOR<CT>& v, CCT x)
{ return PRODXV<CT,CT>(CT(x),v); }

template <class T> inline PRODXV<CT,CT> operator*(
    const GENVECTOR<CT>& v, VCT x)
{ return PRODXV<CT,CT>(CT(x),v); }

// v/x
template <class T> inline PRODXV<T,T> operator/(
    const GENVECTOR<T>& v, T x) 
{ return PRODXV<T,T>(RealType(T)(1)/x,v); }

template <class T> inline PRODXV<CT,T> operator/(
    const GENVECTOR<T>& v, CT x)
{ return PRODXV<CT,T>(T(1)/x,v); }

template <class T> inline PRODXV<CT,T> operator/(
    const GENVECTOR<T>& v, CCT x)
{ return PRODXV<CT,T>(T(1)/CT(x),v); }

template <class T> inline PRODXV<CT,T> operator/(
    const GENVECTOR<T>& v, VCT x)
{ return PRODXV<CT,T>(T(1)/CT(x),v); }

template <class T> inline PRODXV<CT,CT> operator/(
    const GENVECTOR<CT>& v, T x)
{ return PRODXV<CT,CT>(T(1)/x,v); }

template <class T> inline PRODXV<CT,CT> operator/(
    const GENVECTOR<CT>& v, CCT x)
{ return PRODXV<CT,CT>(T(1)/CT(x),v); }

template <class T> inline PRODXV<CT,CT> operator/(
    const GENVECTOR<CT>& v, VCT x)
{ return PRODXV<CT,CT>(T(1)/CT(x),v); }

// -(x*v)
template <class T, class T2> inline PRODXV<T,T2> operator-(
    const PRODXV<T,T2>& pxv)
{ return PRODXV<T,T2>(-pxv.GetX(),pxv.GetV()); }

// x*(x*v)
template <class T, class T2> inline PRODXV<T,T2> operator*(
    const T x, const PRODXV<T,T2>& pxv)
{ return PRODXV<T,T2>(x*pxv.GetX(),pxv.GetV()); }

template <class T, class T2> inline PRODXV<CT,T2> operator*(
    const T x, const PRODXV<CT,T2>& pxv)
{ return PRODXV<CT,T2>(x*pxv.GetX(),pxv.GetV()); }

template <class T> inline PRODXV<CT,T> operator*(
    const CT x, const PRODXV<T,T>& pxv)
{ return PRODXV<CT,T>(x*pxv.GetX(),pxv.GetV()); }

template <class T> inline PRODXV<CT,T> operator*(
    const CCT x, const PRODXV<T,T>& pxv)
{ return PRODXV<CT,T>(CT(x)*pxv.GetX(),pxv.GetV()); }

template <class T> inline PRODXV<CT,T> operator*(
    const VCT x, const PRODXV<T,T>& pxv)
{ return PRODXV<CT,T>(CT(x)*pxv.GetX(),pxv.GetV()); }

template <class T, class T2> inline PRODXV<CT,T2> operator*(
    const CCT x, const PRODXV<CT,T2>& pxv)
{ return PRODXV<CT,T2>(CT(x)*pxv.GetX(),pxv.GetV()); }

template <class T, class T2> inline PRODXV<CT,T2> operator*(
    const VCT x, const PRODXV<CT,T2>& pxv)
{ return PRODXV<CT,T2>(CT(x)*pxv.GetX(),pxv.GetV()); }

// (x*v)*x
template <class T, class T2> inline PRODXV<T,T2> operator*(
    const PRODXV<T,T2>& pxv, const T x)
{ return PRODXV<T,T2>(x*pxv.GetX(),pxv.GetV()); }

template <class T, class T2> inline PRODXV<CT,T2> operator*(
    const PRODXV<CT,T2>& pxv, const T x)
{ return PRODXV<CT,T2>(x*pxv.GetX(),pxv.GetV()); }

template <class T> inline PRODXV<CT,T> operator*(
    const PRODXV<T,T>& pxv, const CT x)
{ return PRODXV<CT,T>(x*pxv.GetX(),pxv.GetV()); }

template <class T> inline PRODXV<CT,T> operator*(
    const PRODXV<T,T>& pxv, const CCT x)
{ return PRODXV<CT,T>(CT(x)*pxv.GetX(),pxv.GetV()); }

template <class T> inline PRODXV<CT,T> operator*(
    const PRODXV<T,T>& pxv, const VCT x)
{ return PRODXV<CT,T>(CT(x)*pxv.GetX(),pxv.GetV()); }

template <class T, class T2> inline PRODXV<CT,T2> operator*(
    const PRODXV<CT,T2>& pxv, const CCT x)
{ return PRODXV<CT,T2>(CT(x)*pxv.GetX(),pxv.GetV()); }

template <class T, class T2> inline PRODXV<CT,T2> operator*(
    const PRODXV<CT,T2>& pxv, const VCT x)
{ return PRODXV<CT,T2>(CT(x)*pxv.GetX(),pxv.GetV()); }

// (x*v)/x
template <class T, class T2> inline PRODXV<T,T> operator/(
    const PRODXV<T,T2>& pxv, const T x)
{ return PRODXV<T,T2>(pxv.GetX()/x,pxv.GetV()); }

template <class T, class T2> inline PRODXV<CT,T2> operator/(
    const PRODXV<CT,T2>& pxv, const T x)
{ return PRODXV<CT,T2>(pxv.GetX()/x,pxv.GetV()); }

template <class T> inline PRODXV<CT,T> operator/(
    const PRODXV<T,T>& pxv, const CT x)
{ return PRODXV<CT,T>(pxv.GetX()/x,pxv.GetV()); }

template <class T> inline PRODXV<CT,T> operator/(
    const PRODXV<T,T>& pxv, const CCT x)
{ return PRODXV<CT,T>(pxv.GetX()/CT(x),pxv.GetV()); }

template <class T> inline PRODXV<CT,T> operator/(
    const PRODXV<T,T>& pxv, const VCT x)
{ return PRODXV<CT,T>(pxv.GetX()/CT(x),pxv.GetV()); }

template <class T, class T2> inline PRODXV<CT,T2> operator/(
    const PRODXV<CT,T2>& pxv, const CCT x)
{ return PRODXV<CT,T2>(pxv.GetX()/CT(x),pxv.GetV()); }

template <class T, class T2> inline PRODXV<CT,T2> operator/(
    const PRODXV<CT,T2>& pxv, const VCT x)
{ return PRODXV<CT,T2>(pxv.GetX()/CT(x),pxv.GetV()); }

