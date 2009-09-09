// Things that need to be #defined on entry:
// (The values for a normal Matrix are given)
//
// PRODXM	ProdXM
// GENMATRIX	GenMatrix

// -m
template <class T> inline PRODXM<T,T> operator-(const GENMATRIX<T>& m)
{ return PRODXM<T,T>(T(-1),m); }

// x*m
template <class T> inline PRODXM<T,T> operator*(
    T x, const GENMATRIX<T>& m) 
{ return PRODXM<T,T>(x,m); }

template <class T> inline PRODXM<CT,T> operator*(
    CT x, const GENMATRIX<T>& m)
{ return PRODXM<CT,T>(x,m); }

template <class T> inline PRODXM<CT,T> operator*(
    CCT x, const GENMATRIX<T>& m)
{ return PRODXM<CT,T>(CT(x),m); }

template <class T> inline PRODXM<CT,T> operator*(
    VCT x, const GENMATRIX<T>& m)
{ return PRODXM<CT,T>(CT(x),m); }

template <class T> inline PRODXM<CT,CT> operator*(
    T x, const GENMATRIX<CT>& m) 
{ return PRODXM<CT,CT>(x,m); }

template <class T> inline PRODXM<CT,CT> operator*(
    CCT x, const GENMATRIX<CT>& m)
{ return PRODXM<CT,CT>(CT(x),m); }

template <class T> inline PRODXM<CT,CT> operator*(
    VCT x, const GENMATRIX<CT>& m)
{ return PRODXM<CT,CT>(CT(x),m); }

// m*x
template <class T> inline PRODXM<T,T> operator*(
    const GENMATRIX<T>& m, T x) 
{ return PRODXM<T,T>(x,m); }

template <class T> inline PRODXM<CT,T> operator*(
    const GENMATRIX<T>& m,CT x)
{ return PRODXM<CT,T>(x,m); }

template <class T> inline PRODXM<CT,T> operator*(
    const GENMATRIX<T>& m, CCT x)
{ return PRODXM<CT,T>(CT(x),m); }

template <class T> inline PRODXM<CT,T> operator*(
    const GENMATRIX<T>& m, VCT x)
{ return PRODXM<CT,T>(CT(x),m); }

template <class T> inline PRODXM<CT,CT> operator*(
    const GENMATRIX<CT>& m, CCT x)
{ return PRODXM<CT,CT>(CT(x),m); }

template <class T> inline PRODXM<CT,CT> operator*(
    const GENMATRIX<CT>& m, VCT x)
{ return PRODXM<CT,CT>(CT(x),m); }

template <class T> inline PRODXM<CT,CT> operator*(
    const GENMATRIX<CT>& m, T x) 
{ return PRODXM<CT,CT>(x,m); }

// m/x
template <class T> inline PRODXM<T,T> operator/(
    const GENMATRIX<T>& m, T x) 
{ return PRODXM<T,T>(RealType(T)(1)/x,m); }

template <class T> inline PRODXM<CT,T> operator/(
    const GENMATRIX<T>& m, CT x)
{ return PRODXM<CT,T>(T(1)/x,m); }

template <class T> inline PRODXM<CT,T> operator/(
    const GENMATRIX<T>& m, CCT x)
{ return PRODXM<CT,T>(T(1)/CT(x),m); }

template <class T> inline PRODXM<CT,T> operator/(
    const GENMATRIX<T>& m, VCT x)
{ return PRODXM<CT,T>(T(1)/CT(x),m); }

template <class T> inline PRODXM<CT,CT> operator/(
    const GENMATRIX<CT>& m, CCT x)
{ return PRODXM<CT,CT>(T(1)/CT(x),m); }

template <class T> inline PRODXM<CT,CT> operator/(
    const GENMATRIX<CT>& m, VCT x)
{ return PRODXM<CT,CT>(T(1)/CT(x),m); }

template <class T> inline PRODXM<CT,CT> operator/(
    const GENMATRIX<CT>& m, T x)
{ return PRODXM<CT,CT>(T(1)/x,m); }

// -(x*m)
template <class T, class T2> inline PRODXM<T,T2> operator-(
    const PRODXM<T,T2>& pxm)
{ return PRODXM<T,T2>(-pxm.GetX(),pxm.GetM()); }

// x*(x*m)
template <class T, class T2> inline PRODXM<T,T2> operator*(
    const T x, const PRODXM<T,T2>& pxm)
{ return PRODXM<T,T2>(x*pxm.GetX(),pxm.GetM()); }

template <class T, class T2> inline PRODXM<CT,T2> operator*(
    const T x, const PRODXM<CT,T2>& pxm)
{ return PRODXM<CT,T2>(x*pxm.GetX(),pxm.GetM()); }

template <class T> inline PRODXM<CT,T> operator*(
    const CT x, const PRODXM<T,T>& pxm)
{ return PRODXM<CT,T>(x*pxm.GetX(),pxm.GetM()); }

template <class T> inline PRODXM<CT,T> operator*(
    const CCT x, const PRODXM<T,T>& pxm)
{ return PRODXM<CT,T>(CT(x)*pxm.GetX(),pxm.GetM()); }

template <class T> inline PRODXM<CT,T> operator*(
    const VCT x, const PRODXM<T,T>& pxm)
{ return PRODXM<CT,T>(CT(x)*pxm.GetX(),pxm.GetM()); }

template <class T, class T2> inline PRODXM<CT,T2> operator*(
    const CCT x, const PRODXM<CT,T2>& pxm)
{ return PRODXM<CT,T2>(CT(x)*pxm.GetX(),pxm.GetM()); }

template <class T, class T2> inline PRODXM<CT,T2> operator*(
    const VCT x, const PRODXM<CT,T2>& pxm)
{ return PRODXM<CT,T2>(CT(x)*pxm.GetX(),pxm.GetM()); }

// (x*m)*x
template <class T, class T2> inline PRODXM<T,T2> operator*(
    const PRODXM<T,T2>& pxm, const T x)
{ return PRODXM<T,T2>(x*pxm.GetX(),pxm.GetM()); }

template <class T, class T2> inline PRODXM<CT,T2> operator*(
    const PRODXM<CT,T2>& pxm, const T x)
{ return PRODXM<CT,T2>(x*pxm.GetX(),pxm.GetM()); }

template <class T> inline PRODXM<CT,T> operator*(
    const PRODXM<T,T>& pxm, const CT x)
{ return PRODXM<CT,T>(x*pxm.GetX(),pxm.GetM()); }

template <class T> inline PRODXM<CT,T> operator*(
    const PRODXM<T,T>& pxm, const CCT x)
{ return PRODXM<CT,T>(CT(x)*pxm.GetX(),pxm.GetM()); }

template <class T> inline PRODXM<CT,T> operator*(
    const PRODXM<T,T>& pxm, const VCT x)
{ return PRODXM<CT,T>(CT(x)*pxm.GetX(),pxm.GetM()); }

template <class T, class T2> inline PRODXM<CT,T2> operator*(
    const PRODXM<CT,T2>& pxm, const CCT x)
{ return PRODXM<CT,T2>(CT(x)*pxm.GetX(),pxm.GetM()); }

template <class T, class T2> inline PRODXM<CT,T2> operator*(
    const PRODXM<CT,T2>& pxm, const VCT x)
{ return PRODXM<CT,T2>(CT(x)*pxm.GetX(),pxm.GetM()); }

// (x*m)/x
template <class T, class T2> inline PRODXM<T,T2> operator/(
    const PRODXM<T,T2>& pxm, const T x)
{ return PRODXM<T,T2>(pxm.GetX()/x,pxm.GetM()); }

template <class T, class T2> inline PRODXM<CT,T2> operator/(
    const PRODXM<CT,T2>& pxm, const T x)
{ return PRODXM<CT,T2>(pxm.GetX()/x,pxm.GetM()); }

template <class T> inline PRODXM<CT,T> operator/(
    const PRODXM<T,T>& pxm, const CT x)
{ return PRODXM<CT,T>(pxm.GetX()/x,pxm.GetM()); }

template <class T> inline PRODXM<CT,T> operator/(
    const PRODXM<T,T>& pxm, const CCT x)
{ return PRODXM<CT,T>(pxm.GetX()/CT(x),pxm.GetM()); }

template <class T> inline PRODXM<CT,T> operator/(
    const PRODXM<T,T>& pxm, const VCT x)
{ return PRODXM<CT,T>(pxm.GetX()/CT(x),pxm.GetM()); }

template <class T, class T2> inline PRODXM<CT,T2> operator/(
    const PRODXM<CT,T2>& pxm, const CCT x)
{ return PRODXM<CT,T2>(pxm.GetX()/CT(x),pxm.GetM()); }

template <class T, class T2> inline PRODXM<CT,T2> operator/(
    const PRODXM<CT,T2>& pxm, const VCT x)
{ return PRODXM<CT,T2>(pxm.GetX()/CT(x),pxm.GetM()); }

