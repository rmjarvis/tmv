// Things that need to be #defined on entry:
// (The values for a normal Matrix are given)
//
// PRODXM	ProdXM
// GENMATRIX2	GenMatrix

// -m
template <class T> inline PRODXM<T,T> operator-(const GENMATRIX2<T>& m)
{ return PRODXM<T,T>(T(-1),m); }

// x*m
template <class T> inline PRODXM<T,T> operator*(
    T x1, const GENMATRIX2<T>& m2) 
{ return PRODXM<T,T>(x1,m2); }

template <class T> inline PRODXM<CT,T> operator*(
    CT x1, const GENMATRIX2<T>& m2)
{ return PRODXM<CT,T>(x1,m2); }

template <class T> inline PRODXM<CT,T> operator*(
    CCT x1, const GENMATRIX2<T>& m2)
{ return PRODXM<CT,T>(CT(x1),m2); }

template <class T> inline PRODXM<CT,T> operator*(
    VCT x1, const GENMATRIX2<T>& m2)
{ return PRODXM<CT,T>(CT(x1),m2); }

template <class T> inline PRODXM<CT,CT> operator*(
    T x1, const GENMATRIX2<CT>& m2) 
{ return PRODXM<CT,CT>(x1,m2); }

template <class T> inline PRODXM<CT,CT> operator*(
    CCT x1, const GENMATRIX2<CT>& m2)
{ return PRODXM<CT,CT>(CT(x1),m2); }

template <class T> inline PRODXM<CT,CT> operator*(
    VCT x1, const GENMATRIX2<CT>& m2)
{ return PRODXM<CT,CT>(CT(x1),m2); }

// m*x
template <class T> inline PRODXM<T,T> operator*(
    const GENMATRIX2<T>& m1, T x2) 
{ return PRODXM<T,T>(x2,m1); }

template <class T> inline PRODXM<CT,T> operator*(
    const GENMATRIX2<T>& m1,CT x2)
{ return PRODXM<CT,T>(x2,m1); }

template <class T> inline PRODXM<CT,T> operator*(
    const GENMATRIX2<T>& m1, CCT x2)
{ return PRODXM<CT,T>(CT(x2),m1); }

template <class T> inline PRODXM<CT,T> operator*(
    const GENMATRIX2<T>& m1, VCT x2)
{ return PRODXM<CT,T>(CT(x2),m1); }

template <class T> inline PRODXM<CT,CT> operator*(
    const GENMATRIX2<CT>& m1, CCT x2)
{ return PRODXM<CT,CT>(CT(x2),m1); }

template <class T> inline PRODXM<CT,CT> operator*(
    const GENMATRIX2<CT>& m1, VCT x2)
{ return PRODXM<CT,CT>(CT(x2),m1); }

template <class T> inline PRODXM<CT,CT> operator*(
    const GENMATRIX2<CT>& m1, T x2) 
{ return PRODXM<CT,CT>(x2,m1); }

// m/x
template <class T> inline PRODXM<T,T> operator/(
    const GENMATRIX2<T>& m1, T x2) 
{ return PRODXM<T,T>(RealType(T)(1)/x2,m1); }

template <class T> inline PRODXM<CT,T> operator/(
    const GENMATRIX2<T>& m1, CT x2)
{ return PRODXM<CT,T>(T(1)/x2,m1); }

template <class T> inline PRODXM<CT,T> operator/(
    const GENMATRIX2<T>& m1, CCT x2)
{ return PRODXM<CT,T>(T(1)/CT(x2),m1); }

template <class T> inline PRODXM<CT,T> operator/(
    const GENMATRIX2<T>& m1, VCT x2)
{ return PRODXM<CT,T>(T(1)/CT(x2),m1); }

template <class T> inline PRODXM<CT,CT> operator/(
    const GENMATRIX2<CT>& m1, CCT x2)
{ return PRODXM<CT,CT>(T(1)/CT(x2),m1); }

template <class T> inline PRODXM<CT,CT> operator/(
    const GENMATRIX2<CT>& m1, VCT x2)
{ return PRODXM<CT,CT>(T(1)/CT(x2),m1); }

template <class T> inline PRODXM<CT,CT> operator/(
    const GENMATRIX2<CT>& m1, T x2)
{ return PRODXM<CT,CT>(T(1)/x2,m1); }

// -(x*m)
template <class T, class T2> inline PRODXM<T,T2> operator-(
    const PRODXM<T,T2>& pxm)
{ return PRODXM<T,T2>(-pxm.GetX1(),pxm.GetM2()); }

// x*(x*m)
template <class T, class T2> inline PRODXM<T,T2> operator*(
    const T x3, const PRODXM<T,T2>& pxm)
{ return PRODXM<T,T2>(x3*pxm.GetX1(),pxm.GetM2()); }

template <class T, class T2> inline PRODXM<CT,T2> operator*(
    const T x3, const PRODXM<CT,T2>& pxm)
{ return PRODXM<CT,T2>(x3*pxm.GetX1(),pxm.GetM2()); }

template <class T> inline PRODXM<CT,T> operator*(
    const CT x3, const PRODXM<T,T>& pxm)
{ return PRODXM<CT,T>(x3*pxm.GetX1(),pxm.GetM2()); }

template <class T> inline PRODXM<CT,T> operator*(
    const CCT x3, const PRODXM<T,T>& pxm)
{ return PRODXM<CT,T>(CT(x3)*pxm.GetX1(),pxm.GetM2()); }

template <class T> inline PRODXM<CT,T> operator*(
    const VCT x3, const PRODXM<T,T>& pxm)
{ return PRODXM<CT,T>(CT(x3)*pxm.GetX1(),pxm.GetM2()); }

template <class T, class T2> inline PRODXM<CT,T2> operator*(
    const CCT x3, const PRODXM<CT,T2>& pxm)
{ return PRODXM<CT,T2>(CT(x3)*pxm.GetX1(),pxm.GetM2()); }

template <class T, class T2> inline PRODXM<CT,T2> operator*(
    const VCT x3, const PRODXM<CT,T2>& pxm)
{ return PRODXM<CT,T2>(CT(x3)*pxm.GetX1(),pxm.GetM2()); }

// (x*m)*x
template <class T, class T2> inline PRODXM<T,T2> operator*(
    const PRODXM<T,T2>& pxm, const T x3)
{ return PRODXM<T,T2>(x3*pxm.GetX1(),pxm.GetM2()); }

template <class T, class T2> inline PRODXM<CT,T2> operator*(
    const PRODXM<CT,T2>& pxm, const T x3)
{ return PRODXM<CT,T2>(x3*pxm.GetX1(),pxm.GetM2()); }

template <class T> inline PRODXM<CT,T> operator*(
    const PRODXM<T,T>& pxm, const CT x3)
{ return PRODXM<CT,T>(x3*pxm.GetX1(),pxm.GetM2()); }

template <class T> inline PRODXM<CT,T> operator*(
    const PRODXM<T,T>& pxm, const CCT x3)
{ return PRODXM<CT,T>(CT(x3)*pxm.GetX1(),pxm.GetM2()); }

template <class T> inline PRODXM<CT,T> operator*(
    const PRODXM<T,T>& pxm, const VCT x3)
{ return PRODXM<CT,T>(CT(x3)*pxm.GetX1(),pxm.GetM2()); }

template <class T, class T2> inline PRODXM<CT,T2> operator*(
    const PRODXM<CT,T2>& pxm, const CCT x3)
{ return PRODXM<CT,T2>(CT(x3)*pxm.GetX1(),pxm.GetM2()); }

template <class T, class T2> inline PRODXM<CT,T2> operator*(
    const PRODXM<CT,T2>& pxm, const VCT x3)
{ return PRODXM<CT,T2>(CT(x3)*pxm.GetX1(),pxm.GetM2()); }

// (x*m)/x
template <class T, class T2> inline PRODXM<T,T2> operator/(
    const PRODXM<T,T2>& pxm, const T x3)
{ return PRODXM<T,T2>(pxm.GetX1()/x3,pxm.GetM2()); }

template <class T, class T2> inline PRODXM<CT,T2> operator/(
    const PRODXM<CT,T2>& pxm, const T x3)
{ return PRODXM<CT,T2>(pxm.GetX1()/x3,pxm.GetM2()); }

template <class T> inline PRODXM<CT,T> operator/(
    const PRODXM<T,T>& pxm, const CT x3)
{ return PRODXM<CT,T>(pxm.GetX1()/x3,pxm.GetM2()); }

template <class T> inline PRODXM<CT,T> operator/(
    const PRODXM<T,T>& pxm, const CCT x3)
{ return PRODXM<CT,T>(pxm.GetX1()/CT(x3),pxm.GetM2()); }

template <class T> inline PRODXM<CT,T> operator/(
    const PRODXM<T,T>& pxm, const VCT x3)
{ return PRODXM<CT,T>(pxm.GetX1()/CT(x3),pxm.GetM2()); }

template <class T, class T2> inline PRODXM<CT,T2> operator/(
    const PRODXM<CT,T2>& pxm, const CCT x3)
{ return PRODXM<CT,T2>(pxm.GetX1()/CT(x3),pxm.GetM2()); }

template <class T, class T2> inline PRODXM<CT,T2> operator/(
    const PRODXM<CT,T2>& pxm, const VCT x3)
{ return PRODXM<CT,T2>(pxm.GetX1()/CT(x3),pxm.GetM2()); }

