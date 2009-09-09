// Things that need to be #defined on entry:
// (The values for a normal Matrix are given)
//
// SUMMX	SumMX
// PRODXM2	ProdXM
// GENMATRIX2	GenMatrix

// m+x
template <class T> inline SUMMX<T,T> operator+(
    const GENMATRIX2<T>& m1, T x2)
{ return SUMMX<T,T>(T(1),m1,x2); }

template <class T> inline SUMMX<CT,T> operator+(
    const GENMATRIX2<T>& m1, CT x2)
{ return SUMMX<CT,T>(CT(1),m1,x2); }

template <class T> inline SUMMX<CT,T> operator+(
    const GENMATRIX2<T>& m1, CCT x2)
{ return SUMMX<CT,T>(CT(1),m1,CT(x2)); }

template <class T> inline SUMMX<CT,T> operator+(
    const GENMATRIX2<T>& m1, VCT x2)
{ return SUMMX<CT,T>(CT(1),m1,CT(x2)); }

template <class T> inline SUMMX<CT,CT> operator+(
    const GENMATRIX2<CT>& m1, T x2)
{ return SUMMX<CT,CT>(CT(1),m1,x2); }

template <class T> inline SUMMX<CT,CT> operator+(
    const GENMATRIX2<CT>& m1, CCT x2)
{ return SUMMX<CT,CT>(CT(1),m1,CT(x2)); }

template <class T> inline SUMMX<CT,CT> operator+(
    const GENMATRIX2<CT>& m1, VCT x2)
{ return SUMMX<CT,CT>(CT(1),m1,CT(x2)); }

// x+m
template <class T> inline SUMMX<T,T> operator+(
    T x1, const GENMATRIX2<T>& m2)
{ return SUMMX<T,T>(T(1),m2,x1); }

template <class T> inline SUMMX<CT,T> operator+(
    CT x1, const GENMATRIX2<T>& m2)
{ return SUMMX<CT,T>(CT(1),m2,x1); }

template <class T> inline SUMMX<CT,T> operator+(
    CCT x1, const GENMATRIX2<T>& m2)
{ return SUMMX<CT,T>(CT(1),m2,CT(x1)); }

template <class T> inline SUMMX<CT,T> operator+(
    VCT x1, const GENMATRIX2<T>& m2)
{ return SUMMX<CT,T>(CT(1),m2,CT(x1)); }

template <class T> inline SUMMX<CT,CT> operator+(
    T x1, const GENMATRIX2<CT>& m2)
{ return SUMMX<CT,CT>(CT(1),m2,x1); }

template <class T> inline SUMMX<CT,CT> operator+(
    CCT x1, const GENMATRIX2<CT>& m2)
{ return SUMMX<CT,CT>(CT(1),m2,CT(x1)); }

template <class T> inline SUMMX<CT,CT> operator+(
    VCT x1, const GENMATRIX2<CT>& m2)
{ return SUMMX<CT,CT>(CT(1),m2,CT(x1)); }

// m-x
template <class T> inline SUMMX<T,T> operator-(
    const GENMATRIX2<T>& m1, T x2)
{ return SUMMX<T,T>(T(1),m1,-x2); }

template <class T> inline SUMMX<CT,T> operator-(
    const GENMATRIX2<T>& m1, CT x2)
{ return SUMMX<CT,T>(CT(1),m1,-x2); }

template <class T> inline SUMMX<CT,T> operator-(
    const GENMATRIX2<T>& m1, CCT x2)
{ return SUMMX<CT,T>(CT(1),m1,-CT(x2)); }

template <class T> inline SUMMX<CT,T> operator-(
    const GENMATRIX2<T>& m1, VCT x2)
{ return SUMMX<CT,T>(CT(1),m1,-CT(x2)); }

template <class T> inline SUMMX<CT,CT> operator-(
    const GENMATRIX2<CT>& m1, T x2)
{ return SUMMX<CT,CT>(CT(1),m1,-x2); }

template <class T> inline SUMMX<CT,CT> operator-(
    const GENMATRIX2<CT>& m1, CCT x2)
{ return SUMMX<CT,CT>(CT(1),m1,-CT(x2)); }

template <class T> inline SUMMX<CT,CT> operator-(
    const GENMATRIX2<CT>& m1, VCT x2)
{ return SUMMX<CT,CT>(CT(1),m1,-CT(x2)); }

// x-m
template <class T> inline SUMMX<T,T> operator-(
    T x1, const GENMATRIX2<T>& m2)
{ return SUMMX<T,T>(T(-1),m2,x1); }

template <class T> inline SUMMX<CT,T> operator-(
    CT x1, const GENMATRIX2<T>& m2)
{ return SUMMX<CT,T>(CT(-1),m2,x1); }

template <class T> inline SUMMX<CT,T> operator-(
    CCT x1, const GENMATRIX2<T>& m2)
{ return SUMMX<CT,T>(CT(-1),m2,CT(x1)); }

template <class T> inline SUMMX<CT,T> operator-(
    VCT x1, const GENMATRIX2<T>& m2)
{ return SUMMX<CT,T>(CT(-1),m2,CT(x1)); }

template <class T> inline SUMMX<CT,CT> operator-(
    T x1, const GENMATRIX2<CT>& m2)
{ return SUMMX<CT,CT>(CT(-1),m2,x1); }

template <class T> inline SUMMX<CT,CT> operator-(
    CCT x1, const GENMATRIX2<CT>& m2)
{ return SUMMX<CT,CT>(CT(-1),m2,CT(x1)); }

template <class T> inline SUMMX<CT,CT> operator-(
    VCT x1, const GENMATRIX2<CT>& m2)
{ return SUMMX<CT,CT>(CT(-1),m2,CT(x1)); }

// x+(x*m)
template <class T, class T2> inline SUMMX<T,T2> operator+(
    const T x3, const PRODXM<T,T2>& pxm)
{ return SUMMX<T,T2>(pxm.GetX1(),pxm.GetM2(),x3); }

template <class T, class T2> inline SUMMX<CT,T2> operator+(
    const T x3, const PRODXM<CT,T2>& pxm)
{ return SUMMX<CT,T2>(pxm.GetX1(),pxm.GetM2(),x3); }

template <class T> inline SUMMX<CT,T> operator+(
    const CT x3, const PRODXM<T,T>& pxm)
{ return SUMMX<CT,T>(pxm.GetX1(),pxm.GetM2(),x3); }

template <class T> inline SUMMX<CT,T> operator+(
    const CCT x3, const PRODXM<T,T>& pxm)
{ return SUMMX<CT,T>(pxm.GetX1(),pxm.GetM2(),CT(x3)); }

template <class T> inline SUMMX<CT,T> operator+(
    const VCT x3, const PRODXM<T,T>& pxm)
{ return SUMMX<CT,T>(pxm.GetX1(),pxm.GetM2(),CT(x3)); }

template <class T, class T2> inline SUMMX<CT,T2> operator+(
    const CCT x3, const PRODXM<CT,T2>& pxm)
{ return SUMMX<CT,T2>(pxm.GetX1(),pxm.GetM2(),CT(x3)); }

template <class T, class T2> inline SUMMX<CT,T2> operator+(
    const VCT x3, const PRODXM<CT,T2>& pxm)
{ return SUMMX<CT,T2>(pxm.GetX1(),pxm.GetM2(),CT(x3)); }

// (x*m)+x
template <class T, class T2> inline SUMMX<T,T2> operator+(
    const PRODXM<T,T2>& pxm, const T x3)
{ return SUMMX<T,T2>(pxm.GetX1(),pxm.GetM2(),x3); }

template <class T, class T2> inline SUMMX<CT,T2> operator+(
    const PRODXM<CT,T2>& pxm, const T x3)
{ return SUMMX<CT,T2>(pxm.GetX1(),pxm.GetM2(),x3); }

template <class T> inline SUMMX<CT,T> operator+(
    const PRODXM<T,T>& pxm, const CT x3)
{ return SUMMX<CT,T>(pxm.GetX1(),pxm.GetM2(),x3); }

template <class T> inline SUMMX<CT,T> operator+(
    const PRODXM<T,T>& pxm, const CCT x3)
{ return SUMMX<CT,T>(pxm.GetX1(),pxm.GetM2(),CT(x3)); }

template <class T> inline SUMMX<CT,T> operator+(
    const PRODXM<T,T>& pxm, const VCT x3)
{ return SUMMX<CT,T>(pxm.GetX1(),pxm.GetM2(),CT(x3)); }

template <class T, class T2> inline SUMMX<CT,T2> operator+(
    const PRODXM<CT,T2>& pxm, const CCT x3)
{ return SUMMX<CT,T2>(pxm.GetX1(),pxm.GetM2(),CT(x3)); }

template <class T, class T2> inline SUMMX<CT,T2> operator+(
    const PRODXM<CT,T2>& pxm, const VCT x3)
{ return SUMMX<CT,T2>(pxm.GetX1(),pxm.GetM2(),CT(x3)); }

// x-(x*m)
template <class T, class T2> inline SUMMX<T,T2> operator-(
    const T x3, const PRODXM<T,T2>& pxm)
{ return SUMMX<T,T2>(-pxm.GetX1(),pxm.GetM2(),x3); }

template <class T, class T2> inline SUMMX<CT,T2> operator-(
    const T x3, const PRODXM<CT,T2>& pxm)
{ return SUMMX<CT,T2>(-pxm.GetX1(),pxm.GetM2(),x3); }

template <class T> inline SUMMX<CT,T> operator-(
    const CT x3, const PRODXM<T,T>& pxm)
{ return SUMMX<CT,T>(-pxm.GetX1(),pxm.GetM2(),x3); }

template <class T> inline SUMMX<CT,T> operator-(
    const CCT x3, const PRODXM<T,T>& pxm)
{ return SUMMX<CT,T>(-pxm.GetX1(),pxm.GetM2(),CT(x3)); }

template <class T> inline SUMMX<CT,T> operator-(
    const VCT x3, const PRODXM<T,T>& pxm)
{ return SUMMX<CT,T>(-pxm.GetX1(),pxm.GetM2(),CT(x3)); }

template <class T, class T2> inline SUMMX<CT,T2> operator-(
    const CCT x3, const PRODXM<CT,T2>& pxm)
{ return SUMMX<CT,T2>(-pxm.GetX1(),pxm.GetM2(),CT(x3)); }

template <class T, class T2> inline SUMMX<CT,T2> operator-(
    const VCT x3, const PRODXM<CT,T2>& pxm)
{ return SUMMX<CT,T2>(-pxm.GetX1(),pxm.GetM2(),CT(x3)); }

// (x*m)-x
template <class T, class T2> inline SUMMX<T,T2> operator-(
    const PRODXM<T,T2>& pxm, const T x3)
{ return SUMMX<T,T2>(pxm.GetX1(),pxm.GetM2(),-x3); }

template <class T, class T2> inline SUMMX<CT,T2> operator-(
    const PRODXM<CT,T2>& pxm, const T x3)
{ return SUMMX<CT,T2>(pxm.GetX1(),pxm.GetM2(),-x3); }

template <class T> inline SUMMX<CT,T> operator-(
    const PRODXM<T,T>& pxm, const CT x3)
{ return SUMMX<CT,T>(pxm.GetX1(),pxm.GetM2(),-x3); }

template <class T> inline SUMMX<CT,T> operator-(
    const PRODXM<T,T>& pxm, const CCT x3)
{ return SUMMX<CT,T>(pxm.GetX1(),pxm.GetM2(),-CT(x3)); }

template <class T> inline SUMMX<CT,T> operator-(
    const PRODXM<T,T>& pxm, const VCT x3)
{ return SUMMX<CT,T>(pxm.GetX1(),pxm.GetM2(),-CT(x3)); }

template <class T, class T2> inline SUMMX<CT,T2> operator-(
    const PRODXM<CT,T2>& pxm, const CCT x3)
{ return SUMMX<CT,T2>(pxm.GetX1(),pxm.GetM2(),-CT(x3)); }

template <class T, class T2> inline SUMMX<CT,T2> operator-(
    const PRODXM<CT,T2>& pxm, const VCT x3)
{ return SUMMX<CT,T2>(pxm.GetX1(),pxm.GetM2(),-CT(x3)); }

