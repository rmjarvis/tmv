// Need to define the following with #define statements.
// (The given definition is for a regular Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX GenMatrix
// #define QUOTXM QuotXM

// x/m
template <class T> inline QUOTXM<T,T> operator/(
    T x, const GENMATRIX<T>& m)
{ return QUOTXM<T,T>(x,m); }

template <class T> inline QUOTXM<CT,T> operator/(
    CT x, const GENMATRIX<T>& m)
{ return QUOTXM<CT,T>(x,m); }

template <class T> inline QUOTXM<CT,T> operator/(
    CCT x, const GENMATRIX<T>& m)
{ return QUOTXM<CT,T>(CT(x),m); }

template <class T> inline QUOTXM<CT,T> operator/(
    VCT x, const GENMATRIX<T>& m)
{ return QUOTXM<CT,T>(CT(x),m); }

template <class T> inline QUOTXM<CT,CT> operator/(
    T x, const GENMATRIX<CT>& m)
{ return QUOTXM<CT,CT>(CT(x),m); }

template <class T> inline QUOTXM<CT,CT> operator/(
    CCT x, const GENMATRIX<CT>& m)
{ return QUOTXM<CT,CT>(CT(x),m); }

template <class T> inline QUOTXM<CT,CT> operator/(
    VCT x, const GENMATRIX<CT>& m)
{ return QUOTXM<CT,CT>(CT(x),m); }

// x%m
template <class T> inline QUOTXM<T,T> operator%(
    T x, const GENMATRIX<T>& m)
{ return QUOTXM<T,T>(x,m); }

template <class T> inline QUOTXM<CT,T> operator%(
    CT x, const GENMATRIX<T>& m)
{ return QUOTXM<CT,T>(x,m); }

template <class T> inline QUOTXM<CT,T> operator%(
    CCT x, const GENMATRIX<T>& m)
{ return QUOTXM<CT,T>(CT(x),m); }

template <class T> inline QUOTXM<CT,T> operator%(
    VCT x, const GENMATRIX<T>& m)
{ return QUOTXM<CT,T>(CT(x),m); }

template <class T> inline QUOTXM<CT,CT> operator%(
    T x, const GENMATRIX<CT>& m)
{ return QUOTXM<CT,CT>(CT(x),m); }

template <class T> inline QUOTXM<CT,CT> operator%(
    CCT x, const GENMATRIX<CT>& m)
{ return QUOTXM<CT,CT>(CT(x),m); }

template <class T> inline QUOTXM<CT,CT> operator%(
    VCT x, const GENMATRIX<CT>& m)
{ return QUOTXM<CT,CT>(CT(x),m); }


// -(x/m)
template <class T, class T2> inline QUOTXM<T,T2> operator-(
    const QUOTXM<T,T2>& qxm)
{ return QUOTXM<T,T2>(-qxm.GetX(),qxm.GetM()); }

// x*(x/m)
template <class T, class T2> inline QUOTXM<T,T2> operator*(
    const T x, const QUOTXM<T,T2>& qxm)
{ return QUOTXM<T,T2>(x*qxm.GetX(),qxm.GetM()); }

template <class T, class T2> inline QUOTXM<CT,T2> operator*(
    const T x, const QUOTXM<CT,T2>& qxm)
{ return QUOTXM<CT,T2>(x*qxm.GetX(),qxm.GetM()); }

template <class T> inline QUOTXM<CT,T> operator*(
    const CT x, const QUOTXM<T,T>& qxm)
{ return QUOTXM<CT,T>(x*qxm.GetX(),qxm.GetM()); }

template <class T> inline QUOTXM<CT,T> operator*(
    const CCT x, const QUOTXM<T,T>& qxm)
{ return QUOTXM<CT,T>(CT(x)*qxm.GetX(),qxm.GetM()); }

template <class T> inline QUOTXM<CT,T> operator*(
    const VCT x, const QUOTXM<T,T>& qxm)
{ return QUOTXM<CT,T>(CT(x)*qxm.GetX(),qxm.GetM()); }

template <class T, class T2> inline QUOTXM<CT,T2> operator*(
    const CCT x, const QUOTXM<CT,T2>& qxm)
{ return QUOTXM<CT,T2>(CT(x)*qxm.GetX(),qxm.GetM()); }

template <class T, class T2> inline QUOTXM<CT,T2> operator*(
    const VCT x, const QUOTXM<CT,T2>& qxm)
{ return QUOTXM<CT,T2>(CT(x)*qxm.GetX(),qxm.GetM()); }

// (x/m)*x
template <class T, class T2> inline QUOTXM<T,T2> operator*(
    const QUOTXM<T,T2>& qxm, const T x)
{ return QUOTXM<T,T2>(x*qxm.GetX(),qxm.GetM()); }

template <class T, class T2> inline QUOTXM<CT,T2> operator*(
    const QUOTXM<CT,T2>& qxm, const T x)
{ return QUOTXM<CT,T2>(x*qxm.GetX(),qxm.GetM()); }

template <class T> inline QUOTXM<CT,T> operator*(
    const QUOTXM<T,T>& qxm, const CT x)
{ return QUOTXM<CT,T>(x*qxm.GetX(),qxm.GetM()); }

template <class T> inline QUOTXM<CT,T> operator*(
    const QUOTXM<T,T>& qxm, const CCT x)
{ return QUOTXM<CT,T>(CT(x)*qxm.GetX(),qxm.GetM()); }

template <class T> inline QUOTXM<CT,T> operator*(
    const QUOTXM<T,T>& qxm, const VCT x)
{ return QUOTXM<CT,T>(CT(x)*qxm.GetX(),qxm.GetM()); }

template <class T, class T2> inline QUOTXM<CT,T2> operator*(
    const QUOTXM<CT,T2>& qxm, const CCT x)
{ return QUOTXM<CT,T2>(CT(x)*qxm.GetX(),qxm.GetM()); }

template <class T, class T2> inline QUOTXM<CT,T2> operator*(
    const QUOTXM<CT,T2>& qxm, const VCT x)
{ return QUOTXM<CT,T2>(CT(x)*qxm.GetX(),qxm.GetM()); }

// (x/m)/x
template <class T, class T2> inline QUOTXM<T,T2> operator/(
    const QUOTXM<T,T2>& qxm, const T x)
{ return QUOTXM<T,T2>(qxm.GetX()/x,qxm.GetM()); }

template <class T, class T2> inline QUOTXM<CT,T2> operator/(
    const QUOTXM<CT,T2>& qxm, const T x)
{ return QUOTXM<CT,T2>(qxm.GetX()/x,qxm.GetM()); }

template <class T> inline QUOTXM<CT,T> operator/(
    const QUOTXM<T,T>& qxm, const CT x)
{ return QUOTXM<CT,T>(qxm.GetX()/x,qxm.GetM()); }

template <class T> inline QUOTXM<CT,T> operator/(
    const QUOTXM<T,T>& qxm, const CCT x)
{ return QUOTXM<CT,T>(qxm.GetX()/CT(x),qxm.GetM()); }

template <class T> inline QUOTXM<CT,T> operator/(
    const QUOTXM<T,T>& qxm, const VCT x)
{ return QUOTXM<CT,T>(qxm.GetX()/CT(x),qxm.GetM()); }

template <class T, class T2> inline QUOTXM<CT,T2> operator/(
    const QUOTXM<CT,T2>& qxm, const CCT x)
{ return QUOTXM<CT,T2>(qxm.GetX()/CT(x),qxm.GetM()); }

template <class T, class T2> inline QUOTXM<CT,T2> operator/(
    const QUOTXM<CT,T2>& qxm, const VCT x)
{ return QUOTXM<CT,T2>(qxm.GetX()/CT(x),qxm.GetM()); }


