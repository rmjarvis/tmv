// Need to define the following with #define statements.
// (The given definition is for a regular Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX2 GenMatrix
// #define QUOTXM QuotXM

// x/m
template <class T> inline QUOTXM<T,T,T> operator/(
    T x1, const GENMATRIX2<T>& m2)
{ return QUOTXM<T,T,T>(x1,m2); }

template <class T> inline QUOTXM<CT,T,CT> operator/(
    T x1, const GENMATRIX2<CT>& m2)
{ return QUOTXM<CT,T,CT>(x1,m2); }

template <class T> inline QUOTXM<CT,CT,T> operator/(
    CT x1, const GENMATRIX2<T>& m2)
{ return QUOTXM<CT,CT,T>(x1,m2); }

template <class T> inline QUOTXM<CT,CT,T> operator/(
    CCT x1, const GENMATRIX2<T>& m2)
{ return QUOTXM<CT,CT,T>(CT(x1),m2); }

template <class T> inline QUOTXM<CT,CT,T> operator/(
    VCT x1, const GENMATRIX2<T>& m2)
{ return QUOTXM<CT,CT,T>(CT(x1),m2); }

template <class T> inline QUOTXM<CT,CT,CT> operator/(
    CCT x1, const GENMATRIX2<CT>& m2)
{ return QUOTXM<CT,CT,CT>(CT(x1),m2); }

template <class T> inline QUOTXM<CT,CT,CT> operator/(
    VCT x1, const GENMATRIX2<CT>& m2)
{ return QUOTXM<CT,CT,CT>(CT(x1),m2); }

// x%m
template <class T> inline QUOTXM<T,T,T> operator%(
    T x1, const GENMATRIX2<T>& m2)
{ return QUOTXM<T,T,T>(x1,m2); }

template <class T> inline QUOTXM<CT,T,CT> operator%(
    T x1, const GENMATRIX2<CT>& m2)
{ return QUOTXM<CT,T,CT>(x1,m2); }

template <class T> inline QUOTXM<CT,CT,T> operator%(
    CT x1, const GENMATRIX2<T>& m2)
{ return QUOTXM<CT,CT,T>(x1,m2); }

template <class T> inline QUOTXM<CT,CT,T> operator%(
    CCT x1, const GENMATRIX2<T>& m2)
{ return QUOTXM<CT,CT,T>(CT(x1),m2); }

template <class T> inline QUOTXM<CT,CT,T> operator%(
    VCT x1, const GENMATRIX2<T>& m2)
{ return QUOTXM<CT,CT,T>(CT(x1),m2); }

template <class T> inline QUOTXM<CT,CT,CT> operator%(
    CCT x1, const GENMATRIX2<CT>& m2)
{ return QUOTXM<CT,CT,CT>(CT(x1),m2); }

template <class T> inline QUOTXM<CT,CT,CT> operator%(
    VCT x1, const GENMATRIX2<CT>& m2)
{ return QUOTXM<CT,CT,CT>(CT(x1),m2); }

