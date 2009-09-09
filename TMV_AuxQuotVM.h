// Need to define the following with #define statements.
// (The given definition is for a regular Matrix/Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENVECTOR1 GenVecotx
// #define GENMATRIX2 GenMatrix
// #define QUOTVM QuotVM
// #define RQUOTVM RQuotVM

template <class T> inline QUOTVM<T,T,T> operator/(
    const GENVECTOR1<T>& m1, const GENMATRIX2<T>& m2)
{ return QUOTVM<T,T,T>(m1,m2); }

template <class T> inline QUOTVM<CT,CT,T> operator/(
    const GENVECTOR1<CT>& m1, const GENMATRIX2<T>& m2)
{ return QUOTVM<CT,CT,T>(m1,m2); }

template <class T> inline QUOTVM<CT,T,CT> operator/(
    const GENVECTOR1<T>& m1, const GENMATRIX2<CT>& m2)
{ return QUOTVM<CT,T,CT>(m1,m2); }

template <class T> inline RQUOTVM<T,T,T> operator%(
    const GENVECTOR1<T>& m1, const GENMATRIX2<T>& m2)
{ return RQUOTVM<T,T,T>(m1,m2); }

template <class T> inline RQUOTVM<CT,CT,T> operator%(
    const GENVECTOR1<CT>& m1, const GENMATRIX2<T>& m2)
{ return RQUOTVM<CT,CT,T>(m1,m2); }

template <class T> inline RQUOTVM<CT,T,CT> operator%(
    const GENVECTOR1<T>& m1, const GENMATRIX2<CT>& m2)
{ return RQUOTVM<CT,T,CT>(m1,m2); }

