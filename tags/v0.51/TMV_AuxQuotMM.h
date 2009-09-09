// Need to define the following with #define statements.
// (The given definition is for a regular Matrix/Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX1 GenMatrix
// #define GENMATRIX2 GenMatrix
// #define QUOTMM QuotMM
// #define RQUOTMM RQuotMM

template <class T> inline QUOTMM<T,T,T> operator/(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{ return QUOTMM<T,T,T>(m1,m2); }

template <class T> inline QUOTMM<CT,CT,T> operator/(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& m2)
{ return QUOTMM<CT,CT,T>(m1,m2); }

template <class T> inline QUOTMM<CT,T,CT> operator/(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& m2)
{ return QUOTMM<CT,T,CT>(m1,m2); }

template <class T> inline RQUOTMM<T,T,T> operator%(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{ return RQUOTMM<T,T,T>(m1,m2); }

template <class T> inline RQUOTMM<CT,CT,T> operator%(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& m2)
{ return RQUOTMM<CT,CT,T>(m1,m2); }

template <class T> inline RQUOTMM<CT,T,CT> operator%(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& m2)
{ return RQUOTMM<CT,T,CT>(m1,m2); }

