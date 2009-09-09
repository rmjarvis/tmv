// Need to define the following with #define statements.
// (The given definition is for a DiagMatrix/Matrix.  Modify as 
// appropriate for the various other matrices.)
// Note: Vector/Matrix uses the same format - no need to write that separately.
//
// #define GENMATRIX1 GenDiagMatrix
// #define GENMATRIX2 GenMatrix
// #define TQUOTMM TransientQuotMM
// #define TRQUOTMM TransientRQuotMM

template <class T> inline TQUOTMM<T,T,ColMajor,T> operator/(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{ return TQUOTMM<T,T,ColMajor,T>(new Matrix<T,ColMajor>(m1),m2); }

template <class T> inline TQUOTMM<CT,CT,ColMajor,T> operator/(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& m2)
{ return TQUOTMM<CT,CT,ColMajor,T>(new Matrix<CT,ColMajor>(m1),m2); }

template <class T> inline TQUOTMM<CT,T,ColMajor,CT> operator/(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& m2)
{ return TQUOTMM<CT,T,ColMajor,CT>(new Matrix<T,ColMajor>(m1),m2); }

template <class T> inline TRQUOTMM<T,T,RowMajor,T> operator%(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{ return TRQUOTMM<T,T,RowMajor,T>(new Matrix<T,RowMajor>(m1),m2); }

template <class T> inline TRQUOTMM<CT,CT,RowMajor,T> operator%(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& m2)
{ return TRQUOTMM<CT,CT,RowMajor,T>(new Matrix<CT,RowMajor>(m1),m2); }

template <class T> inline TRQUOTMM<CT,T,RowMajor,CT> operator%(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& m2)
{ return TRQUOTMM<CT,T,RowMajor,CT>(new Matrix<T,RowMajor>(m1),m2); }

