// Need to define the following with #define statements.
// (The given definition is for a regular Matrix/Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX1 GenMatrix
// #define GENMATRIX2 GenMatrix
// #define QUOTMM QuotMM
// #define RQUOTMM RQuotMM
// #define PRODXM1 ProdXM
// #define PRODXM2 ProdXM
// #define QUOTXM QuotXM

// m/m
template <class T> inline QUOTMM<T,T,T> operator/(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{ return QUOTMM<T,T,T>(T(1),m1,m2); }

template <class T> inline QUOTMM<CT,CT,T> operator/(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& m2)
{ return QUOTMM<CT,CT,T>(CT(1),m1,m2); }

template <class T> inline QUOTMM<CT,T,CT> operator/(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& m2)
{ return QUOTMM<CT,T,CT>(CT(1),m1,m2); }

// m%m
template <class T> inline RQUOTMM<T,T,T> operator%(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{ return RQUOTMM<T,T,T>(T(1),m1,m2); }

template <class T> inline RQUOTMM<CT,CT,T> operator%(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& m2)
{ return RQUOTMM<CT,CT,T>(CT(1),m1,m2); }

template <class T> inline RQUOTMM<CT,T,CT> operator%(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& m2)
{ return RQUOTMM<CT,T,CT>(CT(1),m1,m2); }

// (x*m)%m
template <class T, class T1> inline RQUOTMM<T,T1,T> operator%(
    const PRODXM1<T,T1>& pxm, const GENMATRIX2<T>& m)
{ return RQUOTMM<T,T1,T>(pxm.GetX(),pxm.GetM(),m); }

template <class T> inline RQUOTMM<CT,T,CT> operator%(
    const PRODXM1<T,T>& pxm, const GENMATRIX2<CT>& m)
{ return RQUOTMM<CT,T,CT>(pxm.GetX(),pxm.GetM(),m); }

template <class T, class T1> inline RQUOTMM<CT,T1,T> operator%(
    const PRODXM1<CT,T1>& pxm, const GENMATRIX2<T>& m)
{ return RQUOTMM<CT,T1,T>(pxm.GetX(),pxm.GetM(),m); }

// m/(x*m)
template <class T, class T2> inline QUOTMM<T,T,T2> operator/(
    const GENMATRIX1<T>& m, const PRODXM2<T,T2>& pxm)
{ return QUOTMM<T,T,T2>(RealType(T)(1)/pxm.GetX(),m,pxm.GetM()); }

template <class T> inline QUOTMM<CT,CT,T> operator/(
    const GENMATRIX1<CT>& m, const PRODXM2<T,T>& pxm)
{ return QUOTMM<CT,CT,T>(RealType(T)(1)/pxm.GetX(),m,pxm.GetM()); }

template <class T, class T2> inline QUOTMM<CT,T,T2> operator/(
    const GENMATRIX1<T>& m, const PRODXM2<CT,T2>& pxm)
{ return QUOTMM<CT,T,T2>(RealType(T)(1)/pxm.GetX(),m,pxm.GetM()); }

// m%(x*m)
template <class T, class T2> inline QUOTMM<T,T,T2> operator%(
    const GENMATRIX1<T>& m, const PRODXM2<T,T2>& pxm)
{ return RQUOTMM<T,T,T2>(RealType(T)(1)/pxm.GetX(),m,pxm.GetM()); }

template <class T> inline QUOTMM<CT,CT,T> operator%(
    const GENMATRIX1<CT>& m, const PRODXM2<T,T>& pxm)
{ return RQUOTMM<CT,CT,T>(RealType(T)(1)/pxm.GetX(),m,pxm.GetM()); }

template <class T, class T2> inline QUOTMM<CT,T,T2> operator%(
    const GENMATRIX1<T>& m, const PRODXM2<CT,T2>& pxm)
{ return RQUOTMM<CT,T,T2>(RealType(T)(1)/pxm.GetX(),m,pxm.GetM()); }

// (x*m)/(x*m)
template <class T, class T1, class T2> inline QUOTMM<T,T1,T2> operator/(
    const PRODXM1<T,T1>& pxm1, const PRODXM2<T,T2>& pxm2)
{ return QUOTMM<T,T1,T2>(pxm1.GetX()/pxm2.GetX(),pxm1.GetM(),pxm2.GetM()); }

template <class T, class T1> inline QUOTMM<CT,T1,T> operator/(
    const PRODXM1<CT,T1>& pxm1, const PRODXM2<T,T>& pxm2)
{ return QUOTMM<CT,T1,T>(pxm1.GetX()/pxm2.GetX(),pxm1.GetM(),pxm2.GetM()); }

template <class T, class T2> inline QUOTMM<CT,T,T2> operator/(
    const PRODXM1<T,T>& pxm1, const PRODXM2<CT,T2>& pxm2)
{ return QUOTMM<CT,T,T2>(pxm1.GetX()/pxm2.GetX(),pxm1.GetM(),pxm2.GetM()); }

// (x*m)%(x*m)
template <class T, class T1, class T2> inline QUOTMM<T,T1,T2> operator%(
    const PRODXM1<T,T1>& pxm1, const PRODXM2<T,T2>& pxm2)
{ return RQUOTMM<T,T1,T2>(pxm1.GetX()/pxm2.GetX(),pxm1.GetM(),pxm2.GetM()); }

template <class T, class T1> inline QUOTMM<CT,T1,T> operator%(
    const PRODXM1<CT,T1>& pxm1, const PRODXM2<T,T>& pxm2)
{ return RQUOTMM<CT,T1,T>(pxm1.GetX()/pxm2.GetX(),pxm1.GetM(),pxm2.GetM()); }

template <class T, class T2> inline QUOTMM<CT,T,T2> operator%(
    const PRODXM1<T,T>& pxm1, const PRODXM2<CT,T2>& pxm2)
{ return RQUOTMM<CT,T,T2>(pxm1.GetX()/pxm2.GetX(),pxm1.GetM(),pxm2.GetM()); }

// (x/m)*m
template <class T, class T2> inline QUOTMM<T,T,T2> operator*(
    const QUOTXM<T,T2>& qxm, const GENMATRIX1<T>& m)
{ return QUOTMM<T,T,T2>(qxm.GetX(),m,qxm.GetM()); }

template <class T> inline QUOTMM<CT,CT,T> operator*(
    const QUOTXM<T,T>& qxm, const GENMATRIX1<CT>& m)
{ return QUOTMM<CT,CT,T>(qxm.GetX(),m,qxm.GetM()); }

template <class T, class T2> inline QUOTMM<CT,T,T2> operator*(
    const QUOTXM<CT,T2>& qxm, const GENMATRIX1<T>& m)
{ return QUOTMM<CT,T,T2>(qxm.GetX(),m,qxm.GetM()); }

// m*(x/m)
template <class T, class T2> inline RQUOTMM<T,T,T2> operator*(
    const GENMATRIX1<T>& m, const QUOTXM<T,T2>& qxm)
{ return RQUOTMM<T,T,T2>(qxm.GetX(),m,qxm.GetM()); }

template <class T> inline RQUOTMM<CT,CT,T> operator*(
    const GENMATRIX1<CT>& m, const QUOTXM<T,T>& qxm)
{ return RQUOTMM<CT,CT,T>(qxm.GetX(),m,qxm.GetM()); }

template <class T, class T2> inline RQUOTMM<CT,T,T2> operator*(
    const GENMATRIX1<T>& m, const QUOTXM<CT,T2>& qxm)
{ return RQUOTMM<CT,T,T2>(qxm.GetX(),m,qxm.GetM()); }

// (x/m)*(x*m)
template <class T, class T1, class T2> inline QUOTMM<T,T,T2> operator*(
    const QUOTXM<T,T2>& qxm, const PRODXM1<T,T1>& pxm)
{ return QUOTMM<T,T,T2>(pxm.GetX()*qxm.GetX(),pxm.GetM(),qxm.GetM()); }

template <class T, class T1> inline QUOTMM<CT,CT,T> operator*(
    const QUOTXM<T,T>& qxm, const PRODXM1<CT,T1>& pxm)
{ return QUOTMM<CT,CT,T>(pxm.GetX()*qxm.GetX(),pxm.GetM(),qxm.GetM()); }

template <class T, class T2> inline QUOTMM<CT,T,T2> operator*(
    const QUOTXM<CT,T2>& qxm, const PRODXM1<T,T>& pxm)
{ return QUOTMM<CT,T,T2>(pxm.GetX()*qxm.GetX(),pxm.GetM(),qxm.GetM()); }

// (x*m)*(x/m)
template <class T, class T1, class T2> inline RQUOTMM<T,T1,T2> operator*(
    const PRODXM1<T,T1>& pxm, const QUOTXM<T,T2>& qxm)
{ return RQUOTMM<T,T1,T2>(pxm.GetX()*qxm.GetX(),pxm.GetM(),qxm.GetM()); }

template <class T, class T1> inline RQUOTMM<CT,T1,T> operator*(
    const PRODXM1<CT,T1>& pxm, const QUOTXM<T,T>& qxm)
{ return RQUOTMM<CT,T1,T>(pxm.GetX()*qxm.GetX(),pxm.GetM(),qxm.GetM()); }

template <class T, class T2> inline RQUOTMM<CT,T,T2> operator*(
    const PRODXM1<T,T>& pxm, const QUOTXM<CT,T2>& qxm)
{ return RQUOTMM<CT,T,T2>(pxm.GetX()*qxm.GetX(),pxm.GetM(),qxm.GetM()); }

