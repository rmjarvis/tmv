// Need to define the following with #define statements.
// (The given definition is for a regular Matrix/Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX1 GenMatrix
// #define GENMATRIX2 GenMatrix
// #define RQUOTMM QuotMM
// #define RQUOTMM RQuotMM
// #define PRODXM1 ProdXM
// #define PRODXM2 ProdXM
// #define QUOTXM1 QuotXM
// #define QUOTXM2 QuotXM

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

// x*(x*m/m)
template <class T, class T1, class T2> inline QUOTMM<T,T1,T2> operator*(
    T x, const QUOTMM<T,T1,T2>& qmm)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    CT x, const QUOTMM<CT,T,T>& qmm)
{ return QUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    CCT x, const QUOTMM<CT,T,T>& qmm)
{ return QUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    VCT x, const QUOTMM<CT,T,T>& qmm)
{ return QUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    T x, const QUOTMM<CT,T1,T2>& qmm)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    CCT x, const QUOTMM<CT,T1,T2>& qmm)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    VCT x, const QUOTMM<CT,T1,T2>& qmm)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

// (x*m/m)*x
template <class T, class T1, class T2> inline QUOTMM<T,T1,T2> operator*(
    const QUOTMM<T,T1,T2>& qmm, T x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    const QUOTMM<CT,T,T>& qmm, CT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    const QUOTMM<CT,T,T>& qmm, CCT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    const QUOTMM<CT,T,T>& qmm, VCT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    const QUOTMM<CT,T1,T2>& qmm, T x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    const QUOTMM<CT,T1,T2>& qmm, CCT x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    const QUOTMM<CT,T1,T2>& qmm, VCT x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

// (x*m/m)/x
template <class T, class T1, class T2> inline QUOTMM<T,T1,T2> operator/(
    const QUOTMM<T,T1,T2>& qmm, T x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator/(
    const QUOTMM<CT,T,T>& qmm, CT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator/(
    const QUOTMM<CT,T,T>& qmm, CCT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator/(
    const QUOTMM<CT,T,T>& qmm, VCT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator/(
    const QUOTMM<CT,T1,T2>& qmm, T x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator/(
    const QUOTMM<CT,T1,T2>& qmm, CCT x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator/(
    const QUOTMM<CT,T1,T2>& qmm, VCT x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }

// x*(x*m%m)
template <class T, class T1, class T2> inline RQUOTMM<T,T1,T2> operator*(
    T x, const RQUOTMM<T,T1,T2>& qmm)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    CT x, const RQUOTMM<CT,T,T>& qmm)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    CCT x, const RQUOTMM<CT,T,T>& qmm)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    VCT x, const RQUOTMM<CT,T,T>& qmm)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    T x, const RQUOTMM<CT,T1,T2>& qmm)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    CCT x, const RQUOTMM<CT,T1,T2>& qmm)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    VCT x, const RQUOTMM<CT,T1,T2>& qmm)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

// (x*m%m)*x
template <class T, class T1, class T2> inline RQUOTMM<T,T1,T2> operator*(
    const RQUOTMM<T,T1,T2>& qmm, T x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    const RQUOTMM<CT,T,T>& qmm, CT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    const RQUOTMM<CT,T,T>& qmm, CCT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    const RQUOTMM<CT,T,T>& qmm, VCT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    const RQUOTMM<CT,T1,T2>& qmm, T x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    const RQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    const RQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

// (x*m%m)*x
template <class T, class T1, class T2> inline RQUOTMM<T,T1,T2> operator/(
    const RQUOTMM<T,T1,T2>& qmm, T x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator/(
    const RQUOTMM<CT,T,T>& qmm, CT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator/(
    const RQUOTMM<CT,T,T>& qmm, CCT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator/(
    const RQUOTMM<CT,T,T>& qmm, VCT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator/(
    const RQUOTMM<CT,T1,T2>& qmm, T x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator/(
    const RQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator/(
    const RQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }

// (x*m)/m
template <class T, class T1> inline QUOTMM<T,T1,T> operator/(
    const PRODXM1<T,T1>& pxm, const GENMATRIX2<T>& m)
{ return QUOTMM<T,T1,T>(pxm.GetX(),pxm.GetM(),m); }

template <class T> inline QUOTMM<CT,T,CT> operator/(
    const PRODXM1<T,T>& pxm, const GENMATRIX2<CT>& m)
{ return QUOTMM<CT,T,CT>(pxm.GetX(),pxm.GetM(),m); }

template <class T, class T1> inline QUOTMM<CT,T1,T> operator/(
    const PRODXM1<CT,T1>& pxm, const GENMATRIX2<T>& m)
{ return QUOTMM<CT,T1,T>(pxm.GetX(),m,pxm.GetM(),m); }

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

