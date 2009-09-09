// Need to define the following with #define statements.
// (The given definition is for a regular Matrix/Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENVECTOR GenVecotx
// #define GENMATRIX GenMatrix
// #define QUOTVM QuotVM
// #define RQUOTVM RQuotVM
// #define PRODXV ProdXV
// #define PRODXM ProdXM
// #define QUOTXM QuotXM

// v/m
template <class T> inline QUOTVM<T,T,T> operator/(
    const GENVECTOR<T>& v, const GENMATRIX<T>& m)
{ return QUOTVM<T,T,T>(T(1),v,m); }

template <class T> inline QUOTVM<CT,CT,T> operator/(
    const GENVECTOR<CT>& v, const GENMATRIX<T>& m)
{ return QUOTVM<CT,CT,T>(CT(1),v,m); }

template <class T> inline QUOTVM<CT,T,CT> operator/(
    const GENVECTOR<T>& v, const GENMATRIX<CT>& m)
{ return QUOTVM<CT,T,CT>(CT(1),v,m); }

// v%m
template <class T> inline RQUOTVM<T,T,T> operator%(
    const GENVECTOR<T>& v, const GENMATRIX<T>& m)
{ return RQUOTVM<T,T,T>(T(1),v,m); }

template <class T> inline RQUOTVM<CT,CT,T> operator%(
    const GENVECTOR<CT>& v, const GENMATRIX<T>& m)
{ return RQUOTVM<CT,CT,T>(CT(1),v,m); }

template <class T> inline RQUOTVM<CT,T,CT> operator%(
    const GENVECTOR<T>& v, const GENMATRIX<CT>& m)
{ return RQUOTVM<CT,T,CT>(CT(1),v,m); }

// x*(x*v/m)
template <class T, class T1, class T2> inline QUOTVM<T,T1,T2> operator*(
    T x, const QUOTVM<T,T1,T2>& qvm)
{ return QUOTVM<T,T1,T2>(x*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T> inline QUOTVM<CT,T,T> operator*(
    CT x, const QUOTVM<T,T,T>& qvm)
{ return QUOTVM<CT,T,T>(x*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T> inline QUOTVM<CT,T,T> operator*(
    CCT x, const QUOTVM<T,T,T>& qvm)
{ return QUOTVM<CT,T,T>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T> inline QUOTVM<CT,T,T> operator*(
    VCT x, const QUOTVM<T,T,T>& qvm)
{ return QUOTVM<CT,T,T>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline QUOTVM<CT,T1,T2> operator*(
    T x, const QUOTVM<CT,T1,T2>& qvm)
{ return QUOTVM<T,T1,T2>(x*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline QUOTVM<CT,T1,T2> operator*(
    CCT x, const QUOTVM<CT,T1,T2>& qvm)
{ return QUOTVM<T,T1,T2>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline QUOTVM<CT,T1,T2> operator*(
    VCT x, const QUOTVM<CT,T1,T2>& qvm)
{ return QUOTVM<T,T1,T2>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

// (x*v/m)*x
template <class T, class T1, class T2> inline QUOTVM<T,T1,T2> operator*(
    const QUOTVM<T,T1,T2>& qvm, T x)
{ return QUOTVM<T,T1,T2>(x*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T> inline QUOTVM<CT,T,T> operator*(
    const QUOTVM<T,T,T>& qvm, CT x)
{ return QUOTVM<CT,T,T>(x*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T> inline QUOTVM<CT,T,T> operator*(
    const QUOTVM<T,T,T>& qvm, CCT x)
{ return QUOTVM<CT,T,T>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T> inline QUOTVM<CT,T,T> operator*(
    const QUOTVM<T,T,T>& qvm, VCT x)
{ return QUOTVM<CT,T,T>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline QUOTVM<CT,T1,T2> operator*(
    const QUOTVM<CT,T1,T2>& qvm, T x)
{ return QUOTVM<T,T1,T2>(x*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline QUOTVM<CT,T1,T2> operator*(
    const QUOTVM<CT,T1,T2>& qvm, CCT x)
{ return QUOTVM<T,T1,T2>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline QUOTVM<CT,T1,T2> operator*(
    const QUOTVM<CT,T1,T2>& qvm, VCT x)
{ return QUOTVM<T,T1,T2>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

// (x*v/m)/x
template <class T, class T1, class T2> inline QUOTVM<T,T1,T2> operator/(
    const QUOTVM<T,T1,T2>& qvm, T x)
{ return QUOTVM<T,T1,T2>(qvm.GetX()/x,qvm.GetV(),qvm.GetM()); }

template <class T> inline QUOTVM<CT,T,T> operator/(
    const QUOTVM<T,T,T>& qvm, CT x)
{ return QUOTVM<CT,T,T>(qvm.GetX()/x,qvm.GetV(),qvm.GetM()); }

template <class T> inline QUOTVM<CT,T,T> operator/(
    const QUOTVM<T,T,T>& qvm, CCT x)
{ return QUOTVM<CT,T,T>(qvm.GetX()/CT(x),qvm.GetV(),qvm.GetM()); }

template <class T> inline QUOTVM<CT,T,T> operator/(
    const QUOTVM<T,T,T>& qvm, VCT x)
{ return QUOTVM<CT,T,T>(qvm.GetX()/CT(x),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline QUOTVM<CT,T1,T2> operator/(
    const QUOTVM<CT,T1,T2>& qvm, T x)
{ return QUOTVM<T,T1,T2>(qvm.GetX()/x,qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline QUOTVM<CT,T1,T2> operator/(
    const QUOTVM<CT,T1,T2>& qvm, CCT x)
{ return QUOTVM<T,T1,T2>(qvm.GetX()/CT(x),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline QUOTVM<CT,T1,T2> operator/(
    const QUOTVM<CT,T1,T2>& qvm, VCT x)
{ return QUOTVM<T,T1,T2>(qvm.GetX()/CT(x),qvm.GetV(),qvm.GetM()); }

// x*(x*v%m)
template <class T, class T1, class T2> inline RQUOTVM<T,T1,T2> operator*(
    T x, const RQUOTVM<T,T1,T2>& qvm)
{ return RQUOTVM<T,T1,T2>(x*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T> inline RQUOTVM<CT,T,T> operator*(
    CT x, const RQUOTVM<T,T,T>& qvm)
{ return RQUOTVM<CT,T,T>(x*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T> inline RQUOTVM<CT,T,T> operator*(
    CCT x, const RQUOTVM<T,T,T>& qvm)
{ return RQUOTVM<CT,T,T>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T> inline RQUOTVM<CT,T,T> operator*(
    VCT x, const RQUOTVM<T,T,T>& qvm)
{ return RQUOTVM<CT,T,T>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline RQUOTVM<CT,T1,T2> operator*(
    T x, const RQUOTVM<CT,T1,T2>& qvm)
{ return RQUOTVM<T,T1,T2>(x*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline RQUOTVM<CT,T1,T2> operator*(
    CCT x, const RQUOTVM<CT,T1,T2>& qvm)
{ return RQUOTVM<T,T1,T2>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline RQUOTVM<CT,T1,T2> operator*(
    VCT x, const RQUOTVM<CT,T1,T2>& qvm)
{ return RQUOTVM<T,T1,T2>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

// (x*v%m)*x
template <class T, class T1, class T2> inline RQUOTVM<T,T1,T2> operator*(
    const RQUOTVM<T,T1,T2>& qvm, T x)
{ return RQUOTVM<T,T1,T2>(x*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T> inline RQUOTVM<CT,T,T> operator*(
    const RQUOTVM<T,T,T>& qvm, CT x)
{ return RQUOTVM<CT,T,T>(x*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T> inline RQUOTVM<CT,T,T> operator*(
    const RQUOTVM<T,T,T>& qvm, CCT x)
{ return RQUOTVM<CT,T,T>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T> inline RQUOTVM<CT,T,T> operator*(
    const RQUOTVM<T,T,T>& qvm, VCT x)
{ return RQUOTVM<CT,T,T>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline RQUOTVM<CT,T1,T2> operator*(
    const RQUOTVM<CT,T1,T2>& qvm, T x)
{ return RQUOTVM<T,T1,T2>(x*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline RQUOTVM<CT,T1,T2> operator*(
    const RQUOTVM<CT,T1,T2>& qvm, CCT x)
{ return RQUOTVM<T,T1,T2>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline RQUOTVM<CT,T1,T2> operator*(
    const RQUOTVM<CT,T1,T2>& qvm, VCT x)
{ return RQUOTVM<T,T1,T2>(CT(x)*qvm.GetX(),qvm.GetV(),qvm.GetM()); }

// (x*v%m)/x
template <class T, class T1, class T2> inline RQUOTVM<T,T1,T2> operator/(
    const RQUOTVM<T,T1,T2>& qvm, T x)
{ return RQUOTVM<T,T1,T2>(qvm.GetX()/x,qvm.GetV(),qvm.GetM()); }

template <class T> inline RQUOTVM<CT,T,T> operator/(
    const RQUOTVM<T,T,T>& qvm, CT x)
{ return RQUOTVM<CT,T,T>(qvm.GetX()/x,qvm.GetV(),qvm.GetM()); }

template <class T> inline RQUOTVM<CT,T,T> operator/(
    const RQUOTVM<T,T,T>& qvm, CCT x)
{ return RQUOTVM<CT,T,T>(qvm.GetX()/CT(x),qvm.GetV(),qvm.GetM()); }

template <class T> inline RQUOTVM<CT,T,T> operator/(
    const RQUOTVM<T,T,T>& qvm, VCT x)
{ return RQUOTVM<CT,T,T>(qvm.GetX()/CT(x),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline RQUOTVM<CT,T1,T2> operator/(
    const RQUOTVM<CT,T1,T2>& qvm, T x)
{ return RQUOTVM<T,T1,T2>(qvm.GetX()/x,qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline RQUOTVM<CT,T1,T2> operator/(
    const RQUOTVM<CT,T1,T2>& qvm, CCT x)
{ return RQUOTVM<T,T1,T2>(qvm.GetX()/CT(x),qvm.GetV(),qvm.GetM()); }

template <class T, class T1, class T2> inline RQUOTVM<CT,T1,T2> operator/(
    const RQUOTVM<CT,T1,T2>& qvm, VCT x)
{ return RQUOTVM<T,T1,T2>(qvm.GetX()/CT(x),qvm.GetV(),qvm.GetM()); }

// (x*v)/m
template <class T, class T1> inline QUOTVM<T,T1,T> operator/(
    const PRODXV<T,T1>& pxv, const GENMATRIX<T>& m)
{ return QUOTVM<T,T1,T>(pxv.GetX(),pxv.GetV(),m); }

template <class T> inline QUOTVM<CT,T,CT> operator/(
    const PRODXV<T,T>& pxv, const GENMATRIX<CT>& m)
{ return QUOTVM<CT,T,CT>(pxv.GetX(),pxv.GetV(),m); }

template <class T, class T1> inline QUOTVM<CT,T1,T> operator/(
    const PRODXV<CT,T1>& pxv, const GENMATRIX<T>& m)
{ return QUOTVM<CT,T1,T>(pxv.GetX(),pxv.GetV(),m); }

// (x*v)%m
template <class T, class T1> inline RQUOTVM<T,T1,T> operator%(
    const PRODXV<T,T1>& pxv, const GENMATRIX<T>& m)
{ return RQUOTVM<T,T1,T>(pxv.GetX(),pxv.GetV(),m); }

template <class T> inline RQUOTVM<CT,T,CT> operator%(
    const PRODXV<T,T>& pxv, const GENMATRIX<CT>& m)
{ return RQUOTVM<CT,T,CT>(pxv.GetX(),pxv.GetV(),m); }

template <class T, class T1> inline RQUOTVM<T,T1,T> operator%(
    const PRODXV<CT,T1>& pxv, const GENMATRIX<T>& m)
{ return RQUOTVM<CT,T1,T>(pxv.GetX(),pxv.GetV(),m); }

// v/(x*m)
template <class T, class T2> inline QUOTVM<T,T,T2> operator/(
    const GENVECTOR<T>& v, const PRODXM<T,T2>& pxm)
{ return QUOTVM<T,T,T2>(RealType(T)(1)/pxm.GetX(),v,pxm.GetM()); }

template <class T> inline QUOTVM<CT,CT,T> operator/(
    const GENVECTOR<CT>& v, const PRODXM<T,T>& pxm)
{ return QUOTVM<CT,CT,T>(T(1)/pxm.GetX(),v,pxm.GetM()); }

template <class T, class T2> inline QUOTVM<CT,T,T2> operator/(
    const GENVECTOR<T>& v, const PRODXM<CT,T2>& pxm)
{ return QUOTVM<CT,T,T2>(T(1)/pxm.GetX(),v,pxm.GetM()); }

// v%(x*m)
template <class T, class T2> inline QUOTVM<T,T,T2> operator%(
    const GENVECTOR<T>& v, const PRODXM<T,T2>& pxm)
{ return QUOTVM<T,T,T2>(RealType(T)(1)/pxm.GetX(),v,pxm.GetM()); }

template <class T> inline QUOTVM<CT,CT,T> operator%(
    const GENVECTOR<CT>& v, const PRODXM<T,T>& pxm)
{ return QUOTVM<CT,CT,T>(T(1)/pxm.GetX(),v,pxm.GetM()); }

template <class T, class T2> inline QUOTVM<CT,T,T2> operator%(
    const GENVECTOR<T>& v, const PRODXM<CT,T2>& pxm)
{ return QUOTVM<CT,T,T2>(T(1)/pxm.GetX(),v,pxm.GetM()); }

// (x*v)/(x*m)
template <class T, class T1, class T2> inline QUOTVM<T,T1,T2> operator/(
    const PRODXV<T,T1>& pxv, const PRODXM<T,T2>& pxm)
{ return QUOTVM<T,T1,T2>(pxv.GetX()/pxm.GetX(),pxv.GetV(),pxm.GetM()); }

template <class T, class T1> inline QUOTVM<CT,T1,T> operator/(
    const PRODXV<CT,T1>& pxv, const PRODXM<T,T>& pxm)
{ return QUOTVM<CT,T1,T>(pxv.GetX()/pxm.GetX(),pxv.GetV(),pxm.GetM()); }

template <class T, class T2> inline QUOTVM<CT,T,T2> operator/(
    const PRODXV<T,T>& pxv, const PRODXM<CT,T2>& pxm)
{ return QUOTVM<CT,T,T2>(pxv.GetX()/pxm.GetX(),pxv.GetV(),pxm.GetM()); }

// (x*v)%(x*m)
template <class T, class T1, class T2> inline RQUOTVM<T,T1,T2> operator%(
    const PRODXV<T,T1>& pxv, const PRODXM<T,T2>& pxm)
{ return RQUOTVM<T,T1,T2>(pxv.GetX()/pxm.GetX(),pxv.GetV(),pxm.GetM()); }

template <class T, class T1> inline RQUOTVM<CT,T1,T> operator%(
    const PRODXV<CT,T1>& pxv, const PRODXM<T,T>& pxm)
{ return RQUOTVM<CT,T1,T>(pxv.GetX()/pxm.GetX(),pxv.GetV(),pxm.GetM()); }

template <class T, class T2> inline RQUOTVM<CT,T,T2> operator%(
    const PRODXV<T,T>& pxv, const PRODXM<CT,T2>& pxm)
{ return RQUOTVM<CT,T,T2>(pxv.GetX()/pxm.GetX(),pxv.GetV(),pxm.GetM()); }

// v*(x/m)
template <class T, class T2> inline RQUOTVM<T,T,T2> operator*(
    const GENVECTOR<T>& v, const QUOTXM<T,T2>& pxm)
{ return RQUOTVM<T,T,T2>(pxm.GetX(),v,pxm.GetM()); }

template <class T> inline RQUOTVM<CT,CT,T> operator*(
    const GENVECTOR<CT>& v, const QUOTXM<T,T>& pxm)
{ return RQUOTVM<CT,CT,T>(pxm.GetX(),v,pxm.GetM()); }

template <class T, class T2> inline RQUOTVM<CT,T,T2> operator*(
    const GENVECTOR<T>& v, const QUOTXM<CT,T2>& pxm)
{ return RQUOTVM<CT,T,T2>(pxm.GetX(),v,pxm.GetM()); }

// (x/m)*v
template <class T, class T2> inline QUOTVM<T,T,T2> operator*(
    const QUOTXM<T,T2>& pxm, const GENVECTOR<T>& v) 
{ return QUOTVM<T,T,T2>(pxm.GetX(),v,pxm.GetM()); }

template <class T> inline QUOTVM<CT,CT,T> operator*(
    const QUOTXM<T,T>& pxm, const GENVECTOR<CT>& v) 
{ return QUOTVM<CT,CT,T>(pxm.GetX(),v,pxm.GetM()); }

template <class T, class T2> inline QUOTVM<CT,T,T2> operator*(
    const QUOTXM<CT,T2>& pxm, const GENVECTOR<T>& v) 
{ return QUOTVM<CT,T,T2>(pxm.GetX(),v,pxm.GetM()); }

// (x*v)*(x/m)
template <class T, class T1, class T2> inline RQUOTVM<T,T1,T2> operator*(
    const PRODXV<T,T1>& pxv, const QUOTXM<T,T2>& pxm)
{ return RQUOTVM<T,T1,T2>(pxv.GetX()*pxm.GetX(),pxv.GetV(),pxm.GetM()); }

template <class T, class T1> inline RQUOTVM<CT,T1,T> operator*(
    const PRODXV<CT,T1>& pxv, const QUOTXM<T,T>& pxm)
{ return RQUOTVM<CT,T1,T>(pxv.GetX()*pxm.GetX(),pxv.GetV(),pxm.GetM()); }

template <class T, class T2> inline RQUOTVM<CT,T,T2> operator*(
    const PRODXV<T,T>& pxv, const QUOTXM<CT,T2>& pxm)
{ return RQUOTVM<CT,T,T2>(pxv.GetX()*pxm.GetX(),pxv.GetV(),pxm.GetM()); }

// (x/m)*(x*v)
template <class T, class T1, class T2> inline QUOTVM<T,T1,T2> operator*(
    const QUOTXM<T,T2>& pxm, const PRODXV<T,T1>& pxv) 
{ return QUOTVM<T,T1,T2>(pxv.GetX()*pxm.GetX(),pxv.GetV(),pxm.GetM()); }

template <class T, class T1> inline QUOTVM<CT,T1,T> operator*(
    const QUOTXM<T,T>& pxm, const PRODXV<CT,T1>& pxv) 
{ return QUOTVM<CT,T1,T>(pxv.GetX()*pxm.GetX(),pxv.GetV(),pxm.GetM()); }

template <class T, class T2> inline QUOTVM<CT,T,T2> operator*(
    const QUOTXM<CT,T2>& pxm, const PRODXV<T,T>& pxv) 
{ return QUOTVM<CT,T,T2>(pxv.GetX()*pxm.GetX(),pxv.GetV(),pxm.GetM()); }


