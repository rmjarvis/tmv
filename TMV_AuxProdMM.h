
// Need to define the following with #define statements.
// (The given definition is for a regular Matrix*Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX1 GenMatrix
// #define GENMATRIX2 GenMatrix
// #define PRODMM ProdMM
// #define PRODXM1 ProdXM
// #define PRODXM2 ProdXM

// m*m
template <class T> inline PRODMM<T,T,T> operator*(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& v2)
{ return PRODMM<T,T,T>(T(1),m1,v2); }

template <class T> inline PRODMM<CT,CT,T> operator*(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& v2)
{ return PRODMM<CT,CT,T>(CT(1),m1,v2); }

template <class T> inline PRODMM<CT,T,CT> operator*(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& v2)
{ return PRODMM<CT,T,CT>(CT(1),m1,v2); }

// x*(x*m*m)
template <class T, class T1, class T2> inline PRODMM<T,T1,T2> operator*(
    T x1, const PRODMM<T,T1,T2>& pmv2)
{ return PRODMM<T,T1,T2>(x1*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetM3()); }

template <class T> inline PRODMM<CT,T,T> operator*(
    CT x1, const PRODMM<T,T,T>& pmv2)
{ return PRODMM<CT,T,T>(x1*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetM3()); }

template <class T> inline PRODMM<CT,T,T> operator*(
    CCT x1, const PRODMM<T,T,T>& pmv2)
{ return PRODMM<CT,T,T>(CT(x1)*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetM3()); }

template <class T> inline PRODMM<CT,T,T> operator*(
    VCT x1, const PRODMM<T,T,T>& pmv2)
{ return PRODMM<CT,T,T>(CT(x1)*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetM3()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator*(
    T x1, const PRODMM<CT,T1,T2>& pmv2)
{ return PRODMM<CT,T1,T2>(x1*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetM3()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator*(
    CCT x1, const PRODMM<CT,T1,T2>& pmv2)
{ return PRODMM<CT,T1,T2>(CT(x1)*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetM3()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator*(
    VCT x1, const PRODMM<CT,T1,T2>& pmv2)
{ return PRODMM<CT,T1,T2>(CT(x1)*pmv2.GetX1(),pmv2.GetM2(),pmv2.GetM3()); }

// (x*m*m)*x
template <class T, class T1, class T2> inline PRODMM<T,T1,T2> operator*(
    const PRODMM<T,T1,T2>& pmv1, T x2)
{ return PRODMM<T,T1,T2>(x2*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3()); }

template <class T> inline PRODMM<CT,T,T> operator*(
    const PRODMM<T,T,T>& pmv1, CT x2)
{ return PRODMM<CT,T,T>(x2*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3()); }

template <class T> inline PRODMM<CT,T,T> operator*(
    const PRODMM<T,T,T>& pmv1, CCT x2)
{ return PRODMM<CT,T,T>(CT(x2)*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3()); }

template <class T> inline PRODMM<CT,T,T> operator*(
    const PRODMM<T,T,T>& pmv1, VCT x2)
{ return PRODMM<CT,T,T>(CT(x2)*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator*(
    const PRODMM<CT,T1,T2>& pmv1, T x2)
{ return PRODMM<CT,T1,T2>(x2*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator*(
    const PRODMM<CT,T1,T2>& pmv1, CCT x2)
{ return PRODMM<CT,T1,T2>(CT(x2)*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator*(
    const PRODMM<CT,T1,T2>& pmv1, VCT x2)
{ return PRODMM<CT,T1,T2>(CT(x2)*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3()); }

// m*(x*m)
template <class T, class T2> inline PRODMM<T,T,T2> operator*(
    const GENMATRIX1<T>& m1, const PRODXM2<T,T2>& pxm2)
{ return PRODMM<T,T,T>(pxm2.GetX1(),m1,pxm2.GetM2()); }

template <class T> inline PRODMM<CT,CT,T> operator*(
    const GENMATRIX1<CT>& m1, const PRODXM2<T,T>& pxm2)
{ return PRODMM<CT,CT,T>(pxm2.GetX1(),m1,pxm2.GetM2()); }

template <class T, class T2> inline PRODMM<CT,T,T2> operator*(
    const GENMATRIX1<T>& m1, const PRODXM2<CT,T2>& pxm2)
{ return PRODMM<CT,T,T2>(pxm2.GetX1(),m1,pxm2.GetM2()); }

// (x*m)*m
template <class T, class T1> inline PRODMM<T,T1,T> operator*(
    const PRODXM1<T,T1>& pxm1, const GENMATRIX2<T>& v2)
{ return PRODMM<T,T1,T>(pxm1.GetX1(),pxm1.GetM2(),v2); }

template <class T> inline PRODMM<CT,T,CT> operator*(
    const PRODXM1<T,T>& pxm1, const GENMATRIX2<CT>& v2)
{ return PRODMM<CT,T,CT>(pxm1.GetX1(),pxm1.GetM2(),v2); }

template <class T, class T1> inline PRODMM<CT,T1,T> operator*(
    const PRODXM1<CT,T1>& pxm1, const GENMATRIX2<T>& v2)
{ return PRODMM<CT,T1,T>(pxm1.GetX1(),pxm1.GetM2(),v2); }

// (x*m)*(x*m)
template <class T, class T1, class T2> inline PRODMM<T,T1,T2> operator*(
    const PRODXM1<T,T1>& pxm1, const PRODXM2<T,T2>& pxm2) 
{ return PRODMM<T,T1,T2>(pxm1.GetX1()*pxm2.GetX1(),pxm1.GetM2(),pxm2.GetM2()); }

template <class T, class T1> inline PRODMM<CT,T1,T> operator*(
    const PRODXM1<CT,T1>& pxm1, const PRODXM2<T,T>& pxm2) 
{ return PRODMM<CT,T1,T>(pxm1.GetX1()*pxm2.GetX1(),pxm1.GetM2(),pxm2.GetM2()); }

template <class T, class T2> inline PRODMM<CT,T,T2> operator*(
    const PRODXM1<T,T>& pxm1, const PRODXM2<CT,T2>& pxm2) 
{ return PRODMM<CT,T,T2>(pxm1.GetX1()*pxm2.GetX1(),pxm1.GetM2(),pxm2.GetM2()); }

