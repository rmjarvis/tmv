
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
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{ return PRODMM<T,T,T>(T(1),m1,m2); }

template <class T> inline PRODMM<CT,CT,T> operator*(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& m2)
{ return PRODMM<CT,CT,T>(CT(1),m1,m2); }

template <class T> inline PRODMM<CT,T,CT> operator*(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& m2)
{ return PRODMM<CT,T,CT>(CT(1),m1,m2); }

// m*(x*m)
template <class T, class T2> inline PRODMM<T,T,T2> operator*(
    const GENMATRIX1<T>& m, const PRODXM2<T,T2>& pxm)
{ return PRODMM<T,T,T2>(pxm.GetX(),m,pxm.GetM()); }

template <class T> inline PRODMM<CT,CT,T> operator*(
    const GENMATRIX1<CT>& m, const PRODXM2<T,T>& pxm)
{ return PRODMM<CT,CT,T>(pxm.GetX(),m,pxm.GetM()); }

template <class T, class T2> inline PRODMM<CT,T,T2> operator*(
    const GENMATRIX1<T>& m, const PRODXM2<CT,T2>& pxm)
{ return PRODMM<CT,T,T2>(pxm.GetX(),m,pxm.GetM()); }

// (x*m)*m
template <class T, class T1> inline PRODMM<T,T1,T> operator*(
    const PRODXM1<T,T1>& pxm, const GENMATRIX2<T>& m)
{ return PRODMM<T,T1,T>(pxm.GetX(),pxm.GetM(),m); }

template <class T> inline PRODMM<CT,T,CT> operator*(
    const PRODXM1<T,T>& pxm, const GENMATRIX2<CT>& m)
{ return PRODMM<CT,T,CT>(pxm.GetX(),pxm.GetM(),m); }

template <class T, class T1> inline PRODMM<CT,T1,T> operator*(
    const PRODXM1<CT,T1>& pxm, const GENMATRIX2<T>& m)
{ return PRODMM<CT,T1,T>(pxm.GetX(),pxm.GetM(),m); }

// (x*m)*(x*m)
template <class T, class T1, class T2> inline PRODMM<T,T1,T2> operator*(
    const PRODXM1<T,T1>& pxm1, const PRODXM2<T,T2>& pxm2) 
{ return PRODMM<T,T1,T2>(pxm1.GetX()*pxm2.GetX(),pxm1.GetM(),pxm2.GetM()); }

template <class T, class T1> inline PRODMM<CT,T1,T> operator*(
    const PRODXM1<CT,T1>& pxm1, const PRODXM2<T,T>& pxm2) 
{ return PRODMM<CT,T1,T>(pxm1.GetX()*pxm2.GetX(),pxm1.GetM(),pxm2.GetM()); }

template <class T, class T2> inline PRODMM<CT,T,T2> operator*(
    const PRODXM1<T,T>& pxm1, const PRODXM2<CT,T2>& pxm2) 
{ return PRODMM<CT,T,T2>(pxm1.GetX()*pxm2.GetX(),pxm1.GetM(),pxm2.GetM()); }

