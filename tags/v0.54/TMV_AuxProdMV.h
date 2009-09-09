
// Need to define the following with #define statements.
// (The given definition is for a regular Matrix*Vector.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENVECTOR GenVector
// #define GENMATRIX GenMatrix
// #define PRODMV ProdMV
// #define PRODXV ProdXV
// #define PRODXM ProdXM

// m*v
template <class T> inline PRODMV<T,T,false,T> operator*(
    const GENMATRIX<T>& m, const GENVECTOR<T>& v)
{ return PRODMV<T,T,false,T>(T(1),m,v); }

template <class T> inline PRODMV<CT,CT,false,T> operator*(
    const GENMATRIX<CT>& m, const GENVECTOR<T>& v)
{ return PRODMV<CT,CT,false,T>(CT(1),m,v); }

template <class T> inline PRODMV<CT,T,false,CT> operator*(
    const GENMATRIX<T>& m, const GENVECTOR<CT>& v)
{ return PRODMV<CT,T,false,CT>(CT(1),m,v); }

// v*m
template <class T> inline PRODMV<T,T,true,T> operator*(
    const GENVECTOR<T>& v, const GENMATRIX<T>& m)
{ return PRODMV<T,T,true,T>(T(1),m,v); }

template <class T> inline PRODMV<CT,T,true,CT> operator*(
    const GENVECTOR<CT>& v, const GENMATRIX<T>& m)
{ return PRODMV<CT,T,true,CT>(CT(1),m,v); }

template <class T> inline PRODMV<CT,CT,true,T> operator*(
    const GENVECTOR<T>& v, const GENMATRIX<CT>& m)
{ return PRODMV<CT,CT,true,T>(CT(1),m,v); }

// m*(x*v)
template <class T, class T2> inline PRODMV<T,T,false,T2> operator*(
    const GENMATRIX<T>& m, const PRODXV<T,T2>& pxv)
{ return PRODMV<T,T,false,T2>(pxv.GetX(),m,pxv.GetV()); }

template <class T> inline PRODMV<CT,CT,false,T> operator*(
    const GENMATRIX<CT>& m, const PRODXV<T,T>& pxv)
{ return PRODMV<CT,CT,false,T>(pxv.GetX(),m,pxv.GetV()); }

template <class T, class T2> inline PRODMV<CT,T,false,T2> operator*(
    const GENMATRIX<T>& m, const PRODXV<CT,T2>& pxv)
{ return PRODMV<CT,T,false,T2>(pxv.GetX(),m,pxv.GetV()); }

// (x*v)*m
template <class T, class T2> inline PRODMV<T,T,true,T2> operator*(
    const PRODXV<T,T2>& pxv, const GENMATRIX<T>& m)
{ return PRODMV<T,T,true,T2>(pxv.GetX(),m,pxv.GetV()); }

template <class T> inline PRODMV<CT,CT,true,T> operator*(
    const PRODXV<T,T>& pxv, const GENMATRIX<CT>& m)
{ return PRODMV<CT,CT,true,T>(pxv.GetX(),m,pxv.GetV()); }

template <class T, class T2> inline PRODMV<CT,T,true,T2> operator*(
    const PRODXV<CT,T2>& pxv, const GENMATRIX<T>& m)
{ return PRODMV<CT,T,true,T2>(pxv.GetX(),m,pxv.GetV()); }

// (x*m)*v
template <class T, class T1> inline PRODMV<T,T1,false,T> operator*(
    const PRODXM<T,T1>& pxm, const GENVECTOR<T>& v)
{ return PRODMV<T,T1,false,T>(pxm.GetX(),pxm.GetM(),v); }

template <class T> inline PRODMV<CT,T,false,CT> operator*(
    const PRODXM<T,T>& pxm, const GENVECTOR<CT>& v)
{ return PRODMV<CT,T,false,CT>(pxm.GetX(),pxm.GetM(),v); }

template <class T, class T1> inline PRODMV<CT,T1,false,T> operator*(
    const PRODXM<CT,T1>& pxm, const GENVECTOR<T>& v)
{ return PRODMV<CT,T1,false,T>(pxm.GetX(),pxm.GetM(),v); }

// v*(x*m)
template <class T, class T1> inline PRODMV<T,T1,true,T> operator*(
    const GENVECTOR<T>& v, const PRODXM<T,T1>& pxm)
{ return PRODMV<T,T1,true,T>(pxm.GetX(),pxm.GetM(),v); }

template <class T> inline PRODMV<CT,T,true,CT> operator*(
    const GENVECTOR<CT>& v, const PRODXM<T,T>& pxm)
{ return PRODMV<CT,T,true,CT>(pxm.GetX(),pxm.GetM(),v); }

template <class T, class T1> inline PRODMV<CT,T1,true,T> operator*(
    const GENVECTOR<T>& v, const PRODXM<CT,T1>& pxm)
{ return PRODMV<CT,T1,true,T>(pxm.GetX(),pxm.GetM(),v); }

// (x*m)*(x*v)
template <class T, class T1, class T2> inline PRODMV<T,T1,false,T2> operator*(
    const PRODXM<T,T1>& pxm, const PRODXV<T,T2>& pxv) 
{ 
  return PRODMV<T,T1,false,T2>(pxm.GetX()*pxv.GetX(),pxm.GetM(),pxv.GetV()); 
}

template <class T, class T1> inline PRODMV<CT,T1,false,T> operator*(
    const PRODXM<CT,T1>& pxm, const PRODXV<T,T>& pxv) 
{ 
  return PRODMV<CT,T1,false,T>(pxm.GetX()*pxv.GetX(),pxm.GetM(),pxv.GetV()); 
}

template <class T, class T2> inline PRODMV<CT,T,false,T2> operator*(
    const PRODXM<T,T>& pxm, const PRODXV<CT,T2>& pxv) 
{ 
  return PRODMV<CT,T,false,T2>(pxm.GetX()*pxv.GetX(),pxm.GetM(),pxv.GetV()); 
}

// (x*v)*(x*m)
template <class T, class T1, class T2> inline PRODMV<T,T1,true,T2> operator*(
    const PRODXV<T,T2>& pxv, const PRODXM<T,T1>& pxm)
{ 
  return PRODMV<T,T1,true,T2>(pxv.GetX()*pxm.GetX(),pxm.GetM(),pxv.GetV()); 
}

template <class T, class T2> inline PRODMV<CT,T,true,T2> operator*(
    const PRODXV<CT,T2>& pxv, const PRODXM<T,T>& pxm)
{ 
  return PRODMV<CT,T,true,T2>(pxv.GetX()*pxm.GetX(),pxm.GetM(),pxv.GetV()); 
}

template <class T, class T1> inline PRODMV<CT,T1,true,T> operator*(
    const PRODXV<T,T>& pxv, const PRODXM<CT,T1>& pxm)
{ 
  return PRODMV<CT,T1,true,T>(pxv.GetX()*pxm.GetX(),pxm.GetM(),pxv.GetV()); 
}
