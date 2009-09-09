// Things that need to be #defined on entry:
// (The values for a normal Vector+Vector are given)
//
// SUMVV	SumVV
// GENVECTOR1	GenVector
// GENVECTOR2	GenVector
// PRODXV1	ProdXV
// PRODXV2	ProdXV

// v+v
template <class T> inline SUMVV<T,T,T> operator+(
    const GENVECTOR1<T>& v1, const GENVECTOR2<T>& v2)
{ return SUMVV<T,T,T>(T(1),v1,T(1),v2); }

template <class T> inline SUMVV<CT,CT,T> operator+(
    const GENVECTOR1<CT>& v1,const GENVECTOR2<T>& v2)
{ return SUMVV<CT,CT,T>(CT(1),v1,CT(1),v2); }

template <class T> inline SUMVV<CT,T,CT> operator+(
    const GENVECTOR1<T>& v1,const GENVECTOR2<CT>& v2)
{ return SUMVV<CT,T,CT>(CT(1),v1,CT(1),v2); }

// v-v
template <class T> inline SUMVV<T,T,T> operator-(
    const GENVECTOR1<T>& v1, const GENVECTOR2<T>& v2)
{ return SUMVV<T,T,T>(T(1),v1,T(-1),v2); }

template <class T> inline SUMVV<CT,CT,T> operator-(
    const GENVECTOR1<CT>& v1,const GENVECTOR2<T>& v2)
{ return SUMVV<CT,CT,T>(CT(1),v1,CT(-1),v2); }

template <class T> inline SUMVV<CT,T,CT> operator-(
    const GENVECTOR1<T>& v1,const GENVECTOR2<CT>& v2)
{ return SUMVV<CT,T,CT>(CT(1),v1,CT(-1),v2); }

// (x*v)+v
template <class T> inline SUMVV<T,T,T> operator+(
    const PRODXV1<T,T>& pxv1, const GENVECTOR2<T>& v2)
{ return SUMVV<T,T,T>(pxv1.GetX1(),pxv1.GetV2(),T(1),v2); }

template <class T> inline SUMVV<CT,T,CT> operator+(
    const PRODXV1<T,T>& pxv1, const GENVECTOR2<CT>& v2)
{ return SUMVV<CT,T,CT>(pxv1.GetX1(),pxv1.GetV2(),CT(1),v2); }

template <class T, class T1> inline SUMVV<CT,T1,T> operator+(
    const PRODXV1<CT,T1>& pxv1, const GENVECTOR2<T>& v2)
{ return SUMVV<CT,T1,T>(pxv1.GetX1(),pxv1.GetV2(),CT(1),v2); }

// v+(x*v)
template <class T, class T2> inline SUMVV<T,T,T2> operator+(
    const GENVECTOR1<T>& v1, const PRODXV2<T,T2>& pxv2)
{ return SUMVV<T,T,T2>(T(1),v1,pxv2.GetX1(),pxv2.GetV2()); }

template <class T> inline SUMVV<CT,CT,T> operator+(
    const GENVECTOR1<CT>& v1, const PRODXV2<T,T>& pxv2)
{ return SUMVV<CT,CT,T>(CT(1),v1,pxv2.GetX1(),pxv2.GetV2()); }

template <class T, class T2> inline SUMVV<CT,T,T2> operator+(
    const GENVECTOR1<T>& v1, const PRODXV2<CT,T2>& pxv2)
{ return SUMVV<CT,T,T2>(CT(1),v1,pxv2.GetX1(),pxv2.GetV2()); }

// (x*v)-v
template <class T, class T1> inline SUMVV<T,T1,T> operator-(
    const PRODXV1<T,T1>& pxv1, const GENVECTOR2<T>& v2)
{ return SUMVV<T,T1,T>(pxv1.GetX1(),pxv1.GetV2(),T(-1),v2); }

template <class T> inline SUMVV<CT,T,CT> operator-(
    const PRODXV1<T,T>& pxv1, const GENVECTOR2<CT>& v2)
{ return SUMVV<CT,T,CT>(pxv1.GetX1(),pxv1.GetV2(),CT(-1),v2); }

template <class T, class T1> 
inline SUMVV<CT,T1,T> operator-(
    const PRODXV1<CT,T1>& pxv1, const GENVECTOR2<T>& v2)
{ return SUMVV<CT,T1,T>(pxv1.GetX1(),pxv1.GetV2(),CT(-1),v2); }

// v-(x*v)
template <class T, class T2> inline SUMVV<T,T,T2> operator-(
    const GENVECTOR1<T>& v1, const PRODXV2<T,T2>& pxv2)
{ return SUMVV<T,T,T2>(T(1),v1,-pxv2.GetX1(),pxv2.GetV2()); }

template <class T> inline SUMVV<CT,CT,T> operator-(
    const GENVECTOR1<CT>& v1, const PRODXV2<T,T>& pxv2)
{ return SUMVV<CT,CT,T>(CT(1),v1,-pxv2.GetX1(),pxv2.GetV2()); }

template <class T, class T2> inline SUMVV<CT,T,T2> operator-(
    const GENVECTOR1<T>& v1, const PRODXV2<CT,T2>& pxv2)
{ return SUMVV<CT,T,T2>(CT(1),v1,-pxv2.GetX1(),pxv2.GetV2()); }

// (x*v)+(x*v)
template <class T, class T1, class T2> 
inline SUMVV<T,T1,T2> operator+(
    const PRODXV1<T,T1>& pxv1, const PRODXV2<T,T2>& pxv2)
{ 
  return SUMVV<T,T1,T2>(
      pxv1.GetX1(),pxv1.GetV2(),pxv2.GetX1(),pxv2.GetV2()); 
}

template <class T, class T2> inline SUMVV<CT,T,T2> operator+(
    const PRODXV1<T,T>& pxv1, const PRODXV2<CT,T2>& pxv2)
{ 
  return SUMVV<CT,T,T2>(
      pxv1.GetX1(),pxv1.GetV2(),pxv2.GetX1(),pxv2.GetV2()); 
}

template <class T, class T1> inline SUMVV<CT,T1,T> operator+(
    const PRODXV1<CT,T1>& pxv1, const PRODXV2<T,T>& pxv2)
{ 
  return SUMVV<CT,T1,T>(
      pxv1.GetX1(),pxv1.GetV2(),pxv2.GetX1(),pxv2.GetV2()); 
}

// (x*v)-(x*v)
template <class T, class T1, class T2> 
inline SUMVV<T,T1,T2> operator-(
    const PRODXV1<T,T1>& pxv1, const PRODXV2<T,T2>& pxv2)
{ 
  return SUMVV<T,T1,T2>(
      pxv1.GetX1(),pxv1.GetV2(),-pxv2.GetX1(),pxv2.GetV2()); 
}

template <class T, class T2> inline SUMVV<CT,T,T2> operator-(
    const PRODXV1<T,T>& pxv1, const PRODXV2<CT,T2>& pxv2)
{ 
  return SUMVV<CT,T,T2>(
      pxv1.GetX1(),pxv1.GetV2(),-pxv2.GetX1(),pxv2.GetV2()); 
}

template <class T, class T1> inline SUMVV<CT,T1,T> operator-(
    const PRODXV1<CT,T1>& pxv1, const PRODXV2<T,T>& pxv2)
{ 
  return SUMVV<CT,T1,T>(
      pxv1.GetX1(),pxv1.GetV2(),-pxv2.GetX1(),pxv2.GetV2()); 
}

