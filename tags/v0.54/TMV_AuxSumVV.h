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

// -(x*v+x*v)
template <class T, class T1, class T2> inline SUMVV<T,T1,T2> operator-(
    const SUMVV<T,T1,T2>& svv)
{ return SUMVV<T,T1,T2>(-svv.GetX1(),svv.GetV1(),-svv.GetX2,svv.GetV2); }

// x*(x*v+x*v)
template <class T, class T1, class T2> inline SUMVV<T,T1,T2> operator*(
    const T x, const SUMVV<T,T1,T2>& svv)
{ return SUMVV<T,T1,T2>(svv.GetX1()*x,svv.GetV1(),svv.GetX2()*x,svv.GetV2()); }

template <class T, class T1, class T2> inline SUMVV<CT,T1,T2> operator*(
    const T x, const SUMVV<CT,T1,T2>& svv)
{ return SUMVV<CT,T1,T2>(svv.GetX1()*x,svv.GetV1(),svv.GetX2()*x,svv.GetV2()); }

template <class T> inline SUMVV<CT,T,T> operator*(
    const CT x, const SUMVV<T,T,T>& svv)
{ return SUMVV<CT,T,T>(svv.GetX1()*x,svv.GetV1(),svv.GetX2()*x,svv.GetV2()); }

template <class T> inline SUMVV<CT,T,T> operator*(
    const CCT x, const SUMVV<T,T,T>& svv)
{ 
  return SUMVV<CT,T,T>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

template <class T> inline SUMVV<CT,T,T> operator*(
    const VCT x, const SUMVV<T,T,T>& svv)
{ 
  return SUMVV<CT,T,T>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

template <class T, class T1, class T2> inline SUMVV<CT,T1,T2> operator*(
    const CCT x, const SUMVV<CT,T1,T2>& svv)
{ 
  return SUMVV<CT,T1,T2>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

template <class T, class T1, class T2> inline SUMVV<CT,T1,T2> operator*(
    const VCT x, const SUMVV<CT,T1,T2>& svv)
{ 
  return SUMVV<CT,T1,T2>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

// (x*v+x*v)*x
template <class T, class T1, class T2> inline SUMVV<T,T1,T2> operator*(
    const SUMVV<T,T1,T2>& svv, const T x)
{ return SUMVV<T,T1,T2>(svv.GetX1()*x,svv.GetV1(),svv.GetX2()*x,svv.GetV2()); }

template <class T, class T1, class T2> inline SUMVV<CT,T1,T2> operator*(
    const SUMVV<CT,T1,T2>& svv, const T x)
{ return SUMVV<CT,T1,T2>(svv.GetX1()*x,svv.GetV1(),svv.GetX2()*x,svv.GetV2()); }

template <class T> inline SUMVV<CT,T,T> operator*(
    const SUMVV<T,T,T>& svv, const CT x)
{ return SUMVV<CT,T,T>(svv.GetX1()*x,svv.GetV1(),svv.GetX2()*x,svv.GetV2()); }

template <class T> inline SUMVV<CT,T,T> operator*(
    const SUMVV<T,T,T>& svv, const CCT x)
{ 
  return SUMVV<CT,T,T>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

template <class T> inline SUMVV<CT,T,T> operator*(
    const SUMVV<T,T,T>& svv, const VCT x)
{ 
  return SUMVV<CT,T,T>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

template <class T, class T1, class T2> inline SUMVV<CT,T1,T2> operator*(
    const SUMVV<CT,T1,T2>& svv, const CCT x)
{ 
  return SUMVV<CT,T1,T2>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

template <class T, class T1, class T2> inline SUMVV<CT,T1,T2> operator*(
    const SUMVV<CT,T1,T2>& svv, const VCT x)
{ 
  return SUMVV<CT,T1,T2>(svv.GetX1()*CT(x),svv.GetV1(),
      svv.GetX2()*CT(x),svv.GetV2()); 
}

// (x*v+x*v)/x
template <class T, class T1, class T2> inline SUMVV<T,T1,T2> operator/(
    const SUMVV<T,T1,T2>& svv, const T x)
{ return SUMVV<T,T1,T2>(svv.GetX1()/x,svv.GetV1(),svv.GetX2()/x,svv.GetV2()); }

template <class T, class T1, class T2> inline SUMVV<CT,T1,T2> operator/(
    const SUMVV<CT,T1,T2>& svv, const T x)
{ return SUMVV<CT,T1,T2>(svv.GetX1()/x,svv.GetV1(),svv.GetX2()/x,svv.GetV2()); }

template <class T> inline SUMVV<CT,T,T> operator/(
    const SUMVV<T,T,T>& svv, const CT x)
{ return SUMVV<CT,T,T>(svv.GetX1()/x,svv.GetV1(),svv.GetX2()/x,svv.GetV2()); }

template <class T> inline SUMVV<CT,T,T> operator/(
    const SUMVV<T,T,T>& svv, const CCT x)
{ 
  return SUMVV<CT,T,T>(svv.GetX1()/CT(x),svv.GetV1(),
      svv.GetX2()/CT(x),svv.GetV2()); 
}

template <class T> inline SUMVV<CT,T,T> operator/(
    const SUMVV<T,T,T>& svv, const VCT x)
{ 
  return SUMVV<CT,T,T>(svv.GetX1()/CT(x),svv.GetV1(),
      svv.GetX2()/CT(x),svv.GetV2()); 
}

template <class T, class T1, class T2> inline SUMVV<CT,T1,T2> operator/(
    const SUMVV<CT,T1,T2>& svv, const CCT x)
{ 
  return SUMVV<CT,T1,T2>(svv.GetX1()/CT(x),svv.GetV1(),
      svv.GetX2()/CT(x),svv.GetV2()); 
}

template <class T, class T1, class T2> inline SUMVV<CT,T1,T2> operator/(
    const SUMVV<CT,T1,T2>& svv, const VCT x)
{ 
  return SUMVV<CT,T1,T2>(svv.GetX1()/CT(x),svv.GetV1(),
      svv.GetX2()/CT(x),svv.GetV2()); 
}

// (x*v)+v
template <class T> inline SUMVV<T,T,T> operator+(
    const PRODXV1<T,T>& pxv, const GENVECTOR2<T>& v)
{ return SUMVV<T,T,T>(pxv.GetX(),pxv.GetV(),T(1),v); }

template <class T> inline SUMVV<CT,T,CT> operator+(
    const PRODXV1<T,T>& pxv, const GENVECTOR2<CT>& v)
{ return SUMVV<CT,T,CT>(pxv.GetX(),pxv.GetV(),CT(1),v); }

template <class T, class T1> inline SUMVV<CT,T1,T> operator+(
    const PRODXV1<CT,T1>& pxv, const GENVECTOR2<T>& v)
{ return SUMVV<CT,T1,T>(pxv.GetX(),pxv.GetV(),CT(1),v); }

// v+(x*v)
template <class T, class T2> inline SUMVV<T,T,T2> operator+(
    const GENVECTOR1<T>& v, const PRODXV2<T,T2>& pxv)
{ return SUMVV<T,T,T2>(T(1),v,pxv.GetX(),pxv.GetV()); }

template <class T> inline SUMVV<CT,CT,T> operator+(
    const GENVECTOR1<CT>& v, const PRODXV2<T,T>& pxv)
{ return SUMVV<CT,CT,T>(CT(1),v,pxv.GetX(),pxv.GetV()); }

template <class T, class T2> inline SUMVV<CT,T,T2> operator+(
    const GENVECTOR1<T>& v, const PRODXV2<CT,T2>& pxv)
{ return SUMVV<CT,T,T2>(CT(1),v,pxv.GetX(),pxv.GetV()); }

// (x*v)-v
template <class T, class T1> inline SUMVV<T,T1,T> operator-(
    const PRODXV1<T,T1>& pxv, const GENVECTOR2<T>& v)
{ return SUMVV<T,T1,T>(pxv.GetX(),pxv.GetV(),T(-1),v); }

template <class T> inline SUMVV<CT,T,CT> operator-(
    const PRODXV1<T,T>& pxv, const GENVECTOR2<CT>& v)
{ return SUMVV<CT,T,CT>(pxv.GetX(),pxv.GetV(),CT(-1),v); }

template <class T, class T1> 
inline SUMVV<CT,T1,T> operator-(
    const PRODXV1<CT,T1>& pxv, const GENVECTOR2<T>& v)
{ return SUMVV<CT,T1,T>(pxv.GetX(),pxv.GetV(),CT(-1),v); }

// v-(x*v)
template <class T, class T2> inline SUMVV<T,T,T2> operator-(
    const GENVECTOR1<T>& v, const PRODXV2<T,T2>& pxv)
{ return SUMVV<T,T,T2>(T(1),v,-pxv.GetX(),pxv.GetV()); }

template <class T> inline SUMVV<CT,CT,T> operator-(
    const GENVECTOR1<CT>& v, const PRODXV2<T,T>& pxv)
{ return SUMVV<CT,CT,T>(CT(1),v,-pxv.GetX(),pxv.GetV()); }

template <class T, class T2> inline SUMVV<CT,T,T2> operator-(
    const GENVECTOR1<T>& v, const PRODXV2<CT,T2>& pxv)
{ return SUMVV<CT,T,T2>(CT(1),v,-pxv.GetX(),pxv.GetV()); }

// (x*v)+(x*v)
template <class T, class T1, class T2> inline SUMVV<T,T1,T2> operator+(
    const PRODXV1<T,T1>& pxv1, const PRODXV2<T,T2>& pxv2)
{ return SUMVV<T,T1,T2>(pxv1.GetX(),pxv1.GetV(),pxv2.GetX(),pxv2.GetV()); }

template <class T, class T2> inline SUMVV<CT,T,T2> operator+(
    const PRODXV1<T,T>& pxv1, const PRODXV2<CT,T2>& pxv2)
{ return SUMVV<CT,T,T2>(pxv1.GetX(),pxv1.GetV(),pxv2.GetX(),pxv2.GetV()); }

template <class T, class T1> inline SUMVV<CT,T1,T> operator+(
    const PRODXV1<CT,T1>& pxv1, const PRODXV2<T,T>& pxv2)
{ return SUMVV<CT,T1,T>(pxv1.GetX(),pxv1.GetV(),pxv2.GetX(),pxv2.GetV()); }

// (x*v)-(x*v)
template <class T, class T1, class T2> inline SUMVV<T,T1,T2> operator-(
    const PRODXV1<T,T1>& pxv1, const PRODXV2<T,T2>& pxv2)
{ return SUMVV<T,T1,T2>(pxv1.GetX(),pxv1.GetV(),-pxv2.GetX(),pxv2.GetV()); }

template <class T, class T2> inline SUMVV<CT,T,T2> operator-(
    const PRODXV1<T,T>& pxv1, const PRODXV2<CT,T2>& pxv2)
{ return SUMVV<CT,T,T2>(pxv1.GetX(),pxv1.GetV(),-pxv2.GetX(),pxv2.GetV()); }

template <class T, class T1> inline SUMVV<CT,T1,T> operator-(
    const PRODXV1<CT,T1>& pxv1, const PRODXV2<T,T>& pxv2)
{ return SUMVV<CT,T1,T>(pxv1.GetX(),pxv1.GetV(),-pxv2.GetX(),pxv2.GetV()); }

