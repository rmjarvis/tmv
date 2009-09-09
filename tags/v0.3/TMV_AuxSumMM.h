// Things that need to be #defined on entry:
// (The values for a normal Matrix+Matrix are given)
//
// SUMMM	SumMM
// GENMATRIX1	GenMatrix
// GENMATRIX2	GenMatrix
// PRODXM1	ProdXM
// PRODXM2	ProdXM

// m+m
template <class T> inline SUMMM<T,T,T> operator+(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{ return SUMMM<T,T,T>(T(1),m1,T(1),m2); }

template <class T> inline SUMMM<CT,CT,T> operator+(
    const GENMATRIX1<CT>& m1,const GENMATRIX2<T>& m2)
{ return SUMMM<CT,CT,T>(CT(1),m1,CT(1),m2); }

template <class T> inline SUMMM<CT,T,CT> operator+(
    const GENMATRIX1<T>& m1,const GENMATRIX2<CT>& m2)
{ return SUMMM<CT,T,CT>(CT(1),m1,CT(1),m2); }

// m-m
template <class T> inline SUMMM<T,T,T> operator-(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{ return SUMMM<T,T,T>(T(1),m1,T(-1),m2); }

template <class T> inline SUMMM<CT,CT,T> operator-(
    const GENMATRIX1<CT>& m1,const GENMATRIX2<T>& m2)
{ return SUMMM<CT,CT,T>(CT(1),m1,CT(-1),m2); }

template <class T> inline SUMMM<CT,T,CT> operator-(
    const GENMATRIX1<T>& m1,const GENMATRIX2<CT>& m2)
{ return SUMMM<CT,T,CT>(CT(1),m1,CT(-1),m2); }

// (x*m)+m
template <class T, class T1> inline SUMMM<T,T1,T> operator+(
    const PRODXM1<T,T1>& pxm1, const GENMATRIX2<T>& m2)
{ return SUMMM<T,T1,T>(pxm1.GetX1(),pxm1.GetM2(),T(1),m2); }

template <class T> inline SUMMM<CT,T,CT> operator+(
    const PRODXM1<T,T>& pxm1, const GENMATRIX2<CT>& m2)
{ return SUMMM<CT,T,CT>(pxm1.GetX1(),pxm1.GetM2(),CT(1),m2); }

template <class T, class T1> inline SUMMM<CT,T1,T> operator+(
    const PRODXM1<CT,T1>& pxm1, const GENMATRIX2<T>& m2)
{ return SUMMM<CT,T1,T>(pxm1.GetX1(),pxm1.GetM2(),CT(1),m2); }

// m+(x*m)
template <class T, class T2> inline SUMMM<T,T,T2> operator+(
    const GENMATRIX1<T>& m1, const PRODXM2<T,T2>& pxm2)
{ return SUMMM<T,T,T2>(T(1),m1,pxm2.GetX1(),pxm2.GetM2()); }

template <class T> inline SUMMM<CT,CT,T> operator+(
    const GENMATRIX1<CT>& m1, const PRODXM2<T,T>& pxm2)
{ return SUMMM<CT,CT,T>(CT(1),m1,pxm2.GetX1(),pxm2.GetM2()); }

template <class T, class T2> inline SUMMM<CT,T,T2> operator+(
    const GENMATRIX1<T>& m1, const PRODXM2<CT,T2>& pxm2)
{ return SUMMM<CT,T,T2>(CT(1),m1,pxm2.GetX1(),pxm2.GetM2()); }

// (x*m)-m
template <class T, class T1> inline SUMMM<T,T1,T> operator-(
    const PRODXM1<T,T1>& pxm1, const GENMATRIX2<T>& m2)
{ return SUMMM<T,T1,T>(pxm1.GetX1(),pxm1.GetM2(),T(-1),m2); }

template <class T> inline SUMMM<CT,T,CT> operator-(
    const PRODXM1<T,T>& pxm1, const GENMATRIX2<CT>& m2)
{ return SUMMM<CT,T,CT>(pxm1.GetX1(),pxm1.GetM2(),CT(-1),m2); }

template <class T, class T1> inline SUMMM<CT,T1,T> operator-(
    const PRODXM1<CT,T1>& pxm1, const GENMATRIX2<T>& m2)
{ return SUMMM<CT,T1,T>(pxm1.GetX1(),pxm1.GetM2(),CT(-1),m2); }

// m-(x*m)
template <class T, class T2> inline SUMMM<T,T,T2> operator-(
    const GENMATRIX1<T>& m1, const PRODXM2<T,T2>& pxm2)
{ return SUMMM<T,T,T2>(T(1),m1,-pxm2.GetX1(),pxm2.GetM2()); }

template <class T> inline SUMMM<CT,CT,T> operator-(
    const GENMATRIX1<CT>& m1, const PRODXM2<T,T>& pxm2)
{ return SUMMM<CT,CT,T>(CT(1),m1,-pxm2.GetX1(),pxm2.GetM2()); }

template <class T, class T2> inline SUMMM<CT,T,T2> operator-(
    const GENMATRIX1<T>& m1, const PRODXM2<CT,T2>& pxm2)
{ return SUMMM<CT,T,T2>(CT(1),m1,-pxm2.GetX1(),pxm2.GetM2()); }

// (x*m)+(x*m)
template <class T, class T1, class T2> inline SUMMM<T,T1,T2> operator+(
    const PRODXM1<T,T1>& pxm1, const PRODXM2<T,T2>& pxm2)
{ return SUMMM<T,T1,T2>(pxm1.GetX1(),pxm1.GetM2(),pxm2.GetX1(),pxm2.GetM2()); }

template <class T, class T2> inline SUMMM<CT,T,T2> operator+(
    const PRODXM1<T,T>& pxm1, const PRODXM2<CT,T2>& pxm2)
{ return SUMMM<CT,T,T2>(pxm1.GetX1(),pxm1.GetM2(),pxm2.GetX1(),pxm2.GetM2()); }

template <class T, class T1> inline SUMMM<CT,T1,T> operator+(
    const PRODXM1<CT,T1>& pxm1, const PRODXM2<T,T>& pxm2)
{ return SUMMM<CT,T1,T>(pxm1.GetX1(),pxm1.GetM2(),pxm2.GetX1(),pxm2.GetM2()); }

// (x*m)-(x*m)
template <class T, class T1, class T2> inline SUMMM<T,T1,T2> operator-(
    const PRODXM1<T,T1>& pxm1, const PRODXM2<T,T2>& pxm2)
{ 
  return SUMMM<T,T1,T2>(pxm1.GetX1(),pxm1.GetM2(),-pxm2.GetX1(),pxm2.GetM2()); 
}

template <class T, class T2> inline SUMMM<CT,T,T2> operator-(
    const PRODXM1<T,T>& pxm1, const PRODXM2<CT,T2>& pxm2)
{ 
  return SUMMM<CT,T,T2>(pxm1.GetX1(),pxm1.GetM2(),-pxm2.GetX1(),pxm2.GetM2()); 
}

template <class T, class T1> inline SUMMM<CT,T1,T> operator-(
    const PRODXM1<CT,T1>& pxm1, const PRODXM2<T,T>& pxm2)
{ 
  return SUMMM<CT,T1,T>(pxm1.GetX1(),pxm1.GetM2(),-pxm2.GetX1(),pxm2.GetM2()); 
}


