
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
template <class T> inline PRODMM<T,T,T,false> operator*(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{ return PRODMM<T,T,T,false>(T(1),m1,m2); }

template <class T> inline PRODMM<CT,CT,T,false> operator*(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& m2)
{ return PRODMM<CT,CT,T,false>(CT(1),m1,m2); }

template <class T> inline PRODMM<CT,T,CT,false> operator*(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& m2)
{ return PRODMM<CT,T,CT,false>(CT(1),m1,m2); }

template <class T> inline PRODMM<T,T,T,true> operator^(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{ return PRODMM<T,T,T,true>(T(1),m1,m2); }

template <class T> inline PRODMM<CT,CT,T,true> operator^(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& m2)
{ return PRODMM<CT,CT,T,true>(CT(1),m1,m2); }

template <class T> inline PRODMM<CT,T,CT,true> operator^(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& m2)
{ return PRODMM<CT,T,CT,true>(CT(1),m1,m2); }

// x*(x*m*m)
template <class T, class T1, class T2, bool trans> 
inline PRODMM<T,T1,T2,trans> operator*(
    T x1, const PRODMM<T,T1,T2,trans>& pmm2)
{ return PRODMM<T,T1,T2,trans>(x1*pmm2.GetX1(),pmm2.GetM2(),pmm2.GetM3()); }

template <class T, bool trans> inline PRODMM<CT,T,T,trans> operator*(
    CT x1, const PRODMM<T,T,T,trans>& pmm2)
{ return PRODMM<CT,T,T,trans>(x1*pmm2.GetX1(),pmm2.GetM2(),pmm2.GetM3()); }

template <class T, bool trans> inline PRODMM<CT,T,T,trans> operator*(
    CCT x1, const PRODMM<T,T,T,trans>& pmm2)
{ return PRODMM<CT,T,T,trans>(CT(x1)*pmm2.GetX1(),pmm2.GetM2(),pmm2.GetM3()); }

template <class T, bool trans> inline PRODMM<CT,T,T,trans> operator*(
    VCT x1, const PRODMM<T,T,T,trans>& pmm2)
{ return PRODMM<CT,T,T,trans>(CT(x1)*pmm2.GetX1(),pmm2.GetM2(),pmm2.GetM3()); }

template <class T, class T1, class T2, bool trans> 
inline PRODMM<CT,T1,T2,trans> operator*(
    T x1, const PRODMM<CT,T1,T2,trans>& pmm2)
{ return PRODMM<CT,T1,T2,trans>(x1*pmm2.GetX1(),pmm2.GetM2(),pmm2.GetM3()); }

template <class T, class T1, class T2, bool trans> 
inline PRODMM<CT,T1,T2,trans> operator*(
    CCT x1, const PRODMM<CT,T1,T2,trans>& pmm2)
{ 
  return PRODMM<CT,T1,T2,trans>(CT(x1)*pmm2.GetX1(),pmm2.GetM2(),pmm2.GetM3());
}

template <class T, class T1, class T2, bool trans> 
inline PRODMM<CT,T1,T2,trans> operator*(
    VCT x1, const PRODMM<CT,T1,T2,trans>& pmm2)
{
  return PRODMM<CT,T1,T2,trans>(CT(x1)*pmm2.GetX1(),pmm2.GetM2(),pmm2.GetM3());
}

// (x*m*m)*x
template <class T, class T1, class T2, bool trans> 
inline PRODMM<T,T1,T2,trans> operator*(
    const PRODMM<T,T1,T2,trans>& pmv1, T x2)
{ return PRODMM<T,T1,T2,trans>(x2*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3()); }

template <class T, bool trans> inline PRODMM<CT,T,T,trans> operator*(
    const PRODMM<T,T,T,trans>& pmv1, CT x2)
{ return PRODMM<CT,T,T,trans>(x2*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3()); }

template <class T, bool trans> inline PRODMM<CT,T,T,trans> operator*(
    const PRODMM<T,T,T,trans>& pmv1, CCT x2)
{ return PRODMM<CT,T,T,trans>(CT(x2)*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3()); }

template <class T, bool trans> inline PRODMM<CT,T,T,trans> operator*(
    const PRODMM<T,T,T,trans>& pmv1, VCT x2)
{ return PRODMM<CT,T,T,trans>(CT(x2)*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3()); }

template <class T, class T1, class T2, bool trans> 
inline PRODMM<CT,T1,T2,trans> operator*(
    const PRODMM<CT,T1,T2,trans>& pmv1, T x2)
{ return PRODMM<CT,T1,T2,trans>(x2*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3()); }

template <class T, class T1, class T2, bool trans> 
inline PRODMM<CT,T1,T2,trans> operator*(
    const PRODMM<CT,T1,T2,trans>& pmv1, CCT x2)
{
  return PRODMM<CT,T1,T2,trans>(CT(x2)*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3());
}

template <class T, class T1, class T2, bool trans> 
inline PRODMM<CT,T1,T2,trans> operator*(
    const PRODMM<CT,T1,T2,trans>& pmv1, VCT x2)
{
  return PRODMM<CT,T1,T2,trans>(CT(x2)*pmv1.GetX1(),pmv1.GetM2(),pmv1.GetM3());
}

// m*(x*m)
template <class T, class T2> inline PRODMM<T,T,T2,false> operator*(
    const GENMATRIX1<T>& m1, const PRODXM2<T,T2>& pxm2)
{ return PRODMM<T,T,T,false>(pxm2.GetX1(),m1,pxm2.GetM2()); }

template <class T> inline PRODMM<CT,CT,T,false> operator*(
    const GENMATRIX1<CT>& m1, const PRODXM2<T,T>& pxm2)
{ return PRODMM<CT,CT,T,false>(pxm2.GetX1(),m1,pxm2.GetM2()); }

template <class T, class T2> inline PRODMM<CT,T,T2,false> operator*(
    const GENMATRIX1<T>& m1, const PRODXM2<CT,T2>& pxm2)
{ return PRODMM<CT,T,T2,false>(pxm2.GetX1(),m1,pxm2.GetM2()); }

template <class T, class T2> inline PRODMM<T,T,T2,true> operator^(
    const GENMATRIX1<T>& m1, const PRODXM2<T,T2>& pxm2)
{ return PRODMM<T,T,T,true>(pxm2.GetX1(),m1,pxm2.GetM2()); }

template <class T> inline PRODMM<CT,CT,T,true> operator^(
    const GENMATRIX1<CT>& m1, const PRODXM2<T,T>& pxm2)
{ return PRODMM<CT,CT,T,true>(pxm2.GetX1(),m1,pxm2.GetM2()); }

template <class T, class T2> inline PRODMM<CT,T,T2,true> operator^(
    const GENMATRIX1<T>& m1, const PRODXM2<CT,T2>& pxm2)
{ return PRODMM<CT,T,T2,true>(pxm2.GetX1(),m1,pxm2.GetM2()); }

// (x*m)*m
template <class T, class T1> inline PRODMM<T,T1,T,false> operator*(
    const PRODXM1<T,T1>& pxm1, const GENMATRIX2<T>& m2)
{ return PRODMM<T,T1,T,false>(pxm1.GetX1(),pxm1.GetM2(),m2); }

template <class T> inline PRODMM<CT,T,CT,false> operator*(
    const PRODXM1<T,T>& pxm1, const GENMATRIX2<CT>& m2)
{ return PRODMM<CT,T,CT,false>(pxm1.GetX1(),pxm1.GetM2(),m2); }

template <class T, class T1> inline PRODMM<CT,T1,T,false> operator*(
    const PRODXM1<CT,T1>& pxm1, const GENMATRIX2<T>& m2)
{ return PRODMM<CT,T1,T,false>(pxm1.GetX1(),pxm1.GetM2(),m2); }

template <class T, class T1> inline PRODMM<T,T1,T,true> operator^(
    const PRODXM1<T,T1>& pxm1, const GENMATRIX2<T>& m2)
{ return PRODMM<T,T1,T,true>(pxm1.GetX1(),pxm1.GetM2(),m2); }

template <class T> inline PRODMM<CT,T,CT,true> operator^(
    const PRODXM1<T,T>& pxm1, const GENMATRIX2<CT>& m2)
{ return PRODMM<CT,T,CT,true>(pxm1.GetX1(),pxm1.GetM2(),m2); }

template <class T, class T1> inline PRODMM<CT,T1,T,true> operator^(
    const PRODXM1<CT,T1>& pxm1, const GENMATRIX2<T>& m2)
{ return PRODMM<CT,T1,T,true>(pxm1.GetX1(),pxm1.GetM2(),m2); }

// (x*m)*(x*m)
template <class T, class T1, class T2> inline PRODMM<T,T1,T2,false> operator*(
    const PRODXM1<T,T1>& pxm1, const PRODXM2<T,T2>& pxm2) 
{ 
  return PRODMM<T,T1,T2,false>(pxm1.GetX1()*pxm2.GetX1(),pxm1.GetM2(),
    pxm2.GetM2()); 
}

template <class T, class T1> inline PRODMM<CT,T1,T,false> operator*(
    const PRODXM1<CT,T1>& pxm1, const PRODXM2<T,T>& pxm2) 
{
  return PRODMM<CT,T1,T,false>(pxm1.GetX1()*pxm2.GetX1(),pxm1.GetM2(),
      pxm2.GetM2()); 
}

template <class T, class T2> inline PRODMM<CT,T,T2,false> operator*(
    const PRODXM1<T,T>& pxm1, const PRODXM2<CT,T2>& pxm2) 
{
  return PRODMM<CT,T,T2,false>(pxm1.GetX1()*pxm2.GetX1(),pxm1.GetM2(),
      pxm2.GetM2()); 
}

template <class T, class T1, class T2> inline PRODMM<T,T1,T2,true> operator^(
    const PRODXM1<T,T1>& pxm1, const PRODXM2<T,T2>& pxm2) 
{ 
  return PRODMM<T,T1,T2,true>(pxm1.GetX1()*pxm2.GetX1(),pxm1.GetM2(),
    pxm2.GetM2()); 
}

template <class T, class T1> inline PRODMM<CT,T1,T,true> operator^(
    const PRODXM1<CT,T1>& pxm1, const PRODXM2<T,T>& pxm2) 
{
  return PRODMM<CT,T1,T,true>(pxm1.GetX1()*pxm2.GetX1(),pxm1.GetM2(),
      pxm2.GetM2()); 
}

template <class T, class T2> inline PRODMM<CT,T,T2,true> operator^(
    const PRODXM1<T,T>& pxm1, const PRODXM2<CT,T2>& pxm2) 
{
  return PRODMM<CT,T,T2,true>(pxm1.GetX1()*pxm2.GetX1(),pxm1.GetM2(),
      pxm2.GetM2()); 
}

