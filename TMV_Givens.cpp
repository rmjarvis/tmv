
#include "TMV.h"
#include "TMV_Givens.h"

namespace tmv {

  // 
  // Givens_Rotate
  //
  template <class T> Givens<T> Givens_Rotate(T& x, T& y)
  {
    // Use a Givens matrix G to rotate the vector so y = 0:
    // G [ x ] = [ r ]
    //   [ y ]   [ 0 ]
    // Also, return the Givens matrix used to do so.
    
    const RealType(T) one(1);
    const RealType(T) zero(0);

    RealType(T) maxabs_x = MAXABS(x);
    RealType(T) maxabs_y = MAXABS(y);
    if (maxabs_y == zero) {
      return Givens<T>(one,T(0));
    } else if (maxabs_x == zero) {
      RealType(T) absy = abs(y);
      T s = SIGN(CONJ(y),absy);
      x = absy; y = zero;
      return Givens<T>(zero,s);
    } else {
      // This test is better than abs(x) > abs(y) for complex,
      // since it saves 2 sqrt's and has roughly the same effect
      // of preventing overflow and reducing rounding errors.
      if (maxabs_x > maxabs_y) {
	if (maxabs_y <= SqrtEpsilon<T>()*maxabs_x) {
	  // Then everything simplifies:
	  // c = 1
	  // s = (y/x)* = xy*/|x|^2
	  // r = f
	  RealType(T) invnormx = one/NORM(x);
	  T s = x*CONJ(y)*invnormx;
	  y = zero;
	  return Givens<T>(one,s);
	} else {
	  // c = 1/sqrt(1+|y/x|^2)
	  // s = (y/x)*/sqrt(1+|y/x|^2)
	  // r = x sqrt(1+|y/x|^2)
	  // We calculate c-1 rather than c, since it will be first order in
	  // y/x, which is nominally small.  This leads to better accuracy
	  // on the Givens multiplies.
	  // c-1 = (1-sqrt(1+n))/sqrt(1+n)
	  //     = -n/sqrt(1+n)/(1+sqrt(1+n)) 
	  //     = -n/(1+n+sqrt(1+n)
	  // Likewise, we get a slightly more accurate calculation of r
	  // if we calculate r-x and add this to x:
	  // r-x = x (sqrt(1+|y/x|^2)-1) = x |y/x|^2 / (1 + sqrt(1+|y/x|^2))
	  T yoverx = y/x;
	  RealType(T) n = NORM(yoverx);
	  RealType(T) sqrtfactor = SQRT(one+n);
	  RealType(T) invsqrtfactor = one/sqrtfactor;
	  T s = CONJ(yoverx)*invsqrtfactor;
	  x += x*(n/(one+sqrtfactor));
	  y = zero;
	  RealType(T) c = one/sqrtfactor;
	  return Givens<T>(c,s);
	}
      } else {
	// As above, we store c-1 rather than c
	// even though it isn't small here.
	T xovery = x/y;
	RealType(T) n = NORM(xovery);
	RealType(T) absxovery = SQRT(n);
	if (n <= Epsilon<T>()) {
	  // c = |x/y|
	  // s = (x/y)/|x/y|
	  // r = x/|x/y|
	  T s = SIGN(xovery,absxovery);
	  x = s*y;
	  y = zero;
	  return Givens<T>(absxovery,s);
	} else {
	  // c = |x/y|/sqrt(1+|x/y|^2)
	  // s = (x/y)/|x/y|/sqrt(1+|x/y|^2)
	  // r = x/|x/y| sqrt(|x|^2+|y|^2)
	  RealType(T) sqrtfactor = SQRT(one+n);
	  RealType(T) invsqrtfactor = one/sqrtfactor;
	  T signxovery = SIGN(xovery,absxovery); // (x/y)/|x/y|
	  T s = signxovery*invsqrtfactor;
	  x = y * signxovery * sqrtfactor; 
	  y = T(0);
	  RealType(T) c = absxovery*invsqrtfactor;
	  return Givens<T>(c,s);
	}
      }
    }
  }

  //
  // GivensMult Scalars
  //
  template <class Tg, class T> void GivensMult(
      RealType(Tg) c, Tg s, T& x, T& y) 
  {
    // [ x' ] = [  c  s  ] [ x ] = [  cx+sy  ]
    // [ y' ]   [ -s* c* ] [ y ]   [ c*y-s*x ]
    T xx = c*x+s*y;
    y = c*y-tmv::CONJ(s)*x;
    x = xx;
  }
  template <class Tg, class T> void GivensMult(
      RealType(Tg) c, Tg s, ConjRef<T> x, ConjRef<T> y)
  {
    T xx = c*x+s*y;
    y = c*y-tmv::CONJ(s)*x;
    x = xx;
  }
  template <class Tg, class T> void GivensMult(
      RealType(Tg) c, Tg s, VarConjRef<T> x, VarConjRef<T> y)
  {
    T xx = c*x+s*y;
    y = c*y-tmv::CONJ(s)*x;
    x = xx;
  }

  //
  // GivensMult vectors
  //

  template <bool c0, class Tg, class T> inline void NonBlasGivensMult(
      RealType(Tg) c, Tg s,
      const VectorView<T>& v1, const VectorView<T>& v2)
  {
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v1.size() > 0);
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);
    TMVAssert(!v1.SameAs(v2));
    TMVAssert(v2.ct() == NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(c0 == v1.isconj());

    T* v1ptr = v1.ptr();
    T* v2ptr = v2.ptr();
    const int step0 = v1.step();
    const int step1 = v2.step();

    if (step0 == 1 && step1 == 1)
      for(size_t i=v1.size();i>0;--i,++v1ptr,++v2ptr) 
	GivensMult(c,s,*v1ptr,*v2ptr); 
    else
      for(size_t i=v1.size();i>0;--i,v1ptr+=step0,v2ptr+=step1) 
	GivensMult(c,s,*v1ptr,*v2ptr); 
  }

#ifdef BLAS
  template <class Tg, class T> inline void BlasGivensMult(
      RealType(Tg) c, Tg s,
      const VectorView<T>& v1, const VectorView<T>& v2)
  { 
    if (v1.isconj()) NonBlasGivensMult<true>(c,s,v1,v2); 
    else NonBlasGivensMult<false>(c,s,v1,v2); 
  }
  template <> inline void BlasGivensMult(
      double c, double s,
      const VectorView<double>& v1, const VectorView<double>& v2)
  { 
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n=v1.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(drot) (BLASV(n),BLASP(v1.ptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2),BLASV(c),BLASV(s)); 
  }
  template <> inline void BlasGivensMult(
      double c, double s,
      const VectorView<complex<double> >& v1,
      const VectorView<complex<double> >& v2)
  { 
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n=v1.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(zdrot) (BLASV(n),BLASP(v1.ptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2),BLASV(c),BLASV(s)); 
  }
#ifndef NOFLOAT
  template <> inline void BlasGivensMult(
      float c, float s,
      const VectorView<float>& v1, const VectorView<float>& v2)
  { 
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n=v1.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(srot) (BLASV(n),BLASP(v1.ptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2),BLASV(c),BLASV(s)); 
  }
  template <> inline void BlasGivensMult(
      float c, float s,
      const VectorView<complex<float> >& v1,
      const VectorView<complex<float> >& v2)
  { 
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n=v1.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(csrot) (BLASV(n),BLASP(v1.ptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2),BLASV(c),BLASV(s)); 
  }
#endif
#ifdef ELAP
  template <> inline void BlasGivensMult(
      double c, complex<double> s,
      const VectorView<complex<double> >& v1,
      const VectorView<complex<double> >& v2)
  {
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n = v1.size();
    int s0 = v1.step();
    int s1 = v2.step();
    LAPNAMEX(zrot) (LAPV(n),LAPP(v1.ptr()),LAPV(s0),LAPP(v2.ptr()),LAPV(s1),
	LAPV(c),LAPP(&s)); 
  }
#ifndef NOFLOAT
  template <> inline void BlasGivensMult(
      float c, complex<float> s,
      const VectorView<complex<float> >& v1,
      const VectorView<complex<float> >& v2)
  {
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n = v1.size();
    int s0 = v1.step();
    int s1 = v2.step();
    LAPNAMEX(crot) (LAPV(n),LAPP(v1.ptr()),LAPV(s0),LAPP(v2.ptr()),LAPV(s1),
	LAPV(c),LAPP(&s)); 
  }
#endif
#endif // ELAP
#endif // BLAS

  template <class Tg, class T> inline void DoGivensMult(
      RealType(Tg) c, Tg s,
      const VectorView<T>& v1, const VectorView<T>& v2)
  {
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v1.size() > 0);
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);
    TMVAssert(!v1.SameAs(v2));
    TMVAssert(v2.ct() == NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);

#ifdef BLAS
    if (v1.step() > 0 && v2.step() > 0 && v1.ct() == NonConj)
      BlasGivensMult(c,s,v1,v2);
    else
#endif
      if (v1.isconj()) NonBlasGivensMult<true>(c,s,v1,v2); 
      else NonBlasGivensMult<false>(c,s,v1,v2); 
  }

  template <class Tg, class T> void GivensMult(
      RealType(Tg) c, Tg s,
      const VectorView<T>& v1, const VectorView<T>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);
    TMVAssert(!v1.SameAs(v2));

    if (v1.size() > 0 && s != Tg(0)) {
      if (ShouldReverse(v1.step(),v2.step()))
	GivensMult(c,s,v1.Reverse(),v2.Reverse());
      else if (v2.isconj()) 
	DoGivensMult(c,CONJ(s),v1.Conjugate(),v2.Conjugate());
      else DoGivensMult(c,s,v1,v2); 
    }
  }

  // 
  // Symmetric Givens Mult
  //
  template <class Tg, class T> void GivensHermMult(
      RealType(Tg) c, Tg s, T& d0, T& d1, T& e0) 
  {
    // [ d0 e0* ] = [  c  s ] [ d0 e0* ] [ c  -s ]
    // [ e0 d1  ]   [ -s* c ] [ e0 d1  ] [ s*  c ]
    // = [ c^2 d0 + 2c Re(s e0) + |s|^2 d1   cs(d1-d0) + c^2 e0* - s^2 e0    ]
    //   [ cs*(d1-d0) + c^2 e0 - s*^2 e0*    c^2 d1 - 2c Re(s e0) + |s|^2 d0 ]
    // (using c^2 = 1-|s|^2:)
    // d0' = d0 + 2c Re(s e0) + |s|^2 (d1-d0)
    // e0' = e0 - 2s* Re(s e0) + c s* (d1-d0)
    // d1' = d1 - 2c Re(s e0) - |s|^2 (d1-d0)

    RealType(Tg) s2 = NORM(s);
    RealType(T) Rese0 = REAL(s*e0);
    T d1md0 = d1-d0;
    T dd = RealType(Tg)(2)*c*Rese0 + s2*d1md0;
    d0 += dd;
    d1 -= dd;
    e0 += CONJ(s)*(c*d1md0 - RealType(Tg)(2)*Rese0);
  }
  template <class Tg, class T> void GivensSymMult(
      RealType(Tg) c, Tg s, T& d0, T& d1, T& e0) 
  {
    // [ d0 e0 ] = [  c  s ] [ d0 e0 ] [ c -s* ]
    // [ e0 d1 ]   [ -s* c ] [ e0 d1 ] [ s  c  ]
    // d0' = d0 + 2c (s e0) + s^2 (d1-d0)
    // e0' = e0 - 2|s|^2 e0 + c (s d1 - s* d0)
    // d1' = d1 - 2c (s e0) - s^2 (d1-d0)

    Tg s2 = s*s;
    T se0 = s*e0;
    T d1md0 = d1-d0;
    T dd = RealType(Tg)(2)*c*se0 + s2*d1md0;
    d0 += dd;
    d1 -= dd;
    if (IsReal(Tg()))
      e0 += s*(c*d1md0 - RealType(Tg)(2)*se0);
    else
      e0 += c*(s*d1-CONJ(s)*d0) - RealType(Tg)(2)*CONJ(s)*se0;
  }

#define InstFile "TMV_Givens.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


