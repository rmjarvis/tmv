
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
    
    const RealType(T) T1(1);
    const RealType(T) T0(0);

    RealType(T) maxabs_x = MAXABS(x);
    RealType(T) maxabs_y = MAXABS(y);
    if (maxabs_y == T0) {
      return Givens<T>(T1,T(0));
    } else if (maxabs_x == T0) {
      RealType(T) absy = abs(y);
      T s = SIGN(CONJ(y),absy);
      x = absy; y = T0;
      return Givens<T>(T0,s);
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
	  RealType(T) invnormx = T1/NORM(x);
	  T s = x*CONJ(y)*invnormx;
	  y = T0;
	  return Givens<T>(T1,s);
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
	  RealType(T) sqrtfactor = SQRT(T1+n);
	  RealType(T) invsqrtfactor = T1/sqrtfactor;
	  T s = CONJ(yoverx)*invsqrtfactor;
	  x += x*(n/(T1+sqrtfactor));
	  y = T0;
	  RealType(T) c = T1/sqrtfactor;
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
	  y = T0;
	  return Givens<T>(absxovery,s);
	} else {
	  // c = |x/y|/sqrt(1+|x/y|^2)
	  // s = (x/y)/|x/y|/sqrt(1+|x/y|^2)
	  // r = x/|x/y| sqrt(|x|^2+|y|^2)
	  RealType(T) sqrtfactor = SQRT(T1+n);
	  RealType(T) invsqrtfactor = T1/sqrtfactor;
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
  template <class T1, class T2> inline void GivensMult(
      RealType(T1) c, T1 s, T2& x, T2& y) 
  {
    // [ x' ] = [  c  s  ] [ x ] = [  cx+sy  ]
    // [ y' ]   [ -s* c* ] [ y ]   [ c*y-s*x ]
    T2 xx = c*x+s*y;
    y = c*y-tmv::CONJ(s)*x;
    x = xx;
  }
  template <class T1, class T2> inline void GivensMult(
      RealType(T1) c, T1 s, ConjRef<T2> x, T2& y)
  {
    T2 xx = c*x+s*y;
    y = c*y-tmv::CONJ(s)*x;
    x = xx;
  }
  template <class T1, class T2> inline void GivensMult(
      RealType(T1) c, T1 s, T2& x, ConjRef<T2> y)
  {
    T2 xx = c*x+s*y;
    y = c*y-tmv::CONJ(s)*x;
    x = xx;
  }
  template <class T1, class T2> inline void GivensMult(
      RealType(T1) c, T1 s, ConjRef<T2> x, ConjRef<T2> y)
  {
    T2 xx = c*x+s*y;
    y = c*y-tmv::CONJ(s)*x;
    x = xx;
  }
  template <class T1, class T2> inline void GivensMult(
      RealType(T1) c, T1 s, VarConjRef<T2> x, T2& y)
  {
    T2 xx = c*x+s*y;
    y = c*y-tmv::CONJ(s)*x;
    x = xx;
  }
  template <class T1, class T2> inline void GivensMult(
      RealType(T1) c, T1 s, T2& x, VarConjRef<T2> y)
  {
    T2 xx = c*x+s*y;
    y = c*y-tmv::CONJ(s)*x;
    x = xx;
  }
  template <class T1, class T2> inline void GivensMult(
      RealType(T1) c, T1 s, VarConjRef<T2> x, VarConjRef<T2> y)
  {
    T2 xx = c*x+s*y;
    y = c*y-tmv::CONJ(s)*x;
    x = xx;
  }

  //
  // GivensMult vectors
  //
  template <class T1, class T2, StepItType S2, ConjItType C2, StepItType S3> 
    inline void NonBlasGivensMult(RealType(T1) c, T1 s,
	VIt<T2,S2,C2> it0, VIt<T2,S3,NonConj> it1, size_t size)
    {
      const VIt<T2,S2,C2> _end = it0+size;
      for(; it0 != _end; ++it0,++it1) GivensMult(c,s,*it0,*it1);
    }
  template <class T1, class T2> inline void DoNonBlasGivensMult(
      RealType(T1) c, T1 s,
      const VectorView<T2>& v0, const VectorView<T2>& v1)
  {
    TMVAssert(v0.size()==v1.size());
    TMVAssert(v1.ct()==NonConj);

    if (v0.step() == 1)
      if (v1.step() == 1)
	if (v0.isconj())
	  NonBlasGivensMult(c,s,VIt<T2,Unit,Conj>(v0.begin()),
	      VIt<T2,Unit,NonConj>(v1.begin()),v0.size());
	else
	  NonBlasGivensMult(c,s,VIt<T2,Unit,NonConj>(v0.begin()),
	      VIt<T2,Unit,NonConj>(v1.begin()),v0.size());
      else
	if (v0.isconj())
	  NonBlasGivensMult(c,s,VIt<T2,Unit,Conj>(v0.begin()),
	      VIt<T2,Step,NonConj>(v1.begin()),v0.size());
	else
	  NonBlasGivensMult(c,s,VIt<T2,Unit,NonConj>(v0.begin()),
	      VIt<T2,Step,NonConj>(v1.begin()),v0.size());
    else
      if (v1.step() == 1)
	if (v0.isconj())
	  NonBlasGivensMult(c,s,VIt<T2,Step,Conj>(v0.begin()),
	      VIt<T2,Unit,NonConj>(v1.begin()),v0.size());
	else
	  NonBlasGivensMult(c,s,VIt<T2,Step,NonConj>(v0.begin()),
	      VIt<T2,Unit,NonConj>(v1.begin()),v0.size());
      else
	if (v0.isconj())
	  NonBlasGivensMult(c,s,VIt<T2,Step,Conj>(v0.begin()),
	      VIt<T2,Step,NonConj>(v1.begin()),v0.size());
	else
	  NonBlasGivensMult(c,s,VIt<T2,Step,NonConj>(v0.begin()),
	      VIt<T2,Step,NonConj>(v1.begin()),v0.size());
  }

  template <class T1, class T2> inline void BlasGivensMult(
      RealType(T1) c, T1 s,
      const VectorView<T2>& v0, const VectorView<T2>& v1)
  { DoNonBlasGivensMult(c,s,v0,v1); }
#ifdef BLAS
  template <> inline void BlasGivensMult(
      double c, double s,
      const VectorView<double>& v0, 
      const VectorView<double>& v1)
  { 
    TMVAssert(v0.size()==v1.size());
    TMVAssert(v0.ct()==NonConj);
    TMVAssert(v1.ct()==NonConj);
    cblas_drot(v0.size(),v0.ptr(),v0.step(),v1.ptr(),v1.step(),c,s); 
  }
  template <> inline void BlasGivensMult(
      double c, double s,
      const VectorView<complex<double> >& v0,
      const VectorView<complex<double> >& v1)
  { 
    TMVAssert(v0.size()==v1.size());
    TMVAssert(v1.ct()==NonConj);
    if (v0.isconj()) DoNonBlasGivensMult(c,s,v0,v1);
    else cblas_zdrot(v0.size(),v0.ptr(),v0.step(),v1.ptr(),v1.step(),c,s); 
  }
#ifndef NOFLOAT
  template <> inline void BlasGivensMult(
      float c, float s,
      const VectorView<float>& v0, const VectorView<float>& v1)
  { 
    TMVAssert(v0.size()==v1.size());
    TMVAssert(v0.ct()==NonConj);
    TMVAssert(v1.ct()==NonConj);
    cblas_srot(v0.size(),v0.ptr(),v0.step(),v1.ptr(),v1.step(),c,s); 
  }
  template <> inline void BlasGivensMult(
      float c, float s,
      const VectorView<complex<float> >& v0,
      const VectorView<complex<float> >& v1)
  { 
    TMVAssert(v0.size()==v1.size());
    TMVAssert(v1.ct()==NonConj);
    if (v0.isconj()) DoNonBlasGivensMult(c,s,v0,v1);
    else cblas_csrot(v0.size(),v0.ptr(),v0.step(),v1.ptr(),v1.step(),c,s); 
  }
#endif
#ifdef LAP
  template <> inline void BlasGivensMult(
      double c, complex<double> s,
      const VectorView<complex<double> >& v0,
      const VectorView<complex<double> >& v1)
  {
    TMVAssert(v0.size()==v1.size());
    TMVAssert(v1.ct()==NonConj);
    if (v0.isconj()) DoNonBlasGivensMult(c,s,v0,v1);
    else {
      int n = v0.size();
      int s0 = v0.step();
      int s1 = v1.step();
      double cc = c;
      zrot(&n,LAP_Complex(v0.ptr()),&s0,LAP_Complex(v1.ptr()),&s1,
	  &cc,LAP_Complex(&s)); 
    }
  }
#ifndef NOFLOAT
  template <> inline void BlasGivensMult(
      float c, complex<float> s,
      const VectorView<complex<float> >& v0,
      const VectorView<complex<float> >& v1)
  {
    TMVAssert(v0.size()==v1.size());
    TMVAssert(v1.ct()==NonConj);
    if (v0.isconj()) DoNonBlasGivensMult(c,s,v0,v1);
    else {
      int n = v0.size();
      int s0 = v0.step();
      int s1 = v1.step();
      float cc = c;
      crot(&n,LAP_Complex(v0.ptr()),&s0,LAP_Complex(v1.ptr()),&s1,
	  &cc,LAP_Complex(&s)); 
    }
  }
#endif
#endif // LAP
#endif // BLAS

  template <class T1, class T2> inline void GivensMult(
      RealType(T1) c, T1 s,
      const VectorView<T2>& v0, const VectorView<T2>& v1)
  { 
    TMVAssert(v0.size() == v1.size());
    if (v1.isconj()) BlasGivensMult(c,CONJ(s),v0.Conjugate(),v1.Conjugate());
    else BlasGivensMult(c,s,v0,v1); 
  }

#define InstFile "TMV_Givens.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


