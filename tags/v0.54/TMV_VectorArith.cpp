
#include "TMV_Blas.h"
#include "TMV_VectorArith.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t TMV_MULTVV_RECURSE_SIZE = TMV_BLOCKSIZE;
#else
  const size_t TMV_MULTVV_RECURSE_SIZE = 64;
#endif

  //
  // MultXV
  // 
  template <class T, class Tx> inline void DoMultXV(
      const Tx x, const VectorView<T>& v)
  {
    TMVAssert(x!=Tx(0));
    TMVAssert(x!=Tx(1)); 
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>0);
    TMVAssert(v.ct() == NonConj);
    
    T* vptr = v.ptr();
    const int s = v.step();
    const size_t N = v.size();
    if (s == 1) {
      const size_t N1 = N/4;
      const size_t N2 = N-4*N1;
      if (N1) for(size_t i=N1;i>0;--i,vptr+=4) {
	*vptr *= x;
	vptr[1] *= x;
	vptr[2] *= x;
	vptr[3] *= x;
      }
      if (N2) for(size_t i=N2;i>0;--i,++vptr) *vptr *= x;
    }
    else
      for(size_t i=N;i>0;--i,vptr+=s) *vptr *= x;
  }

#ifdef BLAS
  template <> inline void DoMultXV(const double x, 
      const VectorView<double>& v)
  {
    TMVAssert(x!=double(0));
    TMVAssert(x!=double(1));
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>0);
    int n=v.size();
    int s=v.step();
    BLASNAME(dscal) (BLASV(n),BLASV(x),BLASP(v.ptr()),BLASV(s));
  }
  template <> inline void DoMultXV(const double x,
      const VectorView<complex<double> >& v)
  {
    TMVAssert(x!=double(0));
    TMVAssert(x!=double(1));
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>0);
    TMVAssert(v.ct() == NonConj);
    int n=v.size();
    int s=v.step();
    BLASNAME(zdscal) (BLASV(n),BLASV(x),BLASP(v.ptr()),BLASV(s));
  }
  template <> inline void DoMultXV(const complex<double> x,
      const VectorView<complex<double> >& v)
  { 
    TMVAssert(x!=double(0));
    TMVAssert(x!=double(1));
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>0);
    TMVAssert(v.ct() == NonConj);
    int n=v.size();
    int s=v.step();
    BLASNAME(zscal) (BLASV(n),BLASP(&x),BLASP(v.ptr()),BLASV(s));
  }
#ifndef NOFLOAT
  template <> inline void DoMultXV(const float x, 
      const VectorView<float>& v)
  {
    TMVAssert(x!=float(0));
    TMVAssert(x!=float(1));
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>0);
    int n=v.size();
    int s=v.step();
    BLASNAME(sscal) (BLASV(n),BLASV(x),BLASP(v.ptr()),BLASV(s));
  }
  template <> inline void DoMultXV(const float x,
      const VectorView<complex<float> >& v)
  {
    TMVAssert(x!=float(0));
    TMVAssert(x!=float(1));
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>0);
    TMVAssert(v.ct() == NonConj);
    int n=v.size();
    int s=v.step();
    BLASNAME(csscal) (BLASV(n),BLASV(x),BLASP(v.ptr()),BLASV(s));
  }
  template <> inline void DoMultXV(const complex<float> x,
      const VectorView<complex<float> >& v)
  { 
    TMVAssert(x!=float(0));
    TMVAssert(x!=float(1));
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>0);
    TMVAssert(v.ct() == NonConj);
    int n=v.size();
    int s=v.step();
    BLASNAME(cscal) (BLASV(n),BLASP(&x),BLASP(v.ptr()),BLASV(s));
  }
#endif
#endif

  template <class T> void MultXV(const T x, 
      const VectorView<T>& v)
  { 
    TMVAssert(v.step() != 0);

#ifdef XDEBUG
    Vector<T> v0 = v;
    Vector<T> vx = v;
    for(size_t i=0;i<vx.size();i++) vx(i) *= x;
#endif

    if (v.size() > 0 && x != T(1)) {
      if (v.step() < 0) MultXV(x,v.Reverse());
      else if (v.isconj()) MultXV(CONJ(x),v.Conjugate());
      else if (x == T(0)) v.Zero();
      else if (IsComplex(T()) && IMAG(x) == RealType(T)(0))
	if (v.step() == 1) DoMultXV(REAL(x),v.Flatten());
	else DoMultXV(REAL(x),v);
      else DoMultXV(x,v); 
    }

#ifdef XDEBUG
    if (Norm(vx-v) > 0.001*max(RealType(T)(1),Norm(v))) {
      cerr<<"MultXV: x = "<<x<<endl;
      cerr<<"v = "<<Type(v)<<"  step "<<v.step()<<"  "<<v0<<endl;
      cerr<<"-> "<<v<<endl;
      cerr<<"vx = "<<vx<<endl;
      cerr<<"Norm(vx-v) = "<<Norm(vx-v)<<endl;
      abort();
    }
#endif
  }

  template <bool c1, class T, class Tx, class T1> inline void DoMultXV(
      const Tx x, const GenVector<T1>& v1, const VectorView<T>& v2)
  {
    TMVAssert(v2.size()==v1.size());
    TMVAssert(v2.size()>0);
    TMVAssert(x!=Tx(0));
    TMVAssert(x!=Tx(1)); 
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);
    TMVAssert(v2.ct() == NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(IsReal(x) || IMAG(x) != RealType(Tx)(0));
    TMVAssert(!(v2.SameAs(v1)));
    TMVAssert(c1 == v1.isconj());

    if (v2.SameAs(v1)) DoMultXV(x,v2);
    else {
      const T1* v1ptr = v1.cptr();
      T* v2ptr = v2.ptr();
      const int s1 = v1.step();
      const int s2 = v2.step();
      const size_t N = v1.size();

      if (s1 == 1 && s2 == 1) {
	const size_t N1 = N/4;
	const size_t N2 = N-4*N1;
	if (N1) for(size_t i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
	  *v2ptr = x * (c1 ? CONJ(*v1ptr) : (*v1ptr));
	  v2ptr[1] = x * (c1 ? CONJ(v1ptr[1]) : v1ptr[1]);
	  v2ptr[2] = x * (c1 ? CONJ(v1ptr[2]) : v1ptr[2]);
	  v2ptr[3] = x * (c1 ? CONJ(v1ptr[3]) : v1ptr[3]);
	}
	if (N2) for(size_t i=N2;i>0;--i,v1ptr++,v2ptr++) 
	  *v2ptr = x * (c1 ? CONJ(*v1ptr) : (*v1ptr));
      }
      else
	for(size_t i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) 
	  *v2ptr = x * (c1 ? CONJ(*v1ptr) : (*v1ptr));
    }
  }

  template <class T, class T1> void MultXV(const T x, 
      const GenVector<T1>& v1, const VectorView<T>& v2)
  { 
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);
    TMVAssert(v1.size() == v2.size());

#ifdef XDEBUG
    Vector<T> vx = v1;
    for(size_t i=0;i<vx.size();i++) vx(i) *= x;
#endif

    if (v2.size() > 0) {

      const int step1 = v1.step();
      const int step2 = v2.step();

      if (ShouldReverse(step1,step2)) MultXV(x,v1.Reverse(),v2.Reverse());
      else if (v2.isconj()) MultXV(CONJ(x),v1.Conjugate(),v2.Conjugate());
      else if (x == T(0)) v2.Zero();
      else if (v2.SameAs(v1)) MultXV(x,v2);
      else if (x == T(1)) DoCopy(v1,v2);
      else if (IsComplex(T()) && IMAG(x)==RealType(T)(0))
	if (IsComplex(T1()) && v2.isconj() == v1.isconj() &&
	    (v1.step()==1 && v2.step()==1))
	  DoMultXV<false>(REAL(x),v1.Flatten(),v2.Flatten());
	else if (v1.isconj()) DoMultXV<true>(REAL(x),v1,v2);
	else DoMultXV<false>(REAL(x),v1,v2);
      else if (v1.isconj()) DoMultXV<true>(x,v1,v2);
      else DoMultXV<false>(x,v1,v2);
    }

#ifdef XDEBUG
    if (Norm(vx-v2) > 0.001*max(RealType(T)(1),Norm(v2))) {
      cerr<<"CopyMultXV: x = "<<x<<endl;
      cerr<<"v2 = "<<Type(v2)<<"  step "<<v2.step()<<endl;
      cerr<<"v1 = "<<Type(v1)<<"  step "<<v1.step()<<endl;
      cerr<<"-> "<<v2<<endl;
      cerr<<"vx = "<<vx<<endl;
      cerr<<"Norm(vx-v2) = "<<Norm(vx-v2)<<endl;
      abort();
    }
#endif
  }

  // 
  // AddVV
  //

  template <bool c1, class T, class Tx, class T1> inline void NonBlasAddVV(
      const Tx x, const GenVector<T1>& v1, const VectorView<T>& v2)
  {
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v2.size()>0);
    TMVAssert(x != Tx(0));
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);
    TMVAssert(v2.ct() == NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(IsReal(x) || IMAG(x) != RealType(Tx)(0));
    TMVAssert(c1 == v1.isconj());

    const T1* v1ptr = v1.cptr();
    T* v2ptr = v2.ptr();
    const int s1 = v1.step();
    const int s2 = v2.step();
    const size_t N = v2.size();

    if (s1 == 1 && s2 == 1) {
      const size_t N1 = N/4;
      const size_t N2 = N-4*N1;
      if (N1) {
	if (x == Tx(1))        
	  for(size_t i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
	    *v2ptr += (c1 ? CONJ(*v1ptr) : (*v1ptr));
	    v2ptr[1] += (c1 ? CONJ(v1ptr[1]) : v1ptr[1]);
	    v2ptr[2] += (c1 ? CONJ(v1ptr[2]) : v1ptr[2]);
	    v2ptr[3] += (c1 ? CONJ(v1ptr[3]) : v1ptr[3]);
	  }
	else if (x == Tx(-1))
	  for(size_t i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
	    *v2ptr -= (c1 ? CONJ(*v1ptr) : (*v1ptr));
	    v2ptr[1] -= (c1 ? CONJ(v1ptr[1]) : v1ptr[1]);
	    v2ptr[2] -= (c1 ? CONJ(v1ptr[2]) : v1ptr[2]);
	    v2ptr[3] -= (c1 ? CONJ(v1ptr[3]) : v1ptr[3]);
	  }
	else
	  for(size_t i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
	    *v2ptr += x * (c1 ? CONJ(*v1ptr) : (*v1ptr));
	    v2ptr[1] += x * (c1 ? CONJ(v1ptr[1]) : v1ptr[1]);
	    v2ptr[2] += x * (c1 ? CONJ(v1ptr[2]) : v1ptr[2]);
	    v2ptr[3] += x * (c1 ? CONJ(v1ptr[3]) : v1ptr[3]);
	  }
      }
      if (N2) {
	if (x == Tx(1))        
	  for(size_t i=N2;i>0;--i,++v1ptr,++v2ptr) 
	    *v2ptr += (c1 ? CONJ(*v1ptr) : (*v1ptr));
	else if (x == Tx(-1))
	  for(size_t i=N2;i>0;--i,++v1ptr,++v2ptr) 
	    *v2ptr -= (c1 ? CONJ(*v1ptr) : (*v1ptr));
	else
	  for(size_t i=N2;i>0;--i,++v1ptr,++v2ptr) 
	    *v2ptr += x * (c1 ? CONJ(*v1ptr) : (*v1ptr));
      }
    } else {
      if (x == Tx(1))        
	for(size_t i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) 
	  *v2ptr += (c1 ? CONJ(*v1ptr) : (*v1ptr));
      else if (x == Tx(-1))
	for(size_t i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) 
	  *v2ptr -= (c1 ? CONJ(*v1ptr) : (*v1ptr));
      else
	for(size_t i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) 
	  *v2ptr += x * (c1 ? CONJ(*v1ptr) : (*v1ptr));
    }
  }

#ifdef BLAS
  template <class T, class Tx, class T1> inline void BlasAddVV(
      const Tx x, const GenVector<T1>& v1, const VectorView<T>& v2)
  { 
    if (v1.isconj()) NonBlasAddVV<true>(x,v1,v2); 
    else NonBlasAddVV<false>(x,v1,v2); 
  }
  template <> inline void BlasAddVV(const double x,
      const GenVector<double>& v1, const VectorView<double>& v2)
  { 
    TMVAssert(v1.size() == v1.size()); 
    TMVAssert(v2.size()>0);
    TMVAssert(x != double(0));
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    TMVAssert(v1.step()>0);
    TMVAssert(v2.step()>0);
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(daxpy) (BLASV(n),BLASV(x),BLASP(v1.cptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2));
  }
  template <> inline void BlasAddVV(const complex<double> x, 
      const GenVector<complex<double> >& v1, 
      const VectorView<complex<double> >& v2)
  { 
    TMVAssert(v1.size() == v1.size()); 
    TMVAssert(v2.size()>0);
    TMVAssert(x != double(0));
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    TMVAssert(v1.step()>0);
    TMVAssert(v2.step()>0);
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(zaxpy) (BLASV(n),BLASP(&x),BLASP(v1.cptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2));
  }
#ifndef NOFLOAT
  template <> inline void BlasAddVV(const float x,
      const GenVector<float>& v1, const VectorView<float>& v2)
  {
    TMVAssert(v1.size() == v1.size()); 
    TMVAssert(v2.size()>0);
    TMVAssert(x != float(0));
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    TMVAssert(v1.step()>0);
    TMVAssert(v2.step()>0);
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(saxpy) (BLASV(n),BLASV(x),BLASP(v1.cptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2));
  }
  template <> inline void BlasAddVV(const complex<float> x, 
      const GenVector<complex<float> >& v1, 
      const VectorView<complex<float> >& v2)
  { 
    TMVAssert(v1.size() == v1.size()); 
    TMVAssert(v2.size()>0);
    TMVAssert(x != float(0));
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    TMVAssert(v1.step()>0);
    TMVAssert(v2.step()>0);
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(caxpy) (BLASV(n),BLASP(&x),BLASP(v1.cptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2));
  }
#endif // FLOAT
#endif // BLAS

  template <class T, class Tx, class T1> inline void DoAddVV(
      const Tx x, const GenVector<T1>& v1, const VectorView<T>& v2)
  {
    TMVAssert(v1.size() == v2.size()); 
    TMVAssert(v2.size() > 0);
    TMVAssert(x != Tx(0));
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);
    TMVAssert(v2.ct() == NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(IsReal(x) || IMAG(x) != RealType(Tx)(0));
    if (v2.SameAs(v1)) DoMultXV(x+Tx(1),v2);
    else
#ifdef BLAS
      if (v1.step() > 0 && v2.step() > 0 && v1.ct() == NonConj)
	BlasAddVV(x,v1,v2);
      else
#endif
	if (v1.isconj()) NonBlasAddVV<true>(x,v1,v2); 
	else NonBlasAddVV<false>(x,v1,v2); 
  }

  template <class T, class T1> void AddVV(
      const T x, const GenVector<T1>& v1, const VectorView<T>& v2)
  { 
    TMVAssert(v1.size() == v2.size()); 
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);

#ifdef XDEBUG
    Vector<T> v0 = v2;
    Vector<T> vx = v2;
    for(size_t i=0;i<v2.size();i++) vx(i) += x * v1(i);
#endif

    if (v2.size() > 0 && x != T(0)) {
      if (ShouldReverse(v1.step(),v2.step())) 
	AddVV(x,v1.Reverse(),v2.Reverse());
      else if (IsComplex(T()) && IMAG(x)==RealType(T)(0))
	if (IsComplex(T1()) && v2.isconj() == v1.isconj() &&
	    v1.step()==1 && v2.step()==1)
	  DoAddVV(REAL(x),v1.Flatten(),v2.Flatten());
	else if (v2.isconj()) DoAddVV(REAL(x),v1.Conjugate(),v2.Conjugate());
	else DoAddVV(REAL(x),v1,v2);
      else 
	if (v2.isconj()) DoAddVV(CONJ(x),v1.Conjugate(),v2.Conjugate());
	else DoAddVV(x,v1,v2);
    }

#ifdef XDEBUG
    // Can't just say Norm(vx-v2) as usual, since that would call this function
    RealType(T) nm(0);
    for(size_t i=0;i<v2.size();i++) nm += NORM(vx(i)-v2(i));
    if (nm > 0.001*abs(x)*Norm(v1)*Norm(v0)) {
      cerr<<"AddVV: x = "<<x<<endl;
      cerr<<"v1 = "<<Type(v1)<<"  step "<<v1.step()<<"  "<<v1<<endl;
      cerr<<"v2 = "<<Type(v2)<<"  step "<<v2.step()<<"  "<<v0<<endl;
      cerr<<"-> "<<v2<<endl;
      cerr<<"vx = "<<vx<<endl;
      cerr<<"Norm(vx-v2) = "<<nm<<endl;
      abort();
    }
#endif
  }

  template <bool c1, class T, class Tx1, class Tx2, class T1> 
    inline void DoAddVV(
	const Tx1 x1, const GenVector<T1>& v1, const Tx2 x2,
	const VectorView<T>& v2)
    {
      TMVAssert(v1.size() == v2.size());
      TMVAssert(v2.size()>0);
      TMVAssert(x1 != Tx1(0));
      TMVAssert(x2 != Tx2(0));
      TMVAssert(x2 != Tx2(1));
      TMVAssert(v1.step() != 0);
      TMVAssert(v2.step() != 0);
      TMVAssert(v2.ct() == NonConj);
      TMVAssert(v2.step() != -1);
      TMVAssert(v1.step() != -1 || v2.step() == 1);
      TMVAssert(v2.step() > 0 || v1.step() == 1);
      TMVAssert(IsReal(x1) || IMAG(x1) != RealType(Tx1)(0));
      TMVAssert(IsReal(x2) || IMAG(x2) != RealType(Tx2)(0));

      if (v2.SameAs(v1)) {
	if (v2.step() < 0) DoMultXV(x1+x2,v2.Reverse());
	else DoMultXV(x1+x2,v2);
      } else {
	const T1* v1ptr = v1.cptr();
	T* v2ptr = v2.ptr();
	const int s1 = v1.step();
	const int s2 = v2.step();
	const size_t N = v2.size();

	if (s1 == 1 && s2 == 1) {
	  const size_t N1 = N/4;
	  const size_t N2 = N-4*N1;
	  if (N1) {
	    if (x1 == Tx1(1))        
	      if (x2 == Tx2(-1))
		for(size_t i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
		  (*v2ptr) = (c1 ? CONJ(*v1ptr) : (*v1ptr)) - (*v2ptr);
		  v2ptr[1] = (c1 ? CONJ(v1ptr[1]) : v1ptr[1]) - v2ptr[1];
		  v2ptr[2] = (c1 ? CONJ(v1ptr[2]) : v1ptr[2]) - v2ptr[2];
		  v2ptr[3] = (c1 ? CONJ(v1ptr[3]) : v1ptr[3]) - v2ptr[3];
		}
	      else
		for(size_t i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
		  (*v2ptr) = x2*(*v2ptr) + (c1 ? CONJ(*v1ptr) : (*v1ptr));
		  v2ptr[1] = x2*v2ptr[1] + (c1 ? CONJ(v1ptr[1]) : v1ptr[1]);
		  v2ptr[2] = x2*v2ptr[2] + (c1 ? CONJ(v1ptr[2]) : v1ptr[2]);
		  v2ptr[3] = x2*v2ptr[3] + (c1 ? CONJ(v1ptr[3]) : v1ptr[3]);
		}
	    else if (x1 == Tx1(-1))
	      if (x2 == Tx2(-1))
		for(size_t i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
		  (*v2ptr) = -(*v2ptr) - (c1 ? CONJ(*v1ptr) : (*v1ptr));
		  v2ptr[1] = -v2ptr[1] - (c1 ? CONJ(v1ptr[1]) : v1ptr[1]);
		  v2ptr[2] = -v2ptr[2] - (c1 ? CONJ(v1ptr[2]) : v1ptr[2]);
		  v2ptr[3] = -v2ptr[3] - (c1 ? CONJ(v1ptr[3]) : v1ptr[3]);
		}
	      else
		for(size_t i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
		  (*v2ptr) = x2*(*v2ptr) - (c1 ? CONJ(*v1ptr) : (*v1ptr));
		  v2ptr[1] = x2*v2ptr[1] - (c1 ? CONJ(v1ptr[1]) : v1ptr[1]);
		  v2ptr[2] = x2*v2ptr[2] - (c1 ? CONJ(v1ptr[2]) : v1ptr[2]);
		  v2ptr[3] = x2*v2ptr[3] - (c1 ? CONJ(v1ptr[3]) : v1ptr[3]);
		}
	    else
	      if (x2 == Tx2(-1))
		for(size_t i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
		  (*v2ptr) = -(*v2ptr) + x1*(c1 ? CONJ(*v1ptr) : (*v1ptr));
		  v2ptr[1] = -v2ptr[1] + x1*(c1 ? CONJ(v1ptr[1]) : v1ptr[1]);
		  v2ptr[2] = -v2ptr[2] + x1*(c1 ? CONJ(v1ptr[2]) : v1ptr[2]);
		  v2ptr[3] = -v2ptr[3] + x1*(c1 ? CONJ(v1ptr[3]) : v1ptr[3]);
		}
	      else
		for(size_t i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
		  (*v2ptr) = x2*(*v2ptr) + x1*(c1 ? CONJ(*v1ptr) : (*v1ptr));
		  v2ptr[1] = x2*v2ptr[1] + x1*(c1 ? CONJ(v1ptr[1]) : v1ptr[1]);
		  v2ptr[2] = x2*v2ptr[2] + x1*(c1 ? CONJ(v1ptr[2]) : v1ptr[2]);
		  v2ptr[3] = x2*v2ptr[3] + x1*(c1 ? CONJ(v1ptr[3]) : v1ptr[3]);
		}
	  }
	  if (N2) {
	    if (x1 == Tx1(1))        
	      if (x2 == Tx2(-1))
		for(size_t i=N2;i>0;--i,++v1ptr,++v2ptr) 
		  (*v2ptr) = (c1 ? CONJ(*v1ptr) : (*v1ptr)) - (*v2ptr);
	      else
		for(size_t i=N2;i>0;--i,++v1ptr,++v2ptr) 
		  (*v2ptr) = x2*(*v2ptr) + (c1 ? CONJ(*v1ptr) : (*v1ptr));
	    else if (x1 == Tx1(-1))
	      if (x2 == Tx2(-1))
		for(size_t i=N2;i>0;--i,++v1ptr,++v2ptr) 
		  (*v2ptr) = -(*v2ptr) - (c1 ? CONJ(*v1ptr) : (*v1ptr));
	      else
		for(size_t i=N2;i>0;--i,++v1ptr,++v2ptr) 
		  (*v2ptr) = x2*(*v2ptr) - (c1 ? CONJ(*v1ptr) : (*v1ptr));
	    else
	      if (x2 == Tx2(-1))
		for(size_t i=N2;i>0;--i,++v1ptr,++v2ptr) 
		  (*v2ptr) = x1*(c1 ? CONJ(*v1ptr) : (*v1ptr)) - (*v2ptr);
	      else
		for(size_t i=N2;i>0;--i,++v1ptr,++v2ptr) 
		  (*v2ptr) = x2*(*v2ptr) + x1*(c1 ? CONJ(*v1ptr) : (*v1ptr));
	  }
	} else {
	  if (x1 == Tx1(1))        
	    if (x2 == Tx2(-1))
	      for(size_t i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) 
		(*v2ptr) = (c1 ? CONJ(*v1ptr) : (*v1ptr)) - (*v2ptr);
	    else
	      for(size_t i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) 
		(*v2ptr) = x2*(*v2ptr) + (c1 ? CONJ(*v1ptr) : (*v1ptr));
	  else if (x1 == Tx1(-1))
	    if (x2 == Tx2(-1))
	      for(size_t i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) 
		(*v2ptr) = -(*v2ptr) - (c1 ? CONJ(*v1ptr) : (*v1ptr));
	    else
	      for(size_t i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) 
		(*v2ptr) = x2*(*v2ptr) - (c1 ? CONJ(*v1ptr) : (*v1ptr));
	  else
	    if (x2 == Tx2(-1))
	      for(size_t i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) 
		(*v2ptr) = x1*(c1 ? CONJ(*v1ptr) : (*v1ptr)) - (*v2ptr);
	    else
	      for(size_t i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) 
		(*v2ptr) = x2*(*v2ptr) + x1*(c1 ? CONJ(*v1ptr) : (*v1ptr));
	}
      }
    }

  template <class T, class Tx1, class Tx2, class T1> inline void DoAddVV1(
      const Tx1 x1, const GenVector<T1>& v1, const Tx2 x2,
      const VectorView<T>& v2)
  {
    if (v1.isconj()) DoAddVV<true>(x1,v1,x2,v2);
    else DoAddVV<false>(x1,v1,x2,v2);
  }

  template <class T, class T1> void AddVV(
      const T x1, const GenVector<T1>& v1,
      const T x2, const VectorView<T>& v2)
  { 
    TMVAssert(v2.size() == v1.size()); 
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);

#ifdef XDEBUG
    Vector<T> v0 = v2;
    Vector<T> vx = v2;
    for(size_t i=0;i<vx.size();i++) {
      vx(i) *= x2;
      vx(i) += x1*v1(i);
    }
#endif

    if (v2.size() > 0) {
      const int step1 = v1.step();
      const int step2 = v2.step();

      if ( (step2 < 0 && (step1 != 1 || step2 == -1)) || 
	  (step1 == -1 && step2 != 1) )
	AddVV(x1,v1.Reverse(),x2,v2.Reverse());
      else if (x1 == T(0)) MultXV(x2,v2);
      else if (x2 == T(0)) MultXV(x1,v1,v2);
      else if (x2 == T(1)) AddVV(x1,v1,v2);
      else if (IsComplex(T()))
	if (IMAG(x1)==RealType(T)(0))
	  if (IMAG(x2)==RealType(T)(0))
	    if (IsComplex(T1()) && v2.isconj() == v1.isconj() &&
		v1.step() == 1 && v2.step() == 1)
	      DoAddVV<false>(REAL(x1),v1.Flatten(),REAL(x2),v2.Flatten());
	    else if (v2.isconj()) 
	      DoAddVV1(REAL(x1),v1.Conjugate(),REAL(x2),v2.Conjugate());
	    else DoAddVV1(REAL(x1),v1,REAL(x2),v2);
	  else
	    if (v2.isconj()) 
	      DoAddVV1(REAL(x1),v1.Conjugate(),x2,v2.Conjugate());
	    else DoAddVV1(REAL(x1),v1,x2,v2);
	else
	  if (IMAG(x2)==RealType(T)(0))
	    if (v2.isconj()) 
	      DoAddVV1(x1,v1.Conjugate(),REAL(x2),v2.Conjugate());
	    else DoAddVV1(x1,v1,REAL(x2),v2);
	  else
	    if (v2.isconj()) 
	      DoAddVV1(x1,v1.Conjugate(),x2,v2.Conjugate());
	    else DoAddVV1(x1,v1,x2,v2);
      else DoAddVV<false>(x1,v1,x2,v2);
    }

#ifdef XDEBUG
    if (Norm(vx-v2) > 0.001*(
	  x1==T(0)?RealType(T)(0):(abs(x1)*Norm(v1)) *
	  x2==T(0)?RealType(T)(0):(abs(x2)*Norm(v0)) )) {
      cerr<<"AddVV: x1 = "<<x1<<"  x2 = "<<x2<<endl;
      cerr<<"v1 = "<<Type(v1)<<"  step "<<v1.step()<<"  "<<v1<<endl;
      cerr<<"v2 = "<<Type(v2)<<"  step "<<v2.step()<<"  "<<v0<<endl;
      cerr<<"-> "<<v2<<endl;
      cerr<<"vx = "<<vx<<endl;
      cerr<<"Norm(vx-v2) = "<<Norm(vx-v2)<<endl;
      abort();
    }
#endif
  }

  //
  // MultVV
  //
  template <bool unit, bool c2, class T, class T2> inline T NonBlasMultVV(
      const GenVector<T>& v1, const GenVector<T2>& v2) 
  {
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.size()>0);
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);
    TMVAssert(v1.ct() == NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(c2 == v2.isconj());

    const T* v1ptr = v1.cptr();
    const T2* v2ptr = v2.cptr();

    const size_t N = v1.size();
    if (N > TMV_MULTVV_RECURSE_SIZE) {
      // This isn't for speed reasons - it's for increased accuracy.
      // For large vectors, the incremental additions can be much smaller
      // than the running sum, so the relative errors can be huge.
      // With the recursive algorithm, the relative error is generally
      // closer to the expected few * epsilon.
      const size_t N1 = N/2;
      return NonBlasMultVV<unit,c2>(v1.SubVector(0,N1),v2.SubVector(0,N1)) +
	NonBlasMultVV<unit,c2>(v1.SubVector(N1,N),v2.SubVector(N1,N));
    } else {
      T res(0);

      if (unit) {
	const size_t N1 = N/4;
	const size_t N2 = N-4*N1;
	if (N1) for(size_t i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
	  res += (*v1ptr) * (c2 ? CONJ(*v2ptr) : (*v2ptr));
	  res += v1ptr[1] * (c2 ? CONJ(v2ptr[1]) : v2ptr[1]);
	  res += v1ptr[2] * (c2 ? CONJ(v2ptr[2]) : v2ptr[2]);
	  res += v1ptr[3] * (c2 ? CONJ(v2ptr[3]) : v2ptr[3]);
	}
	if (N2) for(size_t i=N2;i>0;--i,++v1ptr,++v2ptr) 
	  res += (*v1ptr) * (c2 ? CONJ(*v2ptr) : (*v2ptr));
      } else {
	const int s1 = v1.step();
	const int s2 = v2.step();

	for(size_t i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) 
	  res += (*v1ptr) * (c2 ? CONJ(*v2ptr) : (*v2ptr));
      }
      return res;
    }
  }
#ifdef BLAS
  template <class T, class T2> inline T BlasMultVV(
      const GenVector<T>& v1, const GenVector<T2>& v2) 
  { 
    if (v1.step() == 1 && v2.step() == 1)
      if (v2.isconj()) return NonBlasMultVV<true,true>(v1,v2); 
      else return NonBlasMultVV<true,false>(v1,v2); 
    else
      if (v2.isconj()) return NonBlasMultVV<false,true>(v1,v2); 
      else return NonBlasMultVV<false,false>(v1,v2); 
  }
  template <> inline double BlasMultVV(
      const GenVector<double>& v1, const GenVector<double>& v2) 
  { 
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.size()>0);
    TMVAssert(v1.step()>0);
    TMVAssert(v2.step()>0);
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    return BLASNAME(ddot) (BLASV(n),BLASP(v1.cptr()),BLASV(s1),
	BLASP(v2.cptr()),BLASV(s2));
  }
  template <> inline complex<double> BlasMultVV(
      const GenVector<complex<double> >& v1, 
      const GenVector<complex<double> >& v2) 
  { 
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.size()>0);
    TMVAssert(v1.step()>0);
    TMVAssert(v2.step()>0);
    TMVAssert(v1.ct() == NonConj);
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    complex<double> res;
    if (v2.isconj())
      BLASZDOTSET( res, BLASZDOTNAME(zdotc) (
	    BLASZDOT1(BLASP(&res))
	    BLASV(n),BLASP(v2.cptr()),BLASV(s2),
	    BLASP(v1.cptr()),BLASV(s1)
	    BLASZDOT2(BLASP(&res)) ));
    else
      BLASZDOTSET( res, BLASZDOTNAME(zdotu) (
	    BLASZDOT1(BLASP(&res))
	    BLASV(n),BLASP(v2.cptr()),BLASV(s2),
	    BLASP(v1.cptr()),BLASV(s1)
	    BLASZDOT2(BLASP(&res)) ));
    return res;
  }
#ifndef NOFLOAT
  template <> inline float BlasMultVV(
      const GenVector<float>& v1, const GenVector<float>& v2) 
  { 
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.size()>0);
    TMVAssert(v1.step()>0);
    TMVAssert(v2.step()>0);
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    return BLASNAME(sdot) (BLASV(n),BLASP(v1.cptr()),BLASV(s1),
	BLASP(v2.cptr()),BLASV(s2));
  }
  template <> inline complex<float> BlasMultVV(
      const GenVector<complex<float> >& v1, 
      const GenVector<complex<float> >& v2) 
  { 
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.size()>0);
    TMVAssert(v1.step()>0);
    TMVAssert(v2.step()>0);
    TMVAssert(v1.ct() == NonConj);
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    complex<float> res;
    if (v2.isconj())
      BLASZDOTSET( res, BLASZDOTNAME(cdotc) (
	    BLASZDOT1(BLASP(&res))
	    BLASV(n),BLASP(v2.cptr()),BLASV(s2),
	    BLASP(v1.cptr()),BLASV(s1)
	    BLASZDOT2(BLASP(&res)) ));
    else
      BLASZDOTSET( res, BLASZDOTNAME(cdotu) (
	    BLASZDOT1(BLASP(&res))
	    BLASV(n),BLASP(v2.cptr()),BLASV(s2),
	    BLASP(v1.cptr()),BLASV(s1)
	    BLASZDOT2(BLASP(&res)) ));
    return res;
  }
#endif // FLOAT
#endif // BLAS

  template <class T, class T2> inline T DoMultVV(
      const GenVector<T>& v1, const GenVector<T2>& v2) 
  {
    TMVAssert(v1.size() == v2.size()); 
    TMVAssert(v1.size()>0);
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);
    TMVAssert(v1.ct() == NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);

    if (v2.SameAs(v1.Conjugate())) return v1.NormSq();
    else 
#ifdef BLAS
      if (v1.step() > 0 && v2.step() > 0)
	return BlasMultVV(v1,v2);
      else
#endif
	if (v1.step() == 1 && v2.step() == 1)
	  if (v2.isconj()) return NonBlasMultVV<true,true>(v1,v2); 
	  else return NonBlasMultVV<true,false>(v1,v2); 
	else
	  if (v2.isconj()) return NonBlasMultVV<false,true>(v1,v2); 
	  else return NonBlasMultVV<false,false>(v1,v2); 
  }

  template <class T, class T2> T MultVV(
      const GenVector<T>& v1, const GenVector<T2>& v2) 
  { 
    TMVAssert(v1.size() == v2.size()); 
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);

#ifdef XDEBUG
    T resx(0);
    for(size_t i=0;i<v1.size();i++) {
      resx += v1(i)*v2(i);
    }
#endif

    T res(0);
    if (v1.size() > 0) {
      if (ShouldReverse(v1.step(),v2.step())) 
	if (v1.isconj()) 
	  res = CONJ(DoMultVV(v1.Reverse().Conjugate(),
		v2.Reverse().Conjugate()));
	else 
	  res = DoMultVV(v1.Reverse(),v2.Reverse());
      else 
	if (v1.isconj()) 
	  res = CONJ(DoMultVV(v1.Conjugate(),v2.Conjugate()));
	else 
	  res = DoMultVV(v1,v2);
    }

#ifdef XDEBUG
    if (abs(resx-res) > 0.001*max(RealType(T)(1),Norm(v1)*Norm(v2))) {
      cerr<<"MultVV: \n";
      cerr<<"v1 = "<<Type(v1)<<"  step "<<v1.step()<<"  "<<v1<<endl;
      cerr<<"v2 = "<<Type(v2)<<"  step "<<v2.step()<<"  "<<v2<<endl;
      cerr<<"v1*v2 = "<<resx<<endl;
      cerr<<"res = "<<res<<endl;
      cerr<<"abs(resx-res) = "<<abs(resx-res)<<endl;
      abort();
    }
#endif

    return res;
  }

  //
  // ElementProd
  //
 
  template <bool cx, bool cy, class T, class Ta, class Tx, class Ty>
    inline void DoAddElementProd(const Ta alpha, const GenVector<Tx>& x,
	const GenVector<Ty>& y, const VectorView<T>& z)
    // zi += alpha * xi * yi 
    {
      TMVAssert(z.size() == x.size());
      TMVAssert(z.size() == y.size());
      TMVAssert(alpha != Ta(0));
      TMVAssert(z.size()>0);
      TMVAssert(z.ct() == NonConj);
      TMVAssert(z.step() != -1);
      TMVAssert(x.step() != -1 || z.step() == 1);
      TMVAssert(y.step() != -1 || z.step() == 1);
      TMVAssert(z.step() > 0 || x.step() == 1 || y.step() == 1);
      TMVAssert(cx == x.isconj());
      TMVAssert(cy == y.isconj());

      const Tx* xp = x.cptr();
      const Ty* yp = y.cptr();
      T* zp = z.ptr();
      const int sx = x.step();
      const int sy = y.step();
      const int sz = z.step();
      const size_t N = z.size();

      if (sx == 1 && sy == 1 && sz == 1) {
	const size_t N1 = N/4;
	const size_t N2 = N-4*N1;

	if (N1) {
	  if (alpha == Ta(1))
	    for(size_t i=N1;i>0;--i,xp+=4,yp+=4,zp+=4) {
	      *zp += (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
	      zp[1] += (cx?CONJ(xp[1]):xp[1]) * (cy?CONJ(yp[1]):yp[1]);
	      zp[2] += (cx?CONJ(xp[2]):xp[2]) * (cy?CONJ(yp[2]):yp[2]);
	      zp[3] += (cx?CONJ(xp[3]):xp[3]) * (cy?CONJ(yp[3]):yp[3]);
	    }
	  else if (alpha == Ta(-1))
	    for(size_t i=N1;i>0;--i,xp+=4,yp+=4,zp+=4){
	      *zp -= (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
	      zp[1] -= (cx?CONJ(xp[1]):xp[1]) * (cy?CONJ(yp[1]):yp[1]);
	      zp[2] -= (cx?CONJ(xp[2]):xp[2]) * (cy?CONJ(yp[2]):yp[2]);
	      zp[3] -= (cx?CONJ(xp[3]):xp[3]) * (cy?CONJ(yp[3]):yp[3]);
	    }
	  else 
	    for(size_t i=N1;i>0;--i,xp+=4,yp+=4,zp+=4) {
	      *zp += alpha * (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
	      zp[1] += alpha * (cx?CONJ(xp[1]):xp[1]) * (cy?CONJ(yp[1]):yp[1]);
	      zp[2] += alpha * (cx?CONJ(xp[2]):xp[2]) * (cy?CONJ(yp[2]):yp[2]);
	      zp[3] += alpha * (cx?CONJ(xp[3]):xp[3]) * (cy?CONJ(yp[3]):yp[3]);
	    }
	}
	if (N2) {
	  if (alpha == Ta(1))
	    for(size_t i=N2;i>0;--i,++xp,++yp,++zp) 
	      *zp += (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
	  else if (alpha == Ta(-1))
	    for(size_t i=N2;i>0;--i,++xp,++yp,++zp)
	      *zp -= (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
	  else 
	    for(size_t i=N2;i>0;--i,++xp,++yp,++zp) 
	      *zp += alpha * (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
	}
      }
      else {
	if (alpha == Ta(1))
	  for(size_t i=N;i>0;--i,xp+=sx,yp+=sy,zp+=sz) 
	    *zp += (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
	else if (alpha == Ta(-1))
	  for(size_t i=N;i>0;--i,xp+=sx,yp+=sy,zp+=sz)
	    *zp -= (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
	else 
	  for(size_t i=N;i>0;--i,xp+=sx,yp+=sy,zp+=sz) 
	    *zp += alpha * (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
      }
    }

  template <class T, class Tx, class Ty> void AddElementProd(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const VectorView<T>& z)
  {
    TMVAssert(z.size() == x.size());
    TMVAssert(z.size() == y.size());
#ifdef XDEBUG
    Vector<T> z0 = z;
    Vector<T> zx = z;
    for(size_t i=0;i<zx.size();i++) zx(i) += alpha*x(i)*y(i);
#endif

    if (z.size() > 0 && alpha != T(0)) {
      if (z.isconj()) 
	AddElementProd(CONJ(alpha),x.Conjugate(),y.Conjugate(),
	    z.Conjugate());
      else if ( (z.step()<0 && (z.step()==-1 || (x.step()!=1 && y.step()!=1)))
	  || (z.step()!=1 && (x.step()==-1 || y.step()==-1)) )
	AddElementProd(alpha,x.Reverse(),y.Reverse(),z.Reverse());
      else 
	if (IMAG(alpha) == RealType(T)(0))
	  if (x.isconj())
	    if (y.isconj())
	      DoAddElementProd<true,true>(REAL(alpha),x,y,z);
	    else
	      DoAddElementProd<true,false>(REAL(alpha),x,y,z);
	  else
	    if (y.isconj())
	      DoAddElementProd<false,true>(REAL(alpha),x,y,z);
	    else
	      DoAddElementProd<false,false>(REAL(alpha),x,y,z);
	else
	  if (x.isconj())
	    if (y.isconj())
	      DoAddElementProd<true,true>(alpha,x,y,z);
	    else
	      DoAddElementProd<true,false>(alpha,x,y,z);
	  else
	    if (y.isconj())
	      DoAddElementProd<false,true>(alpha,x,y,z);
	    else
	      DoAddElementProd<false,false>(alpha,x,y,z);
    }

#ifdef XDEBUG
    if (Norm(zx-z) > 0.001*(abs(alpha)*(Norm(x)+Norm(y))+Norm(z0))) {
      cerr<<"AddElProd: alpha = "<<alpha<<endl;
      cerr<<"x = "<<Type(x)<<"  step "<<x.step()<<"  "<<x<<endl;
      cerr<<"y = "<<Type(y)<<"  step "<<y.step()<<"  "<<y<<endl;
      cerr<<"z = "<<Type(z)<<"  step "<<z.step()<<"  "<<z0<<endl;
      cerr<<"-> "<<z<<endl;
      cerr<<"zx = "<<zx<<endl;
      cerr<<"Norm(zx-z) = "<<Norm(zx-z)<<endl;
      abort();
    }
#endif
  }

  template <bool cx, class T, class Ta, class Tx> inline void DoElementProd(
      const Ta alpha, const GenVector<Tx>& x, const VectorView<T>& y)
    // yi = alpha * xi * yi
  {
    TMVAssert(x.size() == y.size());
    TMVAssert(y.size()>0);
    TMVAssert(alpha != Ta(0));
    TMVAssert(y.ct() == NonConj);
    TMVAssert(y.step() != -1);
    TMVAssert(x.step() != -1 || y.step() == 1);
    TMVAssert(y.step() > 0 || x.step() == 1);

    const Tx* xp = x.cptr();
    T* yp = y.ptr();
    const int sx = x.step();
    const int sy = y.step();
    const size_t N = y.size();

    if (sx == 1 && sy == 1) {
      const size_t N1 = N/4;
      const size_t N2 = N-4*N1;
      if (N1) {
	if (alpha == Ta(1))
	  for(size_t i=N1;i>0;--i,xp+=4,yp+=4) {
	    *yp *= cx ? CONJ(*xp) : (*xp);
	    yp[1] *= cx ? CONJ(xp[1]) : xp[1];
	    yp[2] *= cx ? CONJ(xp[2]) : xp[2];
	    yp[3] *= cx ? CONJ(xp[3]) : xp[3];
	  }
	else 
	  for(size_t i=N1;i>0;--i,xp+=4,yp+=4) {
	    *yp *= alpha * (cx ? CONJ(*xp) : (*xp));
	    yp[1] *= alpha * (cx ? CONJ(xp[1]) : xp[1]);
	    yp[2] *= alpha * (cx ? CONJ(xp[2]) : xp[2]);
	    yp[3] *= alpha * (cx ? CONJ(xp[3]) : xp[3]);
	  }
      }
      if (N2) {
	if (alpha == Ta(1))
	  for(size_t i=N2;i>0;--i,++xp,++yp) 
	    *yp *= cx ? CONJ(*xp) : (*xp);
	else 
	  for(size_t i=N2;i>0;--i,++xp,++yp) 
	    *yp *= alpha * (cx ? CONJ(*xp) : (*xp));
      }
    }
    else {
      if (alpha == Ta(1))
	for(size_t i=N;i>0;--i,xp+=sx,yp+=sy) 
	  *yp *= cx ? CONJ(*xp) : (*xp);
      else 
	for(size_t i=N;i>0;--i,xp+=sx,yp+=sy) 
	  *yp *= alpha * (cx ? CONJ(*xp) : (*xp));
    }
  }

  template <class T, class Tx> void ElementProd(const T alpha,
      const GenVector<Tx>& x, const VectorView<T>& y)
  {
    TMVAssert(x.size() == y.size());

#ifdef XDEBUG
    Vector<T> y0 = y;
    Vector<T> yx = y;
    for(size_t i=0;i<yx.size();i++) yx(i) *= alpha*x(i);
#endif

    if (y.size() > 0 && alpha != T(0)) 
      if (y.isconj()) ElementProd(CONJ(alpha),x.Conjugate(),y.Conjugate());
      else if (ShouldReverse(x.step(),y.step()))
	ElementProd(alpha,x.Reverse(),y.Reverse());
      else if (IMAG(alpha) == 0) 
	if (x.isconj()) 
	  DoElementProd<true>(REAL(alpha),x,y);
	else 
	  DoElementProd<false>(REAL(alpha),x,y);
      else
	if (x.isconj()) 
	  DoElementProd<true>(alpha,x,y);
	else 
	  DoElementProd<false>(alpha,x,y);

#ifdef XDEBUG
    if (Norm(yx-y) > 0.001*(abs(alpha)*Norm(x)+Norm(y0))) {
      cerr<<"AddElProd: alpha = "<<alpha<<endl;
      cerr<<"x = "<<Type(x)<<"  step "<<x.step()<<"  "<<x<<endl;
      cerr<<"y = "<<Type(y)<<"  step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"-> "<<y<<endl;
      cerr<<"yx = "<<yx<<endl;
      cerr<<"Norm(yx-y) = "<<Norm(yx-y)<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_VectorArith.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


