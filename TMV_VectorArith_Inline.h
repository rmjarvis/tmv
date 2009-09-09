#include "TMV_VectorArith.h"

namespace tmv {

  //
  // MultXV
  // 
  template <class T1, class T2, StepItType S2> inline void DoMultXV(
      const T1 x, VIt<T2,S2,NonConj> vit, size_t size)
  {
    TMVAssert(x!=T1(0));
    TMVAssert(x!=T1(1));
    TMVAssert(size>0);
    const VIt<T2,S2,NonConj> _end = vit + size;
    for(;vit!=_end;++vit) *vit *= x;
  }
#ifdef BLAS
  template <StepItType S2> inline void DoMultXV(const double x, 
      VIt<double,S2,NonConj> vit, size_t size)
  {
    TMVAssert(x!=double(0));
    TMVAssert(x!=double(1));
    TMVAssert(size>0);
    TMVAssert(vit.step()>0);
    cblas_dscal(size,x,vit.GetP(),vit.step()); 
  }
  template <StepItType S2> inline void DoMultXV(const double x,
      VIt<complex<double>,S2,NonConj> vit, size_t size)
  {
    TMVAssert(x!=double(0));
    TMVAssert(x!=double(1));
    TMVAssert(size>0);
    TMVAssert(vit.step()>0);
    cblas_zdscal(size,x,vit.GetP(),vit.step()); 
  }
  template <StepItType S2> inline void DoMultXV(const complex<double> x,
      VIt<complex<double>,S2,NonConj> vit, size_t size)
  { 
    TMVAssert(x!=double(0));
    TMVAssert(x!=double(1));
    TMVAssert(size>0);
    TMVAssert(vit.step()>0);
    cblas_zscal(size,&x,vit.GetP(),vit.step()); 
  }
#ifndef NOFLOAT
  template <StepItType S2> inline void DoMultXV(const float x,
      VIt<float,S2,NonConj> vit, size_t size)
  { 
    TMVAssert(x!=float(0));
    TMVAssert(x!=float(1));
    TMVAssert(size>0);
    TMVAssert(vit.step()>0);
    cblas_sscal(size,x,vit.GetP(),vit.step()); 
  }
  template <StepItType S2> inline void DoMultXV(const float x,
      VIt<complex<float>,S2,NonConj> vit, size_t size)
  { 
    TMVAssert(x!=float(0));
    TMVAssert(x!=float(1));
    TMVAssert(size>0);
    TMVAssert(vit.step()>0);
    cblas_csscal(size,x,vit.GetP(),vit.step()); 
  }
  template <StepItType S2> inline void DoMultXV(const complex<float> x,
      VIt<complex<float>,S2,NonConj> vit, size_t size)
  { 
    TMVAssert(x!=float(0));
    TMVAssert(x!=float(1));
    TMVAssert(size>0);
    TMVAssert(vit.step()>0);
    cblas_cscal(size,&x,vit.GetP(),vit.step()); 
  }
#endif
#endif

  // 
  // AddVV
  //
  template <class T1, class T2, StepItType S2, ConjItType C2, class T3, StepItType S3>
    inline void NonBlasAddVV(const T1 x1, CVIt<T2,S2,C2> v2it,
	VIt<T3,S3,NonConj> v3it, size_t size)
    {
      TMVAssert(size>0);
      TMVAssert(x1 != T1(0));
      const VIt<T3,S3,NonConj> _end = v3it + size;
      if (x1 == T1(1))        
	for(;v3it!=_end;++v2it,++v3it) (*v3it) += (*v2it);
      else if (x1 == T1(-1))
	for(;v3it!=_end;++v2it,++v3it) (*v3it) -= (*v2it);
      else
	for(;v3it!=_end;++v2it,++v3it) (*v3it) += x1 * (*v2it);
    }
  template <class T1, class T2, StepItType S2, ConjItType C2, class T3, StepItType S3>
    inline void DoAddVV(const T1 x1, CVIt<T2,S2,C2> v2it,
	VIt<T3,S3,NonConj> v3it, size_t size)
    { NonBlasAddVV(x1,v2it,v3it,size); }
#ifdef BLAS
  template <StepItType S2, StepItType S3> inline void DoAddVV(
      const double x, CVIt<double,S2,NonConj> v1it, 
      VIt<double,S3,NonConj> v2it, size_t size)
  { 
    TMVAssert(size>0);
    TMVAssert(x != double(0));
    TMVAssert(v1it.step() > 0);
    TMVAssert(v2it.step() > 0);
    cblas_daxpy(size,x,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step()); 
  }
  template <StepItType S2, StepItType S3> inline void DoAddVV(
      const complex<double> x, CVIt<complex<double>,S2,NonConj> v1it,
      VIt<complex<double>,S3,NonConj> v2it, size_t size)
  {
    TMVAssert(size>0);
    TMVAssert(x != double(0));
    TMVAssert(v1it.step() > 0);
    TMVAssert(v2it.step() > 0);
    cblas_zaxpy(size,&x,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step()); 
  }
  template <StepItType S2, StepItType S3> inline void DoAddVV(
      const complex<double> x, CVIt<double,S2,NonConj> v1it,
      VIt<complex<double>,S3,NonConj> v2it, size_t size)
  {
    TMVAssert(size>0);
    TMVAssert(x != double(0));
    TMVAssert(v1it.step() > 0);
    TMVAssert(v2it.step() > 0);
    if (x.imag() == double(0)) {
      double xr = x.real();
      cblas_daxpy(size,xr,v1it.GetP(),v1it.step(),
	  reinterpret_cast<double*>(v2it.GetP()),2*v2it.step()); 
    } else if (x.real() == double(0)) {
      double xi = x.imag();
      cblas_daxpy(size,xi,v1it.GetP(),v1it.step(),
	  reinterpret_cast<double*>(v2it.GetP())+1,2*v2it.step()); 
    } else NonBlasAddVV(x,v1it,v2it,size);
  }
  template <StepItType S2, StepItType S3> inline void DoAddVV(
      const double x, CVIt<double,S2,NonConj> v1it,
      VIt<complex<double>,S3,NonConj> v2it, size_t size)
  {
    TMVAssert(size>0);
    TMVAssert(x != double(0));
    TMVAssert(v1it.step() > 0);
    TMVAssert(v2it.step() > 0);
    cblas_daxpy(size,x,v1it.GetP(),v1it.step(),
	reinterpret_cast<double*>(v2it.GetP()),2*v2it.step()); 
  }
#ifndef NOFLOAT
  template <StepItType S2, StepItType S3> inline void DoAddVV(
      const float x, CVIt<float,S2,NonConj> v1it,
      VIt<float,S3,NonConj> v2it, size_t size)
  {
    TMVAssert(size>0);
    TMVAssert(x != float(0));
    TMVAssert(v1it.step() > 0);
    TMVAssert(v2it.step() > 0);
    cblas_saxpy(size,x,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step()); 
  }
  template <StepItType S2, StepItType S3> inline void DoAddVV(
      const complex<float> x, CVIt<complex<float>,S2,NonConj> v1it,
      VIt<complex<float>,S3,NonConj> v2it, size_t size)
  { 
    TMVAssert(size>0);
    TMVAssert(x != float(0));
    TMVAssert(v1it.step() > 0);
    TMVAssert(v2it.step() > 0);
    cblas_caxpy(size,&x,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step()); 
  }
  template <StepItType S2, StepItType S3> inline void DoAddVV(
      const complex<float> x, CVIt<float,S2,NonConj> v1it,
      VIt<complex<float>,S3,NonConj> v2it, size_t size)
  {
    TMVAssert(size>0);
    TMVAssert(x != float(0));
    TMVAssert(v1it.step() > 0);
    TMVAssert(v2it.step() > 0);
    if (x.imag() == float(0)) {
      float xr = x.real();
      cblas_saxpy(size,xr,v1it.GetP(),v1it.step(),
	  reinterpret_cast<float*>(v2it.GetP()),2*v2it.step()); 
    } else if (x.real() == float(0)) {
      float xi = x.imag();
      cblas_saxpy(size,xi,v1it.GetP(),v1it.step(),
	  reinterpret_cast<float*>(v2it.GetP())+1,2*v2it.step()); 
    } else NonBlasAddVV(x,v1it,v2it,size);
  }
  template <StepItType S2, StepItType S3> inline void DoAddVV(
      const float x, CVIt<float,S2,NonConj> v1it,
      VIt<complex<float>,S3,NonConj> v2it, size_t size)
  {
    TMVAssert(size>0);
    TMVAssert(x != float(0));
    TMVAssert(v1it.step() > 0);
    TMVAssert(v2it.step() > 0);
    cblas_saxpy(size,x,v1it.GetP(),v1it.step(),
	reinterpret_cast<float*>(v2it.GetP()),2*v2it.step()); 
  }
#endif
#endif

  //
  // MultVV
  //
  template <class T, StepItType S1, class T2, StepItType S2, ConjItType C2> 
    inline void NonBlasMultVV(T& res, CVIt<T,S1,NonConj> v1it, 
	CVIt<T2,S2,C2> v2it, size_t size)
    {
      TMVAssert(size>0);
      res=T(0);
      const CVIt<T,S1,NonConj> _end = v1it + size;
      for(;v1it!=_end;++v1it,++v2it) res += (*v1it) * (*v2it);
    }
  template <class T, StepItType S1, class T2, StepItType S2, ConjItType C2> 
    inline void DoMultVV(T& res, CVIt<T,S1,NonConj> v1it, CVIt<T2,S2,C2> v2it,
	size_t size)
    { NonBlasMultVV(res,v1it,v2it,size); }
#ifdef BLAS
  template <StepItType S1, StepItType S2> inline void DoMultVV(
      double& res, CVIt<double,S1,NonConj> v1it,
      CVIt<double,S2,NonConj> v2it, size_t size)
  { 
    TMVAssert(size>0);
    TMVAssert(v1it.step()>0);
    TMVAssert(v2it.step()>0);
    res = cblas_ddot(size,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step()); 
  }
  template <StepItType S1, StepItType S2> inline void DoMultVV(
      complex<double>& res, CVIt<double,S1,NonConj> v1it,
      CVIt<double,S2,NonConj> v2it, size_t size)
  { 
    TMVAssert(size>0);
    TMVAssert(v1it.step()>0);
    TMVAssert(v2it.step()>0);
    res = cblas_ddot(size,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step()); 
  }
  template <StepItType S1, StepItType S2> inline void DoMultVV(
      complex<double>& res, CVIt<complex<double>,S1,NonConj> v1it, 
      CVIt<complex<double>,S2,NonConj> v2it, size_t size)
  { 
    TMVAssert(size>0);
    TMVAssert(v1it.step()>0);
    TMVAssert(v2it.step()>0);
    cblas_zdotu_sub(size,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step(),&res);
  }
  template <StepItType S1, StepItType S2> inline void DoMultVV(
      complex<double>& res, CVIt<complex<double>,S1,NonConj> v1it, 
      CVIt<complex<double>,S2,Conj> v2it, size_t size)
  { 
    TMVAssert(size>0);
    TMVAssert(v1it.step()>0);
    TMVAssert(v2it.step()>0);
    cblas_zdotc_sub(size,v2it.GetP(),v2it.step(),v1it.GetP(),v1it.step(),&res);
  }
#ifndef NOFLOAT
  template <StepItType S1, StepItType S2> inline void DoMultVV(
      float& res, CVIt<float,S1,NonConj> v1it,
      CVIt<float,S2,NonConj> v2it, size_t size)
  { 
    TMVAssert(size>0);
    TMVAssert(v1it.step()>0);
    TMVAssert(v2it.step()>0);
    res = cblas_sdot(size,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step()); 
  }
  template <StepItType S1, StepItType S2> inline void DoMultVV(
      complex<float>& res, CVIt<float,S1,NonConj> v1it,
      CVIt<float,S2,NonConj> v2it, size_t size)
  { 
    TMVAssert(size>0);
    TMVAssert(v1it.step()>0);
    TMVAssert(v2it.step()>0);
    res = cblas_sdot(size,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step()); 
  }
  template <StepItType S1, StepItType S2> inline void DoMultVV(
      complex<float>& res, CVIt<complex<float>,S1,NonConj> v1it, 
      CVIt<complex<float>,S2,NonConj> v2it, size_t size)
  { 
    TMVAssert(size>0);
    TMVAssert(v1it.step()>0);
    TMVAssert(v2it.step()>0);
    cblas_cdotu_sub(size,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step(),&res);
  }
  template <StepItType S1, StepItType S2> inline void DoMultVV(
      complex<float>& res, CVIt<complex<float>,S1,NonConj> v1it, 
      CVIt<complex<float>,S2,Conj> v2it, size_t size)
  { 
    TMVAssert(size>0);
    TMVAssert(v1it.step()>0);
    TMVAssert(v2it.step()>0);
    cblas_cdotc_sub(size,v2it.GetP(),v2it.step(),v1it.GetP(),v1it.step(),&res);
  }
#endif
#endif // BLAS

  // Remove Conj on v1
  template <class T, StepItType S1, class T2, StepItType S2> 
    inline void DoMultVV(
	T& res, CVIt<T,S1,Conj> v1it, CVIt<T2,S2,Conj> v2it, size_t size)
    { 
      DoMultVV(res,CVIt<T,S1,NonConj>(v1it.GetP(),v1it.step()),
	  CVIt<T2,S2,NonConj>(v2it.GetP(),v2it.step()),size);
      res = CONJ(res);
    }
  template <class T, StepItType S1, StepItType S2> inline void DoMultVV(
      T& res, CVIt<T,S1,Conj> v1it, CVIt<T,S2,NonConj> v2it, size_t size)
  { DoMultVV(res,v2it,v1it,size); }
  template <class T, StepItType S1, StepItType S2> inline void DoMultVV(
      complex<T>& res, CVIt<complex<T>,S1,Conj> v1it,
      CVIt<T,S2,NonConj> v2it, size_t size)
  { 
    DoMultVV(res,CVIt<complex<T>,S1,NonConj>(v1it.GetP(),v1it.step()),
	v2it,size);
    res = CONJ(res);
  }

  // real * complex
  template <class T, StepItType S1, ConjItType C1, StepItType S2, ConjItType C2>
    inline void DoMultVV(complex<T>& res, CVIt<T,S1,C1> v1it,
	CVIt<complex<T>,S2,C2> v2it, size_t size)
    { DoMultVV(res,v2it,v1it,size); }

  // real * real -> complex
  template <class T, StepItType S1, ConjItType C1, StepItType S2, ConjItType C2>
    inline void DoMultVV(complex<T>& res, CVIt<T,S1,C1> v1it,
	CVIt<T,S2,C2> v2it, size_t size)
    { T x; DoMultVV(x,v1it,v1it,size); res = x; }

  template <class T1, class T2, StepItType S2, ConjItType C2, class T3, StepItType S3, ConjItType C3, class T4, StepItType S4>
    inline void DoAddElementProd(const T1 alpha, CVIt<T2,S2,C2> xit,
	CVIt<T3,S3,C3> yit, VIt<T4,S4,NonConj> zit, size_t size)
    // zi += alpha * xi * yi 
    {
      TMVAssert(alpha != T1(0));
      TMVAssert(size>0);
      const VIt<T4,S4,NonConj> _end = zit + size;

      if (alpha == T1(1))
	for(;zit!=_end;++xit,++yit,++zit) *zit += (*xit) * (*yit);
      else if (alpha == T1(-1))
	for(;zit!=_end;++xit,++yit,++zit) *zit -= (*xit) * (*yit);
      else 
	for(;zit!=_end;++xit,++yit,++zit) *zit += alpha * (*xit) * (*yit);
    }

  template <class T1, class T2, StepItType S2, ConjItType C2, class T3, StepItType S3>
    inline void DoElementProd(const T1 alpha, CVIt<T2,S2,C2> xit, VIt<T3,S3,NonConj> yit, size_t size)
    // yi = xi * yi
    {
      TMVAssert(size>0);
      const VIt<T3,S3,NonConj> _end = yit + size;

      if (alpha == T1(1))
	for(;yit!=_end;++xit,++yit) *yit *= (*xit);
      else if (alpha == T1(-1))
	for(;yit!=_end;++xit,++yit) *yit *= -(*xit);
      else 
	for(;yit!=_end;++xit,++yit) *yit *= alpha * (*xit);
    }

} // namespace mv


