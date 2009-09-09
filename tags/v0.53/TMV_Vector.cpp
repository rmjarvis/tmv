
#include "TMV_Vector.h"
#include "TMV_VectorArith.h"

namespace tmv {

  //
  // Access
  //
 
  template <class T> T GenVector<T>::cref(size_t i) const
  {
    TMVAssert(i < size());
    const T* vi;
    if (step() == 1) vi = cptr() + i;
    else if (step() == 0) vi = cptr();
    else vi = cptr() + int(i)*step();
    return (itsct == Conj) ? CONJ(*vi) : *vi;
  }

  template <class T, IndexStyle I> RefType(T) VectorView<T,I>::ref(size_t i) const
  {
    TMVAssert(i < size());
    T* vi;
    if (step() == 1) vi = ptr() + i;
    else if (step() == 0) vi = ptr();
    else vi = ptr() + int(i)*step();
#ifdef TMVFLDEBUG
    TMVAssert(vi >= first);
    TMVAssert(vi < last);
#endif
    return REF(vi,ct());
  }

  //
  // OKSubVector
  //
  template <class T> bool GenVector<T>::OKSubVector(
      int i1, int i2, int istep) const
  {
    if (i1==i2) return true;  // no elements
    bool ok = true;
    if (istep == 0) {
      ok = false;
      cout<<"istep ("<<istep<<") cannot be 0\n";
    }
    if (i1 < 0 || i1 >= int(size())) {
      ok = false;
      cout<<"first element ("<<i1<<") must be in 0 -- "<<size()-1<<endl;
    }
    if (i2-istep < 0 || i2-istep >= int(size())) {
      ok = false;
      cout<<"last element ("<<i2-istep<<") must be in 0 -- "<<size()-1<<endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      cout<<"range ("<<i2-i1<<") must be multiple of istep ("<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      cout<<"n elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
    }
    return ok;
  }

  template <class T> bool ConstVectorView<T,FortranStyle>::OKSubVector(
      int i1, int i2, int istep) const
  {
    if (i1==i2) return true;  // no elements
    bool ok = true;
    if (istep == 0) {
      ok = false;
      cout<<"istep ("<<istep<<") cannot be 0\n";
    }
    if (i1 < 1 || i1 > int(this->size())) {
      ok = false;
      cout<<"first element ("<<i1<<") must be in 1 -- "<<this->size()<<endl;
    }
    if (i2 < 1 || i2 > int(this->size())) {
      ok = false;
      cout<<"last element ("<<i2<<") must be in 1 -- "<<this->size()<<endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      cout<<"range ("<<i2-i1<<") must be multiple of istep ("<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      cout<<"n elements ("<<(i2-i1)/istep+1<<") must be positive\n";
    }
    return ok;
  }

  //
  // Norm2
  //
  template <class T> RealType(T) GenVector<T>::NormSq() const
  {
    const int s = step();
    if (s < 0) return Reverse().NormSq();
    else if (IsComplex(T()) && s == 1) return Flatten().NormSq();
    else if (s == 0) return RealType(T)(size()) * NORM(*cptr());
    else {
      RealType(T) sum(0);
      const T* p = cptr();

      for(size_t i=size();i>0;--i,p+=s) sum += NORM(*p);
      return sum;
    }
  }
  template <class T> inline RealType(T) NonBlasNorm2(const GenVector<T>& v)
  { return tmv::SQRT(NormSq(v)); }
#ifdef BLAS
  template <class T> inline RealType(T) BlasNorm2(const GenVector<T>& v)
  { return NonBlasNorm2(v); }
  template <> inline double BlasNorm2(const GenVector<double>& v)
  {
    return cblas_dnrm2(v.size(),v.cptr(),v.step()); 
  }
  template <> inline double BlasNorm2(const GenVector<complex<double> >& v)
  { 
    return cblas_dznrm2(v.size(),v.cptr(),v.step()); 
  }
#ifndef NOFLOAT
  template <> inline float BlasNorm2(const GenVector<float>& v)
  {
    return cblas_snrm2(v.size(),v.cptr(),v.step()); 
  }
  template <> inline float BlasNorm2(const GenVector<complex<float> >& v)
  { 
    return cblas_scnrm2(v.size(),v.cptr(),v.step()); 
  }
#endif
#endif
  template <class T> RealType(T) GenVector<T>::Norm2() const
  { 
    if (step() < 0)
      return Reverse().Norm2();
    else 
#ifdef BLAS
      return BlasNorm2(*this); 
#else
      return NonBlasNorm2(*this); 
#endif
  }

  //
  // SumElements
  //
  template <class T> T GenVector<T>::SumElements() const
  {
    const int s = step();
    if (s < 0) return Reverse().SumElements();
    else if (s == 0) return RealType(T)(size()) * (*cptr());
    else {
      T sum(0);
      const T* p = cptr();

      if (s==1)
	for(size_t i=size();i>0;--i,++p) sum += *p;
      else
	for(size_t i=size();i>0;--i,p+=s) sum += *p;

      return isconj() ? CONJ(sum) : sum;
    }
  }

  template <class T> RealType(T) DoSumAbsElements(const GenVector<T>& v)
  {
    TMVAssert(v.step() >= 0)
    const int s = v.step();
    RealType(T) sum(0);
    const T* p = v.cptr();
    if (s == 1) 
      for(size_t i=v.size();i>0;--i,++p) sum += abs(*p);
    else if (s == 0)
      sum = RealType(T)(v.size()) * abs(*p);
    else 
      for(size_t i=v.size();i>0;--i,p+=s) sum += abs(*p);
    return sum;
  }
#ifdef BLAS
  template <> inline double DoSumAbsElements(const GenVector<double>& v)
  { 
    TMVAssert(v.step()>=0);
    return cblas_dasum(v.size(),v.cptr(),v.step()); 
  }
#ifndef NOFLOAT
  template <> inline float DoSumAbsElements(const GenVector<float>& v)
  { 
    TMVAssert(v.step()>=0);
    return cblas_sasum(v.size(),v.cptr(),v.step()); 
  }
#endif
#ifdef LAP
  template <> inline double DoSumAbsElements(
      const GenVector<complex<double> >& v)
  { 
    TMVAssert(v.step()>=0);
    int n = v.size();
    int s = v.step();
    return dzsum1(&n,LAP_Complex(v.cptr()),&s); 
  }
#ifndef NOFLOAT
  template <> inline float DoSumAbsElements(
      const GenVector<complex<float> >& v)
  { 
    TMVAssert(v.step()>=0);
    int n = v.size();
    int s = v.step();
    return scsum1(&n,LAP_Complex(v.cptr()),&s); 
  }
#endif
#endif // LAP
#endif // BLAS

  template <class T> RealType(T) GenVector<T>::SumAbsElements() const
  {
    if (step() > 0) return DoSumAbsElements(*this); 
    else if (step() == 0) return size() * abs(*cptr());
    else return DoSumAbsElements(Reverse()); 
  }


  //
  // Find Min/Max Element
  //

  template <class T> T FindMinElement(const GenVector<T>& v, size_t& imin)
  {
    TMVAssert(v.size() > 0);
    TMVAssert(v.step() > 0);

    const T* p = v.cptr();
    const int s = v.step();
    T min = *p;
    imin = 0;
    size_t i=1;
    if (s == 1) {
      ++p;
      for(size_t k=v.size()-1;k>0; --k,++p,++i) {
	if (tmv::REAL(*p) < tmv::REAL(min)) {
	  min = *p;
	  imin = i;
	}
      }
    }
    else {
      p += s;
      for(size_t k=v.size()-1;k>0; --k,p+=s,++i) {
	if (tmv::REAL(*p) < tmv::REAL(min)) {
	  min = *p;
	  imin = i;
	}
      }
    }
    return v.isconj() ? CONJ(min) : min;
  }
  template <class T> T FindMaxElement(const GenVector<T>& v, size_t& imax)
  {
    TMVAssert(v.size() > 0);
    TMVAssert(v.step() > 0);

    const T* p = v.cptr();
    const int s = v.step();
    T max = *p;
    imax = 0;
    size_t i=1;
    if (s == 1) {
      ++p;
      for(size_t k=v.size()-1;k>0; --k,++p,++i) {
	if (tmv::REAL(*p) > tmv::REAL(max)) {
	  max = *p;
	  imax = i;
	}
      }
    }
    else {
      p += s;
      for(size_t k=v.size()-1;k>0; --k,p+=s,++i) {
	if (tmv::REAL(*p) > tmv::REAL(max)) {
	  max = *p;
	  imax = i;
	}
      }
    }
    return v.isconj() ? CONJ(max) : max;
  }
  template <class T> RealType(T) FindMaxAbsElement(
      const GenVector<T>& v, size_t& imax)
  {
    TMVAssert(v.size() > 0);
    TMVAssert(v.step() > 0);

    const T* p = v.cptr();
    const int s = v.step();
    RealType(T) max = IsReal(T()) ? abs(*p) : NORM(*p);
    imax = 0;
    size_t i=1;
    if (s == 1) {
      ++p;
      for(size_t k=v.size()-1;k>0; --k,++p,++i) {
	RealType(T) absval = IsReal(T()) ? abs(*p) : NORM(*p);
	if (absval > max) { 
	  max = absval; 
	  imax = i;
	}
      }
    }
    else {
      p += s;
      for(size_t k=v.size()-1;k>0; --k,p+=s,++i) {
	RealType(T) absval = IsReal(T()) ? abs(*p) : NORM(*p);
	if (absval > max) { 
	  max = absval; 
	  imax = i;
	}
      }
    }
    return IsReal(T()) ? max : SQRT(max);
  }
  template <class T> RealType(T) FindMinAbsElement(
      const GenVector<T>& v, size_t& imin)
  {
    TMVAssert(v.size() > 0);
    TMVAssert(v.step() > 0);

    const T* p = v.cptr();
    const int s = v.step();
    RealType(T) min = IsReal(T()) ? abs(*p) : NORM(*p);
    imin = 0;
    size_t i=1;
    if (s == 1) {
      ++p;
      for(size_t k=v.size()-1;k>0; --k,++p,++i) {
	RealType(T) absval = IsReal(T()) ? abs(*p) : NORM(*p);
	if (absval < min) {
	  min = absval; 
	  imin = i;
	}
      }
    } else {
      p += s;
      for(size_t k=v.size()-1;k>0; --k,p+=s,++i) {
	RealType(T) absval = IsReal(T()) ? abs(*p) : NORM(*p);
	if (absval < min) {
	  min = absval; 
	  imin = i;
	}
      }
    }
    return IsReal(T()) ? min : SQRT(min);
  }
#ifdef BLAS
  template <> double FindMaxAbsElement(
      const GenVector<double>& v, size_t& imax)
  {
    imax = cblas_idamax(v.size(),v.cptr(),v.step()); 
    TMVAssert(imax < v.size());
    return abs(v[imax]);
  }
#ifdef MKL
  // (This is an extension that MKL includes, but not ATLAS)
  template <> double FindMinAbsElement(
      const GenVector<double>& v, size_t& imin)
  {
    imin = cblas_idamin(v.size(),v.cptr(),v.step()); 
    TMVAssert(imin < v.size());
    return abs(v[imin]);
  }
#endif
#ifndef NOFLOAT
  template <> float FindMaxAbsElement(
      const GenVector<float>& v, size_t& imax)
  {
    imax = cblas_isamax(v.size(),v.cptr(),v.step()); 
    TMVAssert(imax < v.size());
    return abs(v[imax]);
  }
#ifdef MKL
  // (This is an extension that MKL includes, but not ATLAS)
  template <> float FindMinAbsElement(
      const GenVector<float>& v, size_t& imin)
  {
    TMVAssert(v.ct() == NonConj);
    imin = cblas_isamin(v.size(),v.cptr(),v.step()); 
    TMVAssert(imin < v.size());
    return abs(v[imin]);
  }
#endif
#endif
#endif

  template <class T> T GenVector<T>::MinElement(size_t* iminout) const
  {
    TMVAssert(size()>0);
    if (step() < 0) {
      T min = Reverse().MinElement(iminout);
      if (iminout) *iminout = size()-1-(*iminout);
      return min;
    } else if (step() == 0) {
      if (iminout) *iminout = 0;
      return *cptr();
    } else {
      size_t imin;
      T min = FindMinElement(*this,imin);
      TMVAssert(imin < size());
      if (iminout) *iminout = imin;
      return min;
    }
  }
  template <class T> T GenVector<T>::MaxElement(size_t* imaxout) const
  {
    TMVAssert(size()>0);
    if (step() < 0) {
      T max = Reverse().MaxElement(imaxout);
      if (imaxout) *imaxout = size()-1-(*imaxout);
      return max;
    } else if (step() == 0) {
      if (imaxout) *imaxout = 0;
      return *cptr();
    } else {
      size_t imax;
      T max = FindMaxElement(*this,imax);
      TMVAssert(imax < size());
      if (imaxout) *imaxout = imax;
      return max;
    }
  }

  template <class T> RealType(T) GenVector<T>::MinAbsElement(
      size_t* iminout) const
  {
    TMVAssert(size()>0);
    if (step() < 0) {
      RealType(T) min = Reverse().MinAbsElement(iminout);
      if (iminout) *iminout = size()-1-(*iminout);
      return min;
    } else if (step() == 0) {
      if (iminout) *iminout = 0;
      return abs(*cptr());
    } else {
      size_t imin;
      RealType(T) min = FindMinAbsElement(*this,imin);
      TMVAssert(imin < size());
      if (iminout) *iminout = imin;
      return min;
    }
  }
  template <class T> RealType(T) GenVector<T>::MaxAbsElement( 
      size_t* imaxout) const
  {
    TMVAssert(size()>0);
    if (step() < 0) {
      RealType(T) max = Reverse().MaxAbsElement(imaxout);
      if (imaxout) *imaxout = size()-1-(*imaxout);
      return max;
    } else if (step() == 0) {
      if (imaxout) *imaxout = 0;
      return abs(*cptr());
    } else {
      size_t imax;
      RealType(T) max = FindMaxAbsElement(*this,imax);
      TMVAssert(imax < size());
      if (imaxout) *imaxout = imax;
      return max;
    }
  }

  //
  // Other Modifying Functions:
  //   Clip
  //   SetAllTo
  //   AddToAll
  //   ConjugateSelf
  //   ReverseSelf
  //   Sort
  //

  template <class T, IndexStyle I> const VectorView<T,I>& VectorView<T,I>::Clip(
      RealType(T) thresh) const
  {
    const int s = step();
    if (s < 0) Reverse().Clip(thresh);
    else if (s == 0) {
      if (IsReal(T())) {
	if (abs(*ptr()) < thresh) *ptr() = T(0);
      } else {
	RealType(T) threshsq = thresh*thresh;
	if (NORM(*ptr()) < threshsq) *ptr() = T(0);
      }
    } else {
      T* p = ptr();

      if (IsReal(T())) {
	if (s == 1) {
	  for(size_t i=size();i>0;--i,++p) 
	    if (abs(*p) < thresh) *p = T(0);
	} else {
	  for(size_t i=size();i>0;--i,p+=s) 
	    if (abs(*p) < thresh) *p = T(0);
	}
      } else {
	RealType(T) threshsq = thresh*thresh;
	if (s == 1) {
	  for(size_t i=size();i>0;--i,++p) 
	    if (NORM(*p) < threshsq) *p = T(0);
	} else {
	  for(size_t i=size();i>0;--i,p+=s) 
	    if (NORM(*p) < threshsq) *p = T(0);
	}
      }
    }
    return *this; 
  }

  template <class T, IndexStyle I> const VectorView<T,I>& VectorView<T,I>::SetAllTo(T x) const
  {
    const int s = step();
    if (s < 0) Reverse().SetAllTo(x);
    else if (s == 0) *ptr() = x;
    else {
      T* p = ptr();

      if (isconj()) 
	if (s == 1) 
	  for(size_t i=size();i>0;--i,++p) *p = CONJ(x); 
	else 
	  for(size_t i=size();i>0;--i,p+=s) *p = CONJ(x); 
      else
	if (s == 1) 
	  for(size_t i=size();i>0;--i,++p) *p = x; 
	else 
	  for(size_t i=size();i>0;--i,p+=s) *p = x; 
    }
    return *this; 
  }

  template <class T, IndexStyle I> const VectorView<T,I>& VectorView<T,I>::AddToAll(T x) const
  {
    const int s = step();
    if (s < 0) Reverse().AddToAll(x);
    else if (s == 0) *ptr() += x;
    else {
      T* p = ptr();

      if (isconj()) 
	if (s == 1) 
	  for(size_t i=size();i>0;--i,++p) *p += CONJ(x); 
	else 
	  for(size_t i=size();i>0;--i,p+=s) *p += CONJ(x); 
      else 
	if (s == 1) 
	  for(size_t i=size();i>0;--i,++p) *p += x; 
	else 
	  for(size_t i=size();i>0;--i,p+=s) *p += x; 
    }
    return *this; 
  }

  template <class T> void NonLapConjugate(const VectorView<complex<T>,CStyle>& v)
  {
    TMVAssert(v.step() > 0);

    T* p = v.Imag().ptr();

    if (v.step() == 1)
      for(size_t i=v.size();i>0;--i,p+=2) *p = -(*p);
    else if (v.step() == 0) *p = -*p;
    else {
      const int s = 2*v.step();
      for(size_t i=v.size();i>0;--i,p+=s) *p = -(*p);
    }
  }
  template <class T> inline void NonLapConjugate(const VectorView<T,CStyle>& v) {}
#ifdef LAP
  template <class T> inline void LapConjugate(const VectorView<T,CStyle>& v)
  { NonLapConjugate(v); }
  template <> inline void LapConjugate(const VectorView<complex<double>,CStyle>& v)
  { 
    int n = v.size();
    int s = v.step();
    zlacgv(&n,LAP_Complex(v.ptr()),&s); 
  }
#ifndef NOFLOAT
  template <> inline void LapConjugate(const VectorView<complex<float>,CStyle>& v)
  {
    int n = v.size();
    int s = v.step();
    clacgv(&n,LAP_Complex(v.ptr()),&s); 
  }
#endif
#endif
  template <class T, IndexStyle I> 
    const VectorView<T,I>& VectorView<T,I>::ConjugateSelf() const
    {
      if (step() < 0) Reverse().ConjugateSelf();
      else {
#ifdef LAP
	LapConjugate(*this);
#else
	NonLapConjugate(*this);
#endif
      }
      return *this; 
    }

  template <class T, IndexStyle I> const VectorView<T,I>& VectorView<T,I>::Permute(
      const size_t* p, size_t i1, size_t i2) const
  { 
    TMVAssert(i2 <= size());
    TMVAssert(i1 <= i2);
    for(size_t i=i1;i<i2;++i) Swap(i,p[i]);
    return *this; 
  }

  template <class T, IndexStyle I> 
    const VectorView<T,I>& VectorView<T,I>::ReversePermute(
	const size_t* p, size_t i1, size_t i2) const
    { 
      TMVAssert(i2 <= size());
      TMVAssert(i1 <= i2);
      for(size_t i=i2;i>i1;) { --i; Swap(i,p[i]); }
      return *this; 
    }

  template <class T, IndexStyle I> 
    const VectorView<T,I>& VectorView<T,I>::ReverseSelf() const
    {
      const int s = step();
      if (s < 0) Reverse().ReverseSelf();
      else if (s > 0) {
	T* p1 = ptr();
	if (s == 1) 
	  for(T* p2=p1+size()-1;p2>p1;++p1,--p2) swap(*p1,*p2);
	else
	  for(T* p2=p1+s*(size()-1);p2>p1;p1+=s,p2-=s) swap(*p1,*p2);
      }
      return *this;
    }

  template <class T> class VTIndex
  {

    private :

      RealType(T) itsvalue;
      size_t itsi;

    public :

      VTIndex() : itsvalue(RealType(T)(0)), itsi(0) {}

      VTIndex(T val, size_t i, ADType ad, COMPType comp) : itsi(i)
      {
	bool neg = ad==DESCEND;
	switch(comp) {
	  case REAL_COMP : itsvalue = neg ? -REAL(val) : REAL(val); break;
	  case ABS_COMP : itsvalue = neg ? -abs(val) : abs(val); break;
	  case IMAG_COMP : itsvalue = neg ? -IMAG(val) : IMAG(val); break;
	  case ARG_COMP : itsvalue = neg ? -ARG(val) : ARG(val); break;
	  default : TMVAssert(false);
	}
      }

      // Use default copy, op=, destructor

      size_t GetI() const { return itsi; }
      size_t GetVal() const { return itsvalue; }
      bool operator<(const VTIndex& rhs) const
      { return itsvalue < rhs.itsvalue; }

  };

  template <class T> class Compare
  {
    public:

      Compare(ADType _ad, COMPType _comp) : ad(_ad), comp(_comp) {}

      bool operator()(T x, T y) const {
	if (ad == ASCEND) {
	  switch(comp) {
	    case REAL_COMP : return REAL(x) < REAL(y);
	    case ABS_COMP : return abs(x) < abs(y);
	    case IMAG_COMP : return IMAG(x) < IMAG(y);
	    case ARG_COMP : return ARG(x) < ARG(y);
	    default : TMVAssert(false); return false;
	  }
	} else {
	  switch(comp) {
	    case REAL_COMP : return REAL(x) > REAL(y);
	    case ABS_COMP : return abs(x) > abs(y);
	    case IMAG_COMP : return IMAG(x) > IMAG(y);
	    case ARG_COMP : return ARG(x) > ARG(y);
	    default : TMVAssert(false); return false;
	  }
	} 
      }

    private:

      const ADType ad;
      const COMPType comp;
  };


  template <class T, IndexStyle I> const VectorView<T,I>& VectorView<T,I>::Sort(
      size_t* P, ADType ad, COMPType comp) const
  {
    if (P) {
      vector<VTIndex<T> > newindex(size());
      for(size_t i=0;i<size();++i) {
	newindex[i] = VTIndex<T>(ref(i),i,ad,comp);
      }
      std::sort(newindex.begin(),newindex.end());

      // newindex[i]=j means value at original j location needs to go to i.
      vector<size_t> currindex(size());
      vector<size_t> origindex(size());
      for(size_t i=0;i<size();++i) {
	currindex[i] = i;
	origindex[i] = i;
      }
      // currindex[i]=j means value at original i location is currently at j.
      // origindex[j]=i means value at original i location is currently at j.
      for(size_t i=0;i<size();++i) {
	size_t ip = currindex[newindex[i].GetI()];
	P[i] = ip;
	if (i != ip) {
	  Swap(i,ip);
	  size_t origi = origindex[i];
	  size_t origip = origindex[ip];
	  currindex[origi] = ip;
	  currindex[origip] = i;
	  origindex[i] = origip;
	  origindex[ip] = origi;
	}
      }
    } else {
      const Compare<T> cc(ad,comp);
      std::sort(begin(),end(),cc);
    }
    return *this;
  }

  //
  // Copy Vectors
  //

  template <bool c1, class T> void NonBlasCopy(
      const GenVector<T>& v1, const VectorView<T>& v2)
  {
    const T* v1ptr = v1.cptr();
    T* v2ptr = v2.ptr();
    const int step1 = v1.step();
    const int step2 = v2.step();

    if (step1 == 1 && step2 == 1)
      if (c1)
	for(size_t i=v2.size();i>0;--i,++v1ptr,++v2ptr)
	  *v2ptr = CONJ(*v1ptr);
      else
	memmove(v2ptr,v1ptr,v2.size()*sizeof(T)); 
    else
      for(size_t i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2)
	*v2ptr = (c1 ? CONJ(*v1ptr) : (*v1ptr));
  }
#ifdef BLAS
  template <class T> inline void BlasCopy(
      const GenVector<T>& v1, const VectorView<T>& v2)
  { NonBlasCopy<false>(v1,v2); }
  template <> inline void BlasCopy(
      const GenVector<double>& v1, const VectorView<double>& v2)
  { cblas_dcopy(v2.size(),v1.cptr(),v1.step(),v2.ptr(),v2.step()); }
  template <> inline void BlasCopy(
      const GenVector<complex<double> >& v1,
      const VectorView<complex<double> >& v2)
  { cblas_zcopy(v2.size(),v1.cptr(),v1.step(),v2.ptr(),v2.step()); }
#ifndef NOFLOAT
  template <> inline void BlasCopy(
      const GenVector<float>& v1, const VectorView<float>& v2)
  { cblas_scopy(v2.size(),v1.cptr(),v1.step(),v2.ptr(),v2.step()); }
  template <> inline void BlasCopy(
      const GenVector<complex<float> >& v1,
      const VectorView<complex<float> >& v2)
  { cblas_ccopy(v2.size(),v1.cptr(),v1.step(),v2.ptr(),v2.step()); }
#endif
#endif
  template <bool c1, class T> void DoCopySameType(
      const GenVector<T>& v1, const VectorView<T>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v2.size()>0);
    TMVAssert(v1.step()!=0);
    TMVAssert(v2.step()!=0);
    TMVAssert(v2.ct()==NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(!v2.SameAs(v1));
    TMVAssert(c1 == v1.isconj());

#ifdef BLAS
    if (!c1 && v1.step() > 0 && v2.step() > 0) BlasCopy(v1,v2);
    else 
#endif
      NonBlasCopy<c1>(v1,v2); 
  }

  //
  // Swap Vectors
  //

  template <class T> inline void conjswap(T& x, T& y)
  {
    T temp = x;
    x = CONJ(y);
    y = CONJ(temp);
  }

  template <class T> void NonBlasSwap(
      const VectorView<T>& v1, const VectorView<T>& v2)
  {
    T* v1ptr = v1.ptr();
    T* v2ptr = v2.ptr();
    const int step1 = v1.step();
    const int step2 = v2.step();

    if (v1.isconj())
      if (step1 == 1 && step2 == 1)
	for(size_t i=v2.size();i>0;--i,++v1ptr,++v2ptr) 
	  conjswap(*v1ptr,*v2ptr); 
      else
	for(size_t i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2) 
	  conjswap(*v1ptr,*v2ptr); 
    else
      if (step1 == 1 && step2 == 1)
	for(size_t i=v2.size();i>0;--i,++v1ptr,++v2ptr) 
	  swap(*v1ptr,*v2ptr); 
      else
	for(size_t i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2) 
	  swap(*v1ptr,*v2ptr); 
  }
#ifdef BLAS
  template <class T> inline void BlasSwap(
      const VectorView<T>& v1, const VectorView<T>& v2)
  { NonBlasSwap(v1,v2); }
  template <> inline void BlasSwap(
      const VectorView<double>& v1, const VectorView<double>& v2)
  { cblas_dswap(v2.size(),v1.ptr(),v1.step(),v2.ptr(),v2.step()); }
  template <> inline void BlasSwap(
      const VectorView<complex<double> >& v1, 
      const VectorView<complex<double> >& v2)
  { cblas_zswap(v2.size(),v1.ptr(),v1.step(),v2.ptr(),v2.step()); }
#ifndef NOFLOAT
  template <> inline void BlasSwap(
      const VectorView<float>& v1, const VectorView<float>& v2)
  { cblas_sswap(v2.size(),v1.ptr(),v1.step(),v2.ptr(),v2.step()); }
  template <> inline void BlasSwap(
      const VectorView<complex<float> >& v1, 
      const VectorView<complex<float> >& v2)
  { cblas_cswap(v2.size(),v1.ptr(),v1.step(),v2.ptr(),v2.step()); }
#endif
#endif
  template <class T> inline void DoSwap(
      const VectorView<T>& v1, const VectorView<T>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v2.size()>0);
    TMVAssert(v1.step()!=0);
    TMVAssert(v2.step()!=0);
    TMVAssert(v2.ct()==NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(!v2.SameAs(v1));

#ifdef BLAS
    if (v1.step() > 0 && v2.step() > 0) BlasSwap(v1,v2);
    else 
#endif
      NonBlasSwap(v1,v2); 
  }

  template <class T> void Swap(
      const VectorView<T>& v1, const VectorView<T>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);
    TMVAssert(!v2.SameAs(v1) || v1.isconj() == v2.isconj());

    if (v1.size() > 0 && !v2.SameAs(v1)) {
      if (ShouldReverse(v1.step(),v2.step()))
	Swap(v1.Reverse(),v2.Reverse());
      else if (v2.isconj()) DoSwap(v1.Conjugate(),v2.Conjugate());
      else DoSwap(v1,v2);
    }
  }

  //
  // op ==
  //
  template <class T1, class T2> bool operator==(
      const GenVector<T1>& v1, const GenVector<T2>& v2)
  {
    TMVAssert(v1.size() == v2.size());

    if (v2.SameAs(v1)) return true;
    else {
      const T1* v1ptr = v1.cptr();
      const T2* v2ptr = v2.cptr();
      const int step1 = v1.step();
      const int step2 = v2.step();

      if (v1.isconj() == v2.isconj()) {
	if (step1 == 1 && step2 == 1) {
	  for(size_t i=v2.size();i>0;--i,++v1ptr,++v2ptr)
	    if ( *v1ptr != *v2ptr ) return false;
	} else {
	  for(size_t i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2)
	    if ( *v1ptr != *v2ptr ) return false;
	}
      } else {
	if (step1 == 1 && step2 == 1) {
	  for(size_t i=v2.size();i>0;--i,++v1ptr,++v2ptr)
	    if ( *v1ptr != CONJ(*v2ptr) ) return false;
	} else {
	  for(size_t i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2)
	    if ( *v1ptr != CONJ(*v2ptr) ) return false;
	}
      }
    }
    return true;
  }

  //
  // I/O
  //
  
  template <class T> void GenVector<T>::Write(ostream& fout) const
  {
    const T* p = cptr();
    const int s = step();

    fout << size() << " (";
    if (isconj()) 
      if (s == 1) 
	for(size_t i=size();i>0;--i,++p) fout << " " << CONJ(*p) << " ";
      else 
	for(size_t i=size();i>0;--i,p+=s) fout << " " << CONJ(*p) << " ";
    else 
      if (step() == 1) 
	for(size_t i=size();i>0;--i,++p) fout << " " << *p << " ";
      else
	for(size_t i=size();i>0;--i,p+=s) fout << " " << *p << " ";
    fout << " )";
  }

  template <class T> void GenVector<T>::Write(ostream& fout,
      RealType(T) thresh) const
  {
    const T* p = cptr();
    const int s = step();

    fout << size() << " (";
    if (isconj()) 
      if (s == 1) 
	for(size_t i=size();i>0;--i,++p) 
	  fout << " " << (abs(*p) < thresh ? T(0) : CONJ(*p)) << " ";
      else 
	for(size_t i=size();i>0;--i,p+=s) 
	  fout << " " << (abs(*p) < thresh ? T(0) : CONJ(*p)) << " ";
    else 
      if (step() == 1) 
	for(size_t i=size();i>0;--i,++p) 
	  fout << " " << (abs(*p) < thresh ? T(0) : *p) << " ";
      else
	for(size_t i=size();i>0;--i,p+=s) 
	  fout << " " << (abs(*p) < thresh ? T(0) : *p) << " ";
    fout << " )";
  }

  template <class T, IndexStyle I> void VectorView<T,I>::Read(istream& fin) const
  {
    T* p = ptr();
    const int s = step();

    char paren;
    fin >> paren;
    if (!fin || paren != '(') tmv_error("reading ( in Vector::Read");
    if (s == 1) 
      for(size_t i=size();i>0;--i,++p) fin >> *p;
    else
      for(size_t i=size();i>0;--i,++p) fin >> *p;
    if (!fin) tmv_error("reading value in Vector::Read");
    fin >> paren;
    if (paren != ')') tmv_error("reading ) in Vector::Read");
    if (isconj()) ConjugateSelf();
  }

  template <class T> istream& operator>>(istream& fin, const VectorView<T>& v)
  {
    size_t n;
    fin >> n;
    if (!fin) tmv_error("reading size in Vector::Read");
    if (n != v.size()) tmv_error("size does not match in Vector::Read");
    v.Read(fin);
    return fin;
  }

  template <class T, IndexStyle I> istream& operator>>(istream& fin, 
      Vector<T,I>*& v)
  {
    size_t n;
    fin >> n;
    if (!fin) tmv_error("reading size in Vector::Read");
    v = new Vector<T,I>(n);
    v->View().Read(fin);
    return fin;
  }


#define InstFile "TMV_Vector.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


