
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
    else vi = cptr() + int(i)*step();
    return (itsct == Conj) ? CONJ(*vi) : *vi;
  }

  template <class T> RefType(T) VectorView<T>::ref(size_t i) const
  {
    TMVAssert(i < size());
    T* vi;
    if (step() == 1) vi = ptr() + i;
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

  //
  // Norm2
  //
  template <class T> RealType(T) GenVector<T>::NormSq() const
  {
    if (isconj()) return Conjugate().NormSq();
    else {
      RealType(T) sum(0);
      if (step() == 1) {
	CVIt<T,Unit,NonConj> it = begin();
	for(size_t i=size();i>0;--i,++it) sum += tmv::NORM(*it);
      } else {
	CVIt<T,Step,NonConj> it = begin();
	for(size_t i=size();i>0;--i,++it) sum += tmv::NORM(*it);
      }
      return sum;
    }
  }
  template <class T> inline RealType(T) DoNorm2(const GenVector<T>& v)
  { return tmv::SQRT(NormSq(v)); }
#ifdef BLAS
  template <> inline double DoNorm2(const GenVector<double>& v)
  {
    if (v.step() > 0)
      return cblas_dnrm2(v.size(),v.cptr(),v.step()); 
    else {
      ConstVectorView<double> v2 = v.Reverse();
      return cblas_dnrm2(v2.size(),v2.cptr(),v2.step()); 
    }
  }
  template <> inline double DoNorm2(const GenVector<complex<double> >& v)
  { 
    if (v.step() > 0)
      return cblas_dznrm2(v.size(),v.cptr(),v.step()); 
    else {
      ConstVectorView<complex<double> > v2 = v.Reverse();
      return cblas_dznrm2(v2.size(),v2.cptr(),v2.step()); 
    }
  }
#ifndef NOFLOAT
  template <> inline float DoNorm2(const GenVector<float>& v)
  {
    if (v.step() > 0)
      return cblas_snrm2(v.size(),v.cptr(),v.step()); 
    else {
      ConstVectorView<float> v2 = v.Reverse();
      return cblas_snrm2(v2.size(),v2.cptr(),v2.step()); 
    }
  }
  template <> inline float DoNorm2(const GenVector<complex<float> >& v)
  { 
    if (v.step() > 0)
      return cblas_scnrm2(v.size(),v.cptr(),v.step()); 
    else {
      ConstVectorView<complex<float> > v2 = v.Reverse();
      return cblas_scnrm2(v2.size(),v2.cptr(),v2.step()); 
    }
  }
#endif
#endif
  template <class T> RealType(T) GenVector<T>::Norm2() const
  { return DoNorm2(*this); }

  //
  // SumElements
  //
  template <class T> T GenVector<T>::SumElements() const
  {
    if (isconj()) return CONJ(Conjugate().SumElements());
    else {
      T sum(0);
      if (step() == 1) {
	CVIt<T,Unit,NonConj> it = begin();
	for(size_t i=size();i>0;--i,++it) sum += *it;
      } else {
	CVIt<T,Step,NonConj> it = begin();
	for(size_t i=size();i>0;--i,++it) sum += *it;
      }
      return sum;
    }
  }

  template <class T> RealType(T) DoSumAbsElements(const GenVector<T>& v)
  {
    TMVAssert(v.step() > 0)
    if (v.isconj()) return DoSumAbsElements(v.Conjugate());
    else {
      RealType(T) sum(0);
      if (v.step() == 1) {
	CVIt<T,Unit,NonConj> it = v.begin();
	for(size_t i=v.size();i>0;--i,++it) sum += abs(*it);
      } else {
	CVIt<T,Step,NonConj> it = v.begin();
	for(size_t i=v.size();i>0;--i,++it) sum += abs(*it);
      }
      return sum;
    }
  }
#ifdef BLAS
  template <> inline double DoSumAbsElements(const GenVector<double>& v)
  { 
    TMVAssert(v.step()>0);
    return cblas_dasum(v.size(),v.cptr(),v.step()); 
  }
#ifndef NOFLOAT
  template <> inline float DoSumAbsElements(const GenVector<float>& v)
  { 
    TMVAssert(v.step()>0);
    return cblas_sasum(v.size(),v.cptr(),v.step()); 
  }
#endif
#ifdef LAP
  template <> inline double DoSumAbsElements(
      const GenVector<complex<double> >& v)
  { 
    TMVAssert(v.step()>0);
    int n = v.size();
    int s = v.step();
    return dzsum1(&n,LAP_Complex(v.cptr()),&s); 
  }
#ifndef NOFLOAT
  template <> inline float DoSumAbsElements(
      const GenVector<complex<float> >& v)
  { 
    TMVAssert(v.step()>0);
    int n = v.size();
    int s = v.step();
    return scsum1(&n,LAP_Complex(v.cptr()),&s); 
  }
#endif
#endif // LAP
#endif // BLAS

  template <class T> RealType(T) GenVector<T>::SumAbsElements() const
  {
    TMVAssert(step() != 0);
    if (step() > 0) return DoSumAbsElements(*this); 
    else return DoSumAbsElements(Reverse()); 
  }


  //
  // Find Min/Max Element
  //
  template <StepItType S, class T> T FindMinElement(
      CVIt<T,S,NonConj> it, size_t size, size_t& imin)
  {
    T min = *it;
    imin = 0;
    ++it;
    size_t i=1;
    for(size_t k=size-1;k>0; --k,++it,++i) {
      if (tmv::REAL(*it) < tmv::REAL(min)) {
	min = *it;
	imin = i;
      }
    }
    return min;
  }
  template <StepItType S, class T> T FindMaxElement(
      CVIt<T,S,NonConj> it, size_t size, size_t& imax)
  {
    T max = *it;
    imax = 0;
    ++it;
    size_t i=1;
    for(size_t k=size-1;k>0; --k,++it,++i) {
      if (tmv::REAL(*it) > tmv::REAL(max)) {
	max = *it;
	imax = i;
      }
    }
    return max;
  }
  template <StepItType S, class T> RealType(T) FindMaxAbsElement(
      CVIt<T,S,NonConj> it, size_t size, size_t& imax)
  {
    RealType(T) max = IsReal(T()) ? abs(*it) : NORM(*it);
    imax = 0;
    ++it;
    size_t i=1;
    for(size_t k=size-1;k>0; --k,++it,++i) {
      RealType(T) absval = IsReal(T()) ? abs(*it) : NORM(*it);
      if (absval > max) { 
	max = absval; 
	imax = i;
      }
    }
    return IsReal(T()) ? max : SQRT(max);
  }
  template <StepItType S, class T> RealType(T) FindMinAbsElement(
      CVIt<T,S,NonConj> it, size_t size, size_t& imin)
  {
    RealType(T) min = IsReal(T()) ? abs(*it) : NORM(*it);
    imin = 0;
    ++it;
    size_t i=1;
    for(size_t k=size-1;k>0; --k,++it,++i) {
      RealType(T) absval = IsReal(T()) ? abs(*it) : NORM(*it);
      if (absval < min) {
	min = absval; 
	imin = i;
      }
    }
    return IsReal(T()) ? min : SQRT(min);
  }
#ifdef BLAS
  template <StepItType S> double FindMaxAbsElement(
      CVIt<double,S,NonConj> v, size_t size, size_t& imax)
  {
    imax = cblas_idamax(size,v.GetP(),v.step()); 
    TMVAssert(imax < size);
    return abs(v[imax]);
  }
  /* (This is an extension that MKL includes, but not ATLAS)
  template <StepItType S> double FindMinAbsElement(
      CVIt<double,S,NonConj> v, size_t size, size_t& imin)
  {
    imin = cblas_idamin(size,v.GetP(),v.step()); 
    TMVAssert(imin < size);
    return abs(v[imin]);
  } */
#ifndef NOFLOAT
  template <StepItType S> float FindMaxAbsElement(
      CVIt<float,S,NonConj> v, size_t size, size_t& imax)
  {
    imax = cblas_isamax(size,v.GetP(),v.step()); 
    TMVAssert(imax < size);
    return abs(v[imax]);
  }
  /* (This is an extension that MKL includes, but not ATLAS)
  template <StepItType S> float FindMinAbsElement(
      CVIt<float,S,NonConj> v, size_t size, size_t& imin)
  {
    imin = cblas_isamin(size,v.GetP(),v.step()); 
    TMVAssert(imin < size);
    return abs(v[imin]);
  } */
#endif
#endif

  template <class T> T GenVector<T>::MinElement(size_t* iminout) const
  {
    TMVAssert(size()>0);
    if (isconj()) return CONJ(Conjugate().MinElement(iminout));
    else {
      size_t imin;
      T min;
      if (step() == 1)
	min = FindMinElement(CVIt<T,Unit,NonConj>(begin()),size(),imin);
      else
	min = FindMinElement(CVIt<T,Step,NonConj>(begin()),size(),imin);
      TMVAssert(imin < size());
      if (iminout) *iminout = imin;
      return min;
    }
  }
  template <class T> T GenVector<T>::MaxElement(size_t* imaxout) const
  {
    TMVAssert(size()>0);
    if (isconj()) return CONJ(Conjugate().MaxElement(imaxout));
    else {
      size_t imax;
      T max;
      if (step() == 1)
	max = FindMaxElement(CVIt<T,Unit,NonConj>(begin()),size(),imax);
      else
	max = FindMaxElement(CVIt<T,Step,NonConj>(begin()),size(),imax);
      TMVAssert(imax < size());
      if (imaxout) *imaxout = imax;
      return max;
    }
  }
  template <class T> RealType(T) GenVector<T>::MinAbsElement(
      size_t* iminout) const
  {
    TMVAssert(size()>0);
    if (isconj()) return Conjugate().MinAbsElement(iminout);
    else {
      size_t imin;
      RealType(T) min;
      if (step() == 1)
	min = FindMinAbsElement(CVIt<T,Unit,NonConj>(begin()),size(),imin);
#ifdef BLAS
      else if (step() < 0) {
	min = FindMinAbsElement(CVIt<T,Step,NonConj>(rbegin()),size(),imin);
	imin = size()-1-imin;
      }
#endif
      else
	min = FindMinAbsElement(CVIt<T,Step,NonConj>(begin()),size(),imin);
      TMVAssert(imin < size());
      if (iminout) *iminout = imin;
      return min;
    }
  }
  template <class T> RealType(T) GenVector<T>::MaxAbsElement( 
      size_t* imaxout) const
  {
    TMVAssert(size()>0);
    if (isconj()) return Conjugate().MaxAbsElement(imaxout);
    else {
      size_t imax;
      RealType(T) max;
      if (step() == 1)
	max = FindMaxAbsElement(CVIt<T,Unit,NonConj>(begin()),size(),imax);
#ifdef BLAS
      else if (step() < 0) {
	max = FindMaxAbsElement(CVIt<T,Step,NonConj>(rbegin()),size(),imax);
	imax = size()-1-imax;
      }
#endif
      else
	max = FindMaxAbsElement(CVIt<T,Step,NonConj>(begin()),size(),imax);
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

  template <class T> const VectorView<T>& VectorView<T>::Clip(
      RealType(T) thresh) const
  {
    if (IsReal(T())) {
      if (isconj()) {
	if (step() == 1) {
	  VIt<T,Unit,Conj> it = begin();
	  for(size_t i=size();i>0;--i,++it) 
	    if (abs(*it) < thresh) *it = T(0);
	} else {
	  VIt<T,Step,Conj> it = begin();
	  for(size_t i=size();i>0;--i,++it) 
	    if (abs(*it) < thresh) *it = T(0);
	}
      } else {
	if (step() == 1) {
	  VIt<T,Unit,NonConj> it = begin();
	  for(size_t i=size();i>0;--i,++it) 
	    if (abs(*it) < thresh) *it = T(0);
	} else {
	  VIt<T,Step,NonConj> it = begin();
	  for(size_t i=size();i>0;--i,++it) 
	    if (abs(*it) < thresh) *it = T(0);
	}
      }
    } else {
      RealType(T) threshsq = thresh*thresh;
      if (isconj()) {
	if (step() == 1) {
	  VIt<T,Unit,Conj> it = begin();
	  for(size_t i=size();i>0;--i,++it) 
	    if (NORM(*it) < threshsq) *it = T(0);
	} else {
	  VIt<T,Step,Conj> it = begin();
	  for(size_t i=size();i>0;--i,++it) 
	    if (NORM(*it) < threshsq) *it = T(0);
	}
      } else {
	if (step() == 1) {
	  VIt<T,Unit,NonConj> it = begin();
	  for(size_t i=size();i>0;--i,++it) 
	    if (NORM(*it) < threshsq) *it = T(0);
	} else {
	  VIt<T,Step,NonConj> it = begin();
	  for(size_t i=size();i>0;--i,++it) 
	    if (NORM(*it) < threshsq) *it = T(0);
	}
      }
    }
    return *this; 
  }

  template <class T> const VectorView<T>& VectorView<T>::SetAllTo(T x) const
  {
    if (isconj()) {
      if (step() == 1) std::fill(ptr(),ptr()+size(),CONJ(x));
      else {
	VIt<T,Step,Conj> it = begin();
	for(size_t i=size();i>0;--i,++it) *it = x; 
      }
    } else {
      if (step() == 1) std::fill(ptr(),ptr()+size(),x);
      else {
	VIt<T,Step,NonConj> it = begin();
	for(size_t i=size();i>0;--i,++it) *it = x; 
      }
    }
    return *this; 
  }

  template <class T> const VectorView<T>& VectorView<T>::AddToAll(T x) const
  {
    if (isconj()) {
      if (step() == 1) {
	VIt<T,Unit,Conj> it = begin();
	for(size_t i=size();i>0;--i,++it) *it += x; 
      } else {
	VIt<T,Step,Conj> it = begin();
	for(size_t i=size();i>0;--i,++it) *it += x; 
      }
    } else {
      if (step() == 1) {
	VIt<T,Unit,NonConj> it = begin();
	for(size_t i=size();i>0;--i,++it) *it += x; 
      } else {
	VIt<T,Step,NonConj> it = begin();
	for(size_t i=size();i>0;--i,++it) *it += x; 
      }
    }
    return *this; 
  }

  template <class T> void NonLapConjugate(const VectorView<complex<T> >& v)
  {
    if (v.isconj()) {
      if (v.step() == 1) {
	VIt<complex<T>,Unit,Conj> it = v.begin();
	for(size_t i=v.size();i>0;--i,++it) *it = tmv::CONJ(*it); 
      } else {
	VIt<complex<T>,Step,Conj> it = v.begin();
	for(size_t i=v.size();i>0;--i,++it) *it = tmv::CONJ(*it); 
      }
    } else {
      if (v.step() == 1) {
	VIt<complex<T>,Unit,NonConj> it = v.begin();
	for(size_t i=v.size();i>0;--i,++it) *it = tmv::CONJ(*it); 
      } else {
	VIt<complex<T>,Step,NonConj> it = v.begin();
	for(size_t i=v.size();i>0;--i,++it) *it = tmv::CONJ(*it); 
      }
    }
  }
  template <class T> inline void LapConjugate(const VectorView<T>& v) {}
  template <class T> inline void LapConjugate(const VectorView<complex<T> >& v)
  { NonLapConjugate(v); }
#ifdef LAP
  template <> inline void LapConjugate(const VectorView<complex<double> >& v)
  { 
    int n = v.size();
    int s = v.step();
    zlacgv(&n,LAP_Complex(v.ptr()),&s); 
  }
#ifndef NOFLOAT
  template <> inline void LapConjugate(const VectorView<complex<float> >& v)
  {
    int n = v.size();
    int s = v.step();
    clacgv(&n,LAP_Complex(v.ptr()),&s); 
  }
#endif
#endif
  template <class T> const VectorView<T>& VectorView<T>::ConjugateSelf() const
  { LapConjugate(*this); return *this; }

  template <class T> const VectorView<T>& VectorView<T>::Swap(
      size_t i1, size_t i2) const
  { 
    TMVAssert(i1 < size() && i2 < size());
    if (i1 != i2) swap(ref(i1),ref(i2)); 
    return *this; 
  }

  template <class T> const VectorView<T>& VectorView<T>::Permute(
      const size_t* p, size_t i1, size_t i2) const
  { 
    TMVAssert(i2 <= size());
    TMVAssert(i1 <= i2);
    for(size_t i=i1;i<i2;i++) Swap(i,p[i]);
    return *this; 
  }

  template <class T> const VectorView<T>& VectorView<T>::ReversePermute(
      const size_t* p, size_t i1, size_t i2) const
  { 
    TMVAssert(i2 <= size());
    TMVAssert(i1 <= i2);
    for(size_t i=i2;i>i1;) { i--; Swap(i,p[i]); }
    return *this; 
  }

  template <class T> const VectorView<T>& VectorView<T>::ReverseSelf() const
  {
    if (isconj()) {
      if (step() == 1) {
	VIt<T,Unit,Conj> it2=begin()+size()-1;
	for(VIt<T,Unit,Conj> it1=begin();it1<it2;++it1,--it2) swap(*it1,*it2);
      } else {
	VIt<T,Step,Conj> it2=begin()+size()-1;
	for(VIt<T,Step,Conj> it1=begin();it1<it2;++it1,--it2) swap(*it1,*it2);
      }
    } else {
      if (step() == 1) {
	VIt<T,Unit,NonConj> it2=begin()+size()-1;
	for(VIt<T,Unit,NonConj> it1=begin();it1<it2;++it1,--it2) 
	  swap(*it1,*it2);
      } else {
	VIt<T,Step,NonConj> it2=begin()+size()-1;
	for(VIt<T,Step,NonConj> it1=begin();it1<it2;++it1,--it2) 
	  swap(*it1,*it2);
      }
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


  template <class T> const VectorView<T>& VectorView<T>::Sort(
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

  template <class T, StepItType S1, ConjItType C1, StepItType S2> 
    inline void NonBlasCopy(CVIt<T,S1,C1> v1it, VIt<T,S2,NonConj> v2it,
	size_t size)
    { for(size_t i=size;i>0;--i,++v1it,++v2it) *v2it = *v1it; }
  template <class T> inline void NonBlasCopy(CVIt<T,Unit,NonConj> v1it, 
      VIt<T,Unit,NonConj> v2it, size_t size)
  { memmove(v2it.GetP(),v1it.GetP(),size*sizeof(T)); }
#ifdef BLAS
  template <class T, StepItType S1, ConjItType C1, StepItType S2> 
    inline void BlasCopy(CVIt<T,S1,C1> v1it, VIt<T,S2,NonConj>& v2it,
	size_t size)
    { NonBlasCopy(v1it,v2it,size); }
  template <StepItType S1, StepItType S2> inline void BlasCopy(
      CVIt<double,S1,NonConj> v1it, VIt<double,S2,NonConj>& v2it, size_t size)
  { cblas_dcopy(size,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step()); }
  template <StepItType S1, StepItType S2> inline void BlasCopy(
      CVIt<complex<double>,S1,NonConj>& v1it,
      VIt<complex<double>,S2,NonConj>& v2it, size_t size) 
  { cblas_zcopy(size,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step()); }

#ifndef NOFLOAT
  template <StepItType S1, StepItType S2> inline void BlasCopy(
      CVIt<float,S1,NonConj> v1it, VIt<float,S2,NonConj>& v2it, size_t size)
  { cblas_scopy(size,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step()); }
  template <StepItType S1, StepItType S2> inline void BlasCopy(
      CVIt<complex<float>,S1,NonConj>& v1it,
      VIt<complex<float>,S2,NonConj>& v2it, size_t size) 
  { cblas_ccopy(size,v1it.GetP(),v1it.step(),v2it.GetP(),v2it.step()); }
#endif
#endif
  template <class T, StepItType S1, ConjItType C1, StepItType S2> 
    inline void DoCopy(CVIt<T,S1,C1> v1it, VIt<T,S2,NonConj> v2it,
	size_t size)
    { 
#ifdef BLAS
      if (v1it.step() > 0 && v2it.step() > 0) BlasCopy(v1it,v2it,size);
      else 
#endif
	NonBlasCopy(v1it,v2it,size); 
    }

  template <class T> inline void CallCopy(
      const GenVector<T>& v1, const VectorView<T>& v2)
  {
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v2.ct()==NonConj);
    if (v1.isconj())
      if (v1.step() == 1)
	if (v2.step() == 1) 
	  DoCopy(CVIt<T,Unit,Conj>(v1.begin()),
	      VIt<T,Unit,NonConj>(v2.begin()),v1.size()); 
	else 
	  DoCopy(CVIt<T,Unit,Conj>(v1.begin()),
	      VIt<T,Step,NonConj>(v2.begin()),v1.size()); 
      else
	if (v2.step() == 1) 
	  DoCopy(CVIt<T,Step,Conj>(v1.begin()),
	      VIt<T,Unit,NonConj>(v2.begin()),v1.size()); 
	else
	  DoCopy(CVIt<T,Step,Conj>(v1.begin()),
	      VIt<T,Step,NonConj>(v2.begin()),v1.size()); 
    else
      if (v1.step() == 1)
	if (v2.step() == 1) 
	  DoCopy(CVIt<T,Unit,NonConj>(v1.begin()),
	      VIt<T,Unit,NonConj>(v2.begin()),v1.size()); 
	else
	  DoCopy(CVIt<T,Unit,NonConj>(v1.begin()),
	      VIt<T,Step,NonConj>(v2.begin()),v1.size()); 
      else
	if (v2.step() == 1)
	  DoCopy(CVIt<T,Step,NonConj>(v1.begin()),
	      VIt<T,Unit,NonConj>(v2.begin()),v1.size()); 
	else
	  DoCopy(CVIt<T,Step,NonConj>(v1.begin()),
	      VIt<T,Step,NonConj>(v2.begin()),v1.size()); 
  }
  template <class T> inline void Copy(
      const GenVector<T>& v1, const VectorView<T>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    if (v1.size() > 0)
      if (v2.isconj()) Copy(v1.Conjugate(),v2.Conjugate());
      else if ((v1.step()<0 && v2.step()<0) ||
	  v2.step()==-1 || (v1.step()==-1 && v2.step()!=1))
	CallCopy(v1.Reverse(),v2.Reverse());
      else CallCopy(v1,v2);
  }

  //
  // Swap Vectors
  //

  template <class T, StepItType S1, ConjItType C1, StepItType S2> 
    inline void NonBlasSwap(
      VIt<T,S1,C1> v1it, VIt<T,S2,NonConj> v2it, size_t size)
  { for(size_t i=size;i>0;--i,++v1it,++v2it) swap(*v1it,*v2it); }
#ifdef BLAS
  template <class T, StepItType S1, ConjItType C1, StepItType S2> 
    inline void BlasSwap(
	VIt<T,S1,C1> v1it, VIt<T,S2,NonConj> v2it, size_t size)
    { NonBlasSwap(v1it,v2it,size); }
  template <StepItType S1, StepItType S2> inline void BlasSwap(
      VIt<double,S1,NonConj> v1, VIt<double,S2,NonConj> v2, size_t size) 
  { cblas_dswap(size,v1.GetP(),v1.step(),v2.GetP(),v2.step()); }
  template <StepItType S1, StepItType S2> inline void BlasSwap(
      VIt<complex<double>,S1,NonConj>& v1, 
      VIt<complex<double>,S2,NonConj>& v2, size_t size) 
  { cblas_zswap(size,v1.GetP(),v1.step(),v2.GetP(),v2.step()); }
#ifndef NOFLOAT
  template <StepItType S1, StepItType S2> inline void BlasSwap(
      VIt<float,S1,NonConj> v1, VIt<float,S2,NonConj> v2, size_t size) 
  { cblas_sswap(size,v1.GetP(),v1.step(),v2.GetP(),v2.step()); }
  template <StepItType S1, StepItType S2> inline void BlasSwap(
      VIt<complex<float>,S1,NonConj> v1, 
      VIt<complex<float>,S2,NonConj> v2, size_t size) 
  { cblas_cswap(size,v1.GetP(),v1.step(),v2.GetP(),v2.step()); }
#endif
#endif
  template <class T, StepItType S1, ConjItType C1, StepItType S2> 
    inline void DoSwap(VIt<T,S1,C1> v1it, VIt<T,S2,NonConj> v2it, size_t size)
    { 
#ifdef BLAS
      if (v1it.step() > 0 && v2it.step() > 0) BlasSwap(v1it,v2it,size);
      else 
#endif
	NonBlasSwap(v1it,v2it,size); 
    }
  template <class T> inline void Swap(
      const VectorView<T>& v1, const VectorView<T>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    if (v1.size() > 0)
      if (v2.isconj()) Swap(v1.Conjugate(),v2.Conjugate());
      else if (v1.step()<0 && v2.step()<0) Swap(v1.Reverse(),v2.Reverse());
      else
	if (v1.isconj())
	  if (v1.step() == 1)
	    if (v2.step() == 1) 
	      DoSwap(VIt<T,Unit,Conj>(v1.begin()),
		  VIt<T,Unit,NonConj>(v2.begin()),v1.size()); 
	    else 
	      DoSwap(VIt<T,Unit,Conj>(v1.begin()),
		  VIt<T,Step,NonConj>(v2.begin()),v1.size()); 
	  else
	    if (v2.step() == 1) 
	      DoSwap(VIt<T,Step,Conj>(v1.begin()),
		  VIt<T,Unit,NonConj>(v2.begin()),v1.size()); 
	    else 
	      DoSwap(VIt<T,Step,Conj>(v1.begin()),
		  VIt<T,Step,NonConj>(v2.begin()),v1.size()); 
	else
	  if (v1.step() == 1)
	    if (v2.step() == 1) 
	      DoSwap(VIt<T,Unit,NonConj>(v1.begin()),
		  VIt<T,Unit,NonConj>(v2.begin()),v1.size()); 
	    else
	      DoSwap(VIt<T,Unit,NonConj>(v1.begin()),
		  VIt<T,Step,NonConj>(v2.begin()),v1.size()); 
	  else
	    if (v2.step() == 1)
	      DoSwap(VIt<T,Step,NonConj>(v1.begin()),
		  VIt<T,Unit,NonConj>(v2.begin()),v1.size()); 
	    else 
	      DoSwap(VIt<T,Step,NonConj>(v1.begin()),
		  VIt<T,Step,NonConj>(v2.begin()),v1.size()); 
  }

  template <class T> inline void Swap(
      const VectorView<T>& v1, Vector<T>& v2) 
  { 
    TMVAssert(v1.size()==v2.size());
    if (v1.size() > 0) 
      if (v1.isconj())
	if (v1.step() == 1)
	  DoSwap(VIt<T,Unit,Conj>(v1.begin()),v2.begin(),v1.size()); 
	else
	  DoSwap(VIt<T,Step,Conj>(v1.begin()),v2.begin(),v1.size()); 
      else
	if (v1.step() == 1)
	  DoSwap(VIt<T,Unit,NonConj>(v1.begin()),v2.begin(),v1.size());
	else
	  DoSwap(VIt<T,Step,NonConj>(v1.begin()),v2.begin(),v1.size()); 
  }

  template <class T> inline void Swap(Vector<T>& v1, Vector<T>& v2) 
  { 
    TMVAssert(v1.size()==v2.size());
    if (v1.size()>0) DoSwap(v1.begin(),v2.begin(),v1.size()); 
  }

  //
  // op ==
  //
  template <class T> bool operator==(
      const GenVector<T>& v1, const GenVector<T>& v2)
  {
    if (v1.size() != v2.size()) return false;
    else if (v1.SameAs(v2)) return true;
    else if (v1.isconj()) return (v1.Conjugate() == v2.Conjugate());
    else if (v1.step() == 1) {
      CVIt<T,Unit,NonConj> v1it=v1.begin();
      if (v2.isconj()) {
	if (v2.step() == 1) {
	  CVIt<T,Unit,Conj> v2it=v2.begin();
	  for(size_t i=v2.size();i>0;--i,++v1it,++v2it) 
	    if ( (*v1it) != (*v2it) ) return false;
	} 
	else {
	  CVIt<T,Step,Conj> v2it=v2.begin();
	  for(size_t i=v2.size();i>0;--i,++v1it,++v2it) 
	    if ( (*v1it) != (*v2it) ) return false;
	}
      } else {
	if (v2.step() == 1) {
	  CVIt<T,Unit,NonConj> v2it=v2.begin();
	  for(size_t i=v2.size();i>0;--i,++v1it,++v2it) 
	    if ( (*v1it) != (*v2it) ) return false;
	} 
	else {
	  CVIt<T,Step,NonConj> v2it=v2.begin();
	  for(size_t i=v2.size();i>0;--i,++v1it,++v2it) 
	    if ( (*v1it) != (*v2it) ) return false;
	}
      }
    } else {
      CVIt<T,Step,NonConj> v1it=v1.begin();
      if (v2.isconj()) {
	if (v2.step() == 1) {
	  CVIt<T,Unit,Conj> v2it=v2.begin();
	  for(size_t i=v2.size();i>0;--i,++v1it,++v2it) 
	    if ( (*v1it) != (*v2it) ) return false;
	} 
	else {
	  CVIt<T,Step,Conj> v2it=v2.begin();
	  for(size_t i=v2.size();i>0;--i,++v1it,++v2it) 
	    if ( (*v1it) != (*v2it) ) return false;
	}
      } else {
	if (v2.step() == 1) {
	  CVIt<T,Unit,NonConj> v2it=v2.begin();
	  for(size_t i=v2.size();i>0;--i,++v1it,++v2it) 
	    if ( (*v1it) != (*v2it) ) return false;
	} 
	else {
	  CVIt<T,Step,NonConj> v2it=v2.begin();
	  for(size_t i=v2.size();i>0;--i,++v1it,++v2it) 
	    if ( (*v1it) != (*v2it) ) return false;
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
    fout << size() << " (";
    if (isconj()) {
      if (step() == 1) {
	CVIt<T,Unit,Conj> it = begin();
	for(size_t i=size();i>0;--i,++it) fout << " " << *it << " ";
      } else {
	CVIt<T,Step,Conj> it = begin();
	for(size_t i=size();i>0;--i,++it) fout << " " << *it << " ";
      }
    } else {
      if (step() == 1) {
	CVIt<T,Unit,NonConj> it = begin();
	for(size_t i=size();i>0;--i,++it) fout << " " << *it << " ";
      } else {
	CVIt<T,Step,NonConj> it = begin();
	for(size_t i=size();i>0;--i,++it) fout << " " << *it << " ";
      }
    }
    fout << " )";
  }

  template <class T> void VectorView<T>::Read(istream& fin) const
  {
    char paren;
    fin >> paren;
    if (!fin || paren != '(') tmv_error("reading ( in Vector::Read");
    if (isconj()) {
      if (step() == 1) {
	VIt<T,Unit,Conj> it = begin();
	for(size_t i=size();i>0;--i,++it) fin >> *it;
      } else {
	VIt<T,Step,Conj> it = begin();
	for(size_t i=size();i>0;--i,++it) fin >> *it;
      }
    } else {
      if (step() == 1) {
	VIt<T,Unit,NonConj> it = begin();
	for(size_t i=size();i>0;--i,++it) fin >> *it;
      } else {
	VIt<T,Step,NonConj> it = begin();
	for(size_t i=size();i>0;--i,++it) fin >> *it;
      }
    }
    if (!fin) tmv_error("reading value in Vector::Read");
    fin >> paren;
    if (paren != ')') tmv_error("reading ) in Vector::Read");
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

  template <class T> istream& operator>>(istream& fin, Vector<T>* v)
  {
    size_t n;
    fin >> n;
    if (!fin) tmv_error("reading size in Vector::Read");
    v = new Vector<T>(n);
    v->View().Read(fin);
    return fin;
  }


#define InstFile "TMV_Vector.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


