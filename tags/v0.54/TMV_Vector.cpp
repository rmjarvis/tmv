
#include "TMV_Blas.h"
#include "TMV_Vector.h"
#include "TMV_VectorArith.h"
#include "TMV_VectorAux.h"

namespace tmv {

  bool FALSE = false;
  int Lap_info = 0;
  char Blas_ch_N = 'N';
  char Blas_ch_C = 'C';
  char Blas_ch_T = 'T';
  char Blas_ch_L = 'L';
  char Blas_ch_R = 'R';
  char Blas_ch_U = 'U';

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
    int n=v.size();
    int s=v.step();
    return BLASNAME(dnrm2) (BLASV(n),BLASP(v.cptr()),BLASV(s));
  }
  template <> inline double BlasNorm2(const GenVector<complex<double> >& v)
  { 
    int n=v.size();
    int s=v.step();
    return BLASNAME(dznrm2) (BLASV(n),BLASP(v.cptr()),BLASV(s));
  }
#ifndef NOFLOAT
  template <> inline float BlasNorm2(const GenVector<float>& v)
  {
    int n=v.size();
    int s=v.step();
    return BLASNAME(snrm2) (BLASV(n),BLASP(v.cptr()),BLASV(s));
  }
  template <> inline float BlasNorm2(const GenVector<complex<float> >& v)
  { 
    int n=v.size();
    int s=v.step();
    return BLASNAME(scnrm2) (BLASV(n),BLASP(v.cptr()),BLASV(s));
  }
#endif
#endif
  template <class T> RealType(T) GenVector<T>::Norm2() const
  { 
    if (step() < 0)
      return Reverse().Norm2();
    else  {
#ifdef BLAS
      return BlasNorm2(*this); 
#else
      return NonBlasNorm2(*this); 
#endif
    }
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

  template <class T> inline RealType(T) DoSumAbsElements(const GenVector<T>& v)
  {
    TMVAssert(v.step() >= 0);
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
    int n=v.size();
    int s=v.step();
    return BLASNAME(dasum) (BLASV(n),BLASP(v.cptr()),BLASV(s));
  }
#ifndef NOFLOAT
  template <> inline float DoSumAbsElements(const GenVector<float>& v)
  { 
    TMVAssert(v.step()>=0);
    int n=v.size();
    int s=v.step();
    return BLASNAME(sasum) (BLASV(n),BLASP(v.cptr()),BLASV(s));
  }
#endif
#ifdef XLAP
  template <> inline double DoSumAbsElements(
      const GenVector<complex<double> >& v)
  { 
    TMVAssert(v.step()>=0);
    int n = v.size();
    int s = v.step();
    return LAPNAMEX(dzsum1) (LAPV(n),LAPP(v.cptr()),LAPV(s)); 
  }
#ifndef NOFLOAT
  template <> inline float DoSumAbsElements(
      const GenVector<complex<float> >& v)
  { 
    TMVAssert(v.step()>=0);
    int n = v.size();
    int s = v.step();
    return LAPNAMEX(scsum1) (LAPV(n),LAPP(v.cptr()),LAPV(s)); 
  }
#endif
#endif // XLAP
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

  template <class T> inline T FindMinElement(const GenVector<T>& v, size_t& imin)
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
  template <class T> inline T FindMaxElement(const GenVector<T>& v, size_t& imax)
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
  template <class T> inline RealType(T) FindMaxAbsElement(
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
  template <class T> inline RealType(T) FindMinAbsElement(
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
  template <> inline double FindMaxAbsElement(
      const GenVector<double>& v, size_t& imax)
  {
    int n=v.size();
    int s=v.step();
    imax = BLASNAME(idamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
    --imax;
#endif
    TMVAssert(imax < v.size());
    return abs(v[imax]);
  }
#ifdef BLASAMIN
  template <> inline double FindMinAbsElement(
      const GenVector<double>& v, size_t& imin)
  {
    int n=v.size();
    int s=v.step();
    imin = BLASNAME(idamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
    --imin;
#endif
    TMVAssert(imin < v.size());
    return abs(v[imin]);
  }
#endif // BLASAMIN
#ifndef NOFLOAT
  template <> inline float FindMaxAbsElement(
      const GenVector<float>& v, size_t& imax)
  {
    int n=v.size();
    int s=v.step();
    imax = BLASNAME(isamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
    --imax;
#endif
    TMVAssert(imax < v.size());
    return abs(v[imax]);
  }
#ifdef BLASAMIN
  template <> inline float FindMinAbsElement(
      const GenVector<float>& v, size_t& imin)
  {
    int n=v.size();
    int s=v.step();
    imin = BLASNAME(isamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
    --imin;
#endif
    TMVAssert(imin < v.size());
    return abs(v[imin]);
  }
#endif // BLASAMIN
#endif // FLOAT
#endif // BLAS

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

  template <class T> inline void NonLapConjugate(const VectorView<complex<T>,CStyle>& v)
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
  template <class T> inline void NonLapConjugate(const VectorView<T,CStyle>& ) {}
#ifdef XLAP
  template <class T> inline void LapConjugate(const VectorView<T,CStyle>& v)
  { NonLapConjugate(v); }
  template <> inline void LapConjugate(const VectorView<complex<double>,CStyle>& v)
  { 
    int n = v.size();
    int s = v.step();
    LAPNAMEX(zlacgv) (LAPV(n),LAPP(v.ptr()),LAPV(s)); 
  }
#ifndef NOFLOAT
  template <> inline void LapConjugate(const VectorView<complex<float>,CStyle>& v)
  {
    int n = v.size();
    int s = v.step();
    LAPNAMEX(clacgv) (LAPV(n),LAPP(v.ptr()),LAPV(s)); 
  }
#endif
#endif // XLAP
  template <class T, IndexStyle I> 
    const VectorView<T,I>& VectorView<T,I>::ConjugateSelf() const
    {
      if (step() < 0) Reverse().ConjugateSelf();
      else {
#ifdef XLAP
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
	  default : TMVAssert(FALSE);
	}
      }

      // Use default copy, op=, destructor

      size_t GetI() const { return itsi; }
      size_t GetVal() const { return itsvalue; }
      bool operator<(const VTIndex& rhs) const
      { return itsvalue < rhs.itsvalue; }
      operator size_t() const { return itsi; }

  };

  template <class T> inline ostream& operator<<(ostream& os, VTIndex<T>& i)
  { os << i.GetVal(); return os; }

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
	    default : TMVAssert(FALSE); return false;
	  }
	} else {
	  switch(comp) {
	    case REAL_COMP : return REAL(x) > REAL(y);
	    case ABS_COMP : return abs(x) > abs(y);
	    case IMAG_COMP : return IMAG(x) > IMAG(y);
	    case ARG_COMP : return ARG(x) > ARG(y);
	    default : TMVAssert(FALSE); return false;
	  }
	} 
      }

    private:

      const ADType ad;
      const COMPType comp;
  };

#ifdef NOSTLSORT
  template <class IT, class COMP>
    inline void Sort3(IT x1, IT x2, IT x3, const COMP& comp)
    {
      if (comp(*x3,*x1)) {
	if (comp(*x1,*x2)) { // x3 < x1 < x2
	  iterator_traits<IT>::value_type temp = *x1;
	  *x1 = *x3;
	  *x3 = *x2;
	  *x2 = temp;
	} else if (comp(*x2,*x3)) { // x2 < x3 < x1
	  iterator_traits<IT>::value_type temp = *x1;
	  *x1 = *x2;
	  *x2 = *x3;
	  *x3 = temp;
	} else { // x3 <= x2 <= x1 and x3 < x1
	  swap(*x1,*x3);
	}
      } else {
	if (comp(*x2,*x1)) { // x2 < x1 <= x3
	  swap(*x1,*x2);
	} else if (comp(*x3,*x2)) { // x1 <= x3 < x2
	  swap(*x2,*x3);
	} else { // x1 <= x2 <= x3
	  // nothing to do
	}
      }
    }

  inline void* Address(void* x) { return x; }
  template <class T> inline void* Address(VIter<T> x) { return x.GetP(); }

  template <class IT, class COMP>
    inline void QuickSort(IT begin, IT end, const COMP& comp)
  {
    TMVAssert(end-begin >= 0);
    if (end-begin <= 3) {
      if (end-begin == 3) { // 3 elements
	Sort3(begin,begin+1,begin+2,comp);
      } else if (end-begin == 2) { // 2 elements
        if (comp(*(begin+1),*begin)) swap(*begin,*(begin+1));
      } // else 0 or 1 element
      return;
    } else {
      IT mid = begin + (end-begin)/2;
      Sort3(begin,mid,end-1,comp);
      iterator_traits<IT>::value_type pivot = *mid;
      swap(*mid,*(end-2));
      IT left = begin+1;
      IT right = end-3;
      while (left < right) {
	while (!comp(*right,pivot) && left < right) --right;
	while (!comp(pivot,*left) && left < right) ++left;
	if (left < right) swap(*left,*right);
      }
      TMVAssert(left == right);
      if (comp(*left,pivot)) ++left;
      swap(*left,*(end-2));
      QuickSort(begin,left,comp);
      QuickSort(left+1,end,comp);
    }
  }

  template <class IT>
    inline void QuickSort(const IT& begin, const IT& end)
    { QuickSort(begin,end,std::less<iterator_traits<IT>::value_type>()); }
#endif

  template <class T, IndexStyle I> const VectorView<T,I>& VectorView<T,I>::Sort(
      size_t* P, ADType ad, COMPType comp) const
  {
    if (P) {
      std::vector<VTIndex<T> > newindex(size());
      for(size_t i=0;i<size();++i) {
	newindex[i] = VTIndex<T>(ref(i),i,ad,comp);
      }
#ifdef NOSTLSORT
      QuickSort(newindex.begin(),newindex.end());
#else
      std::sort(newindex.begin(),newindex.end());
#endif

      ConvertIndexToPermute(size(),newindex,P);
      Permute(P);

    } else {
      const Compare<T> cc(ad,comp);
#ifdef NOSTLSORT
      QuickSort(begin(),end(),cc);
#else
      std::sort(begin(),end(),cc);
#endif
    }
    return *this;
  }

  //
  // Copy Vectors
  //

  template <bool c1, class T> inline void NonBlasCopy(
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
  { 
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(dcopy) (BLASV(n),BLASP(v1.cptr()),BLASV(s1)
	,BLASP(v2.ptr()),BLASV(s2));
  }
  template <> inline void BlasCopy(
      const GenVector<complex<double> >& v1,
      const VectorView<complex<double> >& v2)
  { 
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(zcopy) (BLASV(n),BLASP(v1.cptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2));
  }
#ifndef NOFLOAT
  template <> inline void BlasCopy(
      const GenVector<float>& v1, const VectorView<float>& v2)
  {
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(scopy) (BLASV(n),BLASP(v1.cptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2));
  }
  template <> inline void BlasCopy(
      const GenVector<complex<float> >& v1,
      const VectorView<complex<float> >& v2)
  { 
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(ccopy) (BLASV(n),BLASP(v1.cptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2));
  }
#endif // FLOAT
#endif // BLAS

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

  template <class T> inline void NonBlasSwap(
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
  { 
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(dswap) (BLASV(n),BLASP(v1.ptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2));
  }
  template <> inline void BlasSwap(
      const VectorView<complex<double> >& v1, 
      const VectorView<complex<double> >& v2)
  { 
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(zswap) (BLASV(n),BLASP(v1.ptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2));
  }
#ifndef NOFLOAT
  template <> inline void BlasSwap(
      const VectorView<float>& v1, const VectorView<float>& v2)
  { 
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(sswap) (BLASV(n),BLASP(v1.ptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2));
  }
  template <> inline void BlasSwap(
      const VectorView<complex<float> >& v1, 
      const VectorView<complex<float> >& v2)
  { 
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    BLASNAME(cswap) (BLASV(n),BLASP(v1.ptr()),BLASV(s1),
	BLASP(v2.ptr()),BLASV(s2));
  }
#endif
#endif // BLAS

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

  template <class T> void VectorReadError<T>::Write(ostream& os) const throw()
  {
    os<<"TMV Read Error: Reading istream input for Vector\n";
    if (exp != got) {
      os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
    }
    if (v.get() && s != v->size()) {
      os<<"Wrong size: expected "<<v->size()<<", got "<<s<<".\n";
    }
    if (!is) {
      if (iseof) {
	os<<"Input stream reached end-of-file prematurely.\n";
      } else if (isbad) {
	os<<"Input stream is corrupted.\n";
      } else {
	os<<"Input stream cannot read next character.\n";
      }
    }
    if (v.get()) {
      os<<"The portion of the Vector which was successfully read is: \n";
      os<<"( ";
      for(size_t ii=0;ii<i;++ii)
	os<<' '<<(*v)(ii)<<' ';
      os<<" )\n";
    }
  }

  template <class T, IndexStyle I> void VectorView<T,I>::Read(istream& is) const
  {
    T* p = ptr();
    const int s = step();

    char paren;
    is >> paren;
    if (!is || paren != '(') 
      throw VectorReadError<T>(0,*this,is,'(',is?paren:'(');
    if (s == 1) 
      for(size_t i=size();i>0;--i,++p) {
	is >> *p;
	if (!is) 
	  throw VectorReadError<T>(size()-i,*this,is);
      }
    else
      for(size_t i=size();i>0;--i,p+=s) {
	is >> *p;
	if (!is) 
	  throw VectorReadError<T>(size()-i,*this,is);
      }
    is >> paren;
    if (!is || paren != ')') 
      throw VectorReadError<T>(size(),*this,is,')',is?paren:')');
    if (isconj()) ConjugateSelf();
  }

  template <class T> istream& operator>>(istream& is, const VectorView<T>& v)
  {
    size_t n;
    is >> n;
    if (!is) 
      throw VectorReadError<T>(is);
    if (n != v.size()) 
      throw VectorReadError<T>(v,is,n);
    v.Read(is);
    return is;
  }

  template <class T, IndexStyle I> istream& operator>>(istream& is, 
      auto_ptr<Vector<T,I> >& v)
  {
    size_t n;
    is >> n;
    if (!is) 
      throw VectorReadError<T>(is);
    v.reset(new Vector<T,I>(n));
    v->View().Read(is);
    return is;
  }


#define InstFile "TMV_Vector.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


