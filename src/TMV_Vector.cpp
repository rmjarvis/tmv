///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TMV_Blas.h"
#include <iostream>
#include <algorithm>
#include <limits>
#include <iostream>
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VIt.h"
#include "TMV_ConvertIndex.h"
#include "portable_platform.h"

#ifndef NOTHROW
#include "TMV_VectorRE.h"
#endif

namespace tmv {

#define RT RealType(T)

  // First the things from TMV_Base.h
  template <class T> RT Epsilon()
  { return std::numeric_limits<RT>::epsilon(); }

  template <class T> RT SqrtEpsilon()
  { return SQRT(Epsilon<T>()); }

  bool FALSE = false;
  std::ostream* warn_out = &std::cout;

  // Next one thing from TMV_ListInit.h

  ListInitClass ListInit;

  // And also some things from TMV_Blas.h
  int Lap_info = 0;
  char Blas_ch_N = 'N';
  char Blas_ch_C = 'C';
  char Blas_ch_T = 'T';
  char Blas_ch_L = 'L';
  char Blas_ch_R = 'R';
  char Blas_ch_U = 'U';

  void LAP_Results(const char* fn)
  {
    if (Lap_info < 0) {
#ifdef NOTHROW
      std::cerr<<"info < 0 returned by LAPACK function "<<fn<<std::endl;
      exit(1);
#else
      throw Error("info < 0 returned by LAPACK function ",fn);
#endif
    }
  }

  void LAP_Results(const int lwork_opt, const int m, const int n,
      const int lwork, const char* fn)
  {
    LAP_Results(fn);
    if (lwork_opt > lwork) { 
      std::ostringstream s;
      s << "LAPACK function " << fn << 
      " requested more workspace than provided";
      TMV_Warning(s.str());
      s.str(std::string());
      s<<"for matrix with m,n = "<<m<<','<<n<<std::endl;
      TMV_Warning(s.str());
      s.str(std::string());
      s<<"Given: "<<lwork<<", requested "<<lwork_opt<<std::endl;
      TMV_Warning(s.str());
    }
  }

  //
  // Access
  //

  template <class T> T GenVector<T>::cref(int i) const
  {
    const T* vi = cptr() + int(i)*step();
    return (ct() == Conj) ? CONJ(*vi) : *vi;
  }

  template <class T, IndexStyle I> RefType(T) VectorView<T,I>::ref(int i) const
  {
    T* vi = ptr() + int(i)*step();
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
      std::cout<<"istep ("<<istep<<") cannot be 0\n";
    }
    if (i1 < 0 || i1 >= int(size())) {
      ok = false;
      std::cout<<"first element ("<<i1<<") must be in 0 -- "<<size()-1<<std::endl;
    }
    if (i2-istep < 0 || i2-istep >= int(size())) {
      ok = false;
      std::cout<<"last element ("<<i2-istep<<") must be in 0 -- "<<size()-1<<std::endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      std::cout<<"range ("<<i2-i1<<") must be multiple of istep ("<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      std::cout<<"n elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
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
      std::cout<<"istep ("<<istep<<") cannot be 0\n";
    }
    if (i1 < 1 || i1 > int(this->size())) {
      ok = false;
      std::cout<<"first element ("<<i1<<") must be in 1 -- "<<this->size()<<std::endl;
    }
    if (i2 < 1 || i2 > int(this->size())) {
      ok = false;
      std::cout<<"last element ("<<i2<<") must be in 1 -- "<<this->size()<<std::endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      std::cout<<"range ("<<i2-i1<<") must be multiple of istep ("<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      std::cout<<"n elements ("<<(i2-i1)/istep+1<<") must be positive\n";
    }
    return ok;
  }

  //
  // Norm2
  //
  template <class T> RT GenVector<T>::NormSq(const RT scale) const
  {
    if (size() == 0) return RT(0);

    const int s = step();
    if (s == 1) {
      if (IsComplex(T())) return Flatten().NormSq(scale);
      else {
        RT sum(0);
        const T* p = cptr();
        if (scale == RT(1)) 
          for(int i=size();i>0;--i,++p) sum += NORM(*p);
        else
          for(int i=size();i>0;--i,++p) sum += NORM(scale* (*p));
        return sum;
      }
    }
    else if (s < 0) return Reverse().NormSq(scale);
    else if (s == 0) return RT(size()) * NORM(scale*(*cptr()));
    else {
      RT sum(0);
      const T* p = cptr();
      if (scale == RT(1))
        for(int i=size();i>0;--i,p+=s) sum += NORM(*p);
      else
        for(int i=size();i>0;--i,p+=s) sum += NORM(scale* (*p));
      return sum;
    }
  }
  template <class T> static RT DoNorm2(const GenVector<T>& v)
  { 
    //std::cout<<"DoNorm2: v = "<<v<<std::endl;
    const RT eps = Epsilon<T>();
    const RT inveps = RT(1)/eps;

    RT vmax = v.MaxAbsElement();
    //std::cout<<"vmax = "<<vmax<<std::endl;
    if (vmax == RT(0)) return RT(0);
    else if (vmax * vmax * eps == RT(0)) {
      //std::cout<<"vmax^2 * eps = "<<vmax*vmax*eps<<std::endl;
      // Then we need to rescale, since underflow has caused rounding errors
      // Epsilon is a pure power of 2, so no rounding errors from rescaling.
      RT scale = inveps;
      vmax *= scale;
      while (vmax < eps) { scale *= inveps; vmax *= inveps; }
      return SQRT(v.NormSq(scale))/scale;
    } else if (RT(1) / (vmax) == RT(0)) {
      //std::cout<<"1/vmax = "<<RT(1)/vmax<<std::endl;
      // Then vmax is already inf, so no hope of making more accurate.
      // (And specifically, next section would cause infinite loop.)
      return vmax;
    } else if (RT(1) / (vmax*vmax) == RT(0)) {
      //std::cout<<"1/vmax^2 = "<<RT(1)/(vmax*vmax)<<std::endl;
      // Then we have overflow, so we need to rescale:
      RT scale = eps;
      vmax *= scale;
      while (vmax > RT(1)) { scale *= eps; vmax *= eps; }
      return SQRT(v.NormSq(scale))/scale;
    } 
    return SQRT(v.NormSq());

  }
#ifdef BLAS
#ifndef BLASNORETURN
#ifdef INST_DOUBLE
  template <> double DoNorm2(const GenVector<double>& v)
  {
    int n=v.size();
    int s=v.step();
    const double* vp = v.cptr();
    return BLASNAME(dnrm2) (BLASV(n),BLASP(vp),BLASV(s));
  }
  template <> double DoNorm2(const GenVector<std::complex<double> >& v)
  { 
    int n=v.size();
    int s=v.step();
    const std::complex<double>* vp = v.cptr();
    return BLASNAME(dznrm2) (BLASV(n),BLASP(vp),BLASV(s));
  }
#endif
#ifdef INST_FLOAT
  template <> float DoNorm2(const GenVector<float>& v)
  {
    int n=v.size();
    int s=v.step();
    const float* vp = v.cptr();
    return BLASNAME(snrm2) (BLASV(n),BLASP(vp),BLASV(s));
  }
  template <> float DoNorm2(const GenVector<std::complex<float> >& v)
  { 
    int n=v.size();
    int s=v.step();
    const std::complex<float>* vp = v.cptr();
    return BLASNAME(scnrm2) (BLASV(n),BLASP(vp),BLASV(s));
  }
#endif
#endif // BLASNORETURN
#endif // BLAS
  template <class T> RT GenVector<T>::Norm2() const
  { 
    if (size() == 0) return RT(0);
    else if (step() < 0) return DoNorm2(Reverse());
    else  return DoNorm2(*this);
  }

  //
  // SumElements
  //
  template <class T> T GenVector<T>::SumElements() const
  {
    if (size() == 0) return T(0);

    const int s = step();
    if (s == 1) {
      T sum(0);
      const T* p = cptr();
      for(int i=size();i>0;--i,++p) sum += *p;
      return this->isconj() ? CONJ(sum) : sum;
    }
    else if (s < 0) return Reverse().SumElements();
    else if (s == 0) return RT(size()) * (*cptr());
    else {
      T sum(0);
      const T* p = cptr();
      for(int i=size();i>0;--i,p+=s) sum += *p;
      return this->isconj() ? CONJ(sum) : sum;
    }
  }

  template <class T> static RT DoSumAbsElements(const GenVector<T>& v)
  {
    TMVAssert(v.step() > 0);
    const int s = v.step();
    RT sum(0);
    const T* p = v.cptr();
    if (s == 1) 
      for(int i=v.size();i>0;--i,++p) sum += ABS(*p);
    else 
      for(int i=v.size();i>0;--i,p+=s) sum += ABS(*p);
    return sum;
  }
#ifdef BLAS
#ifndef BLASNORETURN
#ifdef INST_DOUBLE
  template <> double DoSumAbsElements(const GenVector<double>& v)
  { 
    TMVAssert(v.step()>=0);
    int n=v.size();
    int s=v.step();
    return BLASNAME(dasum) (BLASV(n),BLASP(v.cptr()),BLASV(s));
  }
#endif
#ifdef INST_FLOAT
  template <> float DoSumAbsElements(const GenVector<float>& v)
  { 
    TMVAssert(v.step()>=0);
    int n=v.size();
    int s=v.step();
    return BLASNAME(sasum) (BLASV(n),BLASP(v.cptr()),BLASV(s));
  }
#endif
#ifdef XLAP
#ifdef INST_DOUBLE
  template <> double DoSumAbsElements(
      const GenVector<std::complex<double> >& v)
  { 
    TMVAssert(v.step()>=0);
    int n = v.size();
    int s = v.step();
    return LAPNAME(dzsum1) (LAPV(n),LAPP(v.cptr()),LAPV(s)); 
  }
#endif
#ifdef INST_FLOAT
  template <> float DoSumAbsElements(
      const GenVector<std::complex<float> >& v)
  { 
    TMVAssert(v.step()>=0);
    int n = v.size();
    int s = v.step();
    return LAPNAME(scsum1) (LAPV(n),LAPP(v.cptr()),LAPV(s)); 
  }
#endif
#endif // XLAP
#endif // BLASNORETURN
#endif // BLAS

  template <class T> RT GenVector<T>::SumAbsElements() const
  {
    if (size() == 0) return RT(0);
    else if (step() > 0) return DoSumAbsElements(*this); 
    else if (step() < 0) return DoSumAbsElements(Reverse()); 
    // If s == 0, the BLAS standard is to return 0 for the sum,
    // rather than the correct value.  Weird.
    // The non-blas version works correctly, but just do it here anyway.
    else return size() * ABS(*cptr());
  }

  //
  // Find Min/Max Element
  //

  template <class T> static T FindMinElement(
      const GenVector<T>& v, int& imin)
  {
    TMVAssert(v.size() > 0);
    TMVAssert(v.step() > 0);

    const T* p = v.cptr();
    const int s = v.step();
    T min = *p;
    imin = 0;
    int i=1;
    if (s == 1) {
      ++p;
      for(int k=v.size()-1;k>0; --k,++p,++i) {
        if (REAL(*p) < REAL(min)) {
          min = *p;
          imin = i;
        }
      }
    }
    else {
      p += s;
      for(int k=v.size()-1;k>0; --k,p+=s,++i) {
        if (REAL(*p) < REAL(min)) {
          min = *p;
          imin = i;
        }
      }
    }
    return v.isconj() ? CONJ(min) : min;
  }
  template <class T> static T FindMaxElement(const GenVector<T>& v, int& imax)
  {
    TMVAssert(v.size() > 0);
    TMVAssert(v.step() > 0);

    const T* p = v.cptr();
    const int s = v.step();
    T max = *p;
    imax = 0;
    int i=1;
    if (s == 1) {
      ++p;
      for(int k=v.size()-1;k>0; --k,++p,++i) {
        if (REAL(*p) > REAL(max)) {
          max = *p;
          imax = i;
        }
      }
    }
    else {
      p += s;
      for(int k=v.size()-1;k>0; --k,p+=s,++i) {
        if (REAL(*p) > REAL(max)) {
          max = *p;
          imax = i;
        }
      }
    }
    return v.isconj() ? CONJ(max) : max;
  }
  template <class T> static RT FindMaxAbsElement(
      const GenVector<T>& v, int& imax)
  {
    TMVAssert(v.size() > 0);
    TMVAssert(v.step() > 0);

    const T* p = v.cptr();
    const int s = v.step();
    RT max = ABS(*p);
    imax = 0;
    int i=1;
    if (s == 1) {
      ++p;
      for(int k=v.size()-1;k>0; --k,++p,++i) {
        RT absval = ABS(*p);
        if (absval > max) { 
          max = absval; 
          imax = i;
        }
      }
    }
    else {
      p += s;
      for(int k=v.size()-1;k>0; --k,p+=s,++i) {
        RT absval = ABS(*p);
        if (absval > max) { 
          max = absval; 
          imax = i;
        }
      }
    }
    return max;
  }
  template <class T> static RT FindMinAbsElement(
      const GenVector<T>& v, int& imin)
  {
    TMVAssert(v.size() > 0);
    TMVAssert(v.step() > 0);

    const T* p = v.cptr();
    const int s = v.step();
    RT min = ABS(*p);
    imin = 0;
    int i=1;
    if (s == 1) {
      ++p;
      for(int k=v.size()-1;k>0; --k,++p,++i) {
        RT absval = ABS(*p);
        if (absval < min) {
          min = absval; 
          imin = i;
        }
      }
    } else {
      p += s;
      for(int k=v.size()-1;k>0; --k,p+=s,++i) {
        RT absval = ABS(*p);
        if (absval < min) {
          min = absval; 
          imin = i;
        }
      }
    }
    return IsReal(T()) ? min : SQRT(min);
  }
#ifdef BLAS
  // These return values seem to work, so I don't guard this segment 
  // with BLASNORETURN
#ifdef INST_DOUBLE
  template <> double FindMaxAbsElement(
      const GenVector<double>& v, int& imax)
  {
    int n=v.size();
    int s=v.step();
    imax = BLASNAME(idamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
    --imax;
#endif
    TMVAssert(imax < int(v.size()));
    return ABS(v[imax]);
  }
#ifdef BLASIDAMIN
  template <> double FindMinAbsElement(
      const GenVector<double>& v, int& imin)
  {
    int n=v.size();
    int s=v.step();
    imin = BLASNAME(idamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
    --imin;
#endif
    TMVAssert(imin < int(v.size()));
    return ABS(v[imin]);
  }
#endif // BLASIDAMIN
#endif
#ifdef INST_FLOAT
  template <> float FindMaxAbsElement(
      const GenVector<float>& v, int& imax)
  {
    int n=v.size();
    int s=v.step();
    imax = BLASNAME(isamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
    --imax;
#endif
    TMVAssert(imax < int(v.size()));
    return ABS(v[imax]);
  }
#ifdef BLASIDAMIN
  template <> float FindMinAbsElement(
      const GenVector<float>& v, int& imin)
  {
    int n=v.size();
    int s=v.step();
    imin = BLASNAME(isamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
    --imin;
#endif
    TMVAssert(imin < int(v.size()));
    return ABS(v[imin]);
  }
#endif // BLASAMIN
#endif // FLOAT
#endif // BLAS

  template <class T> T GenVector<T>::MinElement(int* iminout) const
  {
    if (size() == 0) {
      if (iminout) *iminout = -1;
      return T(0);
    }
    if (step() > 0) {
      int imin;
      T min = FindMinElement(*this,imin);
      TMVAssert(imin < int(size()));
      if (iminout) *iminout = imin;
      return min;
    } else if (step() == 0) {
      if (iminout) *iminout = 0;
      return *cptr();
    } else {
      T min = Reverse().MinElement(iminout);
      if (iminout) *iminout = size()-1-(*iminout);
      return min;
    }
  }
  template <class T> T GenVector<T>::MaxElement(int* imaxout) const
  {
    if (size() == 0) {
      if (imaxout) *imaxout = -1;
      return T(0);
    }
    if (step() > 0) {
      int imax;
      T max = FindMaxElement(*this,imax);
      TMVAssert(imax < int(size()));
      if (imaxout) *imaxout = imax;
      return max;
    } else if (step() == 0) {
      if (imaxout) *imaxout = 0;
      return *cptr();
    } else {
      T max = Reverse().MaxElement(imaxout);
      if (imaxout) *imaxout = size()-1-(*imaxout);
      return max;
    }
  }

  template <class T> RT GenVector<T>::MinAbsElement(
      int* iminout) const
  {
    if (size() == 0) {
      if (iminout) *iminout = -1;
      return RT(0);
    }
    if (step() > 0) {
      int imin;
      RT min = FindMinAbsElement(*this,imin);
      TMVAssert(imin < int(size()));
      if (iminout) *iminout = imin;
      return min;
    } else if (step() == 0) {
      if (iminout) *iminout = 0;
      return ABS(*cptr());
    } else {
      RT min = Reverse().MinAbsElement(iminout);
      if (iminout) *iminout = size()-1-(*iminout);
      return min;
    }
  }
  template <class T> RT GenVector<T>::MaxAbsElement( 
      int* imaxout) const
  {
    if (size() == 0) {
      if (imaxout) *imaxout = -1;
      return RT(0);
    }
    if (step() > 0) {
      int imax;
      RT max = FindMaxAbsElement(*this,imax);
      TMVAssert(imax < int(size()));
      if (imaxout) *imaxout = imax;
      return max;
    } else if (step() == 0) {
      if (imaxout) *imaxout = 0;
      return ABS(*cptr());
    } else {
      RT max = Reverse().MaxAbsElement(imaxout);
      if (imaxout) *imaxout = size()-1-(*imaxout);
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

  template <class T, IndexStyle I> 
  const VectorView<T,I>& VectorView<T,I>::Zero() const
  {
    if (step() == 1) std::fill_n(ptr(),size(),T(0));
    else SetAllTo(0);
    return *this;
  }

  template <class T, IndexStyle I> Vector<T,I>& Vector<T,I>::Zero() 
  {
    std::fill_n(ptr(),size(),T(0));
    return *this;
  }

  template <class T, IndexStyle I> 
  const VectorView<T,I>& VectorView<T,I>::Clip(RT thresh) const
  {
    const int s = step();
    if (s < 0) Reverse().Clip(thresh);
    else if (s == 0) {
      if (ABS(*ptr()) < thresh) *ptr() = T(0);
    } else {
      T* p = ptr();

      if (s == 1) {
        for(int i=size();i>0;--i,++p) 
          if (ABS(*p) < thresh) *p = T(0);
      } else {
        for(int i=size();i>0;--i,p+=s) 
          if (ABS(*p) < thresh) *p = T(0);
      }
    }
    return *this; 
  }

  template <class T, IndexStyle I> Vector<T,I>& Vector<T,I>::Clip(
      RT thresh)
  { View().Clip(thresh); return *this; }

  template <class T, IndexStyle I> 
  const VectorView<T,I>& VectorView<T,I>::SetAllTo(T x) const
  {
    const int s = step();
    if (s < 0) Reverse().SetAllTo(x);
    else if (s == 0) {
#ifdef TMVFLDEBUG
      TMVAssert(ptr() >= first);
      TMVAssert(ptr() < last);
#endif
      *ptr() = x;
    }
    else if (s == 1) {
      std::fill_n(ptr(),size(),x);
    } else {
      T* p = ptr();
      if (this->isconj()) 
        for(int i=size();i>0;--i,p+=s) {
#ifdef TMVFLDEBUG
          TMVAssert(p >= first);
          TMVAssert(p < last);
#endif
          *p = CONJ(x); 
        }
      else
        for(int i=size();i>0;--i,p+=s) {
#ifdef TMVFLDEBUG
          TMVAssert(p >= first);
          TMVAssert(p < last);
#endif
          *p = x; 
        }
    }
    return *this; 
  }

  template <class T, IndexStyle I> Vector<T,I>& Vector<T,I>::SetAllTo(T x) 
  {
    std::fill_n(ptr(),size(),x);
    return *this;
  }

  template <class T, IndexStyle I> 
  const VectorView<T,I>& VectorView<T,I>::AddToAll(T x) const
  {
    const int s = step();
    if (s < 0) Reverse().AddToAll(x);
    else if (s == 0) *ptr() += x;
    else {
      T* p = ptr();

      if (this->isconj()) 
        if (s == 1) 
          for(int i=size();i>0;--i,++p) *p += CONJ(x); 
        else 
          for(int i=size();i>0;--i,p+=s) *p += CONJ(x); 
      else 
        if (s == 1) 
          for(int i=size();i>0;--i,++p) *p += x; 
        else 
          for(int i=size();i>0;--i,p+=s) *p += x; 
    }
    return *this; 
  }

  template <class T, IndexStyle I> Vector<T,I>& Vector<T,I>::AddToAll(T x)
  {
    T* p = ptr();
    for(int i=size();i>0;--i,++p) *p += x;
    return *this;
  }

  template <class T> static void NonLapConjugate(
      const VectorView<std::complex<T>,CStyle>& v)
  {
    TMVAssert(v.step() > 0);

    T* p = v.Real().ptr();
    p++;

    if (v.step() == 1)
      for(int i=v.size();i>0;--i,p+=2) *p = -(*p);
    else if (v.step() == 0) *p = -*p;
    else {
      const int s = 2*v.step();
      for(int i=v.size();i>0;--i,p+=s) *p = -(*p);
    }
  }
  template <class T> static inline void NonLapConjugate(
      const VectorView<T,CStyle>& ) 
  {}
#ifdef ELAP
  template <class T> static inline void LapConjugate(const VectorView<T,CStyle>& v)
  { NonLapConjugate(v); }
#ifdef INST_DOUBLE
  template <> void LapConjugate(
      const VectorView<std::complex<double>,CStyle>& v)
  { 
    int n = v.size();
    int s = v.step();
    LAPNAME(zlacgv) (LAPV(n),LAPP(v.ptr()),LAPV(s)); 
  }
#endif
#ifdef INST_FLOAT
  template <> void LapConjugate(
      const VectorView<std::complex<float>,CStyle>& v)
  {
    int n = v.size();
    int s = v.step();
    LAPNAME(clacgv) (LAPV(n),LAPP(v.ptr()),LAPV(s)); 
  }
#endif
#endif // ELAP
  template <class T, IndexStyle I> 
  const VectorView<T,I>& VectorView<T,I>::ConjugateSelf() const
  {
    if (step() < 0) Reverse().ConjugateSelf();
    else {
#ifdef ELAP
      LapConjugate(*this);
#else
      NonLapConjugate(*this);
#endif
    }
    return *this; 
  }

  template <class T, IndexStyle I> Vector<T,I>& Vector<T,I>::ConjugateSelf()
  { 
    View().ConjugateSelf(); 
    return *this; 
  }


  template <class T, IndexStyle I> 
  const VectorView<T,I>& VectorView<T,I>::DoMakeBasis(int i) const
  { 
    TMVAssert(I==CStyle);
    TMVAssert(i < int(size()));
    Zero();
    ref(i) = T(1);
    return *this; 
  }

  template <class T, IndexStyle I> Vector<T,I>& Vector<T,I>::DoMakeBasis(
      int i)
  {
    if (I == CStyle) { TMVAssert(i<int(size())); }
    else { TMVAssert(i>0 && i<=int(size())); }

    const int ix = (I==CStyle ? i : i-1);
    Zero();
    ref(ix) = T(1);
    return *this;
  }

  template <class T, IndexStyle I> 
  const VectorView<T,I>& VectorView<T,I>::DoSwap(int i1, int i2) const
  {
    TMVAssert(i1 < int(size()));
    TMVAssert(i2 < int(size()));
    TMVAssert(I==CStyle);
    if (i1 != i2) {
      const int s = step();
      if (s == 1) 
        __TMV_SWAP(*(ptr()+i1),*(ptr()+i2));
      else
        __TMV_SWAP(*(ptr()+i1*s),*(ptr()+i2*s));
    }
    return *this;
  }

  template <class T, IndexStyle I> Vector<T,I>& Vector<T,I>::DoSwap(
      int i1, int i2)
  {
    TMVAssert(i1 < int(size()));
    TMVAssert(i2 < int(size()));
    if (i1 != i2)  {
      if (I==CStyle) __TMV_SWAP(*(ptr()+i1),*(ptr()+i2));
      else __TMV_SWAP(*(ptr()+i1-1),*(ptr()+i2-1));
    }
    return *this;
  }

  template <class T, IndexStyle I> 
  const VectorView<T,I>& VectorView<T,I>::DoPermute(
      const int* p, int i1, int i2) const
  { 
    TMVAssert(i2 <= int(size()));
    TMVAssert(i1 <= i2);
    TMVAssert(I==CStyle);
    for(int i=i1;i<i2;++i) DoSwap(i,p[i]);
    return *this; 
  }

  template <class T, IndexStyle I> 
  const VectorView<T,I>& VectorView<T,I>::DoReversePermute(
      const int* p, int i1, int i2) const
  { 
    TMVAssert(i2 <= int(size()));
    TMVAssert(i1 <= i2);
    TMVAssert(I==CStyle);
    for(int i=i2;i>i1;) { --i; DoSwap(i,p[i]); }
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
        for(T* p2=p1+size()-1;p2>p1;++p1,--p2) __TMV_SWAP(*p1,*p2);
      else
        for(T* p2=p1+s*(size()-1);p2>p1;p1+=s,p2-=s) __TMV_SWAP(*p1,*p2);
    }
    return *this;
  }

  template <class T> class VTIndex
  {

  private :

    RT itsvalue;
    int itsi;

  public :

    VTIndex() : itsvalue(RT(0)), itsi(0) {}

    VTIndex(T val, int i, ADType ad, COMPType comp) : 
      itsvalue(0), itsi(i)
    {
      bool neg = ad==DESCEND;
      switch(comp) {
        case REAL_COMP : itsvalue = neg ? -REAL(val) : REAL(val); break;
        case ABS_COMP : itsvalue = neg ? -ABS(val) : ABS(val); break;
        case IMAG_COMP : itsvalue = neg ? -IMAG(val) : IMAG(val); break;
        case ARG_COMP : itsvalue = neg ? -ARG(val) : ARG(val); break;
        default : TMVAssert2(FALSE);
      }
    }

    // Use default copy, op=, destructor

    int GetI() const { return itsi; }
    RT GetVal() const { return itsvalue; }
    bool operator<(const VTIndex& rhs) const
    { return itsvalue < rhs.itsvalue; }
    operator int() const { return itsi; }

  };

  template <class T> static inline std::ostream& operator<<(
      std::ostream& os, VTIndex<T>& i)
  { os << i.GetVal(); return os; }

  template <class T> class Compare
  {
  public:

    Compare(ADType _ad, COMPType _comp) : ad(_ad), comp(_comp) {}

    bool operator()(const T& x, const T& y) const {
      if (ad == ASCEND) {
        switch(comp) {
          case REAL_COMP : return REAL(x) < REAL(y);
          case ABS_COMP : return ABS(x) < ABS(y);
          case IMAG_COMP : return IMAG(x) < IMAG(y);
          case ARG_COMP : return ARG(x) < ARG(y);
          default : TMVAssert2(FALSE); return false;
        }
      } else {
        switch(comp) {
          case REAL_COMP : return REAL(x) > REAL(y);
          case ABS_COMP : return ABS(x) > ABS(y);
          case IMAG_COMP : return IMAG(x) > IMAG(y);
          case ARG_COMP : return ARG(x) > ARG(y);
          default : TMVAssert2(FALSE); return false;
        }
      } 
    }

  private:

    const ADType ad;
    const COMPType comp;
  };

#ifdef NOSTL
  template <class IT, class COMP>
  static void Sort3(IT x1, IT x2, IT x3, const COMP& comp)
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
        __TMV_SWAP(*x1,*x3);
      }
    } else {
      if (comp(*x2,*x1)) { // x2 < x1 <= x3
        __TMV_SWAP(*x1,*x2);
      } else if (comp(*x3,*x2)) { // x1 <= x3 < x2
        __TMV_SWAP(*x2,*x3);
      } else { // x1 <= x2 <= x3
        // nothing to do
      }
    }
  }

  static inline void* Address(void* x) { return x; }
  template <class T> static inline void* Address(VIter<T> x) { return x.GetP(); }

  template <class IT, class COMP>
  static void QuickSort(IT begin, IT end, const COMP& comp)
  {
    TMVAssert(end-begin >= 0);
    if (end-begin <= 3) {
      if (end-begin == 3) { // 3 elements
        Sort3(begin,begin+1,begin+2,comp);
      } else if (end-begin == 2) { // 2 elements
        if (comp(*(begin+1),*begin)) __TMV_SWAP(*begin,*(begin+1));
      } // else 0 or 1 element
      return;
    } else {
      IT mid = begin + (end-begin)/2;
      Sort3(begin,mid,end-1,comp);
      iterator_traits<IT>::value_type pivot = *mid;
      __TMV_SWAP(*mid,*(end-2));
      IT left = begin+1;
      IT right = end-3;
      while (left < right) {
        while (!comp(*right,pivot) && left < right) --right;
        while (!comp(pivot,*left) && left < right) ++left;
        if (left < right) __TMV_SWAP(*left,*right);
      }
      TMVAssert(left == right);
      if (comp(*left,pivot)) ++left;
      __TMV_SWAP(*left,*(end-2));
      QuickSort(begin,left,comp);
      QuickSort(left+1,end,comp);
    }
  }

  template <class IT>
  static inline void QuickSort(const IT& begin, const IT& end)
  { QuickSort(begin,end,std::less<iterator_traits<IT>::value_type>()); }
#endif

  template <class T, IndexStyle I> const VectorView<T,I>& VectorView<T,I>::Sort(
      int* P, ADType ad, COMPType comp) const
  {
    if (P) {
      std::vector<VTIndex<T> > newindex(size());
      const int N = size();
      for(int i=0;i<N;++i) newindex[i] = VTIndex<T>(ref(i),i,ad,comp);
#ifdef NOSTL
      QuickSort(newindex.begin(),newindex.end());
#else
      std::sort(newindex.begin(),newindex.end());
#endif
      ConvertIndexToPermute(size(),newindex,P);
      Permute(P);
    } else {
      // Swap ad necessary accroding the the conj status of the vector:
      if (this->isconj() && (comp==IMAG_COMP || comp==ARG_COMP)) {
        if (ad == ASCEND) ad = DESCEND;
        else ad = ASCEND;
      }
      const Compare<T> cc(ad,comp);
#ifdef NOSTL
      QuickSort(ptr(),ptr()+size(),cc);
#else
      std::sort(ptr(),ptr()+size(),cc);
#endif
    }
    return *this;
  }

  //
  // Special Constructors
  //

  template <class T, IndexStyle I> Vector<T,I> DoBasisVector(
      size_t n, int i)
  {
    if (I == CStyle) { TMVAssert(i<int(n)); }
    else { TMVAssert(i>0 && i<=int(n)); }
    Vector<T,I> temp(n,T(0));
    temp(i) = T(1);
    return temp;
  }

  //
  // Copy Vectors
  //

  template <class T>  void DoCopySameType(
      const GenVector<T>& v1, const VectorView<T>& v2)
  {
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v2.size()>0);
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(!v2.SameAs(v1));
    const T* v1ptr = v1.cptr();
    T* v2ptr = v2.ptr();
    const int step1 = v1.step();
    const int step2 = v2.step();

    if (step1 == 1 && step2 == 1) {
      std::copy(v1ptr,v1ptr+v2.size(),v2ptr);
    }
    else {
      for(int i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2) {
#ifdef TMVFLDEBUG
        TMVAssert(v2ptr >= v2.first);
        TMVAssert(v2ptr < v2.last);
#endif
        *v2ptr = *v1ptr;
      }
    }
  }
#ifdef BLAS
#ifdef INST_DOUBLE
  template <> void DoCopySameType(
      const GenVector<double>& v1, const VectorView<double>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v2.size()>0);
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(!v2.SameAs(v1));
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    const double* v1p = v1.cptr();
    if (s1 < 0) v1p += (n-1)*s1;
    double* v2p = v2.ptr();
    if (s2 < 0) v2p += (n-1)*s2;
    BLASNAME(dcopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
  }
  template <> void DoCopySameType(
      const GenVector<std::complex<double> >& v1,
      const VectorView<std::complex<double> >& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v2.size()>0);
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(!v2.SameAs(v1));
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    const std::complex<double>* v1p = v1.cptr();
    if (s1 < 0) v1p += (n-1)*s1;
    std::complex<double>* v2p = v2.ptr();
    if (s2 < 0) v2p += (n-1)*s2;
    BLASNAME(zcopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
  }
#endif
#ifdef INST_FLOAT
  template <> void DoCopySameType(
      const GenVector<float>& v1, const VectorView<float>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v2.size()>0);
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(!v2.SameAs(v1));
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    const float* v1p = v1.cptr();
    if (s1 < 0) v1p += (n-1)*s1;
    float* v2p = v2.ptr();
    if (s2 < 0) v2p += (n-1)*s2;
    BLASNAME(scopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
  }
  template <> void DoCopySameType(
      const GenVector<std::complex<float> >& v1,
      const VectorView<std::complex<float> >& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v2.size()>0);
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(!v2.SameAs(v1));
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    const std::complex<float>* v1p = v1.cptr();
    if (s1 < 0) v1p += (n-1)*s1;
    std::complex<float>* v2p = v2.ptr();
    if (s2 < 0) v2p += (n-1)*s2;
    BLASNAME(ccopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
  }
#endif // FLOAT
#endif // BLAS

  //
  // Swap Vectors
  //

  template <class T> static void conjswap(T& x, T& y)
  {
    T temp = x;
    x = CONJ(y);
    y = CONJ(temp);
  }

  template <class T> static void DoSwap(
      const VectorView<T>& v1, const VectorView<T>& v2)
  {
    T* v1ptr = v1.ptr();
    T* v2ptr = v2.ptr();
    const int step1 = v1.step();
    const int step2 = v2.step();

    if (v1.isconj())
      if (step1 == 1 && step2 == 1)
        for(int i=v2.size();i>0;--i,++v1ptr,++v2ptr) 
          conjswap(*v1ptr,*v2ptr); 
      else
        for(int i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2) 
          conjswap(*v1ptr,*v2ptr); 
    else
      if (step1 == 1 && step2 == 1)
        for(int i=v2.size();i>0;--i,++v1ptr,++v2ptr) 
          __TMV_SWAP(*v1ptr,*v2ptr); 
      else
        for(int i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2) 
          __TMV_SWAP(*v1ptr,*v2ptr); 
  }
#ifdef BLAS
#ifdef INST_DOUBLE
  template <> void DoSwap(
      const VectorView<double>& v1, const VectorView<double>& v2)
  { 
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    double* v1p = v1.ptr();
    if (s1 < 0) v1p += (n-1)*s1;
    double* v2p = v2.ptr();
    if (s2 < 0) v2p += (n-1)*s2;
    BLASNAME(dswap) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
  }
  template <> void DoSwap(
      const VectorView<std::complex<double> >& v1, 
      const VectorView<std::complex<double> >& v2)
  { 
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    std::complex<double>* v1p = v1.ptr();
    if (s1 < 0) v1p += (n-1)*s1;
    std::complex<double>* v2p = v2.ptr();
    if (s2 < 0) v2p += (n-1)*s2;
    BLASNAME(zswap) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
  }
#endif
#ifdef INST_FLOAT
  template <> void DoSwap(
      const VectorView<float>& v1, const VectorView<float>& v2)
  { 
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    float* v1p = v1.ptr();
    if (s1 < 0) v1p += (n-1)*s1;
    float* v2p = v2.ptr();
    if (s2 < 0) v2p += (n-1)*s2;
    BLASNAME(sswap) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
  }
  template <> void DoSwap(
      const VectorView<std::complex<float> >& v1, 
      const VectorView<std::complex<float> >& v2)
  { 
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    std::complex<float>* v1p = v1.ptr();
    if (s1 < 0) v1p += (n-1)*s1;
    std::complex<float>* v2p = v2.ptr();
    if (s2 < 0) v2p += (n-1)*s2;
    BLASNAME(cswap) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
  }
#endif
#endif // BLAS

  template <class T> void Swap(
      const VectorView<T>& v1, const VectorView<T>& v2)
  { 
    TMVAssert2(v1.size() == v2.size());
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
    if (v1.size() != v2.size()) return false;
    else if (v2.SameAs(v1)) return true;
    else {
      const T1* v1ptr = v1.cptr();
      const T2* v2ptr = v2.cptr();
      const int step1 = v1.step();
      const int step2 = v2.step();

      if (v1.isconj() == v2.isconj()) {
        if (step1 == 1 && step2 == 1) {
          for(int i=v2.size();i>0;--i,++v1ptr,++v2ptr)
            if ( *v1ptr != *v2ptr ) return false;
        } else {
          for(int i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2)
            if ( *v1ptr != *v2ptr ) return false;
        }
      } else {
        if (step1 == 1 && step2 == 1) {
          for(int i=v2.size();i>0;--i,++v1ptr,++v2ptr)
            if ( *v1ptr != CONJ(*v2ptr) ) return false;
        } else {
          for(int i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2)
            if ( *v1ptr != CONJ(*v2ptr) ) return false;
        }
      }
      return true;
    }
  }

  //
  // I/O
  //

  // This bit is to workaround a bug in pgCC that was fixed in version 7.
  // I don't know if versions earlier than 6.1 had the bug, but 
  // I apply the workaround to all version before 7.
  template <class T> inline T Value(T x) { return x; }
#ifdef PLATFORM_COMPILER_PGI
#if PLATFORM_COMPILER_VERSION < 0x070000
  inline double Value(long double x) { return double(x); }
  inline std::complex<double> Value(std::complex<long double> x) 
  { return std::complex<double>(x); }
#endif
#endif

  template <class T> void GenVector<T>::Write(std::ostream& os) const
  {
    const T* p = cptr();
    const int s = step();

    os << size() << " (";
    if (this->isconj()) 
      if (s == 1) 
        for(int i=size();i>0;--i,++p) os << " " << Value(CONJ(*p)) << " ";
      else 
        for(int i=size();i>0;--i,p+=s) os << " " << Value(CONJ(*p)) << " ";
    else 
      if (step() == 1) 
        for(int i=size();i>0;--i,++p) os << " " << Value(*p) << " ";
      else
        for(int i=size();i>0;--i,p+=s) os << " " << Value(*p) << " ";
    os << ")";
  }

  template <class T> void GenVector<T>::Write(std::ostream& os,
      RT thresh) const
  {
    const T* p = cptr();
    const int s = step();

    os << size() << " (";
    if (this->isconj()) 
      if (s == 1) 
        for(int i=size();i>0;--i,++p) 
          os << " " << Value((ABS(*p) < thresh ? T(0) : CONJ(*p))) << " ";
      else 
        for(int i=size();i>0;--i,p+=s) 
          os << " " << Value((ABS(*p) < thresh ? T(0) : CONJ(*p))) << " ";
    else 
      if (step() == 1) 
        for(int i=size();i>0;--i,++p) 
          os << " " << Value((ABS(*p) < thresh ? T(0) : *p)) << " ";
      else
        for(int i=size();i>0;--i,p+=s) 
          os << " " << Value((ABS(*p) < thresh ? T(0) : *p)) << " ";
    os << ")";
  }

#ifndef NOTHROW
  template <class T> void VectorReadError<T>::Write(std::ostream& os) const throw()
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
      os<<"(";
      for(int ii=0;ii<i;++ii)
        os<<' '<<(*v)(ii)<<' ';
      os<<")\n";
    }
  }
#endif

  template <class T, IndexStyle I> void VectorView<T,I>::Read(
      std::istream& is) const
  {
    T* p = ptr();
    const int s = step();

    char paren;
    is >> paren;
    if (!is || paren != '(') 
#ifdef NOTHROW
    { std::cerr<<"Vector ReadError: "<<paren<<" != (\n"; exit(1); }
#else
    throw VectorReadError<T>(0,*this,is,'(',is?paren:'(');
#endif
    if (s == 1) 
      for(int i=size();i>0;--i,++p) {
        is >> *p;
        if (!is) 
#ifdef NOTHROW
        { std::cerr<<"Vector ReadError: !is \n"; exit(1); }
#else
        throw VectorReadError<T>(size()-i,*this,is);
#endif
      }
    else
      for(int i=size();i>0;--i,p+=s) {
        is >> *p;
        if (!is) 
#ifdef NOTHROW
        { std::cerr<<"Vector ReadError: !is \n"; exit(1); }
#else
        throw VectorReadError<T>(size()-i,*this,is);
#endif
      }
    is >> paren;
    if (!is || paren != ')') 
#ifdef NOTHROW
    { std::cerr<<"Vector ReadError: "<<paren<<" != )\n"; exit(1); }
#else
    throw VectorReadError<T>(size(),*this,is,')',is?paren:')');
#endif
    if (this->isconj()) ConjugateSelf();
  }

  template <class T> std::istream& operator>>(
      std::istream& is, const VectorView<T>& v)
  {
    size_t n;
    is >> n;
    if (!is) 
#ifdef NOTHROW
    { std::cerr<<"Vector ReadError: !is \n"; exit(1); }
#else
    throw VectorReadError<T>(is);
#endif
    if (n != v.size()) 
#ifdef NOTHROW
    { std::cerr<<"Vector ReadError: Wrong size \n"; exit(1); }
#else
    throw VectorReadError<T>(v,is,n);
#endif
    v.Read(is);
    return is;
  }

  template <class T, IndexStyle I> std::istream& operator>>(
      std::istream& is, auto_ptr<Vector<T,I> >& v)
  {
    size_t n;
    is >> n;
    if (!is) 
#ifdef NOTHROW
    { std::cerr<<"Vector ReadError: !is \n"; exit(1); }
#else
    throw VectorReadError<T>(is);
#endif
    v.reset(new Vector<T,I>(n));
    v->View().Read(is);
    return is;
  }

#ifdef BLAS
#define INST_SKIP_BLAS
#endif

#define InstFile "TMV_Vector.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


