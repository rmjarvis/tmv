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
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_MinMax.h"
#include "tmv/TMV_NormV.h"
#include "tmv/TMV_SortV.h"
#include "tmv/TMV_SwapV.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_VectorIO.h"

namespace tmv {

  const int XX = UNKNOWN;

  // First warn_out from TMV_Base.h

  std::ostream* warn_out = &std::cout;

  // And also some things from TMV_Blas.h
  int Lap_info = 0;

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
  // Copy Vectors
  //

  template <class T, bool C> 
  void InstCopy(const ConstVectorView<T,XX,C>& v1, VectorView<T> v2)
  {
    TMVAssert(v1.cptr() != v2.cptr());
    if (v2.step() == 1) 
    {
      VectorView<T,1> v2u = v2.UnitView();
      if (v1.step() == 1) InlineCopy(v1.UnitView(),v2u);
      else InlineCopy(v1,v2u);
    } 
    else 
      if (v1.step() == 1) InlineCopy(v1.UnitView(),v2);
      else InlineCopy(v1,v2);
  }

#ifdef BLAS
#define TMV_INST_SKIP_BLAS
#ifdef TMV_INST_DOUBLE
  template <> 
  void InstCopy(const ConstVectorView<double>& v1, VectorView<double> v2)
  {
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v1.cptr() != v2.cptr());
    int n=v2.size();
    if (n == 0) return;
    int s1=v1.step();
    int s2=v2.step();
    const double* v1p = v1.cptr();
    if (s1 < 0) v1p += (n-1)*s1;
    double* v2p = v2.ptr();
    if (s2 < 0) v2p += (n-1)*s2;
    BLASNAME(dcopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
  }
  template <> 
  void InstCopy(const ConstVectorView<std::complex<double> >& v1,
      VectorView<std::complex<double> > v2)
  {
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v1.cptr() != v2.cptr());
    int n=v2.size();
    if (n == 0) return;
    int s1=v1.step();
    int s2=v2.step();
    const std::complex<double>* v1p = v1.cptr();
    if (s1 < 0) v1p += (n-1)*s1;
    std::complex<double>* v2p = v2.ptr();
    if (s2 < 0) v2p += (n-1)*s2;
    BLASNAME(zcopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
  }
  template <> 
  void InstCopy(const ConstVectorView<std::complex<double>,XX,true>& v1, 
      VectorView<std::complex<double> > v2)
  {
    InstCopy(v1.Conjugate(),v2);
    v2.ConjugateSelf();
  }
#endif
#ifdef TMV_INST_FLOAT
  template <> 
  void InstCopy(const ConstVectorView<float>& v1, VectorView<float> v2)
  {
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v1.cptr() != v2.cptr());
    int n=v2.size();
    if (n == 0) return;
    int s1=v1.step();
    int s2=v2.step();
    const float* v1p = v1.cptr();
    if (s1 < 0) v1p += (n-1)*s1;
    float* v2p = v2.ptr();
    if (s2 < 0) v2p += (n-1)*s2;
    BLASNAME(scopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
  }
  template <> 
  void InstCopy(const ConstVectorView<std::complex<float> >& v1,
      VectorView<std::complex<float> > v2)
  {
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v1.cptr() != v2.cptr());
    int n=v2.size();
    if (n == 0) return;
    int s1=v1.step();
    int s2=v2.step();
    const std::complex<float>* v1p = v1.cptr();
    if (s1 < 0) v1p += (n-1)*s1;
    std::complex<float>* v2p = v2.ptr();
    if (s2 < 0) v2p += (n-1)*s2;
    BLASNAME(ccopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
  }
  template <> 
  void InstCopy(const ConstVectorView<std::complex<float>,XX,true>& v1, 
      VectorView<std::complex<float> > v2)
  {
    InstCopy(v1.Conjugate(),v2);
    v2.ConjugateSelf();
  }
#endif // FLOAT
#endif // BLAS

  //
  // Swap Vectors
  //

  template <class T, bool C> 
  void InstSwap(VectorView<T,XX,C> v1, VectorView<T> v2)
  {
    if (v1.step() == 1) 
    {
      VectorView<T,1,C> v1u = v1.UnitView();
      if (v2.step() == 1) 
      {
        VectorView<T,1> v2u = v2.UnitView();
        InlineSwap(v1u,v2u);
      }
      else InlineSwap(v1u,v2);
    }
    else 
    {
      if (v2.step() == 1) 
      {
        VectorView<T,1> v2u = v2.UnitView();
        InlineSwap(v1,v2u);
      }
      else InlineSwap(v1,v2);
    }
  }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
  template <> 
  void InstSwap(VectorView<double> v1, VectorView<double> v2)
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
  template <> 
  void InstSwap(VectorView<std::complex<double> > v1, 
      VectorView<std::complex<double> > v2)
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
  template <> 
  void InstSwap(VectorView<std::complex<double>,XX,true> v1, 
      VectorView<std::complex<double> > v2)
  {
    InstSwap(v1.Conjugate(),v2);
    v1.ConjugateSelf();
    v2.ConjugateSelf();
  }
#endif
#ifdef TMV_INST_FLOAT
  template <> 
  void InstSwap(VectorView<float> v1, VectorView<float> v2)
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
  template <> 
  void InstSwap(VectorView<std::complex<float> > v1, 
      VectorView<std::complex<float> > v2)
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
  template <> 
  void InstSwap(VectorView<std::complex<float>,XX,true> v1, 
      VectorView<std::complex<float> > v2)
  {
    InstSwap(v1.Conjugate(),v2);
    v1.ConjugateSelf();
    v2.ConjugateSelf();
  }
#endif
#endif // BLAS


  //
  // ReverseSelf
  //

  template <class T> 
  void InstReverseSelf(VectorView<T> v)
  {
    if (v.step() == 1) 
    {
      VectorView<T,1> vu = v.UnitView();
      InlineReverseSelf(vu);
    }
    else InlineReverseSelf(v);
  }


  //
  // ConjugateSelf
  //

  template <class T> 
  void InstConjugateSelf(VectorView<T> v)
  {
    if (v.step() == 1) 
    {
      VectorView<T,1> vu = v.UnitView();
      InlineConjugateSelf(vu);
    }
    else InlineConjugateSelf(v);
  }

#ifdef ELAP
#ifdef TMV_INST_DOUBLE
  template <> 
  void InstConjugateSelf(VectorView<std::complex<double> > v)
  { 
    int n = v.size();
    int s = v.step();
    std::complex<double>* vp = v.ptr();
    if (s < 0) vp += (n-1)*s;
    LAPNAME(zlacgv) (LAPV(n),LAPP(vp),LAPV(s)); 
  }
#endif
#ifdef TMV_INST_FLOAT
  template <> 
  void InstConjugateSelf(VectorView<std::complex<float> > v)
  {
    int n = v.size();
    int s = v.step();
    std::complex<float>* vp = v.ptr();
    if (s < 0) vp += (n-1)*s;
    LAPNAME(clacgv) (LAPV(n),LAPP(vp),LAPV(s)); 
  }
#endif
#endif // ELAP


  //
  // SumElements
  //

  template <class T> 
  T InstSumElements(const ConstVectorView<T>& v)
  {
    return (v.step() == 1) ? 
      InlineSumElements(v.UnitView()) : 
      InlineSumElements(v);
  }


  //
  // SumAbsElements
  //

  template <class T> 
  RealType(T) InstSumAbsElements(const ConstVectorView<T>& v)
  {
    if (v.step() == 1) return InlineSumAbsElements(v.UnitView());
    else return InlineSumAbsElements(v);
  }

  template <class T> 
  RealType(T) InstSumAbs2Elements(const ConstVectorView<T>& v)
  {
    if (v.step() == 1) return InlineSumAbs2Elements(v.UnitView());
    else return InlineSumAbs2Elements(v);
  }

#ifdef BLAS
#ifndef BLASNORETURN
#define BLAS_SumAbsElements // So we know not to do this in the .inst file
#ifdef TMV_INST_DOUBLE
  template <> 
  double InstSumAbs2Elements(const ConstVectorView<double>& v)
  {
    int n = v.size();
    if (n == 0) return double(0);
    int s = v.step();
    // If s == 0, the BLAS standard is to return 0 for the sum,
    // rather than the correct value.  Weird.
    // The non-blas version works correctly, but just do it here anyway.
    if (s == 0) return n*TMV_ABS(v[0]);
    const double* vp = v.cptr();
    if (s < 0) vp += (n-1)*s;
    return BLASNAME(dasum) (BLASV(n),BLASP(vp),BLASV(s));
  }
  template <> 
  double InstSumAbsElements(const ConstVectorView<double>& v)
  { return InstSumAbs2Elements(v); }
  template <> 
  double InstSumAbs2Elements(const ConstVectorView<std::complex<double> >& v)
  {
    int n = v.size();
    if (n == 0) return double(0);
    int s = v.step();
    if (s == 0) return n*TMV_ABS(v[0]);
    const std::complex<double>* vp = v.cptr();
    if (s < 0) vp += (n-1)*s;
    return BLASNAME(dzasum) (BLASV(n),BLASP(vp),BLASV(s));
  }
#endif
#ifdef TMV_INST_FLOAT
  template <> 
  float InstSumAbs2Elements(const ConstVectorView<float>& v)
  { 
    int n = v.size();
    if (n == 0) return float(0);
    int s = v.step();
    if (s == 0) return n*TMV_ABS(v[0]);
    const float* vp = v.cptr();
    if (s < 0) vp += (n-1)*s;
    return BLASNAME(sasum) (BLASV(n),BLASP(vp),BLASV(s));
  }
  template <> 
  float InstSumAbsElements(const ConstVectorView<float>& v)
  { return InstSumAbs2Elements(v); }
  template <> 
  float InstSumAbs2Elements(const ConstVectorView<std::complex<float> >& v)
  {
    int n = v.size();
    if (n == 0) return float(0);
    int s = v.step();
    if (s == 0) return n*TMV_ABS(v[0]);
    const std::complex<float>* vp = v.cptr();
    if (s < 0) vp += (n-1)*s;
    return BLASNAME(scasum) (BLASV(n),BLASP(vp),BLASV(s));
  }
#endif
#ifdef XLAP
#define BLAS_SumAbsElements_complex
#ifdef TMV_INST_DOUBLE
  template <> 
  double InstSumAbsElements(
      const ConstVectorView<std::complex<double> >& v)
  { 
    int n = v.size();
    if (n == 0) return double(0);
    int s = v.step();
    if (s == 0) return n*TMV_ABS(v[0]);
    const std::complex<double>* vp = v.cptr();
    if (s < 0) vp += (n-1)*s;
    return LAPNAME(dzsum1) (LAPV(n),LAPP(vp),LAPV(s)); 
  }
#endif
#ifdef TMV_INST_FLOAT
  template <> 
  float InstSumAbsElements(const ConstVectorView<std::complex<float> >& v)
  { 
    int n = v.size();
    if (n == 0) return float(0);
    int s = v.step();
    if (s == 0) return n*TMV_ABS(v[0]);
    const std::complex<float>* vp = v.cptr();
    if (s < 0) vp += (n-1)*s;
    return LAPNAME(scsum1) (LAPV(n),LAPP(vp),LAPV(s)); 
  }
#endif
#endif // XLAP
#endif // BLASNORETURN
#endif // BLAS


  //
  // MaxElement
  //

  template <class T> 
  T InstMaxElement(const ConstVectorView<T>& v, int*const imax)
  {
    if (v.step() == 1) return InlineMaxElement(v.UnitView(),imax);
    else return InlineMaxElement(v,imax);
  }

  //
  // MaxAbsElement
  //

  template <class T> 
  RealType(T) InstMaxAbsElement(const ConstVectorView<T>& v, int*const imax)
  {
    if (v.step() == 1) return InlineMaxAbsElement(v.UnitView(),imax);
    else return InlineMaxAbsElement(v,imax);
  }

  template <class T> 
  RealType(T) InstMaxAbs2Element(const ConstVectorView<T>& v, int*const imax)
  {
    if (v.step() == 1) return InlineMaxAbs2Element(v.UnitView(),imax);
    else return InlineMaxAbs2Element(v,imax);
  }

#ifdef BLAS
#define BLAS_MaxAbs
  // These return values seem to work, so I don't guard this segment 
  // with BLASNORETURN
#ifdef TMV_INST_DOUBLE
  template <> 
  double InstMaxAbs2Element(const ConstVectorView<double>& v, int*const imax)
  {
    int n=v.size();
    int s=v.step();
    if (s == 0) {
      if (imax) *imax = 0;
      return TMV_ABS(v[0]);
    } else if (s < 0) {
      double max = InstMaxAbs2Element(v.Reverse(),imax);
      if (imax) *imax = n-1-(*imax);
      return max;
    } else {
      int i = BLASNAME(idamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
      --i;
#endif
      if (imax) *imax = i;
      return TMV_ABS(v[i]);
    }
  }
  template <> 
  double InstMaxAbsElement(const ConstVectorView<double>& v, int*const imax)
  { return InstMaxAbs2Element(v,imax); }
  template <> 
  double InstMaxAbs2Element(
      const ConstVectorView<std::complex<double> >& v, int*const imax)
  {
    int n=v.size();
    int s=v.step();
    if (s == 0) {
      if (imax) *imax = 0;
      return TMV_ABS(v[0]);
    } else if (s < 0) {
      double max = InstMaxAbs2Element(v.Reverse(),imax);
      if (imax) *imax = n-1-(*imax);
      return max;
    } else {
      int i = BLASNAME(izamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
      --i;
#endif
      if (imax) *imax = i;
      return TMV_ABS(v[i]);
    }
  }
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
  template <> 
  float InstMaxAbs2Element(const ConstVectorView<float>& v, int*const imax)
  {
    int n=v.size();
    int s=v.step();
    if (s == 0) {
      if (imax) *imax = 0;
      return TMV_ABS(v[0]);
    } else if (s < 0) {
      float max = InstMaxAbs2Element(v.Reverse(),imax);
      if (imax) *imax = n-1-(*imax);
      return max;
    } else {
      int i = BLASNAME(isamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
      --i;
#endif
      if (imax) *imax = i;
      return TMV_ABS(v[i]);
    }
  }
  template <> 
  float InstMaxAbsElement(const ConstVectorView<float>& v, int*const imax)
  { return InstMaxAbs2Element(v,imax); }
  template <> 
  float InstMaxAbs2Element(
      const ConstVectorView<std::complex<float> >& v, int*const imax)
  {
    int n=v.size();
    int s=v.step();
    if (s == 0) {
      if (imax) *imax = 0;
      return TMV_ABS(v[0]);
    } else if (s < 0) {
      float max = InstMaxAbs2Element(v.Reverse(),imax);
      if (imax) *imax = n-1-(*imax);
      return max;
    } else {
      int i = BLASNAME(icamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
      --i;
#endif
      if (imax) *imax = i;
      return TMV_ABS(v[i]);
    }
  }
#endif // FLOAT
#endif // BLAS


  //
  // MinElement
  //

  template <class T> 
  T InstMinElement(const ConstVectorView<T>& v, int*const imin)
  {
    if (v.step() == 1) return InlineMinElement(v.UnitView(),imin);
    else return InlineMinElement(v,imin);
  }


  //
  // MinAbsElement
  //

  template <class T> 
  RealType(T) InstMinAbsElement(const ConstVectorView<T>& v, int*const imin)
  {
    if (v.step() == 1) return InlineMinAbsElement(v.UnitView(),imin);
    else return InlineMinAbsElement(v,imin);
  }

  template <class T> 
  RealType(T) InstMinAbs2Element(const ConstVectorView<T>& v, int*const imin)
  {
    if (v.step() == 1) return InlineMinAbs2Element(v.UnitView(),imin);
    else return InlineMinAbs2Element(v,imin);
  }

#ifdef BLAS
#ifdef BLASIDAMIN
#ifdef TMV_INST_DOUBLE
  template <> 
  double InstMinAbs2Element(const ConstVectorView<double>& v, int*const imin)
  {
    int n=v.size();
    int s=v.step();
    if (s == 0) {
      if (imin) *imin = 0;
      return TMV_ABS(v[0]);
    } else if (s < 0) {
      double min = v.Reverse().MinAbsElement(imin);
      if (imin) *imin = n-1-(*imin);
      return min;
    } else {
      int i = BLASNAME(idamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
      --i;
#endif
      if (imin) *imin = i;
      return TMV_ABS(v[i]);
    }
  }
  template <> 
  double InstMinAbsElement(const ConstVectorView<double>& v, int*const imin)
  { InstMinAbs2Element(v,imin); }
  template <> 
  double InstMinAbs2Element(
      const ConstVectorView<std::complex<double> >& v, int*const imin)
  {
    int n=v.size();
    int s=v.step();
    if (s == 0) {
      if (imin) *imin = 0;
      return TMV_ABS(v[0]);
    } else if (s < 0) {
      double min = v.Reverse().MinAbsElement(imin);
      if (imin) *imin = n-1-(*imin);
      return min;
    } else {
      int i = BLASNAME(izamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
      --i;
#endif
      if (imin) *imin = i;
      return TMV_ABS(v[i]);
    }
  }
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
  template <> 
  float InstMinAbs2Element(const ConstVectorView<float>& v, int*const imin)
  {
    int n=v.size();
    int s=v.step();
    if (s == 0) {
      if (imin) *imin = 0;
      return TMV_ABS(v[0]);
    } else if (s < 0) {
      float min = v.Reverse().MinAbsElement(imin);
      if (imin) *imin = n-1-(*imin);
      return min;
    } else {
      int i = BLASNAME(isamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
      --i;
#endif
      if (imin) *imin = i;
      return TMV_ABS(v[i]);
    }
  }
  template <> 
  float InstMinAbsElement(const ConstVectorView<float>& v, int*const imin)
  { InstMinAbs2Element(v,imin); }
  template <> 
  float InstMinAbs2Element(
      const ConstVectorView<std::complex<float> >& v, int*const imin)
  {
    int n=v.size();
    int s=v.step();
    if (s == 0) {
      if (imin) *imin = 0;
      return TMV_ABS(v[0]);
    } else if (s < 0) {
      float min = v.Reverse().MinAbsElement(imin);
      if (imin) *imin = n-1-(*imin);
      return min;
    } else {
      int i = BLASNAME(icamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
      --i;
#endif
      if (imin) *imin = i;
      return TMV_ABS(v[i]);
    }
  }
#endif // FLOAT
#endif // BLASIDAMIN
#endif // BLAS


  //
  // NormSq
  //

  template <class T> 
  RealType(T) InstNormSq(const ConstVectorView<T>& v)
  {
    if (v.step() == 1) return InlineNormSq(v.UnitView());
    else return InlineNormSq(v);
  }

  template <class T> 
  RealType(T) InstNormSq(const ConstVectorView<T>& v, const RealType(T) scale)
  {
    if (v.step() == 1) return InlineNormSq(v.UnitView(),scale);
    else return InlineNormSq(v,scale);
  }


  //
  // Norm2
  //

  template <class T> 
  RealType(T) InstNorm2(const ConstVectorView<T>& v)
  {
    if (v.step() == 1) return InlineNorm2(v.UnitView());
    else return InlineNorm2(v);
  }

#ifdef BLAS
#ifndef BLASNORETURN
#define BLAS_Norm2
#ifdef TMV_INST_DOUBLE
  template <> 
  double InstNorm2(const ConstVectorView<double>& v)
  {
    int n=v.size();
    int s=v.step();
    const double* vp = v.cptr();
    if (s < 0) vp += (n-1)*s;
    return BLASNAME(dnrm2) (BLASV(n),BLASP(vp),BLASV(s));
  }
  template <> 
  double InstNorm2(const ConstVectorView<std::complex<double> >& v)
  { 
    int n=v.size();
    int s=v.step();
    const std::complex<double>* vp = v.cptr();
    if (s < 0) vp += (n-1)*s;
    return BLASNAME(dznrm2) (BLASV(n),BLASP(vp),BLASV(s));
  }
#endif
#ifdef TMV_INST_FLOAT
  template <> 
  float InstNorm2(const ConstVectorView<float>& v)
  {
    int n=v.size();
    int s=v.step();
    const float* vp = v.cptr();
    if (s < 0) vp += (n-1)*s;
    return BLASNAME(snrm2) (BLASV(n),BLASP(vp),BLASV(s));
  }
  template <> 
  float InstNorm2(const ConstVectorView<std::complex<float> >& v)
  { 
    int n=v.size();
    int s=v.step();
    const std::complex<float>* vp = v.cptr();
    if (s < 0) vp += (n-1)*s;
    return BLASNAME(scnrm2) (BLASV(n),BLASP(vp),BLASV(s));
  }
#endif
#endif // BLASNORETURN
#endif // BLAS


  //
  // Sort
  //

  template <class T> 
  void InstSort(VectorView<T> v, ADType ad, COMPType comp)
  { 
    if (v.step() == 1) 
    {
      VectorView<T,1> vu = v.UnitView();
      InlineSort(vu,ad,comp);
    }
    else InlineSort(v,ad,comp); 
  }
  template <class T> 
  void InstSort(VectorView<T> v, int*const P, ADType ad, COMPType comp)
  { 
    if (v.step() == 1) 
    {
      VectorView<T,1> vu = v.UnitView();
      InlineSort(vu,P,ad,comp);
    }
    else InlineSort(v,P,ad,comp); 
  }


  //
  // I/O
  //

  template <class T, bool C> 
  void InstWrite(std::ostream& os, const ConstVectorView<T,XX,C>& v)
  { 
    if (v.step() == 1) InlineWrite(os,v.UnitView()); 
    else InlineWrite(os,v); 
  }

  template <class T, bool C> 
  void InstWrite(std::ostream& os, const ConstVectorView<T,XX,C>& v,
      RealType(T) thresh)
  {
    if (v.step() == 1) InlineWrite(os,v.UnitView(),thresh); 
    else InlineWrite(os,v,thresh); 
  }

  template <class T, bool C> 
  void InstRead(std::istream& is, VectorView<T,XX,C> v)
  {
    if (v.step() == 1)
    {
      VectorView<T,1,C> vu = v.UnitView();
      InlineRead(is,vu); 
    }
    else InlineRead(is,v); 
  }

#define InstFile "TMV_Vector.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


