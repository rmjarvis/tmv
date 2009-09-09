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
#include "tmv/TMV_MultXV.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SwapV.h"

namespace tmv {

  const int XX = UNKNOWN;

  //
  // MultXV: v *= x
  // 

  template <class T, class V> 
  static void DoMultXV(const T x, V& v)
  { 
    if (x == T(-1))
      InlineMultXV(Scaling<-1,T>(x),v); 
    else if (x == T(0))
      v.Zero();
    else if (x != T(1))
      InlineMultXV(Scaling<0,T>(x),v); 
  }

  template <class T, class V> 
  static void DoMultXV(const std::complex<T> x, V& v)
  {
    if (imag(x) == T(0)) 
    {
      if (real(x) == T(-1))
        InlineMultXV(Scaling<-1,T>(real(x)),v);
      else if (real(x) == T(0))
        v.Zero();
      else if (real(x) != T(1))
        InlineMultXV(Scaling<0,T>(real(x)),v);
    }
    else InlineMultXV(Scaling<0,std::complex<T> >(x),v); 
  }

  template <class T> 
  void InstMultXV(const T x, VectorView<T> v)
  { 
    if (v.step() == 1) {
      VectorView<T,1> vu = v.UnitView();
      DoMultXV(x,vu);
    }
    else DoMultXV(x,v); 
  }

#ifdef BLAS
#define TMV_INST_SKIP_BLAS
#ifdef TMV_INST_DOUBLE
  template <> 
  void InstMultXV(const double x, VectorView<double> v)
  {
    int n=v.size();
    if (n==0) return;
    TMVAssert(v.size()>0);
    int s=v.step();
    if (s == 0) { v[0] *= x; return; }
    double* vp = v.ptr();
    if (s < 0) vp += (n-1)*s;
    BLASNAME(dscal) (BLASV(n),BLASV(x),BLASP(vp),BLASV(s));
  }
  template <> 
  void InstMultXV(const std::complex<double> x,
      VectorView<std::complex<double> > v)
  {
    if (imag(x) == double(0)) {
      int n=v.size();
      if (n==0) return;
      TMVAssert(v.size()>0);
      int s=v.step();
      double xr = real(x);
      if (s == 0) { v[0] *= x; return; }
      std::complex<double>* vp = v.ptr();
      if (s < 0) vp += (n-1)*s;
      BLASNAME(zdscal) (BLASV(n),BLASV(xr),BLASP(vp),BLASV(s));
    } else {
      int n=v.size();
      if (n==0) return;
      TMVAssert(v.size()>0);
      int s=v.step();
      if (s == 0) { v[0] *= x; return; }
      std::complex<double>* vp = v.ptr();
      if (s < 0) vp += (n-1)*s;
      BLASNAME(zscal) (BLASV(n),BLASP(&x),BLASP(vp),BLASV(s));
    }
  }
#endif
#ifdef TMV_INST_FLOAT
  template <> 
  void InstMultXV(const float x, VectorView<float> v)
  {
    int n=v.size();
    if (n==0) return;
    TMVAssert(v.size()>0);
    int s=v.step();
    if (s == 0) { v[0] *= x; return; }
    float* vp = v.ptr();
    if (s < 0) vp += (n-1)*s;
    BLASNAME(sscal) (BLASV(n),BLASV(x),BLASP(vp),BLASV(s));
  }
  template <> 
  void InstMultXV(const std::complex<float> x,
      VectorView<std::complex<float> > v)
  {
    if (imag(x) == float(0)) {
      int n=v.size();
      if (n==0) return;
      TMVAssert(v.size()>0);
      int s=v.step();
      float xr = real(x);
      if (s == 0) { v[0] *= x; return; }
      std::complex<float>* vp = v.ptr();
      if (s < 0) vp += (n-1)*s;
      BLASNAME(csscal) (BLASV(n),BLASV(xr),BLASP(vp),BLASV(s));
    } else {
      int n=v.size();
      if (n==0) return;
      TMVAssert(v.size()>0);
      int s=v.step();
      if (s == 0) { v[0] *= x; return; }
      std::complex<float>* vp = v.ptr();
      if (s < 0) vp += (n-1)*s;
      BLASNAME(cscal) (BLASV(n),BLASP(&x),BLASP(vp),BLASV(s));
    }
  }
#endif
#endif

  template <class T, class V1, class V2>
  static void DoMultXV(const T x, const V1& v1, V2& v2)
  { 
    if (x == T(-1))
      InlineMultXV(Scaling<-1,T>(x),v1,v2); 
    else if (x == T(1))
      v2 = v1;
    else
      InlineMultXV(Scaling<0,T>(x),v1,v2); 
  }

  template <class T, class V1, class V2>
  static void DoMultXV(const std::complex<T> x, const V1& v1, V2& v2)
  {
    if (imag(x) == T(0)) {
      if (x == T(-1))
        InlineMultXV(Scaling<-1,T>(real(x)),v1,v2);
      else if (x == T(1))
        v2 = v1;
      else
        InlineMultXV(Scaling<0,T>(real(x)),v1,v2);
    }
    else InlineMultXV(Scaling<0,std::complex<T> >(x),v1,v2); 
  }

  template <class T1, bool C1, class T2>
  void InstMultXV(
      const T2 x, const ConstVectorView<T1,XX,C1>& v1, VectorView<T2> v2)
  {
    if (v2.step() == 1) 
    {
      VectorView<T2,1> v2unit = v2.UnitView();
      if (v1.step() == 1) 
        DoMultXV(x,v1.UnitView(),v2unit);
      else
        DoMultXV(x,v1,v2unit);
    } 
    else 
      if (v1.step() == 1)
        DoMultXV(x,v1.UnitView(),v2); 
      else
        DoMultXV(x,v1,v2); 
  }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
  template <> 
  void InstMultXV(const double x, 
      const ConstVectorView<double>& v1, VectorView<double> v2)
  { InstMultXV(x,v2=v1); }
#ifdef TMV_INST_MIX
  template <> 
  void InstMultXV(const std::complex<double> x,
      const ConstVectorView<double>& v1,
      VectorView<std::complex<double> > v2)
  { InstMultXV(x,v2=v1); }
#endif
  template <> 
  void InstMultXV(const std::complex<double> x,
      const ConstVectorView<std::complex<double> >& v1,
      VectorView<std::complex<double> > v2)
  { InstMultXV(x,v2=v1); }
  template <> 
  void InstMultXV(const std::complex<double> x,
      const ConstVectorView<std::complex<double>,UNKNOWN,true>& v1,
      VectorView<std::complex<double> > v2)
  { InstMultXV(x,v2=v1); }
#endif
#ifdef TMV_INST_FLOAT
  template <> 
  void InstMultXV(const float x, 
      const ConstVectorView<float>& v1, VectorView<float> v2)
  { InstMultXV(x,v2=v1); }
#ifdef TMV_INST_MIX
  template <> 
  void InstMultXV(const std::complex<float> x,
      const ConstVectorView<float>& v1,
      VectorView<std::complex<float> > v2)
  { InstMultXV(x,v2=v1); }
#endif
  template <> 
  void InstMultXV(const std::complex<float> x,
      const ConstVectorView<std::complex<float> >& v1,
      VectorView<std::complex<float> > v2)
  { InstMultXV(x,v2=v1); }
  template <> 
  void InstMultXV(const std::complex<float> x,
      const ConstVectorView<std::complex<float>,UNKNOWN,true>& v1,
      VectorView<std::complex<float> > v2)
  { InstMultXV(x,v2=v1); }
#endif
#endif


#define InstFile "TMV_MultXV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


