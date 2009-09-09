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
#include "tmv/TMV_AddVV.h"
#include "tmv/TMV_Vector.h"

namespace tmv {

  // 
  // AddVV
  //

  //
  // v2 += x * v1
  //

  template <class T, class V1, class V2> 
  static void DoAddVV(const T x, const V1& v1, V2& v2)
  {
    if (V2::vstep != 1 && v2.step() == 1) 
    {
      typename V2::unitview_type v2u = v2.UnitView();
      DoAddVV(x,v1,v2u);
    } 
    else if (V1::vstep != 1 && v1.step() == 1) 
      DoAddVV(x,v1.UnitView(),v2);
    else if (x == T(1))
      InlineAddVV(Scaling<1,T>(x),v1,v2);
    else if (x == T(-1))
      InlineAddVV(Scaling<-1,T>(x),v1,v2);
    else if (x != T(0))
      InlineAddVV(Scaling<0,T>(x),v1,v2);
  }

  template <class T, class V1, class V2> 
  static void DoAddVV(const std::complex<T> x, const V1& v1, V2& v2)
  {
    if (imag(x) == T(0)) 
    {
      if (real(x) == T(1))
        InlineAddVV(Scaling<1,T>(real(x)),v1,v2);
      else if (real(x) == T(-1))
        InlineAddVV(Scaling<-1,T>(real(x)),v1,v2);
      else if (real(x) != T(0))
        InlineAddVV(Scaling<0,T>(real(x)),v1,v2);
    }
    else
      InlineAddVV(Scaling<0,std::complex<T> >(x),v1,v2);
  }

  template <class T1, bool C1, class T2>
  void InstAddVV(const T2 x,
      const ConstVectorView<T1,UNKNOWN,C1>& v1, VectorView<T2> v2)
  {
    if (v2.step() == 1) 
    {
      VectorView<T2,1> v2u = v2;
      if (v1.step() == 1) 
        DoAddVV(x,v1.UnitView(),v2u);
      else
        DoAddVV(x,v1,v2u);
    }
    else 
    {
      if (v1.step() == 1) 
        DoAddVV(x,v1.UnitView(),v2);
      else
        DoAddVV(x,v1,v2);
    }
  }

#ifdef BLAS
#define TMV_INST_SKIP_BLAS
#ifdef TMV_INST_DOUBLE
  template <> void InstAddVV(const double x,
      const ConstVectorView<double>& v1, VectorView<double> v2)
  { 
    TMVAssert(v1.size() == v1.size()); 
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    const double* v1p = v1.cptr();
    if (s1<0) v1p += (n-1)*s1;
    double* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*s2;
    BLASNAME(daxpy) (BLASV(n),BLASV(x),BLASP(v1p),BLASV(s1),
        BLASP(v2p),BLASV(s2));
  }
  template <> void InstAddVV(const std::complex<double> x, 
      const ConstVectorView<std::complex<double> >& v1, 
      VectorView<std::complex<double> > v2)
  {
    TMVAssert(v1.size() == v1.size()); 
    if (imag(x) == 0.) {
      if (v1.step() == 1 && v2.step() == 1) {
        InstAddVV(real(x),v1.Flatten().XView(),v2.Flatten().XView());
      }
      else {
        InstAddVV(real(x),v1.Real(),v2.Real());
        InstAddVV(real(x),v1.Imag(),v2.Imag());
      }
    } else {
      int n=v2.size();
      int s1=v1.step();
      int s2=v2.step();
      const std::complex<double>* v1p = v1.cptr();
      if (s1<0) v1p += (n-1)*s1;
      std::complex<double>* v2p = v2.ptr();
      if (s2<0) v2p += (n-1)*s2;
      BLASNAME(zaxpy) (BLASV(n),BLASP(&x),BLASP(v1p),BLASV(s1),
          BLASP(v2p),BLASV(s2));
    }
  }
#ifdef TMV_INST_MIX
  template <> void InstAddVV(const std::complex<double> x, 
      const ConstVectorView<double>& v1, 
      VectorView<std::complex<double> > v2)
  {
    TMVAssert(v1.size() == v1.size()); 
    double xr = real(x);
    if (xr != 0.) InstAddVV(xr,v1,v2.Real());
    double xi = imag(x);
    if (xi != 0.) InstAddVV(xi,v1,v2.Imag());
  }
#endif
#endif
#ifdef TMV_INST_FLOAT
  template <> void InstAddVV(const float x,
      const ConstVectorView<float>& v1, VectorView<float> v2)
  {
    TMVAssert(v1.size() == v1.size()); 
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    const float* v1p = v1.cptr();
    if (s1<0) v1p += (n-1)*s1;
    float* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*s2;
    BLASNAME(saxpy) (BLASV(n),BLASV(x),BLASP(v1p),BLASV(s1),
        BLASP(v2p),BLASV(s2));
  }
  template <> void InstAddVV(const std::complex<float> x, 
      const ConstVectorView<std::complex<float> >& v1, 
      VectorView<std::complex<float> > v2)
  {
    TMVAssert(v1.size() == v1.size()); 
    if (imag(x) == 0.F) {
      if (v1.step() == 1 && v2.step() == 1) 
        InstAddVV(real(x),v1.Flatten().XView(),v2.Flatten().XView());
      else {
        InstAddVV(real(x),v1.Real(),v2.Real());
        InstAddVV(real(x),v1.Imag(),v2.Imag());
      }
    } else {
      int n=v2.size();
      int s1=v1.step();
      int s2=v2.step();
      const std::complex<float>* v1p = v1.cptr();
      if (s1<0) v1p += (n-1)*s1;
      std::complex<float>* v2p = v2.ptr();
      if (s2<0) v2p += (n-1)*s2;
      BLASNAME(caxpy) (BLASV(n),BLASP(&x),BLASP(v1p),BLASV(s1),
          BLASP(v2p),BLASV(s2));
    }
  }
#ifdef TMV_INST_MIX
  template <> void InstAddVV(const std::complex<float> x, 
      const ConstVectorView<float>& v1, 
      VectorView<std::complex<float> > v2)
  {
    TMVAssert(v1.size() == v1.size()); 
    float xr = real(x);
    if (xr != 0.F) InstAddVV(xr,v1,v2.Real());
    float xi = imag(x);
    if (xi != 0.F) InstAddVV(xi,v1,v2.Imag());
  }
#endif
#endif
#endif // BLAS

#define InstFile "TMV_AddVV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


