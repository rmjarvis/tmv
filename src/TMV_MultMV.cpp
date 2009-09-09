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
#include "tmv/TMV_MultMV.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_ProdXV.h"
#include "tmv/TMV_SumVV.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_ProdXM.h"

#ifdef BLAS
#include "TMV_MultMV_Blas.h"
#endif

namespace tmv {

  // 
  //
  // MultMV
  //

  template <bool add, class T, class M1, class V2, class V3>
  static void DoMultMV(const T x, const M1& m1, const V2& v2, V3& v3)
  {
    // Check for non-unit step and x != 1, and do the necessary copies here,
    // rather than in the InlineMultMV function.  
    // This is faster to compile, since it keeps the InlineMultMV
    // algo path to the ones that have vstep == 1.

    typedef RealType(T) RT;
    const Scaling<1,RT> one;

    if (x == RT(0)) 
    { Maybe<!add>::zero(v3); }
    else if (v2.step() != 1) 
    {
      Vector<T> v2c = x*v2;
      if (v3.step() != 1) 
      {
        Vector<T> v3c(v3.size());
        VectorView<T,1> v3u = v3c.UnitView();
        InlineMultMV<false>(one,m1,v2c.UnitView(),v3u);
        Maybe<add>::add(v3,v3c);
      }
      else 
      {
        VectorView<T,1> v3u = v3.UnitView();
        InlineMultMV<add>(one,m1,v2c.UnitView(),v3u);
      }
    } 
    else if (v3.step() != 1)
    {
      Vector<T> v3c(v3.size());
      VectorView<T,1> v3u = v3c.UnitView();
      InlineMultMV<false>(one,m1,v2.UnitView(),v3u);
      Maybe<add>::add(v3,x*v3c);
    }
    else 
    {
      VectorView<T,1> v3u = v3.UnitView();
      if (x == RT(1))
        InlineMultMV<add>(one,m1,v2.UnitView(),v3u);
      else if (x == RT(-1))
        InlineMultMV<add>(Scaling<-1,RT>(),m1,v2.UnitView(),v3u);
      else if (TMV_IMAG(x) == RT(0))
        InlineMultMV<add>(Scaling<0,RT>(TMV_REAL(x)),m1,v2.UnitView(),v3u);
      else 
        InlineMultMV<add>(Scaling<0,T>(x),m1,v2.UnitView(),v3u);
    }
  }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
  template <bool add>
  void DoMultMV(const double x,
      const ConstMatrixView<double,1>& m1,
      const ConstVectorView<double>& v2, VectorView<double> v3)
  { BlasMultMV(x,m1,v2,add?1:0,v3); }

  template <bool add, class T1, bool C1, class T2, bool C2>
  void DoMultMV(const std::complex<double>  x,
      const ConstMatrixView<T1,1,UNKNOWN,C1>& m1,
      const ConstVectorView<T2,UNKNOWN,C2>& v2,
      VectorView<std::complex<double> > v3)
  { BlasMultMV(x,m1,v2,add?1:0,v3); }

  template <bool add>
  void DoMultMV(const double x,
      const ConstMatrixView<double,UNKNOWN,1>& m1,
      const ConstVectorView<double>& v2, VectorView<double> v3)
  { BlasMultMV(x,m1,v2,add?1:0,v3); }

  template <bool add, class T1, bool C1, class T2, bool C2>
  void DoMultMV(const std::complex<double>  x,
      const ConstMatrixView<T1,UNKNOWN,1,C1>& m1,
      const ConstVectorView<T2,UNKNOWN,C2>& v2,
      VectorView<std::complex<double> > v3)
  { BlasMultMV(x,m1,v2,add?1:0,v3); }
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
  template <bool add>
  void DoMultMV(const float x,
      const ConstMatrixView<float,UNKNOWN,1>& m1,
      const ConstVectorView<float>& v2, VectorView<float> v3)
  { BlasMultMV(x,m1,v2,add?1:0,v3); }

  template <bool add, class T1, bool C1, class T2, bool C2>
  void DoMultMV(const std::complex<float>  x,
      const ConstMatrixView<T1,UNKNOWN,1,C1>& m1,
      const ConstVectorView<T2,UNKNOWN,C2>& v2,
      VectorView<std::complex<float> > v3)
  { BlasMultMV(x,m1,v2,add?1:0,v3); }

  template <bool add>
  void DoMultMV(const float x,
      const ConstMatrixView<float,1>& m1,
      const ConstVectorView<float>& v2, VectorView<float> v3)
  { BlasMultMV(x,m1,v2,add?1:0,v3); }

  template <bool add, class T1, bool C1, class T2, bool C2>
  void DoMultMV(const std::complex<float>  x,
      const ConstMatrixView<T1,1,UNKNOWN,C1>& m1,
      const ConstVectorView<T2,UNKNOWN,C2>& v2,
      VectorView<std::complex<float> > v3)
  { BlasMultMV(x,m1,v2,add?1:0,v3); }
#endif // FLOAT
#endif // BLAS

  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstMultMV(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
      const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3)
  {
    if (m1.isrm())
      DoMultMV<false>(x,m1.RMView(),v2,v3);
    else if (m1.iscm())
      DoMultMV<false>(x,m1.CMView(),v2,v3);
    else 
      DoMultMV<false>(T3(1),(x*m1).calc().View(),v2,v3);
  }
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstAddMultMV(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
      const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3)
  {
    if (m1.isrm())
      DoMultMV<true>(x,m1.RMView(),v2,v3);
    else if (m1.iscm())
      DoMultMV<true>(x,m1.CMView(),v2,v3);
    else 
      DoMultMV<true>(T3(1),(x*m1).calc().View(),v2,v3);
  }


#define InstFile "TMV_MultMV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


