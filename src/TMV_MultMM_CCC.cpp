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

#include "tmv/TMV_MultMM.h"
#include "TMV_MultMM_Blas.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_SumMM.h"

namespace tmv {

  template <bool add, class T1, bool C1, class T2, bool C2, class T3>
  void DoMultMM(const T3 x,
      const ConstMatrixView<T1,1,UNKNOWN,C1>& m1,
      const ConstMatrixView<T2,1,UNKNOWN,C2>& m2, MatrixView<T3,1> m3)
  {
    //std::cout<<"Start DoMultMM: x = "<<x<<std::endl;
    //std::cout<<"m1 = "<<TypeText(m1)<<std::endl;
    //std::cout<<"m2 = "<<TypeText(m2)<<std::endl;
    //std::cout<<"m3 = "<<TypeText(m3)<<std::endl;
    typedef RealType(T3) RT;
    if (x == RT(0))
      Maybe<!add>::zero(m3); 
    else if (x == RT(1))
      InlineMultMM<add>(Scaling<1,RT>(),m1,m2,m3);
    else if (x == RT(-1))
      InlineMultMM<add>(Scaling<-1,RT>(),m1,m2,m3);
    else if (TMV_IMAG(x) == RT(0))
      InlineMultMM<add>(Scaling<0,RT>(TMV_REAL(x)),m1,m2,m3);
    else
      InlineMultMM<add>(Scaling<0,T3>(x),m1,m2,m3);
  }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
  template <bool add>
  static void DoMultMM(const double x,
      const ConstMatrixView<double,1>& m1,
      const ConstMatrixView<double,1>& m2, MatrixView<double,1> m3)
  { BlasMultMM(x,m1,m2,add?1:0,m3); }
  template <bool add, class T1, bool C1, class T2, bool C2>
  static void DoMultMM(const std::complex<double> x,
      const ConstMatrixView<T1,1,UNKNOWN,C1>& m1,
      const ConstMatrixView<T2,1,UNKNOWN,C2>& m2,
      MatrixView<std::complex<double>,1> m3)
  { BlasMultMM(x,m1,m2,add?1:0,m3); }
#endif // TMV_INST_DOUBLE
#ifdef TMV_INST_FLOAT
  template <bool add>
  static void DoMultMM(const float x,
      const ConstMatrixView<float,1>& m1,
      const ConstMatrixView<float,1>& m2, MatrixView<float,1> m3)
  { BlasMultMM(x,m1,m2,add?1:0,m3); }
  template <bool add, class T1, bool C1, class T2, bool C2>
  static void DoMultMM(const std::complex<float> x,
      const ConstMatrixView<T1,1,UNKNOWN,C1>& m1,
      const ConstMatrixView<T2,1,UNKNOWN,C2>& m2,
      MatrixView<std::complex<float>,1> m3)
  { BlasMultMM(x,m1,m2,add?1:0,m3); }
#endif // TMV_INST_FLOAT
#endif // BLAS

  template <bool add, class T1, bool C1, class T2, bool C2, class T3>
  void DoInstMultMM(const T3 x,
      const ConstMatrixView<T1,1,UNKNOWN,C1>& m1,
      const ConstMatrixView<T2,1,UNKNOWN,C2>& m2, MatrixView<T3,1> m3)
  { DoMultMM<add>(x,m1,m2,m3); }

#define InstFile "TMV_MultMM_CCC.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


