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

#include "tmv/TMV_AddUU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SumVV.h"

namespace tmv {

  // 
  // AddMM
  //

  template <class T, class M1, class M2> 
  static void DoAddMM(const T x, const M1& m1, M2& m2)
  {
    if (x == T(1))
      InlineAddMM(Scaling<1,T>(1),m1,m2);
    else if (x == T(-1))
      InlineAddMM(Scaling<-1,T>(-1),m1,m2);
    else if (x != T(0))
      InlineAddMM(Scaling<0,T>(x),m1,m2);
  }

  template <class T, class M1, class M2> 
  static void DoAddMM(const std::complex<T> x, const M1& m1, M2& m2)
  {
    if (imag(x) == T(0)) 
    {
      if (real(x) == T(1))
        InlineAddMM(Scaling<1,T>(1),m1,m2);
      else if (real(x) == T(-1))
        InlineAddMM(Scaling<-1,T>(-1),m1,m2);
      else if (real(x) != T(0))
        InlineAddMM(Scaling<0,T>(real(x)),m1,m2);
    }
    else 
      InlineAddMM(Scaling<0,std::complex<T> >(x),m1,m2);
  }

  template <class T1, DiagType D1, bool C1, class T2>
  void InstAddMM(const T2 x,
      const ConstUpperTriMatrixView<T1,D1,UNKNOWN,UNKNOWN,C1>& m1,
      UpperTriMatrixView<T2,NonUnitDiag> m2)
  {
    if (m2.isrm())
    {
      UpperTriMatrixView<T2,NonUnitDiag,UNKNOWN,1> m2rm = m2;
      if (m1.isrm())
        DoAddMM(x,m1.RMView(),m2rm);
      else if (m2.iscm())
        DoAddMM(x,m1.CMView(),m2rm);
      else
        DoAddMM(x,m1,m2rm);
    }
    else if (m2.iscm())
    {
      UpperTriMatrixView<T2,NonUnitDiag,1> m2cm = m2;
      if (m1.isrm())
        DoAddMM(x,m1.RMView(),m2cm);
      else if (m2.iscm())
        DoAddMM(x,m1.CMView(),m2cm);
      else
        DoAddMM(x,m1,m2cm);
    }
    else 
    {
      if (m1.isrm())
        DoAddMM(x,m1.RMView(),m2);
      else if (m2.iscm())
        DoAddMM(x,m1.CMView(),m2);
      else
        DoAddMM(x,m1,m2);
    }
  }

  template <class T1, DiagType D1, bool C1, class T2>
  void InstAddMM(const T2 x,
      const ConstLowerTriMatrixView<T1,D1,UNKNOWN,UNKNOWN,C1>& m1,
      LowerTriMatrixView<T2,NonUnitDiag> m2)
  { InstAddMM(x,m1.Transpose(),m2.Transpose()); }

#define InstFile "TMV_AddUU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


