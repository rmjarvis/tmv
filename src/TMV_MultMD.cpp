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

//#define PRINTALGO_MD
//#include <iostream>
//#include "tmv/TMV_DiagMatrixIO.h"
//#include "tmv/TMV_MatrixIO.h"


#include "tmv/TMV_MultMD.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_AddMM.h"
#include "tmv/TMV_SumMM.h"

namespace tmv {

  // 
  //
  // MultMD
  //

  template <bool add, class T, class M1, class M2, class M3>
  static void DoMultMD(const T x, const M1& m1, const M2& m2, M3& m3)
  { 
    //std::cout<<"DoMultMD add = "<<add<<std::endl;
    //std::cout<<"m1 = "<<TypeText(m1)<<std::endl;
    //std::cout<<"m2 = "<<TypeText(m2)<<std::endl;
    //std::cout<<"m3 = "<<TypeText(m3)<<std::endl;
    if (x == T(1))
      InlineMultMM<add>(Scaling<1,T>(x),m1,m2,m3); 
    else if (x == T(-1))
      InlineMultMM<add>(Scaling<-1,T>(x),m1,m2,m3); 
    else if (x == T(0))
    { if (!add) m3.Zero(); }
    else
      InlineMultMM<add>(Scaling<0,T>(x),m1,m2,m3); 
  }

  template <bool add, class T, class M1, class M2, class M3>
  static void DoMultMD(const std::complex<T> x, 
      const M1& m1, const M2& m2, M3& m3)
  {
    //std::cout<<"DoMultMD add = "<<add<<std::endl;
    //std::cout<<"m1 = "<<TypeText(m1)<<std::endl;
    //std::cout<<"m2 = "<<TypeText(m2)<<std::endl;
    //std::cout<<"m3 = "<<TypeText(m3)<<std::endl;
    if (imag(x) == T(0)) 
    {
      if (real(x) == T(1))
        InlineMultMM<add>(Scaling<1,T>(real(x)),m1,m2,m3); 
      else if (real(x) == T(-1))
        InlineMultMM<add>(Scaling<-1,T>(real(x)),m1,m2,m3); 
      else if (real(x) == T(0))
      { if (!add) m3.Zero(); }
      else
        InlineMultMM<add>(Scaling<0,T>(real(x)),m1,m2,m3); 
    }
    else
      InlineMultMM<add>(Scaling<0,std::complex<T> >(x),m1,m2,m3); 
  }

  template <class T1, bool C1, class T2, bool C2, class T3>  
  void InstMultMM(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
      const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2, MatrixView<T3> m3)
  {
    //std::cout<<"InstMultMD\n";
    //std::cout<<"m1 = "<<TypeText(m1)<<" = "<<m1<<std::endl;
    //std::cout<<"m2 = "<<TypeText(m2)<<" = "<<m2<<std::endl;
    //std::cout<<"m3 = "<<TypeText(m3)<<" = "<<m3<<std::endl;
    if (m2.step() == 1) {
      if (m1.isrm() && m3.isrm()) {
        MatrixView<T3,UNKNOWN,1> m3rm = m3.RMView();
        DoMultMD<false>(x,m1.RMView(),m2.CMView(),m3rm);
      } else if (m3.iscm()) {
        MatrixView<T3,1> m3cm = m3.CMView();
        if (m1.iscm())
          DoMultMD<false>(x,m1.CMView(),m2.CMView(),m3cm);
        else
          DoMultMD<false>(x,m1,m2.CMView(),m3cm);
      } else {
        if (m1.iscm())
          DoMultMD<false>(x,m1.CMView(),m2.CMView(),m3);
        else
          DoMultMD<false>(x,m1,m2.CMView(),m3);
      }
    } else {
      InstMultMM(T3(1),m1,(x*m2).calc().XView(),m3);
    }
    //std::cout<<"m3 => "<<m3<<std::endl;
  }

  template <class T1, bool C1, class T2, bool C2, class T3>  
  void InstAddMultMM(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
      const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2, MatrixView<T3> m3)
  {
    //std::cout<<"InstAddMultMD\n";
    //std::cout<<"m1 = "<<TypeText(m1)<<" = "<<m1<<std::endl;
    //std::cout<<"m2 = "<<TypeText(m2)<<" = "<<m2<<std::endl;
    //std::cout<<"m3 = "<<TypeText(m3)<<" = "<<m3<<std::endl;
    if (m2.step() == 1) {
      if (m1.isrm() && m3.isrm()) {
        MatrixView<T3,UNKNOWN,1> m3rm = m3.RMView();
        DoMultMD<true>(x,m1.RMView(),m2.CMView(),m3rm);
      } else if (m3.iscm()) {
        MatrixView<T3,1> m3cm = m3.CMView();
        if (m1.iscm())
          DoMultMD<true>(x,m1.CMView(),m2.CMView(),m3cm);
        else
          DoMultMD<true>(x,m1,m2.CMView(),m3cm);
      } else {
        if (m1.iscm())
          DoMultMD<true>(x,m1.CMView(),m2.CMView(),m3);
        else
          DoMultMD<true>(x,m1,m2.CMView(),m3);
      }
    } else {
      InstAddMultMM(T3(1),m1,(x*m2).calc().XView(),m3);
    }
    //std::cout<<"m3 => "<<m3<<std::endl;
  }


#define InstFile "TMV_MultMD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


