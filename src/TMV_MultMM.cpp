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

  // Defined in TMV_MultMM_CCC.cpp
  template <bool add, class T1, bool C1, class T2, bool C2, class T3>
  void DoInstMultMM(const T3 x,
      const ConstMatrixView<T1,1,UNKNOWN,C1>& m1,
      const ConstMatrixView<T2,1,UNKNOWN,C2>& m2, MatrixView<T3,1> m3);

  // Defined in TMV_MultMM_CRC.cpp
  template <bool add, class T1, bool C1, class T2, bool C2, class T3>
  void DoInstMultMM(const T3 x,
      const ConstMatrixView<T1,1,UNKNOWN,C1>& m1,
      const ConstMatrixView<T2,UNKNOWN,1,C2>& m2, MatrixView<T3,1> m3);

  // Defined in TMV_MultMM_RCC.cpp
  template <bool add, class T1, bool C1, class T2, bool C2, class T3>
  void DoInstMultMM(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,1,C1>& m1,
      const ConstMatrixView<T2,1,UNKNOWN,C2>& m2, MatrixView<T3,1> m3);

  // Defined in TMV_MultMM_RRC.cpp
  template <bool add, class T1, bool C1, class T2, bool C2, class T3>
  void DoInstMultMM(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,1,C1>& m1,
      const ConstMatrixView<T2,UNKNOWN,1,C2>& m2, MatrixView<T3,1> m3);

  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstMultMM(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
      const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3)
  {
    if (m3.isrm() && !m3.iscm())
      InstMultMM(x,m2.Transpose(),m1.Transpose(),m3.Transpose());
    else if (!m3.iscm())
    {
      Matrix<T3,ColMajor> m3c(m3.colsize(),m3.rowsize());
      InstMultMM(T3(1),m1,m2,m3c.XView());
      m3 = x * m3c;
    }
    else 
    {
      MatrixView<T3,1> m3cm = m3.CMView();
      if (m1.isrm())
      {
        if (m2.isrm())
          DoInstMultMM<false>(x,m1.RMView(),m2.RMView(),m3cm);
        else if (m2.iscm())
          DoInstMultMM<false>(x,m1.RMView(),m2.CMView(),m3cm);
        else 
        {
          const Matrix<T3,ColMajor> m2c = x*m2;
          DoInstMultMM<false>(T3(1),m1.RMView(),m2c.View(),m3cm);
        }
      }
      else if (m1.iscm())
      {
        if (m2.isrm())
          DoInstMultMM<false>(x,m1.CMView(),m2.RMView(),m3cm);
        else if (m2.iscm())
          DoInstMultMM<false>(x,m1.CMView(),m2.CMView(),m3cm);
        else 
        {
          const Matrix<T3,ColMajor> m2c = x*m2;
          DoInstMultMM<false>(T3(1),m1.CMView(),m2c.View(),m3cm);
        }
      }
      else 
      {
        if (m2.isrm())
        {
          const Matrix<T3,ColMajor> m1c = x*m1;
          DoInstMultMM<false>(T3(1),m1c.View(),m2.RMView(),m3cm);
        }
        else if (m2.iscm())
        {
          const Matrix<T3,RowMajor> m1c = x*m1;
          DoInstMultMM<false>(T3(1),m1c.View(),m2.CMView(),m3cm);
        }
        else 
        {
          const Matrix<T3,RowMajor> m1c = x*m1;
          const Matrix<T2,ColMajor> m2c = m2;
          DoInstMultMM<false>(T3(1),m1c.View(),m2c.View(),m3cm);
        }
      }
    }
  }

  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstAddMultMM(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
      const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3)
  {
    if (m3.isrm() && !m3.iscm())
      InstAddMultMM(x,m2.Transpose(),m1.Transpose(),m3.Transpose());
    else if (!m3.iscm())
    {
      Matrix<T3,ColMajor> m3c(m3.colsize(),m3.rowsize());
      InstMultMM(T3(1),m1,m2,m3c.XView());
      m3 += x * m3c;
    }
    else 
    {
      MatrixView<T3,1> m3cm = m3.CMView();
      if (m1.isrm())
      {
        if (m2.isrm())
          DoInstMultMM<true>(x,m1.RMView(),m2.RMView(),m3cm);
        else if (m2.iscm())
          DoInstMultMM<true>(x,m1.RMView(),m2.CMView(),m3cm);
        else 
        {
          const Matrix<T3,ColMajor> m2c = x*m2;
          DoInstMultMM<true>(T3(1),m1.RMView(),m2c.View(),m3cm);
        }
      }
      else if (m1.iscm())
      {
        if (m2.isrm())
          DoInstMultMM<true>(x,m1.CMView(),m2.RMView(),m3cm);
        else if (m2.iscm())
          DoInstMultMM<true>(x,m1.CMView(),m2.CMView(),m3cm);
        else 
        {
          const Matrix<T3,ColMajor> m2c = x*m2;
          DoInstMultMM<true>(T3(1),m1.CMView(),m2c.View(),m3cm);
        }
      }
      else 
      {
        if (m2.isrm())
        {
          const Matrix<T3,ColMajor> m1c = x*m1;
          DoInstMultMM<true>(T3(1),m1c.View(),m2.RMView(),m3cm);
        }
        else if (m2.iscm())
        {
          const Matrix<T3,RowMajor> m1c = x*m1;
          DoInstMultMM<true>(T3(1),m1c.View(),m2.CMView(),m3cm);
        }
        else 
        {
          const Matrix<T3,RowMajor> m1c = x*m1;
          const Matrix<T2,ColMajor> m2c = m2;
          DoInstMultMM<true>(T3(1),m1c.View(),m2c.View(),m3cm);
        }
      }
    }
  }

#define InstFile "TMV_MultMM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


