///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
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


#ifndef TMV_TriBandArithFunc_H
#define TMV_TriBandArithFunc_H

#define CT std::complex<T>

#include "TMV_BaseTriMatrix.h"
#include "TMV_BaseBandMatrix.h"
#include "TMV_BandMatrixArithFunc.h"

namespace tmv {

  template <class T, class T2> inline void AddMM(
      const T x2, const GenUpperTriMatrix<T2>& m2,
      const BandMatrixView<T>& m1)
  {
    if (m2.isunit()) {
      if (m2.size() > 1)
	AddMM(x2,BandMatrixViewOf(m2.OffDiag()),m1.Diags(1,m2.size()));
      m1.diag().AddToAll(x2);
    } else {
      AddMM(x2,BandMatrixViewOf(m2),m1);
    }
  }

  template <class T, class T2> inline void AddMM(
      const T x2, const GenLowerTriMatrix<T2>& m2,
      const BandMatrixView<T>& m1)
  {
    if (m2.isunit()) {
      if (m2.size() > 1)
	AddMM(x2,BandMatrixViewOf(m2.OffDiag()),m1.Diags(-int(m2.size())+1,0));
      m1.diag().AddToAll(x2);
    } else {
      AddMM(x2,BandMatrixViewOf(m2),m1);
    }
  }

  template <bool add, class T, class T1, class T2> inline void MultMM(
      const T x, const GenUpperTriMatrix<T1>& m1,
      const GenBandMatrix<T2>& m2, const BandMatrixView<T>& m0)
  { 
    if (m1.isunit()) {
      UpperTriMatrix<T1,NonUnitDiag,RowMajor> m1x = m1;
      MultMM<add>(x,BandMatrixViewOf(m1x),m2,m0);
    } else {
      MultMM<add>(x,BandMatrixViewOf(m1),m2,m0);
    }
  }
  template <bool add, class T, class T1, class T2> inline void MultMM(
      const T x, const GenBandMatrix<T1>& m1,
      const GenUpperTriMatrix<T2>& m2, const BandMatrixView<T>& m0)
  { 
    if (m2.isunit()) {
      UpperTriMatrix<T2,NonUnitDiag,RowMajor> m2x = m2;
      MultMM<add>(x,m1,BandMatrixViewOf(m2x),m0);
    } else {
      MultMM<add>(x,m1,BandMatrixViewOf(m2),m0);
    }
  }
  template <bool add, class T, class T1, class T2> inline void MultMM(
      const T x, const GenLowerTriMatrix<T1>& m1,
      const GenBandMatrix<T2>& m2, const BandMatrixView<T>& m0)
  { 
    if (m1.isunit()) {
      LowerTriMatrix<T1,NonUnitDiag,RowMajor> m1x = m1;
      MultMM<add>(x,BandMatrixViewOf(m1x),m2,m0);
    } else {
      MultMM<add>(x,BandMatrixViewOf(m1),m2,m0);
    }
  }
  template <bool add, class T, class T1, class T2> inline void MultMM(
      const T x, const GenBandMatrix<T1>& m1,
      const GenLowerTriMatrix<T2>& m2, const BandMatrixView<T>& m0)
  { 
    if (m2.isunit()) {
      LowerTriMatrix<T2,NonUnitDiag,RowMajor> m2x = m2;
      MultMM<add>(x,m1,BandMatrixViewOf(m2x),m0);
    } else {
      MultMM<add>(x,m1,BandMatrixViewOf(m2),m0);
    }
  }

  template <class T, class T2> inline void AddMM(
      const T x2, const GenUpperTriMatrix<T2>& m2,
      const BandMatrixView<CT>& m1)
  { AddMM(CT(x2),m2,m1); }
  template <class T, class T2> inline void AddMM(
      const T x2, const GenLowerTriMatrix<T2>& m2,
      const BandMatrixView<CT>& m1)
  { AddMM(CT(x2),m2,m1); }

  template <class T, class T2> inline void AddMM(
      const CT , const GenUpperTriMatrix<T2>& ,
      const BandMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class T2> inline void AddMM(
      const CT , const GenLowerTriMatrix<T2>& ,
      const BandMatrixView<T>& )
  { TMVAssert(FALSE); }

  template <bool add, class T, class T1, class T2> inline void MultMM(
      const T x, const GenUpperTriMatrix<T1>& m1,
      const GenBandMatrix<T2>& m2, const BandMatrixView<CT>& m0)
  { MultMM<add>(CT(x),m1,m2,m0); }
  template <bool add, class T, class T1, class T2> inline void MultMM(
      const T x, const GenBandMatrix<T1>& m1,
      const GenUpperTriMatrix<T2>& m2, const BandMatrixView<CT>& m0)
  { MultMM<add>(CT(x),m1,m2,m0); }
  template <bool add, class T, class T1, class T2> inline void MultMM(
      const T x, const GenLowerTriMatrix<T1>& m1,
      const GenBandMatrix<T2>& m2, const BandMatrixView<CT>& m0)
  { MultMM<add>(CT(x),m1,m2,m0); }
  template <bool add, class T, class T1, class T2> inline void MultMM(
      const T x, const GenBandMatrix<T1>& m1,
      const GenLowerTriMatrix<T2>& m2, const BandMatrixView<CT>& m0)
  { MultMM<add>(CT(x),m1,m2,m0); }

  template <bool add, class T, class T1, class T2> inline void MultMM(
      const CT , const GenUpperTriMatrix<T1>& ,
      const GenBandMatrix<T2>& , const BandMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <bool add, class T, class T1, class T2> inline void MultMM(
      const CT , const GenBandMatrix<T1>& ,
      const GenUpperTriMatrix<T2>& , const BandMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <bool add, class T, class T1, class T2> inline void MultMM(
      const CT , const GenLowerTriMatrix<T1>& ,
      const GenBandMatrix<T2>& , const BandMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <bool add, class T, class T1, class T2> inline void MultMM(
      const CT , const GenBandMatrix<T1>& ,
      const GenLowerTriMatrix<T2>& , const BandMatrixView<T>& )
  { TMVAssert(FALSE); }

}

#undef CT

#endif
