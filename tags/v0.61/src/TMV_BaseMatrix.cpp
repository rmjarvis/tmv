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



#include "TMV_BaseMatrix.h"
#include "TMV_Vector.h"
#include "TMV_Matrix.h"
#include "TMV_Divider.h"
#include "TMV_DivImpl.h"

#include <iostream>

namespace tmv {

  template <class T> DivHelper<T>::~DivHelper() 
  { if (pdiv) { delete pdiv; pdiv=0; } }

  template <class T> void DivHelper<T>::SetupDiv() const
  { if (!(pdiv)) pdiv = new DivImpl(GetMatrix()); }

  template <class T> T DivHelper<T>::DoDet() const
  {
    TMVAssert(colsize() == rowsize());
    SetDiv();
    T det = pdiv->itsdiv->Det();
    DoneDiv();
    return det;
  }

  template <class T> RealType(T) DivHelper<T>::DoLogDet(T* sign) const
  {
    TMVAssert(colsize() == rowsize());
    SetDiv();
    RealType(T) logdet = pdiv->itsdiv->LogDet(sign);
    DoneDiv();
    return logdet;
  }

  template <class T, class T1> static inline void DoInverse1(
      const Divider<T>& div, const MatrixView<T1>& minv)
  { div.Inverse(minv); }
  template <class T> static inline void DoInverse1(
      const Divider<std::complex<T> >& , const MatrixView<T>& )
  { TMVAssert(FALSE); }

  template <class T> template <class T1> void DivHelper<T>::DoInverse(
      const MatrixView<T1>& minv) const
  {
    TMVAssert(minv.colsize() == rowsize());
    TMVAssert(minv.rowsize() == colsize());
    SetDiv();
    DoInverse1(*pdiv->itsdiv,minv);
    DoneDiv();
  }

  template <class T> void DivHelper<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == MIN(rowsize(),colsize()));
    TMVAssert(minv.rowsize() == MIN(rowsize(),colsize()));
    SetDiv();
    pdiv->itsdiv->InverseATA(minv);
    DoneDiv();
  }

  template <class T> bool DivHelper<T>::DoSingular() const
  {
    SetDiv();
    bool s = pdiv->itsdiv->Singular();
    DoneDiv();
    return s;
  }

  template <class T> RealType(T) DivHelper<T>::DoNorm2() const
  {
    SetupDiv();
    TMVAssert(DivIsSet());
    TMVAssert(pdiv->itsdt == SV);
    return pdiv->itsdiv->Norm2();
  }

  template <class T> RealType(T) DivHelper<T>::DoCondition() const
  {
    SetupDiv();
    TMVAssert(DivIsSet());
    TMVAssert(pdiv->itsdt == SV);
    return pdiv->itsdiv->Condition();
  }

  template <class T, class T1> static inline void DoRDivEq1(
      const Divider<T>& div, const MatrixView<T1>& m0)
  { div.RDivEq(m0); }

  template <class T, class T1> static inline void DoLDivEq1(
      const Divider<T>& div, const MatrixView<T1>& m0)
  { div.LDivEq(m0); }

  template <class T, class T1, class T0> static inline void DoLDiv1(
      const Divider<T>& div,
      const GenMatrix<T1>& m1, const MatrixView<T0>& m0)
  { div.LDiv(m1,m0); }

  template <class T, class T1, class T0> static inline void DoRDiv1(
      const Divider<T>& div,
      const GenMatrix<T1>& m1, const MatrixView<T0>& m0)
  { div.RDiv(m1,m0); }

  template <class T> static inline void DoRDivEq1(
      const Divider<std::complex<T> >& , const MatrixView<T>& )
  { TMVAssert(FALSE); }

  template <class T> static inline void DoLDivEq1(
      const Divider<std::complex<T> >& , const MatrixView<T>& )
  { TMVAssert(FALSE); }

  template <class T, class T1> static inline void DoLDiv1(
      const Divider<std::complex<T> >& ,
      const GenMatrix<T1>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }

  template <class T, class T1> static inline void DoRDiv1(
      const Divider<std::complex<T> >& ,
      const GenMatrix<T1>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }

  template <class T> static inline void DoLDiv1(
      const Divider<T>& ,
      const GenMatrix<std::complex<T> >& , const MatrixView<T>& )
  { TMVAssert(FALSE); }

  template <class T> static inline void DoRDiv1(
      const Divider<T>& ,
      const GenMatrix<std::complex<T> >& , const MatrixView<T>& )
  { TMVAssert(FALSE); }

  template <class T> template <class T1> void DivHelper<T>::DoLDivEq(
      const VectorView<T1>& v) const
  {
    TMVAssert(colsize() == rowsize());
    TMVAssert(colsize() == v.size());
    SetDiv();
    DoLDivEq1(*pdiv->itsdiv,ColVectorViewOf(v));
    DoneDiv();
  }

  template <class T> template <class T1> void DivHelper<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(colsize() == rowsize());
    TMVAssert(colsize() == m.colsize());
    SetDiv();
    DoLDivEq1(*pdiv->itsdiv,m);
    DoneDiv();
  }

  template <class T> template <class T1> void DivHelper<T>::DoRDivEq(
      const VectorView<T1>& v) const
  {
    TMVAssert(colsize() == rowsize());
    TMVAssert(colsize() == v.size());
    SetDiv();
    DoRDivEq1(*pdiv->itsdiv,RowVectorViewOf(v));
    DoneDiv();
  }
		
  template <class T> template <class T1> void DivHelper<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(colsize() == rowsize());
    TMVAssert(colsize() == m.rowsize());
    SetDiv();
    DoRDivEq1(*pdiv->itsdiv,m);
    DoneDiv();
  }

  template <class T> template <class T1, class T0> 
    void DivHelper<T>::DoLDiv(
	const GenVector<T1>& v1, const VectorView<T0>& v0) const
    {
      TMVAssert(rowsize() == v0.size());
      TMVAssert(colsize() == v1.size());
      SetDiv();
      DoLDiv1(*pdiv->itsdiv,ColVectorViewOf(v1),ColVectorViewOf(v0));
      DoneDiv();
    }

  template <class T> template <class T1, class T0> 
    void DivHelper<T>::DoLDiv(
	const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
    {
      TMVAssert(rowsize() == m0.colsize());
      TMVAssert(colsize() == m1.colsize());
      TMVAssert(m1.rowsize() == m0.rowsize());
      SetDiv();
      DoLDiv1(*pdiv->itsdiv,m1,m0);
      DoneDiv();
    }

  template <class T> template <class T1, class T0> 
    void DivHelper<T>::DoRDiv(
	const GenVector<T1>& v1, const VectorView<T0>& v0) const
    {
      TMVAssert(rowsize() == v1.size());
      TMVAssert(colsize() == v0.size());
      SetDiv();
      DoRDiv1(*pdiv->itsdiv,RowVectorViewOf(v1),RowVectorViewOf(v0));
      DoneDiv();
    }

  template <class T> template <class T1, class T0> 
    void DivHelper<T>::DoRDiv(
	const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
    {
      TMVAssert(rowsize() == m1.rowsize());
      TMVAssert(colsize() == m0.rowsize());
      TMVAssert(m1.colsize() == m0.colsize());
      SetDiv();
      DoRDiv1(*pdiv->itsdiv,m1,m0);
      DoneDiv();
    }

  template <class T> void DivHelper<T>::DivideInPlace() const
  {
    SetupDiv();
    SaveDiv();
    pdiv->inplace = true;
  }

  template <class T> bool DivHelper<T>::IsDivInPlace() const
  {
    SetupDiv();
    return pdiv->inplace;
  }

  template <class T> void DivHelper<T>::SaveDiv() const
  { 
    SetupDiv();
    pdiv->cache = true; 
  }

  template <class T> void DivHelper<T>::DivideUsing(DivType dt) const
  { 
    SetupDiv();
    if (dt != pdiv->itsdt) UnSetDiv();
    pdiv->itsdt = dt;
  }

  template <class T> void DivHelper<T>::SetDiv() const
  { 
    SetupDiv();
    if (!pdiv->itsdiv) {
      if (pdiv->itsdt == XXX) 
	pdiv->itsdt = (colsize() == rowsize()) ? LU : QR;
      NewDivider();
    }
  }

  template <class T> void DivHelper<T>::UnSetDiv() const
  { 
    SetupDiv();
    if (pdiv->itsdiv) {
      delete pdiv->itsdiv;
      pdiv->itsdiv = 0;
    }
  }

  template <class T> void DivHelper<T>::ReSetDiv() const
  {
    UnSetDiv(); 
    SetDiv(); 
  }

  template <class T> bool DivHelper<T>::DivIsSet() const 
  {
    SetupDiv();
    return pdiv->itsdiv; 
  }

  template <class T> void DivHelper<T>::DoneDiv() const
  { 
    TMVAssert(pdiv);
    if (!pdiv->cache) UnSetDiv(); 
  }

  template <class T> const Divider<T>* DivHelper<T>::GetDiv() const 
  {
    SetupDiv();
    return pdiv->itsdiv; 
  }

  template <class T> void DivHelper<T>::SetDiv(Divider<T>* d) const 
  {
    SetupDiv();
    if (pdiv->itsdiv) delete pdiv->itsdiv;
    pdiv->itsdiv = d; 
  }

  template <class T> DivType DivHelper<T>::GetDivType() const 
  {
    SetupDiv();
    return pdiv->itsdt; 
  }

  template <class T> bool DivHelper<T>::CheckDecomp(std::ostream* fout) const
  {
    TMVAssert(pdiv);
    TMVAssert(pdiv->itsdiv);
    return pdiv->itsdiv->CheckDecomp(GetMatrix(),fout);
  }

  template <class T> bool DivHelper<T>::CheckDecomp(
      const BaseMatrix<T>& m2, std::ostream* fout) const
  {
    TMVAssert(pdiv);
    TMVAssert(pdiv->itsdiv);
    return pdiv->itsdiv->CheckDecomp(m2,fout);
  }

#define InstFile "TMV_BaseMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


