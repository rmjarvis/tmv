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



#include "TMV_TriDiv.h"
#include "TMV_TriMatrix.h"
#include "TMV_Matrix.h"
#include "TMV_Vector.h"
#include "TMV_DiagMatrix.h"
#include "TMV_TriMatrixArith.h"

namespace tmv {

  template <class T> QuotXU<T,T> GenUpperTriMatrix<T>::QInverse() const
  { return QuotXU<T,T>(T(1),*this); }

  template <class T> QuotXL<T,T> GenLowerTriMatrix<T>::QInverse() const
  { return QuotXL<T,T>(T(1),*this); }

  template <class T> template <class T2> void GenUpperTriMatrix<T>::DoLDivEq(
      const VectorView<T2>& v) const
  { Tri_LDivEq(*this,v); }

  template <class T> template <class T2> void GenUpperTriMatrix<T>::DoLDivEq(
      const MatrixView<T2>& m) const
  { Tri_LDivEq(*this,m); }

  template <class T> template <class T2> void GenUpperTriMatrix<T>::DoLDivEq(
      const UpperTriMatrixView<T2>& m) const
  { Tri_LDivEq(*this,m); }

  template <class T> template <class T1, class T2> 
    void GenUpperTriMatrix<T>::DoLDiv(
	const GenVector<T1>& v1, const VectorView<T2>& v2) const
    {
      if (SameStorage(*this,v1)) {
	Vector<T2> temp = v1;
	Tri_LDivEq(*this,temp.View());
	v2 = temp;
      } else {
	Tri_LDivEq(*this,v2=v1);
      }
    }

  template <class T> template <class T1, class T2> 
    void GenUpperTriMatrix<T>::DoLDiv(
	const GenMatrix<T1>& m1, const MatrixView<T2>& m2) const
    {
      if (SameStorage(*this,m2)) {
	if (m2.isrm()) {
	  Matrix<T2,RowMajor> temp = m1;
	  Tri_LDivEq(*this,temp.View());
	  m2 = temp;
	} else {
	  Matrix<T2,ColMajor> temp = m1;
	  Tri_LDivEq(*this,temp.View());
	  m2 = temp;
	}
      } else {
	Tri_LDivEq(*this,m2=m1);
      }
    }

  template <class T> template <class T1, class T2> 
    void GenUpperTriMatrix<T>::DoLDiv(
	const GenUpperTriMatrix<T1>& m1, const UpperTriMatrixView<T2>& m2) const
    {
      if (SameStorage(*this,m2)) {
	if (m2.isrm()) {
	  UpperTriMatrix<T2,NonUnitDiag,RowMajor> temp = m1;
	  Tri_LDivEq(*this,temp.View());
	  m2 = temp;
	} else {
	  UpperTriMatrix<T2,NonUnitDiag,ColMajor> temp = m1;
	  Tri_LDivEq(*this,temp.View());
	  m2 = temp;
	}
      } else {
	Tri_LDivEq(*this,m2=m1);
      }
    }

  template <class T> template <class T2> void GenLowerTriMatrix<T>::DoLDivEq(
      const VectorView<T2>& v) const
  { Tri_LDivEq(*this,v); }

  template <class T> template <class T2> void GenLowerTriMatrix<T>::DoLDivEq(
      const MatrixView<T2>& m) const
  { Tri_LDivEq(*this,m); }

  template <class T> template <class T2> void GenLowerTriMatrix<T>::DoLDivEq(
      const LowerTriMatrixView<T2>& m) const
  { Tri_LDivEq(*this,m); }

  template <class T> template <class T1, class T2> 
    void GenLowerTriMatrix<T>::DoLDiv(
	const GenVector<T1>& v1, const VectorView<T2>& v2) const
    {
      if (SameStorage(*this,v2)) {
	Vector<T2> temp = v1;
	Tri_LDivEq(*this,temp.View());
	v2 = temp;
      } else {
	Tri_LDivEq(*this,v2=v1);
      }
    }

  template <class T> template <class T1, class T2> 
    void GenLowerTriMatrix<T>::DoLDiv(
	const GenMatrix<T1>& m1, const MatrixView<T2>& m2) const
    {
      if (SameStorage(*this,m2)) {
	if (m2.isrm()) {
	  Matrix<T2,RowMajor> temp = m1;
	  Tri_LDivEq(*this,temp.View());
	  m2 = temp;
	} else {
	  Matrix<T2,ColMajor> temp = m1;
	  Tri_LDivEq(*this,temp.View());
	  m2 = temp;
	}
      } else {
	Tri_LDivEq(*this,m2=m1);
      }
    }

  template <class T> template <class T1, class T2> 
    void GenLowerTriMatrix<T>::DoLDiv(
	const GenLowerTriMatrix<T1>& m1, const LowerTriMatrixView<T2>& m2) const
    {
      if (SameStorage(*this,m2)) {
	if (m2.isrm()) {
	  LowerTriMatrix<T2,NonUnitDiag,RowMajor> temp = m1;
	  Tri_LDivEq(*this,temp.View());
	  m2 = temp;
	} else {
	  LowerTriMatrix<T2,NonUnitDiag,ColMajor> temp = m1;
	  Tri_LDivEq(*this,temp.View());
	  m2 = temp;
	}
      } else {
	Tri_LDivEq(*this,m2=m1);
      }
    }

  template <class T> T GenUpperTriMatrix<T>::Det() const
  { return isunit() ? T(1) : DiagMatrixViewOf(diag()).Det(); }                  

  template <class T> RealType(T) GenUpperTriMatrix<T>::LogDet(T* s) const
  { 
    if (isunit()) {
      if (s) *s = T(1);
      return RealType(T)(0);
    } else {
      return DiagMatrixViewOf(diag()).LogDet(s); 
    }                  
  }

  template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::InvertSelf() const
    {
      Tri_Inverse(*this);
      return *this;
    }

  template <class T> template <class T1> void GenUpperTriMatrix<T>::DoInverse(
      const MatrixView<T1>& minv) const
  {
    bool ss = SameStorage(*this,minv);

    if (!ss) minv.Zero();

    UpperTriMatrixView<T1> U = UpperTriMatrixViewOf(minv);
    U = *this;
    U.InvertSelf();

    if (ss && minv.rowsize() > 1)
      LowerTriMatrixViewOf(minv).OffDiag().Zero();
  }

  template <class T> template <class T1> void GenUpperTriMatrix<T>::DoInverse(
      const UpperTriMatrixView<T1>& minv) const
  {
    minv = *this;
    minv.InvertSelf();
  }

  template <class T> void GenUpperTriMatrix<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == size());
    TMVAssert(minv.rowsize() == size());

    UpperTriMatrixView<T> U = UpperTriMatrixViewOf(minv);
    U = *this;
    U.InvertSelf();
    minv = U * U.Adjoint();
  }

  template <class T> void GenLowerTriMatrix<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == size());
    TMVAssert(minv.rowsize() == size());

    LowerTriMatrixView<T> L = LowerTriMatrixViewOf(minv);
    L = *this;
    L.InvertSelf();
    minv = L * L.Adjoint();
  }

#define InstFile "TMV_TriDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


