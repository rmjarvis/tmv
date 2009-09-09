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


#ifndef TMV_BaseTriMatrix_H
#define TMV_BaseTriMatrix_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

  template <class T> 
  class GenUpperTriMatrix;

  template <class T> 
  class GenLowerTriMatrix;

  template <class T, IndexStyle I=CStyle> 
  class ConstUpperTriMatrixView;

  template <class T, IndexStyle I=CStyle> 
  class ConstLowerTriMatrixView;

  template <class T, IndexStyle I=CStyle> 
  class UpperTriMatrixView;

  template <class T, IndexStyle I=CStyle> 
  class LowerTriMatrixView;

  template <class T, DiagType D=NonUnitDiag, StorageType S=ColMajor, IndexStyle I=CStyle>
  class UpperTriMatrix;

  template <class T, DiagType D=NonUnitDiag, StorageType S=ColMajor, IndexStyle I=CStyle>
  class LowerTriMatrix;

  template <class T> 
  class UpperTriDiv;

  template <class T> 
  class LowerTriDiv;

  template <class T, class Tm> 
  class QuotXU;

  template <class T, class Tm> 
  class QuotXL;

  template <class T1, class T2> 
  void Copy(
      const GenUpperTriMatrix<T1>& m1, const UpperTriMatrixView<T2>& m2);

  template <class T> 
  struct AssignableToUpperTriMatrix :
    virtual public AssignableToMatrix<T>
  {
    virtual size_t size() const = 0;
    virtual DiagType dt() const = 0;
    virtual void AssignToU(const UpperTriMatrixView<RealType(T)>& m) const = 0;
    virtual void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m) const = 0;
    virtual inline ~AssignableToUpperTriMatrix() {}
  };

  template <class T> 
  struct AssignableToLowerTriMatrix :
    virtual public AssignableToMatrix<T>
  {
    virtual size_t size() const = 0;
    virtual DiagType dt() const = 0;
    virtual void AssignToL(const LowerTriMatrixView<RealType(T)>& m) const = 0;
    virtual void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m) const = 0;
    virtual inline ~AssignableToLowerTriMatrix() {}
  };

  template <class T, DiagType D, StorageType S, IndexStyle I>
  inline std::string TypeText(const UpperTriMatrix<T,D,S,I>& )
  {
    return std::string("UpperTriMatrix<")+TypeText(T())+","
    + Text(D)+","+Text(S)+","+Text(I)+">";
  }
  template <class T, DiagType D, StorageType S, IndexStyle I>
  inline std::string TypeText(const LowerTriMatrix<T,D,S,I>& )
  {
    return std::string("LowerTriMatrix<")+TypeText(T())+","
    + Text(D)+","+Text(S)+","+Text(I)+">";
  }

  template <class T> 
  inline std::string TypeText(const GenUpperTriMatrix<T>& m)
  {
    return std::string("GenUpperTriMatrix<")+TypeText(T())+","
    + Text(m.dt())+","+Text(m.ct())+","+Text(m.stor())+">";
  }
  template <class T> 
  inline std::string TypeText(const GenLowerTriMatrix<T>& m)
  {
    return std::string("GenLowerTriMatrix<")+TypeText(T())+","
    + Text(m.dt())+","+Text(m.ct())+","+Text(m.stor())+">";
  }

  template <class T, IndexStyle I> 
  inline std::string TypeText(
      const ConstUpperTriMatrixView<T,I>& m)
  {
    return std::string("ConstUpperTriMatrixView<")+TypeText(T())+","
    + Text(m.dt())+","+Text(m.ct())+","+Text(m.stor())+","+Text(I)+">";
  }
  template <class T, IndexStyle I> 
  inline std::string TypeText(
      const ConstLowerTriMatrixView<T,I>& m)
  {
    return std::string("ConstLowerTriMatrixView<")+TypeText(T())+","
    + Text(m.dt())+","+Text(m.ct())+","+Text(m.stor())+","+Text(I)+">";
  }

  template <class T, IndexStyle I> 
  inline std::string TypeText(
      const UpperTriMatrixView<T,I>& m)
  {
    return std::string("UpperTriMatrixView<")+TypeText(T())+","
    + Text(m.dt())+","+Text(m.ct())+","+Text(m.stor())+","+Text(I)+">";
  }
  template <class T, IndexStyle I> 
  inline std::string TypeText(
      const LowerTriMatrixView<T,I>& m)
  {
    return std::string("LowerTriMatrixView<")+TypeText(T())+","
    + Text(m.dt())+","+Text(m.ct())+","+Text(m.stor())+","+Text(I)+">";
  }

} // namespace tmv

#endif
