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


#ifndef TMV_BaseDiagMatrix_H
#define TMV_BaseDiagMatrix_H

#include "TMV_BaseMatrix.h"

namespace tmv {

  template <class T> class GenDiagMatrix;
  template <class T, IndexStyle I=CStyle> class ConstDiagMatrixView;
  template <class T, IndexStyle I=CStyle> class DiagMatrixView;
  template <class T, IndexStyle I=CStyle> class DiagMatrix;

  template <class T, class Tm> class QuotXD;

  template <class T> struct AssignableToDiagMatrix :
    virtual public AssignableToMatrix<T>
  {
    virtual size_t size() const = 0;
    virtual void AssignToD(const DiagMatrixView<RealType(T)>& m) const = 0;
    virtual void AssignToD(const DiagMatrixView<ComplexType(T)>& m) const = 0;
    virtual inline ~AssignableToDiagMatrix() {}
  };

#ifdef XDEBUG
  template <class T, IndexStyle I> inline std::string Type(
      const DiagMatrix<T,I>& )
  { return std::string("DiagMatrix<")+Type(T())+","+Text(I)+">"; }
  template <class T> inline std::string Type(const GenDiagMatrix<T>& m)
  {
    return std::string("GenDiagMatrix<")+Type(T())+","+Text(m.diag().ct())+">";
  }
  template <class T, IndexStyle I> inline std::string Type(
      const ConstDiagMatrixView<T,I>& m)
  {
    return std::string("ConstDiagMatrixView<")+Type(T())+","+Text(I)+","+
      Text(m.diag().ct())+">";
  }
  template <class T, IndexStyle I> inline std::string Type(
      const DiagMatrixView<T,I>& m)
  {
    return std::string("DiagMatrixView<")+Type(T())+","+Text(I)+","+
      Text(m.diag().ct())+">";
  }
#endif

} // namespace tmv

#endif
