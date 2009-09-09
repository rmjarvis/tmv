///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef TMV_BaseSymMatrix_H
#define TMV_BaseSymMatrix_H

#include "TMV_BaseMatrix.h"

namespace tmv {

  template <class T> class GenSymMatrix;
  template <class T, IndexStyle I=CStyle> class ConstSymMatrixView;
  template <class T, IndexStyle I=CStyle> class SymMatrixView;
  template <class T, UpLoType U=Upper, StorageType S=RowMajor, IndexStyle I=CStyle> 
    class SymMatrix;
  template <class T, UpLoType U=Upper, StorageType S=RowMajor, IndexStyle I=CStyle> 
    class HermMatrix;
  template <class T> class SymDivider;

  template <class T> class SymLUDiv;
  template <class T> class HermCHDiv;
  template <class T> class HermSVDiv;
  template <class T> class SymSVDiv;
  template <class T, class Tm> class QuotXS;

  template <class T> struct AssignableToSymMatrix :
    virtual public AssignableToMatrix<T>
  {
    virtual size_t size() const = 0;
    virtual SymType sym() const = 0;
    inline bool issym() const 
    { return IsReal(T()) || sym() == Sym; }
    inline bool isherm() const 
    { return IsReal(T()) || sym() == Herm; }
    virtual void AssignToS(const SymMatrixView<RealType(T)>& m) const = 0;
    virtual void AssignToS(const SymMatrixView<ComplexType(T)>& m) const = 0;
    virtual inline ~AssignableToSymMatrix() {}
  };

  template <class T, UpLoType U, StorageType S, IndexStyle I>
    inline std::string Type(const SymMatrix<T,U,S,I>& )
    {
      return std::string("SymMatrix<")+Type(T())+","+Text(U)+","+Text(S)+
	","+Text(I)+">";
    }

  template <class T, UpLoType U, StorageType S, IndexStyle I>
    inline std::string Type(const HermMatrix<T,U,S,I>& )
    {
      return std::string("HermMatrix<")+Type(T())+","+Text(U)+","+Text(S)+
	","+Text(I)+">";
    }

  template <class T> inline std::string Type(const GenSymMatrix<T>& m)
  {
    return std::string("GenSymMatrix<")+Type(T())+","+Text(m.sym())+","
      +Text(m.uplo())+","+Text(m.ct())+","+Text(m.stor())+">";
  }

  template <class T, IndexStyle I> inline std::string Type(
      const ConstSymMatrixView<T,I>& m)
  {
    return std::string("ConstSymMatrixView<")+Type(T())+","+Text(m.sym())+
      ","+Text(m.uplo())+ ","+Text(m.stor())+ ","+Text(I)+
      ","+Text(m.ct())+ ">";
  }

  template <class T, IndexStyle I> inline std::string Type(
      const SymMatrixView<T,I>& m)
  {
    return std::string("SymMatrixView<")+Type(T())+","+Text(m.sym())+
      ","+Text(m.uplo())+ ","+Text(m.stor())+ ","+Text(I)+
      ","+Text(m.ct())+ ">";
  }

}

#endif
