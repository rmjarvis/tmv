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


#ifndef TMV_BaseSymBandMatrix_H
#define TMV_BaseSymBandMatrix_H

#include "TMV_BaseSymMatrix.h"
#include "TMV_BaseBandMatrix.h"

namespace tmv {

  template <class T> class GenSymBandMatrix;
  template <class T, IndexStyle I=CStyle> class ConstSymBandMatrixView;
  template <class T, IndexStyle I=CStyle> class SymBandMatrixView;
  template <class T, UpLoType U=Upper, StorageType S=RowMajor, IndexStyle I=CStyle>
    class SymBandMatrix;
  template <class T, UpLoType U=Upper, StorageType S=RowMajor, IndexStyle I=CStyle>
    class HermBandMatrix;

  template <class T> class HermBandCHDiv;
  template <class T> class HermBandSVDiv;
  template <class T> class SymBandSVDiv;
  template <class T, class Tm> class QuotXsB;

  template <class T, class T1, class T2> class ProdBB;

  template <class T> inline StorageType BaseStorOf(
      const GenSymBandMatrix<T>& m)
  { return m.stor() != NoMajor ? m.stor() : DiagMajor; }

  template <class T> struct AssignableToSymBandMatrix :
    virtual public AssignableToSymMatrix<T>,
    virtual public AssignableToBandMatrix<T>
  {
    virtual void AssignTosB(const SymBandMatrixView<RealType(T)>& m) const = 0;
    virtual void AssignTosB(const SymBandMatrixView<ComplexType(T)>& m) const = 0;
    virtual inline ~AssignableToSymBandMatrix() {}
  };


  template <class T, UpLoType U, StorageType S, IndexStyle I>
    inline std::string Type(const SymBandMatrix<T,U,S,I>& )
    {
      return std::string("SymBandMatrix<")+Type(T())+","+Text(U)+","+Text(S)+
	","+Text(I)+">";
    }

  template <class T, UpLoType U, StorageType S, IndexStyle I>
    inline std::string Type(const HermBandMatrix<T,U,S,I>& )
    {
      return std::string("HermBandMatrix<")+Type(T())+","+Text(U)+","+Text(S)+
	","+Text(I)+">";
    }

  template <class T> inline std::string Type(const GenSymBandMatrix<T>& m)
  {
    return std::string("GenSymBandMatrix<")+Type(T())+","+Text(m.sym())+","
      +Text(m.uplo())+","+Text(m.ct())+","+Text(m.stor())+">";
  }

  template <class T, IndexStyle I> inline std::string Type(
      const ConstSymBandMatrixView<T,I>& m)
  {
    return std::string("ConstSymBandMatrixView<")+Type(T())+","+Text(m.sym())+
      ","+Text(m.uplo())+ ","+Text(m.stor())+ ","+Text(I)+
      ","+Text(m.ct())+ ">";
  }

  template <class T, IndexStyle I> inline std::string Type(
      const SymBandMatrixView<T,I>& m)
  {
    return std::string("SymBandMatrixView<")+Type(T())+","+Text(m.sym())+
      ","+Text(m.uplo())+ ","+Text(m.stor())+ ","+Text(I)+
      ","+Text(m.ct())+ ">";
  }

}

#endif
