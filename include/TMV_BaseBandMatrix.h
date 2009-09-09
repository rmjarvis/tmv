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


#ifndef TMV_BaseBandMatrix_H
#define TMV_BaseBandMatrix_H

#include "TMV_BaseMatrix.h"

namespace tmv {

  template <class T> class GenBandMatrix;
  template <class T, IndexStyle I=CStyle> class ConstBandMatrixView;
  template <class T, IndexStyle I=CStyle> class BandMatrixView;
  template <class T, StorageType S=RowMajor, IndexStyle I=CStyle> 
    class BandMatrix;

  template <class T> class BandLUDiv;
  template <class T> class BandSVDiv;
  template <class T> class BandQRDiv;
  template <class T, class Tm> class QuotXB;

  template <class T> inline StorageType BaseStorOf(
      const GenBandMatrix<T>& m)
  {
    return (m.stor()==RowMajor || m.stor()==ColMajor || m.stor()==DiagMajor) ?
      m.stor() : DiagMajor;
  }

  size_t BandStorageLength(StorageType s, size_t cs, size_t rs,
      int lo, int hi);

  template <class T1, class T2> void Copy(
      const GenBandMatrix<T1>& m1, const BandMatrixView<T2>& m2);

  template <class T> struct AssignableToBandMatrix :
    virtual public AssignableToMatrix<T>
  {
    virtual int nlo() const = 0;
    virtual int nhi() const = 0;
    virtual void AssignToB(const BandMatrixView<RealType(T)>& m) const = 0;
    virtual void AssignToB(const BandMatrixView<ComplexType(T)>& m) const = 0;
    virtual inline ~AssignableToBandMatrix() {}
  };

  template <class T, StorageType S, IndexStyle I> inline std::string Type(
      const BandMatrix<T,S,I>& )
  { return std::string("BandMatrix<")+Type(T())+","+Text(I)+","+Text(S)+">"; }

  template <class T> inline std::string Type(const GenBandMatrix<T>& m)
  {
    return std::string("GenBandMatrix<")+Type(T())+","+Text(m.ct())+
      ","+Text(m.stor())+">";
  }
  template <class T, IndexStyle I> inline std::string Type(
      const ConstBandMatrixView<T,I>& m)
  {
    return std::string("ConstBandMatrixView<")+Type(T())+","+Text(I)+
      ","+Text(m.ct())+","+Text(m.stor())+">";
  }
  template <class T, IndexStyle I> inline std::string Type(
      const BandMatrixView<T,I>& m)
  {
    return std::string("BandMatrixView<")+Type(T())+","+Text(I)+
      ","+Text(m.ct())+","+Text(m.stor())+">";
  }

} // namespace tmv

#endif
