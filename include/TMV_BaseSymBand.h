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


#ifndef TMV_BaseBand_H
#define TMV_BaseBand_H

#include "TMV_BaseMatrix.h"

namespace tmv {

  template <class T> class GenSymBandMatrix;
  template <class T, IndexStyle I=CStyle> class ConstSymBandMatrixView;
  template <class T, IndexStyle I=CStyle> class SymBandMatrixView;
  template <class T, UpLoType U=Upper, StorageType S=RowMajor, IndexStyle I=CStyle>
    class SymBandMatrix;
  template <class T, UpLoType U=Upper, StorageType S=RowMajor, IndexStyle I=CStyle>
    class HermBandMatrix;
  template <class T> class SymBandMatrixComposite;
  template <class T, class Tm> class SumsBX;
  template <class T, class Tm> class ProdXsB;
  template <class T, class T1, class T2> class SumsBsB;
  template <class T> class SymDivider;

  template <class T> class SymLUDiv;
  template <class T> class HermCHDiv;
  template <class T> class HermSVDiv;
  template <class T> class SymSVDiv;
  template <class T, class Tm> class QuotXsB;


} // namespace tmv

#endif
