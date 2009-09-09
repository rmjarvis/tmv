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


#define RT RealType(T)
#define CT ComplexType(T)

#define DefDivEq(T) \
  inline void LDivEq(const MatrixView<T>& m) const \
{ \
  TMVAssert(colsize() == rowsize()); \
  TMVAssert(m.colsize() == colsize()); \
  DoLDivEq(m);  \
} \
inline void RDivEq(const MatrixView<T>& m) const \
{ \
  TMVAssert(colsize() == rowsize()); \
  TMVAssert(m.rowsize() == rowsize()); \
  DoRDivEq(m); \
} \
inline void Inverse(const MatrixView<T>& minv) const \
{ \
  TMVAssert(minv.rowsize() == colsize()); \
  TMVAssert(minv.colsize() == rowsize()); \
  DoInverse(minv); \
} \


DefDivEq(RT);
DefDivEq(CT);
#undef DefDivEq

#define DefDiv(T1,T2) \
  inline void LDiv(const GenMatrix<T1>& m1, \
      const MatrixView<T2>& m0) const \
{ \
  TMVAssert(m0.rowsize() == m1.rowsize()); \
  TMVAssert(m0.colsize() == rowsize()); \
  TMVAssert(m1.colsize() == colsize()); \
  DoLDiv(m1,m0); \
} \
inline void RDiv(const GenMatrix<T1>& m1, \
    const MatrixView<T2>& m0) const \
{ \
  TMVAssert(m0.colsize() == m1.colsize()); \
  TMVAssert(m1.rowsize() == rowsize()); \
  TMVAssert(m0.rowsize() == colsize()); \
  DoRDiv(m1,m0); \
}
DefDiv(RT,RT);
DefDiv(RT,CT);
DefDiv(CT,CT);
#undef DefDiv

#undef RT
#undef CT

inline void InverseATA(const MatrixView<T>& minv) const
{
  if (colsize() < rowsize()) {
    TMVAssert(minv.rowsize() == colsize());
    TMVAssert(minv.colsize() == colsize());
#ifdef TMVDEBUG
#ifndef NWARN
    std::cout << "Warning InverseATA called for matrix with more columns than rows.\n";
    std::cout << "The result is really (AAt)^-1\n";
#endif
#endif
  } else {
    TMVAssert(minv.rowsize() == rowsize());
    TMVAssert(minv.colsize() == rowsize());
  }
  DoInverseATA(minv);
}

