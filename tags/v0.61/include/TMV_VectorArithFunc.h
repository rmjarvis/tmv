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


#ifndef TMV_VectorArithFunc_H
#define TMV_VectorArithFunc_H

#define CT std::complex<T>

#include "TMV_BaseVector.h"

namespace tmv {

  // v *= x
  template <class T> void MultXV(const T x, const VectorView<T>& v2);
  // v2 = x * v1
  template <class T, class T1> void MultXV(const T x,
      const GenVector<T1>& v1, const VectorView<T>& v2);

  // v2 += x * v1
  template <class T, class T1> void AddVV(const T x, 
      const GenVector<T1>& v1, const VectorView<T>& v2);
  // v3 = x1 * v1 + x2 * v2
  template <class T, class T1, class T2> void AddVV(
      const T x1, const GenVector<T1>& v1,
      const T x2, const GenVector<T2>& v2, const VectorView<T>& v3);

  // v1 * v2 (dot product)
  // Note: the return type is the type of the first vector
  // This is important for mixing complex and real vectors
  template <class T, class T2> T MultVV(
      const GenVector<T>& v1, const GenVector<T2>& v2);

  template <class T, class Tx> void ElementProd(
      const T alpha, const GenVector<Tx>& x, const VectorView<T>& y);
  template <class T, class Tx, class Ty> void AddElementProd(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const VectorView<T>& z);

  template <class T> class VectorComposite :
    public GenVector<T>
  {
    public:

      inline VectorComposite() : itsv(0) {}
      inline VectorComposite(const VectorComposite<T>&) : itsv(0) {}
      virtual inline ~VectorComposite() {}

      const T* cptr() const;
      inline int step() const { return 1; }
      inline ConjType ct() const { return NonConj; }
      inline bool isconj() const { return false; }

    private:
      mutable auto_array<T> itsv;
  };

  // Specialize allowed complex combinations:
  template <class T> inline void MultXV(const T x, const VectorView<CT>& v2)
  { MultXV(CT(x),v2); }
  template <class T, class T1> inline void MultXV(const T x,
      const GenVector<T1>& v1, const VectorView<CT>& v2)
  { MultXV(CT(x),v1,v2); }

  template <class T, class T1> inline void AddVV(const T x, 
      const GenVector<T1>& v1, const VectorView<CT>& v2)
  { AddVV(CT(x),v1,v2); }
  template <class T, class T1, class T2> inline void AddVV(
      const T x1, const GenVector<T1>& v1,
      const T x2, const GenVector<T2>& v2, const VectorView<CT>& v3)
  { AddVV(CT(x1),v1,CT(x2),v2,v3); }
  template <class T> inline void AddVV(
      const CT x1, const GenVector<CT>& v1,
      const CT x2, const GenVector<T>& v2, const VectorView<CT>& v3)
  { AddVV(x2,v2,x1,v1,v3); }

  template <class T> inline CT MultVV(
      const GenVector<T>& v1, const GenVector<CT>& v2)
  { return MultVV(v2,v1); }

  template <class T, class Tx> inline void ElementProd(
      const T alpha, const GenVector<Tx>& x, const VectorView<CT>& y)
  { ElementProd(CT(alpha),x,y); }
  template <class T, class Tx, class Ty> inline void AddElementProd(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const VectorView<CT>& z)
  { AddElementProd(CT(alpha),x,y,z); }
  template <class T> inline void AddElementProd(
      const CT alpha, const GenVector<CT>& x, const GenVector<T>& y,
      const VectorView<CT>& z)
  { AddElementProd(alpha,y,x,z); }

  // Specialize disallowed complex combinations:
  template <class T> inline void MultXV(const CT , const VectorView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta> inline void MultXV(const CT ,
      const GenVector<Ta>& , const VectorView<T>& )
  { TMVAssert(FALSE); }

  template <class T, class T1> inline void AddVV(const CT , 
      const GenVector<T1>& , const VectorView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb> inline void AddVV(
      const CT , const GenVector<Ta>& ,
      const CT , const GenVector<Tb>& , const VectorView<T>& )
  { TMVAssert(FALSE); }

} // namespace tmv

#undef CT

#endif 
