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


#ifndef TMV_SmallVectorArithFunc_H
#define TMV_SmallVectorArithFunc_H

#define CT std::complex<T>

namespace tmv {

  // v *= x
  template <class T, class Tv, size_t N, int S, bool C> inline void MultXV(
      const T x, const SmallVectorView<Tv,N,S,C>& v)
  {
    typename SmallVectorView<Tv,N,S,C>::iterator vi = v.begin();
    for(size_t i=N;i>0;--i,++vi) *vi *= x;
  }
  template <class T, size_t N, int S, bool C> inline void MultXV(
      const CT , const SmallVectorView<T,N,S,C>& )
  { TMVAssert(FALSE); }
  template <size_t N, class T, class Tv> inline void MultXV(
      const T x, Tv* vp, int S)
  { if (x != T(1)) for(size_t i=N;i>0;--i,vp+=S) *vp *= x; }
  template <size_t N, class T> inline void MultXV(const CT , T* , int )
  { TMVAssert(FALSE); }

  // v2 = x * v1
  template <class T, class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2> 
    inline void MultXV(
	const T x, const GenSmallVector<T1,N,S1,C1>& v1, 
	const SmallVectorView<T2,N,S2,C2>& v2)
    {
      typename GenSmallVector<T1,N,S1,C1>::const_iterator v1i = v1.begin();
      typename SmallVectorView<T2,N,S2,C2>::iterator v2i = v2.begin();
      for(size_t i=N;i>0;--i,++v1i,++v2i) *v2i = x * (*v1i);
    }
  template <class T, class T1, size_t N, int S1, bool C1, int S2, bool C2> 
    inline void MultXV(
	const CT , const GenSmallVector<T1,N,S1,C1>& , 
	const SmallVectorView<T,N,S2,C2>& )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, size_t N, int S1, bool C1> 
    inline void MultXV(
	const T x, const GenSmallVector<T1,N,S1,C1>& v1, 
	T2* v2p, const int S2)
    {
      typename GenSmallVector<T1,N,S1,C1>::const_iterator v1i = v1.begin();
      for(size_t i=N;i>0;--i,++v1i,v2p+=S2) *v2p = x * (*v1i);
    }
  template <class T, class T1, size_t N, int S1, bool C1> 
    inline void MultXV(
	const CT , const GenSmallVector<T1,N,S1,C1>& , 
	T* , const int )
    { TMVAssert(FALSE); }

  // v2 += x * v1
  template <class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2> 
    inline void AddVV_1(const GenSmallVector<T1,N,S1,C1>& v1, 
	const SmallVectorView<T2,N,S2,C2>& v2)
    {
      typename GenSmallVector<T1,N,S1,C1>::const_iterator v1i = v1.begin();
      typename SmallVectorView<T2,N,S2,C2>::iterator v2i = v2.begin();
      for(size_t i=N;i>0;--i,++v1i,++v2i) *v2i += *v1i;
    }
  template <class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2> 
    inline void AddVV_m1(const GenSmallVector<T1,N,S1,C1>& v1, 
	const SmallVectorView<T2,N,S2,C2>& v2)
    {
      typename GenSmallVector<T1,N,S1,C1>::const_iterator v1i = v1.begin();
      typename SmallVectorView<T2,N,S2,C2>::iterator v2i = v2.begin();
      for(size_t i=N;i>0;--i,++v1i,++v2i) *v2i -= *v1i;
    }
  template <class T, class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2> 
    inline void AddVV(const T x, const GenSmallVector<T1,N,S1,C1>& v1, 
	const SmallVectorView<T2,N,S2,C2>& v2)
    {
      typename GenSmallVector<T1,N,S1,C1>::const_iterator v1i = v1.begin();
      typename SmallVectorView<T2,N,S2,C2>::iterator v2i = v2.begin();
      for(size_t i=N;i>0;--i,++v1i,++v2i) *v2i += x * (*v1i);
    }
  template <class T, class T1, size_t N, int S1, bool C1, int S2, bool C2> 
    inline void AddVV(const CT , const GenSmallVector<T1,N,S1,C1>& , 
	const SmallVectorView<T,N,S2,C2>& )
    { TMVAssert(FALSE); }
  template <class T1, class T2, size_t N, int S1, bool C1> 
    inline void AddVV_1(
	const GenSmallVector<T1,N,S1,C1>& v1, T2* v2p, const int S2)
    {
      typename GenSmallVector<T1,N,S1,C1>::const_iterator v1i = v1.begin();
      for(size_t i=N;i>0;--i,++v1i,v2p+=S2) *v2p += *v1i;
    }
  template <class T1, class T2, size_t N, int S1, bool C1> 
    inline void AddVV_m1(
	const GenSmallVector<T1,N,S1,C1>& v1, 
	T2* v2p, const int S2)
    {
      typename GenSmallVector<T1,N,S1,C1>::const_iterator v1i = v1.begin();
      for(size_t i=N;i>0;--i,++v1i,v2p+=S2) *v2p -= *v1i;
    }
  template <class T, class T1, class T2, size_t N, int S1, bool C1> 
    inline void AddVV(
	const T x, const GenSmallVector<T1,N,S1,C1>& v1, 
	T2* v2p, const int S2)
    {
      typename GenSmallVector<T1,N,S1,C1>::const_iterator v1i = v1.begin();
      for(size_t i=N;i>0;--i,++v1i,v2p+=S2) *v2p += x * (*v1i);
    }
  template <class T, class T1, size_t N, int S1, bool C1> 
    inline void AddVV(
	const CT , const GenSmallVector<T1,N,S1,C1>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T1, class T2, size_t N, int S2, bool C2> 
    inline void AddVV_1(
	const T1* v1p, const int S1, const SmallVectorView<T2,N,S2,C2>& v2)
    {
      typename SmallVectorView<T2,N,S2,C2>::iterator v2i = v2.begin();
      for(size_t i=N;i>0;--i,++v2i,v1p+=S1) *v2i += *v1p;
    }
  template <class T1, class T2, size_t N, int S2, bool C2> 
    inline void AddVV_m1(
	const T1* v1p, const int S1, const SmallVectorView<T2,N,S2,C2>& v2) 
    {
      typename SmallVectorView<T2,N,S2,C2>::iterator v2i = v2.begin();
      for(size_t i=N;i>0;--i,++v2i,v1p+=S1) *v2i -= *v1p;
    }
  template <class T, class T1, class T2, size_t N, int S2, bool C2> 
    inline void AddVV(
	const T x, const T1* v1p, const int S1,
	const SmallVectorView<T2,N,S2,C2>& v2)
    {
      typename SmallVectorView<T2,N,S2,C2>::iterator v2i = v2.begin();
      for(size_t i=N;i>0;--i,++v2i,v1p+=S1) *v2i += x * (*v1p);
    }
  template <class T, class T1, size_t N, int S1, bool C1> 
    inline void AddVV(
	const CT , const T* , const int , const SmallVectorView<T1,N,S1,C1>& )
    { TMVAssert(FALSE); }

  // v3 = x1 * v1 + x2 * v2
  template <class T1, class T2, class T3, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_1(
	const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	const SmallVectorView<T3,N,S3,C3>& v3)
    { Copy(v1,v3); AddVV_1(v2,v3); }
  template <class T1, class T2, class T3, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_m1(
	const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	const SmallVectorView<T3,N,S3,C3>& v3)
    { Copy(v1,v3); AddVV_m1(v2,v3); }
  template <class T, class T1, class T2, class T3, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_x(
	const GenSmallVector<T1,N,S1,C1>& v1,
	const T x2, const GenSmallVector<T2,N,S2,C2>& v2, 
	const SmallVectorView<T3,N,S3,C3>& v3)
    { MultXV(x2,v2,v3); AddVV_1(v1,v3); }
  template <class T, class T1, class T2, class T3, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_x_1(
	const T x1, const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	const SmallVectorView<T3,N,S3,C3>& v3)
    { MultXM(x1,v1,v3); AddVV_1(v2,v3); }
  template <class T, class T1, class T2, class T3, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_x_m1(
	const T x1, const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	const SmallVectorView<T3,N,S3,C3>& v3)
    { MultXM(x1,v1,v3); AddVV_m1(v2,v3); }
  template <class T, class T1, class T2, class T3, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV(
	const T x1, const GenSmallVector<T1,N,S1,C1>& v1,
	const T x2, const GenSmallVector<T2,N,S2,C2>& v2, 
	const SmallVectorView<T3,N,S3,C3>& v3)
    { MultXM(x1,v1,v3); AddVV(x2,v2,v3); }
  template <class T, class T2, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_1(
	const GenSmallVector<CT,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_1(
	const GenSmallVector<T,N,S1,C1>& ,
	const GenSmallVector<CT,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, class T2, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_m1(
	const GenSmallVector<CT,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_m1(
	const GenSmallVector<T,N,S1,C1>& ,
	const GenSmallVector<CT,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_x(
	const GenSmallVector<T1,N,S1,C1>& ,
	const CT , const GenSmallVector<T2,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_x_1(
	const CT , const GenSmallVector<T1,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_x_m1(
	const CT , const GenSmallVector<T1,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV(
	const CT , const GenSmallVector<T1,N,S1,C1>& ,
	const CT , const GenSmallVector<T2,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T1, class T2, class T3, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_1(
	const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	T3* v3p, const int S3)
    { Copy(v1,v3p,S3); AddVV_1(v2,v3p,S3); }
  template <class T1, class T2, class T3, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_m1(
	const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	T3* v3p, const int S3)
    { Copy(v1,v3p,S3); AddVV_m1(v2,v3p,S3); }
  template <class T, class T1, class T2, class T3, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_x(
	const GenSmallVector<T1,N,S1,C1>& v1,
	const T x2, const GenSmallVector<T2,N,S2,C2>& v2, 
	T3* v3p, const int S3)
    { MultXV(x2,v2,v3p,S3); AddVV_1(v1,v3p,S3); }
  template <class T, class T1, class T2, class T3, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_x_1(
	const T x1, const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	T3* v3p, const int S3)
    { MultXV(x1,v1,v3p,S3); AddVV_1(v2,v3p,S3); }
  template <class T, class T1, class T2, class T3, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_x_m1(
	const T x1, const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	T3* v3p, const int S3)
    { MultXV(x1,v1,v3p,S3); AddVV_m1(v2,v3p,S3); }
  template <class T, class T1, class T2, class T3, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV(
	const T x1, const GenSmallVector<T1,N,S1,C1>& v1,
	const T x2, const GenSmallVector<T2,N,S2,C2>& v2, 
	T3* v3p, const int S3)
    { MultXV(x1,v1,v3p,S3); AddVV(x2,v2,v3p,S3); }
  template <class T, class T2, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_1(
	const GenSmallVector<CT,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_1(
	const GenSmallVector<T,N,S1,C1>& ,
	const GenSmallVector<CT,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, class T2, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_m1(
	const GenSmallVector<CT,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_m1(
	const GenSmallVector<T,N,S1,C1>& ,
	const GenSmallVector<CT,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_x(
	const GenSmallVector<T1,N,S1,C1>& ,
	const CT , const GenSmallVector<T2,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_x_1(
	const CT , const GenSmallVector<T1,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_x_m1(
	const CT , const GenSmallVector<T1,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV(
	const CT , const GenSmallVector<T1,N,S1,C1>& ,
	const CT , const GenSmallVector<T2,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }

  // v1 * v2 (dot product)
  template <class T, class T2, size_t N, int S1, int S2, bool C1, bool C2>
    inline T MultVV(
	const GenSmallVector<T,N,S1,C1>& v1, 
	const GenSmallVector<T2,N,S2,C2>& v2)
    {
      typename GenSmallVector<T,N,S1,C1>::const_iterator v1i = v1.begin();
      typename GenSmallVector<T2,N,S2,C2>::const_iterator v2i = v2.begin();
      T res(0);
      for(size_t i=N;i>0;--i,++v1i,++v2i) res += (*v1i)*(*v2i);
      return res;
    }
  template <class T, size_t N, int S1, int S2, bool C1, bool C2>
    inline CT MultVV(
	const GenSmallVector<T,N,S1,C1>& v1, 
	const GenSmallVector<CT,N,S2,C2>& v2)
    { return MultVV(v2,v1); }

  template <class T, class T1, class T2, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void ElementProd(
	const T alpha, const GenSmallVector<T1,N,S1,C1>& v1, 
	const SmallVectorView<T2,N,S2,C2>& v2)
    {
      typename GenSmallVector<T1,N,S1,C1>::const_iterator v1i = v1.begin();
      typename SmallVectorView<T2,N,S2,C2>::iterator v2i = v2.begin();
      if (alpha == T(1))
	for(size_t i=N;i>0;--i,++v1i,++v2i) 
	  *v2i = (*v1i) * (*v2i);
      else if (alpha == T(-1))
	for(size_t i=N;i>0;--i,++v1i,++v2i) 
	  *v2i = (-*v1i) * (*v2i);
      else
	for(size_t i=N;i>0;--i,++v1i,++v2i) 
	  *v2i = alpha * (*v1i) * (*v2i);
    }
  template <class T, class T1, size_t N, int S1, int S2, bool C1, bool C2> 
    inline void ElementProd(
	const CT , const GenSmallVector<T1,N,S1,C1>& , 
	const SmallVectorView<T,N,S2,C2>& )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, class T3, size_t N, int S1, int S2, int S3, bool C1, bool C2, bool C3> 
    inline void AddElementProd(
	const T alpha, const GenSmallVector<T1,N,S1,C1>& v1, 
	const GenSmallVector<T2,N,S2,C2>& v2,
	const SmallVectorView<T3,N,S3,C3>& v3)
    {
      typename GenSmallVector<T1,N,S1,C1>::const_iterator v1i = v1.begin();
      typename GenSmallVector<T2,N,S2,C2>::const_iterator v2i = v2.begin();
      typename SmallVectorView<T3,N,S3,C3>::iterator v3i = v3.begin();
      if (alpha == T(1))
	for(size_t i=N;i>0;--i,++v1i,++v2i,++v3i) 
	  *v3i += (*v1i) * (*v2i);
      else if (alpha == T(-1))
	for(size_t i=N;i>0;--i,++v1i,++v2i,++v3i) 
	  *v3i -= (*v1i) * (*v2i);
      else
	for(size_t i=N;i>0;--i,++v1i,++v2i,++v3i) 
	  *v3i += alpha * (*v1i) * (*v2i);
    }
  template <class T, class T1, class T2, size_t N, int S1, int S2, int S3, bool C1, bool C2, bool C3> 
    inline void AddElementProd(
	const CT , const GenSmallVector<T1,N,S1,C1>& , 
	const GenSmallVector<T2,N,S2,C2>& ,
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }

} // namespace tmv

#undef CT

#endif 
