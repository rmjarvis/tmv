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


#ifndef TMV_SmallVectorArithFunc_H
#define TMV_SmallVectorArithFunc_H

#define CT std::complex<T>

namespace tmv {

  // v *= x
  template <class T, class Tv, int N, int S, bool C> inline void MultXV(
      const T x, const SmallVectorView<Tv,N,S,C>& v)
  { for(int i=0;i<N;++i) v[i] *= x; }
  template <class T, int N, int S, bool C> inline void MultXV(
      const CT , const SmallVectorView<T,N,S,C>& )
  { TMVAssert(FALSE); }
  template <int N, class T, class Tv> inline void MultXV(
      const T x, Tv* vp, int S)
  { if (x != T(1)) for(int i=0,ii=0;i<N;++i,ii+=S) vp[ii] *= x; }
  template <int N, class T> inline void MultXV(const CT , T* , int )
  { TMVAssert(FALSE); }

  // v2 = x * v1
  template <class T, class T1, class T2, int N, int S1, bool C1, int S2, bool C2> 
    inline void MultXV(
	const T x, const GenSmallVector<T1,N,S1,C1>& v1, 
	const SmallVectorView<T2,N,S2,C2>& v2)
    { for(int i=0;i<N;++i) v2[i] = x * v1[i]; }
  template <class T, class T1, int N, int S1, bool C1, int S2, bool C2> 
    inline void MultXV(
	const CT , const GenSmallVector<T1,N,S1,C1>& , 
	const SmallVectorView<T,N,S2,C2>& )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, int N, int S1, bool C1> 
    inline void MultXV(
	const T x, const GenSmallVector<T1,N,S1,C1>& v1, 
	T2* v2p, const int S2)
    { for(int i=0,ii=0;i<N;++i,ii+=S2) v2p[ii] = x * v1[i]; }
  template <class T, class T1, int N, int S1, bool C1> 
    inline void MultXV(
	const CT , const GenSmallVector<T1,N,S1,C1>& , 
	T* , const int )
    { TMVAssert(FALSE); }

  // v2 += x * v1
  template <class T1, class T2, int N, int S1, bool C1, int S2, bool C2> 
    inline void AddVV_1(const GenSmallVector<T1,N,S1,C1>& v1, 
	const SmallVectorView<T2,N,S2,C2>& v2)
    { for(int i=0;i<N;++i) v2[i] += v1[i]; }
  template <class T1, class T2, int N, int S1, bool C1, int S2, bool C2> 
    inline void AddVV_m1(const GenSmallVector<T1,N,S1,C1>& v1, 
	const SmallVectorView<T2,N,S2,C2>& v2)
    { for(int i=0;i<N;++i) v2[i] -= v1[i]; }
  template <class T, class T1, class T2, int N, int S1, bool C1, int S2, bool C2> 
    inline void AddVV(const T x, const GenSmallVector<T1,N,S1,C1>& v1, 
	const SmallVectorView<T2,N,S2,C2>& v2)
    { for(int i=0;i<N;++i) v2[i] += x*v1[i]; }
  template <class T, class T1, int N, int S1, bool C1, int S2, bool C2> 
    inline void AddVV(const CT , const GenSmallVector<T1,N,S1,C1>& , 
	const SmallVectorView<T,N,S2,C2>& )
    { TMVAssert(FALSE); }
  template <class T1, class T2, int N, int S1, bool C1> 
    inline void AddVV_1(
	const GenSmallVector<T1,N,S1,C1>& v1, T2* v2p, const int S2)
    { for(int i=0,ii=0;i<N;++i,ii+=S2) v2p[ii] += v1[i]; }
  template <class T1, class T2, int N, int S1, bool C1> 
    inline void AddVV_m1(
	const GenSmallVector<T1,N,S1,C1>& v1, 
	T2* v2p, const int S2)
    { for(int i=0,ii=0;i<N;++i,ii+=S2) v2p[ii] -= v1[i]; }
  template <class T, class T1, class T2, int N, int S1, bool C1> 
    inline void AddVV(
	const T x, const GenSmallVector<T1,N,S1,C1>& v1, 
	T2* v2p, const int S2)
    { for(int i=0,ii=0;i<N;++i,ii+=S2) v2p[ii] += x*v1[i]; }
  template <class T, class T1, int N, int S1, bool C1> 
    inline void AddVV(
	const CT , const GenSmallVector<T1,N,S1,C1>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T1, class T2, int N, int S2, bool C2> 
    inline void AddVV_1(
	const T1* v1p, const int S1, const SmallVectorView<T2,N,S2,C2>& v2)
    { for(int i=0,ii=0;i<N;++i,ii+=S1) v2[i] += v1p[ii]; }
  template <class T1, class T2, int N, int S2, bool C2> 
    inline void AddVV_m1(
	const T1* v1p, const int S1, const SmallVectorView<T2,N,S2,C2>& v2) 
    { for(int i=0,ii=0;i<N;++i,ii+=S1) v2[i] -= v1p[ii]; }
  template <class T, class T1, class T2, int N, int S2, bool C2> 
    inline void AddVV(
	const T x, const T1* v1p, const int S1,
	const SmallVectorView<T2,N,S2,C2>& v2)
    { for(int i=0,ii=0;i<N;++i,ii+=S1) v2[i] += x*v1p[ii]; }
  template <class T, class T1, int N, int S1, bool C1> 
    inline void AddVV(
	const CT , const T* , const int , const SmallVectorView<T1,N,S1,C1>& )
    { TMVAssert(FALSE); }

  // v3 = x1 * v1 + x2 * v2
  template <class T1, class T2, class T3, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_1(
	const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	const SmallVectorView<T3,N,S3,C3>& v3)
    { Copy(v1,v3); AddVV_1(v2,v3); }
  template <class T1, class T2, class T3, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_m1(
	const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	const SmallVectorView<T3,N,S3,C3>& v3)
    { Copy(v1,v3); AddVV_m1(v2,v3); }
  template <class T, class T1, class T2, class T3, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_x(
	const GenSmallVector<T1,N,S1,C1>& v1,
	const T x2, const GenSmallVector<T2,N,S2,C2>& v2, 
	const SmallVectorView<T3,N,S3,C3>& v3)
    { MultXV(x2,v2,v3); AddVV_1(v1,v3); }
  template <class T, class T1, class T2, class T3, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_x_1(
	const T x1, const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	const SmallVectorView<T3,N,S3,C3>& v3)
    { MultXM(x1,v1,v3); AddVV_1(v2,v3); }
  template <class T, class T1, class T2, class T3, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_x_m1(
	const T x1, const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	const SmallVectorView<T3,N,S3,C3>& v3)
    { MultXM(x1,v1,v3); AddVV_m1(v2,v3); }
  template <class T, class T1, class T2, class T3, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV(
	const T x1, const GenSmallVector<T1,N,S1,C1>& v1,
	const T x2, const GenSmallVector<T2,N,S2,C2>& v2, 
	const SmallVectorView<T3,N,S3,C3>& v3)
    { MultXM(x1,v1,v3); AddVV(x2,v2,v3); }
  template <class T, class T2, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_1(
	const GenSmallVector<CT,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_1(
	const GenSmallVector<T,N,S1,C1>& ,
	const GenSmallVector<CT,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, class T2, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_m1(
	const GenSmallVector<CT,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_m1(
	const GenSmallVector<T,N,S1,C1>& ,
	const GenSmallVector<CT,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_1_x(
	const GenSmallVector<T1,N,S1,C1>& ,
	const CT , const GenSmallVector<T2,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_x_1(
	const CT , const GenSmallVector<T1,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV_x_m1(
	const CT , const GenSmallVector<T1,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, int N, int S1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline void AddVV(
	const CT , const GenSmallVector<T1,N,S1,C1>& ,
	const CT , const GenSmallVector<T2,N,S2,C2>& , 
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }
  template <class T1, class T2, class T3, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_1(
	const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	T3* v3p, const int S3)
    { Copy(v1,v3p,S3); AddVV_1(v2,v3p,S3); }
  template <class T1, class T2, class T3, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_m1(
	const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	T3* v3p, const int S3)
    { Copy(v1,v3p,S3); AddVV_m1(v2,v3p,S3); }
  template <class T, class T1, class T2, class T3, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_x(
	const GenSmallVector<T1,N,S1,C1>& v1,
	const T x2, const GenSmallVector<T2,N,S2,C2>& v2, 
	T3* v3p, const int S3)
    { MultXV(x2,v2,v3p,S3); AddVV_1(v1,v3p,S3); }
  template <class T, class T1, class T2, class T3, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_x_1(
	const T x1, const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	T3* v3p, const int S3)
    { MultXV(x1,v1,v3p,S3); AddVV_1(v2,v3p,S3); }
  template <class T, class T1, class T2, class T3, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_x_m1(
	const T x1, const GenSmallVector<T1,N,S1,C1>& v1,
	const GenSmallVector<T2,N,S2,C2>& v2, 
	T3* v3p, const int S3)
    { MultXV(x1,v1,v3p,S3); AddVV_m1(v2,v3p,S3); }
  template <class T, class T1, class T2, class T3, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV(
	const T x1, const GenSmallVector<T1,N,S1,C1>& v1,
	const T x2, const GenSmallVector<T2,N,S2,C2>& v2, 
	T3* v3p, const int S3)
    { MultXV(x1,v1,v3p,S3); AddVV(x2,v2,v3p,S3); }
  template <class T, class T2, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_1(
	const GenSmallVector<CT,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_1(
	const GenSmallVector<T,N,S1,C1>& ,
	const GenSmallVector<CT,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, class T2, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_m1(
	const GenSmallVector<CT,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_m1(
	const GenSmallVector<T,N,S1,C1>& ,
	const GenSmallVector<CT,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_1_x(
	const GenSmallVector<T1,N,S1,C1>& ,
	const CT , const GenSmallVector<T2,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_x_1(
	const CT , const GenSmallVector<T1,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV_x_m1(
	const CT , const GenSmallVector<T1,N,S1,C1>& ,
	const GenSmallVector<T2,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, int N, int S1, int S2, bool C1, bool C2> 
    inline void AddVV(
	const CT , const GenSmallVector<T1,N,S1,C1>& ,
	const CT , const GenSmallVector<T2,N,S2,C2>& , T* , const int )
    { TMVAssert(FALSE); }

  // v1 * v2 (dot product)
  template <class T, class T2, int N, int S1, int S2, bool C1, bool C2>
    inline T MultVV(
	const GenSmallVector<T,N,S1,C1>& v1, 
	const GenSmallVector<T2,N,S2,C2>& v2)
    {
      T res(0);
      for(int i=0;i<N;++i) res += v1[i]*v2[i];
      return res;
    }
  template <class T, int N, int S1, int S2, bool C1, bool C2>
    inline CT MultVV(
	const GenSmallVector<T,N,S1,C1>& v1, 
	const GenSmallVector<CT,N,S2,C2>& v2)
    { return MultVV(v2,v1); }

  template <class T, class T1, class T2, int N, int S1, int S2, bool C1, bool C2> 
    inline void ElementProd(
	const T alpha, const GenSmallVector<T1,N,S1,C1>& v1, 
	const SmallVectorView<T2,N,S2,C2>& v2)
    {
      if (alpha == T(1))
	for(int i=0;i<N;++i) v2[i] = v1[i] * v2[i];
      else if (alpha == T(-1))
	for(int i=0;i<N;++i) v2[i] = -v1[i] * v2[i];
      else
	for(int i=0;i<N;++i) v2[i] = alpha * v1[i] * v2[i];
    }
  template <class T, class T1, int N, int S1, int S2, bool C1, bool C2> 
    inline void ElementProd(
	const CT , const GenSmallVector<T1,N,S1,C1>& , 
	const SmallVectorView<T,N,S2,C2>& )
    { TMVAssert(FALSE); }
  template <class T, class T1, class T2, class T3, int N, int S1, int S2, int S3, bool C1, bool C2, bool C3> 
    inline void AddElementProd(
	const T alpha, const GenSmallVector<T1,N,S1,C1>& v1, 
	const GenSmallVector<T2,N,S2,C2>& v2,
	const SmallVectorView<T3,N,S3,C3>& v3)
    {
      if (alpha == T(1))
	for(int i=0;i<N;++i) v3[i] += v1[i] * v2[i];
      else if (alpha == T(-1))
	for(int i=0;i<N;++i) v3[i] -= v1[i] * v2[i];
      else
	for(int i=0;i<N;++i) v3[i] += alpha * v1[i] * v2[i];
    }
  template <class T, class T1, class T2, int N, int S1, int S2, int S3, bool C1, bool C2, bool C3> 
    inline void AddElementProd(
	const CT , const GenSmallVector<T1,N,S1,C1>& , 
	const GenSmallVector<T2,N,S2,C2>& ,
	const SmallVectorView<T,N,S3,C3>& )
    { TMVAssert(FALSE); }

} // namespace tmv

#undef CT

#endif 
