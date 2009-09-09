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



#include "TMV_Blas.h"
#include "TMV_Givens.h"
#include "TMV_Vector.h"
#include "TMV_VIt.h"

namespace tmv {

#define RT RealType(T)

  // 
  // Givens_Rotate
  //
  template <class T> Givens<T> Givens_Rotate(T& x, T& y)
  {
    // Use a Givens matrix G to rotate the vector so y = 0:
    // G [ x ] = [ r ]
    //   [ y ]   [ 0 ]
    // Also, return the Givens matrix used to do so.
    
    const RT eps = Epsilon<T>();

    RT maxabs_x = MAXABS(x);
    RT maxabs_y = MAXABS(y);
    if (maxabs_y == RT(0) || maxabs_y*maxabs_y*eps == RT(0)) {
      y = RT(0);
      return Givens<T>(RT(1),T(0));
    } else if (maxabs_x == RT(0) || maxabs_x*maxabs_x*eps == RT(0)) {
      x = RT(0);
      RT absy = ABS(y);
      T s = SIGN(CONJ(y),absy);
      x = absy; y = RT(0);
      return Givens<T>(RT(0),s);
    } else {
      // This test is better than abs(x) > abs(y) for complex,
      // since it saves 2 sqrt's and has roughly the same effect
      // of preventing overflow and reducing rounding errors.
      if (maxabs_x > maxabs_y) {
	if (maxabs_y <= SqrtEpsilon<T>()*maxabs_x) {
	  // Then everything simplifies:
	  // c = 1
	  // s = (y/x)* = xy*/|x|^2
	  // r = f
	  RT invnormx = RT(1)/NORM(x);
	  T s = x*CONJ(y)*invnormx;
	  y = RT(0);
	  return Givens<T>(RT(1),s);
	} else {
	  // c = 1/sqrt(1+|y/x|^2)
	  // s = (y/x)*/sqrt(1+|y/x|^2)
	  // r = x sqrt(1+|y/x|^2)
	  // We calculate c-1 rather than c, since it will be first order in
	  // y/x, which is nominally small.  This leads to better accuracy
	  // on the Givens multiplies.
	  // c-1 = (1-sqrt(1+n))/sqrt(1+n)
	  //     = -n/sqrt(1+n)/(1+sqrt(1+n)) 
	  //     = -n/(1+n+sqrt(1+n)
	  // Likewise, we get a slightly more accurate calculation of r
	  // if we calculate r-x and add this to x:
	  // r-x = x (sqrt(1+|y/x|^2)-1) = x |y/x|^2 / (1 + sqrt(1+|y/x|^2))
	  T yoverx = y/x;
	  RT n = NORM(yoverx);
	  RT sqrtfactor = SQRT(RT(1)+n);
	  RT invsqrtfactor = RT(1)/sqrtfactor;
	  T s = CONJ(yoverx)*invsqrtfactor;
	  x += x*(n/(RT(1)+sqrtfactor));
	  y = RT(0);
	  RT c = RT(1)/sqrtfactor;
	  return Givens<T>(c,s);
	}
      } else {
	// As above, we store c-1 rather than c
	// even though it isn't small here.
	T xovery = x/y;
	RT n = NORM(xovery);
	RT absxovery = SQRT(n);
	if (n <= Epsilon<T>()) {
	  // c = |x/y|
	  // s = (x/y)/|x/y|
	  // r = x/|x/y|
	  T s = SIGN(xovery,absxovery);
	  x = s*y;
	  y = RT(0);
	  return Givens<T>(absxovery,s);
	} else {
	  // c = |x/y|/sqrt(1+|x/y|^2)
	  // s = (x/y)/|x/y|/sqrt(1+|x/y|^2)
	  // r = x/|x/y| sqrt(|x|^2+|y|^2)
	  RT sqrtfactor = SQRT(RT(1)+n);
	  RT invsqrtfactor = RT(1)/sqrtfactor;
	  T signxovery = SIGN(xovery,absxovery); // (x/y)/|x/y|
	  T s = signxovery*invsqrtfactor;
	  x = y * signxovery * sqrtfactor; 
	  y = T(0);
	  RT c = absxovery*invsqrtfactor;
	  return Givens<T>(c,s);
	}
      }
    }
  }

  //
  // GivensMult Scalars
  //
  template <class T, class Tx> void GivensMult(
      RT c, T s, Tx& x, Tx& y) 
  {
    // [ x' ] = [  c  s  ] [ x ] = [  cx+sy  ]
    // [ y' ]   [ -s* c* ] [ y ]   [ c*y-s*x ]
    Tx xx = c*x+s*y;
    y = c*y-tmv::CONJ(s)*x;
    x = xx;
  }
  template <class T, class Tx> void GivensMult(
      RT c, T s, ConjRef<Tx> x, ConjRef<Tx> y)
  {
    Tx xx = c*x+s*y;
    y = c*y-tmv::CONJ(s)*x;
    x = xx;
  }
  template <class T, class Tx> void GivensMult(
      RT c, T s, VarConjRef<Tx> x, VarConjRef<Tx> y)
  {
    Tx xx = c*x+s*y;
    y = c*y-tmv::CONJ(s)*x;
    x = xx;
  }

  //
  // GivensMult vectors
  //

  template <bool c0, class T, class Tx> static void NonBlasGivensMult(
      RT c, T s,
      const VectorView<Tx>& v1, const VectorView<Tx>& v2)
  {
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v1.size() > 0);
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);
    TMVAssert(!v1.SameAs(v2));
    TMVAssert(v2.ct() == NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(c0 == v1.isconj());

    Tx* v1ptr = v1.ptr();
    Tx* v2ptr = v2.ptr();
    const int step0 = v1.step();
    const int step1 = v2.step();

    if (step0 == 1 && step1 == 1)
      for(int i=v1.size();i>0;--i,++v1ptr,++v2ptr) {
#ifdef TMVFLDEBUG
	TMVAssert(v1ptr >= v1.first);
	TMVAssert(v1ptr < v1.last);
	TMVAssert(v2ptr >= v2.first);
	TMVAssert(v2ptr < v2.last);
#endif
	GivensMult(c,s,*v1ptr,*v2ptr); 
      }
    else
      for(int i=v1.size();i>0;--i,v1ptr+=step0,v2ptr+=step1) {
#ifdef TMVFLDEBUG
	TMVAssert(v1ptr >= v1.first);
	TMVAssert(v1ptr < v1.last);
	TMVAssert(v2ptr >= v2.first);
	TMVAssert(v2ptr < v2.last);
#endif
	GivensMult(c,s,*v1ptr,*v2ptr); 
      }
  }

  template <class T, class Tx> static void DoGivensMult(
      RT c, T s,
      const VectorView<Tx>& v1, const VectorView<Tx>& v2)
  {
    if (v1.isconj()) NonBlasGivensMult<true>(c,s,v1,v2); 
    else NonBlasGivensMult<false>(c,s,v1,v2); 
  }

#ifdef BLAS
#ifdef INST_DOUBLE
  template <> void DoGivensMult(
      double c, double s,
      const VectorView<double>& v1, const VectorView<double>& v2)
  { 
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n=v1.size();
    int s1=v1.step();
    int s2=v2.step();
    double* v1p = v1.ptr();
    if (s1<0) v1p += (n-1)*s1;
    double* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*s2;
    BLASNAME(drot) (BLASV(n),BLASP(v1p),BLASV(s1),
	BLASP(v2p),BLASV(s2),BLASV(c),BLASV(s)); 
  }
  template <> void DoGivensMult(
      double c, double s,
      const VectorView<std::complex<double> >& v1,
      const VectorView<std::complex<double> >& v2)
  { 
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n=v1.size();
    int s1=v1.step();
    int s2=v2.step();
    std::complex<double>* v1p = v1.ptr();
    if (s1<0) v1p += (n-1)*s1;
    std::complex<double>* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*s2;
    BLASNAME(zdrot) (BLASV(n),BLASP(v1p),BLASV(s1),
	BLASP(v2p),BLASV(s2),BLASV(c),BLASV(s)); 
  }
#ifdef ELAP
  template <> void DoGivensMult(
      double c, std::complex<double> s,
      const VectorView<std::complex<double> >& v1,
      const VectorView<std::complex<double> >& v2)
  {
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n = v1.size();
    int s1 = v1.step();
    int s2 = v2.step();
    std::complex<double>* v1p = v1.ptr();
    if (s1<0) v1p += (n-1)*s1;
    std::complex<double>* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*s2;
    LAPNAMEX(zrot) (LAPV(n),LAPP(v1p),LAPV(s1),LAPP(v2p),LAPV(s2),
	LAPV(c),LAPP(&s)); 
  }
#endif // ELAP
#endif // DOUBLE
#ifdef INST_FLOAT
  template <> void DoGivensMult(
      float c, float s,
      const VectorView<float>& v1, const VectorView<float>& v2)
  { 
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n=v1.size();
    int s1=v1.step();
    int s2=v2.step();
    float* v1p = v1.ptr();
    if (s1<0) v1p += (n-1)*s1;
    float* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*s2;
    BLASNAME(srot) (BLASV(n),BLASP(v1p),BLASV(s1),
	BLASP(v2p),BLASV(s2),BLASV(c),BLASV(s)); 
  }
  template <> void DoGivensMult(
      float c, float s,
      const VectorView<std::complex<float> >& v1,
      const VectorView<std::complex<float> >& v2)
  { 
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n=v1.size();
    int s1=v1.step();
    int s2=v2.step();
    std::complex<float>* v1p = v1.ptr();
    if (s1<0) v1p += (n-1)*s1;
    std::complex<float>* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*s2;
    BLASNAME(csrot) (BLASV(n),BLASP(v1p),BLASV(s1),
	BLASP(v2p),BLASV(s2),BLASV(c),BLASV(s)); 
  }
#ifdef ELAP
  template <> void DoGivensMult(
      float c, std::complex<float> s,
      const VectorView<std::complex<float> >& v1,
      const VectorView<std::complex<float> >& v2)
  {
    TMVAssert(v1.size()==v2.size());
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n = v1.size();
    int s1 = v1.step();
    int s2 = v2.step();
    std::complex<float>* v1p = v1.ptr();
    if (s1<0) v1p += (n-1)*s1;
    std::complex<float>* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*s2;
    LAPNAMEX(crot) (LAPV(n),LAPP(v1p),LAPV(s1),LAPP(v2p),LAPV(s2),
	LAPV(c),LAPP(&s)); 
  }
#endif // ELAP
#endif // DOUBLE
#endif // BLAS

  template <class T, class Tx> void GivensMult(
      RT c, T s,
      const VectorView<Tx>& v1, const VectorView<Tx>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    TMVAssert(!v1.SameAs(v2));

    if (v1.size() > 0 && s != T(0)) {
      if (ShouldReverse(v1.step(),v2.step()))
	GivensMult(c,s,v1.Reverse(),v2.Reverse());
      else if (v2.isconj()) 
	GivensMult(c,CONJ(s),v1.Conjugate(),v2.Conjugate());
#ifdef BLAS
      else if (v1.isconj()) {
	Vector<Tx> v1x = v1;
	DoGivensMult(c,s,v1x.View(),v2);
	v1 = v1x;
      }
#endif
      else DoGivensMult(c,s,v1,v2); 
    }
  }

  // 
  // Symmetric Givens Mult
  //
  template <class T, class Tx> void GivensHermMult(
      RT c, T s, Tx& d0, Tx& d1, Tx& e0) 
  {
    // [ d0 e0* ] = [  c  s ] [ d0 e0* ] [ c  -s ]
    // [ e0 d1  ]   [ -s* c ] [ e0 d1  ] [ s*  c ]
    // = [ c^2 d0 + 2c Re(s e0) + |s|^2 d1   cs(d1-d0) + c^2 e0* - s^2 e0    ]
    //   [ cs*(d1-d0) + c^2 e0 - s*^2 e0*    c^2 d1 - 2c Re(s e0) + |s|^2 d0 ]
    // (using c^2 = 1-|s|^2:)
    // d0' = d0 + 2c Re(s e0) + |s|^2 (d1-d0)
    // e0' = e0 - 2s* Re(s e0) + c s* (d1-d0)
    // d1' = d1 - 2c Re(s e0) - |s|^2 (d1-d0)

    RT s2 = NORM(s);
    RealType(Tx) Rese0 = REAL(s*e0);
    Tx d1md0 = d1-d0;
    Tx dd = RT(2)*c*Rese0 + s2*d1md0;
    d0 += dd;
    d1 -= dd;
    e0 += CONJ(s)*(c*d1md0 - RT(2)*Rese0);
  }
  template <class T, class Tx> void GivensSymMult(
      RT c, T s, Tx& d0, Tx& d1, Tx& e0) 
  {
    // [ d0 e0 ] = [  c  s ] [ d0 e0 ] [ c -s* ]
    // [ e0 d1 ]   [ -s* c ] [ e0 d1 ] [ s  c  ]
    // d0' = d0 + 2c (s e0) + s^2 (d1-d0)
    // e0' = e0 - 2|s|^2 e0 + c (s d1 - s* d0)
    // d1' = d1 - 2c (s e0) - s^2 (d1-d0)

    T s2 = s*s;
    Tx se0 = s*e0;
    Tx d1md0 = d1-d0;
    Tx dd = RT(2)*c*se0 + s2*d1md0;
    d0 += dd;
    d1 -= dd;
    if (IsReal(T()))
      e0 += s*(c*d1md0 - RT(2)*se0);
    else
      e0 += c*(s*d1-CONJ(s)*d0) - RT(2)*CONJ(s)*se0;
  }

#undef RT

#define InstFile "TMV_Givens.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


