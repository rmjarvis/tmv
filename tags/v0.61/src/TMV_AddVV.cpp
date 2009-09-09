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
#include "TMV_VectorArithFunc.h"
#include "TMV_Vector.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_VIt.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

  // 
  // AddVV
  //

  template <bool c1, class T, class Tx, class T1> static void NonBlasAddVV(
      const Tx x, const GenVector<T1>& v1, const VectorView<T>& v2)
  {
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v2.size()>0);
    TMVAssert(x != Tx(0));
    TMVAssert(v2.ct() == NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(IsReal(x) || IMAG(x) != RealType(Tx)(0));
    TMVAssert(c1 == v1.isconj());

    const T1* v1ptr = v1.cptr();
    T* v2ptr = v2.ptr();
    const int s1 = v1.step();
    const int s2 = v2.step();
    const int N = v2.size();

    if (s1 == 1 && s2 == 1) {
      const int N1 = N/4;
      const int N2 = N-4*N1;
      if (N1) {
	if (x == Tx(1))        
	  for(int i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
#ifdef TMVFLDEBUG
	    TMVAssert(v2ptr >= v2.first);
	    TMVAssert(v2ptr+3 < v2.last);
#endif
	    *v2ptr += (c1 ? CONJ(*v1ptr) : (*v1ptr));
	    v2ptr[1] += (c1 ? CONJ(v1ptr[1]) : v1ptr[1]);
	    v2ptr[2] += (c1 ? CONJ(v1ptr[2]) : v1ptr[2]);
	    v2ptr[3] += (c1 ? CONJ(v1ptr[3]) : v1ptr[3]);
	  }
	else if (x == Tx(-1))
	  for(int i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
#ifdef TMVFLDEBUG
	    TMVAssert(v2ptr >= v2.first);
	    TMVAssert(v2ptr+3 < v2.last);
#endif
	    *v2ptr -= (c1 ? CONJ(*v1ptr) : (*v1ptr));
	    v2ptr[1] -= (c1 ? CONJ(v1ptr[1]) : v1ptr[1]);
	    v2ptr[2] -= (c1 ? CONJ(v1ptr[2]) : v1ptr[2]);
	    v2ptr[3] -= (c1 ? CONJ(v1ptr[3]) : v1ptr[3]);
	  }
	else
	  for(int i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
#ifdef TMVFLDEBUG
	    TMVAssert(v2ptr >= v2.first);
	    TMVAssert(v2ptr+3 < v2.last);
#endif
	    *v2ptr += x * (c1 ? CONJ(*v1ptr) : (*v1ptr));
	    v2ptr[1] += x * (c1 ? CONJ(v1ptr[1]) : v1ptr[1]);
	    v2ptr[2] += x * (c1 ? CONJ(v1ptr[2]) : v1ptr[2]);
	    v2ptr[3] += x * (c1 ? CONJ(v1ptr[3]) : v1ptr[3]);
	  }
      }
      if (N2) {
	if (x == Tx(1))        
	  for(int i=N2;i>0;--i,++v1ptr,++v2ptr) {
#ifdef TMVFLDEBUG
	    TMVAssert(v2ptr >= v2.first);
	    TMVAssert(v2ptr < v2.last);
#endif
	    *v2ptr += (c1 ? CONJ(*v1ptr) : (*v1ptr));
	  }
	else if (x == Tx(-1))
	  for(int i=N2;i>0;--i,++v1ptr,++v2ptr) {
#ifdef TMVFLDEBUG
	    TMVAssert(v2ptr >= v2.first);
	    TMVAssert(v2ptr < v2.last);
#endif
	    *v2ptr -= (c1 ? CONJ(*v1ptr) : (*v1ptr));
	  }
	else
	  for(int i=N2;i>0;--i,++v1ptr,++v2ptr) {
#ifdef TMVFLDEBUG
	    TMVAssert(v2ptr >= v2.first);
	    TMVAssert(v2ptr < v2.last);
#endif
	    *v2ptr += x * (c1 ? CONJ(*v1ptr) : (*v1ptr));
	  }
      }
    } else {
      if (x == Tx(1))        
	for(int i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) {
#ifdef TMVFLDEBUG
	  TMVAssert(v2ptr >= v2.first);
	  TMVAssert(v2ptr < v2.last);
#endif
	  *v2ptr += (c1 ? CONJ(*v1ptr) : (*v1ptr));
	}
      else if (x == Tx(-1))
	for(int i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) {
#ifdef TMVFLDEBUG
	  TMVAssert(v2ptr >= v2.first);
	  TMVAssert(v2ptr < v2.last);
#endif
	  *v2ptr -= (c1 ? CONJ(*v1ptr) : (*v1ptr));
	}
      else
	for(int i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2) {
#ifdef TMVFLDEBUG
	  TMVAssert(v2ptr >= v2.first);
	  TMVAssert(v2ptr < v2.last);
#endif
	  *v2ptr += x * (c1 ? CONJ(*v1ptr) : (*v1ptr));
	}
    }
  }

  template <class T, class T1> static void DoAddVV(
      const T x, const GenVector<T1>& v1, const VectorView<T>& v2)
  {
    if (IMAG(x) == RealType(T)(0))
      if (v1.isconj()) NonBlasAddVV<true>(REAL(x),v1,v2); 
      else NonBlasAddVV<false>(REAL(x),v1,v2); 
    else
      if (v1.isconj()) NonBlasAddVV<true>(x,v1,v2); 
      else NonBlasAddVV<false>(x,v1,v2); 
  }

#ifdef BLAS
#ifdef INST_DOUBLE
  template <> void DoAddVV(const double x,
      const GenVector<double>& v1, const VectorView<double>& v2)
  { 
    TMVAssert(v1.size() == v1.size()); 
    TMVAssert(v2.size()>0);
    TMVAssert(x != double(0));
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    const double* v1p = v1.cptr();
    if (s1<0) v1p += (n-1)*s1;
    double* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*s2;
    BLASNAME(daxpy) (BLASV(n),BLASV(x),BLASP(v1p),BLASV(s1),
	BLASP(v2p),BLASV(s2));
  }
  template <> void DoAddVV(const std::complex<double> x, 
      const GenVector<std::complex<double> >& v1, 
      const VectorView<std::complex<double> >& v2)
  { 
    TMVAssert(v1.size() == v1.size()); 
    TMVAssert(v2.size()>0);
    TMVAssert(x != double(0));
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    const std::complex<double>* v1p = v1.cptr();
    if (s1<0) v1p += (n-1)*s1;
    std::complex<double>* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*s2;
    BLASNAME(zaxpy) (BLASV(n),BLASP(&x),BLASP(v1p),BLASV(s1),
	BLASP(v2p),BLASV(s2));
  }
  template <> void DoAddVV(const std::complex<double> x, 
      const GenVector<double>& v1, 
      const VectorView<std::complex<double> >& v2)
  { 
    TMVAssert(v1.size() == v1.size()); 
    TMVAssert(v2.size()>0);
    TMVAssert(x != double(0));
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    double xr = REAL(x);
    double xi = IMAG(x);
    int n=v2.size();
    int s1=v1.step();
    int s2=2*v2.step();
    const double* v1p = v1.cptr();
    if (s1<0) v1p += (n-1)*s1;
    std::complex<double>* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*v2.step();
    if (xr != double(0))
      BLASNAME(daxpy) (BLASV(n),BLASV(xr),BLASP(v1p),BLASV(s1),
	  BLASP((double*)v2p),BLASV(s2));
    if (xi != double(0))
      BLASNAME(daxpy) (BLASV(n),BLASV(xi),BLASP(v1p),BLASV(s1),
	  BLASP((double*)v2p+1),BLASV(s2));
  }
#endif
#ifdef INST_FLOAT
  template <> void DoAddVV(const float x,
      const GenVector<float>& v1, const VectorView<float>& v2)
  { 
    TMVAssert(v1.size() == v1.size()); 
    TMVAssert(v2.size()>0);
    TMVAssert(x != float(0));
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    const float* v1p = v1.cptr();
    if (s1<0) v1p += (n-1)*s1;
    float* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*s2;
    BLASNAME(saxpy) (BLASV(n),BLASV(x),BLASP(v1p),BLASV(s1),
	BLASP(v2p),BLASV(s2));
  }
  template <> void DoAddVV(const std::complex<float> x, 
      const GenVector<std::complex<float> >& v1, 
      const VectorView<std::complex<float> >& v2)
  { 
    TMVAssert(v1.size() == v1.size()); 
    TMVAssert(v2.size()>0);
    TMVAssert(x != float(0));
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    int n=v2.size();
    int s1=v1.step();
    int s2=v2.step();
    const std::complex<float>* v1p = v1.cptr();
    if (s1<0) v1p += (n-1)*s1;
    std::complex<float>* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*s2;
    BLASNAME(caxpy) (BLASV(n),BLASP(&x),BLASP(v1p),BLASV(s1),
	BLASP(v2p),BLASV(s2));
  }
  template <> void DoAddVV(const std::complex<float> x, 
      const GenVector<float>& v1, 
      const VectorView<std::complex<float> >& v2)
  { 
    TMVAssert(v1.size() == v1.size()); 
    TMVAssert(v2.size()>0);
    TMVAssert(x != float(0));
    TMVAssert(v1.ct()==NonConj);
    TMVAssert(v2.ct()==NonConj);
    float xr = REAL(x);
    float xi = IMAG(x);
    int n=v2.size();
    int s1=v1.step();
    int s2=2*v2.step();
    const float* v1p = v1.cptr();
    if (s1<0) v1p += (n-1)*s1;
    std::complex<float>* v2p = v2.ptr();
    if (s2<0) v2p += (n-1)*v2.step();
    if (xr != float(0))
      BLASNAME(saxpy) (BLASV(n),BLASV(xr),BLASP(v1p),BLASV(s1),
	  BLASP((float*)v2p),BLASV(s2));
    if (xi != float(0))
      BLASNAME(saxpy) (BLASV(n),BLASV(xi),BLASP(v1p),BLASV(s1),
	  BLASP((float*)v2p+1),BLASV(s2));
  }
#endif
#endif // BLAS

  template <class T, class T1> void AddVV(
      const T x, const GenVector<T1>& v1, const VectorView<T>& v2)
  { 
    TMVAssert(v1.size() == v2.size()); 
    TMVAssert(v2.step() != 0 || v1.step() == 0 || v2.size() <= 1);

#ifdef XDEBUG
    //std::cout<<"Start AddVV: x = "<<x<<endl;
    //std::cout<<"v1 = "<<Type(v1)<<"  "<<v1<<endl;
    //std::cout<<"v2 = "<<Type(v2)<<"  "<<v2<<endl;
    Vector<T> v0 = v2;
    Vector<T> vx = v2;
    for(int i=0;i<int(v2.size());i++) vx(i) += x * v1(i);
#endif

    if (v2.size() > 0 && x != T(0)) {
      if (v2.SameAs(v1)) {
	if (x == T(-1)) v2.Zero();
	else MultXV(x+T(1),v2);
      }
      else if (v2.isconj()) AddVV(CONJ(x),v1.Conjugate(),v2.Conjugate());
      else if (ShouldReverse(v1.step(),v2.step())) 
	AddVV(x,v1.Reverse(),v2.Reverse());
#ifdef BLAS
      else if (v1.isconj())
	DoAddVV(x,Vector<T1>(v1),v2);
#endif
      else
	DoAddVV(x,v1,v2);
    }

#ifdef XDEBUG
    Vector<T> diff(v2.size());
    for(int i=0;i<int(v2.size());i++) diff(i) = vx(i) - v2(i);
    if (Norm(diff) > 0.001*ABS(x)*(Norm(v1)+Norm(v0))) {
      cerr<<"AddVV: x = "<<x<<endl;
      cerr<<"v1 = "<<Type(v1)<<"  step "<<v1.step()<<"  "<<v1<<endl;
      cerr<<"v2 = "<<Type(v2)<<"  step "<<v2.step()<<"  "<<v0<<endl;
      cerr<<"-> "<<v2<<endl;
      cerr<<"vx = "<<vx<<endl;
      cerr<<"Norm(vx-v2) = "<<Norm(diff)<<endl;
      abort();
    }
#endif
  }

  template <class T, class T1, class T2> void AddVV(
      const T x1, const GenVector<T1>& v1,
      const T x2, const GenVector<T2>& v2, const VectorView<T>& v3)
  { 
    TMVAssert(v2.size() == v1.size()); 
    TMVAssert(v3.size() == v1.size()); 
    TMVAssert(v3.step() != 0 || (v1.step() == 0 && v2.step() == 0) || 
	v3.size() <= 1);

#ifdef XDEBUG
    //cout<<"Start AddVV: x1 = "<<x1<<", x2 = "<<x2<<endl;
    //cout<<"v1 = "<<Type(v1)<<"  "<<v1<<endl;
    //cout<<"v2 = "<<Type(v2)<<"  "<<v2<<endl;
    //cout<<"v3 = "<<Type(v3)<<"  "<<v3<<endl;
    Vector<T> v10 = v1;
    Vector<T> v20 = v2;
    Vector<T> vx = v2;
    for(int i=0;i<int(vx.size());i++) {
      vx(i) *= x2;
      vx(i) += x1*v1(i);
    }
#endif

    if (v3.size() > 0) {
      if (SameStorage(v1,v3)) {
	if (SameStorage(v2,v3)) {
	  Vector<T> temp(v3.size());
	  MultXV(x2,v2,temp.View());
	  AddVV(x1,v1,temp.View());
	  v3 = temp;
	} else {
	  MultXV(x1,v1,v3);
	  AddVV(x2,v2,v3);
	}
      } else {
	MultXV(x2,v2,v3);
	AddVV(x1,v1,v3);
      }
    }

#ifdef XDEBUG
    Vector<T> diff(v3.size());
    for(int i=0;i<int(v3.size());i++) diff(i) = vx(i) - v3(i);
    if (Norm(diff) > 0.001*(
	  x1==T(0)?RealType(T)(0):(ABS(x1)*Norm(v10)) +
	  x2==T(0)?RealType(T)(0):(ABS(x2)*Norm(v20)) )) {
      cerr<<"AddVV: x1 = "<<x1<<"  x2 = "<<x2<<endl;
      cerr<<"v1 = "<<Type(v1)<<"  step "<<v1.step()<<"  "<<v10<<endl;
      cerr<<"v2 = "<<Type(v2)<<"  step "<<v2.step()<<"  "<<v20<<endl;
      cerr<<"-> "<<v3<<endl;
      cerr<<"vx = "<<vx<<endl;
      cerr<<"Norm(vx-v3) = "<<Norm(diff)<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_AddVV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


