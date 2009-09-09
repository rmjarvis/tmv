///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
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


//#define XDEBUG


#include "TMV_Blas.h"
#include "tmv/TMV_VectorArithFunc.h"
#include "tmv/TMV_Vector.h"

#ifdef XDEBUG
#include "tmv/TMV_VIt.h"
#include "tmv/TMV_VectorArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

  //
  // MultXV
  // 
  template <class T, class Tx> static void DoMultXV(
      const Tx x, const VectorView<T>& v)
  {
    //cout<<"DoMultXV: x = "<<x<<endl;
    TMVAssert(x!=Tx(0));
    TMVAssert(x!=Tx(1)); 
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>0);
    TMVAssert(v.ct() == NonConj);

    T* vptr = v.ptr();
    const int s = v.step();
    const int N = v.size();
    if (s == 1) {
      const int N1 = N/4;
      const int N2 = N-4*N1;
      if (N1) for(int i=N1;i>0;--i,vptr+=4) {
#ifdef TMVFLDEBUG
        TMVAssert(vptr >= v.first);
        TMVAssert(vptr+3 < v.last);
#endif
        *vptr *= x;
        vptr[1] *= x;
        vptr[2] *= x;
        vptr[3] *= x;
      }
      if (N2) for(int i=N2;i>0;--i,++vptr) {
#ifdef TMVFLDEBUG
        TMVAssert(vptr >= v.first);
        TMVAssert(vptr < v.last);
#endif
        *vptr *= x;
      }
    }
    else
      for(int i=N;i>0;--i,vptr+=s) {
#ifdef TMVFLDEBUG
        TMVAssert(vptr >= v.first);
        TMVAssert(vptr < v.last);
#endif
        *vptr *= x;
      }
  }

#ifdef BLAS
#ifdef INST_DOUBLE
  template <> void DoMultXV(const double x, 
      const VectorView<double>& v)
  {
    TMVAssert(x!=0.);
    TMVAssert(x!=1.);
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>=0);
    int n=v.size();
    int s=v.step();
    BLASNAME(dscal) (BLASV(n),BLASV(x),BLASP(v.ptr()),BLASV(s));
  }
  template <> void DoMultXV(const double x,
      const VectorView<std::complex<double> >& v)
  {
    TMVAssert(x!=0.);
    TMVAssert(x!=1.);
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>=0);
    TMVAssert(v.ct() == NonConj);
    int n=v.size();
    int s=v.step();
    BLASNAME(zdscal) (BLASV(n),BLASV(x),BLASP(v.ptr()),BLASV(s));
  }
  template <> void DoMultXV(const std::complex<double> x,
      const VectorView<std::complex<double> >& v)
  { 
    TMVAssert(x!=0.);
    TMVAssert(x!=1.);
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>=0);
    TMVAssert(v.ct() == NonConj);
    int n=v.size();
    int s=v.step();
    BLASNAME(zscal) (BLASV(n),BLASP(&x),BLASP(v.ptr()),BLASV(s));
  }
#endif
#ifdef INST_FLOAT
  template <> void DoMultXV(const float x, 
      const VectorView<float>& v)
  {
    TMVAssert(x!=0.F);
    TMVAssert(x!=1.F);
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>=0);
    int n=v.size();
    int s=v.step();
    BLASNAME(sscal) (BLASV(n),BLASV(x),BLASP(v.ptr()),BLASV(s));
  }
  template <> void DoMultXV(const float x,
      const VectorView<std::complex<float> >& v)
  {
    TMVAssert(x!=0.F);
    TMVAssert(x!=1.F);
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>=0);
    TMVAssert(v.ct() == NonConj);
    int n=v.size();
    int s=v.step();
    BLASNAME(csscal) (BLASV(n),BLASV(x),BLASP(v.ptr()),BLASV(s));
  }
  template <> void DoMultXV(const std::complex<float> x,
      const VectorView<std::complex<float> >& v)
  { 
    TMVAssert(x!=0.F);
    TMVAssert(x!=1.F);
    TMVAssert(v.size()>0);
    TMVAssert(v.step()>=0);
    TMVAssert(v.ct() == NonConj);
    int n=v.size();
    int s=v.step();
    BLASNAME(cscal) (BLASV(n),BLASP(&x),BLASP(v.ptr()),BLASV(s));
  }
#endif
#endif

  template <class T> void MultXV(const T x, 
      const VectorView<T>& v)
  { 
#ifdef XDEBUG
    //cout<<"MultXV: x = "<<x<<endl;
    Vector<T> v0 = v;
    Vector<T> vx = v;
    for(int i=0;i<int(vx.size());i++) vx(i) *= x;
#endif

    if (v.size() > 0 && x != T(1)) {
      if (v.step() < 0) MultXV(x,v.Reverse());
      else if (v.isconj()) MultXV(CONJ(x),v.Conjugate());
      else if (x == T(0)) v.Zero();
      else if (IsComplex(T()) && IMAG(x) == RealType(T)(0))
        if (v.step() == 1) DoMultXV(REAL(x),v.Flatten());
        else DoMultXV(REAL(x),v);
      else DoMultXV(x,v); 
    }

#ifdef XDEBUG
    Vector<T> diff(v.size());
    for(int i=0;i<int(v.size());i++) diff(i) = vx(i) - v(i);
    if (Norm(diff) > 0.001*MAX(RealType(T)(1),Norm(v))) {
      cerr<<"MultXV: x = "<<x<<endl;
      cerr<<"v = "<<TypeText(v)<<"  step "<<v.step()<<"  "<<v0<<endl;
      cerr<<"-> "<<v<<endl;
      cerr<<"vx = "<<vx<<endl;
      cerr<<"Norm(vx-v) = "<<Norm(vx-v)<<endl;
      abort();
    }
#endif
  }

  template <bool c1, class T, class Tx, class T1> static void DoMultXV(
      const Tx x, const GenVector<T1>& v1, const VectorView<T>& v2)
  {
    TMVAssert(v2.size()==v1.size());
    TMVAssert(v2.size()>0);
    TMVAssert(x!=Tx(0));
    TMVAssert(x!=Tx(1)); 
    TMVAssert(v1.step() != 0);
    TMVAssert(v2.step() != 0);
    TMVAssert(v2.ct() == NonConj);
    TMVAssert(v2.step() != -1);
    TMVAssert(v1.step() != -1 || v2.step() == 1);
    TMVAssert(v2.step() > 0 || v1.step() == 1);
    TMVAssert(IsReal(x) || IMAG(x) != RealType(Tx)(0));
    TMVAssert(!(v2.SameAs(v1)));
    TMVAssert(c1 == v1.isconj());

    const T1* v1ptr = v1.cptr();
    T* v2ptr = v2.ptr();
    const int s1 = v1.step();
    const int s2 = v2.step();
    const int N = v1.size();

    if (s1 == 1 && s2 == 1) {
      const int N1 = N/4;
      const int N2 = N-4*N1;
      if (N1) for(int i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
#ifdef TMVFLDEBUG
        TMVAssert(v2ptr >= v2.first);
        TMVAssert(v2ptr+3 < v2.last);
#endif
        *v2ptr = x * (c1 ? CONJ(*v1ptr) : (*v1ptr));
        v2ptr[1] = x * (c1 ? CONJ(v1ptr[1]) : v1ptr[1]);
        v2ptr[2] = x * (c1 ? CONJ(v1ptr[2]) : v1ptr[2]);
        v2ptr[3] = x * (c1 ? CONJ(v1ptr[3]) : v1ptr[3]);
      }
      if (N2) for(int i=N2;i>0;--i,v1ptr++,v2ptr++)  {
#ifdef TMVFLDEBUG
        TMVAssert(v2ptr >= v2.first);
        TMVAssert(v2ptr < v2.last);
#endif
        *v2ptr = x * (c1 ? CONJ(*v1ptr) : (*v1ptr));
      }
    }
    else
      for(int i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2)  {
#ifdef TMVFLDEBUG
        TMVAssert(v2ptr >= v2.first);
        TMVAssert(v2ptr < v2.last);
#endif
        *v2ptr = x * (c1 ? CONJ(*v1ptr) : (*v1ptr));
      }
  }

  template <class T, class T1> void MultXV(const T x, 
      const GenVector<T1>& v1, const VectorView<T>& v2)
  { 
    TMVAssert(v2.step() != 0 || v1.step() == 0 || v2.size() <= 1);
    TMVAssert(v1.size() == v2.size());

#ifdef XDEBUG
    Vector<T> vx = v1;
    for(int i=0;i<int(vx.size());i++) vx(i) *= x;
#endif

    if (v2.size() > 0) {
      if (v2.isconj()) MultXV(CONJ(x),v1.Conjugate(),v2.Conjugate());
      else if (v2.size() == 1) {
#ifdef TMVFLDEBUG
        TMVAssert(v2.ptr() >= v2.first);
        TMVAssert(v2.ptr() < v2.last);
#endif
        *v2.ptr() = x * (v1.isconj() ? CONJ(*v1.cptr()) : (*v1.cptr()));
      }
      else if (ShouldReverse(v1.step(),v2.step())) 
        MultXV(x,v1.Reverse(),v2.Reverse());
      else if (x == T(0)) v2.Zero();
      else if (x == T(1)) v2 = v1;
      else if (v1.step() == 0) v2.SetAllTo(x *
          (v1.isconj() ? CONJ(*v1.cptr()) : (*v1.cptr())));
      else if (v2.SameAs(v1)) MultXV(x,v2);
      else if (IsComplex(T()) && IMAG(x)==RealType(T)(0))
        if (IsComplex(T1()) && v2.isconj() == v1.isconj() &&
            (v1.step()==1 && v2.step()==1))
          DoMultXV<false>(REAL(x),v1.Flatten(),v2.Flatten());
        else if (v1.isconj()) DoMultXV<true>(REAL(x),v1,v2);
        else DoMultXV<false>(REAL(x),v1,v2);
      else if (v1.isconj()) DoMultXV<true>(x,v1,v2);
      else DoMultXV<false>(x,v1,v2);
    }

#ifdef XDEBUG
    Vector<T> diff(v2.size());
    for(int i=0;i<int(v2.size());i++) diff(i) = vx(i) - v2(i);
    if (Norm(diff) > 0.001*MAX(RealType(T)(1),Norm(v2))) {
      cerr<<"CopyMultXV: x = "<<x<<endl;
      cerr<<"v2 = "<<TypeText(v2)<<"  step "<<v2.step()<<endl;
      cerr<<"v1 = "<<TypeText(v1)<<"  step "<<v1.step()<<endl;
      cerr<<"-> "<<v2<<endl;
      cerr<<"vx = "<<vx<<endl;
      cerr<<"Norm(vx-v2) = "<<Norm(vx-v2)<<endl;
      abort();
    }
#endif
  }

  //
  // ElementProd
  //

  template <bool cx, bool cy, class T, class Ta, class Tx, class Ty>
  static void DoAddElementProd(const Ta alpha, const GenVector<Tx>& x,
      const GenVector<Ty>& y, const VectorView<T>& z)
  // zi += alpha * xi * yi 
  {
    TMVAssert(z.size() == x.size());
    TMVAssert(z.size() == y.size());
    TMVAssert(alpha != Ta(0));
    TMVAssert(z.size()>0);
    TMVAssert(z.ct() == NonConj);
    TMVAssert(z.step() != -1 || (x.step() == 1 && y.step() == 1));
    TMVAssert(x.step() != -1 || y.step() == 1 || z.step() == 1);
    TMVAssert(y.step() != -1 || x.step() == 1 || z.step() == 1);
    TMVAssert(z.step() >= 0 || x.step() == 1 || y.step() == 1);
    TMVAssert(cx == x.isconj());
    TMVAssert(cy == y.isconj());

    const Tx* xp = x.cptr();
    const Ty* yp = y.cptr();
    T* zp = z.ptr();
    const int sx = x.step();
    const int sy = y.step();
    const int sz = z.step();
    const int N = z.size();

    if (sx == 1 && sy == 1 && sz == 1) {
      const int N1 = N/4;
      const int N2 = N-4*N1;

      if (N1) {
        if (alpha == Ta(1))
          for(int i=N1;i>0;--i,xp+=4,yp+=4,zp+=4) {
#ifdef TMVFLDEBUG
            TMVAssert(zp >= z.first);
            TMVAssert(zp+3 < z.last);
#endif
            *zp += (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
            zp[1] += (cx?CONJ(xp[1]):xp[1]) * (cy?CONJ(yp[1]):yp[1]);
            zp[2] += (cx?CONJ(xp[2]):xp[2]) * (cy?CONJ(yp[2]):yp[2]);
            zp[3] += (cx?CONJ(xp[3]):xp[3]) * (cy?CONJ(yp[3]):yp[3]);
          }
        else if (alpha == Ta(-1))
          for(int i=N1;i>0;--i,xp+=4,yp+=4,zp+=4){
#ifdef TMVFLDEBUG
            TMVAssert(zp >= z.first);
            TMVAssert(zp+3 < z.last);
#endif
            *zp -= (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
            zp[1] -= (cx?CONJ(xp[1]):xp[1]) * (cy?CONJ(yp[1]):yp[1]);
            zp[2] -= (cx?CONJ(xp[2]):xp[2]) * (cy?CONJ(yp[2]):yp[2]);
            zp[3] -= (cx?CONJ(xp[3]):xp[3]) * (cy?CONJ(yp[3]):yp[3]);
          }
        else 
          for(int i=N1;i>0;--i,xp+=4,yp+=4,zp+=4) {
#ifdef TMVFLDEBUG
            TMVAssert(zp >= z.first);
            TMVAssert(zp+3 < z.last);
#endif
            *zp += alpha * (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
            zp[1] += alpha * (cx?CONJ(xp[1]):xp[1]) * (cy?CONJ(yp[1]):yp[1]);
            zp[2] += alpha * (cx?CONJ(xp[2]):xp[2]) * (cy?CONJ(yp[2]):yp[2]);
            zp[3] += alpha * (cx?CONJ(xp[3]):xp[3]) * (cy?CONJ(yp[3]):yp[3]);
          }
      }
      if (N2) {
        if (alpha == Ta(1))
          for(int i=N2;i>0;--i,++xp,++yp,++zp)  {
#ifdef TMVFLDEBUG
            TMVAssert(zp >= z.first);
            TMVAssert(zp < z.last);
#endif
            *zp += (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
          }
        else if (alpha == Ta(-1))
          for(int i=N2;i>0;--i,++xp,++yp,++zp) {
#ifdef TMVFLDEBUG
            TMVAssert(zp >= z.first);
            TMVAssert(zp < z.last);
#endif
            *zp -= (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
          }
        else 
          for(int i=N2;i>0;--i,++xp,++yp,++zp)  {
#ifdef TMVFLDEBUG
            TMVAssert(zp >= z.first);
            TMVAssert(zp < z.last);
#endif
            *zp += alpha * (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
          }
      }
    }
    else {
      if (alpha == Ta(1))
        for(int i=N;i>0;--i,xp+=sx,yp+=sy,zp+=sz) {
#ifdef TMVFLDEBUG
          TMVAssert(zp >= z.first);
          TMVAssert(zp < z.last);
#endif
          *zp += (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
        }
      else if (alpha == Ta(-1))
        for(int i=N;i>0;--i,xp+=sx,yp+=sy,zp+=sz) {
#ifdef TMVFLDEBUG
          TMVAssert(zp >= z.first);
          TMVAssert(zp < z.last);
#endif
          *zp -= (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
        }
      else 
        for(int i=N;i>0;--i,xp+=sx,yp+=sy,zp+=sz) {
#ifdef TMVFLDEBUG
          TMVAssert(zp >= z.first);
          TMVAssert(zp < z.last);
#endif
          *zp += alpha * (cx?CONJ(*xp):(*xp)) * (cy?CONJ(*yp):(*yp));
        }
    }
  }

  template <class T, class Tx, class Ty> void AddElementProd(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const VectorView<T>& z)
  {
    TMVAssert(z.size() == x.size());
    TMVAssert(z.size() == y.size());
#ifdef XDEBUG
    //cout<<"Start AddElProd:\n";
    //cout<<"steps = "<<x.step()<<"  "<<y.step()<<"  "<<z.step()<<endl;
    Vector<T> z0 = z;
    Vector<T> zx = z;
    for(int i=0;i<int(zx.size());i++) zx(i) += alpha*x(i)*y(i);
#endif

    if (z.size() > 0 && alpha != T(0)) {
      if (z.isconj()) 
        AddElementProd(CONJ(alpha),x.Conjugate(),y.Conjugate(),
            z.Conjugate());
      else if (
          (z.step()==-1 && !(x.step()==1 && y.step()==1)) ||
          (z.step()!=1 && (x.step()==-1 || (x.step()!=1 && y.step()==-1))) ||
          (z.step()<0 && !(x.step()==1 || y.step()==1)))
        AddElementProd(alpha,x.Reverse(),y.Reverse(),z.Reverse());
      else 
        if (IMAG(alpha) == RealType(T)(0))
          if (x.isconj())
            if (y.isconj())
              DoAddElementProd<true,true>(REAL(alpha),x,y,z);
            else
              DoAddElementProd<true,false>(REAL(alpha),x,y,z);
          else
            if (y.isconj())
              DoAddElementProd<false,true>(REAL(alpha),x,y,z);
            else
              DoAddElementProd<false,false>(REAL(alpha),x,y,z);
        else
          if (x.isconj())
            if (y.isconj())
              DoAddElementProd<true,true>(alpha,x,y,z);
            else
              DoAddElementProd<true,false>(alpha,x,y,z);
          else
            if (y.isconj())
              DoAddElementProd<false,true>(alpha,x,y,z);
            else
              DoAddElementProd<false,false>(alpha,x,y,z);
    }

#ifdef XDEBUG
    if (Norm(zx-z) > 0.001*(ABS(alpha)*(Norm(x)+Norm(y))+Norm(z0))) {
      cerr<<"AddElProd: alpha = "<<alpha<<endl;
      cerr<<"x = "<<TypeText(x)<<"  step "<<x.step()<<"  "<<x<<endl;
      cerr<<"y = "<<TypeText(y)<<"  step "<<y.step()<<"  "<<y<<endl;
      cerr<<"z = "<<TypeText(z)<<"  step "<<z.step()<<"  "<<z0<<endl;
      cerr<<"-> "<<z<<endl;
      cerr<<"zx = "<<zx<<endl;
      cerr<<"Norm(zx-z) = "<<Norm(zx-z)<<endl;
      abort();
    }
#endif
  }

  template <bool cx, class T, class Ta, class Tx> static void DoElementProd(
      const Ta alpha, const GenVector<Tx>& x, const VectorView<T>& y)
  // yi = alpha * xi * yi
  {
    TMVAssert(x.size() == y.size());
    TMVAssert(y.size()>0);
    TMVAssert(alpha != Ta(0));
    TMVAssert(y.ct() == NonConj);
    TMVAssert(y.step() != 0 || x.step() == 0 || y.size() <= 1);
    TMVAssert(y.step() != -1);
    TMVAssert(x.step() != -1 || y.step() == 1);
    TMVAssert(y.step() > 0 || x.step() == 1);

    const Tx* xp = x.cptr();
    T* yp = y.ptr();
    const int sx = x.step();
    const int sy = y.step();
    const int N = y.size();

    if (sx == 1 && sy == 1) {
      const int N1 = N/4;
      const int N2 = N-4*N1;
      if (N1) {
        if (alpha == Ta(1))
          for(int i=N1;i>0;--i,xp+=4,yp+=4) {
#ifdef TMVFLDEBUG
            TMVAssert(yp >= y.first);
            TMVAssert(yp+3 < y.last);
#endif
            *yp *= cx ? CONJ(*xp) : (*xp);
            yp[1] *= cx ? CONJ(xp[1]) : xp[1];
            yp[2] *= cx ? CONJ(xp[2]) : xp[2];
            yp[3] *= cx ? CONJ(xp[3]) : xp[3];
          }
        else 
          for(int i=N1;i>0;--i,xp+=4,yp+=4) {
#ifdef TMVFLDEBUG
            TMVAssert(yp >= y.first);
            TMVAssert(yp+3 < y.last);
#endif
            *yp *= alpha * (cx ? CONJ(*xp) : (*xp));
            yp[1] *= alpha * (cx ? CONJ(xp[1]) : xp[1]);
            yp[2] *= alpha * (cx ? CONJ(xp[2]) : xp[2]);
            yp[3] *= alpha * (cx ? CONJ(xp[3]) : xp[3]);
          }
      }
      if (N2) {
        if (alpha == Ta(1))
          for(int i=N2;i>0;--i,++xp,++yp) {
#ifdef TMVFLDEBUG
            TMVAssert(yp >= y.first);
            TMVAssert(yp < y.last);
#endif
            *yp *= cx ? CONJ(*xp) : (*xp);
          }
        else 
          for(int i=N2;i>0;--i,++xp,++yp) {
#ifdef TMVFLDEBUG
            TMVAssert(yp >= y.first);
            TMVAssert(yp < y.last);
#endif
            *yp *= alpha * (cx ? CONJ(*xp) : (*xp));
          }
      }
    }
    else {
      if (alpha == Ta(1))
        for(int i=N;i>0;--i,xp+=sx,yp+=sy) {
#ifdef TMVFLDEBUG
          TMVAssert(yp >= y.first);
          TMVAssert(yp < y.last);
#endif
          *yp *= cx ? CONJ(*xp) : (*xp);
        }
      else 
        for(int i=N;i>0;--i,xp+=sx,yp+=sy) {
#ifdef TMVFLDEBUG
          TMVAssert(yp >= y.first);
          TMVAssert(yp < y.last);
#endif
          *yp *= alpha * (cx ? CONJ(*xp) : (*xp));
        }
    }
  }

  template <class T, class Tx> void ElementProd(const T alpha,
      const GenVector<Tx>& x, const VectorView<T>& y)
  {
    TMVAssert(x.size() == y.size());
    TMVAssert(y.step() != 0 || x.step() == 0 || y.size() <= 1);

#ifdef XDEBUG
    Vector<T> y0 = y;
    Vector<T> yx = y;
    for(int i=0;i<int(yx.size());i++) yx(i) *= alpha*x(i);
#endif

    if (y.size() > 0 && alpha != T(0)) {
      if (y.isconj()) ElementProd(CONJ(alpha),x.Conjugate(),y.Conjugate());
      else if (y.size() == 1) {
#ifdef TMVFLDEBUG
        TMVAssert(y.ptr() >= y.first);
        TMVAssert(y.ptr() < y.last);
#endif
        *y.ptr() *= alpha * (x.isconj() ? CONJ(*x.cptr()) : (*x.cptr()));
      } else if (ShouldReverse(x.step(),y.step())) {
        ElementProd(alpha,x.Reverse(),y.Reverse());
      } else if (IMAG(alpha) == 0) {
        if (x.isconj()) 
          DoElementProd<true>(REAL(alpha),x,y);
        else 
          DoElementProd<false>(REAL(alpha),x,y);
      } else {
        if (x.isconj()) 
          DoElementProd<true>(alpha,x,y);
        else 
          DoElementProd<false>(alpha,x,y);
      }
    }

#ifdef XDEBUG
    if (Norm(yx-y) > 0.001*(ABS(alpha)*Norm(x)+Norm(y0))) {
      cerr<<"AddElProd: alpha = "<<alpha<<endl;
      cerr<<"x = "<<TypeText(x)<<"  step "<<x.step()<<"  "<<x<<endl;
      cerr<<"y = "<<TypeText(y)<<"  step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"-> "<<y<<endl;
      cerr<<"yx = "<<yx<<endl;
      cerr<<"Norm(yx-y) = "<<Norm(yx-y)<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_MultXV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


