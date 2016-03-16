///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
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
    template <class T, class Tx> 
    static void DoMultXV(const Tx x, VectorView<T> v)
    {
        //cout<<"DoMultXV: x = "<<x<<endl;
        TMVAssert(x!=Tx(0));
        TMVAssert(x!=Tx(1)); 
        TMVAssert(v.size()>0);
        TMVAssert(v.step()>=0);
        TMVAssert(v.ct() == NonConj);

        T* vptr = v.ptr();
        const ptrdiff_t s = v.step();
        const ptrdiff_t N = v.size();
        if (s == 1) {
            const ptrdiff_t N1 = N/4;
            const ptrdiff_t N2 = N-4*N1;
            if (N1) for(ptrdiff_t i=N1;i>0;--i,vptr+=4) {
#ifdef TMVFLDEBUG
                TMVAssert(vptr >= v._first);
                TMVAssert(vptr+3 < v._last);
#endif
                *vptr *= x;
                vptr[1] *= x;
                vptr[2] *= x;
                vptr[3] *= x;
            }
            if (N2) for(ptrdiff_t i=N2;i>0;--i,++vptr) {
#ifdef TMVFLDEBUG
                TMVAssert(vptr >= v._first);
                TMVAssert(vptr < v._last);
#endif
                *vptr *= x;
            }
        } else {
            for(ptrdiff_t i=N;i>0;--i,vptr+=s) {
#ifdef TMVFLDEBUG
                TMVAssert(vptr >= v._first);
                TMVAssert(vptr < v._last);
#endif
                *vptr *= x;
            }
        }
    }

#ifdef BLAS
#ifdef INST_DOUBLE
    template <> 
    void DoMultXV(const double x, VectorView<double> v)
    {
        int n=v.size();
        int s=v.step();
        BLASNAME(dscal) (BLASV(n),BLASV(x),BLASP(v.ptr()),BLASV(s));
    }
    template <> 
    void DoMultXV(const double x, VectorView<std::complex<double> > v)
    {
        int n=v.size();
        int s=v.step();
        BLASNAME(zdscal) (BLASV(n),BLASV(x),BLASP(v.ptr()),BLASV(s));
    }
    template <> 
    void DoMultXV(
        const std::complex<double> x,
        VectorView<std::complex<double> > v)
    { 
        int n=v.size();
        int s=v.step();
        BLASNAME(zscal) (BLASV(n),BLASP(&x),BLASP(v.ptr()),BLASV(s));
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void DoMultXV(const float x, VectorView<float> v)
    {
        int n=v.size();
        int s=v.step();
        BLASNAME(sscal) (BLASV(n),BLASV(x),BLASP(v.ptr()),BLASV(s));
    }
    template <> 
    void DoMultXV(const float x, VectorView<std::complex<float> > v)
    {
        int n=v.size();
        int s=v.step();
        BLASNAME(csscal) (BLASV(n),BLASV(x),BLASP(v.ptr()),BLASV(s));
    }
    template <> 
    void DoMultXV(
        const std::complex<float> x,
        VectorView<std::complex<float> > v)
    { 
        int n=v.size();
        int s=v.step();
        BLASNAME(cscal) (BLASV(n),BLASP(&x),BLASP(v.ptr()),BLASV(s));
    }
#endif
#endif

    template <class T> 
    void MultXV(const T x, VectorView<T> v)
    { 
#ifdef XDEBUG
        Vector<T> v0 = v;
        Vector<T> vx = v;
        for(ptrdiff_t i=0;i<vx.size();i++) vx(i) *= x;
#endif

        if (v.size() > 0 && x != T(1)) {
            if (v.step() < 0) MultXV(x,v.reverse());
            else if (v.isconj()) MultXV(TMV_CONJ(x),v.conjugate());
            else if (x == T(0)) v.setZero();
            else if (isComplex(T()) && TMV_IMAG(x) == TMV_RealType(T)(0))
                if (v.step() == 1) DoMultXV(TMV_REAL(x),v.flatten());
                else DoMultXV(TMV_REAL(x),v);
            else DoMultXV(x,v); 
        }

#ifdef XDEBUG
        Vector<T> diff(v.size());
        for(ptrdiff_t i=0;i<v.size();i++) diff(i) = vx(i) - v(i);
        if (!(Norm(diff) <= 0.001*TMV_MAX(TMV_RealType(T)(1),Norm(v)))) {
            cerr<<"MultXV: x = "<<x<<endl;
            cerr<<"v = "<<TMV_Text(v)<<"  step "<<v.step()<<"  "<<v0<<endl;
            cerr<<"-> "<<v<<endl;
            cerr<<"vx = "<<vx<<endl;
            cerr<<"Norm(vx-v) = "<<Norm(vx-v)<<endl;
            abort();
        }
#endif
    }

    template <bool c1, class T, class Tx, class T1> 
    static void DoMultXV(
        const Tx x, const GenVector<T1>& v1, VectorView<T> v2)
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
        TMVAssert(isReal(x) || TMV_IMAG(x) != TMV_RealType(Tx)(0));
        TMVAssert(!(v2.isSameAs(v1)));
        TMVAssert(c1 == v1.isconj());

        const T1* v1ptr = v1.cptr();
        T* v2ptr = v2.ptr();
        const ptrdiff_t s1 = v1.step();
        const ptrdiff_t s2 = v2.step();
        const ptrdiff_t N = v1.size();

        if (s1 == 1 && s2 == 1) {
            const ptrdiff_t N1 = N/4;
            const ptrdiff_t N2 = N-4*N1;
            if (N1) for(ptrdiff_t i=N1;i>0;--i,v1ptr+=4,v2ptr+=4) {
#ifdef TMVFLDEBUG
                TMVAssert(v2ptr >= v2._first);
                TMVAssert(v2ptr+3 < v2._last);
#endif
                *v2ptr = x * (c1 ? TMV_CONJ(*v1ptr) : (*v1ptr));
                v2ptr[1] = x * (c1 ? TMV_CONJ(v1ptr[1]) : v1ptr[1]);
                v2ptr[2] = x * (c1 ? TMV_CONJ(v1ptr[2]) : v1ptr[2]);
                v2ptr[3] = x * (c1 ? TMV_CONJ(v1ptr[3]) : v1ptr[3]);
            }
            if (N2) for(ptrdiff_t i=N2;i>0;--i,v1ptr++,v2ptr++)  {
#ifdef TMVFLDEBUG
                TMVAssert(v2ptr >= v2._first);
                TMVAssert(v2ptr < v2._last);
#endif
                *v2ptr = x * (c1 ? TMV_CONJ(*v1ptr) : (*v1ptr));
            }
        } else
            for(ptrdiff_t i=N;i>0;--i,v1ptr+=s1,v2ptr+=s2)  {
#ifdef TMVFLDEBUG
                TMVAssert(v2ptr >= v2._first);
                TMVAssert(v2ptr < v2._last);
#endif
                *v2ptr = x * (c1 ? TMV_CONJ(*v1ptr) : (*v1ptr));
            }
    }

    template <class T, class T1> 
    void MultXV(const T x, const GenVector<T1>& v1, VectorView<T> v2)
    { 
        TMVAssert(v2.step() != 0 || v1.step() == 0 || v2.size() <= 1);
        TMVAssert(v1.size() == v2.size());

#ifdef XDEBUG
        Vector<T> vx = v1;
        for(ptrdiff_t i=0;i<vx.size();i++) vx(i) *= x;
#endif

        if (v2.size() > 0) {
            if (v2.isconj()) {
                MultXV(TMV_CONJ(x),v1.conjugate(),v2.conjugate());
            } else if (v2.size() == 1) {
#ifdef TMVFLDEBUG
                TMVAssert(v2.ptr() >= v2._first);
                TMVAssert(v2.ptr() < v2._last);
#endif
                *v2.ptr() = x * 
                    ( v1.isconj() ?
                      TMV_CONJ(*v1.cptr()) :
                      (*v1.cptr())
                    );
            } else if (shouldReverse(v1.step(),v2.step())) {
                MultXV(x,v1.reverse(),v2.reverse());
            } else if (x == T(0)) {
                v2.setZero();
            } else if (x == T(1)) {
                v2 = v1;
            } else if (v1.step() == 0) {
                v2.setAllTo(x * (
                        v1.isconj() ?
                        TMV_CONJ(*v1.cptr()) :
                        (*v1.cptr()))
                    );
            } else if (v2.isSameAs(v1)) {
                MultXV(x,v2);
            } else if (isComplex(T()) && TMV_IMAG(x)==TMV_RealType(T)(0)) {
                if (isComplex(T1()) && v2.isconj() == v1.isconj() &&
                    (v1.step()==1 && v2.step()==1))
                    DoMultXV<false>(TMV_REAL(x),v1.flatten(),v2.flatten());
                else if (v1.isconj()) DoMultXV<true>(TMV_REAL(x),v1,v2);
                else DoMultXV<false>(TMV_REAL(x),v1,v2);
            } else if (v1.isconj()) {
                DoMultXV<true>(x,v1,v2);
            } else {
                DoMultXV<false>(x,v1,v2);
            }
        }

#ifdef XDEBUG
        Vector<T> diff(v2.size());
        for(ptrdiff_t i=0;i<v2.size();i++) diff(i) = vx(i) - v2(i);
        if (!(Norm(diff) <= 0.001*TMV_MAX(TMV_RealType(T)(1),Norm(v2)))) {
            cerr<<"CopyMultXV: x = "<<x<<endl;
            cerr<<"v2 = "<<TMV_Text(v2)<<"  step "<<v2.step()<<endl;
            cerr<<"v1 = "<<TMV_Text(v1)<<"  step "<<v1.step()<<endl;
            cerr<<"-> "<<v2<<endl;
            cerr<<"vx = "<<vx<<endl;
            cerr<<"Norm(vx-v2) = "<<Norm(vx-v2)<<endl;
            abort();
        }
#endif
    }

    //
    // ElemMult
    //

    template <bool yn>
    struct Maybe  // yn = true
    {
        template <class T1, class T2>
        static void add(T1& a, const T2& b) { a += b; }
        template <class T>
        static T conj(const T& a) { return TMV_CONJ(a); }
    };
    template <>
    struct Maybe<false>  // yn = false
    {
        template <class T1, class T2>
        static void add(T1& a, const T2& b) { a = b; }
        template <class T>
        static const T& conj(const T& a) { return a; }
    };

    template <bool add, bool cx, bool cy, class T, class Ta, class Tx, class Ty>
    static void DoElemMultVV(
        const Ta alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, VectorView<T> z)
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
        const ptrdiff_t sx = x.step();
        const ptrdiff_t sy = y.step();
        const ptrdiff_t sz = z.step();
        const ptrdiff_t N = z.size();

        if (sx == 1 && sy == 1 && sz == 1) {
            const ptrdiff_t N1 = N/4;
            const ptrdiff_t N2 = N-4*N1;

            if (N1) {
                if (alpha == Ta(1)) {
                    for(ptrdiff_t i=N1;i>0;--i,xp+=4,yp+=4,zp+=4) {
#ifdef TMVFLDEBUG
                        TMVAssert(zp >= z._first);
                        TMVAssert(zp+3 < z._last);
#endif
                        Maybe<add>::add(
                            *zp,
                            Maybe<cx>::conj(*xp) * Maybe<cy>::conj(*yp));
                        Maybe<add>::add(
                            zp[1],
                            Maybe<cx>::conj(xp[1]) * Maybe<cy>::conj(yp[1]));
                        Maybe<add>::add(
                            zp[2],
                            Maybe<cx>::conj(xp[2]) * Maybe<cy>::conj(yp[2]));
                        Maybe<add>::add(
                            zp[3],
                            Maybe<cx>::conj(xp[3]) * Maybe<cy>::conj(yp[3]));
                    }
                } else {
                    for(ptrdiff_t i=N1;i>0;--i,xp+=4,yp+=4,zp+=4) {
#ifdef TMVFLDEBUG
                        TMVAssert(zp >= z._first);
                        TMVAssert(zp+3 < z._last);
#endif
                        Maybe<add>::add(
                            *zp, alpha * 
                            Maybe<cx>::conj(*xp) * Maybe<cy>::conj(*yp));
                        Maybe<add>::add(
                            zp[1], alpha * 
                            Maybe<cx>::conj(xp[1]) * Maybe<cy>::conj(yp[1]));
                        Maybe<add>::add(
                            zp[2], alpha * 
                            Maybe<cx>::conj(xp[2]) * Maybe<cy>::conj(yp[2]));
                        Maybe<add>::add(
                            zp[3], alpha * 
                            Maybe<cx>::conj(xp[3]) * Maybe<cy>::conj(yp[3]));
                    }
                }
            }
            if (N2) {
                if (alpha == Ta(1)) {
                    for(ptrdiff_t i=N2;i>0;--i,++xp,++yp,++zp)  {
#ifdef TMVFLDEBUG
                        TMVAssert(zp >= z._first);
                        TMVAssert(zp < z._last);
#endif
                        Maybe<add>::add(
                            *zp,
                            Maybe<cx>::conj(*xp) * Maybe<cy>::conj(*yp));
                    }
                } else {
                    for(ptrdiff_t i=N2;i>0;--i,++xp,++yp,++zp)  {
#ifdef TMVFLDEBUG
                        TMVAssert(zp >= z._first);
                        TMVAssert(zp < z._last);
#endif
                        Maybe<add>::add(
                            *zp, alpha * 
                            Maybe<cx>::conj(*xp) * Maybe<cy>::conj(*yp));
                    }
                }
            }
        } else {
            if (alpha == Ta(1)) {
                for(ptrdiff_t i=N;i>0;--i,xp+=sx,yp+=sy,zp+=sz) {
#ifdef TMVFLDEBUG
                    TMVAssert(zp >= z._first);
                    TMVAssert(zp < z._last);
#endif
                    Maybe<add>::add(
                        *zp, Maybe<cx>::conj(*xp) * Maybe<cy>::conj(*yp));
                }
            } else {
                for(ptrdiff_t i=N;i>0;--i,xp+=sx,yp+=sy,zp+=sz) {
#ifdef TMVFLDEBUG
                    TMVAssert(zp >= z._first);
                    TMVAssert(zp < z._last);
#endif
                    Maybe<add>::add(
                        *zp, alpha * 
                        Maybe<cx>::conj(*xp) * Maybe<cy>::conj(*yp));
                }
            }
        }
    }

    template <bool add, class T, class Tx, class Ty> 
    void ElemMultVV(
        const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
        VectorView<T> z)
    {
        typedef typename Traits<T>::real_type RT;
        TMVAssert(z.size() == x.size());
        TMVAssert(z.size() == y.size());
#ifdef XDEBUG
        Vector<T> z0 = z;
        Vector<T> zx(z.size());
        for(ptrdiff_t i=0;i<z.size();i++) zx(i) = alpha*x(i)*y(i);
        if (add) zx += z;
        //cerr<<"ElemMultVV: alpha = "<<alpha<<", add = "<<add<<endl;
        //cerr<<"x = "<<TMV_Text(x)<<"  step "<<x.step()<<"  "<<x<<endl;
        //cerr<<"y = "<<TMV_Text(y)<<"  step "<<y.step()<<"  "<<y<<endl;
        //cerr<<"z = "<<TMV_Text(z)<<"  step "<<z.step()<<"  "<<z0<<endl;
        //cerr<<"zx = "<<zx<<endl;
#endif

        if (z.size() > 0 && alpha != T(0)) {
            if (z.isconj()) {
                ElemMultVV<add>(
                    TMV_CONJ(alpha),x.conjugate(),y.conjugate(),z.conjugate());
            } else if (
                (z.step()==-1 && !(x.step()==1 && y.step()==1)) ||
                ( z.step()!=1 && 
                  (x.step()==-1 || (x.step()!=1 && y.step()==-1))) ||
                (z.step()<0 && !(x.step()==1 || y.step()==1))) {
                ElemMultVV<add>(alpha,x.reverse(),y.reverse(),z.reverse());
            } else if (SameStorage(x,z) && x.step() > z.step()) {
                if (add)
                    ElemMultVV<add>(alpha,Vector<Tx>(x),y,z);
                else if (SameStorage(y,z)) {
                    Vector<Tx> xx(x);
                    z = y;
                    ElemMultVV<false>(alpha,xx,z,z);
                } else {
                    z = x;
                    ElemMultVV<false>(alpha,z,y,z);
                }
            } else if (SameStorage(y,z) && y.step() > z.step()) {
                if (add)
                    ElemMultVV<add>(alpha,x,Vector<Ty>(y),z);
                else if (SameStorage(x,z)) {
                    Vector<Ty> yy(y);
                    z = x;
                    ElemMultVV<false>(alpha,z,yy,z);
                } else {
                    z = y;
                    ElemMultVV<false>(alpha,x,z,z);
                }
            } else {
                if (TMV_IMAG(alpha) == RT(0)) {
                    const RT ralpha = TMV_REAL(alpha);
                    if (x.isconj())
                        if (y.isconj())
                            DoElemMultVV<add,true,true>(ralpha,x,y,z);
                        else
                            DoElemMultVV<add,true,false>(ralpha,x,y,z);
                    else
                        if (y.isconj())
                            DoElemMultVV<add,false,true>(ralpha,x,y,z);
                        else
                            DoElemMultVV<add,false,false>(ralpha,x,y,z);
                } else {
                    if (x.isconj())
                        if (y.isconj())
                            DoElemMultVV<add,true,true>(alpha,x,y,z);
                        else
                            DoElemMultVV<add,true,false>(alpha,x,y,z);
                    else
                        if (y.isconj())
                            DoElemMultVV<add,false,true>(alpha,x,y,z);
                        else
                            DoElemMultVV<add,false,false>(alpha,x,y,z);
                }
            }
        }

#ifdef XDEBUG
        if (!(Norm(zx-z) <= 
              0.001*(TMV_ABS(alpha)*(Norm(x)+Norm(y))+(add?Norm(z0):RT(0))))) {
            cerr<<"ElemMultVV: alpha = "<<alpha<<", add = "<<add<<endl;
            cerr<<"x = "<<TMV_Text(x)<<"  step "<<x.step()<<"  "<<x<<endl;
            cerr<<"y = "<<TMV_Text(y)<<"  step "<<y.step()<<"  "<<y<<endl;
            cerr<<"z = "<<TMV_Text(z)<<"  step "<<z.step()<<"  "<<z0<<endl;
            cerr<<"-> "<<z<<endl;
            cerr<<"zx = "<<zx<<endl;
            cerr<<"Norm(zx-z) = "<<Norm(zx-z)<<endl;
            abort();
        }
#endif
    }


#define InstFile "TMV_MultXV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


