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
#include "tmv/TMV_MatrixArithFunc.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"

// CBLAS trick of using RowMajor when we have a 
// case of x.conjugate() ^ y doesn't seem to be working with MKL 10.2.2.
// I haven't been able to figure out why.  (e.g. Is it a bug in the MKL
// code, or am I doing something wrong?)  So for now, just disable it.
#ifdef CBLAS
#undef CBLAS
#endif

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    // 
    // Rank1Update
    //

    template <bool cx, bool cy, bool cm, bool add, class T, class Tx, class Ty> 
    static void ColRank1Update(
        const GenVector<Tx>& x, const GenVector<Ty>& y, MatrixView<T> A)
    {
        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.step() == 1);
        TMVAssert(cm == A.iscm());
        TMVAssert(cx == x.isconj());
        TMVAssert(cy == y.isconj());

        const Ty* yj = y.cptr();
        const Tx*const xptr = x.cptr();
        T* Acolj = A.ptr();
        const ptrdiff_t sj = A.stepj();
        const ptrdiff_t si = (cm ? 1 : A.stepi());
        const ptrdiff_t ys = y.step();
        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();

        for (ptrdiff_t j=N; j>0; --j,yj+=ys,Acolj+=sj) {
            if (*yj!=Ty(0)) {
                T* Aij = Acolj;
                const Tx* xi = xptr;
                for (ptrdiff_t i=M; i>0; --i,++xi,(cm?++Aij:Aij+=si)) {
                    const T temp = 
                        (cx ? TMV_CONJ(*xi) : *xi) *
                        (cy ? TMV_CONJ(*yj) : *yj);
#ifdef TMVFLDEBUG
                    TMVAssert(Aij >= A._first);
                    TMVAssert(Aij < A._last);
#endif
                    if (add) *Aij += temp;
                    else *Aij = temp;
                }
            } else if (!add) {
                T* Aij = Acolj;
                if (cm) std::fill_n(Aij,M,T(0));
                else for (ptrdiff_t i=M; i>0; --i,Aij+=si) {
#ifdef TMVFLDEBUG
                    TMVAssert(Aij >= A._first);
                    TMVAssert(Aij < A._last);
#endif
                    *Aij = T(0);
                }
            }
        }
    }

    template <bool cx, bool add, class T, class Tx, class Ty> 
    static void UnitARank1Update(
        const GenVector<Tx>& x,
        const GenVector<Ty>& y, MatrixView<T> A)
    {
        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.step() == 1);
        TMVAssert(cx == x.isconj());

        if (A.iscm()) 
            if (y.isconj())
                ColRank1Update<cx,true,true,add>(x,y,A);
            else
                ColRank1Update<cx,false,true,add>(x,y,A);
        else
            if (y.isconj())
                ColRank1Update<cx,true,false,add>(x,y,A);
            else
                ColRank1Update<cx,false,false,add>(x,y,A);
    }

    template <bool add, class T, class Tx, class Ty> 
    static void NonBlasRank1Update(
        const T alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, MatrixView<T> A)
    {
        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());
        TMVAssert(alpha != T(0));
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(A.ct() == NonConj);

        if (x.step() != 1 || alpha != T(1)) {
            if (x.step() == 1 && y.size() < x.size()) {
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                    Vector<Ty> yy = TMV_REAL(alpha)*y;
                    if (x.isconj()) UnitARank1Update<true,add>(x,yy,A);
                    else UnitARank1Update<false,add>(x,yy,A);
                } else {
                    Vector<T> yy = alpha*y;
                    if (x.isconj()) UnitARank1Update<true,add>(x,yy,A);
                    else UnitARank1Update<false,add>(x,yy,A);
                }
            } else {
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                    Vector<Tx> xx = TMV_REAL(alpha)*x;
                    UnitARank1Update<false,add>(xx,y,A);
                } else {
                    Vector<T> xx = alpha*x;
                    UnitARank1Update<false,add>(xx,y,A);
                }
            }
        } else {
            if (x.isconj()) UnitARank1Update<true,add>(x,y,A);
            else UnitARank1Update<false,add>(x,y,A);
        }
    }

#ifdef BLAS
    template <class T, class Tx, class Ty> 
    static inline void BlasRank1Update(
        const T alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, MatrixView<T> A)
    { NonBlasRank1Update<true>(alpha,x,y,A); }
#ifdef INST_DOUBLE
    template <> 
    void BlasRank1Update(
        const double alpha, const GenVector<double>& x,
        const GenVector<double>& y, MatrixView<double> A)
    {
        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());
        TMVAssert(alpha != 0.);
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(y.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.iscm());

        int m = A.colsize();
        int n = A.rowsize();
        int xs = x.step();
        int ys = y.step();
        const double* xp = x.cptr();
        if (xs < 0) xp += (m-1)*xs;
        const double* yp = y.cptr();
        if (ys < 0) yp += (n-1)*ys;
        int lda = A.stepj();
        if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
        if (lda < m) { TMVAssert(n == 1); lda = m; }
        BLASNAME(dger) (
            BLASCM BLASV(m),BLASV(n),BLASV(alpha),
            BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
            BLASP(A.ptr()),BLASV(lda));
    }
    template <> 
    void BlasRank1Update(
        const std::complex<double> alpha,
        const GenVector<std::complex<double> >& x, 
        const GenVector<std::complex<double> >& y,
        MatrixView<std::complex<double> > A)
    {
        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());
        TMVAssert(alpha != 0.);
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj || y.ct() == NonConj);
        TMVAssert(A.iscm());

        int m = A.colsize();
        int n = A.rowsize();
        int xs = x.step();
        int ys = y.step();
        const std::complex<double>* xp = x.cptr();
        if (xs < 0) xp += (m-1)*xs;
        const std::complex<double>* yp = y.cptr();
        if (ys < 0) yp += (n-1)*ys;
        int lda = A.stepj();
        if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
        if (lda < m) { TMVAssert(n == 1); lda = m; }
        if (x.isconj()) {
#ifdef CBLAS
            BLASNAME(zgerc) (
                BLASRM BLASV(n),BLASV(m),BLASP(&alpha),
                BLASP(yp),BLASV(ys),BLASP(xp),BLASV(xs),
                BLASP(A.ptr()),BLASV(lda));
#else
            Vector<std::complex<double> > xx = alpha*x;
            xs = 1;
            xp = xx.cptr();
            std::complex<double> alpha2(1);
            BLASNAME(zgeru) (
                BLASCM BLASV(m),BLASV(n),BLASP(&alpha2),
                BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
                BLASP(A.ptr()),BLASV(lda));
#endif
        }
        else if (y.isconj())
            BLASNAME(zgerc) (
                BLASCM BLASV(m),BLASV(n),BLASP(&alpha),
                BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
                BLASP(A.ptr()),BLASV(lda));
        else
            BLASNAME(zgeru) (
                BLASCM BLASV(m),BLASV(n),BLASP(&alpha),
                BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
                BLASP(A.ptr()),BLASV(lda));
    }
    template <> 
    void BlasRank1Update(
        const std::complex<double> alpha,
        const GenVector<double>& x, 
        const GenVector<std::complex<double> >& y,
        MatrixView<std::complex<double> > A)
    {
        // A += a * x ^ y
        // (Ar + I Ai) += (ar + I ai) * x ^ (yr + I yi)
        // Ar += ar * x ^ yr - ai * x ^ yi
        // Ai += ai * x ^ yr + ar * x ^ yi
        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());
        TMVAssert(alpha != 0.);
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.iscm());

        int m = 2*A.colsize();
        int n = A.rowsize();
        int xs = 1;
        int ys = 2*y.step();
        const double* yp = (const double*) y.cptr();
        if (ys < 0) yp += (n-1)*ys;
        int lda = 2*A.stepj();
        if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
        if (lda < m) { TMVAssert(n == 1); lda = m; }
        Vector<double> xx(2*x.size());
        xx.subVector(0,xx.size(),2) = TMV_REAL(alpha)*x;
        xx.subVector(1,xx.size()+1,2) = TMV_IMAG(alpha)*x;
        const double* xp = xx.cptr();
        double xalpha(1);
        BLASNAME(dger) (
            BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
            BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
            BLASP((double*)A.ptr()),BLASV(lda));
        if (y.isconj()) {
            xx.subVector(0,xx.size(),2) = TMV_IMAG(alpha)*x;
            xx.subVector(1,xx.size()+1,2) = -TMV_REAL(alpha)*x;
        } else {
            xx.subVector(0,xx.size(),2) = -TMV_IMAG(alpha)*x;
            xx.subVector(1,xx.size()+1,2) = TMV_REAL(alpha)*x;
        }
        BLASNAME(dger) (
            BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
            BLASP(xp),BLASV(xs),BLASP(yp+1),BLASV(ys),
            BLASP((double*)A.ptr()),BLASV(lda));
    }
    template <> 
    void BlasRank1Update(
        const std::complex<double> alpha,
        const GenVector<std::complex<double> >& x,
        const GenVector<double>& y, 
        MatrixView<std::complex<double> > A)
    {
        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());
        TMVAssert(alpha != 0.);
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(y.ct() == NonConj);
        TMVAssert(A.iscm());

        int m = 2*A.colsize();
        int n = A.rowsize();
        int xs = 1;
        int ys = y.step();
        const double* yp = y.cptr();
        if (ys < 0) yp += (n-1)*ys;
        int lda = 2*A.stepj();
        if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
        if (lda < m) { TMVAssert(n == 1); lda = m; }
        if (x.step() == 1 && !x.isconj() && TMV_IMAG(alpha) == 0.) {
            const double* xp = (double*) x.cptr();
            double xalpha(TMV_REAL(alpha));
            BLASNAME(dger) (
                BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
                BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
                BLASP((double*)A.ptr()),BLASV(lda));
        } else {
            Vector<std::complex<double> > xx = alpha*x;
            const double* xp = (double*) xx.cptr();
            double xalpha(1);
            BLASNAME(dger) (
                BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
                BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
                BLASP((double*)A.ptr()),BLASV(lda));
        } 
    }
    template <> 
    void BlasRank1Update(
        const std::complex<double> alpha,
        const GenVector<double>& x, const GenVector<double>& y, 
        MatrixView<std::complex<double> > A)
    {
        // A += a * x ^ y
        // (Ar + I Ai) += (ar + I ai) * x ^ y
        // Ar += ar * x ^ y
        // Ai += ai * x ^ y
        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());
        TMVAssert(alpha != 0.);
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(y.ct() == NonConj);
        TMVAssert(A.iscm());

        int m = 2*A.colsize();
        int n = A.rowsize();
        int xs = 1;
        int ys = y.step();
        const double* yp = y.cptr();
        if (ys < 0) yp += (n-1)*ys;
        int lda = 2*A.stepj();
        if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
        if (lda < m) { TMVAssert(n == 1); lda = m; }
        Vector<double> xx(2*x.size());
        xx.subVector(0,xx.size(),2) = TMV_REAL(alpha)*x;
        xx.subVector(1,xx.size()+1,2) = TMV_IMAG(alpha)*x;
        const double* xp = xx.cptr();
        double xalpha(1);
        BLASNAME(dger) (
            BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
            BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
            BLASP((double*)A.ptr()),BLASV(lda));
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void BlasRank1Update(
        const float alpha, const GenVector<float>& x,
        const GenVector<float>& y, MatrixView<float> A)
    {
        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());
        TMVAssert(alpha != 0.F);
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(y.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.iscm());

        int m = A.colsize();
        int n = A.rowsize();
        int xs = x.step();
        int ys = y.step();
        const float* xp = x.cptr();
        if (xs < 0) xp += (m-1)*xs;
        const float* yp = y.cptr();
        if (ys < 0) yp += (n-1)*ys;
        int lda = A.stepj();
        if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
        if (lda < m) { TMVAssert(n == 1); lda = m; }
        BLASNAME(sger) (
            BLASCM BLASV(m),BLASV(n),BLASV(alpha),
            BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
            BLASP(A.ptr()),BLASV(lda));
    }
    template <> 
    void BlasRank1Update(
        const std::complex<float> alpha,
        const GenVector<std::complex<float> >& x, 
        const GenVector<std::complex<float> >& y,
        MatrixView<std::complex<float> > A)
    {
        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());
        TMVAssert(alpha != 0.F);
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj || y.ct() == NonConj);
        TMVAssert(A.iscm());

        int m = A.colsize();
        int n = A.rowsize();
        int xs = x.step();
        int ys = y.step();
        const std::complex<float>* xp = x.cptr();
        if (xs < 0) xp += (m-1)*xs;
        const std::complex<float>* yp = y.cptr();
        if (ys < 0) yp += (n-1)*ys;
        int lda = A.stepj();
        if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
        if (lda < m) { TMVAssert(n == 1); lda = m; }
        if (x.isconj()) {
#ifdef CBLAS
            BLASNAME(cgerc) (
                BLASRM BLASV(n),BLASV(m),BLASP(&alpha),
                BLASP(yp),BLASV(ys),BLASP(xp),BLASV(xs),
                BLASP(A.ptr()),BLASV(lda));
#else
            Vector<std::complex<float> > xx = alpha*x;
            xs = 1;
            xp = xx.cptr();
            std::complex<float> alpha2(1);
            BLASNAME(cgeru) (
                BLASCM BLASV(m),BLASV(n),BLASP(&alpha2),
                BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
                BLASP(A.ptr()),BLASV(lda));
#endif
        }
        else if (y.isconj())
            BLASNAME(cgerc) (
                BLASCM BLASV(m),BLASV(n),BLASP(&alpha),
                BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
                BLASP(A.ptr()),BLASV(lda));
        else
            BLASNAME(cgeru) (
                BLASCM BLASV(m),BLASV(n),BLASP(&alpha),
                BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
                BLASP(A.ptr()),BLASV(lda));
    }
    template <> 
    void BlasRank1Update(
        const std::complex<float> alpha,
        const GenVector<float>& x, 
        const GenVector<std::complex<float> >& y,
        MatrixView<std::complex<float> > A)
    {
        // A += a * x ^ y
        // (Ar + I Ai) += (ar + I ai) * x ^ (yr + I yi)
        // Ar += ar * x ^ yr - ai * x ^ yi
        // Ai += ai * x ^ yr + ar * x ^ yi
        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());
        TMVAssert(alpha != 0.F);
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.iscm());

        int m = 2*A.colsize();
        int n = A.rowsize();
        int xs = 1;
        int ys = 2*y.step();
        const float* yp = (float*) y.cptr();
        if (ys < 0) yp += (n-1)*ys;
        int lda = 2*A.stepj();
        if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
        if (lda < m) { TMVAssert(n == 1); lda = m; }
        Vector<float> xx(2*x.size());
        xx.subVector(0,xx.size(),2) = TMV_REAL(alpha)*x;
        xx.subVector(1,xx.size()+1,2) = TMV_IMAG(alpha)*x;
        const float* xp = xx.cptr();
        float xalpha(1);
        BLASNAME(sger) (
            BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
            BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
            BLASP((float*)A.ptr()),BLASV(lda));
        if (y.isconj()) {
            xx.subVector(0,xx.size(),2) = TMV_IMAG(alpha)*x;
            xx.subVector(1,xx.size()+1,2) = -TMV_REAL(alpha)*x;
        } else {
            xx.subVector(0,xx.size(),2) = -TMV_IMAG(alpha)*x;
            xx.subVector(1,xx.size()+1,2) = TMV_REAL(alpha)*x;
        }
        BLASNAME(sger) (
            BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
            BLASP(xp),BLASV(xs),BLASP(yp+1),BLASV(ys),
            BLASP((float*)A.ptr()),BLASV(lda));
    }
    template <> 
    void BlasRank1Update(
        const std::complex<float> alpha,
        const GenVector<std::complex<float> >& x,
        const GenVector<float>& y, 
        MatrixView<std::complex<float> > A)
    {
        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());
        TMVAssert(alpha != 0.F);
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(y.ct() == NonConj);
        TMVAssert(A.iscm());

        int m = 2*A.colsize();
        int n = A.rowsize();
        int xs = 1;
        int ys = y.step();
        const float* yp = y.cptr();
        if (ys < 0) yp += (n-1)*ys;
        int lda = 2*A.stepj();
        if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
        if (lda < m) { TMVAssert(n == 1); lda = m; }
        if (x.step() == 1 && !x.isconj() && TMV_IMAG(alpha) == 0.F) {
            const float* xp = (float*) x.cptr();
            float xalpha(TMV_REAL(alpha));
            BLASNAME(sger) (
                BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
                BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
                BLASP((float*)A.ptr()),BLASV(lda));
        } else {
            Vector<std::complex<float> > xx = alpha*x;
            const float* xp = (float*) xx.cptr();
            float xalpha(1);
            BLASNAME(sger) (
                BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
                BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
                BLASP((float*)A.ptr()),BLASV(lda));
        } 
    }
    template <> 
    void BlasRank1Update(
        const std::complex<float> alpha,
        const GenVector<float>& x, const GenVector<float>& y, 
        MatrixView<std::complex<float> > A)
    {
        // A += a * x ^ y
        // (Ar + I Ai) += (ar + I ai) * x ^ y
        // Ar += ar * x ^ y
        // Ai += ai * x ^ y
        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());
        TMVAssert(alpha != 0.F);
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(y.ct() == NonConj);
        TMVAssert(A.iscm());

        int m = 2*A.colsize();
        int n = A.rowsize();
        int xs = 1;
        int ys = y.step();
        const float* yp = y.cptr();
        if (ys < 0) yp += (n-1)*ys;
        int lda = 2*A.stepj();
        if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
        if (lda < m) { TMVAssert(n == 1); lda = m; }
        Vector<float> xx(2*x.size());
        xx.subVector(0,xx.size(),2) = TMV_REAL(alpha)*x;
        xx.subVector(1,xx.size()+1,2) = TMV_IMAG(alpha)*x;
        const float* xp = xx.cptr();
        float xalpha(1);
        BLASNAME(sger) (
            BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
            BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
            BLASP((float*)A.ptr()),BLASV(lda));
    }
#endif
#endif // BLAS

    template <bool add, class T, class Tx, class Ty> 
    void Rank1Update(
        const T alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, MatrixView<T> A)
        // A (+)= alpha * x * yT
    {
#ifdef XDEBUG
        //cout<<"Rank1Update: alpha = "<<alpha<<endl;
        //cout<<"add = "<<add<<endl;
        //cout<<"x = "<<TMV_Text(x)<<"  "<<x<<endl;
        //cout<<"y = "<<TMV_Text(y)<<"  "<<y<<endl;
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
        Vector<Tx> x0 = x;
        Vector<Ty> y0 = y;
        Matrix<T> A0 = A;
        Matrix<T> A2 = A;
        for(ptrdiff_t i=0;i<x.size();i++) for(ptrdiff_t j=0;j<y.size();j++) 
            if (add)
                A2(i,j) += alpha*x0(i)*y0(j);
            else
                A2(i,j) = alpha*x0(i)*y0(j);
#endif

        TMVAssert(A.colsize() == x.size());
        TMVAssert(A.rowsize() == y.size());

        if (A.colsize() > 0 && A.rowsize() > 0) {
            if (alpha == T(0)) {
                if (!add) A.setZero();
            } else if (A.isconj()) 
                Rank1Update<add>(
                    TMV_CONJ(alpha),x.conjugate(),y.conjugate(),
                    A.conjugate());
            else if (!BlasIsCM(A) && BlasIsRM(A)) 
                Rank1Update<add>(alpha,y,x,A.transpose());
#ifdef BLAS
            else if (!BlasIsCM(A)) {
                Matrix<T,ColMajor> A2(A);
                Rank1Update<add>(alpha,x,y,A2.view());
                A = A2;
            } else if (x.step() != 1 || SameStorage(x,A)) {
                // Most BLAS implementations do fine with the x.step() != 1.
                // However, some implementations seem to propagate nan's from
                // the temporary memory they create to do the unit-1 
                // calculation.
                // So to make sure they don't have to make a temporary, I just
                // do it here for them.
                if (y.step() != 1 || SameStorage(y,A)) {
                    if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                        if (x.size() <= y.size()) {
                            Vector<Tx> xx = TMV_REAL(alpha)*x;
                            Vector<Ty> yy = y;
                            if (!add) A.setZero();
                            BlasRank1Update(T(1),xx,yy,A);
                        } else {
                            Vector<Tx> xx = x;
                            Vector<Ty> yy = TMV_REAL(alpha)*y;
                            if (!add) A.setZero();
                            BlasRank1Update(T(1),xx,yy,A);
                        }
                    } else {
                        if (x.size() <= y.size()) {
                            Vector<T> xx = alpha*x;
                            Vector<Ty> yy = y;
                            if (!add) A.setZero();
                            BlasRank1Update(T(1),xx,yy,A);
                        } else {
                            Vector<Tx> xx = x;
                            Vector<T> yy = alpha*y;
                            if (!add) A.setZero();
                            BlasRank1Update(T(1),xx,yy,A);
                        }
                    }
                } else {
                    if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                        Vector<Tx> xx = TMV_REAL(alpha)*x;
                        if (!add) A.setZero();
                        BlasRank1Update(T(1),xx,y,A);
                    } else {
                        Vector<T> xx = alpha*x;
                        if (!add) A.setZero();
                        BlasRank1Update(T(1),xx,y,A);
                    }
                }
            } else {
                if (y.step() != 1 || SameStorage(A,y)) {
                    if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                        Vector<Ty> yy = TMV_REAL(alpha)*y;
                        if (!add) A.setZero();
                        BlasRank1Update(T(1),x,yy,A);
                    } else {
                        Vector<T> yy = alpha*y;
                        if (!add) A.setZero();
                        BlasRank1Update(T(1),x,yy,A);
                    }
                } else {
                    if (!add) A.setZero();
                    if (x.isconj() && y.isconj()) {
                        if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                            if (x.size() <= y.size()) {
                                Vector<Tx> xx = TMV_REAL(alpha)*x;
                                BlasRank1Update(T(1),xx,y,A);
                            } else {
                                Vector<Ty> yy = TMV_REAL(alpha)*y;
                                BlasRank1Update(T(1),x,yy,A);
                            }
                        } else {
                            if (x.size() <= y.size()) {
                                Vector<T> xx = alpha*x;
                                BlasRank1Update(T(1),xx,y,A);
                            } else {
                                Vector<T> yy = alpha*y;
                                BlasRank1Update(T(1),x,yy,A);
                            }
                        }
                    } else {
                        BlasRank1Update(alpha,x,y,A);
                    }
                }
            }
#else
            else NonBlasRank1Update<add>(alpha,x,y,A);
#endif
        }

#ifdef XDEBUG
        //cout<<"Done Rank1Update: A->"<<A<<endl;
        if (!(Norm(A-A2) <= 0.001*(TMV_ABS(alpha)*Norm(x0)*Norm(y0)+
                                  (add?Norm(A0):TMV_RealType(T)(0))))) {
            cerr<<"Rank1Update: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"x = "<<TMV_Text(x)<<"  step = "<<x.step()<<"  "<<x0<<endl;
            cerr<<"y = "<<TMV_Text(y)<<"  step = "<<y.step()<<"  "<<y0<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"-> A = "<<A<<endl;
            cerr<<"A2 = "<<A2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_Rank1_VVM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


