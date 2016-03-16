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
#include "tmv/TMV_BandMatrixArithFunc.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArithFunc.h"
#ifdef BLAS
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_BandMatrixArith.h"
#endif
#include <iostream>

// CBLAS trick of using RowMajor with ConjTrans when we have a 
// case of A.conjugate() * x doesn't seem to be working with MKL 10.2.2.
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
    // BandMatrixComposite
    //

    template <class T> 
    ptrdiff_t BandMatrixComposite<T>::ls() const 
    {
        return BandStorageLength(
            ColMajor,this->colsize(),this->rowsize(),
            this->nlo(),this->nhi());
    }

    template <class T> 
    ConstVectorView<T> BandMatrixComposite<T>::constLinearView() const 
    {
        cptr(); // This makes the instantiation, but we don't need the result.
        return ConstVectorView<T>(itsm1.get(),ls(),1,NonConj);
    }

    template <class T> 
    const T* BandMatrixComposite<T>::cptr() const
    {
        if (!itsm1.get()) {
            ptrdiff_t cs = this->colsize();
            ptrdiff_t rs = this->rowsize();
            ptrdiff_t lo = this->nlo();
            ptrdiff_t hi = this->nhi();
            ptrdiff_t len = ls();
            itsm1.resize(len);
            ptrdiff_t si = stepi();
            ptrdiff_t sj = stepj();
            ptrdiff_t ds = si + sj;
            itsm = this->isdm() ? itsm1.get()-lo*si : itsm1.get();
            this->assignToB(BandMatrixView<T>(
                    itsm,cs,rs,lo,hi,si,sj,ds,NonConj,len
                    TMV_FIRSTLAST1(itsm1.get(),itsm1.get()+len) ) );
        }
        return itsm;
    }

    template <class T> 
    ptrdiff_t BandMatrixComposite<T>::stepi() const
    { return 1; }

    template <class T> 
    ptrdiff_t BandMatrixComposite<T>::stepj() const
    { return this->nlo()+this->nhi(); }

    template <class T> 
    ptrdiff_t BandMatrixComposite<T>::diagstep() const
    { return this->nlo()+this->nhi()+1; }


    //
    // MultMV
    //

    template <bool add, bool cx, bool ca, bool rm,  class T, class Ta, class Tx>
    static void RowMultMV(
        const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y)
    {
        TMVAssert(A.rowsize()==x.size());
        TMVAssert(A.colsize()==y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(y.ct() == NonConj);
        TMVAssert(x.step()==1);
        TMVAssert(y.step()==1);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(cx == x.isconj());
        TMVAssert(ca == A.isconj());
        TMVAssert(rm == A.isrm());

        const ptrdiff_t si = A.stepi();
        const ptrdiff_t sj = (rm ? 1 : A.stepj());
        const ptrdiff_t ds = A.diagstep();
        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();

        const Ta* Aij1 = A.cptr();
        const Tx* xj1 = x.cptr();
        T* yi = y.ptr();

        ptrdiff_t k=A.nlo();
        ptrdiff_t i=0;
        ptrdiff_t j1=0;
        ptrdiff_t j2=A.nhi()+1;
        ptrdiff_t len = j2; // len = j2-j1
        for(; i<M; ++i, ++yi) {
#ifdef TMVFLDEBUG
            TMVAssert(yi >= y._first);
            TMVAssert(yi < y._last);
#endif
            if (!add) *yi = T(0);

            // *yi += A.row(i,j1,j2) * x.subVector(j1,j2);
            const Ta* Aij = Aij1;
            const Tx* xj = xj1;
            for(ptrdiff_t j=len;j>0;--j,++xj,(rm?++Aij:Aij+=sj)) {
#ifdef TMVFLDEBUG
                TMVAssert(yi >= y._first);
                TMVAssert(yi < y._last);
#endif
                *yi += (cx ? TMV_CONJ(*xj) : *xj) * (ca ? TMV_CONJ(*Aij) : *Aij);
            }

            if (k>0) { --k; ++len; Aij1+=si; }
            else { ++j1; ++xj1; Aij1+=ds; }
            if (j2<N) ++j2;
            else { --len; if (j1==N) { ++i, ++yi; break; } }
        }
        if (!add) for(;i<M; ++i, ++yi) {
#ifdef TMVFLDEBUG
            TMVAssert(yi >= y._first);
            TMVAssert(yi < y._last);
#endif
            *yi = T(0);
        }
    }

    template <bool add, bool cx, bool ca, bool cm, class T, class Ta, class Tx>
    static void ColMultMV(
        const GenBandMatrix<Ta>& A, const GenVector<Tx>& x, VectorView<T> y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(y.ct() == NonConj);
        TMVAssert(x.step()==1);
        TMVAssert(y.step()==1);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(cx == x.isconj());
        TMVAssert(ca == A.isconj());
        TMVAssert(cm == A.iscm());

        const ptrdiff_t N = A.rowsize();
        const ptrdiff_t M = A.colsize();
        const ptrdiff_t si = (cm ? 1 : A.stepi());
        const ptrdiff_t sj = A.stepj();
        const ptrdiff_t ds = A.diagstep();

        const Ta* Ai1j = A.cptr();
        const Tx* xj = x.cptr();
        T* yi1 = y.ptr();

        ptrdiff_t k=A.nhi();
        ptrdiff_t i1=0;
        ptrdiff_t i2=A.nlo()+1;
        ptrdiff_t len = i2; // = i2-i1

        if (!add) y.setZero();

        for(ptrdiff_t j=N; j>0; --j,++xj) {
            if (*xj != Tx(0)) {
                // y.subVector(i1,i2) += *xj * A.col(j,i1,i2);
                const Ta* Aij = Ai1j;
                T* yi = yi1;
                for(ptrdiff_t i=len;i>0;--i,++yi,(cm?++Aij:Aij+=si)) {
#ifdef TMVFLDEBUG
                    TMVAssert(yi >= y._first);
                    TMVAssert(yi < y._last);
#endif
                    *yi += 
                        (cx ? TMV_CONJ(*xj) : *xj) *
                        (ca ? TMV_CONJ(*Aij) : *Aij);
                }
            }
            if (k>0) { --k; Ai1j+=sj; ++len; }
            else { ++i1; ++yi1; Ai1j+=ds; }
            if (i2<M) ++i2; 
            else { --len; if (i1==M) break; }
        }
    }

    template <bool add, bool cx, bool ca, bool dm, class T, class Ta, class Tx>
    static void DiagMultMV(
        const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(y.ct() == NonConj);
        TMVAssert(x.step()==1);
        TMVAssert(y.step()==1);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(cx == x.isconj());
        TMVAssert(ca == A.isconj());
        TMVAssert(dm == A.isdm());

        const ptrdiff_t si = A.stepi();
        const ptrdiff_t sj = A.stepj();
        const ptrdiff_t ds = A.diagstep();
        const ptrdiff_t lo = A.nlo();
        const ptrdiff_t hi = A.nhi();
        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();

        const Ta* Ai1j1 = A.cptr() + lo*si;
        const Tx* xj1 = x.cptr();
        T* yi1 = y.ptr() + lo;

        ptrdiff_t j2=TMV_MIN(M-lo,N);
        ptrdiff_t len=j2; // == j2-j1

        if (!add) y.setZero();

        for(ptrdiff_t k=-A.nlo(); k<=hi; ++k) {
            // y.subVector(i1,i2) += DiagMatrixViewOf(A.diag(k)) * 
            //     x.subVector(j1,j2);
            const Ta* Aij = Ai1j1;
            const Tx* xj = xj1;
            T* yi = yi1;
            for(ptrdiff_t i=len;i>0;--i,++yi,++xj,(dm?++Aij:Aij+=ds)) {
#ifdef TMVFLDEBUG
                TMVAssert(yi >= y._first);
                TMVAssert(yi < y._last);
#endif
                *yi += 
                    (cx ? TMV_CONJ(*xj) : *xj) *
                    (ca ? TMV_CONJ(*Aij) : *Aij);
            }
            if (k<0) { --yi1; ++len; Ai1j1-=si; } 
            else { ++xj1; Ai1j1+=sj; }
            if (j2 < N) ++j2; else --len; 
        }
    }

    template <bool add, bool cx, class T, class Ta, class Tx> 
    static void UnitAMultMV1(
        const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(y.ct() == NonConj);
        TMVAssert(x.step()==1);
        TMVAssert(y.step()==1);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(cx == x.isconj());

        if (A.isrm())
            if (A.isconj())
                RowMultMV<add,cx,true,true>(A,x,y);
            else
                RowMultMV<add,cx,false,true>(A,x,y);
        else if (A.iscm())
            if (A.isconj())
                ColMultMV<add,cx,true,true>(A,x,y);
            else
                ColMultMV<add,cx,false,true>(A,x,y);
        else if (A.isdm())
            if (A.isconj())
                DiagMultMV<add,cx,true,true>(A,x,y);
            else
                DiagMultMV<add,cx,false,true>(A,x,y);
        else
            if (A.isconj())
                DiagMultMV<add,cx,true,false>(A,x,y);
            else
                DiagMultMV<add,cx,false,false>(A,x,y);
    }

    template <bool add, bool cx, class T, class Ta, class Tx> 
    static void UnitAMultMV(
        const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y)
    {
        // Check for 0's in beginning or end of x:
        // y += [ A1 A2 A3 ] [ 0 ]  -->  y += A2 x
        //                   [ x ]
        //                   [ 0 ]

        const ptrdiff_t N = x.size(); // = A.rowsize()
        ptrdiff_t j2 = N;
        for(const Tx* x2=x.cptr()+N-1; j2>0 && *x2==Tx(0); --j2,--x2);
        if (j2 == 0) {
            if (!add) y.setZero();
            return;
        }
        ptrdiff_t j1 = 0;
        for(const Tx* x1=x.cptr(); *x1==Tx(0); ++j1,++x1);
        if (j1 == 0 && j2 == N) UnitAMultMV1<add,cx>(A,x,y);
        else {
            const ptrdiff_t hi = A.nhi();
            const ptrdiff_t lo = A.nlo();
            const ptrdiff_t M = y.size(); // = A.colsize()
            // This next bit is copied from the BandMatrix colRange function
            ptrdiff_t i1 = j1 > hi ? j1-hi : 0;
            ptrdiff_t i2 = TMV_MIN(j2+lo,M);
            ptrdiff_t newhi = j1 < hi ? hi-j1 : 0;
            ptrdiff_t newlo = lo+hi-newhi;
            ptrdiff_t newM = i2-i1;
            ptrdiff_t newN = j2-j1;
            TMVAssert(newM > 0);
            TMVAssert(newN > 0);
            if (newhi >= newN) newhi = newN-1;
            if (newlo >= newM) newlo = newM-1;
            TMVAssert(A.hasSubBandMatrix(i1,i2,j1,j2,newlo,newhi,1,1));
            const Ta* p = A.cptr()+i1*A.stepi()+j1*A.stepj();
            ConstBandMatrixView<Ta> Acols(
                p,newM,newN,newlo,newhi,
                A.stepi(),A.stepj(),A.diagstep(),A.ct());
            UnitAMultMV1<add,cx>(Acols,x.subVector(j1,j2),y.subVector(i1,i2));
            if (!add) {
                y.subVector(0,i1).setZero();
                y.subVector(i2,M).setZero();
            }
        }
    }

    template <bool add, class T, class Ta, class Tx> 
    static void NonBlasMultMV(
        const T alpha, const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y)
    // y (+)= alpha * A * x
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(alpha != T(0));
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(y.ct() == NonConj);

#ifdef XDEBUG
        cout<<"NonBlasMultMV: A = "<<A<<endl;
        Vector<T> y0 = y;
        Vector<Tx> x0 = x;
        Matrix<Ta> A0 = A;
        Vector<T> y2 = alpha*A0*x0;
        if (add) y2 += y0;
#endif

        if (x.step() != 1 || SameStorage(x,y)) {
            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                Vector<Tx> xx = TMV_REAL(alpha) * x;
                if (y.step() != 1) {
                    Vector<T> yy(y.size());
                    UnitAMultMV<false,false>(A,xx,yy.view());
                    if (add) y += yy;
                    else y = yy;
                } else {
                    UnitAMultMV<add,false>(A,xx,y);
                }
            } else {
                Vector<T> xx = alpha * x;
                if (y.step() != 1) {
                    Vector<T> yy(y.size());
                    UnitAMultMV<false,false>(A,xx,yy.view());
                    if (add) y += yy;
                    else y = yy;
                } else {
                    UnitAMultMV<add,false>(A,xx,y);
                }
            }
        } else if (y.step() != 1 || alpha != TMV_RealType(T)(1)) {
            Vector<T> yy(y.size());
            if (x.isconj())
                UnitAMultMV<false,true>(A,x,yy.view());
            else
                UnitAMultMV<false,false>(A,x,yy.view());
            if (add) y += alpha * yy;
            else y = alpha * yy;
        } else {
            TMVAssert(alpha == T(1));
            TMVAssert(y.step() == 1);
            TMVAssert(x.step() == 1);
            TMVAssert(!SameStorage(x,y));
            if (x.isconj())
                UnitAMultMV<add,true>(A,x,y);
            else
                UnitAMultMV<add,false>(A,x,y);
        }

#ifdef XDEBUG
        if (!(Norm(y2-y) <= 0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(x0)+
                                (add?Norm(y0):TMV_RealType(T)(0))))) {
            cerr<<"NonBlas MultMV: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"  "<<A0<<endl;
            cerr<<"x = "<<TMV_Text(x)<<"  "<<x.cptr()<<
                " step "<<x.step()<<"  "<<x0<<endl;
            cerr<<"y = "<<TMV_Text(y)<<"  "<<y.cptr()<<
                " step "<<y.step()<<"  "<<y0<<endl;
            cerr<<"--> y = "<<y<<endl;
            cerr<<"y2 = "<<y2<<endl;
            abort();
        }
#endif
    }

#ifdef BLAS
    template <class T, class Ta, class Tx> 
    static inline void BlasMultMV(
        const T alpha, const GenBandMatrix<Ta>& A,
        const GenVector<Tx>& x, int beta, VectorView<T> y)
    { 
        if (beta == 1) NonBlasMultMV<true>(alpha,A,x,y); 
        else NonBlasMultMV<false>(alpha,A,x,y); 
    }
#ifdef INST_DOUBLE
    template <> 
    void BlasMultMV(
        const double alpha,
        const GenBandMatrix<double>& A, const GenVector<double>& x,
        int beta, VectorView<double> y)
    {
        int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
        int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
        int lo = BlasIsCM(A) ? A.nlo() : A.nhi();
        int hi = BlasIsCM(A) ? A.nhi() : A.nlo();
        int ds = A.diagstep();
        int xs = x.step();
        int ys = y.step();
        const double* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        double* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        if (beta == 0) y.setZero();
        double xbeta(1);
        BLASNAME(dgbmv) (
            BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
            BLASV(alpha),BLASP(A.cptr()-hi),BLASV(ds),
            BLASP(xp),BLASV(xs),BLASV(xbeta),
            BLASP(yp),BLASV(ys) BLAS1); 
    }
    template <> 
    void BlasMultMV(
        const std::complex<double> alpha,
        const GenBandMatrix<std::complex<double> >& A,
        const GenVector<std::complex<double> >& x,
        int beta, VectorView<std::complex<double> > y)
    {
        if (x.isconj()
#ifndef CBLAS
            && !(A.isconj() && BlasIsCM(A))
#endif
        ) {
            Vector<std::complex<double> > xx = alpha*x;
            return BlasMultMV(std::complex<double>(1),A,xx,beta,y);
        }

        int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
        int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
        int lo = BlasIsCM(A) ? A.nlo() : A.nhi();
        int hi = BlasIsCM(A) ? A.nhi() : A.nlo();
        int ds = A.diagstep();
        int xs = x.step();
        int ys = y.step();
        const std::complex<double>* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        std::complex<double>* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        if (beta == 0) y.setZero();
        std::complex<double> xbeta(1);

        if (A.isconj() && BlasIsCM(A)) {
#ifdef CBLAS
            TMV_SWAP(m,n);
            TMV_SWAP(lo,hi);
            BLASNAME(zgbmv) (
                BLASRM BLASCH_CT, BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASP(&alpha),BLASP(A.cptr()-lo),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),
                BLASP(yp),BLASV(ys) BLAS1); 
#else
            std::complex<double> ca = TMV_CONJ(alpha);
            if (x.isconj()) {
                y.conjugateSelf();
                BLASNAME(zgbmv) (
                    BLASCM BLASCH_NT, BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                    BLASP(&ca),BLASP(A.cptr()-hi),BLASV(ds),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),
                    BLASP(yp),BLASV(ys) BLAS1); 
                y.conjugateSelf();
            } else {
                Vector<std::complex<double> > xx=ca*x.conjugate();
                ca = std::complex<double>(1);
                xs = 1;
                xp = xx.cptr();
                y.conjugateSelf();
                BLASNAME(zgbmv) (
                    BLASCM BLASCH_NT, BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                    BLASP(&ca),BLASP(A.cptr()-hi),BLASV(ds),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),
                    BLASP(yp),BLASV(ys) BLAS1); 
                y.conjugateSelf();
            }
#endif
        } else {
            BLASNAME(zgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASP(&alpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),
                BLASP(yp),BLASV(ys) BLAS1); 
        }
    }
    template <> 
    void BlasMultMV(
        const std::complex<double> alpha,
        const GenBandMatrix<std::complex<double> >& A,
        const GenVector<double>& x,
        int beta, VectorView<std::complex<double> > y)
    { BlasMultMV(alpha,A,Vector<std::complex<double> >(x),beta,y); }
    template <> 
    void BlasMultMV(
        const std::complex<double> alpha,
        const GenBandMatrix<double>& A,
        const GenVector<std::complex<double> >& x,
        int beta, VectorView<std::complex<double> > y)
    {
        if (beta == 0) {
            int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
            int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
            int lo = BlasIsCM(A) ? A.nlo() : A.nhi();
            int hi = BlasIsCM(A) ? A.nhi() : A.nlo();
            int ds = A.diagstep();
            int xs = 2*x.step();
            int ys = 2*y.step();
            const double* xp = (const double*) x.cptr();
            if (xs < 0) xp += (x.size()-1)*xs;
            double* yp = (double*) y.ptr();
            if (ys < 0) yp += (y.size()-1)*ys;
            double xalpha(1);
            y.setZero();
            double xbeta(1);
            BLASNAME(dgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(xalpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1); 
            BLASNAME(dgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(xalpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp+1),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1); 
            if (x.isconj()) y.conjugateSelf();
            y *= alpha;
        } else if (TMV_IMAG(alpha) == 0. && !x.isconj()) {
            int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
            int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
            int lo = BlasIsCM(A) ? A.nlo() : A.nhi();
            int hi = BlasIsCM(A) ? A.nhi() : A.nlo();
            int ds = A.diagstep();
            int xs = 2*x.step();
            int ys = 2*y.step();
            const double* xp = (const double*) x.cptr();
            if (xs < 0) xp += (x.size()-1)*xs;
            double* yp = (double*)  y.ptr();
            if (ys < 0) yp += (y.size()-1)*ys;
            double xalpha(TMV_REAL(alpha));
            double xbeta(beta);
            BLASNAME(dgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(xalpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1); 
            BLASNAME(dgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(xalpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp+1),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1); 
        } else {
            Vector<std::complex<double> > xx = alpha*x;
            BlasMultMV(std::complex<double>(1),A,xx,1,y);
        }
    }
    template <> 
    void BlasMultMV(
        const std::complex<double> alpha,
        const GenBandMatrix<double>& A,
        const GenVector<double>& x,
        int beta, VectorView<std::complex<double> > y)
    {
        int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
        int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
        int lo = BlasIsCM(A) ? A.nlo() : A.nhi();
        int hi = BlasIsCM(A) ? A.nhi() : A.nlo();
        int ds = A.diagstep();
        int xs = x.step();
        int ys = 2*y.step();
        const double* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        double* yp = (double*) y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        double ar(TMV_REAL(alpha));
        double ai(TMV_IMAG(alpha));
        if (beta == 0) y.setZero();
        double xbeta(1);
        if (ar != 0.) {
            BLASNAME(dgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(ar),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1); 
        }
        if (ai != 0.) {
            BLASNAME(dgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(ai),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1); 
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void BlasMultMV(
        const float alpha,
        const GenBandMatrix<float>& A, const GenVector<float>& x,
        int beta, VectorView<float> y)
    {
        int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
        int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
        int lo = BlasIsCM(A) ? A.nlo() : A.nhi();
        int hi = BlasIsCM(A) ? A.nhi() : A.nlo();
        int ds = A.diagstep();
        int xs = x.step();
        int ys = y.step();
        const float* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        float* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        if (beta == 0) y.setZero();
        float xbeta(1);
        BLASNAME(sgbmv) (
            BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
            BLASV(alpha),BLASP(A.cptr()-hi),BLASV(ds),
            BLASP(xp),BLASV(xs),BLASV(xbeta),
            BLASP(yp),BLASV(ys) BLAS1); 
    }
    template <> 
    void BlasMultMV(
        const std::complex<float> alpha,
        const GenBandMatrix<std::complex<float> >& A,
        const GenVector<std::complex<float> >& x,
        int beta, VectorView<std::complex<float> > y)
    {
        if (x.isconj()
#ifndef CBLAS
            && !(A.isconj() && BlasIsCM(A))
#endif
        ) {
            Vector<std::complex<float> > xx = alpha*x;
            return BlasMultMV(std::complex<float>(1),A,xx,beta,y);
        }

        int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
        int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
        int lo = BlasIsCM(A) ? A.nlo() : A.nhi();
        int hi = BlasIsCM(A) ? A.nhi() : A.nlo();
        int ds = A.diagstep();
        int xs = x.step();
        int ys = y.step();
        const std::complex<float>* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        std::complex<float>* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        if (beta == 0) y.setZero();
        std::complex<float> xbeta(1);

        if (A.isconj() && BlasIsCM(A)) {
#ifdef CBLAS
            TMV_SWAP(m,n);
            TMV_SWAP(lo,hi);
            BLASNAME(cgbmv) (
                BLASRM BLASCH_CT, BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASP(&alpha),BLASP(A.cptr()-lo),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),
                BLASP(yp),BLASV(ys) BLAS1); 
#else
            std::complex<float> ca = TMV_CONJ(alpha);
            if (x.isconj()) {
                y.conjugateSelf();
                BLASNAME(cgbmv) (
                    BLASCM BLASCH_NT, BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                    BLASP(&ca),BLASP(A.cptr()-hi),BLASV(ds),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),
                    BLASP(yp),BLASV(ys) BLAS1); 
                y.conjugateSelf();
            } else {
                Vector<std::complex<float> > xx=ca*x.conjugate();
                ca = std::complex<float>(1);
                xs = 1;
                xp = xx.cptr();
                y.conjugateSelf();
                BLASNAME(cgbmv) (
                    BLASCM BLASCH_NT, BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                    BLASP(&ca),BLASP(A.cptr()-hi),BLASV(ds),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),
                    BLASP(yp),BLASV(ys) BLAS1); 
                y.conjugateSelf();
            }
#endif
        } else {
            BLASNAME(cgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASP(&alpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),
                BLASP(yp),BLASV(ys) BLAS1); 
        }
    }
    template <> 
    void BlasMultMV(
        const std::complex<float> alpha,
        const GenBandMatrix<std::complex<float> >& A,
        const GenVector<float>& x,
        int beta, VectorView<std::complex<float> > y)
    { BlasMultMV(alpha,A,Vector<std::complex<float> >(x),beta,y); }
    template <> 
    void BlasMultMV(
        const std::complex<float> alpha,
        const GenBandMatrix<float>& A,
        const GenVector<std::complex<float> >& x,
        int beta, VectorView<std::complex<float> > y)
    {
        if (beta == 0) {
            int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
            int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
            int lo = BlasIsCM(A) ? A.nlo() : A.nhi();
            int hi = BlasIsCM(A) ? A.nhi() : A.nlo();
            int ds = A.diagstep();
            int xs = 2*x.step();
            int ys = 2*y.step();
            const float* xp = (const float*) x.cptr();
            if (xs < 0) xp += (x.size()-1)*xs;
            float* yp = (float*) y.ptr();
            if (ys < 0) yp += (y.size()-1)*ys;
            float xalpha(1);
            y.setZero();
            float xbeta(1);
            BLASNAME(sgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(xalpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1); 
            BLASNAME(sgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(xalpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp+1),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1); 
            if (x.isconj()) y.conjugateSelf();
            y *= alpha;
        } else if (TMV_IMAG(alpha) == 0.F && !x.isconj()) {
            int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
            int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
            int lo = BlasIsCM(A) ? A.nlo() : A.nhi();
            int hi = BlasIsCM(A) ? A.nhi() : A.nlo();
            int ds = A.diagstep();
            int xs = 2*x.step();
            int ys = 2*y.step();
            const float* xp = (const float*) x.cptr();
            if (xs < 0) xp += (x.size()-1)*xs;
            float* yp = (float*) y.ptr();
            if (ys < 0) yp += (y.size()-1)*ys;
            float xalpha(TMV_REAL(alpha));
            float xbeta(1);
            BLASNAME(sgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(xalpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1); 
            BLASNAME(sgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(xalpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp+1),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1); 
        } else {
            Vector<std::complex<float> > xx = alpha*x;
            BlasMultMV(std::complex<float>(1),A,xx,1,y);
        }
    }
    template <> 
    void BlasMultMV(
        const std::complex<float> alpha,
        const GenBandMatrix<float>& A,
        const GenVector<float>& x,
        int beta, VectorView<std::complex<float> > y)
    {
        int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
        int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
        int lo = BlasIsCM(A) ? A.nlo() : A.nhi();
        int hi = BlasIsCM(A) ? A.nhi() : A.nlo();
        int ds = A.diagstep();
        int xs = x.step();
        int ys = 2*y.step();
        const float* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        float* yp = (float*) y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        float ar(TMV_REAL(alpha));
        float ai(TMV_IMAG(alpha));
        if (beta == 0) y.setZero();
        float xbeta(1);
        if (ar != 0.F) {
            BLASNAME(sgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(ar),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1); 
        }
        if (ai != 0.F) {
            BLASNAME(sgbmv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(ai),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1); 
        }
    }
#endif
#endif // BLAS

    template <bool add, class T, class Ta, class Tx> 
    static void DoMultMV(
        const T alpha, const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y)
    {
        //cout<<"Start DoMultMV\n";
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"  "<<A<<endl;
        //cout<<"x = "<<TMV_Text(x)<<"  "<<x.cptr()<<"  step "<<x.step()<<"  "<<x<<endl;
        //cout<<"y = "<<TMV_Text(y)<<"  "<<y.cptr()<<"  step "<<y.step()<<"  "<<y<<endl;
        //cout<<"alpha = "<<alpha<<", add = "<<add<<endl;
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(alpha != T(0));
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);

        if (y.isconj()) {
            DoMultMV<add>(
                TMV_CONJ(alpha),A.conjugate(),x.conjugate(),y.conjugate());
        } else {
#ifdef BLAS
            if (x.step() == 0) {
                if (x.size() <= 1) 
                    DoMultMV<add>(
                        alpha,A,
                        ConstVectorView<Tx>(x.cptr(),x.size(),1,x.ct()),y);
                else 
                    DoMultMV<add>(alpha,A,Vector<Tx>(x),y);
            } else if (y.step() == 0) {
                TMVAssert(y.size() <= 1);
                DoMultMV<add>(
                    alpha,A,x,VectorView<T>(y.ptr(),y.size(),1,y.ct()));
            } else if (BlasIsRM(A) || BlasIsCM(A)) {
                if (!SameStorage(A,y)) {
                    if (!SameStorage(x,y) && !SameStorage(A,x)) {
                        BlasMultMV(alpha,A,x,add?1:0,y);
                    } else {
                        Vector<T> xx = alpha*x;
                        BlasMultMV(T(1),A,xx,add?1:0,y);
                    }
                } else {
                    Vector<T> yy(y.size(),T(0));
                    if (!SameStorage(A,x)) {
                        BlasMultMV(T(1),A,x,0,yy.view());
                        if (add) y += alpha*yy;
                        else y = alpha*yy;
                    } else {
                        Vector<T> xx = alpha*x;
                        BlasMultMV(T(1),A,xx,0,yy.view());
                        if (add) y += yy;
                        else y = yy;
                    }
                }
            } else if ((A.isrm() && A.stepi() < A.nlo()+A.nhi()) ||
                       (A.iscm() && A.stepj() < A.nlo()+A.nhi())) {
                if (SameStorage(A,y)) {
                    Vector<T> yy(y.size(),T(0));
                    DoMultMV<false>(T(1),A,x,yy.view());
                    if (add) y += alpha*yy;
                    else y = alpha*yy;
                } else if (SameStorage(x,y)) {
                    DoMultMV<add>(T(1),A,Vector<T>(alpha*x),y);
                } else if (A.nlo()+1 == A.colsize()) {
                    if (A.nhi()+1 == A.rowsize()) {
                        ConstMatrixView<Ta> A1 =
                            A.subMatrix(0,A.colsize(),0,A.rowsize());
                        MultMV<add>(alpha,A1,x,y);
                    } else {
                        ConstMatrixView<Ta> A1 =
                            A.subMatrix(0,A.colsize(),0,A.nhi());
                        MultMV<add>(alpha,A1,x.subVector(0,A.nhi()),y);
                        ConstBandMatrixView<Ta> A2 = A.colRange(A.nhi(),A.rowsize());
                        BlasMultMV(alpha,A2,x.subVector(A.nhi(),A.rowsize()),1,y);
                    }
                } else {
                    TMVAssert(A.nlo()>0);
                    if (A.nhi()+1 == A.rowsize()) {
                        ConstMatrixView<Ta> A1 =
                            A.subMatrix(0,A.nlo(),0,A.rowsize());
                        MultMV<add>(alpha,A1,x,y.subVector(0,A.nlo()));
                    } else {
                        ConstBandMatrixView<Ta> A1 = A.rowRange(0,A.nlo());
                        BlasMultMV(alpha,A1,x.subVector(0,A1.rowsize()),
                                   add?1:0,y.subVector(0,A.nlo()));
                    }
                    ConstBandMatrixView<Ta> A2 = A.rowRange(A.nlo(),A.colsize());
                    BlasMultMV(alpha,A2,x,add?1:0,y.subVector(A.nlo(),A.colsize()));
                }
            } else {
                if (TMV_IMAG(alpha) == T(0)) {
                    BandMatrix<Ta,RowMajor> A2 = TMV_REAL(alpha)*A;
                    DoMultMV<add>(T(1),A2,x,y);
                } else {
                    BandMatrix<T,RowMajor> A2 = alpha*A;
                    DoMultMV<add>(T(1),A2,x,y);
                }
            }
#else
            NonBlasMultMV<add>(alpha,A,x,y);
#endif
        }
        //std::cout<<"Done DoMultMV\n";
    }

    //
    // MultEqMV
    //

    template <bool rm, bool ca, class T, class Ta> 
    static void DoRowUpperMultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.colsize() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step()==1);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(rm == A.isrm());
        TMVAssert(ca == A.isconj());

        const ptrdiff_t N = x.size();
        const ptrdiff_t sj = (rm ? 1 : A.stepj());
        const ptrdiff_t ds = A.diagstep();

        T* xi = x.ptr();
        const Ta* Aii = A.cptr();
        ptrdiff_t j2=A.nhi()+1;
        ptrdiff_t len = j2-1;

        for(; len>0; ++xi,Aii+=ds) {
            // i = 0..N-2
            // x(i) = A.row(i,i,j2) * x.subVector(i,j2);
#ifdef TMVFLDEBUG
            TMVAssert(xi >= x._first);
            TMVAssert(xi < x._last);
#endif
            *xi *= (ca ? TMV_CONJ(*Aii) : *Aii);
            const T* xj = xi+1;
            const Ta* Aij = Aii+sj;
            for(ptrdiff_t j=len;j>0;--j,++xj,(rm?++Aij:Aij+=sj)) {
#ifdef TMVFLDEBUG
                TMVAssert(xi >= x._first);
                TMVAssert(xi < x._last);
#endif
                *xi += (*xj) * (ca ? TMV_CONJ(*Aij) : *Aij);
            }

            if (j2<N) ++j2;
            else --len;
        }
#ifdef TMVFLDEBUG
        TMVAssert(xi >= x._first);
        TMVAssert(xi < x._last);
#endif
        *xi *= (ca ? TMV_CONJ(*Aii) : *Aii);
    }

    template <bool rm, class T, class Ta> 
    static inline void RowUpperMultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    { 
        if (A.isconj())
            DoRowUpperMultEqMV<rm,true>(A,x);
        else
            DoRowUpperMultEqMV<rm,false>(A,x);
    }

    template <bool cm, bool ca, class T, class Ta> 
    static void DoColUpperMultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.colsize() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step()==1);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(cm == A.iscm());
        TMVAssert(ca == A.isconj());

        const ptrdiff_t N = x.size();
        const ptrdiff_t si = cm ? 1 : A.stepi();
        const ptrdiff_t sj = A.stepj();
        const ptrdiff_t ds = A.diagstep();

        const Ta* Ai1j = A.cptr();
        T* xi1 = x.ptr();

        *xi1 *= (ca ? TMV_CONJ(*Ai1j) : *Ai1j);
        Ai1j += sj;
        const T* xj = x.cptr()+1;

        ptrdiff_t k=A.nhi()-1;
        ptrdiff_t len = 1;
        for(ptrdiff_t j=1; j<N; ++j,++xj) {
            if (*xj != T(0)) {
                // j = 1..N-1
                // x.subVector(i1,j) += x(j) * A.col(j,i1,j);
                const Ta* Aij = Ai1j;
                T* xi = xi1;
                for(ptrdiff_t i=len;i>0;--i,++xi,(cm?++Aij:Aij+=si)) {
#ifdef TMVFLDEBUG
                    TMVAssert(xi >= x._first);
                    TMVAssert(xi < x._last);
#endif
                    *xi += *xj * (ca ? TMV_CONJ(*Aij) : *Aij);
                }
                // Now Aij == Ajj, xi == xj
                // so this next statement is really *xj *= *Ajj
#ifdef TMVFLDEBUG
                TMVAssert(xi >= x._first);
                TMVAssert(xi < x._last);
#endif
                *xi *= (ca ? TMV_CONJ(*Aij) : *Aij);
            }

            if (k>0) { --k; Ai1j+=sj; ++len; }
            else { ++xi1; Ai1j+=ds; }
        }
    }

    template <bool cm, class T, class Ta> 
    static inline void ColUpperMultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    { 
        if (A.isconj())
            DoColUpperMultEqMV<cm,true>(A,x);
        else
            DoColUpperMultEqMV<cm,false>(A,x);
    }

    template <bool rm, bool ca, class T, class Ta> 
    static void DoRowLowerMultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.colsize() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step()==1);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(rm == A.isrm());
        TMVAssert(ca == A.isconj());

        const ptrdiff_t N = x.size();
        const ptrdiff_t si = A.stepi();
        const ptrdiff_t sj = (rm ? 1 : A.stepj());
        const ptrdiff_t ds = A.diagstep();

        ptrdiff_t j1 = N-1-A.nlo();
        const T* xj1 = x.cptr() + j1;
        T* xi = x.ptr() + N-1;
        const Ta* Aii = A.cptr() + (N-1)*ds;
        const Ta* Aij1 = Aii - A.nlo()*sj;
        ptrdiff_t len = A.nlo();

        for(; len>0; --xi,Aii-=ds) {
            // i = N-1..1
            // x(i) = A.row(i,j1,i+1) * x.subVector(j1,i+1);
#ifdef TMVFLDEBUG
            TMVAssert(xi >= x._first);
            TMVAssert(xi < x._last);
#endif
            *xi *= (ca ? TMV_CONJ(*Aii) : *Aii);
            const Ta* Aij = Aij1;
            const T* xj = xj1;
            for(ptrdiff_t j=len;j>0;--j,++xj,(rm?++Aij:Aij+=sj)) {
#ifdef TMVFLDEBUG
                TMVAssert(xi >= x._first);
                TMVAssert(xi < x._last);
#endif
                *xi += *xj * (ca ? TMV_CONJ(*Aij) : *Aij);
            }

            if (j1>0) { --j1; Aij1-=ds; --xj1; }
            else { --len; Aij1-=si; }
        }
#ifdef TMVFLDEBUG
        TMVAssert(xi >= x._first);
        TMVAssert(xi < x._last);
#endif
        *xi *= (ca ? TMV_CONJ(*Aii) : *Aii);
    }

    template <bool rm, class T, class Ta> 
    static inline void RowLowerMultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    { 
        if (A.isconj())
            DoRowLowerMultEqMV<rm,true>(A,x);
        else
            DoRowLowerMultEqMV<rm,false>(A,x);
    }

    template <bool cm, bool ca, class T, class Ta> 
    static void DoColLowerMultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.colsize() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(cm == A.iscm());
        TMVAssert(ca == A.isconj());

        const ptrdiff_t N = x.size();
        const ptrdiff_t si = cm ? 1 : A.stepi();
        const ptrdiff_t ds = A.diagstep();

        T* xj = x.ptr() + N-1;
        const Ta* Ajj = A.cptr()+(N-1)*ds;

#ifdef TMVFLDEBUG
        TMVAssert(xj >= x._first);
        TMVAssert(xj < x._last);
#endif
        *xj *= (ca ? TMV_CONJ(*Ajj) : *Ajj);
        --xj;
        Ajj -= ds;

        ptrdiff_t k=A.nlo()-1;
        for(ptrdiff_t j=N-1,len=1;j>0;--j,--xj,Ajj-=ds) {
            if (*xj!=T(0)) {
                // Actual j = N-2..0
                // x.subVector(j+1,N) += *xj * A.col(j,j+1,N);
                T* xi = xj+1;
                const Ta* Aij = Ajj+si;
                for (ptrdiff_t i=len;i>0;--i,++xi,(cm?++Aij:Aij+=si)) {
#ifdef TMVFLDEBUG
                    TMVAssert(xi >= x._first);
                    TMVAssert(xi < x._last);
#endif
                    *xi += *xj * (ca ? TMV_CONJ(*Aij) : *Aij);
                }
#ifdef TMVFLDEBUG
                TMVAssert(xj >= x._first);
                TMVAssert(xj < x._last);
#endif
                *xj *= (ca ? TMV_CONJ(*Ajj) : *Ajj);

            }
            if (k>0) { --k; ++len; }
        }
    }

    template <bool cm, class T, class Ta> 
    static inline void ColLowerMultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    { 
        if (A.isconj())
            DoColLowerMultEqMV<cm,true>(A,x);
        else
            DoColLowerMultEqMV<cm,false>(A,x);
    }

    template <class T, class Ta> 
    static inline void DoUpperMultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    // x = A * x
    {
        if (A.isrm()) RowUpperMultEqMV<true>(A,x);
        else if (A.iscm()) ColUpperMultEqMV<true>(A,x);
        else RowUpperMultEqMV<false>(A,x);
    }

    template <class T, class Ta> 
    static inline void DoLowerMultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    {
        if (A.isrm()) RowLowerMultEqMV<true>(A,x);
        else if (A.iscm() && !SameStorage(A,x))
            ColLowerMultEqMV<true>(A,x);
        else RowLowerMultEqMV<false>(A,x);
    }

    template <class T, class Ta> 
    static void NonBlasUpperMultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.colsize() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(x.ct() == NonConj);

        //     [ A11 A12  0  ] [ 0  ]   [ A12 x2 ]
        // x = [  0  A22 A23 ] [ x2 ] = [ A22 x2 ]
        //     [  0   0  A33 ] [ 0  ]   [   0    ]

        const ptrdiff_t N = x.size(); // = A.size()
        ptrdiff_t j2 = N;
        for(const T* x2=x.cptr()+N-1; j2>0 && *x2==T(0); --j2,--x2);
        if (j2 == 0) return;
        ptrdiff_t j1 = 0;
        for(const T* x1=x.cptr(); *x1==T(0); ++j1,++x1);
        if (j1 == 0 && j2 == N) DoUpperMultEqMV(A,x);
        else {
            TMVAssert(j1 < j2);
            const Ta* p22 = A.cptr() + j1*A.diagstep();
            const ptrdiff_t N22 = j2-j1;
            VectorView<T> x2 = x.subVector(j1,j2);
            if (N22 > A.nhi()) {
                ConstBandMatrixView<Ta> A22(
                    p22,N22,N22,0,A.nhi(),
                    A.stepi(),A.stepj(),A.diagstep(),A.ct());
                if (j1 > 0) {
                    const ptrdiff_t jx = j1+A.nhi();
                    if (j1 < A.nhi()) {
                        const Ta* p12 = A.cptr() + j1*A.stepj();
                        ConstBandMatrixView<Ta> A12(
                            p12,j1,A.nhi(),j1-1,A.nhi()-j1,
                            A.stepi(),A.stepj(),A.diagstep(),A.ct());
                        UnitAMultMV1<false,false>(
                            A12,x.subVector(j1,jx),x.subVector(0,j1));
                    } else {
                        const Ta* p12 = p22 - A.nhi()*A.stepi();
                        ConstBandMatrixView<Ta> A12(
                            p12,A.nhi(),A.nhi(),A.nhi()-1,0,
                            A.stepi(),A.stepj(),A.diagstep(),A.ct());
                        VectorView<T> x1x = x.subVector(j1-A.nhi(),j1);
                        x1x = x.subVector(j1,jx);
                        DoLowerMultEqMV(A12,x1x);
                    }
                }
                DoUpperMultEqMV(A22,x2);
            } else {
                ConstBandMatrixView<Ta> A22(
                    p22,N22,N22,0,N22-1,
                    A.stepi(),A.stepj(),A.diagstep(),A.ct());
                if (j1 > 0) {
                    const ptrdiff_t M12 = (j1 < A.nhi()) ? j1 : A.nhi();
                    const Ta* p12 = p22 - M12*A.stepi();
                    ptrdiff_t newhi = A.nhi()-M12;
                    if (newhi >= N22) newhi = N22-1;
                    ConstBandMatrixView<Ta> A12(
                        p12,M12,N22,M12-1,newhi,
                        A.stepi(),A.stepj(),A.diagstep(),A.ct());
                    UnitAMultMV1<false,false>(A12,x2,x.subVector(j1-M12,j1));
                } 
                DoUpperMultEqMV(A22,x2);
            }
        }
    }

    template <class T, class Ta> 
    static void NonBlasLowerMultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    // x = A * x
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.colsize() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(x.ct() == NonConj);

        //     [ A11  0   0  ] [ 0  ]   [   0    ]
        // x = [ A21 A22  0  ] [ x2 ] = [ A22 x2 ]
        //     [  0  A32 A33 ] [ 0  ]   [ A32 x2 ]

        const ptrdiff_t N = x.size(); // = A.size()
        ptrdiff_t j2 = N;
        for(const T* x2=x.cptr()+N-1; j2>0 && *x2==T(0); --j2,--x2);
        if (j2 == 0) return;
        ptrdiff_t j1 = 0;
        for(const T* x1=x.cptr(); *x1==T(0); ++j1,++x1);
        if (j1 == 0 && j2 == N) DoLowerMultEqMV(A,x);
        else {
            TMVAssert(j1 < j2);
            const Ta* p22 = A.cptr() + j1*A.diagstep();
            const ptrdiff_t N22 = j2-j1;
            VectorView<T> x2 = x.subVector(j1,j2);
            if (N22 > A.nlo()) {
                ConstBandMatrixView<Ta> A22(
                    p22,N22,N22,A.nlo(),0,
                    A.stepi(),A.stepj(),A.diagstep(),A.ct());
                if (j2 < N) {
                    const ptrdiff_t jx = j2-A.nlo();
                    const Ta* p32 = A.cptr() +
                        j2*A.diagstep() - A.nlo()*A.stepj();
                    if (j2+A.nlo() > N) {
                        ConstBandMatrixView<Ta> A32(
                            p32,N-j2,A.nlo(),0,A.nlo()-1,
                            A.stepi(),A.stepj(),A.diagstep(),A.ct());
                        UnitAMultMV1<false,false>(
                            A32,x.subVector(jx,j2),x.subVector(j2,N));
                    } else {
                        ConstBandMatrixView<Ta> A32(
                            p32,A.nlo(),A.nlo(),0,A.nlo()-1,
                            A.stepi(),A.stepj(),A.diagstep(),A.ct());
                        VectorView<T> x3x = x.subVector(j2,j2+A.nlo());
                        x3x = x.subVector(jx,j2);
                        DoUpperMultEqMV(A32,x3x);
                    }
                }
                DoLowerMultEqMV(A22,x2);
            } else {
                ConstBandMatrixView<Ta> A22(
                    p22,N22,N22,N22-1,0,
                    A.stepi(),A.stepj(),A.diagstep(),A.ct());
                if (j2 < N) {
                    const Ta* p32 = p22 + N22*A.stepi();
                    const ptrdiff_t M32 = (j2+A.nlo() > N) ? N-j2 : A.nlo();
                    ptrdiff_t newlo = A.nlo()-N22;
                    if (newlo >= M32) newlo = M32-1;
                    ConstBandMatrixView<Ta> A32(
                        p32,M32,N22,newlo,N22-1,
                        A.stepi(),A.stepj(),A.diagstep(),A.ct());
                    UnitAMultMV1<false,false>(A32,x2,x.subVector(j2,j2+M32));
                } 
                DoLowerMultEqMV(A22,x2);
            }
        }
    }

#ifdef BLAS
    template <class T, class Ta> 
    static inline void BlasMultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    { 
        if (A.nlo() == 0) NonBlasUpperMultEqMV(A,x);
        else NonBlasLowerMultEqMV(A,x);
    }
#ifdef INST_DOUBLE
    template <> 
    void BlasMultEqMV( 
        const GenBandMatrix<double>& A, VectorView<double> x)
    {
        bool up = A.nlo()==0;
        int n=A.colsize();
        int lohi = up ? A.nhi() : A.nlo();
        int aoffset = 
            up && BlasIsCM(A) ? A.nhi() : 
            !up && !BlasIsCM(A) ? A.nlo() : 0;
        int ds = A.diagstep();
        int xs = x.step();
        BLASNAME(dtbmv) (
            BLASCM BlasIsCM(A) == up?BLASCH_UP:BLASCH_LO, 
            BlasIsCM(A) ? BLASCH_NT : BLASCH_T, BLASCH_NU,
            BLASV(n),BLASV(lohi),BLASP(A.cptr()-aoffset),BLASV(ds),
            BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasMultEqMV(
        const GenBandMatrix<std::complex<double> >& A,
        VectorView<std::complex<double> > x)
    {
        bool up = A.nlo()==0;
        int n=A.colsize();
        int lohi = up ? A.nhi() : A.nlo();
        int aoffset = 
            up && BlasIsCM(A) ? A.nhi() : 
            !up && !BlasIsCM(A) ? A.nlo() : 0;
        int ds = A.diagstep();
        int xs = x.step();
        if (BlasIsCM(A) && A.isconj()) {
#ifdef CBLAS
            BLASNAME(ztbmv) (
                BLASRM up ? BLASCH_LO : BLASCH_UP, BLASCH_CT, BLASCH_NU,
                BLASV(n),BLASV(lohi),BLASP(A.cptr()-A.nhi()),
                BLASV(ds),BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
#else
            x.conjugateSelf();
            BLASNAME(ztbmv) (
                BLASCM up?BLASCH_UP:BLASCH_LO, BLASCH_NT, BLASCH_NU,
                BLASV(n),BLASV(lohi),BLASP(A.cptr()-aoffset),BLASV(ds),
                BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
            x.conjugateSelf();
#endif
        } else {
            BLASNAME(ztbmv) (
                BLASCM BlasIsCM(A) == up?BLASCH_UP:BLASCH_LO, 
                BlasIsCM(A) ? BLASCH_NT : A.isconj() ? BLASCH_CT : BLASCH_T,
                BLASCH_NU,
                BLASV(n),BLASV(lohi),BLASP(A.cptr()-aoffset),BLASV(ds),
                BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasMultEqMV( 
        const GenBandMatrix<double>& A,
        VectorView<std::complex<double> > x)
    {
        bool up = A.nlo()==0;
        int n=A.colsize();
        int lohi = up ? A.nhi() : A.nlo();
        int aoffset = 
            up && BlasIsCM(A) ? A.nhi() : 
            !up && !BlasIsCM(A) ? A.nlo() : 0;
        int ds = A.diagstep();
        int xs = 2*x.step();
        BLASNAME(dtbmv) (
            BLASCM BlasIsCM(A) == up?BLASCH_UP:BLASCH_LO, 
            BlasIsCM(A) ? BLASCH_NT : BLASCH_T, BLASCH_NU,
            BLASV(n),BLASV(lohi),BLASP(A.cptr()-aoffset),BLASV(ds),
            BLASP((double*)x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
        BLASNAME(dtbmv) (
            BLASCM BlasIsCM(A) == up?BLASCH_UP:BLASCH_LO, 
            BlasIsCM(A) ? BLASCH_NT : BLASCH_T, BLASCH_NU,
            BLASV(n),BLASV(lohi),BLASP(A.cptr()-aoffset),BLASV(ds),
            BLASP((double*)x.ptr()+1),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void BlasMultEqMV( 
        const GenBandMatrix<float>& A, VectorView<float> x)
    {
        bool up = A.nlo()==0;
        int n=A.colsize();
        int lohi = up ? A.nhi() : A.nlo();
        int aoffset = 
            up && BlasIsCM(A) ? A.nhi() : 
            !up && !BlasIsCM(A) ? A.nlo() : 0;
        int ds = A.diagstep();
        int xs = x.step();
        BLASNAME(stbmv) (
            BLASCM BlasIsCM(A) == up?BLASCH_UP:BLASCH_LO, 
            BlasIsCM(A) ? BLASCH_NT : BLASCH_T, BLASCH_NU,
            BLASV(n),BLASV(lohi),BLASP(A.cptr()-aoffset),BLASV(ds),
            BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasMultEqMV(
        const GenBandMatrix<std::complex<float> >& A,
        VectorView<std::complex<float> > x)
    {
        bool up = A.nlo()==0;
        int n=A.colsize();
        int lohi = up ? A.nhi() : A.nlo();
        int aoffset = 
            up && BlasIsCM(A) ? A.nhi() : 
            !up && !BlasIsCM(A) ? A.nlo() : 0;
        int ds = A.diagstep();
        int xs = x.step();
        //cout<<"Before ctbmv\n";
        //cout<<"up = "<<up<<std::endl;
        //cout<<"n = "<<n<<std::endl;
        //cout<<"lohi = "<<lohi<<std::endl;
        //cout<<"aoffset = "<<aoffset<<std::endl;
        //cout<<"ds = "<<ds<<std::endl;
        //cout<<"xs = "<<xs<<std::endl;
        if (BlasIsCM(A) && A.isconj()) {
            //cout<<"cm && conj\n";
#ifdef CBLAS
            BLASNAME(ctbmv) (
                BLASRM up ? BLASCH_LO : BLASCH_UP, BLASCH_CT, BLASCH_NU,
                BLASV(n),BLASV(lohi),BLASP(A.cptr()-A.nhi()),
                BLASV(ds),BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
#else
            x.conjugateSelf();
            BLASNAME(ctbmv) (
                BLASCM up?BLASCH_UP:BLASCH_LO, BLASCH_NT, BLASCH_NU,
                BLASV(n),BLASV(lohi),BLASP(A.cptr()-aoffset),BLASV(ds),
                BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
            x.conjugateSelf();
#endif
        } else {
            //cout<<"!cm || !conj\n";
            BLASNAME(ctbmv) (
                BLASCM BlasIsCM(A) == up?BLASCH_UP:BLASCH_LO, 
                BlasIsCM(A) ? BLASCH_NT : A.isconj() ? BLASCH_CT : BLASCH_T, 
                BLASCH_NU,
                BLASV(n),BLASV(lohi),BLASP(A.cptr()-aoffset),BLASV(ds),
                BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasMultEqMV( 
        const GenBandMatrix<float>& A,
        VectorView<std::complex<float> > x)
    {
        bool up = A.nlo()==0;
        int n=A.colsize();
        int lohi = up ? A.nhi() : A.nlo();
        int aoffset = 
            up && BlasIsCM(A) ? A.nhi() : 
            !up && !BlasIsCM(A) ? A.nlo() : 0;
        int ds = A.diagstep();
        int xs = 2*x.step();
        BLASNAME(stbmv) (
            BLASCM BlasIsCM(A) == up?BLASCH_UP:BLASCH_LO, 
            BlasIsCM(A) ? BLASCH_NT : BLASCH_T, BLASCH_NU,
            BLASV(n),BLASV(lohi),BLASP(A.cptr()-aoffset),BLASV(ds),
            BLASP((float*)x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
        BLASNAME(stbmv) (
            BLASCM BlasIsCM(A) == up?BLASCH_UP:BLASCH_LO, 
            BlasIsCM(A) ? BLASCH_NT : BLASCH_T, BLASCH_NU,
            BLASV(n),BLASV(lohi),BLASP(A.cptr()-aoffset),BLASV(ds),
            BLASP((float*)x.ptr()+1),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
#endif 
#endif // BLAS

    template <class T, class Ta> 
    static void MultEqMV(
        const GenBandMatrix<Ta>& A, VectorView<T> x)
    {
#ifdef XDEBUG
        cout<<"Start MultEqMV\n";
        Vector<T> x0 = x;
        Matrix<Ta> A0 = A;
        Vector<T> x2 = A0 * x0;
#endif
        TMVAssert(A.isSquare());
        TMVAssert(A.colsize() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(A.nlo() == 0 || A.nhi() == 0);

        if (x.isconj()) MultEqMV(A.conjugate(),x.conjugate());
        else {
#ifdef BLAS
            if ( !(BlasIsRM(A) || BlasIsCM(A)) ) {
                //cout<<"!(rm or cm)\n";
                BandMatrix<Ta,ColMajor> AA = A;
                BlasMultEqMV(AA,x);
            } else if (SameStorage(A,x) || x.step() != 1) {
                //cout<<"copy x\n";
                Vector<T> xx = x;
                BlasMultEqMV(A,xx.view());
                x = xx;
            } else {
                //cout<<"normal\n";
                BlasMultEqMV(A,x);
            }
#else
            if (A.nlo() == 0) NonBlasUpperMultEqMV(A,x);
            else NonBlasLowerMultEqMV(A,x);
#endif
        }

#ifdef XDEBUG
        cout<<"Done MultEqMV\n";
        if (!(Norm(x-x2) <= 0.001*(Norm(A0)*Norm(x0)))) {
            cerr<<"MultEqMV: \n";
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"x = "<<TMV_Text(x)<<" step "<<x.step()<<"  "<<x0<<endl;
            cerr<<"-> x = "<<x<<endl;
            cerr<<"x2 = "<<x2<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tx> 
    void MultMV(
        const T alpha, const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y)
    // y (+)= alpha * A * x
    { 
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
#ifdef XDEBUG
        cout<<"Start Band: MultMV\n";
        cout<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"  "<<A<<endl;
        cout<<"x = "<<TMV_Text(x)<<"  "<<x.cptr()<<"  step "<<x.step()<<"  "<<x<<endl;
        cout<<"y = "<<TMV_Text(y)<<"  "<<y.cptr()<<"  step "<<y.step()<<"  "<<y<<endl;
        cout<<"alpha = "<<alpha<<", add = "<<add<<endl;
        Vector<T> y0 = y;
        cout<<"y0 = "<<y0<<std::endl;
        Vector<Tx> x0 = x;
        cout<<"x0 = "<<x0<<std::endl;
        Matrix<Ta> A0 = A;
        cout<<"A0 = "<<A0<<std::endl;
        Vector<T> y2 = alpha*A0*x0;
        cout<<"y2 = "<<y2<<std::endl;
        if (add) y2 += y0;
        cout<<"y2 => "<<y2<<std::endl;
#endif

        if (y.size() > 0) {
            if (x.size()==0 || alpha==T(0)) {
                if (!add) y.setZero();
            } else if (A.rowsize() > A.colsize()+A.nhi()) {
                MultMV<add>(
                    alpha,A.colRange(0,A.colsize()+A.nhi()),
                    x.subVector(0,A.colsize()+A.nhi()),y);
            } else if (A.colsize() > A.rowsize()+A.nlo()) {
                MultMV<add>(
                    alpha,A.rowRange(0,A.rowsize()+A.nlo()),
                    x,y.subVector(0,A.rowsize()+A.nlo()));
                if (!add) 
                    y.subVector(A.rowsize()+A.nlo(),A.colsize()).setZero();
            } else if (A.isSquare() && (A.nlo() == 0 || A.nhi() == 0)) {
                if (A.nlo() == 0 && A.nhi() == 0)
                    MultMV<add>(alpha,DiagMatrixViewOf(A.diag()),x,y);
                else if (!add && y.step() == 1) {
                    y = alpha * x;
                    MultEqMV(A,y);
                } else {
                    Vector<T> xx = alpha*x;
                    MultEqMV(A,xx.view());
                    if (add) y += xx;
                    else y = xx;
                }
            } else {
                if (SameStorage(y,A)) {
                    Vector<T> yy(y.size());
                    DoMultMV<false>(T(1),A,x,yy.view());
                    if (add) y += alpha*yy;
                    else y = alpha*yy;
                } else {
                    DoMultMV<add>(alpha,A,x,y);
                }
            }
        }
#ifdef XDEBUG
        cout<<"->y = "<<y<<endl;
        if (!(Norm(y2-y) <= 0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(x0)+
                                (add?Norm(y0):TMV_RealType(T)(0))))) {
            cerr<<"MultMV: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A.cptr()<<"  "<<A0<<endl;
            cerr<<"x = "<<TMV_Text(x)<<"  "<<x.cptr()<<" step "<<x.step()<<"  "<<x0<<endl;
            cerr<<"y = "<<TMV_Text(y)<<"  "<<y.cptr()<<" step "<<y.step()<<"  "<<y0<<endl;
            cerr<<"--> y = "<<y<<endl;
            cerr<<"y2 = "<<y2<<endl;
            cerr<<"Norm(diff) = "<<Norm(y2-y)<<endl;
            cerr<<"abs(alpha)*Norm(A0)*Norm(x0) = "<<TMV_ABS(alpha)*Norm(A0)*Norm(x0)<<endl;
            cerr<<"Norm(y0) = "<<Norm(y0)<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_MultBV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


