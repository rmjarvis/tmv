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
#include "TMV_MultMV.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_VIt.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

// CBLAS trick of using RowMajor with ConjTrans when we have a 
// case of A.conjugate() * x doesn't seem to be working with MKL 10.2.2.
// I haven't been able to figure out why.  (e.g. Is it a bug in the MKL
// code, or am I doing something wrong?)  So for now, just disable it.
#ifdef CBLAS
#undef CBLAS
#endif

namespace tmv {

    template <class T> const T* MatrixComposite<T>::cptr() const
    {
        if (!itsm.get()) {
            ptrdiff_t len = this->colsize()*this->rowsize();
            itsm.resize(len);
            MatrixView<T>(itsm.get(),this->colsize(),this->rowsize(),
                          stepi(),stepj(),NonConj,len 
                          TMV_FIRSTLAST1(itsm.get(),itsm.get()+len) ) = *this;
        }
        return itsm.get();
    }

    template <class T> ptrdiff_t MatrixComposite<T>::stepi() const 
    { return 1; }

    template <class T> ptrdiff_t MatrixComposite<T>::stepj() const 
    { return this->colsize(); }

    template <class T> ptrdiff_t MatrixComposite<T>::ls() const 
    { return this->rowsize() * this->colsize(); }

    // 
    //
    // MultMV
    //

    // These routines are designed to work even if y has the same storage
    // as either x or the first row/column of A.

    template <bool add, bool cx, bool ca, bool rm, class T, class Ta, class Tx>
    static void RowMultMV(
        const GenMatrix<Ta>& A, const GenVector<Tx>& x, VectorView<T> y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(y.ct()==NonConj);
        TMVAssert(x.step() == 1);
        TMVAssert(y.step() == 1);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(cx == x.isconj());
        TMVAssert(ca == A.isconj());

        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();
        const ptrdiff_t si = A.stepi();
        const ptrdiff_t sj = (rm ? 1 : A.stepj());

        const Ta* Ai0 = A.cptr();
        const Tx*const x0 = x.cptr();
        T* yi = y.ptr();

        for(ptrdiff_t i=M; i>0; --i,++yi,Ai0+=si) {
            // *yi += A.row(i) * x

            const Ta* Aij = Ai0;
            const Tx* xj = x0;
            register T temp(0);
            for(ptrdiff_t j=N; j>0; --j,++xj,(rm?++Aij:Aij+=sj))
                temp += 
                    (cx ? TMV_CONJ(*xj) : *xj) *
                    (ca ? TMV_CONJ(*Aij) : *Aij);

#ifdef TMVFLDEBUG
            TMVAssert(yi >= y._first);
            TMVAssert(yi < y._last);
#endif
            if (add) *yi += temp;
            else *yi = temp;
        }
    }

    template <bool add, bool cx, bool ca, bool cm, class T, class Ta, class Tx> 
    static void ColMultMV(
        const GenMatrix<Ta>& A, const GenVector<Tx>& x, VectorView<T> y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(y.ct()==NonConj);
        TMVAssert(x.step() == 1);
        TMVAssert(y.step() == 1);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(cx == x.isconj());
        TMVAssert(ca == A.isconj());
        TMVAssert(cm == A.iscm());

        const ptrdiff_t M = A.colsize();
        ptrdiff_t N = A.rowsize();
        const ptrdiff_t si = (cm ? 1 : A.stepi());
        const ptrdiff_t sj = A.stepj();

        const Ta* A0j = A.cptr();
        const Tx* xj = x.cptr();
        T*const y0 = y.ptr();

        if (!add) {
            if (*xj == Tx(0)) {
                y.setZero();
            } else {
                const Ta* Aij = A0j;
                T* yi = y0;
                const Tx xjval = (cx ? TMV_CONJ(*xj) : *xj);
                for(ptrdiff_t i=M; i>0; --i,++yi,(cm?++Aij:Aij+=si)) {
#ifdef TMVFLDEBUG
                    TMVAssert(yi >= y._first);
                    TMVAssert(yi < y._last);
#endif
                    *yi = xjval * (ca ? TMV_CONJ(*Aij) : *Aij);
                }
            }
            ++xj; A0j+=sj; --N;
        }

        for(; N>0; --N,++xj,A0j+=sj) {
            // y += *xj * A.col(j)
            if (*xj != Tx(0)) {
                const Ta* Aij = A0j;
                T* yi = y0;
                const Tx xjval = (cx ? TMV_CONJ(*xj) : *xj);
                for(ptrdiff_t i=M; i>0; --i,++yi,(cm?++Aij:Aij+=si)) {
#ifdef TMVFLDEBUG
                    TMVAssert(yi >= y._first);
                    TMVAssert(yi < y._last);
#endif
                    *yi += xjval * (ca ? TMV_CONJ(*Aij) : *Aij);
                }
            }
        }
    }

    template <bool add, bool cx, class T, class Ta, class Tx> 
    void UnitAMultMV1(
        const GenMatrix<Ta>& A, const GenVector<Tx>& x, VectorView<T> y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(y.ct() == NonConj);
        TMVAssert(x.step() == 1);
        TMVAssert(y.step() == 1);
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
        else if ( A.rowsize() >= A.colsize() )
            if (A.isconj())
                RowMultMV<add,cx,true,false>(A,x,y);
            else
                RowMultMV<add,cx,false,false>(A,x,y);
        else 
            if (A.isconj())
                ColMultMV<add,cx,true,false>(A,x,y);
            else
                ColMultMV<add,cx,false,false>(A,x,y);
    }

    template <bool add, bool cx, class T, class Ta, class Tx> 
    static void UnitAMultMV(
        const GenMatrix<Ta>& A, const GenVector<Tx>& x, VectorView<T> y)
    {
#ifdef XDEBUG
        cout<<"Start UnitAMultMV: \n";
        cout<<"add = "<<add<<endl;
        cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
        cout<<"x = "<<TMV_Text(x)<<" step "<<x.step()<<"  "<<x<<endl;
        cout<<"y = "<<TMV_Text(y)<<" step "<<y.step()<<"  "<<y<<endl;
        Vector<Tx> x0 = x;
        Vector<T> y0 = y;
        Matrix<Ta> A0 = A;
        Vector<T> y2 = y;
        for(ptrdiff_t i=0;i<y.size();i++) {
            if (add)
                y2(i) += (A.row(i) * x0);
            else
                y2(i) = (A.row(i) * x0);
        }
        cout<<"y2 = "<<y2<<endl;
#endif
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
        for(const Tx* x1=x.cptr(); *x1==Tx(0); ++j1,++x1) {}
        TMVAssert(j1 !=j2);
        if (j1 == 0 && j2 == N) UnitAMultMV1<add,cx>(A,x,y);
        else UnitAMultMV1<add,cx>(A.colRange(j1,j2),x.subVector(j1,j2),y);

#ifdef XDEBUG
        cout<<"y => "<<y<<endl;
        if (!(Norm(y-y2) <=
              0.001*(Norm(A0)*Norm(x0)+
                     (add?Norm(y0):TMV_RealType(T)(0))))) {
            cerr<<"MultMV: \n";
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A);
            if (A.rowsize() < 30 && A.colsize() < 30) cerr<<"  "<<A0;
            else cerr<<"  "<<A.colsize()<<" x "<<A.rowsize();
            cerr<<endl<<"x = "<<TMV_Text(x)<<" step "<<x.step();
            if (x.size() < 30) cerr<<"  "<<x0;
            cerr<<endl<<"y = "<<TMV_Text(y)<<" step "<<y.step();
            if (add && y.size() < 30) cerr<<"  "<<y0;
            cerr<<endl<<"Aptr = "<<A.cptr();
            cerr<<", xptr = "<<x.cptr()<<", yptr = "<<y.cptr()<<endl;
            if (y.size() < 200) {
                cerr<<"--> y = "<<y<<endl;
                cerr<<"y2 = "<<y2<<endl;
            } else {
                ptrdiff_t imax;
                (y-y2).maxAbsElement(&imax);
                cerr<<"y("<<imax<<") = "<<y(imax)<<endl;
                cerr<<"y2("<<imax<<") = "<<y2(imax)<<endl;
            }
            cerr<<"Norm(A0) = "<<Norm(A0)<<endl;
            cerr<<"Norm(x0) = "<<Norm(x0)<<endl;
            if (add) cerr<<"Norm(y0) = "<<Norm(y0)<<endl;
            cerr<<"|A0|*|x0|+?|y0| = "<<
                Norm(A0)*Norm(x0)+
                (add?Norm(y0):TMV_RealType(T)(0))<<endl;
            cerr<<"Norm(y-y2) = "<<Norm(y-y2)<<endl;
            cerr<<"NormInf(y-y2) = "<<NormInf(y-y2)<<endl;
            cerr<<"Norm1(y-y2) = "<<Norm1(y-y2)<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tx> static void NonBlasMultMV(
        const T alpha, const GenMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y)
        // y (+)= alpha * A * x
    {
#ifdef XDEBUG
        cout<<"Start MultMV: alpha = "<<alpha<<endl;
        cout<<"add = "<<add<<endl;
        cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
        cout<<"x = "<<TMV_Text(x)<<" step "<<x.step()<<"  "<<x<<endl;
        cout<<"y = "<<TMV_Text(y)<<" step "<<y.step()<<"  "<<y<<endl;
        Vector<Tx> x0 = x;
        Vector<T> y0 = y;
        Matrix<Ta> A0 = A;
        Vector<T> y2 = y;
        for(ptrdiff_t i=0;i<y.size();i++) {
            if (add)
                y2(i) += alpha * (A.row(i) * x0);
            else
                y2(i) = alpha * (A.row(i) * x0);
        }
        cout<<"y2 = "<<y2<<endl;
#endif
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(alpha != T(0));
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(y.ct() == NonConj);

        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();

        if (x.step() != 1 || SameStorage(x,y) ||
            (alpha != TMV_RealType(T)(1) && y.step() == 1 && M/4 >= N)) {
            // This last check is taken from the ATLAS version of this code.
            // Apparently M = 4N is the dividing line between applying alpha
            // here versus at the end when adding Ax to y
            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                Vector<Tx> xx = TMV_REAL(alpha)*x;
                if (y.step()!=1) {
                    Vector<T> yy(y.size());
                    UnitAMultMV<false,false>(A,xx,yy.view());
                    if (add) y += yy;
                    else y = yy;
                } 
                else 
                    UnitAMultMV<add,false>(A,xx,y);
            } else {
                Vector<T> xx = alpha*x;
                if (y.step() != 1) {
                    Vector<T> yy(y.size());
                    UnitAMultMV<false,false>(A,xx,yy.view());
                    if (add) y += yy;
                    else y = yy;
                } 
                else 
                    UnitAMultMV<add,false>(A,xx,y);
            }
        } else if (y.step() != 1 || alpha != TMV_RealType(T)(1)) {
            Vector<T> yy(y.size());
            if (x.isconj())
                UnitAMultMV<false,true>(A,x,yy.view());
            else
                UnitAMultMV<false,false>(A,x,yy.view());
            if (add) y += alpha*yy;
            else y = alpha*yy;
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
        cout<<"y => "<<y<<endl;
        if (!(Norm(y-y2) <=
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(x0)+
                     (add?Norm(y0):TMV_RealType(T)(0))))) {
            cerr<<"MultMV: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A);
            if (A.rowsize() < 30 && A.colsize() < 30) cerr<<"  "<<A0;
            else cerr<<"  "<<A.colsize()<<" x "<<A.rowsize();
            cerr<<endl<<"x = "<<TMV_Text(x)<<" step "<<x.step();
            if (x.size() < 30) cerr<<"  "<<x0;
            cerr<<endl<<"y = "<<TMV_Text(y)<<" step "<<y.step();
            if (add && y.size() < 30) cerr<<"  "<<y0;
            cerr<<endl<<"Aptr = "<<A.cptr();
            cerr<<", xptr = "<<x.cptr()<<", yptr = "<<y.cptr()<<endl;
            if (y.size() < 200) {
                cerr<<"--> y = "<<y<<endl;
                cerr<<"y2 = "<<y2<<endl;
            } else {
                ptrdiff_t imax;
                (y-y2).maxAbsElement(&imax);
                cerr<<"y("<<imax<<") = "<<y(imax)<<endl;
                cerr<<"y2("<<imax<<") = "<<y2(imax)<<endl;
            }
            cerr<<"Norm(A0) = "<<Norm(A0)<<endl;
            cerr<<"Norm(x0) = "<<Norm(x0)<<endl;
            if (add) cerr<<"Norm(y0) = "<<Norm(y0)<<endl;
            cerr<<"|alpha|*|A0|*|x0|+?|y0| = "<<
                TMV_ABS(alpha)*Norm(A0)*Norm(x0)+
                (add?Norm(y0):TMV_RealType(T)(0))<<endl;
            cerr<<"Norm(y-y2) = "<<Norm(y-y2)<<endl;
            cerr<<"NormInf(y-y2) = "<<NormInf(y-y2)<<endl;
            cerr<<"Norm1(y-y2) = "<<Norm1(y-y2)<<endl;
            abort();
        }
#endif
    }

#ifdef BLAS
    template <class T, class Ta, class Tx> static inline void BlasMultMV(
        const T alpha, const GenMatrix<Ta>& A,
        const GenVector<Tx>& x, const int beta, VectorView<T> y)
    { 
        if (beta == 0) NonBlasMultMV<false>(alpha,A,x,y); 
        else NonBlasMultMV<true>(alpha,A,x,y); 
    }
#ifdef INST_DOUBLE
    template <> void BlasMultMV(
        const double alpha, const GenMatrix<double>& A,
        const GenVector<double>& x, const int beta, VectorView<double> y)
    {
        int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
        int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
        int lda = BlasIsCM(A) ? A.stepj() : A.stepi();
        if (lda < m) { TMVAssert(n==1); lda = m; }
        int xs = x.step();
        int ys = y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const double* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        double* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        // Some BLAS implementations seem to have trouble if the 
        // input y has a nan in it.
        // They propagate the nan into the  output.
        // I guess they strictly interpret y = beta*y + alpha*m*x,
        // so if beta = 0, then beta*nan = nan.
        // Anyway, to fix this problem, we always use beta=1, and just
        // zero out y before calling the blas function if beta is 0.
        if (beta == 0) y.setZero();
        double xbeta(1);

#if 0
        std::cout<<"Before dgemv"<<std::endl;
        std::cout<<"A = "<<A<<std::endl;
        std::cout<<"x = "<<x<<std::endl;
        std::cout<<"y = "<<y<<std::endl;
        std::cout<<"m = "<<m<<std::endl;
        std::cout<<"n = "<<n<<std::endl;
        std::cout<<"lda = "<<lda<<std::endl;
        std::cout<<"xs = "<<xs<<std::endl;
        std::cout<<"ys = "<<ys<<std::endl;
        std::cout<<"alpha = "<<alpha<<std::endl;
        std::cout<<"beta = "<<xbeta<<std::endl;
        std::cout<<"aptr = "<<A.cptr()<<std::endl;
        std::cout<<"xp = "<<xp<<std::endl;
        std::cout<<"yp = "<<yp<<std::endl;
        std::cout<<"NT = "<<(A.isrm()?'T':'N')<<std::endl;
        if (A.isrm()) {
            std::cout<<"x.size = "<<x.size()<<std::endl;
            std::cout<<"x = ";
            for(int i=0;i<m;i++) std::cout<<*(xp+i*xs)<<" ";
            std::cout<<std::endl;
            std::cout<<"y.size = "<<y.size()<<std::endl;
            std::cout<<"y = ";
            for(int i=0;i<n;i++) std::cout<<*(yp+i*ys)<<" ";
            std::cout<<std::endl;
            std::cout<<"A.size = "<<A.colsize()<<','<<A.rowsize()<<std::endl;
            std::cout<<"A = ";
            for(int i=0;i<n*lda;i++) std::cout<<*(A.cptr()+i)<<" ";
            std::cout<<std::endl;
        }
#endif
        BLASNAME(dgemv) (
            BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASV(xbeta),BLASP(yp),BLASV(ys) BLAS1);
#if 0
        std::cout<<"After dgemv"<<std::endl;
#endif
    }
    template <> void BlasMultMV(
        const std::complex<double> alpha,
        const GenMatrix<std::complex<double> >& A,
        const GenVector<std::complex<double> >& x,
        const int beta, VectorView<std::complex<double> > y)
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
        int lda = BlasIsCM(A) ? A.stepj() : A.stepi();
        if (lda < m) { TMVAssert(n==1); lda = m; }
        int xs = x.step();
        int ys = y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const std::complex<double>* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        std::complex<double>* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        if (beta == 0) y.setZero();
        std::complex<double> xbeta(1);
#if 0
        std::cout<<"Before zgemv"<<std::endl;
        std::cout<<"A = "<<A<<std::endl;
        std::cout<<"x = "<<x<<std::endl;
        std::cout<<"y = "<<y<<std::endl;
        std::cout<<"m = "<<m<<std::endl;
        std::cout<<"n = "<<n<<std::endl;
        std::cout<<"lda = "<<lda<<std::endl;
        std::cout<<"xs = "<<xs<<std::endl;
        std::cout<<"ys = "<<ys<<std::endl;
        std::cout<<"alpha = "<<alpha<<std::endl;
        std::cout<<"beta = "<<xbeta<<std::endl;
        std::cout<<"aptr = "<<A.cptr()<<std::endl;
        std::cout<<"xp = "<<xp<<std::endl;
        std::cout<<"yp = "<<yp<<std::endl;
#endif
        if (A.isconj() && BlasIsCM(A)) {
#ifdef CBLAS
            TMV_SWAP(m,n);
            BLASNAME(zgemv) (
                BLASRM BLASCH_CT,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys) BLAS1);
#else
            std::complex<double> ca = TMV_CONJ(alpha);
            if (x.isconj()) {
                y.conjugateSelf();
                BLASNAME(zgemv) (
                    BLASCM BLASCH_NT,
                    BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                    BLAS1);
                y.conjugateSelf();
            } else {
                Vector<std::complex<double> > xx = ca*x.conjugate();
                ca = std::complex<double>(1);
                xs = 1;
                xp = xx.cptr();
                y.conjugateSelf();
                BLASNAME(zgemv) (
                    BLASCM BLASCH_NT,
                    BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                    BLAS1);
                y.conjugateSelf();
            }
#endif
        } else {
            BLASNAME(zgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                BLAS1);
        }
    }
    template <> void BlasMultMV(
        const std::complex<double> alpha,
        const GenMatrix<std::complex<double> >& A,
        const GenVector<double>& x,
        const int beta, VectorView<std::complex<double> > y)
    {
        if (BlasIsCM(A)) {
            if (y.step() != 1) {
                Vector<std::complex<double> > yy(y.size());
                BlasMultMV(std::complex<double>(1),A,x,0,yy.view());
                if (beta == 0) y = alpha*yy;
                else y += alpha*yy;
            } else {
                if (beta == 0) {
                    int m = 2*A.colsize();
                    int n = A.rowsize();
                    int lda = 2*A.stepj();
                    if (lda < m) { TMVAssert(n==1); lda = m; }
                    int xs = x.step();
                    int ys = 1;
                    const double* xp = x.cptr();
                    if (xs < 0) xp += (x.size()-1)*xs;
                    double* yp = (double*) y.ptr();
                    double xalpha(1);
                    y.setZero();
                    double xbeta(1);
                    BLASNAME(dgemv) (
                        BLASCM BLASCH_NT,
                        BLASV(m),BLASV(n),BLASV(xalpha),
                        BLASP((double*)A.cptr()),BLASV(lda),
                        BLASP(xp),BLASV(xs),BLASV(xbeta),
                        BLASP(yp),BLASV(ys) BLAS1);
                    if (A.isconj()) y.conjugateSelf();
                    y *= alpha;
                } else if (A.isconj()) {
                    Vector<std::complex<double> > yy(y.size());
                    BlasMultMV(
                        std::complex<double>(1),A.conjugate(),x,0,yy.view());
                    y += alpha*yy.conjugate();
                } else if (TMV_IMAG(alpha) == 0.) {
                    int m = 2*A.colsize();
                    int n = A.rowsize();
                    int lda = 2*A.stepj();
                    if (lda < m) { TMVAssert(n==1); lda = m; }
                    int xs = x.step();
                    int ys = 1;
                    const double* xp = x.cptr();
                    if (xs < 0) xp += (x.size()-1)*xs;
                    double* yp = (double*) y.ptr();
                    if (ys < 0) yp += (y.size()-1)*ys;
                    double xalpha(TMV_REAL(alpha));
                    double xbeta(1);
                    BLASNAME(dgemv) (
                        BLASCM BLASCH_NT,
                        BLASV(m),BLASV(n),BLASV(xalpha),
                        BLASP((double*)A.cptr()),BLASV(lda),
                        BLASP(xp),BLASV(xs),BLASV(xbeta),
                        BLASP(yp),BLASV(ys) BLAS1);
                } else {
                    Vector<std::complex<double> > yy(y.size());
                    BlasMultMV(std::complex<double>(1),A,x,0,yy.view());
                    y += alpha*yy;
                }
            }
        } else { // A.isrm
            BlasMultMV(alpha,A,Vector<std::complex<double> >(x),beta,y);
        }
    }
    template <> void BlasMultMV(
        const std::complex<double> alpha,
        const GenMatrix<double>& A,
        const GenVector<std::complex<double> >& x,
        const int beta, VectorView<std::complex<double> > y)
    {
        if (beta == 0) {
            int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
            int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
            int lda = BlasIsCM(A) ? A.stepj() : A.stepi();
            if (lda < m) { TMVAssert(n==1); lda = m; }
            int xs = 2*x.step();
            int ys = 2*y.step();
            if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
            if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
            const double* xp = (const double*) x.cptr();
            if (xs < 0) xp += (x.size()-1)*xs;
            double* yp = (double*) y.ptr();
            if (ys < 0) yp += (y.size()-1)*ys;
            double xalpha(1);
            y.setZero();
            double xbeta(1);
            BLASNAME(dgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(dgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp+1),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
            if (x.isconj()) y.conjugateSelf();
            y *= alpha;
        } else if (TMV_IMAG(alpha) == 0. && !x.isconj()) {
            int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
            int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
            int lda = BlasIsCM(A) ? A.stepj() : A.stepi();
            if (lda < m) { TMVAssert(n==1); lda = m; }
            int xs = 2*x.step();
            int ys = 2*y.step();
            if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
            if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
            const double* xp = (const double*) x.cptr();
            if (xs < 0) xp += (x.size()-1)*xs;
            double* yp = (double*) y.ptr();
            if (ys < 0) yp += (y.size()-1)*ys;
            double xalpha(TMV_REAL(alpha));
            double xbeta(1);
            BLASNAME(dgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(dgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp+1),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        } else {
            Vector<std::complex<double> > xx = alpha*x;
            BlasMultMV(std::complex<double>(1),A,xx,1,y);
        }
    }
    template <> void BlasMultMV(
        const std::complex<double> alpha,
        const GenMatrix<double>& A,
        const GenVector<double>& x,
        const int beta, VectorView<std::complex<double> > y)
    {
        int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
        int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
        int lda = BlasIsCM(A) ? A.stepj() : A.stepi();
        if (lda < m) { TMVAssert(n==1); lda = m; }
        int xs = x.step();
        int ys = 2*y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const double* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        double* yp = (double*) y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        double ar(TMV_REAL(alpha));
        double ai(TMV_IMAG(alpha));
        double xbeta(beta);
        if (beta == 0) y.setZero();
        if (ar != 0.) {
            BLASNAME(dgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(ar),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
        }
        if (ai != 0.) {
            BLASNAME(dgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(ai),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        }
    }
#endif
#ifdef INST_FLOAT
    template <> void BlasMultMV(
        const float alpha, const GenMatrix<float>& A,
        const GenVector<float>& x, const int beta, VectorView<float> y)
    {
        int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
        int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
        int lda = BlasIsCM(A) ? A.stepj() : A.stepi();
        if (lda < m) { TMVAssert(n==1); lda = m; }
        int xs = x.step();
        int ys = y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const float* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        float* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        if (beta == 0) y.setZero();
        float xbeta(1);
#if 0
        std::cout<<"Before sgemv:\n";
        std::cout<<"A = "<<A<<std::endl;
        std::cout<<"x = "<<x<<std::endl;
        std::cout<<"y = "<<y<<std::endl;
        std::cout<<"m,n = "<<m<<','<<n<<std::endl;
        std::cout<<"alpha,beta = "<<alpha<<','<<xbeta<<std::endl;
        std::cout<<"A.cptr = "<<A.cptr()<<std::endl;
        std::cout<<"xp = "<<xp<<std::endl;
        std::cout<<"yp = "<<yp<<std::endl;
        std::cout<<"lda,xs,ys = "<<lda<<','<<xs<<','<<ys<<std::endl;
        std::cout<<"cm = "<<BlasIsCM(A)<<std::endl;
#endif

        BLASNAME(sgemv) (
            BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASV(xbeta),BLASP(yp),BLASV(ys)
            BLAS1);
#if 0
        std::cout<<"After sgemv:"<<std::endl;
        std::cout<<"y -> "<<y<<std::endl;
#endif
    }
    template <> void BlasMultMV(
        const std::complex<float> alpha,
        const GenMatrix<std::complex<float> >& A,
        const GenVector<std::complex<float> >& x,
        const int beta, VectorView<std::complex<float> > y)
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
        int lda = BlasIsCM(A) ? A.stepj() : A.stepi();
        if (lda < m) { TMVAssert(n==1); lda = m; }
        int xs = x.step();
        int ys = y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const std::complex<float>* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        std::complex<float>* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        if (beta == 0) y.setZero();
        std::complex<float> xbeta(1);
#if 0
        std::cout<<"Before cgemv:\n";
        std::cout<<"A = "<<A<<std::endl;
        std::cout<<"x = "<<x<<std::endl;
        std::cout<<"y = "<<y<<std::endl;
        std::cout<<"m,n = "<<m<<','<<n<<std::endl;
        std::cout<<"alpha,beta = "<<alpha<<','<<xbeta<<std::endl;
        std::cout<<"A.cptr = "<<A.cptr()<<std::endl;
        std::cout<<"xp = "<<xp<<std::endl;
        std::cout<<"yp = "<<yp<<std::endl;
        std::cout<<"lda,xs,ys = "<<lda<<','<<xs<<','<<ys<<std::endl;
        std::cout<<"conj = "<<A.isconj()<<std::endl;
        std::cout<<"cm = "<<A.iscm()<<std::endl;
#endif
        if (A.isconj() && BlasIsCM(A)) {
#ifdef CBLAS
            TMV_SWAP(m,n);
            BLASNAME(cgemv) (
                BLASRM BLASCH_CT,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys) BLAS1);
#else
            std::complex<float> ca = TMV_CONJ(alpha);
            if (x.isconj()) {
                y.conjugateSelf();
                BLASNAME(cgemv) (
                    BLASCM BLASCH_NT,
                    BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                    BLAS1);
                y.conjugateSelf();
            } else {
                Vector<std::complex<float> > xx = ca*x.conjugate();
                ca = std::complex<float>(1);
                xs = 1;
                xp = xx.cptr();
                y.conjugateSelf();
                BLASNAME(cgemv) (
                    BLASCM BLASCH_NT,
                    BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                    BLAS1);
                y.conjugateSelf();
            }
#endif
        } else {
            BLASNAME(cgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                BLAS1);
        }
#if 0
        std::cout<<"After cgemv:"<<std::endl;
        std::cout<<"y -> "<<y<<std::endl;
#endif
    }
    template <> void BlasMultMV(
        const std::complex<float> alpha,
        const GenMatrix<std::complex<float> >& A,
        const GenVector<float>& x,
        const int beta, VectorView<std::complex<float> > y)
    {
        if (BlasIsCM(A)) {
            if (y.step() != 1) {
                Vector<std::complex<float> > yy(y.size());
                BlasMultMV(std::complex<float>(1),A,x,0,yy.view());
                if (beta == 0) y = alpha*yy;
                else y += alpha*yy;
            } else {
                if (beta == 0) {
                    int m = 2*A.colsize();
                    int n = A.rowsize();
                    int lda = 2*A.stepj();
                    if (lda < m) { TMVAssert(n==1); lda = m; }
                    int xs = x.step();
                    int ys = 1;
                    const float* xp = x.cptr();
                    if (xs < 0) xp += (x.size()-1)*xs;
                    float* yp = (float*) y.ptr();
                    float xalpha(1);
                    y.setZero();
                    float xbeta(1);
                    BLASNAME(sgemv) (
                        BLASCM BLASCH_NT,
                        BLASV(m),BLASV(n),BLASV(xalpha),
                        BLASP((float*)A.cptr()),BLASV(lda),
                        BLASP(xp),BLASV(xs),BLASV(xbeta),
                        BLASP(yp),BLASV(ys) BLAS1);
                    if (A.isconj()) y.conjugateSelf();
                    y *= alpha;
                } else if (A.isconj()) {
                    Vector<std::complex<float> > yy(y.size());
                    BlasMultMV(std::complex<float>(1),A.conjugate(),x,0,yy.view());
                    y += alpha*yy.conjugate();
                } else if (TMV_IMAG(alpha) == 0.F) {
                    int m = 2*A.colsize();
                    int n = A.rowsize();
                    int lda = 2*A.stepj();
                    if (lda < m) { TMVAssert(n==1); lda = m; }
                    int xs = x.step();
                    int ys = 1;
                    const float* xp = x.cptr();
                    if (xs < 0) xp += (x.size()-1)*xs;
                    float* yp = (float*) y.ptr();
                    float xalpha(TMV_REAL(alpha));
                    float xbeta(1);
                    BLASNAME(sgemv) (
                        BLASCM BLASCH_NT,
                        BLASV(m),BLASV(n),BLASV(xalpha),
                        BLASP((float*)A.cptr()),BLASV(lda),
                        BLASP(xp),BLASV(xs),BLASV(xbeta),
                        BLASP(yp),BLASV(ys) BLAS1);
                } else {
                    Vector<std::complex<float> > yy(y.size());
                    BlasMultMV(std::complex<float>(1),A,x,0,yy.view());
                    y += alpha*yy;
                }
            } 
        } else { // A.isrm
            BlasMultMV(alpha,A,Vector<std::complex<float> >(x),beta,y);
        }
    }
    template <> void BlasMultMV(
        const std::complex<float> alpha,
        const GenMatrix<float>& A,
        const GenVector<std::complex<float> >& x,
        const int beta, VectorView<std::complex<float> > y)
    {
        if (beta == 0) {
            int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
            int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
            int lda = BlasIsCM(A) ? A.stepj() : A.stepi();
            if (lda < m) { TMVAssert(n==1); lda = m; }
            int xs = 2*x.step();
            int ys = 2*y.step();
            if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
            if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
            const float* xp = (const float*) x.cptr();
            if (xs < 0) xp += (x.size()-1)*xs;
            float* yp = (float*) y.ptr();
            if (ys < 0) yp += (y.size()-1)*ys;
            float xalpha(1);
            y.setZero();
            float xbeta(1);
            BLASNAME(sgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(sgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp+1),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
            if (x.isconj()) y.conjugateSelf();
            y *= alpha;
        } else if (TMV_IMAG(alpha) == 0.F && !x.isconj()) {
            int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
            int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
            int lda = BlasIsCM(A) ? A.stepj() : A.stepi();
            if (lda < m) { TMVAssert(n==1); lda = m; }
            int xs = 2*x.step();
            int ys = 2*y.step();
            if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
            if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
            const float* xp = (const float*) x.cptr();
            if (xs < 0) xp += (x.size()-1)*xs;
            float* yp = (float*) y.ptr();
            if (ys < 0) yp += (y.size()-1)*ys;
            float xalpha(TMV_REAL(alpha));
            float xbeta(1);
            BLASNAME(sgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(sgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp+1),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        } else {
            Vector<std::complex<float> > xx = alpha*x;
            BlasMultMV(std::complex<float>(1),A,xx,1,y);
        }
    }
    template <> void BlasMultMV(
        const std::complex<float> alpha,
        const GenMatrix<float>& A,
        const GenVector<float>& x,
        const int beta, VectorView<std::complex<float> > y)
    {
        int m = BlasIsCM(A) ? A.colsize() : A.rowsize();
        int n = BlasIsCM(A) ? A.rowsize() : A.colsize();
        int lda = BlasIsCM(A) ? A.stepj() : A.stepi();
        if (lda < m) { TMVAssert(n==1); lda = m; }
        int xs = x.step();
        int ys = 2*y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const float* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        float* yp = (float*) y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        float ar(TMV_REAL(alpha));
        float ai(TMV_IMAG(alpha));
        if (beta == 0) y.setZero();
        float xbeta(1);
        if (ar != 0.F) {
            BLASNAME(sgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(ar),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
        }
        if (ai != 0.F) {
            BLASNAME(sgemv) (
                BLASCM BlasIsCM(A)?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(ai),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        }
    }
#endif
#endif // BLAS

    template <bool add, class T, class Ta, class Tx> 
    static void DoMultMV(
        const T alpha, const GenMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<T> y)
    {
#ifdef XDEBUG
        std::cout<<"Start DoMultMV\n";
        std::cout<<"alpha = "<<alpha<<std::endl;
        std::cout<<"add = "<<add<<std::endl;
        std::cout<<"A = "<<A<<std::endl;
        std::cout<<"x = "<<x<<std::endl;
        if (add) std::cout<<"y = "<<y<<std::endl;
#endif

        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(alpha != T(0));
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(y.ct() == NonConj);

#ifdef BLAS
        if (x.step() == 0) {
            if (x.size() <= 1) 
                DoMultMV<add>(
                    alpha,A,ConstVectorView<Tx>(x.cptr(),x.size(),1,x.ct()),y);
            else 
                DoMultMV<add>(alpha,A,Vector<Tx>(x),y);
        } else if (y.step() == 0) {
            TMVAssert(y.size() <= 1);
            DoMultMV<add>(alpha,A,x,VectorView<T>(y.ptr(),y.size(),1,y.ct()));
#if 1
        } else if (y.step() != 1) {
            // Most BLAS implementations do fine with the y.step() < 0.
            // And in fact, they _should_ do ok even if the step is negative.
            // However, some implementations seem to propagate nan's from
            // the temporary memory they create to do the unit-1 calculation.
            // So to make sure they don't have to make a temporary, I just
            // do it here for them.
#else
        } else if (y.step() < 0) {
#endif
            Vector<T> yy(y.size());
            DoMultMV<false>(T(1),A,x,yy.view());
            if (add) y += alpha*yy;
            else y = alpha*yy;
#if 1
        } else if (x.step() != 1) {
            // I don't think the non-unit step is a problem, but just to be
            // sure...
#else
        } else if (x.step() < 0) {
#endif
            Vector<T> xx = alpha*x;
            DoMultMV<add>(T(1),A,xx,y);
        } else if (BlasIsCM(A) || BlasIsRM(A)) {
            if (!SameStorage(A,y)) {
                if (!SameStorage(x,y) && !SameStorage(A,x)) {
                    BlasMultMV(alpha,A,x,add?1:0,y);
                } else {
                    Vector<T> xx = alpha*x;
                    BlasMultMV(T(1),A,xx,add?1:0,y);
                }
            } else {
                Vector<T> yy(y.size());
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
        } else {
            if (TMV_IMAG(alpha) == T(0)) {
                Matrix<Ta,RowMajor> A2 = TMV_REAL(alpha)*A;
                DoMultMV<add>(T(1),A2,x,y);
            } else {
                Matrix<T,RowMajor> A2 = alpha*A;
                DoMultMV<add>(T(1),A2,x,y);
            }
        }
#else
        NonBlasMultMV<add>(alpha,A,x,y);
#endif
#ifdef XDEBUG
        std::cout<<"y => "<<y<<std::endl;
#endif
    }

    template <bool add, class T, class Ta, class Tx> void MultMV(
        const T alpha, const GenMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y)
        // y (+)= alpha * A * x
    { 
#ifdef XDEBUG
        cout<<"Start MultMV: alpha = "<<alpha<<endl;
        cout<<"add = "<<add<<endl;
        cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
        cout<<"x = "<<TMV_Text(x)<<" step "<<x.step()<<"  "<<x<<endl;
        cout<<"y = "<<TMV_Text(y)<<" step "<<y.step()<<endl;
        if (add) cout<<"y = "<<y<<endl;
        cout<<"ptrs = "<<A.cptr()<<"  "<<x.cptr()<<"  "<<y.cptr()<<std::endl;
        Vector<Tx> x0 = x;
        Vector<T> y0 = y;
        Matrix<Ta> A0 = A;
        Vector<T> y2 = y;
        for(ptrdiff_t i=0;i<y.size();i++) {
            if (add) y2(i) += alpha * (A.row(i) * x0);
            else y2(i) = alpha * (A.row(i) * x0);
        }
        cout<<"y2 = "<<y2<<endl;
#endif
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());

        if (y.size() > 0) {
            if (x.size()==0 || alpha==T(0)) {
                if (!add) y.setZero();
            } else if (y.isconj()) {
                DoMultMV<add>(
                    TMV_CONJ(alpha),A.conjugate(),x.conjugate(),y.conjugate());
            } else {
                DoMultMV<add>(alpha,A,x,y);
            }
        }

#ifdef XDEBUG
        cout<<"y => "<<y<<endl;
        if (!(Norm(y-y2) <= 
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(x0)+
                     (add?Norm(y0):TMV_RealType(T)(0))))) {
            cerr<<"MultMV: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A);
            if (A.rowsize() < 30 && A.colsize() < 30) cerr<<"  "<<A0;
            else cerr<<"  "<<A.colsize()<<" x "<<A.rowsize();
            cerr<<endl<<"x = "<<TMV_Text(x)<<" step "<<x.step();
            if (x.size() < 30) cerr<<"  "<<x0;
            cerr<<endl<<"y = "<<TMV_Text(y)<<" step "<<y.step();
            if (add && y.size() < 30) cerr<<"  "<<y0;
            cerr<<endl<<"Aptr = "<<A.cptr();
            cerr<<", xptr = "<<x.cptr()<<", yptr = "<<y.cptr()<<endl;
            if (y.size() < 200) {
                cerr<<"--> y = "<<y<<endl;
                cerr<<"y2 = "<<y2<<endl;
            } else {
                ptrdiff_t imax;
                (y-y2).maxAbsElement(&imax);
                cerr<<"y("<<imax<<") = "<<y(imax)<<endl;
                cerr<<"y2("<<imax<<") = "<<y2(imax)<<endl;
            }
            cerr<<"Norm(A0) = "<<Norm(A0)<<endl;
            cerr<<"Norm(x0) = "<<Norm(x0)<<endl;
            if (add) cerr<<"Norm(y0) = "<<Norm(y0)<<endl;
            cerr<<"|alpha|*|A0|*|x0|+?|y0| = "<<
                TMV_ABS(alpha)*Norm(A0)*Norm(x0)+
                (add?Norm(y0):TMV_RealType(T)(0))<<endl;
            cerr<<"Norm(y-y2) = "<<Norm(y-y2)<<endl;
            cerr<<"NormInf(y-y2) = "<<NormInf(y-y2)<<endl;
            cerr<<"Norm1(y-y2) = "<<Norm1(y-y2)<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_MultMV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


