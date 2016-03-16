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
#include "tmv/TMV_SymMatrixArithFunc.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    // 
    // Rank2Update
    //

    template <bool ha, bool add, class T, class Tx, class Ty> 
    static void UpperRank2Update(
        const GenVector<Tx>& x, const GenVector<Ty>& y,
        SymMatrixView<T> A)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(A.size() == y.size());
        TMVAssert(A.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.uplo() == Upper);
        TMVAssert(x.step()==1);
        TMVAssert(y.step()==1);
        TMVAssert(A.iscm());
        TMVAssert(x.ct() == NonConj);
        TMVAssert(y.ct() == NonConj);
        TMVAssert(ha == !A.issym());
#ifdef XDEBUG
        TMV_RealType(T) normA = Norm(A);
        TMV_RealType(T) normx = Norm(x);
        TMV_RealType(T) normy = Norm(y);
        TMV_RealType(T) eps = TMV_RealType(T)(2)*normx*normy;
        if (add) eps += normA;
        eps *= A.size() * TMV_Epsilon<T>();
#endif

        //std::cout<<"Start Upper"<<std::endl;
        //std::cout<<"x = "<<x<<std::endl;
        //std::cout<<"x.step = "<<x.step()<<std::endl;
        //std::cout<<"y = "<<x<<std::endl;
        //std::cout<<"y.step = "<<y.step()<<std::endl;
        //std::cout<<"A = "<<x<<std::endl;
        //std::cout<<"A.step = "<<A.stepi()<<" "<<A.stepj()<<std::endl;

        const ptrdiff_t sj = A.stepj();
        const ptrdiff_t N = A.size();
        const Tx*const x0 = x.cptr();
        const Ty*const y0 = y.cptr();

        const Tx* xj = x0;
        const Ty* yj = y0;

        T A00;
        if (*xj == Tx(0) || *yj == Ty(0)) A00 = T(0);
        else if (ha) A00 = TMV_RealType(T)(2) * TMV_REAL(*xj * TMV_CONJ(*yj));
        else A00 = TMV_RealType(T)(2) * (*xj) * (*yj);
        ++xj; ++yj;
        T* Acolj = A.ptr()+sj;

        for (ptrdiff_t j=1;j<N;++j,++xj,++yj,Acolj+=sj) {
            // A.col(j,0,j+1) += ax * y.subVector(0,j+1) + ay * x.subVector(0,j+1);
            if (*xj != Tx(0)) {
                T* Aij = Acolj;
                const Ty* yi = y0;
                if (*yj != Ty(0)) {
                    const Tx* xi = x0;
                    for(ptrdiff_t i=j+1;i>0;--i,++yi,++xi,++Aij) {
                        T temp = *xi * (ha ? TMV_CONJ(*yj) : *yj);
                        temp += *yi * (ha ? TMV_CONJ(*xj) : *xj);
#ifdef TMVFLDEBUG
                        TMVAssert(Aij >= A._first);
                        TMVAssert(Aij < A._last);
#endif
                        if (add) *Aij += temp;
                        else *Aij = temp;
                    }
                } else {
                    for(ptrdiff_t i=j+1;i>0;--i,++yi,++Aij) {
                        const T temp = *yi * (ha ? TMV_CONJ(*xj) : *xj);
#ifdef TMVFLDEBUG
                        TMVAssert(Aij >= A._first);
                        TMVAssert(Aij < A._last);
#endif
                        if (add) *Aij += temp;
                        else *Aij = temp;
                    }
                }
            } else if (*yj != Ty(0)) {
                T* Aij = Acolj;
                const Tx* xi = x0;
                for(ptrdiff_t i=j+1;i>0;--i,++xi,++Aij) {
                    const T temp = *xi * (ha ? TMV_CONJ(*yj) : *yj);
#ifdef TMVFLDEBUG
                    TMVAssert(Aij >= A._first);
                    TMVAssert(Aij < A._last);
#endif
                    if (add) *Aij += temp;
                    else *Aij = temp;
                }
            } else if (!add) {
                T* Aij = Acolj;
                std::fill_n(Aij,j+1,T(0));
            }
        }
#ifdef TMVFLDEBUG
        TMVAssert(A.ptr() >= A._first);
        TMVAssert(A.ptr() < A._last);
#endif
        if (add) *A.ptr() += A00;
        else *A.ptr() = A00;
        if (ha && isComplex(T())) {
#ifdef XDEBUG
            TMVAssert(NormInf(A.diag().imagPart()) <= eps);
#endif
            A.diag().imagPart().setZero();
        }
        //std::cout<<"Done Upper"<<std::endl;
    }

    template <bool ha, bool add, class T, class Tx, class Ty> 
    static void LowerRank2Update(
        const GenVector<Tx>& x, const GenVector<Ty>& y,
        SymMatrixView<T> A)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(A.size() == y.size());
        TMVAssert(A.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.iscm());
        TMVAssert(x.step()==1);
        TMVAssert(y.step()==1);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(y.ct() == NonConj);
        TMVAssert(ha == !A.issym());
#ifdef XDEBUG
        TMV_RealType(T) normA = Norm(A);
        TMV_RealType(T) normx = Norm(x);
        TMV_RealType(T) normy = Norm(y);
        TMV_RealType(T) eps = TMV_RealType(T)(2)*normx*normy;
        if (add) eps += normA;
        eps *= A.size() * TMV_Epsilon<T>();
#endif

        //std::cout<<"Start Lower"<<std::endl;
        const ptrdiff_t ds = A.stepj() + 1;
        const ptrdiff_t N = A.size();
        const Tx* xj = x.cptr()+N-1;
        const Ty* yj = y.cptr()+N-1;
        T* Ajj = A.ptr()+(N-1)*ds;

        for (ptrdiff_t jj=N,Nmj=1;jj>0;--jj,++Nmj,--xj,--yj,Ajj-=ds) {
            // Nmj = N-j
            // A.col(j,j,N) += (A.isherm() ? TMV_CONJ(*yj) : *yj) * x.subVector(j,N);
            // A.col(j,j,N) += (A.isherm() ? TMV_CONJ(*xj) : *xj) * y.subVector(j,N);
            if (*yj!=Ty(0)) {
                T* Aij = Ajj;
                const Tx* xi = xj;
                if (*xj!=Tx(0)) {
                    const Ty* yi = yj;
                    for(ptrdiff_t i=Nmj;i>0;--i,++xi,++yi,++Aij) {
                        T temp = *xi * (ha ? TMV_CONJ(*yj) : *yj);
                        temp += *yi * (ha ? TMV_CONJ(*xj) : *xj);
#ifdef TMVFLDEBUG
                        TMVAssert(Aij >= A._first);
                        TMVAssert(Aij < A._last);
#endif
                        if (add) *Aij += temp;
                        else *Aij = temp;
                    }
                } else {
                    for(ptrdiff_t i=Nmj;i>0;--i,++xi,++Aij) {
                        const T temp = *xi * (ha ? TMV_CONJ(*yj) : *yj);
#ifdef TMVFLDEBUG
                        TMVAssert(Aij >= A._first);
                        TMVAssert(Aij < A._last);
#endif
                        if (add) *Aij += temp;
                        else *Aij = temp;
                    }
                }
            } else if (*xj!=Tx(0)) {
                T* Aij = Ajj;
                const Ty* yi = yj;
                for(ptrdiff_t i=Nmj;i>0;--i,++yi,++Aij) {
                    const T temp = *yi * (ha ? TMV_CONJ(*xj) : *xj);
#ifdef TMVFLDEBUG
                    TMVAssert(Aij >= A._first);
                    TMVAssert(Aij < A._last);
#endif
                    if (add) *Aij += temp;
                    else *Aij = temp;
                }
            } else if (!add) {
                T* Aij = Ajj;
                std::fill_n(Aij,Nmj,T(0));
            }
        }
        if (ha && isComplex(T())) {
#ifdef XDEBUG
            TMVAssert(NormInf(A.diag().imagPart()) <= eps);
#endif
            A.diag().imagPart().setZero();
        }
        //std::cout<<"Done Lower"<<std::endl;
    }

    template <bool add, class T, class Tx, class Ty> 
    struct UnitARank2Update
    {
        static void F(
            const GenVector<Tx>& x, const GenVector<Ty>& y, SymMatrixView<T> A)
        {
            TMVAssert(A.size() == x.size());
            TMVAssert(A.size() == y.size());
            TMVAssert(A.size() > 0);
            TMVAssert(A.ct() == NonConj);
            TMVAssert(A.iscm());
            TMVAssert(x.step() == 1);
            TMVAssert(y.step() == 1);
            TMVAssert(x.ct() == NonConj);
            TMVAssert(y.ct() == NonConj);

            if (A.isupper()) UpperRank2Update<false,add>(x,y,A);
            else LowerRank2Update<false,add>(x,y,A);
        }
    };

    template <bool add, class T, class Tx, class Ty> 
    struct UnitARank2Update<add,std::complex<T>,Tx,Ty>
    {
        static void F(
            const GenVector<Tx>& x, const GenVector<Ty>& y,
            SymMatrixView<std::complex<T> > A)
        {
            TMVAssert(A.size() == x.size());
            TMVAssert(A.size() == y.size());
            TMVAssert(A.size() > 0);
            TMVAssert(A.ct() == NonConj);
            TMVAssert(A.iscm());
            TMVAssert(x.step() == 1);
            TMVAssert(y.step() == 1);
            TMVAssert(x.ct() == NonConj);
            TMVAssert(y.ct() == NonConj);

            if (A.isherm())
                if (A.isupper()) UpperRank2Update<true,add>(x,y,A);
                else LowerRank2Update<true,add>(x,y,A);
            else
                if (A.isupper()) UpperRank2Update<false,add>(x,y,A);
                else LowerRank2Update<false,add>(x,y,A);
        }
    };


    template <bool add, class T, class Tx, class Ty> 
    static void NonBlasRank2Update(
        const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
        SymMatrixView<T> A)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(A.size() == y.size());
        TMVAssert(A.iscm());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(alpha != T(0));
        TMVAssert(A.size() > 0);
        //std::cout<<"start NonBlas Rank2Update"<<std::endl;

        if (x.step() != 1 || x.isconj()) {
            // Copy x to new storage
            if (y.step() != 1 || y.isconj()) {
                // Copy x and y to new storage
                if (x.size() <= y.size()) {
                    if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                        Vector<Tx> xx = TMV_REAL(alpha)*x;
                        Vector<Ty> yy = y;
                        UnitARank2Update<add,T,Tx,Ty>::F(xx,yy,A);
                    } else {
                        Vector<T> xx = alpha*x;
                        Vector<Ty> yy = y;
                        UnitARank2Update<add,T,T,Ty>::F(xx,yy,A);
                    }
                } else {
                    if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                        Vector<Tx> xx = x;
                        Vector<Ty> yy = TMV_REAL(alpha)*y;
                        UnitARank2Update<add,T,Tx,Ty>::F(xx,yy,A);
                    } else {
                        Vector<Tx> xx = x;
                        Vector<T> yy = alpha*y;
                        UnitARank2Update<add,T,Tx,T>::F(xx,yy,A);
                    }
                }
            } else {
                // Copy only x to new storage
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                    Vector<Tx> xx = TMV_REAL(alpha)*x;
                    UnitARank2Update<add,T,Tx,Ty>::F(xx,y,A);
                } else {
                    Vector<T> xx = alpha*x;
                    UnitARank2Update<add,T,T,Ty>::F(xx,y,A);
                }
            }
        } else if (y.step() != 1 || y.isconj()) {
            // Copy only y to new storage
            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                Vector<Ty> yy = TMV_REAL(alpha)*y;
                UnitARank2Update<add,T,Tx,Ty>::F(x,yy,A);
            } else {
                Vector<T> yy = alpha*y;
                UnitARank2Update<add,T,Tx,T>::F(x,yy,A);
            }
        } else if (alpha != T(1)) {
            // Copy something to new storage to incorporate alpha
            if (x.size() <= y.size()) {
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                    Vector<Tx> xx = TMV_REAL(alpha)*x;
                    UnitARank2Update<add,T,Tx,Ty>::F(xx,y,A);
                } else {
                    Vector<T> xx = alpha*x;
                    UnitARank2Update<add,T,T,Ty>::F(xx,y,A);
                }
            } else {
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                    Vector<Ty> yy = TMV_REAL(alpha)*y;
                    UnitARank2Update<add,T,Tx,Ty>::F(x,yy,A);
                } else {
                    Vector<T> yy = alpha*y;
                    UnitARank2Update<add,T,Tx,T>::F(x,yy,A);
                }
            }
        } else {
            UnitARank2Update<add,T,Tx,Ty>::F(x,y,A);
        }
        //std::cout<<"Done NonBlas"<<std::endl;
    }

#ifdef BLAS
    template <class T, class Tx, class Ty> 
    static inline void BlasRank2Update(
        const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
        SymMatrixView<T> A)
    { NonBlasRank2Update<true>(alpha,x,y,A); }
#ifdef INST_DOUBLE
    template <> 
    void BlasRank2Update(
        const double alpha, const GenVector<double>& x,
        const GenVector<double>& y, SymMatrixView<double> A)
    {
        int n=A.size();
        int xs=x.step();
        int ys=y.step();
        const double* xp = x.cptr();
        if (xs < 0) xp += (n-1)*xs;
        const double* yp = y.cptr();
        if (ys < 0) yp += (n-1)*ys;
        int lda=A.stepj();
        BLASNAME(dsyr2) (
            BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
            BLASV(n),BLASV(alpha),BLASP(xp),BLASV(xs),
            BLASP(yp),BLASV(ys),BLASP(A.ptr()),BLASV(lda) BLAS1);
    }
    template <> 
    void BlasRank2Update(
        const std::complex<double> alpha,
        const GenVector<std::complex<double> >& x, 
        const GenVector<std::complex<double> >& y, 
        SymMatrixView<std::complex<double> > A)
    {
        if (A.issym() && (x.step() != 1 || y.step() != 1)) {
            if (x.step() != 1) {
                Vector<std::complex<double> > xx = x;
                if (y.step() != 1) {
                    Vector<std::complex<double> > yy = y;
                    return BlasRank2Update(alpha,xx,yy,A);
                } 
                else return BlasRank2Update(alpha,xx,y,A);
            } else {
                Vector<std::complex<double> > yy = y;
                return BlasRank2Update(alpha,x,yy,A);
            }
        } else {
            int n=A.size();
            int xs=x.step();
            int ys=y.step();
            const std::complex<double>* xp = x.cptr();
            if (xs < 0) xp += (n-1)*xs;
            const std::complex<double>* yp = y.cptr();
            if (ys < 0) yp += (n-1)*ys;
            int lda=A.stepj();
            if (A.isherm()) {
                BLASNAME(zher2) (
                    BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
                    BLASV(n),BLASP(&alpha),BLASP(xp),BLASV(xs),
                    BLASP(yp),BLASV(ys),BLASP(A.ptr()),BLASV(lda) BLAS1);
            } else {
                int k=1;
                std::complex<double> beta(1);
                BLASNAME(zsyr2k) (
                    BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
                    BLASCH_NT,BLASV(n),BLASV(k),BLASP(&alpha),
                    BLASP(xp),BLASV(n),BLASP(yp),BLASV(n),
                    BLASP(&beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
            }
        }
    }
    template <> 
    void BlasRank2Update(
        const std::complex<double> alpha,
        const GenVector<std::complex<double> >& x, 
        const GenVector<double>& y, 
        SymMatrixView<std::complex<double> > A)
    {
        SymMatrix<double> A1(A.size(),0.);
        BlasRank2Update(1.,x.realPart(),y,A1.view());
        A += alpha*A1;
        A1.setZero();
        BlasRank2Update(1.,x.imagPart(),y,A1.view());
        A += std::complex<double>(0,1)*alpha*A1;
    }
    template <> 
    void BlasRank2Update(
        const std::complex<double> alpha,
        const GenVector<double>& x, 
        const GenVector<std::complex<double> >& y, 
        SymMatrixView<std::complex<double> > A)
    {
        SymMatrix<double> A1(A.size(),0.);
        BlasRank2Update(1.,x,y.realPart(),A1.view());
        A += alpha*A1;
        A1.setZero();
        BlasRank2Update(1.,x,y.imagPart(),A1.view());
        A += std::complex<double>(0,1)*alpha*A1;
    }
    template <> 
    void BlasRank2Update(
        const std::complex<double> alpha,
        const GenVector<double>& x, const GenVector<double>& y, 
        SymMatrixView<std::complex<double> > A)
    {
        SymMatrix<double> A1(A.size(),0.);
        BlasRank2Update(1.,x,y,A1.view());
        A += alpha*A1;
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void BlasRank2Update(
        const float alpha, const GenVector<float>& x,
        const GenVector<float>& y, SymMatrixView<float> A)
    {
        int n=A.size();
        int xs=x.step();
        int ys=y.step();
        const float* xp = x.cptr();
        if (xs < 0) xp += (n-1)*xs;
        const float* yp = y.cptr();
        if (ys < 0) yp += (n-1)*ys;
        int lda=A.stepj();
        BLASNAME(ssyr2) (
            BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
            BLASV(n),BLASV(alpha),BLASP(xp),BLASV(xs),
            BLASP(yp),BLASV(ys),BLASP(A.ptr()),BLASV(lda) BLAS1);
    }
    template <> 
    void BlasRank2Update(
        const std::complex<float> alpha,
        const GenVector<std::complex<float> >& x, 
        const GenVector<std::complex<float> >& y, 
        SymMatrixView<std::complex<float> > A)
    {
        if (A.issym() && (x.step() != 1 || y.step() != 1)) {
            if (x.step() != 1) {
                Vector<std::complex<float> > xx = x;
                if (y.step() != 1) {
                    Vector<std::complex<float> > yy = y;
                    return BlasRank2Update(alpha,xx,yy,A);
                } 
                else return BlasRank2Update(alpha,xx,y,A);
            } else {
                Vector<std::complex<float> > yy = y;
                return BlasRank2Update(alpha,x,yy,A);
            }
        } else {
            int n=A.size();
            int xs=x.step();
            int ys=y.step();
            const std::complex<float>* xp = x.cptr();
            if (xs < 0) xp += (n-1)*xs;
            const std::complex<float>* yp = y.cptr();
            if (ys < 0) yp += (n-1)*ys;
            int lda=A.stepj();
            if (A.isherm()) {
                BLASNAME(cher2) (
                    BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
                    BLASV(n),BLASP(&alpha),BLASP(xp),BLASV(xs),
                    BLASP(yp),BLASV(ys),BLASP(A.ptr()),BLASV(lda) BLAS1);
            } else {
                int k=1;
                std::complex<float> beta(1);
                BLASNAME(csyr2k) (
                    BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
                    BLASCH_NT,BLASV(n),BLASV(k),BLASP(&alpha),
                    BLASP(xp),BLASV(n),BLASP(yp),BLASV(n),
                    BLASP(&beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
            }
        }
    }
    template <> 
    void BlasRank2Update(
        const std::complex<float> alpha,
        const GenVector<std::complex<float> >& x, 
        const GenVector<float>& y, 
        SymMatrixView<std::complex<float> > A)
    {
        SymMatrix<float> A1(A.size(),0.F);
        BlasRank2Update(1.F,x.realPart(),y,A1.view());
        A += alpha*A1;
        A1.setZero();
        BlasRank2Update(1.F,x.imagPart(),y,A1.view());
        A += std::complex<float>(0,1)*alpha*A1;
    }
    template <> 
    void BlasRank2Update(
        const std::complex<float> alpha,
        const GenVector<float>& x, 
        const GenVector<std::complex<float> >& y, 
        SymMatrixView<std::complex<float> > A)
    {
        SymMatrix<float> A1(A.size(),0.F);
        BlasRank2Update(1.F,x,y.realPart(),A1.view());
        A += alpha*A1;
        A1.setZero();
        BlasRank2Update(1.F,x,y.imagPart(),A1.view());
        A += std::complex<float>(0,1)*alpha*A1;
    }
    template <> 
    void BlasRank2Update(
        const std::complex<float> alpha,
        const GenVector<float>& x, const GenVector<float>& y, 
        SymMatrixView<std::complex<float> > A)
    {
        SymMatrix<float> A1(A.size(),0.F);
        BlasRank2Update(1.F,x,y,A1.view());
        A += alpha*A1;
    }
#endif 
#endif // BLAS

    template <bool add, class T, class Tx, class Ty> 
    void Rank2Update(
        const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
        SymMatrixView<T> A)
    // if A is sym:  A (+)= alpha * (x ^ y + y ^ x)
    // if A is herm: A (+)= alpha * x ^ y* + conj(alpha) * y ^ x*
    {
#ifdef XDEBUG
        Vector<Tx> x0 = x;
        Vector<Ty> y0 = y;
        Matrix<T> A0 = A;
        Matrix<T> A2 = A;
        if (A.isherm()) {
            if (add) A2 += (alpha*x^y.conjugate());
            else A2 = (alpha*x^y.conjugate());
            A2 += (TMV_CONJ(alpha)*y^x.conjugate());
        }
        else {
            if (add) A2 += alpha*(x^y);
            else A2 = alpha*(x^y);
            A2 += alpha*(y^x);
        }
        //cout<<"Start Rank2Update: alpha = "<<alpha<<endl;
        //cout<<"add = "<<add<<endl;
        //cout<<"x = "<<TMV_Text(x)<<"  step = "<<x.step()<<"  "<<x0<<endl;
        //cout<<"y = "<<TMV_Text(y)<<"  step = "<<y.step()<<"  "<<y0<<endl;
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
#endif

        TMVAssert(A.size() == x.size());
        TMVAssert(A.size() == y.size());
        if (alpha != T(0) && A.size() > 0) {
            if (A.isconj()) 
                Rank2Update<add>(
                    TMV_CONJ(alpha),x.conjugate(),y.conjugate(),A.conjugate());
            else if (!A.iscm() && A.isrm())
                if (A.isherm()) Rank2Update<add>(alpha,x,y,A.adjoint());
                else Rank2Update<add>(alpha,x,y,A.transpose());
            else if (!(A.iscm() 
#ifdef BLAS
                       && A.stepj()>0
#endif
            )) {
                if (A.isherm()) {
                    HermMatrix<T,Lower|ColMajor> AA(A.size());
                    Rank2Update<false>(alpha,x,y,AA.view());
                    if (add) A += AA;
                    else A = AA;
                } else {
                    SymMatrix<T,Lower|ColMajor> AA(A.size());
                    Rank2Update<false>(alpha,x,y,AA.view());
                    if (add) A += AA;
                    else A = AA;
                }
            } else {
#ifdef BLAS
                if (x.isconj() || x.step()!=1 || SameStorage(x,A)) {
                    if (y.isconj() || y.step()!=1 || SameStorage(y,A)) {
                        if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                            Vector<Tx> xx = TMV_REAL(alpha)*x;
                            Vector<Ty> yy = y;
                            if (!add) A.setZero();
                            BlasRank2Update(T(1),xx,yy,A);
                        } else {
                            Vector<T> xx = alpha*x;
                            Vector<Ty> yy = y;
                            if (!add) A.setZero();
                            BlasRank2Update(T(1),xx,yy,A);
                        }
                    } else {
                        if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                            Vector<Tx> xx = TMV_REAL(alpha)*x;
                            if (!add) A.setZero();
                            BlasRank2Update(T(1),xx,y,A);
                        } else {
                            Vector<T> xx = alpha*x;
                            if (!add) A.setZero();
                            BlasRank2Update(T(1),xx,y,A);
                        }
                    }
                } else {
                    if (y.isconj() || y.step()!=1 || SameStorage(y,A)) {
                        if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                            Vector<Ty> yy = TMV_REAL(alpha)*y;
                            if (!add) A.setZero();
                            BlasRank2Update(T(1),x,yy,A);
                        } else {
                            Vector<T> yy = TMV_CONJ(alpha)*y;
                            if (!add) A.setZero();
                            BlasRank2Update(T(1),x,yy,A);
                        }
                    } else {
                        if (!add) A.setZero();
                        BlasRank2Update(alpha,x,y,A);
                    }
                }
#else
                NonBlasRank2Update<add>(alpha,x,y,A);
#endif
            }
        }

#ifdef XDEBUG
        TMVAssert(A.isHermOK());
        //cout<<"Done Rank2Update: \n";
        //cout<<"Norm(A-A2) = "<<Norm(A-A2)<<endl;
        if (!(Norm(A-A2) < 0.001*(TMV_ABS(alpha)*Norm(x0)*Norm(y0)+Norm(A0)))) {
            cerr<<"Rank2Update: alpha = "<<alpha<<endl;
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

#define InstFile "TMV_Rank2_VVS.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


