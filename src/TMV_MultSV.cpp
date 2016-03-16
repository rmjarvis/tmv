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
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_MatrixArithFunc.h"
#ifdef BLAS
#include "tmv/TMV_SymMatrixArith.h"
#endif

#ifdef XDEBUG
#include <iostream>
#include "tmv/TMV_MatrixArith.h"
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    template <class T> 
    const T* SymMatrixComposite<T>::cptr() const
    {
        if (!itsm.get()) {
            ptrdiff_t s = this->size();
            ptrdiff_t len = s*s;
            itsm.resize(len);
            this->assignToS(SymMatrixView<T>(
                    itsm.get(),s,stepi(),stepj(),Sym,uplo(),NonConj 
                    TMV_FIRSTLAST1(itsm.get(),itsm.get()+len) ));
        }
        return itsm.get();
    }

    template <class T> 
    ptrdiff_t SymMatrixComposite<T>::stepi() const
    { return 1; }

    template <class T> 
    ptrdiff_t SymMatrixComposite<T>::stepj() const
    { return this->size(); }

    //
    // MultMV
    //

    template <bool add, class T, class Ta, class Tx> 
    static void DoUnitAMultMV(
        const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y)
    {
        const ptrdiff_t N = A.size();
        if (add) y += A.lowerTri() * x;
        else y = A.lowerTri() * x;

        if (N > 1)
            y.subVector(0,N-1) += A.upperTri().offDiag() * x.subVector(1,N);
    }

    template <bool add, class T, class Ta, class Tx> 
    static void UnitAMultMV(
        const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y)
    {
        // Check for 0's in the beginning or end of x:
        //      [ A11 A12 A13 ] [ 0 ]          [ A12 ]
        // y += [ A21 A22 A23 ] [ x ] --> y += [ A22 ] x
        //      [ A31 A32 A33 ] [ 0 ]          [ A32 ]

        const ptrdiff_t N = x.size(); // == A.size()
        ptrdiff_t j2 = N;
        for(const Tx* x2=x.cptr()+N-1; j2>0 && *x2==Tx(0); --j2,--x2);
        if (j2 == 0) {
            if (!add) y.setZero();
            return;
        }
        ptrdiff_t j1 = 0;
        for(const Tx* x1=x.cptr(); *x1==Tx(0); ++j1,++x1);
        if (j1 == 0 && j2 == N) DoUnitAMultMV<add>(A,x,y);
        else {
            if (j1 > 0)
                MultMV<add>(T(1),A.subMatrix(0,j1,j1,j2),x.subVector(j1,j2),
                            y.subVector(0,j1));
            TMVAssert(j1 != j2);
            DoUnitAMultMV<add>(A.subSymMatrix(j1,j2),x.subVector(j1,j2),
                               y.subVector(j1,j2));
            if (j2 < N)
                MultMV<add>(T(1),A.subMatrix(j2,N,j1,j2),x.subVector(j1,j2),
                            y.subVector(j2,N));
        }
    }

    template <bool add, class T, class Ta, class Tx> 
    static void NonBlasMultMV(
        const T alpha, const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y)
    // y (+)= alpha * A * x
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(A.size() == y.size());
        TMVAssert(alpha != T(0));
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));

        if (A.uplo() == Upper) 
            if (A.isherm()) NonBlasMultMV<add>(alpha,A.adjoint(),x,y);
            else NonBlasMultMV<add>(alpha,A.transpose(),x,y);
        else if (y.isconj())
            NonBlasMultMV<add>(TMV_CONJ(alpha),A.conjugate(),x.conjugate(),
                               y.conjugate());
        else {
            if (x.step() != 1) {
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
                    Vector<Tx> xx = TMV_REAL(alpha)*x;
                    if (y.step() != 1) {
                        Vector<T> yy(y.size());
                        UnitAMultMV<false>(A,xx,yy.view());
                        if (!add) y = yy;
                        else y += yy;
                    }
                    else 
                        UnitAMultMV<add>(A,xx,y);
                } else {
                    Vector<T> xx = alpha*x;
                    if (y.step()!=1) {
                        Vector<T> yy(y.size());
                        UnitAMultMV<false>(A,xx,yy.view());
                        if (add) y += yy;
                        else y = yy;
                    }
                    else
                        UnitAMultMV<add>(A,xx,y);
                }
            } else if (y.step()!=1 || alpha!=T(1)) {
                Vector<T> yy(y.size());
                UnitAMultMV<false>(A,x,yy.view());
                if (add) y += alpha * yy;
                else y = alpha * yy;
            } else {
                TMVAssert(alpha == T(1));
                TMVAssert(y.step() == 1);
                TMVAssert(x.step() == 1);
                UnitAMultMV<add>(A,x,y);
            }
        }
    }

#ifdef BLAS
    template <class T, class Ta, class Tx> 
    static inline void BlasMultMV(
        const T alpha, const GenSymMatrix<Ta>& A,
        const GenVector<Tx>& x, int beta, VectorView<T> y)
    { 
        if (beta==1) NonBlasMultMV<true>(alpha,A,x,y); 
        else NonBlasMultMV<false>(alpha,A,x,y); 
    }
#ifdef INST_DOUBLE
    template <> 
    void BlasMultMV(
        const double alpha,
        const GenSymMatrix<double>& A, const GenVector<double>& x,
        int beta, VectorView<double> y)
    {
        int n = A.size();
        int lda = A.stepj();
        int xs = x.step();
        int ys = y.step();
        const double* xp = x.cptr();
        if (xs < 0) xp += (n-1)*xs;
        double* yp = y.ptr();
        if (ys < 0) yp += (n-1)*ys;
        if (beta == 0) y.setZero();
        double xbeta(1);
        BLASNAME(dsymv) (
            BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
            BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASV(xbeta),
            BLASP(yp),BLASV(ys) BLAS1);
    }
    template <> 
    void BlasMultMV(
        const std::complex<double> alpha,
        const GenSymMatrix<std::complex<double> >& A,
        const GenVector<std::complex<double> >& x,
        int beta, VectorView<std::complex<double> > y)
    {
        if (A.isherm()) {
            int n = A.size();
            int lda = A.stepj();
            int xs = x.step();
            int ys = y.step();
            const std::complex<double>* xp = x.cptr();
            if (xs < 0) xp += (n-1)*xs;
            std::complex<double>* yp = y.ptr();
            if (ys < 0) yp += (n-1)*ys;
            if (beta == 0) y.setZero();
            std::complex<double> xbeta(1);
            BLASNAME(zhemv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
        } else {
#ifdef ELAP
            int n = A.size();
            int lda = A.stepj();
            int xs = x.step();
            int ys = y.step();
            const std::complex<double>* xp = x.cptr();
            if (xs < 0) xp += (n-1)*xs;
            std::complex<double>* yp = y.ptr();
            if (ys < 0) yp += (n-1)*ys;
            if (beta == 0) y.setZero();
            std::complex<double> xbeta(1);
            LAPNAME(zsymv) (
                LAPCM A.uplo()==Upper ? LAPCH_UP : LAPCH_LO,
                LAPV(n),LAPP(&alpha),LAPP(A.cptr()),LAPV(lda),
                LAPP(xp),LAPV(xs),LAPP(&xbeta),
                LAPP(yp),LAPV(ys) LAP1);
#else
            if (beta==1)
                NonBlasMultMV<true>(alpha,A,x,y);
            else
                NonBlasMultMV<false>(alpha,A,x,y);
#endif
        }
    }
    template <> 
    void BlasMultMV(
        const std::complex<double> alpha,
        const GenSymMatrix<std::complex<double> >& A,
        const GenVector<double>& x,
        int beta, VectorView<std::complex<double> > y)
    { BlasMultMV(alpha,A,Vector<std::complex<double> >(x),beta,y); }
    template <> 
    void BlasMultMV(
        const std::complex<double> alpha,
        const GenSymMatrix<double>& A,
        const GenVector<std::complex<double> >& x,
        int beta, VectorView<std::complex<double> > y)
    {
        if (beta == 0) {
            int n = A.size();
            int lda = A.stepj();
            int xs = 2*x.step();
            int ys = 2*y.step();
            const double* xp = (const double*) x.cptr();
            if (xs < 0) xp += (n-1)*xs;
            double* yp = (double*) y.ptr();
            if (ys < 0) yp += (n-1)*ys;
            double xalpha(1);
            y.setZero();
            double xbeta(1);
            BLASNAME(dsymv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(dsymv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp+1),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
            if (x.isconj()) y.conjugateSelf();
            y *= alpha;
        } else if (TMV_IMAG(alpha) == 0. && !x.isconj()) {
            int n = A.size();
            int lda = A.stepj();
            int xs = 2*x.step();
            int ys = 2*y.step();
            const double* xp = (const double*) x.cptr();
            if (xs < 0) xp += (n-1)*xs;
            double* yp = (double*) y.ptr();
            if (ys < 0) yp += (n-1)*ys;
            double xalpha(TMV_REAL(alpha));
            double xbeta(1);
            BLASNAME(dsymv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(dsymv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
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
        const GenSymMatrix<double>& A,
        const GenVector<double>& x,
        int beta, VectorView<std::complex<double> > y)
    {
        int n = A.size();
        int lda = A.stepj();
        int xs = x.step();
        int ys = 2*y.step();
        const double* xp = x.cptr();
        if (xs < 0) xp += (n-1)*xs;
        double* yp = (double*) y.ptr();
        if (ys < 0) yp += (n-1)*ys;
        double ar(TMV_REAL(alpha));
        double ai(TMV_IMAG(alpha));
        if (beta == 0) y.setZero();
        double xbeta(1);
        if (ar != 0.) {
            BLASNAME(dsymv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASV(ar),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
        }
        if (ai != 0.) {
            BLASNAME(dsymv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASV(ai),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void BlasMultMV(
        const float alpha,
        const GenSymMatrix<float>& A, const GenVector<float>& x,
        int beta, VectorView<float> y)
    {
        int n = A.size();
        int lda = A.stepj();
        int xs = x.step();
        int ys = y.step();
        const float* xp = x.cptr();
        if (xs < 0) xp += (n-1)*xs;
        float* yp = y.ptr();
        if (ys < 0) yp += (n-1)*ys;
        if (beta == 0) y.setZero();
        float xbeta(1);
        BLASNAME(ssymv) (
            BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
            BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASV(xbeta),
            BLASP(yp),BLASV(ys) BLAS1);
    }
    template <> 
    void BlasMultMV(
        const std::complex<float> alpha,
        const GenSymMatrix<std::complex<float> >& A,
        const GenVector<std::complex<float> >& x,
        int beta, VectorView<std::complex<float> > y)
    {
        if (A.isherm()) {
            int n = A.size();
            int lda = A.stepj();
            int xs = x.step();
            int ys = y.step();
            const std::complex<float>* xp = x.cptr();
            if (xs < 0) xp += (n-1)*xs;
            std::complex<float>* yp = y.ptr();
            if (ys < 0) yp += (n-1)*ys;
            if (beta == 0) y.setZero();
            std::complex<float> xbeta(1);
            BLASNAME(chemv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
        } else {
#ifdef ELAP
            int n = A.size();
            int lda = A.stepj();
            int xs = x.step();
            int ys = y.step();
            const std::complex<float>* xp = x.cptr();
            if (xs < 0) xp += (n-1)*xs;
            std::complex<float>* yp = y.ptr();
            if (ys < 0) yp += (n-1)*ys;
            if (beta == 0) y.setZero();
            std::complex<float> xbeta(1);
            LAPNAME(csymv) (
                LAPCM A.uplo()==Upper ? LAPCH_UP : LAPCH_LO,
                LAPV(n),LAPP(&alpha),LAPP(A.cptr()),LAPV(lda),
                LAPP(xp),LAPV(xs),LAPP(&xbeta),
                LAPP(yp),LAPV(ys) LAP1);
#else
            if (beta==1)
                NonBlasMultMV<true>(alpha,A,x,y);
            else
                NonBlasMultMV<false>(alpha,A,x,y);
#endif
        }
    }
    template <> 
    void BlasMultMV(
        const std::complex<float> alpha,
        const GenSymMatrix<std::complex<float> >& A,
        const GenVector<float>& x,
        int beta, VectorView<std::complex<float> > y)
    { BlasMultMV(alpha,A,Vector<std::complex<float> >(x),beta,y); }
    template <> 
    void BlasMultMV(
        const std::complex<float> alpha,
        const GenSymMatrix<float>& A,
        const GenVector<std::complex<float> >& x,
        int beta, VectorView<std::complex<float> > y)
    {
        if (beta == 0) {
            int n = A.size();
            int lda = A.stepj();
            int xs = 2*x.step();
            int ys = 2*y.step();
            const float* xp = (const float*) x.cptr();
            if (xs < 0) xp += (n-1)*xs;
            float* yp = (float*) y.ptr();
            if (ys < 0) yp += (n-1)*ys;
            float xalpha(1);
            y.setZero();
            float xbeta(1);
            BLASNAME(ssymv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(ssymv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp+1),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
            if (x.isconj()) y.conjugateSelf();
            y *= alpha;
        } else if (TMV_IMAG(alpha) == 0.F && !x.isconj()) {
            int n = A.size();
            int lda = A.stepj();
            int xs = 2*x.step();
            int ys = 2*y.step();
            const float* xp = (const float*) x.cptr();
            if (xs < 0) xp += (n-1)*xs;
            float* yp = (float*) y.ptr();
            if (ys < 0) yp += (n-1)*ys;
            float xalpha(TMV_REAL(alpha));
            float xbeta(1);
            BLASNAME(ssymv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(ssymv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
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
        const GenSymMatrix<float>& A,
        const GenVector<float>& x,
        int beta, VectorView<std::complex<float> > y)
    {
        int n = A.size();
        int lda = A.stepj();
        int xs = x.step();
        int ys = 2*y.step();
        const float* xp = x.cptr();
        if (xs < 0) xp += (n-1)*xs;
        float* yp = (float*) y.ptr();
        if (ys < 0) yp += (n-1)*ys;
        float ar(TMV_REAL(alpha));
        float ai(TMV_IMAG(alpha));
        if (beta == 0) y.setZero();
        float xbeta(1);
        if (ar != 0.F) {
            BLASNAME(ssymv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASV(ar),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
        }
        if (ai != 0.F) {
            BLASNAME(ssymv) (
                BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
                BLASV(n),BLASV(ai),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        }
    }
#endif 
#endif // BLAS

    template <bool add, class T, class Ta, class Tx> 
    static void DoMultMV(
        const T alpha, const GenSymMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<T> y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(alpha != T(0));
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));

#ifdef BLAS
        if (A.isconj()) 
            DoMultMV<add>(
                TMV_CONJ(alpha),A.conjugate(),x.conjugate(),y.conjugate());
        else if (!A.iscm() && A.isrm()) 
            if (A.isherm()) DoMultMV<add>(alpha,A.adjoint(),x,y);
            else DoMultMV<add>(alpha,A.transpose(),x,y);
        else if (x.step() == 0) {
            if (x.size() <= 1) 
                DoMultMV<add>(
                    alpha,A,ConstVectorView<Tx>(x.cptr(),x.size(),1,x.ct()),y);
            else 
                DoMultMV<add>(alpha,A,Vector<Tx>(x),y);
        } else if (y.step() == 0) {
            TMVAssert(y.size() <= 1);
            DoMultMV<add>(alpha,A,x,VectorView<T>(y.ptr(),y.size(),1,y.ct()));
        } else if (A.iscm()&&A.stepj()>0) {
            if (!y.isconj() && y.step() != 1) { 
                if (!x.isconj() && x.step() != 1) {
                    BlasMultMV(alpha,A,x,add?1:0,y);
                } else {
                    Vector<T> xx = alpha*x;
                    BlasMultMV(T(1),A,xx,add?1:0,y);
                }
            } else {
                Vector<T> yy(y.size());
                if (!x.isconj() && x.step() != 1) {
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
        } else if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) {
            if (A.isherm()) {
                if (A.uplo() == Upper) {
                    HermMatrix<Ta,Upper|ColMajor> A2 =
                        TMV_REAL(alpha)*A;
                    DoMultMV<add>(T(1),A2,x,y);
                } else {
                    HermMatrix<Ta,Lower|ColMajor> A2 =
                        TMV_REAL(alpha)*A;
                    DoMultMV<add>(T(1),A2,x,y);
                }
            } else {
                if (A.uplo() == Upper) {
                    SymMatrix<Ta,Upper|ColMajor> A2 =
                        TMV_REAL(alpha)*A;
                    DoMultMV<add>(T(1),A2,x,y);
                } else {
                    SymMatrix<Ta,Lower|ColMajor> A2 =
                        TMV_REAL(alpha)*A;
                    DoMultMV<add>(T(1),A2,x,y);
                }
            }
        } else {
            if (!A.issym()) {
                if (A.uplo() == Upper) {
                    HermMatrix<Ta,Upper|ColMajor> A2 = A;
                    DoMultMV<add>(alpha,A2,x,y);
                } else {
                    HermMatrix<Ta,Lower|ColMajor> A2 = A;
                    DoMultMV<add>(alpha,A2,x,y);
                }
            } else {
                if (A.uplo() == Upper) {
                    SymMatrix<T,Upper|ColMajor> A2 = alpha*A;
                    DoMultMV<add>(T(1),A2,x,y);
                } else {
                    SymMatrix<T,Lower|ColMajor> A2 = alpha*A;
                    DoMultMV<add>(T(1),A2,x,y);
                }
            }
        }
#else
        NonBlasMultMV<add>(alpha,A,x,y);
#endif
    }

    template <bool add, class T, class Ta, class Tx> 
    void MultMV(
        const T alpha, const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y)
    // y (+)= alpha * A * x
    { 
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
#ifdef XDEBUG
        cout<<"Start MultMV: alpha = "<<alpha<<endl;
        cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
        cout<<"x = "<<TMV_Text(x)<<" step "<<x.step()<<"  "<<x<<endl;
        cout<<"y = "<<TMV_Text(y)<<" step "<<y.step()<<"  "<<y<<endl;
        Matrix<Ta> A0 = A;
        Vector<Tx> x0 = x;
        Vector<T> y0 = y;
        Vector<T> y2 = alpha*A0*x0;
        if (add) y2 += y0;
#endif

        if (y.size() > 0) {
            if (x.size()==0 || alpha==T(0)) {
                if (!add) y.setZero();
            } else if (SameStorage(x,y)) {
                Vector<T> yy(y.size());
                DoMultMV<false>(T(1),A,x,yy.view());
                if (add) y += alpha*yy;
                else y = alpha*yy;
            } else {
                DoMultMV<add>(alpha,A,x,y);
            } 
        }
#ifdef XDEBUG
        cout<<"--> y = "<<y<<endl;
        cout<<"y2 = "<<y2<<endl;
        if (!(Norm(y-y2) <=
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(x0)+
                     (add?Norm(y0):TMV_RealType(T)(0))))) {
            cerr<<"MultMV: alpha = "<<alpha<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"x = "<<TMV_Text(x)<<" step "<<x.step()<<"  "<<x0<<endl;
            cerr<<"y = "<<TMV_Text(y)<<" step "<<y.step()<<"  "<<y0<<endl;
            cerr<<"--> y = "<<y<<endl;
            cerr<<"y2 = "<<y2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_MultSV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


