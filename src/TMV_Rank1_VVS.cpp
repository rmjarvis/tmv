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
#ifdef BLAS
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_VectorArith.h"
#endif

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_SymMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    // 
    // Rank1Update
    //

    template <bool a1, bool cx, bool ha, bool rm, bool add, class T, class Tx, class Ta> 
    static void RowRank1Update(
        const Ta alpha, const GenVector<Tx>& x, SymMatrixView<T> A)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.uplo() == Lower);
        TMVAssert(a1 == (alpha == Ta(1)));
        TMVAssert(x.step()==1);
        TMVAssert(cx == x.isconj());
        TMVAssert(ha == A.isherm());
#ifdef XDEBUG
        TMV_RealType(T) normA = Norm(A);
        TMV_RealType(T) normx = Norm(x);
        TMV_RealType(T) eps = TMV_ABS(alpha)*normx*normx;
        if (add) eps += normA;
        eps *= A.size() * TMV_Epsilon<T>();
#endif

        const ptrdiff_t si = A.stepi();
        const ptrdiff_t sj = A.stepj();
        const ptrdiff_t N = A.size();
        const Tx* xi = x.cptr();
        const Tx*const x0 = x.cptr();

        T A00;
        if (*xi == Tx(0)) A00 = T(0);
        else if (ha) A00 = alpha * TMV_NORM(*xi);
        else if (cx) A00 = alpha * TMV_CONJ(*xi * *xi);
        else A00 = alpha * (*xi * *xi);
        ++xi;
        T* Arowi = A.ptr()+si;

        for (ptrdiff_t i=1;i<N;++i,++xi,Arowi+=si) {
            if (*xi != Tx(0)) {
                // A.row(i,0,i+1) += ax * x.subVector(0,i+1);
                T* Aij = Arowi;
                const Tx* xj = x0;

                if (a1 || isReal(alpha)) {
                    Tx axi = cx ? TMV_CONJ(*xi) : *xi;
                    if (!a1) axi *= TMV_REAL(alpha);
                    for(ptrdiff_t j=i+1;j>0;--j,++xj,(rm?++Aij:Aij+=sj)) {
                        const T temp = axi * (ha==cx ? *xj : TMV_CONJ(*xj));
#ifdef TMVFLDEBUG
                        TMVAssert(Aij >= A._first);
                        TMVAssert(Aij < A._last);
#endif
                        if (add) *Aij += temp;
                        else *Aij = temp;
                    }
                } else {
                    T axi = alpha * (cx ? TMV_CONJ(*xi) : *xi);
                    for(ptrdiff_t j=i+1;j>0;--j,++xj,(rm?++Aij:Aij+=sj)) {
                        const T temp = axi * (ha==cx ? *xj : TMV_CONJ(*xj));
#ifdef TMVFLDEBUG
                        TMVAssert(Aij >= A._first);
                        TMVAssert(Aij < A._last);
#endif
                        if (add) *Aij += temp;
                        else *Aij = temp;
                    }
                }
            } else if (!add) {
                T* Aij = Arowi;
                if (rm) std::fill_n(Aij,i+1,T(0));
                else for(ptrdiff_t j=i+1;j>0;--j,Aij+=sj) {
#ifdef TMVFLDEBUG
                    TMVAssert(Aij >= A._first);
                    TMVAssert(Aij < A._last);
#endif
                    *Aij = T(0);
                }
            }
        }

        // Do this at end in case A.ptr = x.ptr, so don't mess up x
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
    }

    template <bool a1, bool cx, bool ha, bool cm, bool add, class T, class Tx, class Ta> 
    static void ColRank1Update(
        const Ta alpha, const GenVector<Tx>& x, SymMatrixView<T> A)
    {
        TMVAssert(x.step()==1);
        TMVAssert(x.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.uplo() == Lower);
        TMVAssert(a1 == (alpha == Ta(1)));
        TMVAssert(x.step() == 1);
        TMVAssert(cx == x.isconj());
        TMVAssert(cm == A.iscm());
        TMVAssert(ha == A.isherm());
#ifdef XDEBUG
        TMV_RealType(T) normA = Norm(A);
        TMV_RealType(T) normx = Norm(x);
        TMV_RealType(T) eps = TMV_ABS(alpha)*normx*normx;
        if (add) eps += normA;
        eps *= A.size() * TMV_Epsilon<T>();
#endif

        const ptrdiff_t si = cm ? 1 : A.stepi();
        const ptrdiff_t ds = A.stepj()+si;
        const ptrdiff_t N = A.size();
        const Tx* xj = x.cptr()+N-1;
        T* Ajj = A.ptr()+(N-1)*ds;

        for (ptrdiff_t jj=N,Nmj=1;jj>0;--jj,++Nmj,--xj,Ajj-=ds) {
            if (*xj!=Tx(0)) {
                // Nmj = N-j
                // A.col(j,j,N) += *xj * x.subVector(j,N);
                T* Aij = Ajj;
                const Tx* xi = xj;
                if (a1 || isReal(alpha)) {
                    Tx axj = (ha!=cx) ? TMV_CONJ(*xj) : *xj;
                    if (!a1) axj *= TMV_REAL(alpha);
                    for(ptrdiff_t i=Nmj;i>0;--i,++xi,(cm?++Aij:Aij+=si)) {
                        const T temp = axj * (cx ? TMV_CONJ(*xi) : *xi);
#ifdef TMVFLDEBUG
                        TMVAssert(Aij >= A._first);
                        TMVAssert(Aij < A._last);
#endif
                        if (add) *Aij += temp;
                        else *Aij = temp;
                    }
                } else {
                    T axj = alpha * ((ha!=cx) ? TMV_CONJ(*xj) : *xj);
                    for(ptrdiff_t i=Nmj;i>0;--i,++xi,(cm?++Aij:Aij+=si)) {
                        const T temp = axj * (cx ? TMV_CONJ(*xi) : *xi);
#ifdef TMVFLDEBUG
                        TMVAssert(Aij >= A._first);
                        TMVAssert(Aij < A._last);
#endif
                        if (add) *Aij += temp;
                        else *Aij = temp;
                    }
                }
            } else if (!add) {
                T* Aij = Ajj;
                if (cm) std::fill_n(Aij,Nmj,T(0));
                else for(ptrdiff_t i=Nmj;i>0;--i,Aij+=si) {
#ifdef TMVFLDEBUG
                    TMVAssert(Aij >= A._first);
                    TMVAssert(Aij < A._last);
#endif
                    *Aij = T(0);
                }
            }
        }
        if (ha && isComplex(T())) {
#ifdef XDEBUG
            TMVAssert(NormInf(A.diag().imagPart()) <= eps);
#endif
            A.diag().imagPart().setZero();
        }
    }

    template <bool a1, bool cx, bool add, class T, class Ta, class Tx>
    static void DoRank1Update(
        const Ta alpha, const GenVector<Tx>& x, SymMatrixView<T> A)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.uplo() == Lower);
        TMVAssert(a1 == (alpha == Ta(1)));
        TMVAssert(x.step() == 1);
        TMVAssert(cx == x.isconj());

        if (A.isherm()) 
            if (A.isrm()) RowRank1Update<a1,cx,true,true,add>(alpha,x,A);
            else if (A.iscm()) ColRank1Update<a1,cx,true,true,add>(alpha,x,A);
            else RowRank1Update<a1,cx,true,false,add>(alpha,x,A);
        else
            if (A.isrm()) RowRank1Update<a1,cx,false,true,add>(alpha,x,A);
            else if (A.iscm()) ColRank1Update<a1,cx,false,true,add>(alpha,x,A);
            else RowRank1Update<a1,cx,false,false,add>(alpha,x,A);
    }

    template <bool add, class T, class Tx> 
    static void NonBlasRank1Update(
        const T alpha, const GenVector<Tx>& x, SymMatrixView<T> A)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(alpha != T(0));
        TMVAssert(TMV_IMAG(alpha)==TMV_RealType(T)(0) || !A.isherm());
        TMVAssert(x.size() > 0);

        if (A.uplo() == Upper) 
            return NonBlasRank1Update<add>(
                alpha,x,A.issym()?A.transpose():A.adjoint());
        else if (A.isconj()) 
            return NonBlasRank1Update<add>(
                TMV_CONJ(alpha),x.conjugate(),A.conjugate());
        else if (x.step() != 1) {
            Vector<T> xx = x;
            if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) 
                if (TMV_REAL(alpha) == TMV_RealType(T)(1)) 
                    DoRank1Update<true,false,add>(TMV_REAL(alpha),xx,A);
                else
                    DoRank1Update<false,false,add>(TMV_REAL(alpha),xx,A);
            else 
                DoRank1Update<false,false,add>(alpha,xx,A);
        } else {
            if (x.isconj())
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) 
                    if (TMV_REAL(alpha) == TMV_RealType(T)(1)) 
                        DoRank1Update<true,true,add>(TMV_REAL(alpha),x,A);
                    else
                        DoRank1Update<false,true,add>(TMV_REAL(alpha),x,A);
                else 
                    DoRank1Update<false,true,add>(alpha,x,A);
            else
                if (TMV_IMAG(alpha) == TMV_RealType(T)(0)) 
                    if (TMV_REAL(alpha) == TMV_RealType(T)(1)) 
                        DoRank1Update<true,false,add>(TMV_REAL(alpha),x,A);
                    else
                        DoRank1Update<false,false,add>(TMV_REAL(alpha),x,A);
                else 
                    DoRank1Update<false,false,add>(alpha,x,A);
        }
    }

#ifdef BLAS
    template <class T, class Tx> 
    static inline void BlasRank1Update(
        const T alpha, const GenVector<Tx>& x, SymMatrixView<T> A)
    { NonBlasRank1Update<true>(alpha,x,A); }
#ifdef INST_DOUBLE
    template <> 
    void BlasRank1Update(
        const double alpha, const GenVector<double>& x,
        SymMatrixView<double> A)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(alpha != 0.);
        TMVAssert(x.size() > 0);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.iscm());

        int n=A.size();
        int xs=x.step();
        const double* xp = x.cptr();
        if (xs < 0) xp += (n-1)*xs;
        int lda=A.stepj();
        BLASNAME(dsyr) (
            BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
            BLASV(n),BLASV(alpha),BLASP(xp),BLASV(xs),
            BLASP(A.ptr()),BLASV(lda) BLAS1);
    }
    template <> 
    void BlasRank1Update(
        const std::complex<double> alpha,
        const GenVector<std::complex<double> >& x, 
        SymMatrixView<std::complex<double> > A)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(alpha != 0.);
        TMVAssert(TMV_IMAG(alpha)==0. || !A.isherm());
        TMVAssert(x.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.iscm());

        //std::cout<<"BlasRank1Update:\n";
        //std::cout<<"alpha = "<<alpha<<std::endl;
        //std::cout<<"x = "<<TMV_Text(x)<<"  step "<<x.step()<<std::endl;
        //std::cout<<"A = "<<TMV_Text(A)<<"  U = "<<A.isupper()<<", step "<<A.stepi()<<" "<<A.stepj()<<std::endl;
#ifndef ELAP
        if (A.issym() && x.step() != 1) {
            //std::cout<<"ndef ELAP and not step1, so copy x\n";
            Vector<std::complex<double> > xx = x;
            return BlasRank1Update(alpha,xx,A);
        } else
#endif
            if (A.isherm()) {
                //std::cout<<"A.isherm()\n";
                TMVAssert(TMV_IMAG(alpha)==0.);
                int n=A.size();
                int xs=x.step();
                const std::complex<double>* xp = x.cptr();
                if (xs < 0) xp += (n-1)*xs;
                int lda=A.stepj();
                double ralpha = TMV_REAL(alpha);
                //std::cout<<"Before zher:\n";
                //std::cout<<"n = "<<n<<std::endl;
                //std::cout<<"ralpha = "<<ralpha<<std::endl;
                //std::cout<<"xp = "<<xp<<std::endl;
                //std::cout<<"xs = "<<xs<<std::endl;
                //std::cout<<"Ap = "<<A.ptr()<<std::endl;
                //std::cout<<"lda = "<<lda<<std::endl;
                BLASNAME(zher) (
                    BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
                    BLASV(n),BLASV(ralpha),BLASP(xp),BLASV(xs),
                    BLASP(A.ptr()),BLASV(lda) BLAS1);
            } else {
                //std::cout<<"not A.isherm()\n";
                int n = A.size();
                int lda = A.stepj();
#ifdef ELAP
                int xs = x.step();
                const std::complex<double>* xp = x.cptr();
                if (xs < 0) xp += (n-1)*xs;
                //std::cout<<"Before zsyr:\n";
                //std::cout<<"n = "<<n<<std::endl;
                //std::cout<<"alpha = "<<alpha<<std::endl;
                //std::cout<<"xp = "<<xp<<std::endl;
                //std::cout<<"xs = "<<xs<<std::endl;
                //std::cout<<"Ap = "<<A.ptr()<<std::endl;
                //std::cout<<"lda = "<<lda<<std::endl;
                LAPNAME(zsyr) (
                    LAPCM A.uplo()==Upper ? LAPCH_UP : LAPCH_LO,
                    LAPV(n),LAPP(&alpha),LAPP(xp),LAPV(xs),
                    LAPP(A.ptr()),LAPV(lda) LAP1);
#else
                int k=1;
                std::complex<double> beta(1);
                if (x.step() == 1) {
                    const std::complex<double>* xp = x.cptr();
                    //std::cout<<"Before zsyrk:\n";
                    //std::cout<<"n = "<<n<<std::endl;
                    //std::cout<<"k = "<<k<<std::endl;
                    //std::cout<<"alpha = "<<alpha<<std::endl;
                    //std::cout<<"xp = "<<xp<<std::endl;
                    //std::cout<<"beta = "<<beta<<std::endl;
                    //std::cout<<"Ap = "<<A.ptr()<<std::endl;
                    //std::cout<<"lda = "<<lda<<std::endl;
                    BLASNAME(zsyrk) (
                        BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                        BLASV(n),BLASV(k),BLASP(&alpha),BLASP(xp),BLASV(n),
                        BLASP(&beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
                } else {
                    //std::cout<<"x.step() != 1, so copy alpha*x\n";
                    Vector<std::complex<double> > xx = alpha * x;
                    std::complex<double> xa(1);
                    const std::complex<double>* xp = xx.cptr();
                    //std::cout<<"Before zsyrk:\n";
                    //std::cout<<"n = "<<n<<std::endl;
                    //std::cout<<"k = "<<k<<std::endl;
                    //std::cout<<"xa = "<<xa<<std::endl;
                    //std::cout<<"xp = "<<xp<<std::endl;
                    //std::cout<<"beta = "<<beta<<std::endl;
                    //std::cout<<"Ap = "<<A.ptr()<<std::endl;
                    //std::cout<<"lda = "<<lda<<std::endl;
                    BLASNAME(zsyrk) (
                        BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                        BLASV(n),BLASV(k),BLASP(&xa),BLASP(xp),BLASV(n),
                        BLASP(&beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
                }
#endif
            }
    }
    template <> 
    void BlasRank1Update(
        const std::complex<double> alpha,
        const GenVector<double>& x, 
        SymMatrixView<std::complex<double> > A)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(alpha != 0.);
        TMVAssert(x.size() > 0);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.iscm());

        SymMatrix<double,Lower|ColMajor> A1(A.size(),0.);
        BlasRank1Update(1.,x,A1.view());
        A += alpha*A1;
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void BlasRank1Update(
        const float alpha, const GenVector<float>& x,
        SymMatrixView<float> A)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(alpha != 0.);
        TMVAssert(x.size() > 0);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.iscm());

        int n=A.size();
        int xs=x.step();
        const float* xp = x.cptr();
        if (xs < 0) xp += (n-1)*xs;
        int lda=A.stepj();
        BLASNAME(ssyr) (
            BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
            BLASV(n),BLASV(alpha),BLASP(xp),BLASV(xs),
            BLASP(A.ptr()),BLASV(lda) BLAS1);
    }
    template <> 
    void BlasRank1Update(
        const std::complex<float> alpha,
        const GenVector<std::complex<float> >& x, 
        SymMatrixView<std::complex<float> > A)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(alpha != 0.F);
        TMVAssert(TMV_IMAG(alpha)==0.F || !A.isherm());
        TMVAssert(x.size() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.iscm());

#ifndef ELAP
        if (A.issym() && x.step() != 1) {
            Vector<std::complex<float> > xx = x;
            return BlasRank1Update(alpha,xx,A);
        } else
#endif
            if (A.isherm()) {
                TMVAssert(TMV_IMAG(alpha)==0.F);
                int n=A.size();
                int xs=x.step();
                const std::complex<float>* xp = x.cptr();
                if (xs < 0) xp += (n-1)*xs;
                int lda=A.stepj();
                float ralpha = TMV_REAL(alpha);
                BLASNAME(cher) (
                    BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
                    BLASV(n),BLASV(ralpha),BLASP(xp),BLASV(xs),
                    BLASP(A.ptr()),BLASV(lda) BLAS1);
            } else {
                int n = A.size();
                int lda = A.stepj();
#ifdef ELAP
                int xs = x.step();
                const std::complex<float>* xp = x.cptr();
                if (xs < 0) xp += (n-1)*xs;
                LAPNAME(csyr) (
                    LAPCM A.uplo()==Upper ? LAPCH_UP : LAPCH_LO,
                    LAPV(n),LAPP(&alpha),LAPP(xp),LAPV(xs),
                    LAPP(A.ptr()),LAPV(lda) LAP1);
#else
                int k=1;
                std::complex<float> beta(1);
                if (x.step() == 1) {
                    const std::complex<float>* xp = x.cptr();
                    BLASNAME(csyrk) (
                        BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                        BLASV(n),BLASV(k),BLASP(&alpha),BLASP(xp),BLASV(n),
                        BLASP(&beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
                } else {
                    Vector<std::complex<float> > xx = alpha * x;
                    std::complex<float> xa(1);
                    const std::complex<float>* xp = xx.cptr();
                    BLASNAME(csyrk) (
                        BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                        BLASV(n),BLASV(k),BLASP(&xa),BLASP(xp),BLASV(n),
                        BLASP(&beta),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
                }
#endif
            }
    }
    template <> 
    void BlasRank1Update(
        const std::complex<float> alpha,
        const GenVector<float>& x, 
        SymMatrixView<std::complex<float> > A)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(alpha != 0.F);
        TMVAssert(x.size() > 0);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(A.iscm());

        SymMatrix<float,Lower|ColMajor> A1(A.size(),0.F);
        BlasRank1Update(1.F,x,A1.view());
        A += alpha*A1;
    }
#endif 
#endif // BLAS

    template <bool add, class T, class Tx> 
    void Rank1Update(
        const T alpha, const GenVector<Tx>& x, SymMatrixView<T> A)
    // A (+)= alpha * x * xT
    {
#ifdef XDEBUG
        Vector<Tx> x0 = x;
        Matrix<T> A0 = A;
        Matrix<T> A2 = A;
        if (!add) A2.setZero();
        if (A.isherm())
            A2 += (alpha*x0^x0.conjugate());
        else 
            A2 += (alpha*x0^x0);
        cout<<"Start Rank1Update: alpha = "<<alpha<<endl;
        cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
        cout<<"x = "<<TMV_Text(x)<<"  "<<x<<endl;
        cout<<"herm = "<<A.isherm()<<std::endl;
        cout<<"sym = "<<A.issym()<<std::endl;
        cout<<"x.ptr = "<<x.cptr()<<std::endl;
        cout<<"A.ptr = "<<A.cptr()<<" ... "<<A.cptr()+A.size()*A.size()<<std::endl;
#endif

        TMVAssert(A.size() == x.size());
        TMVAssert(TMV_IMAG(alpha)==TMV_RealType(T)(0) || !A.isherm());
        if (alpha != T(0) && A.size() > 0) {
#ifdef BLAS
            typedef TMV_RealType(T) RT;
            if (!A.iscm() && A.isrm()) {
                //std::cout<<"Transpose A\n";
                return Rank1Update<add>(
                    alpha,x,A.issym()?A.transpose():A.adjoint());
            } else if (A.isconj())  {
                //std::cout<<"Conjugate A\n";
                return Rank1Update<add>(
                    TMV_CONJ(alpha),x.conjugate(),A.conjugate());
            } else if (A.iscm() && A.stepj() >= A.size() && A.stepj() >= 1) {
                //std::cout<<"A.iscm()\n";
                // Most BLAS implementations do fine with the x.step() != 1.
                // However, some implementations seem to propagate nan's from
                // the temporary memory they create to do the unit-1 
                // calculation.
                // So to make sure they don't have to make a temporary, I just
                // do it here for them.
                if (x.step() != 1 || x.isconj() || SameStorage(x,A)) {
                    //std::cout<<"Copy x\n";
                    Vector<Tx> xx = x;
                    if (!add) A.setZero();
                    BlasRank1Update(alpha,xx,A);
                } else {
                    //std::cout<<"No need to copy x\n";
                    if (!add) A.setZero();
                    BlasRank1Update(alpha,x,A);
                }
            } else {
                //std::cout<<"not A.iscm()\n";
                if (A.isherm()) {
                    //std::cout<<"A.isherm()\n";
                    HermMatrix<T,Lower|ColMajor> AA(A.size(),RT(0));
                    Rank1Update<false>(alpha,x,AA.view());
                    if (add) A += AA;
                    else A = AA;
                } else {
                    //std::cout<<"not A.isherm()\n";
                    SymMatrix<T,Lower|ColMajor> AA(A.size(),T(0));
                    Rank1Update<false>(alpha,x,AA.view());
                    if (add) A += AA;
                    else A = AA;
                }
            }
#else
            NonBlasRank1Update<add>(alpha,x,A);
#endif
        }

#ifdef XDEBUG
        TMVAssert(A.isHermOK());
        cout<<"Done Rank1\n";
        cout<<"Norm(A-A2) = "<<Norm(A-A2)<<std::endl;
        if (!(Norm(A-A2) <= 0.001*(TMV_ABS(alpha)*TMV_SQR(Norm(x0))+Norm(A0)))) {
            cerr<<"Rank1Update: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"x = "<<TMV_Text(x)<<"  step = "<<x.step()<<"  "<<x<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"-> A = "<<A<<endl;
            cerr<<"A2 = "<<A2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_Rank1_VVS.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


