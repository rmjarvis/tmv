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
#include "TMV_BandSVDiv.h"
#include "tmv/TMV_BandSVD.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Householder.h"

#ifdef XDEBUG
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_BandMatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

#define RT TMV_RealType(T)

    template <class T> 
    static void MakeBidiagReal(
        Vector<T>& Udiag, Vector<T>& Vtdiag, 
        const GenVector<T>& cD, const GenVector<T>& cE,
        VectorView<T> D, VectorView<T> E, T& )
    {
        TMVAssert(Vtdiag.size() == Udiag.size());
        TMVAssert(cD.size() == D.size());
        TMVAssert(cE.size() == E.size());
        TMVAssert(D.size() == Udiag.size());
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);

        Udiag.setAllTo(T(1));
        Vtdiag.setAllTo(T(1));
        D = cD;
        E = cE;
    }

    template <class T> 
    static void MakeBidiagReal(
        Vector<std::complex<T> >& Udiag, Vector<std::complex<T> >& Vtdiag, 
        const GenVector<std::complex<T> >& cD, 
        const GenVector<std::complex<T> >& cE,
        VectorView<T> D, VectorView<T> E, 
        std::complex<T>& signdet)
    {
        TMVAssert(Vtdiag.size() == Udiag.size());
        TMVAssert(cD.size() == D.size());
        TMVAssert(cE.size() == E.size());
        TMVAssert(D.size() == Udiag.size());
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);

        const ptrdiff_t N = D.size();

        std::complex<T>* Uj = Udiag.ptr();
        std::complex<T>* Vtj = Vtdiag.ptr();
        T* Dj = D.ptr();
        T* Ej = E.ptr();
        Vector<std::complex<T> > xcD = cD;
        Vector<std::complex<T> > xcE = cE;
        const std::complex<T>* cDj = xcD.cptr();
        const std::complex<T>* cEj = xcE.cptr();
#ifdef TMVFLDEBUG
        TMVAssert(Vtj >= Vtdiag._first);
        TMVAssert(Vtj < Vtdiag._last);
#endif
        *Vtj = T(1);
        std::complex<T> newcDj = *cDj;
        for(ptrdiff_t j=0;j<N-1;++j,++Uj,++Dj,++Ej,++cEj) {
#ifdef TMVFLDEBUG
            TMVAssert(Dj >= D._first);
            TMVAssert(Dj < D._last);
            TMVAssert(Ej >= E._first);
            TMVAssert(Ej < E._last);
            TMVAssert(Uj >= Udiag._first);
            TMVAssert(Uj < Udiag._last);
#endif
            *Dj = TMV_ABS(newcDj);
            *Uj = TMV_SIGN(newcDj,*Dj);
            std::complex<T> newcEj = TMV_CONJ(*Uj) * (*cEj);
            *Ej = TMV_ABS(newcEj);
            ++Vtj; // Now Vdiag(j+1)
#ifdef TMVFLDEBUG
            TMVAssert(Vtj >= Vtdiag._first);
            TMVAssert(Vtj < Vtdiag._last);
#endif
            *Vtj = TMV_SIGN(newcEj,*Ej);
            ++cDj; // Now cd(j+1)
            newcDj = TMV_CONJ(*Vtj) * (*cDj);
        }
#ifdef TMVFLDEBUG
        TMVAssert(Dj >= D._first);
        TMVAssert(Dj < D._last);
        TMVAssert(Uj >= Udiag._first);
        TMVAssert(Uj < Udiag._last);
#endif
        *Dj = TMV_ABS(newcDj);
        *Uj = TMV_SIGN(newcDj,*Dj);
        std::complex<T> su, sv;
        DiagMatrixViewOf(Udiag).logDet(&su);
        DiagMatrixViewOf(Vtdiag).logDet(&sv);
        signdet *= su*sv;
    }

    template <class T> 
    static void NonLapBidiagonalize(
        const GenBandMatrix<T>& A,
        MatrixView<T> U, VectorView<RT> D,
        VectorView<RT> E, MatrixView<T> Vt, RT& logdet, T& signdet)
    {
        // Decompose A into U B Vt
        // The Bidiagonal Matrix B is stored as two vectors: D, E
        // D is the diagonal, E is the super-diagonal
        // We use Householder reflections to reduce A to the bidiagonal form:

        TMVAssert(A.rowsize() <= A.colsize());
        TMVAssert(A.rowsize() > 0);
        if (U.cptr()) {
            TMVAssert(U.colsize() == A.colsize());
            TMVAssert(U.rowsize() == A.rowsize());
        } 
        if (Vt.cptr()) {
            TMVAssert(Vt.colsize() == A.rowsize());
            TMVAssert(Vt.rowsize() == A.rowsize());
        }
        TMVAssert(D.size() == A.rowsize());
        TMVAssert(D.size() == E.size()+1);

        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();

        const ptrdiff_t nlo = A.nlo();
        const ptrdiff_t nhi = A.nhi();

        if (nlo == 0 && nhi == 1) {
            Vector<T> Ud(N);
            Vector<T> Vtd(N);
            MakeBidiagReal(Ud,Vtd,A.diag(),A.diag(1),D,E,signdet);
            if (U.cptr()) {
                U.setZero();
                U.diag() = Ud;
            }
            if (Vt.cptr()) {
                Vt.setZero();
                Vt.diag() = Vtd;
            }
        } else if (A.isSquare() && nlo == 1 && nhi == 0) {
            Vector<T> Ud(N);
            Vector<T> Vtd(N);
            MakeBidiagReal(Ud,Vtd,A.diag().reverse(),A.diag(-1).reverse(),
                           D,E,signdet);
            if (U.cptr()) {
                U.setZero();
                U.subVector(N-1,0,-1,1,N) = Ud;
            }
            if (Vt.cptr()) {
                Vt.setZero();
                Vt.subVector(0,N-1,1,-1,N) = Vtd;
            }
        } else {
            auto_ptr<Matrix<T,ColMajor> > UU;
            auto_ptr<MatrixView<T> > U1;
            if (U.cptr()) {
                U = A;
                U1.reset(new MatrixView<T>(U.view()));
            } else {
                UU.reset(new Matrix<T,ColMajor>(A));
                U1.reset(new MatrixView<T>(UU->view()));
            }

            std::vector<ptrdiff_t> vec(N), ver(N-1);
            ptrdiff_t endcol = nlo+1;
            Vector<T> Ubeta(N);
            Vector<T> Vtbeta(N-1);

            T* Ubj = Ubeta.ptr();
            for(ptrdiff_t j=0;j<N-1;++j,++Ubj) {
                vec[j] = endcol;
                ptrdiff_t endrow = TMV_MIN(endcol+nhi,N);
                ver[j] = endrow;
#ifdef TMVFLDEBUG
                TMVAssert(Ubj >= Ubeta._first);
                TMVAssert(Ubj < Ubeta._last);
#endif
                *Ubj = HouseholderReflect(
                    U1->subMatrix(j,endcol,j,endrow),signdet);
                if (endcol < M) endcol = TMV_MIN(endrow+nlo,M);
                Vtbeta(j) = HouseholderReflect(
                    U1->transpose().subMatrix(j+1,endrow,j,endcol),signdet);
            }
            vec[N-1] = endcol;
#ifdef TMVFLDEBUG
            TMVAssert(Ubj >= Ubeta._first);
            TMVAssert(Ubj < Ubeta._last);
#endif
            *Ubj = HouseholderReflect(U1->subMatrix(N-1,endcol,N-1,N),signdet);

            // Now U stores Householder vectors for U in lower diagonal columns (HLi)
            // and Householder vectors for Vt in upper diagonal rows (HRi)
            // except for the bidiagonal which is the bidiagonal we want:
            D = U1->diag().realPart();
            E = U1->diag(1).realPart();
            if (isComplex(T())) {
                TMVAssert(NormInf(U1->diag().imagPart()) == RT(0));
                TMVAssert(NormInf(U1->diag(1).imagPart()) == RT(0));
            }

            if (Vt.cptr()) {
                Vt.setToIdentity();
                for (ptrdiff_t j=N-2;j>=0;--j) {
                    Vt.row(j+1,j+2,ver[j]) = U1->row(j,j+2,ver[j]);
                    HouseholderUnpack(
                        Vt.transpose().subMatrix(j+1,ver[j],j+1,N),Vtbeta(j));
                }
            }

            if (U.cptr()) {
                U.diag().setZero();
                U.diag(1).setZero();
                // Ubj is currently &U(N-1)
                HouseholderUnpack(U.subMatrix(N-1,vec[N-1],N-1,N),*Ubj);
                for (ptrdiff_t j=N-2;j>=0;--j) {
                    U.row(j,j,ver[j]).setZero();
                    HouseholderUnpack(U.subMatrix(j,vec[j],j,N),*(--Ubj));
                }
            }
        }
        if (signdet != T(0)) {
            RT s;
            logdet += DiagMatrixViewOf(D).logDet(&s);
            signdet *= s;
        }
    }

#ifdef LAP
    template <class T> 
    static inline void LapBidiagonalize(
        const GenBandMatrix<T>& A,
        MatrixView<T> U, VectorView<RT> D,
        VectorView<RT> E, MatrixView<T> Vt, RT& logdet, T& signdet)
    { NonLapBidiagonalize(A,U,D,E,Vt,logdet,signdet); }
#ifdef INST_DOUBLE
    template <> 
    void LapBidiagonalize(
        const GenBandMatrix<double>& A, MatrixView<double> U,
        VectorView<double> D, VectorView<double> E,
        MatrixView<double> Vt, double& logdet, double& signdet)
    {
        TMVAssert(A.rowsize() == A.colsize());
        // The Lap routines can do NonSquare matrices, but they want to
        // write out to a square (MxM) U matrix which is larger than
        // what we have stored here.
        TMVAssert(A.rowsize() > 0);
        if (U.cptr()) {
            TMVAssert(U.colsize() == A.colsize());
            TMVAssert(U.rowsize() == A.rowsize());
            TMVAssert(U.iscm());
            TMVAssert(U.ct() == NonConj);
        }
        if (Vt.cptr()) {
            TMVAssert(Vt.colsize() == A.rowsize());
            TMVAssert(Vt.rowsize() == A.rowsize());
            TMVAssert(Vt.iscm());
            TMVAssert(Vt.ct() == NonConj);
        }
        TMVAssert(D.size() == A.rowsize());
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.step()==1);
        TMVAssert(E.step()==1);
        TMVAssert(D.ct() == NonConj);
        TMVAssert(E.ct() == NonConj);

        int m = A.colsize();
        int n = A.rowsize();
        int ncc = 0;
        int kl = A.nlo();
        int ku = A.nhi();
#ifndef LAPNOWORK
        int lwork = 2*TMV_MAX(m,n);
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        // LAP version overwrites original BandMatrix with crap.
        // Hence, copy BandMatrix before running.
        BandMatrix<double,ColMajor> A2 = A;
        if (U.cptr()) U.setZero();
        if (Vt.cptr()) Vt.setZero();
        D.setZero();
        E.setZero();
        int lda = A2.diagstep();
        int ldu = U.stepj();
        int ldv = Vt.stepj();
        char vect = U.cptr() ? Vt.cptr() ? 'B' : 'Q' : Vt.cptr() ? 'P' : 'N';
        double* VV = Vt.ptr();
        double* UU = U.ptr();
        int Lap_info=0;
        LAPNAME(dgbbrd) (
            LAPCM LAPV(vect),LAPV(m),LAPV(n),LAPV(ncc),
            LAPV(kl),LAPV(ku),
            LAPP(A2.cptr()-A.nhi()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(UU),LAPV(ldu),LAPP(VV),LAPV(ldv),
            0,LAPV(n) LAPWK(work.get()) LAPINFO LAP1 );

        if (signdet != 0.) {
            double s;
            logdet += DiagMatrixViewOf(D).logDet(&s);
            signdet *= s;
        }
        LAP_Results(Lap_info,"dgbbrd");
    }
    template <> 
    void LapBidiagonalize(
        const GenBandMatrix<std::complex<double> >& A,
        MatrixView<std::complex<double> > U,
        VectorView<double> D, VectorView<double> E,
        MatrixView<std::complex<double> > Vt, 
        double& logdet, std::complex<double>& signdet)
    {
        TMVAssert(A.rowsize() == A.colsize());
        TMVAssert(A.rowsize() > 0);
        if (U.cptr()) {
            TMVAssert(U.colsize() == A.colsize());
            TMVAssert(U.rowsize() == A.rowsize());
            TMVAssert(U.iscm());
            TMVAssert(U.ct() == NonConj);
        }
        if (Vt.cptr()) {
            TMVAssert(Vt.colsize() == A.rowsize());
            TMVAssert(Vt.rowsize() == A.rowsize());
            TMVAssert(Vt.iscm());
            TMVAssert(Vt.ct() == NonConj);
        }
        TMVAssert(D.size() == A.rowsize());
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.step()==1);
        TMVAssert(E.step()==1);
        TMVAssert(D.ct() == NonConj);
        TMVAssert(E.ct() == NonConj);

        int m = A.colsize();
        int n = A.rowsize();
        int o = 0;
        int kl = A.nlo();
        int ku = A.nhi();
#ifndef LAPNOWORK
        int lwork = TMV_MAX(m,n);
        AlignedArray<std::complex<double> > work(lwork);
        AlignedArray<double> rwork(lwork);
        VectorViewOf(work.get(),lwork).setZero();
        VectorViewOf(rwork.get(),lwork).setZero();
#endif
        BandMatrix<std::complex<double>,ColMajor> A2 = A;
        if (U.cptr()) U.setZero();
        if (Vt.cptr()) Vt.setZero();
        D.setZero();
        E.setZero();
        int lda = A2.diagstep();
        int ldu = U.stepj();
        int ldv = Vt.stepj();
        char vect = U.cptr() ? Vt.cptr() ? 'B' : 'Q' : Vt.cptr() ? 'P' : 'N';
        std::complex<double>* VV = Vt.ptr();
        std::complex<double>* UU = U.ptr();
        int Lap_info=0;
        LAPNAME(zgbbrd) (
            LAPCM LAPV(vect),LAPV(m),LAPV(n),LAPV(o),
            LAPV(kl),LAPV(ku),LAPP(A2.cptr()-A.nhi()),LAPV(lda),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(UU),LAPV(ldu),
            LAPP(VV),LAPV(ldv),0,LAPV(n)
            LAPWK(work.get()) LAPWK(rwork.get()) LAPINFO LAP1);

        // OK, I lied, the output A2 isn't complete crap.  Its diagonal
        // holds a version of D which includes all the complex arguments
        // so it can be used to find the determinant.
        if (signdet != 0.) {
            std::complex<double> s;
            logdet += DiagMatrixViewOf(A2.diag()).logDet(&s);
            signdet *= s;
        }
        LAP_Results(Lap_info,"zgbbrd");
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapBidiagonalize(
        const GenBandMatrix<float>& A, MatrixView<float> U,
        VectorView<float> D, VectorView<float> E,
        MatrixView<float> Vt, float& logdet, float& signdet)
    {
        TMVAssert(A.rowsize() == A.colsize());
        // The Lap routines can do NonSquare matrices, but they want to
        // write out to a square (MxM) U matrix which is larger than
        // what we have stored here.
        TMVAssert(A.rowsize() > 0);
        if (U.cptr()) {
            TMVAssert(U.colsize() == A.colsize());
            TMVAssert(U.rowsize() == A.rowsize());
            TMVAssert(U.iscm());
            TMVAssert(U.ct() == NonConj);
        }
        if (Vt.cptr()) {
            TMVAssert(Vt.colsize() == A.rowsize());
            TMVAssert(Vt.rowsize() == A.rowsize());
            TMVAssert(Vt.iscm());
            TMVAssert(Vt.ct() == NonConj);
        }
        TMVAssert(D.size() == A.rowsize());
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.step()==1);
        TMVAssert(E.step()==1);
        TMVAssert(D.ct() == NonConj);
        TMVAssert(E.ct() == NonConj);

        int m = A.colsize();
        int n = A.rowsize();
        int ncc = 0;
        int kl = A.nlo();
        int ku = A.nhi();
#ifndef LAPNOWORK
        int lwork = 2*TMV_MAX(m,n);
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        BandMatrix<float,ColMajor> A2 = A;
        if (U.cptr()) U.setZero();
        if (Vt.cptr()) Vt.setZero();
        D.setZero();
        E.setZero();
        int lda = A2.diagstep();
        int ldu = U.stepj();
        int ldv = Vt.stepj();
        char vect = U.cptr() ? Vt.cptr() ? 'B' : 'Q' : Vt.cptr() ? 'P' : 'N';
        float* VV = Vt.ptr();
        float* UU = U.ptr();
        int Lap_info=0;
        LAPNAME(sgbbrd) (
            LAPCM LAPV(vect),LAPV(m),LAPV(n),LAPV(ncc),
            LAPV(kl),LAPV(ku),
            LAPP(A2.cptr()-A.nhi()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(UU),LAPV(ldu),LAPP(VV),LAPV(ldv),
            0,LAPV(n) LAPWK(work.get()) LAPINFO LAP1 );

        if (signdet != 0.F) {
            float s;
            logdet += DiagMatrixViewOf(D).logDet(&s);
            signdet *= s;
        }
        LAP_Results(Lap_info,"sgbbrd");
    }
    template <> 
    void LapBidiagonalize(
        const GenBandMatrix<std::complex<float> >& A,
        MatrixView<std::complex<float> > U,
        VectorView<float> D, VectorView<float> E,
        MatrixView<std::complex<float> > Vt, 
        float& logdet, std::complex<float>& signdet)
    {
        TMVAssert(A.rowsize() == A.colsize());
        TMVAssert(A.rowsize() > 0);
        if (U.cptr()) {
            TMVAssert(U.colsize() == A.colsize());
            TMVAssert(U.rowsize() == A.rowsize());
            TMVAssert(U.iscm());
            TMVAssert(U.ct() == NonConj);
        }
        if (Vt.cptr()) {
            TMVAssert(Vt.colsize() == A.rowsize());
            TMVAssert(Vt.rowsize() == A.rowsize());
            TMVAssert(Vt.iscm());
            TMVAssert(Vt.ct() == NonConj);
        }
        TMVAssert(D.size() == A.rowsize());
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.step()==1);
        TMVAssert(E.step()==1);
        TMVAssert(D.ct() == NonConj);
        TMVAssert(E.ct() == NonConj);

        int m = A.colsize();
        int n = A.rowsize();
        int o = 0;
        int kl = A.nlo();
        int ku = A.nhi();
#ifndef LAPNOWORK
        int lwork = TMV_MAX(m,n);
        AlignedArray<std::complex<float> > work(lwork);
        AlignedArray<float> rwork(lwork);
        VectorViewOf(work.get(),lwork).setZero();
        VectorViewOf(rwork.get(),lwork).setZero();
#endif
        BandMatrix<std::complex<float>,ColMajor> A2 = A;
        if (U.cptr()) U.setZero();
        if (Vt.cptr()) Vt.setZero();
        D.setZero();
        E.setZero();
        int lda = A2.diagstep();
        int ldu = U.stepj();
        int ldv = Vt.stepj();
        char vect = U.cptr() ? Vt.cptr() ? 'B' : 'Q' : Vt.cptr() ? 'P' : 'N';
        std::complex<float>* VV = Vt.ptr();
        std::complex<float>* UU = U.ptr();
        int Lap_info=0;
        LAPNAME(cgbbrd) (
            LAPCM LAPV(vect),LAPV(m),LAPV(n),LAPV(o),
            LAPV(kl),LAPV(ku),LAPP(A2.cptr()-A.nhi()),LAPV(lda),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(UU),LAPV(ldu),
            LAPP(VV),LAPV(ldv),0,LAPV(n)
            LAPWK(work.get()) LAPWK(rwork.get()) LAPINFO LAP1);

        if (signdet != 0.F) {
            std::complex<float> s;
            logdet += DiagMatrixViewOf(A2.diag()).logDet(&s);
            signdet *= s;
        }
        LAP_Results(Lap_info,"cgbbrd");
    }
#endif
#endif
    template <class T> 
    static void Bidiagonalize(
        const GenBandMatrix<T>& A, MatrixView<T> U,
        VectorView<RT> D, VectorView<RT> E,
        MatrixView<T> Vt, RT& logdet, T& signdet)
    {
        TMVAssert(A.rowsize() <= A.colsize());
        TMVAssert(A.rowsize() > 0);
        if (U.cptr()) {
            TMVAssert(U.colsize() == A.colsize());
            TMVAssert(U.rowsize() == A.rowsize());
            TMVAssert(U.ct() == NonConj);
        }
        if (Vt.cptr()) {
            TMVAssert(Vt.colsize() == A.rowsize());
            TMVAssert(Vt.rowsize() == A.rowsize());
            TMVAssert(Vt.ct() == NonConj);
        }
        TMVAssert(D.size() == A.rowsize());
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.step()==1);
        TMVAssert(E.step()==1);
        TMVAssert(D.ct() == NonConj);
        TMVAssert(E.ct() == NonConj);

#ifdef XDEBUG
        std::cout<<"Start Band Bidiagonalize:\n";
        std::cout<<"A = "<<A<<std::endl;
        Matrix<T> A0(A);
#ifdef LAP
        BandMatrix<T> A2(A);
        Vector<RT> D2(D);
        Vector<RT> E2(E);
        Matrix<T> U2(A.colsize(),A.rowsize());
        Matrix<T> Vt2(D.size(),D.size());
        RT logdet2(0);
        T signdet2(1);
        NonLapBidiagonalize<T>(
            A2,U2.view(),D2.view(),E2.view(),Vt2.view(),logdet2,signdet2);
#endif
#endif

        if (A.rowsize() > 0) {
            TMVAssert(E.size() == D.size()-1);
#ifdef LAP
            if (A.isSquare() && 
                (!U.cptr() || U.iscm()) && (!Vt.cptr() || Vt.iscm())) 
                LapBidiagonalize(A,U,D,E,Vt,logdet,signdet);
            else 
#endif
                NonLapBidiagonalize(A,U,D,E,Vt,logdet,signdet);
        }
#ifdef XDEBUG
        if (U.cptr() && Vt.cptr()) {
            std::cout<<"Done Band Bidiagonalize:\n";
            std::cout<<"U = "<<U<<std::endl;
            std::cout<<"D = "<<D<<std::endl;
            std::cout<<"E = "<<E<<std::endl;
            std::cout<<"Vt = "<<Vt<<std::endl;
            Matrix<T> UBVt = U*UpperBiDiagMatrix(D,E)*Vt;
            std::cout<<"UBVt = "<<UBVt<<std::endl;
            std::cout<<"Norm(UBVt-A0) = "<<Norm(UBVt-A0)<<std::endl;
            if (!(Norm(UBVt-A0) < 0.001*Norm(A0))) {
                cerr<<"Bidiagonalize:\n";
                cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
                cerr<<"-> D = "<<D<<endl;
                cerr<<"E = "<<E<<endl;
                cerr<<"U = "<<U<<endl;
                cerr<<"Vt = "<<Vt<<endl;
#ifdef LAP
                cerr<<"Nonlap D = "<<D2<<endl;
                cerr<<"Norm(diff) = "<<Norm(D-D2)<<endl;
                cerr<<"Nonlap E = "<<E2<<endl;
                cerr<<"Norm(diff) = "<<Norm(E-E2)<<endl;
                cerr<<"U2 = "<<U2<<endl;
                cerr<<"Vt2 = "<<Vt2<<endl;
#endif
                cerr<<"UBVt = "<<UBVt<<endl;
                cerr<<"Norm(UBVt-A0) = "<<Norm(UBVt-A0)<<endl;
                abort();
            }
        }
#endif
    }

    template <class T> 
    void SV_Decompose(
        const GenBandMatrix<T>& A, 
        MatrixView<T> U, DiagMatrixView<RT> S,
        MatrixView<T> Vt, RT& logdet, T& signdet)
    {
        // Decompose A into U S Vt
        // where S is a diagonal real matrix, and U,V are unitary matrices.
        // All matrices are square N x N
        // The determinant is kept track of in det.
        //
        // Everything is identical to the regular SVD except for the 
        // Bidiagonal Step.

        TMVAssert(A.rowsize() <= A.colsize());
        TMVAssert(A.rowsize() > 0);
        if (U.cptr()) {
            TMVAssert(U.rowsize() == A.rowsize());
            TMVAssert(U.colsize() == A.colsize());
            TMVAssert(U.ct() == NonConj);
        } 
        if (Vt.cptr()) {
            TMVAssert(Vt.rowsize() == A.rowsize());
            TMVAssert(Vt.colsize() == A.rowsize());
            TMVAssert(Vt.ct() == NonConj);
        }
        TMVAssert(S.size() == A.rowsize());
        TMVAssert(S.diag().ct() == NonConj);

        if (A.nlo() == 0 && A.nhi() == 0) {
            if (U.cptr()) U.setToIdentity();
            if (Vt.cptr()) Vt.setToIdentity();

            if (signdet != T(0)) {
                T s;
                logdet += DiagMatrixViewOf(A.diag()).logDet(&s);
                signdet *= s;
            }
            const ptrdiff_t N = A.rowsize();
            const T* Ajj = A.cptr();
            const ptrdiff_t Ads = A.stepi()+A.stepj();
            RT* Sj = S.diag().ptr();
            const ptrdiff_t Ss = S.diag().step();
            T* Ujj = U.ptr();
            const ptrdiff_t Uds = U.stepi()+U.stepj();

            if (A.isconj()) {
                for(ptrdiff_t j=0;j<N;++j,Ajj+=Ads,Sj+=Ss) {
#ifdef TMVFLDEBUG
                    TMVAssert(Sj >= S.diag()._first);
                    TMVAssert(Sj < S.diag()._last);
#endif
                    *Sj = TMV_ABS(*Ajj);
                    if(U.cptr()) {
#ifdef TMVFLDEBUG
                        TMVAssert(Ujj >= U._first);
                        TMVAssert(Ujj < U._last);
#endif
                        *Ujj = TMV_SIGN(TMV_CONJ(*Ajj),*Sj); Ujj += Uds; 
                    }
                }
            } else {
                for(ptrdiff_t j=0;j<N;++j,Ajj+=Ads,Sj+=Ss) {
#ifdef TMVFLDEBUG
                    TMVAssert(Sj >= S.diag()._first);
                    TMVAssert(Sj < S.diag()._last);
#endif
                    *Sj = TMV_ABS(*Ajj);
                    if(U.cptr()) { 
#ifdef TMVFLDEBUG
                        TMVAssert(Ujj >= U._first);
                        TMVAssert(Ujj < U._last);
#endif
                        *Ujj = TMV_SIGN(*Ajj,*Sj); Ujj += Uds; 
                    }
                }
            }
            AlignedArray<ptrdiff_t> sortp(N);
            S.diag().sort(sortp.get(),Descend);
            if (U.cptr()) U.permuteCols(sortp.get());
            if (Vt.cptr()) Vt.permuteRows(sortp.get());
        } else {
            Vector<RT> E(S.size()-1);
            Bidiagonalize(A,U,S.diag(),E.view(),Vt,logdet,signdet);

            SV_DecomposeFromBidiagonal(U,S.diag(),E.view(),Vt);
        }
    }

    template <class T> 
    void SV_Decompose(
        const GenBandMatrix<T>& A,
        MatrixView<T> U, DiagMatrixView<RT> S, MatrixView<T> Vt)
    { 
        if (U.isconj()) {
            if (Vt.isconj()) {
                SV_Decompose(A.conjugate(),U.conjugate(),S,Vt.conjugate());
            } else {
                SV_Decompose(A.conjugate(),U.conjugate(),S,Vt);
                Vt.conjugateSelf();
            }
        } else {
            if (Vt.isconj()) {
                SV_Decompose(A,U,S,Vt.conjugate());
                Vt.conjugateSelf();
            } else {
                RT ld=0; T d=0; SV_Decompose<T>(A,U,S,Vt,ld,d); 
            }
        }
    }

    template <class T> 
    void SV_Decompose(
        const GenBandMatrix<T>& A, MatrixView<T> U, DiagMatrixView<RT> S)
    {
        if (U.isconj()) {
            SV_Decompose(A.conjugate(),U.conjugate(),S);
        } else {
            RT ld=0; T d=0;
            MatrixView<T> Vt(0,0,0,1,1,NonConj);
            SV_Decompose<T>(A,U,S,Vt,ld,d); 
        }
    }

    template <class T> 
    void SV_Decompose(
        const GenBandMatrix<T>& A, DiagMatrixView<RT> S, MatrixView<T> Vt)
    {
        if (Vt.isconj()) {
            SV_Decompose(A.conjugate(),S,Vt.conjugate());
        } else {
            RT ld=0; T d=0;
            MatrixView<T> U(0,0,0,1,1,NonConj);
            SV_Decompose<T>(A,U,S,Vt,ld,d); 
        }
    }

    template <class T> 
    void SV_Decompose(const GenBandMatrix<T>& A, DiagMatrixView<RT> S)
    {
        RT ld=0; T d=0;
        MatrixView<T> U(0,0,0,1,1,NonConj);
        MatrixView<T> Vt(0,0,0,1,1,NonConj);
        SV_Decompose<T>(A,U,S,Vt,ld,d); 
    }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_BandSVDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


