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
#include "TMV_SymSVDiv.h"
#include "tmv/TMV_Householder.h"
#include "tmv/TMV_SymHouseholder.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"

#ifdef XDEBUG
#include "tmv/TMV_DiagMatrixArith.h"
#include "TMV_QRDiv.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

#ifdef TMV_BLOCKSIZE
#define SYM_TRIDIAG_BLOCKSIZE TMV_BLOCKSIZE
#else
#define SYM_TRIDIAG_BLOCKSIZE 64
#endif

namespace tmv {

#define RT TMV_RealType(T)

    template <class T, class Td> 
    static void NonBlockTridiagonalize(
        SymMatrixView<T> A, VectorView<T> beta,
        VectorView<Td> D, VectorView<RT> E, 
        T& signdet, ptrdiff_t j1=0)
    {
        // Decompose A into U T Ut
        // The Tridiagonal Matrix T is stored as two vectors: D, E
        // D is the diagonal, E is the sub-diagonal.
        // If A is herm, E* is the super-diagonal, otherwise E is.
        // However, the Householder reflections make E real, so this
        // distinction is irrelevant.
        // A along with beta hold the U matrix.
        const ptrdiff_t N = A.size();

        TMVAssert(N > 0);
        TMVAssert(beta.size() == N-1);
        TMVAssert(D.size() == N);
        TMVAssert(E.size() == N-1);
        TMVAssert(isReal(Td()) || !A.isherm());
        TMVAssert(A.uplo() == Lower);
        TMVAssert(beta.step() == 1);
        TMVAssert(!beta.isconj());

        // We use Householder reflections to reduce A to the tridiagonal form:
        T* bj = beta.ptr() + j1;
        for(ptrdiff_t j=j1;j<N-1;++j,++bj) {
#ifdef TMVFLDEBUG
            TMVAssert(bj >= beta._first);
            TMVAssert(bj < beta._last);
#endif
            *bj = HouseholderReflect(A.col(j,j+1,N),signdet);
            if (*bj != T(0)) 
                HouseholderLRMult(A.col(j,j+2,N),*bj,A.subSymMatrix(j+1,N));
        }

        // The tridiagonal of A is the tridiagonal we want, so copy it to D,E
        if (isReal(Td())) D = A.diag().realPart();
        else D = A.diag();
        E = A.diag(1).realPart();

        if (!A.isherm()) signdet *= signdet;
    }

    template <class T, class Td> 
    static void BlockTridiagonalize(
        SymMatrixView<T> A, VectorView<T> beta,
        VectorView<Td> D, VectorView<RT> E, T& signdet)
    {
        // The block version here is basically the same as the block
        // version of Bidiagonalization for normal matrices.  We have the 
        // same trouble that storing just the Z matrix (as we do for QR
        // decomposition) doesn't work, since we need to multiply it by
        // parts of the matrix which gets overwritten.
        // 
        // The (simplifying) difference is that the two temporary matrices
        // in the bidiagonalization routine are equal here, so there
        // is only one temporary matrix.

        const ptrdiff_t N = A.size();
        const ptrdiff_t Nm1 = N-1;

        TMVAssert(N > 0);
        TMVAssert(beta.size() == N-1);
        TMVAssert(D.size() == N);
        TMVAssert(E.size() == N-1);
        TMVAssert(isReal(Td()) || !A.isherm());
        TMVAssert(A.uplo() == Lower);
        TMVAssert(beta.step() == 1);
        TMVAssert(!beta.isconj());

        Matrix<T,ColMajor> ZYtmT(N,TMV_MIN(SYM_TRIDIAG_BLOCKSIZE,int(Nm1)));
        // It's easier (for me) to think about the calculations with ZYtm.
        // But the eventual Rank2KUpdate call will have unconjugated matrices
        // if we store (ZYtm)t for Hermitian or (ZYtm)T for symmetric.
        // This is what is stored in the Matrix ZYtmT.

        T* bj = beta.ptr();
        for(ptrdiff_t j1=0;j1<Nm1;) {
            ptrdiff_t j2 = TMV_MIN(Nm1,j1+SYM_TRIDIAG_BLOCKSIZE);

            if (j2 < N) {
                for(ptrdiff_t j=j1,jj=0;j<j2;++j,++jj,++bj) { // jj = j-j1

                    // Update current column:
                    // A(j:N,j) -= Y(j:N,0:j) ZYtm(0:j,j) + (ZYtm)t(j:N,0:j) Yt(0:j,j)
                    VectorView<T> Acolj = A.col(j,j,N);
                    MatrixView<T> Y0 = A.subMatrix(j,N,j1,j);
                    MatrixView<T> Y0T = (A.issym() ? Y0.transpose() : Y0.adjoint());
                    MatrixView<T> ZYtm0T = ZYtmT.subMatrix(j,N,0,jj);
                    MatrixView<T> ZYtm0 = (A.issym() ? ZYtm0T.transpose() : 
                                           ZYtm0T.adjoint());
                    if (jj > 0) {
                        Acolj -= Y0 * ZYtm0.col(0);
                        Acolj -= ZYtm0T * Y0T.col(0);
                    }

                    // Do the Householder reflection
                    const ptrdiff_t Nmj = Acolj.size();
                    VectorView<T> u = Acolj.subVector(1,Nmj);
#ifdef TMVFLDEBUG
                    TMVAssert(bj >= beta._first);
                    TMVAssert(bj < beta._last);
#endif
                    *bj = HouseholderReflect(u,signdet);

                    T* u0 = u.ptr();
                    TMVAssert(TMV_IMAG(*u0) == RT(0));
                    E(j) = TMV_REAL(*u0);
                    // Set top of u to 1 explicitly.  Calculations below are easier this 
                    // way, rather than dealing with it being implicitly 1.
#ifdef TMVFLDEBUG
                    TMVAssert(u0 >= A._first);
                    TMVAssert(u0 < A._last);
#endif
                    *u0 = T(1);

                    // Update ZYtmT
                    // ZYtm(j,j+1:N)  = Z(j,0:j+1) Yt(0:j+1,0:N) m'(0:N,j+1:N)
                    //                = beta Yt(j,0:N) m'(0:N,j+1:N)
                    //                = beta Yt(j,j+1:N) m'(j+1:N,j+1:N)
                    //                = beta ut m'(j+1:N,j+1:N)
                    // ZYtmT(j+1:N,j) = bb m'(j+1:N,j+1:N) uu
                    //                = bb [ A(j+1:N,j+1:N) uu
                    //                     - Y(j+1:N,0:j) ZYtm(0:j,j+1:N) uu
                    //                     - ZYtmT(j+1:N,0:j) YT(0:j,j+1:N) uu ]
                    //                = bb [ A(j+1:N,j+1:N) uu
                    //                       - Y1 ZYtm1 uu - ZYtm1T Y1T uu ]
                    // where {bb = beta*, uu = u} if hermitian
                    // and   {bb = beta, uu = u*} if symmetric

                    MatrixView<T> ZYtm1 = ZYtm0.colRange(1,Nmj);
                    MatrixView<T> ZYtm1T = ZYtm0T.rowRange(1,Nmj);
                    MatrixView<T> Y1 = Y0.rowRange(1,Nmj);
                    MatrixView<T> Y1T = Y0T.colRange(1,Nmj);
                    VectorView<T> ZYtmTj = ZYtmT.col(jj,j+1,N);
                    VectorView<T> temp = ZYtmT.col(jj,j1,j);
                    VectorView<T> uu = A.isherm() ? u : u.conjugate();
                    T bb = A.issym() ? *bj : TMV_CONJ(*bj);

                    ZYtmTj = A.subSymMatrix(j+1,N) * uu;
                    if (jj > 0) {
                        ZYtmTj -= Y1 * (temp = ZYtm1 * uu);
                        ZYtmTj -= ZYtm1T * (temp = Y1T * uu);
                    }
                    ZYtmTj *= bb;

                    // There is an extra complication as well now.
                    // What we really want to do to block update the matrix 
                    // would be:
                    // m = m - YZYtm - mYZtYt + YZYtmYZtYt
                    // However, we can fold that last term into the first two
                    // by making the ZYtm that we store really equal to
                    // ZYtm - 1/2 ZYtmYZtYt
                    VectorView<T> ZYtmj =
                        A.issym() ? ZYtmTj : ZYtmTj.conjugate();
                    ZYtmTj -= RT(0.5)*bb*(ZYtmj*uu)*u;
                }

                // Update the rest of the matrix:
                if (j2 < N)
                    Rank2KUpdate<true>(
                        T(-1),A.subMatrix(j2,N,j1,j2),
                        ZYtmT.subMatrix(j2,N,0,j2-j1),A.subSymMatrix(j2,N));

            } else NonBlockTridiagonalize(A,beta,D,E,signdet,j1);

            j1 = j2;
        }

        // E has already been set, but the diagonal of A needs to be copied to D.
        if (isReal(Td())) D = A.diag().realPart();
        else D = A.diag();

        if (!A.isherm()) signdet *= signdet;
    }

    template <class T, class Td> 
    static inline void NonLapTridiagonalize(
        SymMatrixView<T> A, VectorView<T> beta,
        VectorView<Td> D, VectorView<RT> E, T& signdet)
    {
        if (A.size() >= SYM_TRIDIAG_BLOCKSIZE*2)
            BlockTridiagonalize(A,beta,D,E,signdet);
        else
            NonBlockTridiagonalize(A,beta,D,E,signdet);
    }

#ifdef LAP
    template <class T, class Td> 
    static inline void LapTridiagonalize(
        SymMatrixView<T> A, VectorView<T> beta,
        VectorView<Td> D, VectorView<RT> E, T& signdet)
    { NonLapTridiagonalize(A,beta.subVector(0,A.size()-1),D,E,signdet); }
#ifdef INST_DOUBLE
    template <> 
    void LapTridiagonalize(
        SymMatrixView<double> A, VectorView<double> beta,
        VectorView<double> D, VectorView<double> E, double& )
    {
        TMVAssert(A.iscm());
        TMVAssert(A.uplo() == Lower);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(beta.step() == 1);
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size()+1 == A.size());
        TMVAssert(beta.size() == A.size());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        TMVAssert(beta.ct()==NonConj);

        int n = A.size();
        int ldu = A.stepj();
        beta.setZero();
        D.setZero();
        E.setZero();
        int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<double> work(1);
        work.get()[0] = 0.;
        LAPNAME(dsytrd) (
            LAPCM LAPCH_LO,LAPV(n),
            LAPP(A.ptr()),LAPV(ldu),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
        lwork = int(work[0]);
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(dsytrd) (
            LAPCM LAPCH_LO,LAPV(n),
            LAPP(A.ptr()),LAPV(ldu),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"dsytrd");
#else
        LAP_Results(Lap_info,int(work[0]),n,n,lwork,"dsytrd");
#endif
    }
    template <> 
    void LapTridiagonalize(
        SymMatrixView<std::complex<double> > A,
        VectorView<std::complex<double> > beta,
        VectorView<double> D, VectorView<double> E, 
        std::complex<double>& )
    {
        //std::cout<<"Start LapTridiagonalize:\n";
        //std::cout<<"A = "<<TMV_Text(A)<<std::endl;
        //std::cout<<"beta = "<<TMV_Text(beta)<<"  "<<beta.step()<<std::endl;
        //std::cout<<"D = "<<TMV_Text(D)<<"  "<<D.step()<<std::endl;
        //std::cout<<"E = "<<TMV_Text(E)<<"  "<<E.step()<<std::endl;

        TMVAssert(A.iscm());
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.isherm());
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(beta.step() == 1);
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size()+1 == A.size());
        TMVAssert(beta.size() == A.size());
        // Note: beta.size() should only need to be n-1, but some 
        // LAPACK implementations write to beta[n-1], which can lead to 
        // segmentation faults or at least wrong calculations.
        // So we require beta.size() to be n.
        // I have only found this to be the case for the complex varieties,
        // but I require it for the real ones as well just to be safe.
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        TMVAssert(beta.ct()==NonConj);

        int n = A.size();
        int ldu = A.stepj();
        beta.setZero();
        D.setZero();
        E.setZero();
        int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<std::complex<double> > work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<std::complex<double> > work(1);
        work.get()[0] = 0.;
        LAPNAME(zhetrd) (
            LAPCM LAPCH_LO,LAPV(n),
            LAPP(A.ptr()),LAPV(ldu),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
        lwork = int(TMV_REAL(work[0]));
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        //std::cout<<"Before zhetrd\n";
        LAPNAME(zhetrd) (
            LAPCM LAPCH_LO,LAPV(n),
            LAPP(A.ptr()),LAPV(ldu),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
        //std::cout<<"After zhetrd\n";
        if (!A.isconj()) beta.conjugateSelf();
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"zhetrd");
#else
        LAP_Results(Lap_info,int(TMV_REAL(work[0])),n,n,lwork,"zhetrd");
#endif
        //std::cout<<"Done LapTridiag\n";
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapTridiagonalize(
        SymMatrixView<float> A, VectorView<float> beta,
        VectorView<float> D, VectorView<float> E, float& )
    {
        TMVAssert(A.iscm());
        TMVAssert(A.uplo() == Lower);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(beta.step() == 1);
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size()+1 == A.size());
        TMVAssert(beta.size() == A.size());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        TMVAssert(beta.ct()==NonConj);

        int n = A.size();
        int ldu = A.stepj();
        beta.setZero();
        D.setZero();
        E.setZero();
        int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<float> work(1);
        work.get()[0] = 0.F;
        LAPNAME(ssytrd) (
            LAPCM LAPCH_LO,LAPV(n),
            LAPP(A.ptr()),LAPV(ldu),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
        lwork = int(work[0]);
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(ssytrd) (
            LAPCM LAPCH_LO,LAPV(n),
            LAPP(A.ptr()),LAPV(ldu),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"ssytrd");
#else
        LAP_Results(Lap_info,int(work[0]),n,n,lwork,"ssytrd");
#endif
    }
    template <> 
    void LapTridiagonalize(
        SymMatrixView<std::complex<float> > A,
        VectorView<std::complex<float> > beta,
        VectorView<float> D, VectorView<float> E, 
        std::complex<float>& )
    {
        TMVAssert(A.iscm());
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.isherm());
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(beta.step() == 1);
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size()+1 == A.size());
        TMVAssert(beta.size() == A.size());
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        TMVAssert(beta.ct()==NonConj);

        int n = A.size();
        int ldu = A.stepj();
        beta.setZero();
        D.setZero();
        E.setZero();
        int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<std::complex<float> > work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<std::complex<float> > work(1);
        work.get()[0] = 0.F;
        LAPNAME(chetrd) (
            LAPCM LAPCH_LO,LAPV(n),
            LAPP(A.ptr()),LAPV(ldu),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
        lwork = int(TMV_REAL(work[0]));
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        //std::cout<<"Start LapTridiagonalize:\n";
        //std::cout<<"A = "<<TMV_Text(A)<<std::endl;
        //std::cout<<"beta = "<<TMV_Text(beta)<<"  "<<beta<<std::endl;
        //std::cout<<"D = "<<TMV_Text(D)<<"  "<<D<<std::endl;
        //std::cout<<"E = "<<TMV_Text(E)<<"  "<<E<<std::endl;
        //std::cout<<"n = "<<n<<std::endl;
        //std::cout<<"ldu = "<<ldu<<std::endl;
        //std::cout<<"lwork = "<<lwork<<std::endl;
        //std::cout<<"work = "<<VectorViewOf(work.get(),lwork)<<std::endl;
        //std::cout<<"Before chetrd"<<std::endl;
        LAPNAME(chetrd) (
            LAPCM LAPCH_LO,LAPV(n),
            LAPP(A.ptr()),LAPV(ldu),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
        //std::cout<<"After chetrd"<<std::endl;
        if (!A.isconj()) beta.conjugateSelf();
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"chetrd");
#else
        LAP_Results(Lap_info,int(TMV_REAL(work[0])),n,n,lwork,"chetrd");
#endif
    }
#endif
#endif // LAP

    template <class T, class Td> 
    void Tridiagonalize(
        SymMatrixView<T> A, VectorView<T> beta,
        VectorView<Td> D, VectorView<RT> E, T& signdet)
    {
#ifdef XDEBUG
        cout<<"Start Tridiagonalize: \n";
        cout<<"A = "<<TMV_Text(A)<<endl;
        cout<<"beta = "<<TMV_Text(beta)<<endl;
        cout<<"D = "<<TMV_Text(D)<<endl;
        cout<<"E = "<<TMV_Text(E)<<endl;
        cout<<"signdet = "<<signdet<<endl;
#endif
        TMVAssert(A.size() == D.size());
#ifdef LAP
        TMVAssert(beta.size() == A.size());
#else
        TMVAssert(beta.size() == A.size()-1);
#endif
        TMVAssert(beta.ct()==NonConj);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        TMVAssert(isReal(Td()) || !A.isherm());

        // If Hermitian, then the determinant is unaffected by Householder
        // multiplies, so don't bother keeping track of d.
        T d = A.isherm() ? T(0) : T(1);

        if (A.size() > 0) {
            TMVAssert(E.size()+1 == D.size());
            if (A.uplo() == Upper) {
                if (A.isherm()) Tridiagonalize(A.adjoint(),beta,D,E,d);
                else Tridiagonalize(A.transpose(),beta,D,E,d);
            } else {
#ifdef XDEBUG
                Matrix<T> A0(A);
                Matrix<T,ColMajor> A2(A);
                Vector<T> b2(A.size()-1);
                Vector<Td> D2(D.size());
                Vector<RT> E2(E.size());
                T d2 = d;
                if (A.isherm())
                    NonBlockTridiagonalize(
                        HermMatrixViewOf(A2,Lower),
                        b2.view(),D2.view(),E2.view(),d2);
                else
                    NonBlockTridiagonalize(
                        SymMatrixViewOf(A2,Lower),
                        b2.view(),D2.view(),E2.view(),d2);
#endif // XDEBUG
#ifdef LAP
                // TODO: The RowMajor version of the LAPACK Tridiagonalize 
                // algorithm (Upper ColMajor storage is how they refer to it)
                // stores the Householder vectors in a different way than my
                // GetQFromQR function can handle.  I should really write a 
                // version that can handle that, so I can call LAPACK in
                // that case as well.
                if (A.iscm() && A.isherm()) {
                    TMVAssert(isReal(Td()));
                    if (D.step() != 1) {
                        Vector<Td> Dx(D.size());
                        Tridiagonalize(A,beta,Dx.view(),E,d);
                        D = Dx;
                    } else if (E.step() != 1) {
                        Vector<RT> Ex(E.size());
                        Tridiagonalize(A,beta,D,Ex.view(),d);
                        E = Ex;
                    } else if (beta.step() != 1) {
                        Vector<T> betax(beta.size());
                        Tridiagonalize(A,betax.view(),D,E,d);
                        beta = betax;
                    } else {
                        LapTridiagonalize(A,beta,D,E,d);
                    }
                } else {
                    NonLapTridiagonalize(A,beta.subVector(0,A.size()-1),D,E,d);
                }
#else // LAP
                NonLapTridiagonalize(A,beta,D,E,d);
#endif

#ifdef XDEBUG
                const ptrdiff_t N = A.size();
                Matrix<T> UU(N,N,T(0));
                for(ptrdiff_t j=N-1;j>0;--j) UU.col(j,j,N) = A.col(j-1,j,N);
                UU(0,0) = T(1);
                GetQFromQR(UU.subMatrix(1,N,1,N),beta.subVector(0,N-1));
                Matrix<T> TT(N,N,T(0));
                TT.diag() = D; TT.diag(1) = TT.diag(-1) = E;
                Matrix<T> A3 = UU*TT*(A.isherm() ? UU.adjoint() : UU.transpose());
                if (!(Norm(A3-A0) < 0.001*Norm(A0))) {
                    cerr<<"Tridiagonalize: \n";
                    cerr<<"A0 = "<<TMV_Text(A)<<endl;
                    cerr<<"  "<<A0<<endl;
                    cerr<<"Done: A = "<<A<<endl;
                    cerr<<"beta = "<<beta<<endl;
                    cerr<<"UU = "<<UU<<endl;
                    cerr<<"D = "<<D<<endl;
                    cerr<<"E = "<<E<<endl;
                    cerr<<"TT = "<<TT<<endl;
                    cerr<<"UU * TT * UUt = "<<A3<<endl;
                    cerr<<"Norm(diff) = "<<Norm(A0-A3)<<std::endl;
                    cerr<<"NonBlock versions:\n";
                    if (A.isherm())
                        cerr<<"A2 = "<<HermMatrixViewOf(A2,Lower)<<endl;
                    else
                        cerr<<"A2 = "<<SymMatrixViewOf(A2,Lower)<<endl;
                    cerr<<"cf. A = "<<A<<endl;
                    cerr<<"A2 - A = "<<A2-A<<endl;
                    cerr<<"Norm(A2-A) = "<<Norm(A2-A)<<endl;
                    cerr<<"beta2 = "<<b2<<endl;
                    cerr<<"beta = "<<beta<<endl;
                    cerr<<"Norm(beta2-beta) = "<<
                        Norm(b2-beta.subVector(0,b2.size()))<<endl;
                    cerr<<"D2 = "<<D2<<endl;
                    cerr<<"D = "<<D<<endl;
                    cerr<<"Norm(D2-D) = "<<Norm(D2-D)<<endl;
                    cerr<<"E2 = "<<E2<<endl;
                    cerr<<"E = "<<E<<endl;
                    cerr<<"Norm(E2-E) = "<<Norm(E2-E)<<endl;
                    abort();
                }
#endif // XDEBUG
            }
        }
        if (!A.isherm()) signdet *= d;
    }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymSVDecompose_Tridiag.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


