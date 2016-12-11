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
#include "TMV_SymBandSVDiv.h"
#include "tmv/TMV_SymBandSVD.h"
#include "tmv/TMV_SymBandMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include "TMV_SymSVDiv.h"
#include "TMV_BandSVDiv.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_Householder.h"
#include "tmv/TMV_SymHouseholder.h"

#ifdef NOTHROW
#include <iostream>
#endif

#ifdef XDEBUG
#include "tmv/TMV_DiagMatrixArith.h"
#include "tmv/TMV_SymBandMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif


namespace tmv {

#define RT TMV_RealType(T)

#ifdef TMV_BLOCKSIZE
#define SYM_TRIDIAG_BLOCKSIZE TMV_BLOCKSIZE
#else
#define SYM_TRIDIAG_BLOCKSIZE 64
#endif

    template <class T> 
    static void MakeTridiagReal(
        VectorView<T> Udiag, 
        const GenVector<T>& cD, const GenVector<T>& cE,
        VectorView<T> D, VectorView<T> E)
    {
        TMVAssert(cD.size() == D.size());
        TMVAssert(cE.size() == E.size());
        TMVAssert(D.size() == Udiag.size());
        TMVAssert(E.size() == D.size()-1);

        Udiag.setAllTo(T(1));
        D = cD;
        E = cE;
    } 

    template <class T> 
    static void MakeTridiagReal(
        VectorView<std::complex<T> > Udiag, 
        const GenVector<std::complex<T> >& cD,
        const GenVector<std::complex<T> >& cE,
        VectorView<T> D, VectorView<T> E)
    {
        // The complexity of D determines whether the original 
        // SymBandMatrix was hermitian or symmetric.
        // This one is the hermitian case.
        TMVAssert(cD.size() == D.size());
        TMVAssert(cE.size() == E.size());
        TMVAssert(D.size() == Udiag.size());
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(Udiag.step() == 1);
        TMVAssert(Udiag.ct() == NonConj);
        TMVAssert(cE.ct() == NonConj);

        const ptrdiff_t N = D.size();

        std::complex<T>* Uj = Udiag.ptr();
        *Uj = T(1);
        T* Ej = E.ptr();
        const std::complex<T>* cEj = cE.cptr();
        const ptrdiff_t cEstep = cE.step();

        for(ptrdiff_t j=1;j<N;++j,++Ej,cEj+=cEstep) {
            std::complex<T> xcEj = (*Uj) * (*cEj);
            ++Uj;
#ifdef TMVFLDEBUG
            TMVAssert(Ej >= E._first);
            TMVAssert(Ej < E._last);
            TMVAssert(Uj >= Udiag._first);
            TMVAssert(Uj < Udiag._last);
#endif
            *Ej = TMV_ABS(xcEj);
            *Uj = TMV_SIGN(xcEj,*Ej);
        }
        D = cD.realPart();
    }

    template <class T> 
    static inline void MakeTridiagReal(
        VectorView<std::complex<T> > , 
        const GenVector<std::complex<T> >& ,
        const GenVector<std::complex<T> >& ,
        VectorView<std::complex<T> > , VectorView<T> ) 
    {
        // This one is the symmetric case, which should never get called.
        TMVAssert(TMV_FALSE); 
    }

    template <class T, class Td> 
    static void NonLapTridiagonalize(
        const GenSymBandMatrix<T>& A, MatrixView<T> U,
        VectorView<Td> D, VectorView<RT> E, T& signdet)
    {
        // Decompose A into U T Ut
        // The Tridiagonal Matrix T is stored as two vectors: D, E
        // D is the diagonal, E is the sub-diagonal.
        // If A is herm, E* is the super-diagonal, otherwise E is.
        // However, the Householder reflections make E real, so this
        // distinction is irrelevant.

        TMVAssert(A.nlo() > 0);
        TMVAssert(A.size() > A.nlo());
        TMVAssert(A.uplo() == Lower);
        if (U.cptr()) {
            TMVAssert(U.colsize() == A.colsize());
            TMVAssert(U.rowsize() == A.colsize());
            TMVAssert(U.ct() == NonConj);
        }
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size() == A.size()-1);
        TMVAssert(isReal(Td()) || !A.isherm());

        const ptrdiff_t N = A.size();
        const ptrdiff_t nlo = A.nlo();
        //std::cout<<"Start NonLapTridiagonalize\n";
        //std::cout<<"N,nlo = "<<N<<','<<nlo<<std::endl;

        if (nlo == 1) {
            TMVAssert(A.isherm());
            Vector<T> Ud(N);
            if (A.isconj()) {
                MakeTridiagReal(Ud.view(),A.diag(),A.diag(-1).conjugate(),D,E);
                Ud.conjugateSelf();
            } else {
                MakeTridiagReal(Ud.view(),A.diag(),A.diag(-1),D,E);
            }
            if (U.cptr()) {
                U.setZero();
                U.diag() = Ud;
            }
        } else {
            auto_ptr<SymMatrix<T,Lower|ColMajor> > UUS;
            auto_ptr<HermMatrix<T,Lower|ColMajor> > UUH;
            auto_ptr<SymMatrixView<T> > U1;
            if (U.cptr()) {
                U = A;
                if (A.issym())
                    U1.reset(new SymMatrixView<T>(SymMatrixViewOf(U,Lower)));
                else 
                    U1.reset(new SymMatrixView<T>(HermMatrixViewOf(U,Lower)));
            } else {
                if (A.issym()) {
                    UUS.reset(new SymMatrix<T,Lower|ColMajor>(A));
                    U1.reset(new SymMatrixView<T>(UUS->view()));
                } else {
                    UUH.reset(new HermMatrix<T,Lower|ColMajor>(A));
                    U1.reset(new SymMatrixView<T>(UUH->view()));
                }
            }

            std::vector<ptrdiff_t> vec(N-1);
            Vector<T> Ubeta(N-1);
            ptrdiff_t endcol = nlo+1;
            if (endcol > N) endcol = N;

            T* Ubj = Ubeta.ptr();
            // We use Householder reflections to reduce A to the tridiagonal form:
            for(ptrdiff_t j=0;j<N-1;++j,++Ubj) {
#ifdef TMVFLDEBUG
                TMVAssert(Ubj >= Ubeta._first);
                TMVAssert(Ubj < Ubeta._last);
#endif
                *Ubj = HouseholderReflect(U1->col(j,j+1,endcol),signdet);
                if (endcol < N) {
                    endcol+=nlo;
                    if (endcol > N) endcol = N;
                }
                vec[j] = endcol;
                if (*Ubj != T(0)) 
                    HouseholderLRMult(
                        U1->col(j,j+2,endcol),*Ubj,
                        U1->subSymMatrix(j+1,endcol));
            }

            // The tridiagonal of U1 is the tridiagonal we want, so copy it to D,E
            if (isReal(Td())) D = U1->diag().realPart();
            else D = U1->diag();
            E = U1->diag(-1).realPart();
            if (isComplex(T())) {
                TMVAssert(NormInf(U1->diag(-1).imagPart()) == RT(0));
            }

            if (U.cptr()) {
                U.upperTri().setZero();
                U.diag(-1).setZero();
                for (ptrdiff_t j=N-1;j>0;--j) {
                    endcol = vec[j-1];
                    U.col(j,j+1,endcol) = U.col(j-1,j+1,endcol);
                    HouseholderUnpack(U.subMatrix(j,endcol,j,N),*(--Ubj));
                }
                U.col(0,0,vec[0]).makeBasis(0);
            }

            if (!A.isherm()) signdet *= signdet;
        }
        //std::cout<<"Done tridiag: D = "<<D<<std::endl;
        //std::cout<<"E = "<<E<<std::endl;
    }

#ifdef LAP
    template <class T, class Td> 
    static inline void LapTridiagonalize(
        const GenSymBandMatrix<T>& A, MatrixView<T> U,
        VectorView<Td> D, VectorView<RT> E, T& signdet)
    { NonLapTridiagonalize(A,U,D,E,signdet); }
#ifdef INST_DOUBLE
    template <> 
    void LapTridiagonalize(
        const GenSymBandMatrix<double>& A, MatrixView<double> U,
        VectorView<double> D, VectorView<double> E, double& )
    {
        TMVAssert(A.nlo() > 1);
        TMVAssert(A.size() > A.nlo());
        TMVAssert(A.uplo() == Lower);
        if (U.cptr()) {
            TMVAssert(U.colsize() == A.colsize());
            TMVAssert(U.rowsize() == A.colsize());
            TMVAssert(U.iscm());
            TMVAssert(U.ct() == NonConj);
        }
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size() == A.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);

        int n = A.size();
        int kl = A.nlo();
        SymBandMatrix<double,Lower|ColMajor> A2 = A;
        int lda = A2.diagstep();
        int ldu = U.stepj();
        char vect = U.cptr() ? 'V' : 'N';
        double* UU = U.ptr();
        D.setZero();
        E.setZero();
#ifndef LAPNOWORK
        int lwork = n;
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        int Lap_info=0;
        LAPNAME(dsbtrd) (
            LAPCM LAPV(vect),LAPCH_LO,LAPV(n),LAPV(kl),
            LAPP(A2.cptr()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(UU),LAPV(ldu) LAPWK(work.get()) LAPINFO LAP1 LAP1 );
        LAP_Results(Lap_info,"dsbtrd");
    }
    template <> 
    void LapTridiagonalize(
        const GenSymBandMatrix<std::complex<double> >& A,
        MatrixView<std::complex<double> > U,
        VectorView<double> D, VectorView<double> E, 
        std::complex<double>& )
    {
        TMVAssert(A.isherm());
        TMVAssert(A.nlo() > 1);
        TMVAssert(A.size() > A.nlo());
        TMVAssert(A.uplo() == Lower);
        if (U.cptr()) {
            TMVAssert(U.colsize() == A.colsize());
            TMVAssert(U.rowsize() == A.colsize());
            TMVAssert(U.iscm());
            TMVAssert(U.ct() == NonConj);
        }
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size() == A.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);

        int n = A.size();
        int kl = A.nlo();
        HermBandMatrix<std::complex<double>,Lower|ColMajor> A2 = A;
        int lda = A2.diagstep();
        int ldu = U.stepj();
        char vect = U.cptr() ? 'V' : 'N';
        std::complex<double>* UU = U.ptr();
        D.setZero();
        E.setZero();
#ifndef LAPNOWORK
        int lwork = n;
        AlignedArray<std::complex<double> > work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        int Lap_info=0;
        LAPNAME(zhbtrd) (
            LAPCM LAPV(vect),LAPCH_LO,LAPV(n),LAPV(kl),
            LAPP(A2.cptr()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(UU),LAPV(ldu) LAPWK(work.get()) LAPINFO LAP1 LAP1 );
        LAP_Results(Lap_info,"zhbtrd");
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapTridiagonalize(
        const GenSymBandMatrix<float>& A, MatrixView<float> U,
        VectorView<float> D, VectorView<float> E, float& )
    {
        TMVAssert(A.nlo() > 1);
        TMVAssert(A.size() > A.nlo());
        TMVAssert(A.uplo() == Lower);
        if (U.cptr()) {
            TMVAssert(U.colsize() == A.colsize());
            TMVAssert(U.rowsize() == A.colsize());
            TMVAssert(U.iscm());
            TMVAssert(U.ct() == NonConj);
        }
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size() == A.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);

        int n = A.size();
        int kl = A.nlo();
        SymBandMatrix<float,Lower|ColMajor> A2 = A;
        int lda = A2.diagstep();
        int ldu = U.stepj();
        char vect = U.cptr() ? 'V' : 'N';
        float* UU = U.ptr();
        D.setZero();
        E.setZero();
#ifndef LAPNOWORK
        int lwork = n;
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        int Lap_info=0;
        LAPNAME(ssbtrd) (
            LAPCM LAPV(vect),LAPCH_LO,LAPV(n),LAPV(kl),
            LAPP(A2.cptr()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(UU),LAPV(ldu) LAPWK(work.get()) LAPINFO LAP1 LAP1 );
        LAP_Results(Lap_info,"dsbtrd");
    }
    template <> 
    void LapTridiagonalize(
        const GenSymBandMatrix<std::complex<float> >& A,
        MatrixView<std::complex<float> > U,
        VectorView<float> D, VectorView<float> E, 
        std::complex<float>& )
    {
        TMVAssert(A.isherm());
        TMVAssert(A.nlo() > 1);
        TMVAssert(A.size() > A.nlo());
        TMVAssert(A.uplo() == Lower);
        if (U.cptr()) {
            TMVAssert(U.colsize() == A.colsize());
            TMVAssert(U.rowsize() == A.colsize());
            TMVAssert(U.iscm());
            TMVAssert(U.ct() == NonConj);
        }
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size() == A.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);

        int n = A.size();
        int kl = A.nlo();
        HermBandMatrix<std::complex<float>,Lower|ColMajor> A2 = A;
        int lda = A2.diagstep();
        int ldu = U.stepj();
        char vect = U.cptr() ? 'V' : 'N';
        std::complex<float>* UU = U.ptr();
        D.setZero();
        E.setZero();
#ifndef LAPNOWORK
        int lwork = n;
        AlignedArray<std::complex<float> > work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        int Lap_info=0;
        LAPNAME(chbtrd) (
            LAPCM LAPV(vect),LAPCH_LO,LAPV(n),LAPV(kl),
            LAPP(A2.cptr()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(UU),LAPV(ldu) LAPWK(work.get()) LAPINFO LAP1 LAP1 );
        LAP_Results(Lap_info,"chbtrd");
    }
#endif 
#endif // LAP

    template <class T, class Td> 
    static void Tridiagonalize(
        const GenSymBandMatrix<T>& A, MatrixView<T> U,
        VectorView<Td> D, VectorView<RT> E, T& signdet)
    {
        TMVAssert(A.nlo() > 0);
        TMVAssert(A.size() == D.size());
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        TMVAssert(isReal(Td()) || !A.isherm());
        if (U.cptr()) {
            TMVAssert(U.colsize() == D.size());
            TMVAssert(U.rowsize() == D.size());
        }

#ifdef XDEBUG
        cout<<"Start Tridiagonalize: \n";
        cout<<"A = "<<TMV_Text(A)<<endl;
        if (U.cptr()) cout<<"U = "<<TMV_Text(U)<<endl;
        cout<<"D = "<<TMV_Text(D)<<endl;
        cout<<"E = "<<TMV_Text(E)<<endl;
        cout<<"signdet = "<<signdet<<endl;
        Matrix<T> A0 = A;
#endif

        if (A.size() > 0) {
            TMVAssert(E.size()+1 == D.size());
            if (A.uplo() == Upper) {
                if (A.isherm()) Tridiagonalize(A.adjoint(),U,D,E,signdet);
                else Tridiagonalize(A.transpose(),U,D,E,signdet);
            } else {
#ifdef LAP
                if (A.isherm() && A.nlo() > 1 && (!U.cptr() || U.iscm())) {
                    TMVAssert(isReal(Td()));
                    if (D.step() != 1) {
                        Vector<Td> Dx(D.size());
                        Tridiagonalize(A,U,Dx.view(),E,signdet);
                        D = Dx;
                    } else if (E.step() != 1) {
                        Vector<RT> Ex(E.size());
                        Tridiagonalize(A,U,D,Ex.view(),signdet);
                        E = Ex;
                    } else {
                        LapTridiagonalize(A,U,D,E,signdet);
                    }
                }
                else 
#endif // LAP
                    NonLapTridiagonalize(A,U,D,E,signdet);
            }
        }

#ifdef XDEBUG
        if (U.cptr()) {
            const ptrdiff_t N = A.size();
            Matrix<T> TT(N,N,T(0));
            TT.diag() = D;
            TT.diag(1) = TT.diag(-1) = E;
            Matrix<T> A2 = U*TT*(A.isherm() ? U.adjoint() : U.transpose());
            std::cout<<"After Tridiag: Norm(A2-A0) = "<<Norm(A2-A0)<<std::endl;
            std::cout<<"Norm(A0) = "<<Norm(A0)<<std::endl;
            std::cout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
            std::cout<<"Norm(UUt-1) = "<<Norm(U*U.adjoint()-T(1))<<std::endl;
            if (!(Norm(A2-A0) <= 1.e-6*Norm(A0))) {
                cerr<<"Tridiagonalize: \n";
                cerr<<"A0 = "<<TMV_Text(A)<<"  "<<A0<<endl;
                cerr<<"Done: U = "<<U<<endl;
                cerr<<"D = "<<D<<endl;
                cerr<<"E = "<<E<<endl;
                cerr<<"TT = "<<TT<<endl;
                cerr<<"UU * TT * UUt = "<<A2.maxAbs2Element()*1.e-3<<endl;
                cerr<<"A0 = "<<A0<<endl;
                cerr<<"A2-A0 = "<<A2-A0<<endl;
                cerr<<"Norm(A2-A0) = "<<Norm(A2-A0)<<endl;
                cerr<<"Norm(A0) = "<<Norm(A0)<<endl;
                abort();
            }
        }
#endif // XDEBUG
    }

    template <class T> 
    void UnsortedEigen(
        const GenSymBandMatrix<T>& A, MatrixView<T> U, VectorView<RT> SS)
    {
        TMVAssert(A.size() > 0);
        TMVAssert(A.isherm());
        if (U.cptr()) {
            TMVAssert(U.rowsize() == A.size());
            TMVAssert(U.colsize() == A.size());
            TMVAssert(U.ct() == NonConj);
        }
        TMVAssert(SS.size() == A.size());
        TMVAssert(SS.ct() == NonConj);

#ifdef XDEBUG
        Matrix<T> A0 = A;
#endif
        // Decompose Hermitian A (input as lower tri of U) into U S Ut
        // where S is a diagonal real matrix, and U is a unitary matrix.
        // U,S are N x N
        const ptrdiff_t N = A.size();
        if (N == 0) return;

        if (A.nlo() == 0) {
            SS = A.diag().realPart();
            if (U.cptr()) U.setToIdentity();
        } else {
            // First we reduce A to tridiagonal form: A = U * T * Ut
            // using a series of Householder transformations.
            // The diagonal of the Tridiagonal Matrix T is stored in D.
            // The subdiagonal is stored in E.
            T d2 = 0;
            Vector<RT> E(N-1);
            Tridiagonalize(A,U,SS,E.view(),d2);

            // Then finish the decomposition as for a HermMatrix
            EigenFromTridiagonal(U,SS,E.view());
        }

#ifdef XDEBUG
        if (U.cptr()) {
            Matrix<T> A2 = U * DiagMatrixViewOf(SS) * U.adjoint();
            std::cout<<"After UnsortedEigen: Norm(A2-A0) = "<<Norm(A2-A0)<<std::endl;
            std::cout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
            std::cout<<"Norm(UUt-1) = "<<Norm(U*U.adjoint()-T(1))<<std::endl;
            if (!(Norm(A0-A2) <= 0.0001 * Norm(U) * Norm(SS) * Norm(U))) {
                cerr<<"Unsorted Eigen:\n";
                //cerr<<"A = "<<A0<<endl;
                //cerr<<"U = "<<U<<endl;
                cerr<<"S = "<<SS<<endl;
                cerr<<"A-USUt = "<<Matrix<T>(A0-A2).clip(1.e-5)<<endl;
                cerr<<"Norm(A0-A2) = "<<Norm(A0-A2)<<endl;
                cerr<<"Norm(U) = "<<Norm(U)<<endl;
                cerr<<"Norm(SS) = "<<Norm(SS)<<endl;
                cerr<<"Norm(A0) = "<<Norm(A0)<<endl;
                abort();
            }
        }
#endif
    }

    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        MatrixView<T> U, DiagMatrixView<RT> SS, MatrixView<T> Vt, 
        RT& logdet, T& signdet)
    {
        TMVAssert(A.size() > 0);
        if (U.cptr()) {
            TMVAssert(U.rowsize() == A.size());
            TMVAssert(U.colsize() == A.size());
            TMVAssert(U.ct() == NonConj);
        }
        if (Vt.cptr()) {
            TMVAssert(Vt.rowsize() == A.size());
            TMVAssert(Vt.colsize() == A.size());
            TMVAssert(Vt.ct() == NonConj);
        }
        TMVAssert(SS.size() == A.size());

#ifdef XDEBUG
        Matrix<T> A0 = A;
#endif

        if (A.isherm()) {
            if (Vt.cptr() && !U.cptr()) {
                UnsortedEigen<T>(A,Vt.transpose(),SS.diag());
                Vt.conjugateSelf();
            } else  {
                UnsortedEigen<T>(A,U,SS.diag());
            }
            if (Vt.cptr() && U.cptr()) Vt = U.adjoint();
            if (signdet != T(0)) {
                RT s;
                logdet += SS.logDet(&s);
                signdet *= s;
            }
            if (U.cptr() || Vt.cptr()) {
                AlignedArray<ptrdiff_t> sortp(A.size());
                SS.diag().sort(sortp.get(),Descend,AbsComp);
                if (U.cptr()) U.permuteCols(sortp.get());
                if (Vt.cptr()) Vt.permuteRows(sortp.get());
            } else {
                SS.diag().sort(Descend,AbsComp);
            }
        } else {
            TMVAssert(isComplex(T()));
            // Decompose complex symmetric A (input as lower tri of U) into 
            // U S Vt where S is a diagonal real matrix, and U,Vt are
            // unitary matrices.
            // U,S,Vt are N x N
            // If Vt = 0, then U,Vt are not formed.
            // Only S,det are accurate on return.
            const ptrdiff_t N = A.size();
            if (N == 0) return;

            if (A.nlo() == 0) {
                SV_Decompose(A.diagRange(0,1),U,SS,Vt,logdet,signdet);
            } else if (A.nlo() == 1) {
                BandMatrix<T,ColMajor> B = A;
                SV_Decompose(B,U,SS,Vt,logdet,signdet);
            } else {
                // First we reduce A to tridiagonal form: A = U * T * UT
                // using a series of Householder transformations.
                // The diagonal of the Tridiagonal Matrix T is stored in D.
                // The subdiagonal is stored in E.
                Vector<T> D(N);
                Vector<RT> E(N-1);
                if (Vt.cptr() && !U.cptr()) {
                    Tridiagonalize<T>(
                        A,Vt.transpose(),D.view(),E.view(),signdet);
                } else {
                    Tridiagonalize<T>(A,U,D.view(),E.view(),signdet);
                    if (Vt.cptr()) { TMVAssert(U.cptr()); Vt = U.transpose(); }
                }

                BandMatrix<T,ColMajor> B(N,N,1,1);
                B.diag() = D;
                B.diag(-1) = E;
                B.diag(1) = E;

                if (U.cptr()) {
                    Matrix<T,ColMajor> U1(N,N);
                    if (Vt.cptr()) {
                        Matrix<T,ColMajor> Vt1(N,N);
                        SV_Decompose<T>(B,U1.view(),SS,Vt1.view(),logdet,signdet);
                        Vt = Vt1*Vt;
                    } else {
                        SV_Decompose<T>(B,U1.view(),SS,Vt,logdet,signdet);
                    }
                    U = U*U1;
                } else {
                    if (Vt.cptr()) {
                        Matrix<T,ColMajor> Vt1(N,N);
                        SV_Decompose<T>(B,U,SS,Vt1.view(),logdet,signdet);
                        Vt = Vt1*Vt;
                    } else {
                        SV_Decompose<T>(B,U,SS,Vt,logdet,signdet);
                    }
                }
            }
        }
#ifdef XDEBUG
        if (U.cptr()&&Vt.cptr()) {
            Matrix<T> A2 = U * SS * Vt;
            if (!(Norm(A0-A2) <= 0.0001 * Norm(U) * Norm(SS) * Norm(Vt))) {
                cerr<<"SV_Decompose:\n";
                cerr<<"A = "<<A0<<endl;
                cerr<<"U = "<<U<<endl;
                cerr<<"S = "<<SS<<endl;
                cerr<<"Vt = "<<Vt<<endl;
                cerr<<"USVt = "<<A2<<endl;
                abort();
            }
        }
#endif
    }

    template <class T> 
    void Eigen(
        const GenSymBandMatrix<T>& A, MatrixView<T> U,
        VectorView<RT> SS)
    {
        TMVAssert(SS.size() == A.size());
        TMVAssert(U.colsize() == A.size());
        TMVAssert(U.rowsize() == A.size());

        if (A.isconj()) {
            if (U.isconj()) {
                Eigen(A.conjugate(),U.conjugate(),SS);
            } else {
                Eigen(A.conjugate(),U,SS);
                U.conjugateSelf();
            }
        } else {
            if (U.isconj()) {
                Eigen(A,U.conjugate(),SS);
                U.conjugateSelf();
            } else {
                UnsortedEigen<T>(A,U,SS);
                AlignedArray<ptrdiff_t> sortp(A.size());
                SS.sort(sortp.get(),Ascend);
                U.permuteCols(sortp.get());
            }
        }
    }

    template <class T> 
    void Eigen(const GenSymBandMatrix<T>& A, VectorView<RT> SS)
    {
        TMVAssert(SS.size() == A.size());
        if (A.isconj()) {
            Eigen(A.conjugate(),SS);
        } else {
            MatrixView<T> U(0,0,0,1,1,NonConj);
            UnsortedEigen<T>(A,U,SS);
            SS.sort(Ascend);
        }
    }

    // Decompose A into U S V
    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        MatrixView<T> U, DiagMatrixView<RT> SS,
        MatrixView<T> Vt)
    {
        TMVAssert(U.colsize() == A.size());
        TMVAssert(U.rowsize() == A.size());
        TMVAssert(SS.size() == A.size());
        TMVAssert(Vt.colsize() == A.size());
        TMVAssert(Vt.rowsize() == A.size());

        if (A.isconj()) {
            if (U.isconj()) {
                if (Vt.isconj()) {
                    SV_Decompose(A.conjugate(),U.conjugate(),SS,Vt.conjugate());
                } else {
                    SV_Decompose(A.conjugate(),U.conjugate(),SS,Vt);
                    Vt.conjugateSelf();
                }
            } else {
                if (Vt.isconj()) {
                    SV_Decompose(A.conjugate(),U,SS,Vt.conjugate());
                    U.conjugateSelf();
                } else {
                    SV_Decompose(A.conjugate(),U,SS,Vt);
                    U.conjugateSelf();
                    Vt.conjugateSelf();
                }
            }
        } else {
            if (U.isconj()) {
                if (Vt.isconj()) {
                    SV_Decompose(A,U.conjugate(),SS,Vt.conjugate());
                    U.conjugateSelf();
                    Vt.conjugateSelf();
                } else {
                    SV_Decompose(A,U.conjugate(),SS,Vt);
                    U.conjugateSelf();
                }
            } else {
                if (Vt.isconj()) {
                    SV_Decompose(A,U,SS,Vt.conjugate());
                    Vt.conjugateSelf();
                } else {
                    RT ld(0);
                    T d(0);
                    SV_Decompose<T>(A,U,SS,Vt,ld,d);
                    if (A.isherm()) {
                        // Then S values might be negative:
                        for(ptrdiff_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) {
                            SS(i) = -SS(i);
                            Vt.row(i) = -Vt.row(i);
                        }
                    }
                }
            }
        }
    }

    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        MatrixView<T> U, DiagMatrixView<RT> SS)
    {
        TMVAssert(U.colsize() == A.size());
        TMVAssert(U.rowsize() == A.size());
        TMVAssert(SS.size() == A.size());

        if (A.isconj()) {
            if (U.isconj()) {
                SV_Decompose(A.conjugate(),U.conjugate(),SS);
            } else {
                SV_Decompose(A.conjugate(),U,SS);
                U.conjugateSelf();
            }
        } else {
            if (U.isconj()) {
                SV_Decompose(A,U.conjugate(),SS);
                U.conjugateSelf();
            } else {
                RT ld(0);
                T d(0);
                MatrixView<T> Vt(0,0,0,1,1,NonConj);
                SV_Decompose<T>(A,U,SS,Vt,ld,d);
                if (A.isherm()) {
                    for(ptrdiff_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) 
                        SS(i) = -SS(i);
                }
            }
        }
    }

    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A, DiagMatrixView<RT> SS, MatrixView<T> Vt)
    {
        TMVAssert(SS.size() == A.size());
        TMVAssert(Vt.colsize() == A.size());
        TMVAssert(Vt.rowsize() == A.size());

        if (A.isherm()) SV_Decompose(A,Vt.adjoint(),SS);
        else SV_Decompose(A,Vt.transpose(),SS);
    }

    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A, DiagMatrixView<RT> SS)
    {
        TMVAssert(SS.size() == A.size());

        if (A.isconj()) {
            SV_Decompose(A.conjugate(),SS);
        } else {
            RT ld(0);
            T d(0);
            MatrixView<T> U(0,0,0,1,1,NonConj);
            MatrixView<T> Vt(0,0,0,1,1,NonConj);
            SV_Decompose<T>(A,U,SS,Vt,ld,d);
            if (A.isherm()) 
                for(ptrdiff_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) 
                    SS(i) = -SS(i);
        }
    }

    template <class T> 
    void SquareRoot(const GenSymBandMatrix<T>& A, SymMatrixView<T> SS)
    {
        TMVAssert(A.isherm());
        TMVAssert(SS.isherm());
        // A -> A^1/2
        //
        // There are faster algorithms than this, but for now it works.
        //
        // A = V D Vt
        // A = V D^1/2 Vt
        if (A.isconj()) {
            if (SS.isconj()) {
                SquareRoot(A.conjugate(),SS.conjugate());
            } else {
                SquareRoot(A.conjugate(),SS);
                SS.conjugateSelf();
            }
        } else {
            if (SS.isconj()) {
                SquareRoot(A,SS.conjugate());
                SS.conjugateSelf();
            } else {
                Matrix<T> V(A.size(),A.size());
                DiagMatrix<RT> D(A.size());
                Eigen(A,V.view(),D.diag());
                for(ptrdiff_t i=0;i<A.size();i++) {
                    if (D(i) < RT(0))  {
#ifdef NOTHROW
                        std::cerr<<"Non Posdef HermBandMatrix found in "
                            "SqaureRoot\n"; 
                        exit(1); 
#else
                        throw NonPosDef("in SymBandMatrix SquareRoot");
#endif
                    }
                    D(i) = TMV_SQRT(D(i));
                }
                Matrix<T> DVt = D*V.adjoint();
                SymMultMM<false>(T(1),V,DVt,SS);
            }
        }
    }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymBandSVDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


