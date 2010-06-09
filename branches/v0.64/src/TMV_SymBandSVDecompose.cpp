///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
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
#include "TMV_Householder.h"
#include "TMV_SymHouseholder.h"

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
        const VectorView<T>& Udiag, 
        const GenVector<T>& cD, const GenVector<T>& cE,
        const VectorView<T>& D, const VectorView<T>& E)
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
        const VectorView<std::complex<T> >& Udiag, 
        const GenVector<std::complex<T> >& cD,
        const GenVector<std::complex<T> >& cE,
        const VectorView<T>& D, const VectorView<T>& E)
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

        const int N = D.size();

        std::complex<T>* Uj = Udiag.ptr();
        *Uj = T(1);
        T* Ej = E.ptr();
        const std::complex<T>* cEj = cE.cptr();
        const int cEstep = cE.step();

        for(int j=1;j<N;++j,++Ej,cEj+=cEstep) {
            std::complex<T> xcEj = (*Uj) * (*cEj);
            ++Uj;
#ifdef TMVFLDEBUG
            TMVAssert(Ej >= E.first);
            TMVAssert(Ej < E.last);
            TMVAssert(Uj >= Udiag.first);
            TMVAssert(Uj < Udiag.last);
#endif
            *Ej = TMV_ABS(xcEj);
            *Uj = TMV_SIGN(xcEj,*Ej);
        }
        D = cD.realPart();
    }

    template <class T> 
    static inline void MakeTridiagReal(
        const VectorView<std::complex<T> >& , 
        const GenVector<std::complex<T> >& ,
        const GenVector<std::complex<T> >& ,
        const VectorView<std::complex<T> >& , const VectorView<T>& ) 
    {
        // This one is the symmetric case, which should never get called.
        TMVAssert(TMV_FALSE); 
    }

    template <class T, class Td> 
    static void NonLapTridiagonalize(
        const GenSymBandMatrix<T>& A, MVP<T> U,
        const VectorView<Td>& D, const VectorView<RT>& E, T& signdet)
    {
        // Decompose A into U T Ut
        // The Tridiagonal Matrix T is stored as two vectors: D, E
        // D is the diagonal, E is the sub-diagonal.
        // If A is herm, E* is the super-diagonal, otherwise E is.
        // However, the Householder reflections make E real, so this
        // distinction is irrelevant.

        TMVAssert(A.nlo() > 0);
        TMVAssert(int(A.size()) > A.nlo());
        TMVAssert(A.uplo() == Lower);
        if (U) {
            TMVAssert(U->colsize() == A.colsize());
            TMVAssert(U->rowsize() == A.colsize());
            TMVAssert(U->ct() == NonConj);
        }
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size() == A.size()-1);
        TMVAssert(isReal(Td()) || !A.isherm());

        const int N = A.size();
        const int nlo = A.nlo();

        if (nlo == 1) {
            TMVAssert(A.isherm());
            Vector<T> Ud(N);
            if (A.isconj()) {
                MakeTridiagReal(Ud.view(),A.diag(),A.diag(-1).conjugate(),D,E);
                Ud.conjugateSelf();
            } else {
                MakeTridiagReal(Ud.view(),A.diag(),A.diag(-1),D,E);
            }
            if (U) {
                U->setZero();
                U->diag() = Ud;
            }
        } else {
            auto_ptr<SymMatrix<T,Lower,ColMajor> > UUS(0);
            auto_ptr<HermMatrix<T,Lower,ColMajor> > UUH(0);
            auto_ptr<SymMatrixView<T> > U1(0);
            if (U) {
                *U = A;
                if (A.issym())
                    U1.reset(new SymMatrixView<T>(SymMatrixViewOf(*U,Lower)));
                else
                    U1.reset(new SymMatrixView<T>(HermMatrixViewOf(*U,Lower)));
            } else {
                if (A.issym()) {
                    UUS.reset(new SymMatrix<T,Lower,ColMajor>(A));
                    U1.reset(new SymMatrixView<T>(UUS->view()));
                } else {
                    UUH.reset(new HermMatrix<T,Lower,ColMajor>(A));
                    U1.reset(new SymMatrixView<T>(UUH->view()));
                }
            }

            std::vector<int> vec(N-1);
            Vector<T> Ubeta(N-1);
            int endcol = nlo+1;
            if (endcol > N) endcol = N;

            T* Ubj = Ubeta.ptr();
            // We use Householder reflections to reduce A to the tridiagonal form:
            for(int j=0;j<N-1;++j,++Ubj) {
#ifdef TMVFLDEBUG
                TMVAssert(Ubj >= Ubeta.first);
                TMVAssert(Ubj < Ubeta.last);
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
#ifdef XTEST
            if (isComplex(T())) {
                if (isReal(Td()))
                    TMVAssert(normInf(U1->diag().imagPart()) == RT(0));
                TMVAssert(normInf(U1->diag(-1).imagPart()) == RT(0));
            }
#endif

            if (U) {
                U->upperTri().setZero();
                U->diag(-1).setZero();
                for (int j=N-1;j>0;--j) {
                    endcol = vec[j-1];
                    U->col(j,j+1,endcol) = U->col(j-1,j+1,endcol);
                    HouseholderUnpack(U->subMatrix(j,endcol,j,N),*(--Ubj));
                }
                U->col(0,0,vec[0]).makeBasis(0);
            }

            if (!A.isherm()) signdet *= signdet;
        }
    }

#ifdef LAP
    template <class T, class Td> 
    static inline void LapTridiagonalize(
        const GenSymBandMatrix<T>& A, MVP<T> U,
        const VectorView<Td>& D, const VectorView<RT>& E, T& signdet)
    { NonLapTridiagonalize(A,U,D,E,signdet); }
#ifdef INST_DOUBLE
    template <> 
    void LapTridiagonalize(
        const GenSymBandMatrix<double>& A, MVP<double> U,
        const VectorView<double>& D, const VectorView<double>& E, double& )
    {
        TMVAssert(A.nlo() > 1);
        TMVAssert(int(A.size()) > A.nlo());
        TMVAssert(A.uplo() == Lower);
        if (U) {
            TMVAssert(U->colsize() == A.colsize());
            TMVAssert(U->rowsize() == A.colsize());
            TMVAssert(U->iscm());
            TMVAssert(U->ct() == NonConj);
        }
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size() == A.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);

        int n = A.size();
        int kl = A.nlo();
        SymBandMatrix<double,Lower,ColMajor> A2 = A;
        int lda = A2.diagstep();
        int ldu = U ? U->stepj() : 1;
        char vect = U ? 'V' : 'N';
        double* UU = U ? U->ptr() : 0;
#ifndef LAPNOWORK
        int lwork = n;
        auto_array<double> work(new double[lwork]);
#endif
        LAPNAME(dsbtrd) (
            LAPCM LAPV(vect),LAPCH_LO,LAPV(n),LAPV(kl),
            LAPP(A2.cptr()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(UU),LAPV(ldu) LAPWK(work.get()) LAPINFO LAP1 LAP1 );
        LAP_Results("dsbtrd");
    }
    template <> 
    void LapTridiagonalize(
        const GenSymBandMatrix<std::complex<double> >& A,
        MVP<std::complex<double> > U,
        const VectorView<double>& D, const VectorView<double>& E, 
        std::complex<double>& )
    {
        TMVAssert(A.isherm());
        TMVAssert(A.nlo() > 1);
        TMVAssert(int(A.size()) > A.nlo());
        TMVAssert(A.uplo() == Lower);
        if (U) {
            TMVAssert(U->colsize() == A.colsize());
            TMVAssert(U->rowsize() == A.colsize());
            TMVAssert(U->iscm());
            TMVAssert(U->ct() == NonConj);
        }
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size() == A.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);

        int n = A.size();
        int kl = A.nlo();
        HermBandMatrix<std::complex<double>,Lower,ColMajor> A2 = A;
        int lda = A2.diagstep();
        int ldu = U ? U->stepj() : 1;
        char vect = U ? 'V' : 'N';
        std::complex<double>* UU = U ? U->ptr() : 0;
#ifndef LAPNOWORK
        int lwork = n;
        auto_array<std::complex<double> > work(new std::complex<double>[lwork]);
#endif
        LAPNAME(zhbtrd) (
            LAPCM LAPV(vect),LAPCH_LO,LAPV(n),LAPV(kl),
            LAPP(A2.cptr()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(UU),LAPV(ldu) LAPWK(work.get()) LAPINFO LAP1 LAP1 );
        LAP_Results("zhbtrd");
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapTridiagonalize(
        const GenSymBandMatrix<float>& A, MVP<float> U,
        const VectorView<float>& D, const VectorView<float>& E, float& )
    {
        TMVAssert(A.nlo() > 1);
        TMVAssert(int(A.size()) > A.nlo());
        TMVAssert(A.uplo() == Lower);
        if (U) {
            TMVAssert(U->colsize() == A.colsize());
            TMVAssert(U->rowsize() == A.colsize());
            TMVAssert(U->iscm());
            TMVAssert(U->ct() == NonConj);
        }
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size() == A.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);

        int n = A.size();
        int kl = A.nlo();
        SymBandMatrix<float,Lower,ColMajor> A2 = A;
        int lda = A2.diagstep();
        int ldu = U ? U->stepj() : 1;
        char vect = U ? 'V' : 'N';
        float* UU = U ? U->ptr() : 0;
#ifndef LAPNOWORK
        int lwork = n;
        auto_array<float> work(new float[lwork]);
#endif
        LAPNAME(ssbtrd) (
            LAPCM LAPV(vect),LAPCH_LO,LAPV(n),LAPV(kl),
            LAPP(A2.cptr()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(UU),LAPV(ldu) LAPWK(work.get()) LAPINFO LAP1 LAP1 );
        LAP_Results("dsbtrd");
    }
    template <> 
    void LapTridiagonalize(
        const GenSymBandMatrix<std::complex<float> >& A,
        MVP<std::complex<float> > U,
        const VectorView<float>& D, const VectorView<float>& E, 
        std::complex<float>& )
    {
        TMVAssert(A.isherm());
        TMVAssert(A.nlo() > 1);
        TMVAssert(int(A.size()) > A.nlo());
        TMVAssert(A.uplo() == Lower);
        if (U) {
            TMVAssert(U->colsize() == A.colsize());
            TMVAssert(U->rowsize() == A.colsize());
            TMVAssert(U->iscm());
            TMVAssert(U->ct() == NonConj);
        }
        TMVAssert(D.size() == A.size());
        TMVAssert(E.size() == A.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);

        int n = A.size();
        int kl = A.nlo();
        HermBandMatrix<std::complex<float>,Lower,ColMajor> A2 = A;
        int lda = A2.diagstep();
        int ldu = U ? U->stepj() : 1;
        char vect = U ? 'V' : 'N';
        std::complex<float>* UU = U ? U->ptr() : 0;
#ifndef LAPNOWORK
        int lwork = n;
        auto_array<std::complex<float> > work(new std::complex<float>[lwork]);
#endif
        LAPNAME(chbtrd) (
            LAPCM LAPV(vect),LAPCH_LO,LAPV(n),LAPV(kl),
            LAPP(A2.cptr()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
            LAPP(UU),LAPV(ldu) LAPWK(work.get()) LAPINFO LAP1 LAP1 );
        LAP_Results("chbtrd");
    }
#endif 
#endif // LAP

    template <class T, class Td> 
    static void Tridiagonalize(
        const GenSymBandMatrix<T>& A, MVP<T> U,
        const VectorView<Td>& D, const VectorView<RT>& E, T& signdet)
    {
        TMVAssert(A.nlo() > 0);
        TMVAssert(A.size() == D.size());
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.ct()==NonConj);
        TMVAssert(E.ct()==NonConj);
        TMVAssert(isReal(Td()) || !A.isherm());
        if (U) {
            TMVAssert(U->colsize() == D.size());
            TMVAssert(U->rowsize() == D.size());
        }

#ifdef XDEBUG
        Matrix<T> A0 = A;
#endif

        if (A.size() > 0) {
            TMVAssert(E.size()+1 == D.size());
            if (A.uplo() == Upper) {
                if (A.isherm()) Tridiagonalize(A.adjoint(),U,D,E,signdet);
                else Tridiagonalize(A.transpose(),U,D,E,signdet);
            } else {
#ifdef LAP
                if (A.isherm() && A.nlo() > 1 && (!U || U->iscm())) {
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
        if (U) {
            const int N = A.size();
            Matrix<T> TT(N,N,T(0));
            TT.diag() = D;
            TT.diag(1) = TT.diag(-1) = E;
            Matrix<T> A2 = (*U)*TT*(A.isherm() ? U->adjoint() : U->transpose());
            std::cout<<"After Tridiag: Norm(A2-A0) = "<<Norm(A2-A0)<<std::endl;
            std::cout<<"Norm(UtU-1) = "<<Norm(U->adjoint()*(*U)-T(1))<<std::endl;
            std::cout<<"Norm(UUt-1) = "<<Norm((*U)*U->adjoint()-T(1))<<std::endl;
            if (Norm(A2-A0) > 0.001*Norm(A0)) {
                cerr<<"Tridiagonalize: \n";
                cerr<<"A0 = "<<TMV_Text(A)<<"  "<<A0<<endl;
                cerr<<"Done: U = "<<*U<<endl;
                cerr<<"D = "<<D<<endl;
                cerr<<"E = "<<E<<endl;
                cerr<<"TT = "<<TT<<endl;
                cerr<<"UU * TT * UUt = ";
                A2.write(cerr,RT(1.e-12));
                cerr<<endl;
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
        const GenSymBandMatrix<T>& A, MVP<T> U, const VectorView<RT>& SS)
    {
        TMVAssert(A.size() > 0);
        TMVAssert(A.isherm());
        if (U) {
            TMVAssert(U->rowsize() == A.size());
            TMVAssert(U->colsize() == A.size());
            TMVAssert(U->ct() == NonConj);
        }
        TMVAssert(SS.size() == A.size());
        TMVAssert(SS.ct() == NonConj);

#ifdef XDEBUG
        Matrix<T> A0 = A;
#endif
        // Decompose Hermitian A (input as lower tri of U) into U S Ut
        // where S is a diagonal real matrix, and U is a unitary matrix.
        // U,S are N x N
        const int N = A.size();
        if (N == 0) return;

        if (A.nlo() == 0) {
            SS = A.diag().realPart();
            if (U) U->setToIdentity();
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
        if (U) {
            Matrix<T> A2 = (*U) * DiagMatrixViewOf(SS) * U->adjoint();
            std::cout<<"After UnsortedEigen: Norm(A2-A0) = "<<Norm(A2-A0)<<std::endl;
            std::cout<<"Norm(UtU-1) = "<<Norm(U->adjoint()*(*U)-T(1))<<std::endl;
            std::cout<<"Norm(UUt-1) = "<<Norm((*U)*U->adjoint()-T(1))<<std::endl;
            if (Norm(A0-A2) > 0.0001 * Norm(*U) * Norm(SS) * Norm(*U)) {
                cerr<<"Unsorted Eigen:\n";
                //cerr<<"A = "<<A0<<endl;
                //cerr<<"U = "<<*U<<endl;
                cerr<<"S = "<<SS<<endl;
                cerr<<"A-USUt = "<<Matrix<T>(A0-A2).clip(1.e-5)<<endl;
                cerr<<"Norm(A0-A2) = "<<Norm(A0-A2)<<endl;
                cerr<<"Norm(U) = "<<Norm(*U)<<endl;
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
        MVP<T> U, const DiagMatrixView<RT>& SS, MVP<T> V, 
        RT& logdet, T& signdet)
    {
        TMVAssert(A.size() > 0);
        if (U) {
            TMVAssert(U->rowsize() == A.size());
            TMVAssert(U->colsize() == A.size());
            TMVAssert(U->ct() == NonConj);
        }
        if (V) {
            TMVAssert(V->rowsize() == A.size());
            TMVAssert(V->colsize() == A.size());
            TMVAssert(V->ct() == NonConj);
        }
        TMVAssert(SS.size() == A.size());

#ifdef XDEBUG
        Matrix<T> A0 = A;
#endif

        if (A.isherm()) {
            if (V && !U) {
                UnsortedEigen<T>(A,V->transpose(),SS.diag());
                V->conjugateSelf();
            } else  {
                UnsortedEigen<T>(A,U,SS.diag());
            }
            if (V && U) *V = U->adjoint();
            if (signdet != T(0)) {
                RT s;
                logdet += SS.logDet(&s);
                signdet *= s;
            }
            if (U || V) {
                auto_array<int> sortp(new int[A.size()]);
                SS.diag().sort(sortp.get(),Descend,AbsComp);
                if (U) U->permuteCols(sortp.get());
                if (V) V->permuteRows(sortp.get());
            } else {
                SS.diag().sort(Descend,AbsComp);
            }
        } else {
            TMVAssert(isComplex(T()));
            // Decompose complex symmetric A (input as lower tri of U) into 
            // U S V where S is a diagonal real matrix, and U,V are
            // unitary matrices.
            // U,S,V are N x N
            // If V = 0, then U,V are not formed.
            // Only S,det are accurate on return.
            const int N = A.size();
            if (N == 0) return;

            if (A.nlo() == 0) {
                SV_Decompose(A.diagRange(0,1),U,SS,V,logdet,signdet);
            } else if (A.nlo() == 1) {
                BandMatrix<T,ColMajor> B = A;
                SV_Decompose(B,U,SS,V,logdet,signdet);
            } else {
                // First we reduce A to tridiagonal form: A = U * T * UT
                // using a series of Householder transformations.
                // The diagonal of the Tridiagonal Matrix T is stored in D.
                // The subdiagonal is stored in E.
                Vector<T> D(N);
                Vector<RT> E(N-1);
                if (V && !U) {
                    Tridiagonalize<T>(
                        A,V->transpose(),D.view(),E.view(),signdet);
                } else {
                    Tridiagonalize<T>(A,U,D.view(),E.view(),signdet);
                    if (V) { TMVAssert(U); *V = U->transpose(); }
                }

                BandMatrix<T,ColMajor> B(N,N,1,1);
                B.diag() = D;
                B.diag(-1) = E;
                B.diag(1) = E;

                if (U) {
                    Matrix<T,ColMajor> U1(N,N);
                    if (V) {
                        Matrix<T,ColMajor> V1(N,N);
                        SV_Decompose<T>(B,U1.view(),SS,V1.view(),logdet,signdet);
                        *V = V1*(*V);
                    } else {
                        SV_Decompose<T>(B,U1.view(),SS,0,logdet,signdet);
                    }
                    *U = *U*U1;
                } else {
                    if (V) {
                        Matrix<T,ColMajor> V1(N,N);
                        SV_Decompose<T>(B,0,SS,V1.view(),logdet,signdet);
                        *V = V1*(*V);
                    } else {
                        SV_Decompose<T>(B,0,SS,0,logdet,signdet);
                    }
                }
            }
        }
#ifdef XDEBUG
        if (U&&V) {
            Matrix<T> A2 = (*U) * SS * (*V);
            if (Norm(A0-A2) > 0.0001 * Norm(*U) * Norm(SS) * Norm(*V)) {
                cerr<<"SV_Decompose:\n";
                cerr<<"A = "<<A0<<endl;
                cerr<<"U = "<<U<<endl;
                cerr<<"S = "<<SS<<endl;
                cerr<<"V = "<<*V<<endl;
                cerr<<"USV = "<<A2<<endl;
                abort();
            }
        }
#endif
    }

    template <class T> 
    void Eigen(
        const GenSymBandMatrix<T>& A, const MatrixView<T>& U,
        const VectorView<RT>& SS)
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
                auto_array<int> sortp(new int[A.size()]);
                SS.sort(sortp.get(),Ascend);
                U.permuteCols(sortp.get());
            }
        }
    }

    template <class T> 
    void Eigen(
        const GenSymBandMatrix<T>& A, const VectorView<RT>& SS)
    {
        TMVAssert(SS.size() == A.size());
        if (A.isconj()) {
            Eigen(A.conjugate(),SS);
        } else {
            UnsortedEigen<T>(A,0,SS);
            SS.sort(Ascend);
        }
    }

    // Decompose A into U S V
    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        const MatrixView<T>& U, const DiagMatrixView<RT>& SS,
        const MatrixView<T>& V)
    {
        TMVAssert(U.colsize() == A.size());
        TMVAssert(U.rowsize() == A.size());
        TMVAssert(SS.size() == A.size());
        TMVAssert(V.colsize() == A.size());
        TMVAssert(V.rowsize() == A.size());

        if (A.isconj()) {
            if (U.isconj()) {
                if (V.isconj()) {
                    SV_Decompose(A.conjugate(),U.conjugate(),SS,V.conjugate());
                } else {
                    SV_Decompose(A.conjugate(),U.conjugate(),SS,V);
                    V.conjugateSelf();
                }
            } else {
                if (V.isconj()) {
                    SV_Decompose(A.conjugate(),U,SS,V.conjugate());
                    U.conjugateSelf();
                } else {
                    SV_Decompose(A.conjugate(),U,SS,V);
                    U.conjugateSelf();
                    V.conjugateSelf();
                }
            }
        } else {
            if (U.isconj()) {
                if (V.isconj()) {
                    SV_Decompose(A,U.conjugate(),SS,V.conjugate());
                    U.conjugateSelf();
                    V.conjugateSelf();
                } else {
                    SV_Decompose(A,U.conjugate(),SS,V);
                    U.conjugateSelf();
                }
            } else {
                if (V.isconj()) {
                    SV_Decompose(A,U,SS,V.conjugate());
                    V.conjugateSelf();
                } else {
                    RT ld(0);
                    T d(0);
                    SV_Decompose<T>(A,U,SS,V,ld,d);
                    if (A.isherm()) {
                        // Then S values might be negative:
                        for(size_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) {
                            SS(i) = -SS(i);
                            V.row(i) = -V.row(i);
                        }
                    }
                }
            }
        }
    }

    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        const MatrixView<T>& U, const DiagMatrixView<RT>& SS)
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
                SV_Decompose<T>(A,U,SS,0,ld,d);
                if (A.isherm()) 
                    for(size_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) 
                        SS(i) = -SS(i);
            }
        }
    }

    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        const DiagMatrixView<RT>& SS, const MatrixView<T>& V)
    {
        TMVAssert(SS.size() == A.size());
        TMVAssert(V.colsize() == A.size());
        TMVAssert(V.rowsize() == A.size());

        if (A.isherm()) SV_Decompose(A,V.adjoint(),SS);
        else SV_Decompose(A,V.transpose(),SS);
    }

    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A, const DiagMatrixView<RT>& SS)
    {
        TMVAssert(SS.size() == A.size());

        if (A.isconj()) {
            SV_Decompose(A.conjugate(),SS);
        } else {
            RT ld(0);
            T d(0);
            SV_Decompose<T>(A,0,SS,0,ld,d);
            if (A.isherm()) 
                for(size_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) 
                    SS(i) = -SS(i);
        }
    }

    template <class T> 
    void SquareRoot(const GenSymBandMatrix<T>& A, const SymMatrixView<T>& SS)
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
                for(size_t i=0;i<A.size();i++) {
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

#define InstFile "TMV_SymBandSVDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


