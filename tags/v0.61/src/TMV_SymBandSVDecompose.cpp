///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
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



#include "TMV_Blas.h"
#include "TMV_SymBandSVDiv.h"
#include "TMV_SymBandSVD.h"
#include "TMV_SymBandMatrix.h"
#include "TMV_Matrix.h"
#include "TMV_Vector.h"
#include "TMV_DiagMatrix.h"
#include "TMV_DiagMatrixArith.h"
#include "TMV_SymSVDiv.h"
#include "TMV_BandSVDiv.h"
#include "TMV_MatrixArith.h"
#include "TMV_VectorArith.h"
#include "TMV_Householder.h"
#include "TMV_SymHouseholder.h"
#include <sstream>

//#define XDEBUG
//#define TIME

#ifdef XDEBUG
#include "TMV_DiagMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

#ifdef TIME
#include <sys/time.h>
#endif

namespace tmv {

#define RT RealType(T)

#ifdef TMV_BLOCKSIZE
#define SYM_TRIDIAG_BLOCKSIZE TMV_BLOCKSIZE
#else
#define SYM_TRIDIAG_BLOCKSIZE 64
#endif

  template <class T> static void MakeTridiagReal(
      const VectorView<T>& Udiag, 
      const GenVector<T>& cD, const GenVector<T>& cE,
      const VectorView<T>& D, const VectorView<T>& E)
  {
    TMVAssert(cD.size() == D.size());
    TMVAssert(cE.size() == E.size());
    TMVAssert(D.size() == Udiag.size());
    TMVAssert(E.size() == D.size()-1);

    Udiag.SetAllTo(T(1));
    D = cD;
    E = cE;
  } 

  template <class T> static void MakeTridiagReal(
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
      *Ej = ABS(xcEj);
      *Uj = SIGN(xcEj,*Ej);
    }
    D = cD.Real();
  }

  template <class T> static inline void MakeTridiagReal(
      const VectorView<std::complex<T> >& , 
      const GenVector<std::complex<T> >& ,
      const GenVector<std::complex<T> >& ,
      const VectorView<std::complex<T> >& , const VectorView<T>& ) 
  {
    // This one is the symmetric case, which should never get called.
    TMVAssert(FALSE); 
  }

  template <class T, class Td> static void NonLapTridiagonalize(
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
    TMVAssert(IsReal(Td()) || !A.isherm());

    const int N = A.size();
    const int nlo = A.nlo();

    if (nlo == 1) {
      TMVAssert(A.isherm());
      Vector<T> Ud(N);
      if (A.isconj()) {
	MakeTridiagReal(Ud.View(),A.diag(),A.diag(-1).Conjugate(),D,E);
	Ud.ConjugateSelf();
      }
      else
	MakeTridiagReal(Ud.View(),A.diag(),A.diag(-1),D,E);
      if (U) {
	U->Zero();
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
	  U1.reset(new SymMatrixView<T>(UUS->View()));
	} else {
	  UUH.reset(new HermMatrix<T,Lower,ColMajor>(A));
	  U1.reset(new SymMatrixView<T>(UUH->View()));
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
	*Ubj = Householder_Reflect(U1->col(j,j+1,endcol),signdet);
	if (endcol < N) {
	  endcol+=nlo;
	  if (endcol > N) endcol = N;
	}
	vec[j] = endcol;
	if (*Ubj != T(0)) 
	  Householder_LRMult(U1->col(j,j+2,endcol),*Ubj,
	      U1->SubSymMatrix(j+1,endcol));
      }

      // The tridiagonal of U1 is the tridiagonal we want, so copy it to D,E
      if (IsReal(Td())) D = U1->diag().Real();
      else D = U1->diag();
      E = U1->diag(-1).Real();
#ifdef XTEST
      if (IsComplex(T())) {
	if (IsReal(Td()))
	  TMVAssert(NormInf(U1->diag().Imag()) == RT(0));
	TMVAssert(NormInf(U1->diag(-1).Imag()) == RT(0));
      }
#endif

      if (U) {
	UpperTriMatrixViewOf(*U).Zero();
	U->diag(-1).Zero();
	for (int j=N-1;j>0;--j) {
	  endcol = vec[j-1];
	  U->col(j,j+1,endcol) = U->col(j-1,j+1,endcol);
	  Householder_Unpack(U->SubMatrix(j,endcol,j,N),*(--Ubj));
	}
	U->col(0,0,vec[0]).MakeBasis(0);
      }

      if (!A.isherm()) signdet *= signdet;
    }
  }

#ifdef LAP
  template <class T, class Td> static inline void LapTridiagonalize(
      const GenSymBandMatrix<T>& A, MVP<T> U,
      const VectorView<Td>& D, const VectorView<RT>& E, T& signdet)
  { NonLapTridiagonalize(A,U,D,E,signdet); }
#ifdef INST_DOUBLE
  template <> void LapTridiagonalize(
      const GenSymBandMatrix<double>& A, MVP<double> U,
      const VectorView<double>& D, const VectorView<double>& E, double& )
  {
    TMVAssert(A.nlo() > 1);
    TMVAssert(A.size() > A.nlo());
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
#ifndef LAPNOWORK
    int lwork = n;
    double* work = LAP_DWork(lwork);
#endif
    SymBandMatrix<double,Lower,ColMajor> A2 = A;
    int lda = A2.diagstep();
    int ldu = U ? U->stepj() : 1;
    char vect = U ? 'V' : 'N';
    double* UU = U ? U->ptr() : 0;

    LAPNAME(dsbtrd) (LAPCM LAPV(vect),LAPCH_LO,LAPV(n),LAPV(kl),
	LAPP(A2.cptr()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
	LAPP(UU),LAPV(ldu) LAPWK(work) LAPINFO LAP1 LAP1 );
    LAP_Results("dsbtrd");
  }
  template <> void LapTridiagonalize(
      const GenSymBandMatrix<std::complex<double> >& A,
      MVP<std::complex<double> > U,
      const VectorView<double>& D, const VectorView<double>& E, 
      std::complex<double>& )
  {
    TMVAssert(A.isherm());
    TMVAssert(A.nlo() > 1);
    TMVAssert(A.size() > A.nlo());
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
#ifndef LAPNOWORK
    int lwork = n;
    std::complex<double>* work = LAP_ZWork(lwork);
#endif
    HermBandMatrix<std::complex<double>,Lower,ColMajor> A2 = A;
    int lda = A2.diagstep();
    int ldu = U ? U->stepj() : 1;
    char vect = U ? 'V' : 'N';
    std::complex<double>* UU = U ? U->ptr() : 0;

    LAPNAME(zhbtrd) (LAPCM LAPV(vect),LAPCH_LO,LAPV(n),LAPV(kl),
	LAPP(A2.cptr()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
	LAPP(UU),LAPV(ldu) LAPWK(work) LAPINFO LAP1 LAP1 );
    LAP_Results("zhbtrd");
  }
#endif
#ifdef INST_FLOAT
  template <> void LapTridiagonalize(
      const GenSymBandMatrix<float>& A, MVP<float> U,
      const VectorView<float>& D, const VectorView<float>& E, float& )
  {
    TMVAssert(A.nlo() > 1);
    TMVAssert(A.size() > A.nlo());
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
#ifndef LAPNOWORK
    int lwork = n;
    float* work = LAP_SWork(lwork);
#endif
    SymBandMatrix<float,Lower,ColMajor> A2 = A;
    int lda = A2.diagstep();
    int ldu = U ? U->stepj() : 1;
    char vect = U ? 'V' : 'N';
    float* UU = U ? U->ptr() : 0;

    LAPNAME(ssbtrd) (LAPCM LAPV(vect),LAPCH_LO,LAPV(n),LAPV(kl),
	LAPP(A2.cptr()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
	LAPP(UU),LAPV(ldu) LAPWK(work) LAPINFO LAP1 LAP1 );
    LAP_Results("dsbtrd");
  }
  template <> void LapTridiagonalize(
      const GenSymBandMatrix<std::complex<float> >& A,
      MVP<std::complex<float> > U,
      const VectorView<float>& D, const VectorView<float>& E, 
      std::complex<float>& )
  {
    TMVAssert(A.isherm());
    TMVAssert(A.nlo() > 1);
    TMVAssert(A.size() > A.nlo());
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
#ifndef LAPNOWORK
    int lwork = n;
    std::complex<float>* work = LAP_CWork(lwork);
#endif
    HermBandMatrix<std::complex<float>,Lower,ColMajor> A2 = A;
    int lda = A2.diagstep();
    int ldu = U ? U->stepj() : 1;
    char vect = U ? 'V' : 'N';
    std::complex<float>* UU = U ? U->ptr() : 0;

    LAPNAME(chbtrd) (LAPCM LAPV(vect),LAPCH_LO,LAPV(n),LAPV(kl),
	LAPP(A2.cptr()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
	LAPP(UU),LAPV(ldu) LAPWK(work) LAPINFO LAP1 LAP1 );
    LAP_Results("chbtrd");
  }
#endif 
#endif // LAP

  template <class T, class Td> static void Tridiagonalize(
      const GenSymBandMatrix<T>& A, MVP<T> U,
      const VectorView<Td>& D, const VectorView<RT>& E, T& signdet)
  {
    TMVAssert(A.nlo() > 0);
    TMVAssert(A.size() == D.size());
    TMVAssert(E.size() == D.size()-1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(IsReal(Td()) || !A.isherm());
    if (U) {
      TMVAssert(U->colsize() == D.size());
      TMVAssert(U->rowsize() == D.size());
    }

#ifdef XDEBUG
    Matrix<T> A0 = A;
    cout<<"SymBandTridiagonalize: \n";
    cout<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
#endif

    if (A.size() > 0) {
      TMVAssert(E.size()+1 == D.size());
      if (A.uplo() == Upper) {
	if (A.isherm()) Tridiagonalize(A.Adjoint(),U,D,E,signdet);
	else Tridiagonalize(A.Transpose(),U,D,E,signdet);
      } else {
#ifdef LAP
	if (A.isherm() && A.nlo() > 1 && (!U || U->iscm()) &&
	    D.step() == 1 && E.step() == 1) {
	  TMVAssert(IsReal(Td()));
	  LapTridiagonalize(A,U,D,E,signdet);
	}
	else 
#endif // LAP
	  NonLapTridiagonalize(A,U,D,E,signdet);
      }
    }

#ifdef XDEBUG
    cout<<"Done Tridiag\n";
    cout<<"A => "<<A<<endl;
    cout<<"D => "<<D<<endl;
    cout<<"E => "<<E<<endl;
    if (U) {
      cout<<"U = "<<*U<<endl;
      const int N = A.size();
      Matrix<T> TT(N,N,T(0));
      TT.diag() = D;
      TT.diag(1) = TT.diag(-1) = E;
      cout<<"TT = "<<TT<<endl;
      Matrix<T> A2 = (*U)*TT*(A.isherm() ? U->Adjoint() : U->Transpose());
      cout<<"A2 = "<<A2<<endl;
      if (Norm(A2-A0) > 0.001*Norm(A0)) {
	cerr<<"SymBandTridiagonalize: \n";
	cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
	cerr<<"Done: U = "<<*U<<endl;
	cerr<<"D = "<<D<<endl;
	cerr<<"E = "<<E<<endl;
	cerr<<"TT = "<<TT<<endl;
	cerr<<"UU * TT * UUt = ";
	A2.Write(cerr,RT(1.e-12));
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

  template <class T> void UnsortedEigen(
      const GenSymBandMatrix<T>& A,
      MVP<T> U, const VectorView<RT>& SS)
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
    cout<<"Unsorted Eigen:\n";
    cout<<"A = "<<A0<<endl;
#endif
    // Decompose Hermitian A (input as lower tri of U) into U S Ut
    // where S is a diagonal real matrix, and U is a unitary matrix.
    // U,S are N x N
    const int N = A.size();
    if (N == 0) return;

    if (A.nlo() == 0) {
      SS = A.diag().Real();
      if (U) U->SetToIdentity();
    } else {
      // First we reduce A to tridiagonal form: A = U * T * Ut
      // using a series of Householder transformations.
      // The diagonal of the Tridiagonal Matrix T is stored in D.
      // The subdiagonal is stored in E.
      T d2 = 0;
      Vector<RT> E(N-1);
      Tridiagonalize(A,U,SS,E.View(),d2);

      // Then finish the decomposition as for a HermMatrix
      Eigen_From_Tridiagonal(U,SS,E.View());
    }

#ifdef XDEBUG
    cout<<"Done Unsorted Eigen\n";
    if (U) {
      Matrix<T> A2 = (*U) * DiagMatrixViewOf(SS) * U->Adjoint();
      if (Norm(A0-A2) > 0.0001 * Norm(*U) * Norm(SS) * Norm(*U)) {
	cerr<<"Unsorted Eigen:\n";
	cerr<<"A = "<<A0<<endl;
	cerr<<"U = "<<*U<<endl;
	cerr<<"S = "<<SS<<endl;
	cerr<<"USUt = "<<A2<<endl;
	abort();
      }
    }
#endif
  }

  template <class T> void SV_Decompose(
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
    cout<<"SymBandSV_Decompose:\n";
    cout<<"A = "<<A0<<endl;
#endif
    if (A.isherm()) {
      if (V && !U) {
	UnsortedEigen<T>(A,V->Transpose(),SS.diag());
	V->ConjugateSelf();
      } else  {
	UnsortedEigen<T>(A,U,SS.diag());
      }
      if (V && U) *V = U->Adjoint();
      if (signdet != T(0)) {
	RT s;
	logdet += SS.LogDet(&s);
	signdet *= s;
      }
      if (U || V) {
	auto_array<int> sortp(new int[A.size()]);
	SS.diag().Sort(sortp.get(),DESCEND,ABS_COMP);
	if (U) U->PermuteCols(sortp.get());
	if (V) V->PermuteRows(sortp.get());
      } else {
	SS.diag().Sort(DESCEND,ABS_COMP);
      }
    } else {
      TMVAssert(IsComplex(T()));
      // Decompose complex symmetric A (input as lower tri of U) into U S V
      // where S is a diagonal real matrix, and U,V are unitary matrices.
      // U,S,V are N x N
      // If V = 0, then U,V are not formed.  Only S,det are accurate on return.
      const int N = A.size();
      if (N == 0) return;

      if (A.nlo() == 0) {
	SV_Decompose(A.Diags(0,1),U,SS,V,logdet,signdet);
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
	  Tridiagonalize<T>(A,V->Transpose(),D.View(),E.View(),signdet);
	} else {
	  Tridiagonalize<T>(A,U,D.View(),E.View(),signdet);
	  if (V) { TMVAssert(U); *V = U->Transpose(); }
	}

	BandMatrix<T,ColMajor> B(N,N,1,1);
	B.diag() = D;
	B.diag(-1) = E;
	B.diag(1) = E;

	if (U) {
	  Matrix<T,ColMajor> U1(N,N);
	  if (V) {
	    Matrix<T,ColMajor> V1(N,N);
	    SV_Decompose<T>(B,U1.View(),SS,V1.View(),logdet,signdet);
	    *V = V1*(*V);
	  } else {
	    SV_Decompose<T>(B,U1.View(),SS,0,logdet,signdet);
	  }
	  *U = *U*U1;
	} else {
	  if (V) {
	    Matrix<T,ColMajor> V1(N,N);
	    SV_Decompose<T>(B,0,SS,V1.View(),logdet,signdet);
	    *V = V1*(*V);
	  } else {
	    SV_Decompose<T>(B,0,SS,0,logdet,signdet);
	  }
	}
      }
    }
#ifdef XDEBUG
    std::cout<<"Done SymBandSV_Decompose\n";
    if (U&&V) {
      Matrix<T> A2 = (*U) * SS * (*V);
      if (Norm(A0-A2) > 0.0001 * Norm(*U) * Norm(SS) * Norm(*V)) {
	cerr<<"SymBandSV_Decompose:\n";
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

  template <class T> void Eigen(
      const GenSymBandMatrix<T>& A, const MatrixView<T>& U,
      const VectorView<RT>& SS)
  {
    TMVAssert(SS.size() == A.size());
    TMVAssert(U.colsize() == A.size());
    TMVAssert(U.rowsize() == A.size());

    if (A.isconj()) {
      if (U.isconj()) {
	Eigen(A.Conjugate(),U.Conjugate(),SS);
      } else {
	Eigen(A.Conjugate(),U,SS);
	U.ConjugateSelf();
      }
    } else {
      if (U.isconj()) {
	Eigen(A,U.Conjugate(),SS);
	U.ConjugateSelf();
      } else {
	UnsortedEigen<T>(A,U,SS);
	auto_array<int> sortp(new int[A.size()]);
	SS.Sort(sortp.get(),ASCEND);
	U.PermuteCols(sortp.get());
      }
    }
  }

  template <class T> void Eigen(
      const GenSymBandMatrix<T>& A, const VectorView<RT>& SS)
  {
    TMVAssert(SS.size() == A.size());
    if (A.isconj()) {
      Eigen(A.Conjugate(),SS);
    } else {
      UnsortedEigen<T>(A,0,SS);
      SS.Sort(ASCEND);
    }
  }

  // Decompose A into U S V
  template <class T> void SV_Decompose(
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
	  SV_Decompose(A.Conjugate(),U.Conjugate(),SS,V.Conjugate());
	} else {
	  SV_Decompose(A.Conjugate(),U.Conjugate(),SS,V);
	  V.ConjugateSelf();
	}
      } else {
	if (V.isconj()) {
	  SV_Decompose(A.Conjugate(),U,SS,V.Conjugate());
	  U.ConjugateSelf();
	} else {
	  SV_Decompose(A.Conjugate(),U,SS,V);
	  U.ConjugateSelf();
	  V.ConjugateSelf();
	}
      }
    } else {
      if (U.isconj()) {
	if (V.isconj()) {
	  SV_Decompose(A,U.Conjugate(),SS,V.Conjugate());
	  U.ConjugateSelf();
	  V.ConjugateSelf();
	} else {
	  SV_Decompose(A,U.Conjugate(),SS,V);
	  U.ConjugateSelf();
	}
      } else {
	if (V.isconj()) {
	  SV_Decompose(A,U,SS,V.Conjugate());
	  V.ConjugateSelf();
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

  template <class T> void SV_Decompose(
      const GenSymBandMatrix<T>& A,
      const MatrixView<T>& U, const DiagMatrixView<RT>& SS)
  {
    TMVAssert(U.colsize() == A.size());
    TMVAssert(U.rowsize() == A.size());
    TMVAssert(SS.size() == A.size());

    if (A.isconj()) {
      if (U.isconj()) {
	SV_Decompose(A.Conjugate(),U.Conjugate(),SS);
      } else {
	SV_Decompose(A.Conjugate(),U,SS);
	U.ConjugateSelf();
      }
    } else {
      if (U.isconj()) {
	SV_Decompose(A,U.Conjugate(),SS);
	U.ConjugateSelf();
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

  template <class T> void SV_Decompose(
      const GenSymBandMatrix<T>& A,
      const DiagMatrixView<RT>& SS, const MatrixView<T>& V)
  {
    TMVAssert(SS.size() == A.size());
    TMVAssert(V.colsize() == A.size());
    TMVAssert(V.rowsize() == A.size());

    if (A.isconj()) {
      if (V.isconj()) {
	SV_Decompose(A.Conjugate(),SS,V.Conjugate());
      } else {
	SV_Decompose(A.Conjugate(),SS,V);
	V.ConjugateSelf();
      }
    } else {
      if (V.isconj()) {
	SV_Decompose(A,SS,V.Conjugate());
	V.ConjugateSelf();
      } else {
	RT ld(0);
	T d(0);
	SV_Decompose<T>(A,0,SS,V,ld,d);
	if (A.isherm()) 
	  for(size_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) 
	    SS(i) = -SS(i);
      }
    }
  }

  template <class T> void SV_Decompose(
      const GenSymBandMatrix<T>& A, const DiagMatrixView<RT>& SS)
  {
    TMVAssert(SS.size() == A.size());

    if (A.isconj()) {
      SV_Decompose(A.Conjugate(),SS);
    } else {
      RT ld(0);
      T d(0);
      SV_Decompose<T>(A,0,SS,0,ld,d);
      if (A.isherm()) 
	for(size_t i=0;i<SS.size();i++) if (SS(i) < RT(0)) 
	  SS(i) = -SS(i);
    }
  }

  template <class T> void SquareRoot(const GenSymBandMatrix<T>& A, 
      const SymMatrixView<T>& SS)
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
	SquareRoot(A.Conjugate(),SS.Conjugate());
      } else {
	SquareRoot(A.Conjugate(),SS);
	SS.ConjugateSelf();
      }
    } else {
      if (SS.isconj()) {
	SquareRoot(A,SS.Conjugate());
	SS.ConjugateSelf();
      } else {
	Matrix<T> V(A.size(),A.size());
	DiagMatrix<RT> D(A.size());
	Eigen(A,V.View(),D.diag());
	for(size_t i=0;i<A.size();i++) {
	  if (D(i) < RT(0)) throw NonPosDef("in SymBandMatrix SquareRoot");
	  D(i) = SQRT(D(i));
	}
	Matrix<T> DVt = D*V.Adjoint();
	SymMultMM<false>(T(1),V,DVt,SS);
      }
    }
  }

#undef RT

#define InstFile "TMV_SymBandSVDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


