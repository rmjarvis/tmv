///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_Blas.h"
#include "TMV_Matrix.h"
#include "TMV_SymBandMatrix.h"
#include "TMV_Givens.h"
#include "TMV_DiagMatrix.h"
#include "TMV_BandMatrix.h"
#include "TMV_Householder.h"
#include "TMV_SymHouseholder.h"
#include "TMV_SymBandSVDiv.h"
#include "TMV_SymSVDiv.h"
#include "TMV_BandSVDiv.h"
#include "TMV_QRDiv.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_SymMatrixArith.h"
#include "TMV_SymBandMatrixArith.h"
#include "TMV_BandMatrixArith.h"
#include <sstream>

//#define XDEBUG
//#define TIME

#ifdef XDEBUG
#include "TMV_DiagMatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

#ifdef TIME
#include <sys/time.h>
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_TRIDIAG_BLOCKSIZE TMV_BLOCKSIZE
#else
#define SYM_TRIDIAG_BLOCKSIZE 64
#endif

  template <class T> inline const MatrixView<T>* ZMV() 
  { return (const MatrixView<T>*)(0); }

  template <class T> inline void MakeTridiagReal(
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

  template <class T> inline void MakeTridiagReal(
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

    const size_t N = D.size();

    std::complex<T>* Uj = Udiag.ptr();
    *Uj = T(1);
    T* Ej = E.ptr();
    const std::complex<T>* cEj = cE.cptr();
    const int cEstep = cE.step();

    for(size_t j=1;j<N;++j,++Ej,cEj+=cEstep) {
      std::complex<T> xcEj = (*Uj) * (*cEj);
      ++Uj;
      *Ej = ABS(xcEj);
      *Uj = SIGN(xcEj,*Ej);
    }
    D = cD.Real();
  }

  template <class T> inline void MakeTridiagReal(
      const VectorView<std::complex<T> >& , 
      const GenVector<std::complex<T> >& ,
      const GenVector<std::complex<T> >& ,
      const VectorView<std::complex<T> >& , const VectorView<T>& ) 
  {
    // This one is the symmetric case, which should never get called.
    TMVAssert(FALSE); 
  }

  template <class T, class Td> inline void NonLapTridiagonalize(
      const GenSymBandMatrix<T>& A, const MatrixView<T>* U,
      const VectorView<Td>& D, const VectorView<RealType(T)>& E, T& det)
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
      TMVAssert(U->iscm());
      TMVAssert(U->ct() == NonConj);
    }
    TMVAssert(D.size() == A.size());
    TMVAssert(E.size() == A.size()-1);
    TMVAssert(IsReal(Td()) || !A.isherm());

    const size_t N = A.size();
    const int nlo = A.nlo();

    if (nlo == 1) {
      TMVAssert(A.isherm());
      // MJ: This is slightly inefficient, making a new vector for the
      // Ud and then copying it into U.diag() if necessary (which is most
      // of the time I would imagine).  It makes the coding of 
      // MakeTridiagReal a bit easier, since it doesn't require checking
      // the step or conj of Ud, and also whether to write the Uj value
      // into the Ud parameter.  Anyway, it seems to be a small inefficiency,
      // but I should probably still fix it at some point.
      // Also, the corresponding routine in TMV_BandSVDiv_A.cpp has the same
      // inefficiency.
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

      std::vector<size_t> vec(N-1);
      Vector<T> Ubeta(N-1);
      size_t endcol = nlo+1;
      if (endcol > N) endcol = N;

      T* Ubj = Ubeta.ptr();
      // We use Householder reflections to reduce A to the tridiagonal form:
      for(size_t j=0;j<N-1;++j,++Ubj) {
	*Ubj = Householder_Reflect(U1->col(j,j+1,endcol),det);
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
	  TMVAssert(NormInf(U1->diag().Imag()) == RealType(T)(0));
	TMVAssert(NormInf(U1->diag(-1).Imag()) == RealType(T)(0));
      }
#endif

      if (U) {
	UpperTriMatrixViewOf(*U).Zero();
	U->diag(-1).Zero();
	for (size_t j=N-1;j>0;--j) {
	  endcol = vec[j-1];
	  U->col(j,j+1,endcol) = U->col(j-1,j+1,endcol);
	  Householder_Unpack(U->SubMatrix(j,endcol,j,N),*(--Ubj));
	}
	U->col(0,0,vec[0]).MakeBasis(0);
      }

      if (!A.isherm()) det *= det;
    }
  }

#ifdef LAP
  template <class T, class Td> inline void LapTridiagonalize(
      const GenSymBandMatrix<T>& A, const MatrixView<T>* U,
      const VectorView<Td>& D, const VectorView<RealType(T)>& E, T& det)
  { NonLapTridiagonalize(A,U,D,E,det); }
#ifdef INST_DOUBLE
  template <> inline void LapTridiagonalize(
      const GenSymBandMatrix<double>& A, const MatrixView<double>* U,
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
  template <> inline void LapTridiagonalize(
      const GenSymBandMatrix<std::complex<double> >& A,
      const MatrixView<std::complex<double> >* U,
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
  template <> inline void LapTridiagonalize(
      const GenSymBandMatrix<float>& A, const MatrixView<float>* U,
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
  template <> inline void LapTridiagonalize(
      const GenSymBandMatrix<std::complex<float> >& A,
      const MatrixView<std::complex<float> >* U,
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

  template <class T, class Td> inline void Tridiagonalize(
      const GenSymBandMatrix<T>& A, const MatrixView<T>* U,
      const VectorView<Td>& D, const VectorView<RealType(T)>& E, T& det)
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
#endif

    if (A.size() > 0) {
      TMVAssert(E.size()+1 == D.size());
      if (A.uplo() == Upper) {
	if (A.isherm()) Tridiagonalize(A.Adjoint(),U,D,E,det);
	else Tridiagonalize(A.Transpose(),U,D,E,det);
      } else {
#ifdef LAP
	if (A.isherm() && A.nlo() > 1 && D.step() == 1 && E.step() == 1) {
	  TMVAssert(IsReal(Td()));
	  LapTridiagonalize(A,U,D,E,det);
	}
	else 
#endif // LAP
	  NonLapTridiagonalize(A,U,D,E,det);
      }
    }

#ifdef XDEBUG
    if (U) {
      const size_t N = A.size();
      Matrix<T> TT(N,N,T(0));
      TT.diag() = D;
      TT.diag(1) = TT.diag(-1) = E;
      Matrix<T> A2 = (*U)*TT*(A.isherm() ? U->Adjoint() : U->Transpose());
      if (Norm(A2-A0) > 0.001*Norm(A0)) {
	cerr<<"SymBandTridiagonalize: \n";
	cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
	cerr<<"Done: U = "<<*U<<endl;
	cerr<<"D = "<<D<<endl;
	cerr<<"E = "<<E<<endl;
	cerr<<"TT = "<<TT<<endl;
	cerr<<"UU * TT * UUt = ";
	A2.Write(cerr,1.e-12);
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

  template <class T> void HermBandSV_Decompose(
      const GenSymBandMatrix<T>& A,
      const MatrixView<T>* U, const VectorView<RealType(T)>& S, T& det)
  {
    TMVAssert(A.size() > 0);
    if (U) {
      TMVAssert(U->rowsize() == A.size());
      TMVAssert(U->colsize() == A.size());
      TMVAssert(U->ct() == NonConj);
    }
    TMVAssert(S.size() == A.size());
    TMVAssert(S.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif
    // Decompose Hermitian A (input as lower tri of U) into U S Ut
    // where S is a diagonal real matrix, and U is a unitary matrix.
    // U,S are N x N
    const size_t N = A.size();
    if (N == 0) return;

    // First we reduce A to tridiagonal form: A = U * T * Ut
    // using a series of Householder transformations.
    // The diagonal of the Tridiagonal Matrix T is stored in D.
    // The subdiagonal is stored in E.
    if (A.nlo() == 0) {
      S = A.diag().Real();
      det *= DiagMatrixViewOf(S).Det();
      auto_array<size_t> sortp(new size_t[N]);
      S.Sort(sortp.get(),DESCEND,ABS_COMP);
      if (U) {
	U->SetToIdentity();
	U->PermuteCols(sortp.get());
      }
    } else {
      Vector<RealType(T)> E(N-1);
      det = 0;
      Tridiagonalize(A,U,S,E.View(),det);

      HermSV_Decompose_From_Tridiagonal(U,S,E.View());
      det = DiagMatrixViewOf(S).Det();
    }

#ifdef XDEBUG
    if (U) {
      Matrix<T> A2 = (*U) * DiagMatrixViewOf(S) * U->Adjoint();
      if (Norm(A0-A2) > 0.0001 * Norm(*U) * Norm(S) * Norm(*U)) {
	cerr<<"HermSV_Decompose:\n";
	cerr<<"A = "<<A0<<endl;
	cerr<<"U = "<<*U<<endl;
	cerr<<"S = "<<S<<endl;
	cerr<<"USUt = "<<A2<<endl;
	abort();
      }
    }
#endif
  }

  template <class T> void SymBandSV_Decompose(
      const GenSymBandMatrix<T>& A,
      const MatrixView<T>* U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>* V, T& det)
  {
    TMVAssert(IsComplex(T()));
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
    TMVAssert(S.size() == A.size());
    TMVAssert(S.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif
    // Decompose complex symmetric A (input as lower tri of U) into U S V
    // where S is a diagonal real matrix, and U,V are unitary matrices.
    // U,S,V are N x N
    // If V = 0, then U,V are not formed.  Only S,det are accurate on return.
    // The determinant is returned in det.
    // (Technically, det is multiplied by the determinant, so det should
    // be set to 1 on entry.)
    const size_t N = A.size();
    if (N == 0) return;

    // First we reduce A to tridiagonal form: A = U * T * UT
    // using a series of Householder transformations.
    // The diagonal of the Tridiagonal Matrix T is stored in D.
    // The subdiagonal is stored in E.
    if (A.nlo() == 0) {
      BandSV_Decompose(A.Diags(0,1),U,S,V,det);
    } else {

      if (A.nlo() == 1) {
	BandMatrix<T,ColMajor> B = A;
	BandSV_Decompose(B,U,S,V,det);
      } else {
	Vector<T> D(N);
	Vector<RealType(T)> E(N-1);
	if (V && !U) {
	  MatrixView<T> Vt = V->Transpose();
	  Tridiagonalize(A,&Vt,D.View(),E.View(),det);
	} else {
	  Tridiagonalize(A,U,D.View(),E.View(),det);
	  if (V) {
	    TMVAssert(U);
	    *V = U->Transpose();
	  }
	}

	BandMatrix<T,ColMajor> B(N,N,1,1);
	B.diag() = D;
	B.diag(-1) = E;
	B.diag(1) = E;

	if (U) {
	  Matrix<T,ColMajor> U1(N,N);
	  MatrixView<T> U1v = U1.View();
	  if (V) {
	    Matrix<T,ColMajor> V1(N,N);
	    MatrixView<T> V1v = V1.View();
	    BandSV_Decompose(B,&U1v,S,&V1v,det);
	    *V = V1*(*V);
	  } else {
	    BandSV_Decompose(B,&U1v,S,ZMV<T>(),det);
	  }
	  *U = *U*U1;
	} else {
	  if (V) {
	    Matrix<T,ColMajor> V1(N,N);
	    MatrixView<T> V1v = V1.View();
	    BandSV_Decompose(B,ZMV<T>(),S,&V1v,det);
	    *V = V1*(*V);
	  } else {
	    BandSV_Decompose(B,ZMV<T>(),S,ZMV<T>(),det);
	  }
	}
      }
    }
#ifdef XDEBUG
    if (U&&V) {
      Matrix<T> A2 = (*U) * DiagMatrixViewOf(S) * (*V);
      if (Norm(A0-A2) > 0.0001 * Norm(*U) * Norm(S) * Norm(*V)) {
	cerr<<"SymSV_Decompose:\n";
	cerr<<"A = "<<A0<<endl;
	cerr<<"U = "<<U<<endl;
	cerr<<"S = "<<S<<endl;
	cerr<<"V = "<<*V<<endl;
	cerr<<"USV = "<<A2<<endl;
	abort();
      }
    }
#endif
  }

#define InstFile "TMV_SymBandSVDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


