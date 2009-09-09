
#include "TMV_Band.h"
#include "TMV_SVDiv.h"
#include "TMV_Householder.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

  template <class T> inline void MakeBidiagReal(
      const VectorView<T>& Udiag, const VectorView<T>& Vdiag, 
      const GenVector<T>& cD, const GenVector<T>& cE,
      const VectorView<RealType(T)>& D, const VectorView<RealType(T)>& E,
      T& det)
  {
    TMVAssert(Vdiag.size() == Udiag.size());
    TMVAssert(cD.size() == D.size());
    TMVAssert(cE.size() == E.size());
    TMVAssert(D.size() == Udiag.size());
    TMVAssert(E.size() == D.size()-1);

    const size_t N = D.size();

    if (IsReal(T())) {
      Udiag.SetAllTo(T(1));
      Vdiag.SetAllTo(T(1));
      D = cD;
      E = cE;
    } else {
      Vdiag(0) = T(1);
      T cDj = cD(0);
      for(size_t j=0;j<N-1;j++) {
	D(j) = abs(cDj);
	Udiag(j) = SIGN(cDj,D(j));
	T cEj = CONJ(Udiag(j))*cE(j);
	E(j) = abs(cEj);
	Vdiag(j+1) = SIGN(cEj,E(j));
	cDj = CONJ(Vdiag(j+1))*cD(j+1);
      }
      D(N-1) = abs(cDj);
      Udiag(N-1) = SIGN(cDj,D(N-1));
      det *= DiagMatrixViewOf(Udiag).Det() * DiagMatrixViewOf(Vdiag).Det();
    }
  }

  template <class T> inline void NonLapBidiagonalize(
      const GenBandMatrix<T>& A,
      const MatrixView<T>* U, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, const MatrixView<T>* V, T& det)
  {
    // Decompose A into U B V
    // The Bidiagonal Matrix B is stored as two vectors: D, E
    // D is the diagonal, E is the super-diagonal
    // We use Householder reflections to reduce A to the bidiagonal form:

    TMVAssert(A.rowsize() <= A.colsize());
    TMVAssert(A.rowsize() > 0);
    if (U) {
      TMVAssert(U->colsize() == A.colsize());
      TMVAssert(U->rowsize() == A.rowsize());
    } 
    if (V) {
      TMVAssert(V->colsize() == A.rowsize());
      TMVAssert(V->rowsize() == A.rowsize());
    }
    TMVAssert(D.size() == A.rowsize());
    TMVAssert(D.size() == E.size()+1);

    const size_t M = A.colsize();
    const size_t N = A.rowsize();

    const int nlo = A.nlo();
    const int nhi = A.nhi();

    if (nlo == 0 && nhi == 1) {
      if (U) {
	U->Zero();
	if (V) {
	  V->Zero();
	  MakeBidiagReal(U->diag(),V->diag(),A.diag(),A.diag(1),D,E,det);
	} else {
	  Vector<T> Vd(N);
	  MakeBidiagReal(U->diag(),Vd.View(),A.diag(),A.diag(1),D,E,det);
	}
      } else {
	Vector<T> Ud(N);
	if (V) {
	  V->Zero();
	  MakeBidiagReal(Ud.View(),V->diag(),A.diag(),A.diag(1),D,E,det);
	} else {
	  Vector<T> Vd(N);
	  MakeBidiagReal(Ud.View(),Vd.View(),A.diag(),A.diag(1),D,E,det);
	}
      }
    } else if (A.IsSquare() && nlo == 1 && nhi == 0) {
      if (U) {
	U->Zero();
	VectorView<T> Udiag = U->SubVector(N-1,0,-1,1,N);
	if (V) {
	  V->Zero();
	  VectorView<T> Vdiag = V->SubVector(0,N-1,1,-1,N);
	  MakeBidiagReal(Udiag,Vdiag,A.diag().Reverse(),
	      A.diag(-1).Reverse(),D,E,det);
	} else {
	  Vector<T> Vd(N);
	  MakeBidiagReal(Udiag,Vd.View(),A.diag().Reverse(),
	      A.diag(-1).Reverse(),D,E,det);
	}
      } else {
	Vector<T> Ud(N);
	if (V) {
	  V->Zero();
	  VectorView<T> Vdiag = V->SubVector(0,N-1,1,-1,N);
	  MakeBidiagReal(Ud.View(),Vdiag,A.diag().Reverse(),
	      A.diag(-1).Reverse(),D,E,det);
	} else {
	  Vector<T> Vd(N);
	  MakeBidiagReal(Ud.View(),Vd.View(),A.diag().Reverse(),
	      A.diag(-1).Reverse(),D,E,det);
	}
      }
    } else {
      auto_ptr<Matrix<T,ColMajor> > UU(0);
      auto_ptr<MatrixView<T> > U1(0);
      if (U) {
	*U = A;
	U1.reset(new MatrixView<T>(U->View()));
      } else {
	UU.reset(new Matrix<T,ColMajor>(A));
	U1.reset(new MatrixView<T>(UU->View()));
      }

      vector<size_t> vec(N), ver(N-1);
      size_t endcol = nlo+1;
      Vector<T> Ubeta(N);
      Vector<T> Vbeta(N-1);
      for(size_t j=0;j<N-1;++j) {
	vec[j] = endcol;
	size_t endrow = min(endcol+nhi,N);
	ver[j] = endrow;
	Ubeta(j) = Householder_Reflect(U1->SubMatrix(j,endcol,j,endrow),det);
	if (endcol < M) endcol = min(endrow+nlo,M);
	Vbeta(j) = Householder_Reflect(U1->Transpose().SubMatrix(
	      j+1,endrow,j,endcol),det);
      }
      vec[N-1] = endcol;
      Ubeta(N-1) = Householder_Reflect(U1->SubMatrix(N-1,endcol,N-1,N),det);

      // Now U stores Householder vectors for U in lower diagonal columns (HLi)
      // and Householder vectors for V in upper diagonal rows (HRi)
      // except for the bidiagonal which is the bidiagonal we want:
      if (IsComplex(T())) {
	TMVAssert(NormInf(U1->diag().Imag()) == RealType(T)(0));
	TMVAssert(NormInf(U1->diag(1).Imag()) == RealType(T)(0));
      }
      D = U1->diag().Real();
      E = U1->diag(1).Real();

      if (V) {
	V->SetToIdentity();
	for (int j=N-2;j>=0;--j) {
	  V->row(j+1,j+2,ver[j]) = U1->row(j,j+2,ver[j]);
	  Householder_Unpack(V->Transpose().SubMatrix(j+1,ver[j],j+1,N),
	      Vbeta(j));
	}
      }

      if (!UU.get()) {
	U->diag().Zero();
	U->diag(1).Zero();
	Householder_Unpack(U->SubMatrix(N-1,vec[N-1],N-1,N),Ubeta(N-1));
	for (int j=N-2;j>=0;--j) {
	  U->row(j,j,ver[j]).Zero();
	  Householder_Unpack(U->SubMatrix(j,vec[j],j,N),Ubeta(j));
	}
      }
    }
    det *= DiagMatrixViewOf(D).Det();
  }

#ifdef LAP
  template <class T> inline void LapBidiagonalize(
      const GenBandMatrix<T>& A,
      const MatrixView<T>* U, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, const MatrixView<T>* V,
      T& det)
  { NonLapBidiagonalize(A,U,D,E,V,det); }
  template <> inline void LapBidiagonalize(const GenBandMatrix<double>& A,
      const MatrixView<double>* U, const VectorView<double>& D,
      const VectorView<double>& E, const MatrixView<double>* V, 
      double& det)
  {
    TMVAssert(A.rowsize() == A.colsize());
    // The Lap routines can do NonSquare matrices, but they want to
    // write out to a square (MxM) U matrix which is larger than
    // what we have stored here.
    TMVAssert(A.rowsize() > 0);
    if (U) {
      TMVAssert(U->colsize() == A.colsize());
      TMVAssert(U->rowsize() == A.rowsize());
      TMVAssert(U->iscm());
      TMVAssert(U->ct() == NonConj);
    }
    if (V) {
      TMVAssert(V->colsize() == A.rowsize());
      TMVAssert(V->rowsize() == A.rowsize());
      TMVAssert(V->iscm());
      TMVAssert(V->ct() == NonConj);
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
    int lwork = 2*max(m,n);
    double* work = LAP_DWork(lwork);
#endif
    // LAP version overwrites original BandMatrix with crap.
    // Hence, copy BandMatrix before running.
    BandMatrix<double,ColMajor> A2 = A;
    int lda = A2.diagstep();
    int ldu = U ? U->stepj() : 1;
    int ldv = V ? V->stepj() : 1;
    char vect = U ? V ? 'B' : 'Q' : V ? 'P' : 'N';
    double* VV = V ? V->ptr() : 0;
    double* UU = U ? U->ptr() : 0;

    LAPNAME(dgbbrd) (LAPCM LAPV(vect),LAPV(m),LAPV(n),LAPV(ncc),
	LAPV(kl),LAPV(ku),
	LAPP(A2.cptr()-A.nhi()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
	LAPP(UU),LAPV(ldu),LAPP(VV),LAPV(ldv),
	0,LAPV(n) LAPWK(work) LAPINFO LAP1 );

    det *= DiagMatrixViewOf(D).Det();
    LAP_Results("dgbbrd");
  }
  template <> inline void LapBidiagonalize(
      const GenBandMatrix<complex<double> >& A,
      const MatrixView<complex<double> >* U, const VectorView<double>& D,
      const VectorView<double>& E, const MatrixView<complex<double> >* V, 
      complex<double>& det)
  {
    TMVAssert(A.rowsize() == A.colsize());
    TMVAssert(A.rowsize() > 0);
    if (U) {
      TMVAssert(U->colsize() == A.colsize());
      TMVAssert(U->rowsize() == A.rowsize());
      TMVAssert(U->iscm());
      TMVAssert(U->ct() == NonConj);
    }
    if (V) {
      TMVAssert(V->colsize() == A.rowsize());
      TMVAssert(V->rowsize() == A.rowsize());
      TMVAssert(V->iscm());
      TMVAssert(V->ct() == NonConj);
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
    int lwork = max(m,n);
    complex<double>* work = LAP_ZWork(lwork);
    double* rwork = LAP_DWork(lwork);
#endif
    BandMatrix<complex<double>,ColMajor> A2 = A;
    int lda = A2.diagstep();
    int ldu = U ? U->stepj() : 1;
    int ldv = V ? V->stepj() : 1;
    char vect = U ? V ? 'B' : 'Q' : V ? 'P' : 'N';
    complex<double>* VV = V ? V->ptr() : 0;
    complex<double>* UU = U ? U->ptr() : 0;

    LAPNAME(zgbbrd) (LAPCM LAPV(vect),LAPV(m),LAPV(n),LAPV(o),
	LAPV(kl),LAPV(ku),LAPP(A2.cptr()-A.nhi()),LAPV(lda),
	LAPP(D.ptr()),LAPP(E.ptr()),LAPP(UU),LAPV(ldu),
	LAPP(VV),LAPV(ldv),0,LAPV(n)
	LAPWK(work) LAPWK(rwork) LAPINFO LAP1);

    det *= DiagMatrixViewOf(A2.diag()).Det();
    LAP_Results("zgbbrd");
  }
#ifndef NOFLOAT
  template <> inline void LapBidiagonalize(const GenBandMatrix<float>& A,
      const MatrixView<float>* U, const VectorView<float>& D,
      const VectorView<float>& E, const MatrixView<float>* V, 
      float& det)
  {
    TMVAssert(A.rowsize() == A.colsize());
    // The Lap routines can do NonSquare matrices, but they want to
    // write out to a square (MxM) U matrix which is larger than
    // what we have stored here.
    TMVAssert(A.rowsize() > 0);
    if (U) {
      TMVAssert(U->colsize() == A.colsize());
      TMVAssert(U->rowsize() == A.rowsize());
      TMVAssert(U->iscm());
      TMVAssert(U->ct() == NonConj);
    }
    if (V) {
      TMVAssert(V->colsize() == A.rowsize());
      TMVAssert(V->rowsize() == A.rowsize());
      TMVAssert(V->iscm());
      TMVAssert(V->ct() == NonConj);
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
    int lwork = 2*max(m,n);
    float* work = LAP_SWork(lwork);
#endif
    // LAP version overwrites original BandMatrix with crap.
    // Hence, copy BandMatrix before running.
    BandMatrix<float,ColMajor> A2 = A;
    int lda = A2.diagstep();
    int ldu = U ? U->stepj() : 1;
    int ldv = V ? V->stepj() : 1;
    char vect = U ? V ? 'B' : 'Q' : V ? 'P' : 'N';
    float* VV = V ? V->ptr() : 0;
    float* UU = U ? U->ptr() : 0;

    LAPNAME(sgbbrd) (LAPCM LAPV(vect),LAPV(m),LAPV(n),LAPV(ncc),
	LAPV(kl),LAPV(ku),
	LAPP(A2.cptr()-A.nhi()),LAPV(lda),LAPP(D.ptr()),LAPP(E.ptr()),
	LAPP(UU),LAPV(ldu),LAPP(VV),LAPV(ldv),
	0,LAPV(n) LAPWK(work) LAPINFO LAP1 );

    det *= DiagMatrixViewOf(D).Det();
    LAP_Results("sgbbrd");
  }
  template <> inline void LapBidiagonalize(
      const GenBandMatrix<complex<float> >& A,
      const MatrixView<complex<float> >* U, const VectorView<float>& D,
      const VectorView<float>& E, const MatrixView<complex<float> >* V, 
      complex<float>& det)
  {
    TMVAssert(A.rowsize() == A.colsize());
    TMVAssert(A.rowsize() > 0);
    if (U) {
      TMVAssert(U->colsize() == A.colsize());
      TMVAssert(U->rowsize() == A.rowsize());
      TMVAssert(U->iscm());
      TMVAssert(U->ct() == NonConj);
    }
    if (V) {
      TMVAssert(V->colsize() == A.rowsize());
      TMVAssert(V->rowsize() == A.rowsize());
      TMVAssert(V->iscm());
      TMVAssert(V->ct() == NonConj);
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
    int lwork = max(m,n);
    complex<float>* work = LAP_CWork(lwork);
    float* rwork = LAP_SWork(lwork);
#endif
    BandMatrix<complex<float>,ColMajor> A2 = A;
    int lda = A2.diagstep();
    int ldu = U ? U->stepj() : 1;
    int ldv = V ? V->stepj() : 1;
    char vect = U ? V ? 'B' : 'Q' : V ? 'P' : 'N';
    complex<float>* VV = V ? V->ptr() : 0;
    complex<float>* UU = U ? U->ptr() : 0;

    LAPNAME(cgbbrd) (LAPCM LAPV(vect),LAPV(m),LAPV(n),LAPV(o),
	LAPV(kl),LAPV(ku),LAPP(A2.cptr()-A.nhi()),LAPV(lda),
	LAPP(D.ptr()),LAPP(E.ptr()),LAPP(UU),LAPV(ldu),
	LAPP(VV),LAPV(ldv),0,LAPV(n)
	LAPWK(work) LAPWK(rwork) LAPINFO LAP1);

    det *= DiagMatrixViewOf(A2.diag()).Det();
    LAP_Results("cgbbrd");
  }
#endif
#endif
  template <class T> inline void Bidiagonalize(
      const GenBandMatrix<T>& A,
      const MatrixView<T>* U, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, const MatrixView<T>* V,
      T& det)
  {
    TMVAssert(A.rowsize() <= A.colsize());
    TMVAssert(A.rowsize() > 0);
    if (U) {
      TMVAssert(U->colsize() == A.colsize());
      TMVAssert(U->rowsize() == A.rowsize());
      TMVAssert(U->ct() == NonConj);
    }
    if (V) {
      TMVAssert(V->colsize() == A.rowsize());
      TMVAssert(V->rowsize() == A.rowsize());
      TMVAssert(V->ct() == NonConj);
    }
    TMVAssert(D.size() == A.rowsize());
    TMVAssert(E.size() == D.size()-1);
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);
    TMVAssert(D.ct() == NonConj);
    TMVAssert(E.ct() == NonConj);

#ifdef XDEBUG
    BandMatrix<T> A0 = A;
    BandMatrix<T> A2 = A;
    Vector<RealType(T)> D2 = D;
    Vector<RealType(T)> E2 = E;
    Matrix<T> U2(A.colsize(),A.rowsize());
    MatrixView<T> U2v = U2.View();
    Matrix<T> V2(D.size(),D.size());
    MatrixView<T> V2v = V2.View();
    T det2(1);
    NonLapBidiagonalize(A2,&U2v,D2.View(),E2.View(),&V2v,det2);
#endif

    if (A.rowsize() > 0) {
      TMVAssert(E.size() == D.size()-1);
#ifdef LAP
      if (A.IsSquare() && (!U || U->iscm()) && (!V || V->iscm())) 
	LapBidiagonalize(A,U,D,E,V,det);
      else 
#endif
	NonLapBidiagonalize(A,U,D,E,V,det);
    }
#ifdef XDEBUG
    if (U && V) {
      Matrix<T> UBV = *U*UpperBiDiagMatrix(D,E)*(*V);
      if (Norm(UBV-A0) > 0.001*Norm(A0)) {
	cerr<<"Band Bidiagonalize:\n";
	cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
	cerr<<"-> D = "<<D<<endl;
	cerr<<"Nonlap D = "<<D2<<endl;
	cerr<<"Norm(diff) = "<<Norm(D-D2)<<endl;
	cerr<<"E = "<<E<<endl;
	cerr<<"Nonlap E = "<<E2<<endl;
	cerr<<"Norm(diff) = "<<Norm(E-E2)<<endl;
	cerr<<"U = "<<*U<<endl;
	cerr<<"V = "<<*V<<endl;
	cerr<<"U2 = "<<U2<<endl;
	cerr<<"V2 = "<<V2<<endl;
	cerr<<"A0 = "<<A0<<endl;
	cerr<<"UBV = "<<UBV<<endl;
	cerr<<"Norm(UBV-A0) = "<<Norm(UBV-A0)<<endl;
	abort();
      }
    }
#endif
  }

  template <class T> void BandSV_Decompose(
      const GenBandMatrix<T>& A,
      const MatrixView<T>* U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>* V, T& det)
  {
    // Decompose A into U S V
    // where S is a diagonal real matrix, and U,V are unitary matrices.
    // All matrices are square N x N
    // The determinant is kept track of in det.
    //
    // Everything is identical to the regular SVD except for the 
    // Bidiagonal Step.

    TMVAssert(A.rowsize() <= A.colsize());
    TMVAssert(A.rowsize() > 0);
    if (U) {
      TMVAssert(U->rowsize() == A.rowsize());
      TMVAssert(U->colsize() == A.colsize());
      TMVAssert(U->ct() == NonConj);
    } 
    if (V) {
      TMVAssert(V->rowsize() == A.rowsize());
      TMVAssert(V->colsize() == A.rowsize());
      TMVAssert(V->ct() == NonConj);
    }
    TMVAssert(S.size() == A.rowsize());
    TMVAssert(S.ct() == NonConj);

    if (A.nlo() == 0 && A.nhi() == 0) {
      if (U) U->SetToIdentity();
      if (V) V->SetToIdentity();

      det *= DiagMatrixViewOf(A.diag()).Det();
      const size_t N = A.rowsize();
      for(size_t j=0;j<N;j++) {
	T Ajj = A(j,j);
	S(j) = abs(Ajj);
	if(U) (*U)(j,j) = SIGN(Ajj,S(j));
      }
      auto_array<size_t> sortp(new size_t[N]);
      S.Sort(sortp.get(),DESCEND);
      if (U) U->PermuteCols(sortp.get());
      if (V) V->PermuteRows(sortp.get());
    } else {
      Vector<RealType(T)> E(S.size()-1);
      Bidiagonalize(A,U,S,E.View(),V,det);

      SV_Decompose_From_Bidiagonal(U,S,E.View(),V);
    }
  }

#define InstFile "TMV_BandSVDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


