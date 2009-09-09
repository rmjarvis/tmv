
#include "TMV_Band.h"
#include "TMV_SVDiv.h"
#include "TMV_Householder.h"
#include "TMV_Diag.h"

namespace tmv {

  template <class T> BandSVDiv<T>::BandSVDiv(const GenBandMatrix<T>& A,
      bool StoreU, bool StoreV) :
    istrans(A.colsize() < A.rowsize()),
    U(0), S(min(A.rowsize(),A.colsize())), V(0), det(T(1))
  {
    size_t M = max(A.colsize(),A.rowsize());
    size_t N = min(A.colsize(),A.rowsize());

    if (istrans) swap(StoreU,StoreV);

    if (StoreU) U = new Matrix<T,ColMajor>(M,N);
    if (StoreV) V = new Matrix<T,ColMajor>(N,N);
    MatrixView<T>* Uv = U ? new MatrixView<T>(U->View()) : 0;
    MatrixView<T>* Vv = V ? new MatrixView<T>(V->View()) : 0;

    if (istrans)
      BandSV_Decompose(A.Transpose(),Uv,S.View(),Vv,det);
    else
      BandSV_Decompose(A,Uv,S.View(),Vv,det);

    delete Uv;
    delete Vv;

    Thresh(Epsilon<T>());
  }

  template <class T> BandSVDiv<T>::~BandSVDiv() 
  { if (U) delete U; if (V) delete V; }

  template <class T> bool BandSVDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenBandMatrix<T>* bm = dynamic_cast<const GenBandMatrix<T>*>(&m);
    TMVAssert(bm);
    if (fout) {
      *fout << "M = "<<*bm<<endl;
      *fout << "U = "<<GetU()<<endl;
      *fout << "S = "<<GetS()<<endl;
      *fout << "V = "<<GetV()<<endl;
    }
    Matrix<T> m2 = GetU()*DiagMatrixViewOf(GetS())*GetV();
    RealType(T) nm = Norm(m2-(*bm));
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    if (fout) {
      *fout << "USV = "<<m2<<endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm<<endl;
    }
    return nm < S.size()*Epsilon<T>();
  }

  template <class T> void MakeBidiagReal(
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

  template <class T> void NonLapBidiagonalize(
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
      Matrix<T,ColMajor>* UU = (U==0) ? new Matrix<T,ColMajor>(A) : 0;
      if (U) *U = A;
      else U = new MatrixView<T>(UU->View());

      vector<size_t> vec(N), ver(N-1);
      size_t endcol = nlo+1;
      Vector<T> Ubeta(N);
      Vector<T> Vbeta(N-1);
      for(size_t j=0;j<N-1;++j) {
	vec[j] = endcol;
	size_t endrow = min(endcol+nhi,N);
	ver[j] = endrow;
	Ubeta(j) = Householder_Reflect(U->SubMatrix(j,endcol,j,endrow),det);
	if (endcol < M) endcol = min(endrow+nlo,M);
	Vbeta(j) = Householder_Reflect(U->Transpose().SubMatrix(
	      j+1,endrow,j,endcol),det);
      }
      vec[N-1] = endcol;
      Ubeta(N-1) = Householder_Reflect(U->SubMatrix(N-1,endcol,N-1,N),det);

      // Now U stores Householder vectors for U in lower diagonal columns (HLi)
      // and Householder vectors for V in upper diagonal rows (HRi)
      // except for the bidiagonal which is the bidiagonal we want:
      if (IsComplex(T())) {
	TMVAssert(NormInf(U->diag().Imag()) == RealType(T)(0));
	TMVAssert(NormInf(U->diag(1).Imag()) == RealType(T)(0));
      }
      D = U->diag().Real();
      E = U->diag(1).Real();

      if (V) {
	V->SetToIdentity();
	for (int j=N-2;j>=0;--j) {
	  V->row(j+1,j+2,ver[j]) = U->row(j,j+2,ver[j]);
	  Householder_Unpack(V->Transpose().SubMatrix(j+1,ver[j],j+1,N),
	      Vbeta(j));
	}
      }

      if (UU) { delete UU; delete U; U = 0; }
      else {
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
  template <class T> void LapBidiagonalize(
      const GenBandMatrix<T>& A,
      const MatrixView<T>* U, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, const MatrixView<T>* V,
      T& det)
  { NonLapBidiagonalize(A,U,D,E,V,det); }
  template <> void LapBidiagonalize(const GenBandMatrix<double>& A,
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
    int lwork = 2*max(m,n);
    double* work = LAP_DWork(lwork);
    int info;
    // LAP version overwrites original BandMatrix with crap.
    // Hence, copy BandMatrix before running.
    BandMatrix<double,ColMajor> A2 = A;
    int lda = A2.stepj()+1;
    int ldu = U ? U->stepj() : 0;
    int ldv = V ? V->stepj() : 0;
    if (U) {
      if (V) {
	char vect = 'B';
	dgbbrd(&vect,&m,&n,&ncc,&kl,&ku,const_cast<double*>(A2.cptr()-A.nhi()),
	    &lda,D.ptr(),E.ptr(),U->ptr(),&ldu,V->ptr(),&ldv,0,&n,work,&info);
      } else {
	char vect = 'Q';
	dgbbrd(&vect,&m,&n,&ncc,&kl,&ku,const_cast<double*>(A2.cptr()-A.nhi()),
	    &lda,D.ptr(),E.ptr(),U->ptr(),&ldu,0,&ldv,0,&n,work,&info);
      }
    } else {
      if (V) {
	char vect = 'P';
	dgbbrd(&vect,&m,&n,&ncc,&kl,&ku,const_cast<double*>(A2.cptr()-A.nhi()),
	    &lda,D.ptr(),E.ptr(),0,&ldu,V->ptr(),&ldv,0,&n,work,&info);
      } else {
	char vect = 'N';
	dgbbrd(&vect,&m,&n,&ncc,&kl,&ku,const_cast<double*>(A2.cptr()-A.nhi()),
	    &lda,D.ptr(),E.ptr(),0,&ldu,0,&ldv,0,&n,work,&info);
      }
    }
    det *= DiagMatrixViewOf(D).Det();
    LAP_Results(info,"dgbbrd");
  }
  template <> void LapBidiagonalize(
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
    int lwork = max(m,n);
    complex<double>* work = LAP_ZWork(lwork);
    double* rwork = LAP_DWork(lwork);
    int info;
    BandMatrix<complex<double>,ColMajor> A2 = A;
    int lda = A2.stepj()+1;
    int ldu = U ? U->stepj() : 0;
    int ldv = V ? V->stepj() : 0;
    if (U) {
      if (V) {
	char vect = 'B';
	zgbbrd(&vect,&m,&n,&o,&kl,&ku,LAP_Complex(A2.cptr()-A.nhi()),&lda,
	    D.ptr(),E.ptr(),LAP_Complex(U->ptr()),&ldu,LAP_Complex(V->ptr()),
	    &ldv,0,&n,LAP_Complex(work),rwork,&info);
      } else {
	char vect = 'Q';
	zgbbrd(&vect,&m,&n,&o,&kl,&ku,LAP_Complex(A2.cptr()-A.nhi()),&lda,
	    D.ptr(),E.ptr(),LAP_Complex(U->ptr()),&ldu,0,
	    &ldv,0,&n,LAP_Complex(work),rwork,&info);
      }
    } else {
      if (V) {
	char vect = 'P';
	zgbbrd(&vect,&m,&n,&o,&kl,&ku,LAP_Complex(A2.cptr()-A.nhi()),&lda,
	    D.ptr(),E.ptr(),0,&ldu,LAP_Complex(V->ptr()),
	    &ldv,0,&n,LAP_Complex(work),rwork,&info);
      } else {
	char vect = 'N';
	zgbbrd(&vect,&m,&n,&o,&kl,&ku,LAP_Complex(A2.cptr()-A.nhi()),&lda,
	    D.ptr(),E.ptr(),0,&ldu,0,&ldv,0,&n,LAP_Complex(work),rwork,&info);
      }
    }
    det *= DiagMatrixViewOf(A2.diag()).Det();
    LAP_Results(info,"zgbbrd");
  }
#ifndef NOFLOAT
  template <> void LapBidiagonalize(const GenBandMatrix<float>& A,
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
    int lwork = 2*max(m,n);
    float* work = LAP_SWork(lwork);
    int info;
    // LAP version overwrites original BandMatrix with crap.
    // Hence, copy BandMatrix before running.
    BandMatrix<float,ColMajor> A2 = A;
    int lda = A2.stepj()+1;
    int ldu = U ? U->stepj() : 0;
    int ldv = V ? V->stepj() : 0;
    if (U) {
      if (V) {
	char vect = 'B';
	sgbbrd(&vect,&m,&n,&ncc,&kl,&ku,const_cast<float*>(A2.cptr()-A.nhi()),
	    &lda,D.ptr(),E.ptr(),U->ptr(),&ldu,V->ptr(),&ldv,0,&n,work,&info);
      } else {
	char vect = 'Q';
	sgbbrd(&vect,&m,&n,&ncc,&kl,&ku,const_cast<float*>(A2.cptr()-A.nhi()),
	    &lda,D.ptr(),E.ptr(),U->ptr(),&ldu,0,&ldv,0,&n,work,&info);
      }
    } else {
      if (V) {
	char vect = 'P';
	sgbbrd(&vect,&m,&n,&ncc,&kl,&ku,const_cast<float*>(A2.cptr()-A.nhi()),
	    &lda,D.ptr(),E.ptr(),0,&ldu,V->ptr(),&ldv,0,&n,work,&info);
      } else {
	char vect = 'N';
	sgbbrd(&vect,&m,&n,&ncc,&kl,&ku,const_cast<float*>(A2.cptr()-A.nhi()),
	    &lda,D.ptr(),E.ptr(),0,&ldu,0,&ldv,0,&n,work,&info);
      }
    }
    det *= DiagMatrixViewOf(D).Det();
    LAP_Results(info,"sgbbrd");
  }
  template <> void LapBidiagonalize(
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
    int lwork = max(m,n);
    complex<float>* work = LAP_CWork(lwork);
    float* rwork = LAP_SWork(lwork);
    int info;
    BandMatrix<complex<float>,ColMajor> A2 = A;
    int lda = A2.stepj()+1;
    int ldu = U ? U->stepj() : 0;
    int ldv = V ? V->stepj() : 0;
    if (U) {
      if (V) {
	char vect = 'B';
	cgbbrd(&vect,&m,&n,&o,&kl,&ku,LAP_Complex(A2.cptr()-A.nhi()),&lda,
	    D.ptr(),E.ptr(),LAP_Complex(U->ptr()),&ldu,LAP_Complex(V->ptr()),
	    &ldv,0,&n,LAP_Complex(work),rwork,&info);
      } else {
	char vect = 'Q';
	cgbbrd(&vect,&m,&n,&o,&kl,&ku,LAP_Complex(A2.cptr()-A.nhi()),&lda,
	    D.ptr(),E.ptr(),LAP_Complex(U->ptr()),&ldu,0,
	    &ldv,0,&n,LAP_Complex(work),rwork,&info);
      }
    } else {
      if (V) {
	char vect = 'P';
	cgbbrd(&vect,&m,&n,&o,&kl,&ku,LAP_Complex(A2.cptr()-A.nhi()),&lda,
	    D.ptr(),E.ptr(),0,&ldu,LAP_Complex(V->ptr()),
	    &ldv,0,&n,LAP_Complex(work),rwork,&info);
      } else {
	char vect = 'N';
	cgbbrd(&vect,&m,&n,&o,&kl,&ku,LAP_Complex(A2.cptr()-A.nhi()),&lda,
	    D.ptr(),E.ptr(),0,&ldu,0,&ldv,0,&n,LAP_Complex(work),rwork,&info);
      }
    }
    det *= DiagMatrixViewOf(A2.diag()).Det();
    LAP_Results(info,"cgbbrd");
  }
#endif
#endif
  template <class T> void Bidiagonalize(
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

    if (A.rowsize() > 0) {
      TMVAssert(E.size() == D.size()-1);
#ifdef LAP
      if (A.IsSquare() && (!U || U->iscm()) && (!V || V->iscm())) 
	LapBidiagonalize(A,U,D,E,V,det);
      else 
#endif
	NonLapBidiagonalize(A,U,D,E,V,det);
    }
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

    Vector<RealType(T)> E(S.size()-1);
    Bidiagonalize(A,U,S,E.View(),V,det);

    SV_Decompose_From_Bidiagonal(U,S,E.View(),V);
  }

  template <class T> void BandSVDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  { 
    TMVAssert(U && V);
    if (istrans) {
      Matrix<T,ColMajor> SinvV = V->Conjugate().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      minv = U->Conjugate().Cols(0,kmax) * SinvV;
    } else  {
      Matrix<T,ColMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      minv = V->Adjoint().Cols(0,kmax) * SinvUt;
    }
  }

  template <class T> void BandSVDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(V);
    Matrix<T,ColMajor> SinvV = V->Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    if (istrans)
      minv = SinvV.Transpose() * SinvV.Conjugate();
    else
      minv = SinvV.Adjoint() * SinvV;
  }

  template <class T> void BandSVDiv<T>::Thresh(RealType(T) toler,
      ostream* debugout) const
  {
    TMVAssert(toler < RealType(T)(1));
    RealType(T) thresh = S(0)*toler;
    for(kmax=S.size(); kmax>0 && S(kmax-1)<=thresh; --kmax);
    if(debugout) {
      (*debugout)<<"S = "<<S<<endl;
      (*debugout)<<"Smax = "<<S(0)<<", thresh = "<<thresh<<std::endl;
      (*debugout)<<"kmax = "<<kmax<<" (S.size = "<<S.size()<<")"<<endl;
    }
  }

  template <class T> void BandSVDiv<T>::Top(size_t neigen,
      ostream* debugout) const
  {
    TMVAssert(neigen > 0);
    TMVAssert(neigen <= S.size());
    kmax = neigen;
    if(debugout) {
      (*debugout)<<"S = "<<S<<endl;
      (*debugout)<<"kmax = "<<kmax<<" (S.size = "<<S.size()<<")"<<endl;
    }
  }

#define InstFile "TMV_BandSVDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


