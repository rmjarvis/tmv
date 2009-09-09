
#include "TMV_Band.h"
#include "TMV_SVDiv.h"
#include "TMV_Householder.h"
#include "TMV_Diag.h"

namespace tmv {

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
      const MatrixView<T>& U, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, const MatrixView<T>& V, T& det)
  {
    // Decompose A into U B V
    // The Bidiagonal Matrix B is stored as two vectors: D, E
    // D is the diagonal, E is the super-diagonal
    // We use Householder reflections to reduce A to the bidiagonal form:

    TMVAssert(A.rowsize() <= A.colsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(U.colsize() == A.colsize());
    TMVAssert(U.rowsize() == A.rowsize());
    TMVAssert(V.colsize() == A.rowsize());
    TMVAssert(V.rowsize() == A.rowsize());
    TMVAssert(D.size() == A.rowsize());
    TMVAssert(E.size() == D.size()-1);

    const size_t M = A.colsize();
    const size_t N = A.rowsize();

    const int nlo = A.nlo();
    const int nhi = A.nhi();

    if (nlo == 0 && nhi == 1) {
      U.Zero();
      V.Zero();
      VectorView<T> Udiag = U.diag();
      VectorView<T> Vdiag = V.diag();
      MakeBidiagReal(U.diag(),V.diag(),A.diag(),A.diag(1),D,E,det);
    } else if (A.IsSquare() && nlo == 1 && nhi == 0) {
      U.Zero();
      V.Zero();
      VectorView<T> Udiag = U.SubVector(N-1,0,-1,1,N);
      VectorView<T> Vdiag = V.SubVector(0,N-1,1,-1,N);
      MakeBidiagReal(Udiag,Vdiag,A.diag().Reverse(),A.diag(-1).Reverse(),
	  D,E,det);
    } else {
      U = A;
      vector<size_t> vec(N), ver(N-1);
      size_t endcol = nlo+1;
      Vector<T> Ubeta(N);
      Vector<T> Vbeta(N-1);
      for(size_t j=0;j<N-1;++j) {
	vec[j] = endcol;
	size_t endrow = min(endcol+nhi,N);
	ver[j] = endrow;
	Ubeta(j) = Householder_Reflect(U.SubMatrix(j,endcol,j,endrow),det);
	if (endcol < M) endcol = min(endrow+nlo,M);
	Vbeta(j) = Householder_Reflect(U.Transpose().SubMatrix(j+1,endrow,j,endcol),det);
      }
      vec[N-1] = endcol;
      Ubeta(N-1) = Householder_Reflect(U.SubMatrix(N-1,endcol,N-1,N),det);

      // Now U stores Householder vectors for U in lower diagonal columns (HLi)
      // and Householder vectors for V in upper diagonal rows (HRi)
      // except for the bidiagonal which is the bidiagonal we want:
      if (IsComplex(T())) {
	TMVAssert(NormInf(U.diag().Imag()) == RealType(T)(0));
	TMVAssert(NormInf(U.diag(1).Imag()) == RealType(T)(0));
      }
      D = U.diag().Real();
      E = U.diag(1).Real();
      U.diag().Zero();
      U.diag(1).Zero();

      // B = HLn-1 ... HL1 HL0 A HR0T HR1T ... HRn-2T
      // Using the fact that H = Ht = H^-1, we get:
      // U = HL0 ... HLn-1  and  VT = HR0 ... HRn-2
      V.SetToIdentity();
      Householder_Unpack(U.SubMatrix(N-1,vec[N-1],N-1,N),Ubeta(N-1));
      for (int j=N-2;j>=0;--j) {
	V.row(j+1,j+2,ver[j]) = U.row(j,j+2,ver[j]);
	Householder_Unpack(V.Transpose().SubMatrix(j+1,ver[j],j+1,N),Vbeta(j));
	U.row(j,j+2,ver[j]).Zero();
	Householder_Unpack(U.SubMatrix(j,vec[j],j,N),Ubeta(j));
      }
    }
    det *= DiagMatrixViewOf(D).Det();
  }

#ifdef LAP
  template <class T> void LapBidiagonalize(
      const GenBandMatrix<T>& A,
      const MatrixView<T>& U, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, const MatrixView<T>& V,
      T& det)
  { NonLapBidiagonalize(A,U,D,E,V,det); }
  template <> void LapBidiagonalize(const GenBandMatrix<double>& A,
      const MatrixView<double>& U, const VectorView<double>& D,
      const VectorView<double>& E, const MatrixView<double>& V, 
      double& det)
  {
    TMVAssert(A.rowsize() == A.colsize());
    // The Lap routines can do NonSquare matrices, but they want to
    // write out to a square (MxM) U matrix which is larger than
    // what we have stored here.
    TMVAssert(A.rowsize() > 0);
    TMVAssert(U.colsize() == A.colsize());
    TMVAssert(U.rowsize() == A.rowsize());
    TMVAssert(V.colsize() == A.rowsize());
    TMVAssert(V.rowsize() == A.rowsize());
    TMVAssert(D.size() == A.rowsize());
    TMVAssert(E.size() == D.size()-1);
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);
    TMVAssert(U.ct() == NonConj);
    TMVAssert(D.ct() == NonConj);
    TMVAssert(E.ct() == NonConj);
    TMVAssert(V.ct() == NonConj);

    char vect = 'B';
    int m = A.colsize();
    int n = A.rowsize();
    int ncc = 0;
    int kl = A.nlo();
    int ku = A.nhi();
    int ldu = U.stepj();
    int ldv = V.stepj();
    int lwork = 2*max(m,n);
    double* work = LAP_DWork(lwork);
    int info;
    // LAP version overwrites original BandMatrix with crap.
    // Hence, copy BandMatrix before running.
    BandMatrix<double,ColMajor> A2 = A;
    int lda = A2.stepj()+1;
    dgbbrd(&vect,&m,&n,&ncc,&kl,&ku,const_cast<double*>(A2.cptr()-A.nhi()),&lda,
	D.ptr(),E.ptr(),U.ptr(),&ldu,V.ptr(),&ldv,0,&n,work,&info);
    det *= DiagMatrixViewOf(D).Det();
    LAP_Results(info,"dgbbrd");
  }
  template <> void LapBidiagonalize(
      const GenBandMatrix<complex<double> >& A,
      const MatrixView<complex<double> >& U, const VectorView<double>& D,
      const VectorView<double>& E, const MatrixView<complex<double> >& V, 
      complex<double>& det)
  {
    TMVAssert(A.rowsize() == A.colsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(U.colsize() == A.colsize());
    TMVAssert(U.rowsize() == A.rowsize());
    TMVAssert(V.colsize() == A.rowsize());
    TMVAssert(V.rowsize() == A.rowsize());
    TMVAssert(D.size() == A.rowsize());
    TMVAssert(E.size() == D.size()-1);
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);
    TMVAssert(U.ct() == NonConj);
    TMVAssert(D.ct() == NonConj);
    TMVAssert(E.ct() == NonConj);
    TMVAssert(V.ct() == NonConj);

    char vect = 'B';
    int m = A.colsize();
    int n = A.rowsize();
    int o = 0;
    int kl = A.nlo();
    int ku = A.nhi();
    int ldu = U.stepj();
    int ldv = V.stepj();
    int lwork = max(m,n);
    complex<double>* work = LAP_ZWork(lwork);
    double* rwork = LAP_DWork(lwork);
    int info;
    BandMatrix<complex<double>,ColMajor> A2 = A;
    int lda = A2.stepj()+1;
    zgbbrd(&vect,&m,&n,&o,&kl,&ku,LAP_Complex(A2.cptr()-A.nhi()),&lda,
	D.ptr(),E.ptr(),LAP_Complex(U.ptr()),&ldu,LAP_Complex(V.ptr()),&ldv,
	0,&n,LAP_Complex(work),rwork,&info);
    det *= DiagMatrixViewOf(A2.diag()).Det();
    LAP_Results(info,"zgbbrd");
  }
#ifndef NOFLOAT
  template <> void LapBidiagonalize(const GenBandMatrix<float>& A,
      const MatrixView<float>& U, const VectorView<float>& D,
      const VectorView<float>& E, const MatrixView<float>& V, 
      float& det)
  {
    TMVAssert(A.rowsize() == A.colsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(U.colsize() == A.colsize());
    TMVAssert(U.rowsize() == A.rowsize());
    TMVAssert(V.colsize() == A.rowsize());
    TMVAssert(V.rowsize() == A.rowsize());
    TMVAssert(D.size() == A.rowsize());
    TMVAssert(E.size() == D.size()-1);
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);
    TMVAssert(U.ct() == NonConj);
    TMVAssert(D.ct() == NonConj);
    TMVAssert(E.ct() == NonConj);
    TMVAssert(V.ct() == NonConj);

    char vect = 'B';
    int m = A.colsize();
    int n = A.rowsize();
    int o = 0;
    int kl = A.nlo();
    int ku = A.nhi();
    int ldu = U.stepj();
    int ldv = V.stepj();
    int lwork = 2*max(m,n);
    float* work = LAP_SWork(lwork);
    int info;
    BandMatrix<float,ColMajor> A2 = A;
    int lda = A2.stepj()+1;
    sgbbrd(&vect,&m,&n,&o,&kl,&ku,const_cast<float*>(A2.cptr()-A.nhi()),&lda,
	D.ptr(),E.ptr(),U.ptr(),&ldu,V.ptr(),&ldv,0,&n,work,&info);
    det *= DiagMatrixViewOf(D).Det();
    LAP_Results(info,"sgbbrd");
  }
  template <> void LapBidiagonalize(
      const GenBandMatrix<complex<float> >& A,
      const MatrixView<complex<float> >& U, const VectorView<float>& D,
      const VectorView<float>& E, const MatrixView<complex<float> >& V, 
      complex<float>& det)
  {
    TMVAssert(A.rowsize() == A.colsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(U.colsize() == A.colsize());
    TMVAssert(U.rowsize() == A.rowsize());
    TMVAssert(V.colsize() == A.rowsize());
    TMVAssert(V.rowsize() == A.rowsize());
    TMVAssert(D.size() == A.rowsize());
    TMVAssert(E.size() == D.size()-1);
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);
    TMVAssert(U.ct() == NonConj);
    TMVAssert(D.ct() == NonConj);
    TMVAssert(E.ct() == NonConj);
    TMVAssert(V.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);

    char vect = 'B';
    int m = A.colsize();
    int n = A.rowsize();
    int o = 0;
    int kl = A.nlo();
    int ku = A.nhi();
    int ldu = U.stepj();
    int ldv = V.stepj();
    int lwork = max(m,n);
    complex<float>* work = LAP_CWork(lwork);
    float* rwork = LAP_SWork(lwork);
    int info;
    BandMatrix<complex<float>,ColMajor> A2 = A;
    int lda = A2.stepj()+1;
    cgbbrd(&vect,&m,&n,&o,&kl,&ku,LAP_Complex(A2.cptr()-A.nhi()),&lda,
	D.ptr(),E.ptr(),LAP_Complex(U.ptr()),&ldu,LAP_Complex(V.ptr()),&ldv,
	0,&n,LAP_Complex(work),rwork,&info);
    det *= DiagMatrixViewOf(A2.diag()).Det();
    LAP_Results(info,"cgbbrd");
  }
#endif
#endif
  template <class T> void Bidiagonalize(
      const GenBandMatrix<T>& A,
      const MatrixView<T>& U, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, const MatrixView<T>& V,
      T& det)
  {
    TMVAssert(A.rowsize() <= A.colsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(U.colsize() == A.colsize());
    TMVAssert(U.rowsize() == A.rowsize());
    TMVAssert(V.colsize() == A.rowsize());
    TMVAssert(V.rowsize() == A.rowsize());
    TMVAssert(D.size() == A.rowsize());
    TMVAssert(E.size() == D.size()-1);
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);
    TMVAssert(U.ct() == NonConj);
    TMVAssert(D.ct() == NonConj);
    TMVAssert(E.ct() == NonConj);
    TMVAssert(V.ct() == NonConj);

    if (U.rowsize() > 0) {
      TMVAssert(E.size() == D.size()-1);
#ifdef LAP
      if (A.IsSquare() && U.iscm() && V.iscm()) 
	LapBidiagonalize(A,U,D,E,V,det);
      else 
#endif
	NonLapBidiagonalize(A,U,D,E,V,det);
    }
  }

  template <class T> void BandSV_Decompose(
      const GenBandMatrix<T>& A,
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>& V, T& det)
  {
    // Decompose A into U S V
    // where S is a diagonal real matrix, and U,V are unitary matrices.
    // All matrices are square N x N
    // The determinant is kept track of in det.
    //
    // Everything is identical to the regular SVD except for the 
    // Bidiagonal Step.

    TMVAssert(U.rowsize() == A.rowsize());
    TMVAssert(U.colsize() == A.colsize());
    TMVAssert(V.rowsize() == A.rowsize());
    TMVAssert(V.colsize() == A.rowsize());
    TMVAssert(S.size() == A.rowsize());
    TMVAssert(A.rowsize()>0);
    TMVAssert(U.ct() == NonConj);
    TMVAssert(S.ct() == NonConj);
    TMVAssert(V.ct() == NonConj);

    Vector<RealType(T)> E(S.size()-1);
    Bidiagonalize(A,U,S,E.View(),V,det);

    SV_Decompose_From_Bidiagonal(U,S,E.View(),V);
  }

  template <class T> Matrix<T,ColMajor> BandSVDiv<T>::Inverse() const
  { 
    if (istrans) {
      Matrix<T,ColMajor> SinvV = V.QuickConjugate().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      return U.QuickConjugate().Cols(0,kmax) * SinvV;
    } else  {
      Matrix<T,ColMajor> SinvUt = U.QuickAdjoint().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      return V.QuickAdjoint().Cols(0,kmax) * SinvUt;
    }
  }

  template <class T> Matrix<T,ColMajor> BandSVDiv<T>::DoInverseATA() const
  {
    Matrix<T,ColMajor> SinvV = V.Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    if (istrans)
      return SinvV.QuickTranspose() * SinvV.QuickConjugate();
    else
      return SinvV.QuickAdjoint() * SinvV;
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


