
#include "TMV_Sym.h"
#include "TMV.h"
#include "TMV_Givens.h"
#include "TMV_Diag.h"
#include "TMV_Band.h"
#include "TMV_Householder.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t SYM_TRIDIAG_BLOCKSIZE = TMV_BLOCKSIZE;
#else
  const size_t SYM_TRIDIAG_BLOCKSIZE = 64;
#endif

  template <class T> const MatrixView<T>* ZMV() 
  { return (const MatrixView<T>*)(0); }

  template <class T> HermSVDiv<T>::HermSVDiv(const GenSymMatrix<T>& A,
      bool StoreU) :
    U(0), S(A.size()), det(RealType(T)(0)), calcdet(false)
  {
    size_t N = A.size();

    U = new Matrix<T,ColMajor>(N,N);
    LowerTriMatrixViewOf(*U) = A.LowerTri();
    HermSV_Decompose(U->View(),S.View(),StoreU);
    if (!StoreU) { delete U; U = 0; }

    Thresh(Epsilon<T>());
  }

  template <class T> HermSVDiv<T>::~HermSVDiv() 
  { if (U) delete U; }

  template <class T> SymSVDiv<T>::SymSVDiv(const GenSymMatrix<T>& A, 
      bool StoreUV) :
    U(0), S(A.size()), V(0), det(T(1))
  {
    TMVAssert(IsComplex(T()));
    size_t N = A.size();

    U = new Matrix<T,ColMajor>(N,N);
    LowerTriMatrixViewOf(*U) = A.LowerTri();
    if (StoreUV) {
      V = new Matrix<T,ColMajor>(N,N);
      MatrixView<T> VV = V->View();
      SymSV_Decompose(U->View(),S.View(),&VV,det);
    } else {
      SymSV_Decompose(U->View(),S.View(),ZMV<T>(),det);
      delete U; U = 0;
    }

    Thresh(Epsilon<T>());
  }

  template <class T> SymSVDiv<T>::~SymSVDiv() 
  { if (U) delete U; if (V) delete V; }

  template <class T> bool SymSVDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenSymMatrix<T>* sm = dynamic_cast<const GenSymMatrix<T>*>(&m);
    TMVAssert(sm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*sm)<<"  "<<*sm<<endl;
      *fout << "U = "<<GetU()<<endl;
      *fout << "S = "<<GetS()<<endl;
      *fout << "V = "<<GetV()<<endl;
    }
    Matrix<T> usv = GetU()*DiagMatrixViewOf(GetS())*GetV();
    RealType(T) nm = Norm(usv-*sm);
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    if (fout) {
      *fout << "USV = "<<usv<<endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm;
      *fout <<"  "<<Condition()<<" * "<<Epsilon<T>()<<endl;
    }
    return nm < Condition()*sm->colsize()*Epsilon<T>();
  }

  template <class T> bool HermSVDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenSymMatrix<T>* sm = dynamic_cast<const GenSymMatrix<T>*>(&m);
    TMVAssert(sm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*sm)<<"  "<<*sm<<endl;
      *fout << "U = "<<GetU()<<endl;
      *fout << "S = "<<GetS()<<endl;
      *fout << "V = "<<GetV()<<endl;
    }
    Matrix<T> usv = GetU()*DiagMatrixViewOf(GetS())*GetV();
    RealType(T) nm = Norm(usv-*sm);
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    if (fout) {
      *fout << "USV = "<<usv<<endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm;
      *fout <<"  "<<Condition()<<" * "<<Epsilon<T>()<<endl;
    }
    return nm < Condition()*sm->colsize()*Epsilon<T>();
  }

  template <class T, class Td> void NonBlockSymTridiagonalize(
      const SymMatrixView<T>& A, const VectorView<T>& beta,
      const VectorView<Td>& D, const VectorView<RealType(T)>& E, 
      T& det, size_t j1=0)
  {
    // Decompose A into U T Ut
    // The Tridiagonal Matrix T is stored as two vectors: D, E
    // D is the diagonal, E is the sub-diagonal.
    // If A is herm, E* is the super-diagonal, otherwise E is.
    // However, the Householder reflections make E real, so this
    // distinction is irrelevant.
    // A along with beta hold the U matrix.
    const size_t N = A.size();

    TMVAssert(N > 0);
    TMVAssert(beta.size() == N-1);
    TMVAssert(D.size() == N);
    TMVAssert(E.size() == N-1);
    TMVAssert(IsReal(Td()) || !A.isherm());
    TMVAssert(A.uplo() == Lower);

    // We use Householder reflections to reduce A to the tridiagonal form:
    for(size_t j=j1;j<N-1;++j) {
      beta(j) = Householder_Reflect(A.col(j,j+1,N),det);
      if (beta(j) != T(0)) 
	Householder_LRMult(A.col(j,j+2,N),T(beta(j)),A.SubSymMatrix(j+1,N));
    }

    // The tridiagonal of A is the tridiagonal we want, so copy it to D,E
    if (IsReal(Td())) D = A.diag().Real();
    else D = A.diag();
    E = A.diag(1).Real();

    if (!A.isherm()) det *= det;
  }

  template <class T, class Td> void BlockSymTridiagonalize(
      const SymMatrixView<T>& A, const VectorView<T>& beta,
      const VectorView<Td>& D, const VectorView<RealType(T)>& E, T& det)
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

    const size_t N = A.size();
    const size_t Nm1 = N-1;
    const RealType(T) half(0.5);

    TMVAssert(N > 0);
    TMVAssert(beta.size() == N-1);
    TMVAssert(D.size() == N);
    TMVAssert(E.size() == N-1);
    TMVAssert(IsReal(Td()) || !A.isherm());
    TMVAssert(A.uplo() == Lower);

    Matrix<T,ColMajor> ZYtmT(N,min(SYM_TRIDIAG_BLOCKSIZE,Nm1));
    // It's easier (for me) to think about the calculations with ZYtm.
    // But the eventual Rank2KUpdate call will have unconjugated matrices
    // if we store (ZYtm)t for Hermitian or (ZYtm)T for symmetric.
    // This is what is stored in the Matrix ZYtmT.
    // 
    for(size_t j1=0;j1<Nm1;) {
      size_t j2 = min(Nm1,j1+SYM_TRIDIAG_BLOCKSIZE);

      if (j2 < N) {
	for(size_t j=j1,jj=0;j<j2;++j,++jj) { // jj = j-j1

	  // Update current column:
	  // A(j:N,j) -= Y(j:N,0:j) ZYtm(0:j,j) + (ZYtm)t(j:N,0:j) Yt(0:j,j)
	  VectorView<T> Acolj = A.col(j,j,N);
	  MatrixView<T> Y0 = A.SubMatrix(j,N,j1,j);
	  MatrixView<T> Y0T = (A.issym() ? Y0.Transpose() : Y0.Adjoint());
	  MatrixView<T> ZYtm0T = ZYtmT.SubMatrix(j,N,0,jj);
	  MatrixView<T> ZYtm0 = (A.issym() ? ZYtm0T.Transpose() : 
	      ZYtm0T.Adjoint());
	  if (jj > 0) {
	    Acolj -= Y0 * ZYtm0.col(0);
	    Acolj -= ZYtm0T * Y0T.col(0);
	  }

	  // Do the Householder reflection
	  const size_t Nmj = Acolj.size();
	  VectorView<T> u = Acolj.SubVector(1,Nmj);
	  T b = Householder_Reflect(u,det);
	  beta(j) = b;
	  TMVAssert(IMAG(u(0)) == RealType(T)(0));
	  E(j) = REAL(u(0));
	  // Set top of u to 1 explicitly.  Calculations below are easier this 
	  // way, rather than dealing with it being implicitly 1.
	  u(0) = T(1);

	  // Update ZYtmT
	  // ZYtm(j,j+1:N)  = Z(j,0:j+1) Yt(0:j+1,0:N) m'(0:N,j+1:N)
	  //                = beta Yt(j,0:N) m'(0:N,j+1:N)
	  //                = beta Yt(j,j+1:N) m'(j+1:N,j+1:N)
	  //                = beta ut m'(j+1:N,j+1:N)
	  // ZYtmT(j+1:N,j) = bb m'(j+1:N,j+1:N) uu
	  //                = bb [ A(j+1:N,j+1:N) uu
	  //                       - Y(j+1:N,0:j) ZYtm(0:j,j+1:N) uu
	  //                       - ZYtmT(j+1:N,0:j) YT(0:j,j+1:N) uu ]
	  //                = bb [ A(j+1:N,j+1:N) uu
	  //                       - Y1 ZYtm1 uu - ZYtm1T Y1T uu ]
	  // where {bb = beta*, uu = u} if hermitian
	  // and   {bb = beta, uu = u*} if symmetric

	  MatrixView<T> ZYtm1 = ZYtm0.Cols(1,Nmj);
	  MatrixView<T> ZYtm1T = ZYtm0T.Rows(1,Nmj);
	  MatrixView<T> Y1 = Y0.Rows(1,Nmj);
	  MatrixView<T> Y1T = Y0T.Cols(1,Nmj);
	  VectorView<T> ZYtmTj = ZYtmT.col(jj,j+1,N);
	  VectorView<T> temp = ZYtmT.col(jj,j1,j);
	  VectorView<T> uu = A.isherm() ? u : u.Conjugate();
	  T bb = A.issym() ? b : CONJ(b);

	  ZYtmTj = A.SubSymMatrix(j+1,N) * uu;
	  if (jj > 0) {
	    ZYtmTj -= Y1 * (temp = ZYtm1 * uu);
	    ZYtmTj -= ZYtm1T * (temp = Y1T * uu);
	  }
	  ZYtmTj *= bb;

	  // There is an extra complication as well now.
	  // What we really want to do to block update the matrix would be
	  // m = m - YZYtm - mYZtYt + YZYtmYZtYt
	  // However, we can fold that last term into the first two
	  // by making the ZYtm that we store really equal to
	  // ZYtm - 1/2 ZYtmYZtYt
	  VectorView<T> ZYtmj = A.issym() ? ZYtmTj : ZYtmTj.Conjugate();
	  ZYtmTj -= (half*bb)*(ZYtmj*uu)*u;
	}

	// Update the rest of the matrix:
	if (j2 < N)
	  Rank2KUpdate(T(-1),A.SubMatrix(j2,N,j1,j2),
	      ZYtmT.SubMatrix(j2,N,0,j2-j1),A.SubSymMatrix(j2,N));

      } else NonBlockSymTridiagonalize(A,beta,D,E,det,j1);

      j1 = j2;
    }

    // E has already been set, but the diagonal of A needs to be copied to D.
    if (IsReal(Td())) D = A.diag().Real();
    else D = A.diag();

    if (!A.isherm()) det *= det;
  }

  template <class T, class Td> void NonLapSymTridiagonalize(
      const SymMatrixView<T>& A, const VectorView<T>& beta,
      const VectorView<Td>& D, const VectorView<RealType(T)>& E, T& det)
  {
    if (A.size() >= SYM_TRIDIAG_BLOCKSIZE*2)
      BlockSymTridiagonalize(A,beta,D,E,det);
    else
      NonBlockSymTridiagonalize(A,beta,D,E,det);
  }

#ifdef LAP
  template <class T, class Td> void LapSymTridiagonalize(
      const SymMatrixView<T>& A, const VectorView<T>& beta,
      const VectorView<Td>& D, const VectorView<RealType(T)>& E, T& det)
  { NonLapSymTridiagonalize(A,beta,D,E,det); }
  template <> void LapSymTridiagonalize(
      const SymMatrixView<double>& A, const VectorView<double>& beta,
      const VectorView<double>& D, const VectorView<double>& E, double& det)
  {
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.uplo() == Lower);
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(beta.step() == 1);
    TMVAssert(beta.size() == A.size()-1);
    TMVAssert(D.size() == A.size());
    TMVAssert(E.size()+1 == D.size());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(beta.ct()==NonConj);

    char u = A.iscm() ? 'L' : 'U';
    int n = A.size();
    int ldu = A.iscm() ? A.stepj() : A.stepi();
    int lwork = n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    dsytrd(&u,&n,A.ptr(),&ldu,D.ptr(),E.ptr(),beta.ptr(),
	work,&lwork,&info);
    LAP_Results(info,int(work[0]),n,n,lwork,"dsytrd");
  }
  template <> void LapSymTridiagonalize(
      const SymMatrixView<complex<double> >& A,
      const VectorView<complex<double> >& beta,
      const VectorView<double>& D, const VectorView<double>& E, 
      complex<double>& det)
  {
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isherm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(beta.step() == 1);
    TMVAssert(beta.size() == A.size()-1);
    TMVAssert(D.size() == A.size());
    TMVAssert(E.size()+1 == D.size());
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(beta.ct()==NonConj);

    char u = A.iscm() ? 'L' : 'U';
    int n = A.size();
    int ldu = A.iscm() ? A.stepj() : A.stepi();
    int lwork = n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    zhetrd(&u,&n,LAP_Complex(A.ptr()),&ldu,D.ptr(),E.ptr(),
	LAP_Complex(beta.ptr()),LAP_Complex(work),&lwork,&info);
    if (A.iscm() != A.isconj()) beta.ConjugateSelf();
    LAP_Results(info,int(real(work[0])),n,n,lwork,"zhetrd");
  }
#ifndef NOFLOAT
  template <> void LapSymTridiagonalize(
      const SymMatrixView<float>& A, const VectorView<float>& beta,
      const VectorView<float>& D, const VectorView<float>& E, float& det)
  {
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.uplo() == Lower);
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(beta.step() == 1);
    TMVAssert(beta.size() == A.size()-1);
    TMVAssert(D.size() == A.size());
    TMVAssert(E.size()+1 == D.size());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(beta.ct()==NonConj);

    char u = A.iscm() ? 'L' : 'U';
    int n = A.size();
    int ldu = A.iscm() ? A.stepj() : A.stepi();
    int lwork = n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    ssytrd(&u,&n,A.ptr(),&ldu,D.ptr(),E.ptr(),beta.ptr(),
	work,&lwork,&info);
    LAP_Results(info,int(work[0]),n,n,lwork,"ssytrd");
  }
  template <> void LapSymTridiagonalize(
      const SymMatrixView<complex<float> >& A,
      const VectorView<complex<float> >& beta,
      const VectorView<float>& D, const VectorView<float>& E, 
      complex<float>& det)
  {
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isherm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(beta.step() == 1);
    TMVAssert(beta.size() == A.size()-1);
    TMVAssert(D.size() == A.size());
    TMVAssert(E.size()+1 == D.size());
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(beta.ct()==NonConj);

    char u = A.iscm() ? 'L' : 'U';
    int n = A.size();
    int ldu = A.iscm() ? A.stepj() : A.stepi();
    int lwork = n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    chetrd(&u,&n,LAP_Complex(A.ptr()),&ldu,D.ptr(),E.ptr(),
	LAP_Complex(beta.ptr()),LAP_Complex(work),&lwork,&info);
    if (A.iscm() != A.isconj()) beta.ConjugateSelf();
    LAP_Results(info,int(real(work[0])),n,n,lwork,"chetrd");
  }
#endif // NOFLOAT
#endif // LAP

  template <class T, class Td> void SymTridiagonalize(
      const SymMatrixView<T>& A, const VectorView<T>& beta,
      const VectorView<Td>& D, const VectorView<RealType(T)>& E, T& det)
  {
    TMVAssert(A.size() == D.size());
    TMVAssert(beta.size() == A.size()-1);
    TMVAssert(A.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(IsReal(Td()) || !A.isherm());

    // If Hermitian, then the determinant is unaffected by Householder
    // multiplies, so don't bother keeping track of d.
    T d = A.isherm() ? T(0) : T(1);

    if (A.size() > 0) {
      TMVAssert(E.size()+1 == D.size());
      if (A.uplo() == Upper) {
	if (A.isherm()) SymTridiagonalize(A.Adjoint(),beta,D,E,d);
	else SymTridiagonalize(A.Transpose(),beta,D,E,d);
      } else {
#ifdef XDEBUG
	Matrix<T> A0 = A;
	Matrix<T,ColMajor> A2 = A;
	Vector<T> b2 = beta;
	Vector<Td> D2 = D;
	Vector<RealType(T)> E2 = E;
	T d2 = d;
	if (IsReal(T()) || A.isherm())
	  NonBlockSymTridiagonalize(HermMatrixViewOf(A2,Lower),b2.View(),
	      D2.View(),E2.View(),d2);
	else
	  NonBlockSymTridiagonalize(SymMatrixViewOf(A2,Lower),b2.View(),
	      D2.View(),E2.View(),d2);
#endif // XDEBUG
#ifdef LAP
	if ((IsReal(T()) || A.isherm()) && (A.iscm() || A.isrm()) &&
	    D.step() == 1 && E.step() == 1 && beta.step() == 1) {
	  TMVAssert(IsReal(Td()));
	  LapSymTridiagonalize(A,beta,D,E,d);
	}
	else 
#endif // LAP
	  NonLapSymTridiagonalize(A,beta,D,E,d);
#ifdef XDEBUG
	const size_t N = A.size();
	Matrix<T> AA(N,N,T(0));
	for(size_t j=N-1;j>0;--j) AA.col(j,j,N) = A.col(j-1,j,N);
	AA(0,0) = T(1);
	GetQFromQR(AA.SubMatrix(1,N,1,N),beta);
	Matrix<T> TT(N,N,T(0));
	TT.diag() = D;
	TT.diag(1) = TT.diag(-1) = E;
	Matrix<T> A3 = AA*TT*(A.isherm() ? AA.Adjoint() : AA.Transpose());
	if (Norm(A3-A0) > 0.001*max(RealType(T)(1),Norm(A0))) {
	  cerr<<"SymTridiagonalize: \n";
	  cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
	  cerr<<"Done: A = "<<A<<endl;
	  cerr<<"beta = "<<beta<<endl;
	  cerr<<"UU = "<<AA<<endl;
	  cerr<<"D = "<<D<<endl;
	  cerr<<"E = "<<E<<endl;
	  cerr<<"TT = "<<TT<<endl;
	  cerr<<"UU * TT * UUt = "<<A3<<endl;
	  cerr<<"NonBlock versions:\n";
	  if (IsReal(T()) || A.isherm())
	    cerr<<"A2 = "<<HermMatrixViewOf(A2,Lower)<<endl;
	  else
	    cerr<<"A2 = "<<SymMatrixViewOf(A2,Lower)<<endl;
	  cerr<<"Norm(A2-A) = "<<Norm(A2-A)<<endl;
	  cerr<<"beta2 = "<<b2<<endl;
	  cerr<<"beta = "<<beta<<endl;
	  cerr<<"Norm(beta2-beta) = "<<Norm(b2-beta)<<endl;
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
    if (!A.isherm()) det *= d;
  }

  template <class T> void HermTridiagonalChopSmallElements(
      const VectorView<T>& D, const VectorView<T>& E)
  {
    // This routines sets to 0 any elements in E or D which
    // are essentially 0, given the machine precision:
    // if |E(i)| < Epsilon() * (|D(i)| + |D(i+1)|) then E(i) <- 0
    // if |D(i)| < Epsilon() * |T| then D(i) <- 0
    TMVAssert(IsReal(T()));
    TMVAssert(E.size() == D.size()-1);
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    RealType(T) eps = Epsilon<T>();
    RealType(T) dthresh = eps * SQRT(NormSq(D) + 2*NormSq(E));

    T* Di = D.ptr();
    T* Ei = E.ptr();
    for(size_t k=E.size();k>0;--k,++Di,++Ei) {
      if (abs(*Di) < dthresh) *Di = T(0);
      if (abs(*Ei) < eps * (abs(*Di)+abs(*(Di+1)))) *Ei = T(0);
    }
    if (abs(*Di) < dthresh) *Di = T(0);
  }

  template <class T> T HermTridiagonalTrailingEigenValue(
      const VectorView<T>& D, const VectorView<T>& E)
  {
    // Return the Wilkinson choice for an eigenvalue of T, namely
    // the eigenvalue of the trailing 2x2 block of T which is closer
    // to the last diagonal element of T.
    //
    // Trailing 2x2 block =  (Use i = N-2, j = N-1)
    // [ a  b ] = [ |Di|^2 + |Ei-1|^2      Di Ei      ]
    // [ b* c ]   [     Di* Ei*       |Dj|^2 + |Ei|^2 ]
    // 
    // mu = c - d +- sqrt(d^2+|b|^2), where d = (c-a)/2
    // if d>0 we use +, if d<0 we use -.
    // 
    // For stability when |b| is small, we rearrange this to:
    // mu = c + |b|^2/(d +- sqrt(d^2+|b|^2))
    //    = c +- |b|^2/|d|(1 + sqrt(1+|b|^2/d^2))
    TMVAssert(IsReal(T()));
    const size_t N = D.size();
    TMVAssert(E.size() == N-1);
    TMVAssert(N > 1);

    T a = D(N-2);
    T c = D(N-1);
    T b = E(N-2);
    T bsq = b*b;
    T d = (c-a)/2;
    T dsq = d*d;
    if (dsq < bsq) {
      RealType(T) x = SQRT(dsq+bsq);
      if (d > 0) return c - d + x;
      else return c - d - x;
    } else {
      RealType(T) x = bsq/abs(d)/(1 + SQRT(1+bsq/dsq));
      if (d > 0) return c + x;
      else return c - x;
    }
  }

  template <class T> void ReduceHermTridiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E)
  {
    // Reduce the superdiagonal elements of unreduced HermTridiagonal Matrix T 
    // (given by D,E) while maintaining U B Ut. 
    // Note: the input T must be unreduced - ie. all entries are non-zero.
    const size_t N = D.size();
    TMVAssert(N>0);
    TMVAssert(E.size() == N-1);
    if (U) { TMVAssert(U->rowsize() == N); }
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);

#ifdef XDEBUG
    Matrix<RealType(T)> TT(N,N,RealType(T)(0));
    TT.diag() = D;
    TT.diag(1) = TT.diag(-1) = E;
    const size_t M = U ? U->colsize() : 0;
    Matrix<T> A0(M,M);
    if (U) {
      A0 = *U * TT * U->Adjoint();
    }
#endif

    if (N == 1) return;

    // The reduction is based on the QR algorithm to diagonalize the
    // unreduced symmetric tridiagonal matrix T.
    // The basic idea is as follows:
    // (see Golub and van Loan, chapter 8 for a derivation)
    //
    // if T is a symmetric tridiagonal matrix
    // and mu is (approximately) an eigenvalue of T
    // and the QR decomposition of T - mu I = Q R
    // Then, T' = R Q + mu I will be tridiagonal with the last 
    // subdiagonal element small.
    //
    // Note: T' = Qt (T-muI) Q + muI = Qt T Q.
    // There is a theorem then Q is unique if its first column is 
    // specified.  The first column is given by the Givens step for
    // D(0),E(0).  This step yields (for N = 6):
    //
    // G1 T = [ x x + 0 0 0 ]
    //        [ 0 x x 0 0 0 ]
    //        [ 0 x x x 0 0 ]
    //        [ 0 0 x x x 0 ]
    //        [ 0 0 0 x x x ]
    //        [ 0 0 0 0 x x ]
    //
    // G1 T G1t = [ x x + 0 0 0 ]
    //            [ x x x 0 0 0 ]
    //            [ + x x x 0 0 ]
    //            [ 0 0 x x x 0 ]
    //            [ 0 0 0 x x x ]
    //            [ 0 0 0 0 x x ]
    //
    // The + indicates the element which screws up the tri-diagonality of T'.
    // The rest of the Givens matrices will be the ones that chase this 
    // element down the matrix and away.
    //
    // Wilkinson (1968) suggested that a good choice for mu is
    // the eigenvalue of the trailing 2x2 block of T that is 
    // closer to the trailing diagonal element.
    //
    // At the end of this procedure, E(N-1) should be smaller than it was.
    // Note: This procedure works exactly if N=2.

    RealType(T)* Di = D.ptr();
    RealType(T)* Ei = E.ptr();

    RealType(T) mu = HermTridiagonalTrailingEigenValue(D,E);
    RealType(T) y = *Di - mu;  // = T00 - mu
    RealType(T) x = *Ei;       // = T10
    Givens<RealType(T)> G = Givens_Rotate(y,x);
    for(size_t i=0;;++i,++Di,++Ei) {
      G.SymMult(*Di,*(Di+1),*Ei);
      //RealType(T) Ei2 = *Ei; // = T01
      //G.Mult(*Di,*Ei);
      //G.Mult(Ei2,*(Di+1));
      //G.Mult(*Di,Ei2);
      //G.Mult(*Ei,*(++Di));
      if (U) G.Mult(U->ColPair(i,i+1).Transpose());
      if (i==N-2) break;
      TMVAssert(x==RealType(T)(0));
      G.Mult(x,*(Ei+1));
      G = Givens_Rotate(*Ei,x);
    }
#ifdef XDEBUG
    if (U) {
      TT.diag() = D;
      TT.diag(-1) = TT.diag(1) = E;
      Matrix<T> A2 = *U * TT * U->Adjoint();
      if (Norm(A2-A0) > 0.001*max(RealType(T)(1),Norm(A0))) {
	cerr<<"Decompose from Tridiagonal:\n";
	cerr<<"A0 = "<<A0<<endl;
	cerr<<"A2 = "<<A2<<endl;
	cerr<<"U = "<<*U<<endl;
	cerr<<"TT = "<<TT<<endl;
	abort();
      }
    } 
#endif
  }

  // MJ: Write divide and conquer version
  template <class T> void NonLapHermSV_Decompose_From_Tridiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E)
  {
    TMVAssert(D.size() == E.size()+1);
    if (U) {
      TMVAssert(U->IsSquare());
      TMVAssert(U->rowsize() == D.size());
    }

    const size_t N = D.size();

    // We successively reduce the offdiagonals of T (E) to 0
    // using a sequence of Givens rotations. 
    // The reduction procedure tends to push the values up and left, so it 
    // makes sense to start at the lower right and work back up the matrix.
    // We also set to zero any very small values based on machine precision.
    // Loop invariant: all E(i) with i>=q are 0.
    // Initially q = N-1. (ie. All E(i) are potentially non-zero.)
    // When q = 0, we are done.
    HermTridiagonalChopSmallElements(D,E);
    for(size_t q = N-1; q>0; ) {
      if (E(q-1) == T(0)) --q;
      else {
	size_t p=q-1;
	while (p > 0 && (E(p-1) != T(0))) --p; 
	// Set p such that E(p-1) = 0 and all E(i) with p<=i<q are non-zero.
	if (U) {
	  MatrixView<T> U1 = U->SubMatrix(0,N,p,q+1);
	  ReduceHermTridiagonal(&U1,D.SubVector(p,q+1),E.SubVector(p,q));
	} else {
	  ReduceHermTridiagonal(ZMV<T>(),D.SubVector(p,q+1),E.SubVector(p,q));
	}
	HermTridiagonalChopSmallElements(D,E);
      }
    }
  }

#ifdef LAP 
  template <class T> void LapHermSV_Decompose_From_Tridiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E)
  { NonLapHermSV_Decompose_From_Tridiagonal(U,D,E); }
  template <> void LapHermSV_Decompose_From_Tridiagonal(
      const MatrixView<double>* U, const VectorView<double>& D, 
      const VectorView<double>& E)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->iscm());
      TMVAssert(U->ct()==NonConj);
    }
    int n = D.size();
    int info;
    if (U) {
      if (U->iscm()) {
	char c = 'V';
	int ldu = U->stepj();
	int lgn;
	frexp(double(n),&lgn);
	int lwork = 1+n*(3+2*lgn+3*n);
	double* work = LAP_DWork(lwork);
	int liwork = 6+n*(6+5*lgn);
	int* iwork = LAP_IWork(lwork);
	dstedc(&c,&n,D.ptr(),E.ptr(),U->ptr(),&ldu,work,&lwork,iwork,&liwork,
	    &info);
      } else {
	char c = 'I';
	Matrix<double,ColMajor> U1(n,n);
	U1.SetToIdentity();
	int ldu = U1.stepj();
	int lwork = 1+n*(4+n);
	double* work = LAP_DWork(lwork);
	int liwork = 3+5*n;
	int* iwork = LAP_IWork(lwork);
	dstedc(&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,work,&lwork,iwork,&liwork,
	    &info);
	*U *= U1;
      }
    } else {
      char c = 'N';
      int ldu = n;
      int lwork = 1;
      double* work = LAP_DWork(lwork);
      int liwork = 1;
      int* iwork = LAP_IWork(lwork);
      dstedc(&c,&n,D.ptr(),E.ptr(),0,&ldu,work,&lwork,iwork,&liwork,&info);
    }
    LAP_Results(info,"dstedc");
  }
  template <> void LapHermSV_Decompose_From_Tridiagonal(
      const MatrixView<complex<double> >* U, const VectorView<double>& D, 
      const VectorView<double>& E)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->iscm());
      TMVAssert(U->ct()==NonConj);
    }
    int n = D.size();
    int info;
    if (U) {
      char c = 'I';
      Matrix<double,ColMajor> U1(n,n);
      U1.SetToIdentity();
      int ldu = U1.stepj();
      int lwork = 1+n*(4+n);
      double* work = LAP_DWork(lwork);
      int liwork = 3+5*n;
      int* iwork = LAP_IWork(lwork);
      dstedc(&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,work,&lwork,iwork,&liwork,
	  &info);
      *U *= U1;
    } else {
      char c = 'N';
      int ldu = n;
      int lwork = 1;
      double* work = LAP_DWork(lwork);
      int liwork = 1;
      int* iwork = LAP_IWork(lwork);
      dstedc(&c,&n,D.ptr(),E.ptr(),0,&ldu,work,&lwork,iwork,&liwork,&info);
    }
    LAP_Results(info,"dstedc");
  }
#ifndef NOFLOAT
  template <> void LapHermSV_Decompose_From_Tridiagonal(
      const MatrixView<float>* U, const VectorView<float>& D, 
      const VectorView<float>& E)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->iscm());
      TMVAssert(U->ct()==NonConj);
    }
    int n = D.size();
    int info;
    if (U) {
      if (U->iscm()) {
	char c = 'V';
	int ldu = U->stepj();
	int lgn;
	frexp(double(n),&lgn);
	int lwork = 1+n*(3+2*lgn+3*n);
	float* work = LAP_SWork(lwork);
	int liwork = 6+n*(6+5*lgn);
	int* iwork = LAP_IWork(lwork);
	sstedc(&c,&n,D.ptr(),E.ptr(),U->ptr(),&ldu,work,&lwork,iwork,&liwork,
	    &info);
      } else {
	char c = 'I';
	Matrix<float,ColMajor> U1(n,n);
	U1.SetToIdentity();
	int ldu = U1.stepj();
	int lwork = 1+n*(4+n);
	float* work = LAP_SWork(lwork);
	int liwork = 3+5*n;
	int* iwork = LAP_IWork(lwork);
	sstedc(&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,work,&lwork,iwork,&liwork,
	    &info);
	*U *= U1;
      }
    } else {
      char c = 'N';
      int ldu = n;
      int lwork = 1;
      float* work = LAP_SWork(lwork);
      int liwork = 1;
      int* iwork = LAP_IWork(lwork);
      sstedc(&c,&n,D.ptr(),E.ptr(),0,&ldu,work,&lwork,iwork,&liwork,&info);
    }
    LAP_Results(info,"sstedc");
  }
  template <> void LapHermSV_Decompose_From_Tridiagonal(
      const MatrixView<complex<float> >* U, const VectorView<float>& D, 
      const VectorView<float>& E)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->iscm());
      TMVAssert(U->ct()==NonConj);
    }
    int n = D.size();
    int info;
    if (U) {
      char c = 'I';
      Matrix<float,ColMajor> U1(n,n);
      U1.SetToIdentity();
      int ldu = U1.stepj();
      int lwork = 1+n*(4+n);
      float* work = LAP_SWork(lwork);
      int liwork = 3+5*n;
      int* iwork = LAP_IWork(lwork);
      sstedc(&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,work,&lwork,iwork,&liwork,
	  &info);
      *U *= U1;
    } else {
      char c = 'N';
      int ldu = n;
      int lwork = 1;
      float* work = LAP_SWork(lwork);
      int liwork = 1;
      int* iwork = LAP_IWork(lwork);
      sstedc(&c,&n,D.ptr(),E.ptr(),0,&ldu,work,&lwork,iwork,&liwork,&info);
    }
    LAP_Results(info,"sstedc");
  }
#endif
#endif // LAP

  template <class T> void HermSV_Decompose_From_Tridiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() == U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->ct()==NonConj);
    }

    if (D.size() > 0) {
#ifdef XDEBUG
      Matrix<T> A0(D.size(),D.size());
      Vector<RealType(T)> D0 = D;
      Vector<RealType(T)> E0 = E;
      if (U) {
	Matrix<T> EDE(D.size(),D.size(),T(0));
	EDE.diag(-1) = E;
	EDE.diag(0) = D;
	EDE.diag(1) = E;
	A0 = *U * EDE * U->Adjoint();
      }
#ifdef LAP
      Vector<RealType(T)> D2 = D;
      Vector<RealType(T)> E2 = E;
      Matrix<T> U2(D.size(),D.size());
      if (U) { 
	U2 = *U;
	MatrixView<T> U2V = U2.View();
	NonLapHermSV_Decompose_From_Tridiagonal(&U2V,D2.View(),E2.View());
      } else {
	NonLapHermSV_Decompose_From_Tridiagonal(ZMV<T>(),D2.View(),E2.View());
      }
#endif // LAP
#endif // XDEBUG

#ifdef LAP
      if ((!U || U->iscm() || U->isrm()) && D.step() == 1 && E.step() == 1) 
	LapHermSV_Decompose_From_Tridiagonal(U,D,E);
      else 
#endif
	NonLapHermSV_Decompose_From_Tridiagonal(U,D,E);

      // Now A = U * D * Ut
      // Technically, singular values should be positive, but we allow them
      // to be negative, since these are the eigenvalues of A - no sense
      // killing that.  Also, to make them positive, we'd have to break the
      // V = Ut relationship.  So just keep that in mind later when we use S.

      // Sort output singular values by absolute value:
      size_t sortp[D.size()];
      D.Sort(sortp,DESCEND,ABS_COMP);
      if (U) U->PermuteCols(sortp);
#ifdef XDEBUG
#ifdef LAP
      size_t sortp2[D.size()];
      D2.Sort(sortp2,DESCEND,ABS_COMP);
      U2.PermuteCols(sortp2);
#endif // LAP
#endif // XDEBUG

#ifdef XDEBUG
#ifdef LAP
      if (Norm(D-D2) > 0.001*max(RealType(T)(1),Norm(D))) {
	cerr<<"HermSV_Decompose_From_Tridiagonal:\n";
	cerr<<"D = "<<D0<<endl;
	cerr<<"E = "<<E0<<endl;
	cerr<<"Done: D = "<<D<<endl;
	cerr<<"NonLap version: D = "<<D2<<endl;
	if (U) {
	  cerr<<"U = "<<*U<<endl;
	  cerr<<"UDUt = "<<(*U*DiagMatrixViewOf(D)*U->Adjoint())<<endl;
	  cerr<<"A0 = "<<A0<<endl;
	  cerr<<"NonLap U = "<<U2<<endl;
	  cerr<<"UDUt = "<<(U2*DiagMatrixViewOf(D)*U->Adjoint())<<endl;
	}
	abort();
      }
#else
      if (U) {
	Matrix<T> UDU = *U * DiagMatrixViewOf(D)*U->Adjoint();
	if (Norm(UDU-A0) > 0.001*max(RealType(T)(1),Norm(A0))) {
	  cerr<<"HermSV_Decompose_From_Tridiagonal:\n";
	  cerr<<"D = "<<D0<<endl;
	  cerr<<"E = "<<E0<<endl;
	  cerr<<"Done: D = "<<D<<endl;
	  cerr<<"U = "<<*U<<endl;
	  cerr<<"D = "<<DiagMatrixViewOf(D)<<endl;
	  cerr<<"UD = "<<*U*DiagMatrixViewOf(D)<<endl;
	  cerr<<"DUt = "<<DiagMatrixViewOf(D)*U->Adjoint()<<endl;
	  cerr<<"UDUt = "<<UDU<<endl;
	  cerr<<"A0 = "<<A0<<endl;
	  abort();
	}
      }
#endif // LAP
#endif // XDEBUG
    }
  }

  template <class T> void HermSV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, bool StoreU)
  {
    // Decompose Hermitian A (input as lower tri of U) into U S Ut
    // where S is a diagonal real matrix, and U is a unitary matrix.
    // U,S are N x N
    const size_t N = U.colsize();
    if (N == 0) return;
    TMVAssert(U.IsSquare());
    TMVAssert(S.size() == N);

    // First we reduce A to tridiagonal form: A = U * T * Ut
    // using a series of Householder transformations.
    // The diagonal of the Tridiagonal Matrix T is stored in D.
    // The subdiagonal is stored in E.
    Vector<RealType(T)> E(N-1);
    Vector<T> Ubeta(N-1);
    T det(0);
    SymTridiagonalize(HermMatrixViewOf(U,Lower),Ubeta.View(),S,E.View(),det);
    // Now U stores Householder vectors for U in lower diagonal columns.

    if (StoreU) {
      for(size_t j=N-1;j>0;--j) U.col(j,j,N) = U.col(j-1,j,N);
      U.col(0).Zero();
      U.row(0).Zero();
      U(0,0) = T(1);
      GetQFromQR(U.SubMatrix(1,N,1,N),Ubeta);
    }

    if (StoreU) HermSV_Decompose_From_Tridiagonal(&U,S,E.View());
    else HermSV_Decompose_From_Tridiagonal(ZMV<T>(),S,E.View());
  }

  template <class T> T HermSVDiv<T>::Det() const
  {
    if (!calcdet) {
      det = DiagMatrixViewOf(S).Det();
      calcdet = true;
    }
    return det; 
  }

  template <class T> void HermSVDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  { 
    TMVAssert(U);
    // A^-1 = U S^-1 Ut
    Matrix<T,ColMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    minv = U->Cols(0,kmax) * SinvUt;
  }

  template <class T> void HermSVDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(U);
    // A = U S Ut
    // At = U S Ut
    // AtA = U S^2 Ut
    // (AtA)^-1 = U S^-2 Ut
    //
    Matrix<T,RowMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    minv = SinvUt.Adjoint() * SinvUt;
  }

  template <class T> void HermSVDiv<T>::Thresh(RealType(T) toler,
      ostream* debugout) const
  {
    TMVAssert(toler < RealType(T)(1));
    RealType(T) thresh = abs(S(0))*toler;
    for(kmax=S.size(); kmax>0 && abs(S(kmax-1))<=thresh; --kmax);
    if(debugout) {
      (*debugout)<<"S = "<<S<<endl;
      (*debugout)<<"Smax = "<<abs(S(0))<<", thresh = "<<thresh<<endl;
      (*debugout)<<"kmax = "<<kmax<<" (S.size = "<<S.size()<<")"<<endl;
    }
  }

  template <class T> void HermSVDiv<T>::Top(size_t neigen,
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

  template <class T> void SymSV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>* V, T& det)
  {
    // Decompose complex symmetric A (input as lower tri of U) into U S V
    // where S is a diagonal real matrix, and U,V are unitary matrices.
    // U,S,V are N x N
    // If V = 0, then U,V are not formed.  Only S,det are accurate on return.
    // The determinant is returned in det.
    // (Technically, det is multiplied by the determinant, so det should
    // be set to 1 on entry.)
    const size_t N = U.colsize();
    if (N == 0) return;
    TMVAssert(IsComplex(T()));
    TMVAssert(U.IsSquare());
    TMVAssert(S.size() == N);

    // First we reduce A to tridiagonal form: A = U * T * UT
    // using a series of Householder transformations.
    // The diagonal of the Tridiagonal Matrix T is stored in D.
    // The subdiagonal is stored in E.
    Vector<T> D(N);
    Vector<RealType(T)> E(N-1);
    Vector<T> Ubeta(N-1);
    SymTridiagonalize(SymMatrixViewOf(U,Lower),Ubeta.View(),D.View(),E.View(),
	det);
    // Now U stores Householder vectors for U in lower diagonal columns.

    BandMatrix<T,ColMajor> B(N,N,1,1);
    B.diag() = D;
    B.diag(-1) = E;
    B.diag(1) = E;

    if (V) {
      TMVAssert(V->IsSquare());
      TMVAssert(V->colsize() == N);
      for(size_t j=N-1;j>0;--j) U.col(j,j,N) = U.col(j-1,j,N);
      U.col(0).Zero();
      U.row(0).Zero();
      U(0,0) = T(1);
      GetQFromQR(U.SubMatrix(1,N,1,N),Ubeta);
      *V = U.Transpose();
      Matrix<T,ColMajor> U1(N,N);
      MatrixView<T> U1v = U1.View();
      Matrix<T,ColMajor> V1(N,N);
      MatrixView<T> V1v = V1.View();
      BandSV_Decompose(B,&U1v,S,&V1v,det);
      U = U*U1;
      *V = V1*(*V);
    } else {
      BandSV_Decompose(B,ZMV<T>(),S,ZMV<T>(),det);
    }
  }

  template <class T> void SymSVDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  { 
    TMVAssert(U&&V);
    // A = U S V
    // A^-1 = Vt S^-1 Ut
    Matrix<T,ColMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    minv = V->Adjoint().Cols(0,kmax) * SinvUt;
  }

  template <class T> void SymSVDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(U&&V);
    // A = U S V
    // At = Vt S Ut
    // AtA = Vt S^2 V
    // (AtA)^-1 = Vt S^-2 V
    //
    Matrix<T,RowMajor> SinvV = V->Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    minv = SinvV.Adjoint() * SinvV;
  }

  template <class T> void SymSVDiv<T>::Thresh(RealType(T) toler,
      ostream* debugout) const
  {
    TMVAssert(toler < RealType(T)(1));
    RealType(T) thresh = S(0)*toler;
    for(kmax=S.size(); kmax>0 && abs(S(kmax-1))<=thresh; --kmax);
    if(debugout) {
      (*debugout)<<"S = "<<S<<endl;
      (*debugout)<<"Smax = "<<S(0)<<", thresh = "<<thresh<<endl;
      (*debugout)<<"kmax = "<<kmax<<" (S.size = "<<S.size()<<")"<<endl;
    }
  }

  template <class T> void SymSVDiv<T>::Top(size_t neigen,
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

#define InstFile "TMV_SymSVDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


