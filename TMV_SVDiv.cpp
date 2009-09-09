
#include "TMV.h"
#include "TMV_Householder.h"
#include "TMV_Givens.h"
#include "TMV_Diag.h"
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

  bool RecursiveSV = true;

#ifdef TMV_BLOCKSIZE
  const size_t BIDIAG_BLOCKSIZE = TMV_BLOCKSIZE;
  const size_t DC_LIMIT = TMV_BLOCKSIZE;
#else
  const size_t BIDIAG_BLOCKSIZE = 16;
  const size_t DC_LIMIT = 32;
#endif

  template <class T> const MatrixView<T>* ZMV() 
  { return (const MatrixView<T>*)(0); }

  //
  // Start with class methods 
  // (These mostly just call other functions which are below in this file.)
  //
  
#define APTR inplace ? A.NonConst().ptr() : new T[A.colsize()*A.rowsize()]     
#define UX istrans ? \
  inplace ? A.NonConst().Transpose() : \
  MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor) : \
  inplace ? A.NonConst().View() : \
  MatrixViewOf(Aptr,A.colsize(),A.rowsize(), \
      A.IsSquare() ? BaseStorOf(A) : ColMajor)

  template <class T> SVFDiv<T>::SVFDiv(const GenMatrix<T>& A,
      bool _inplace, bool StoreU, bool StoreV) :
    istrans(A.colsize() < A.rowsize()), inplace(_inplace), Aptr(APTR),
    U(new MatrixView<T>(UX)), S(U->rowsize()), V(0), det(T(1))
  {
    if (istrans) {
      if (inplace) { TMVAssert(A.Transpose() == *U); }
      else *U = A.Transpose();
    } else {
      if (inplace) { TMVAssert(A == *U); }
      else *U = A;
    }

    if (istrans) swap(StoreU,StoreV);

    if (StoreV) {
      V = new Matrix<T,ColMajor>(U->rowsize(),U->rowsize());
      SV_Decompose(*U,S.View(),V->View(),det,StoreU);
    }
    else SV_Decompose(*U,S.View(),det,StoreU);

    if (!StoreU && !inplace) { delete [] Aptr; Aptr = 0; }
    if (!StoreU) { delete U; U = 0; }

    // Set kmax for actual 0 elements (to within machine precision).
    // Any further cut in the number of singular values to use
    // should be done by the user.
    Thresh(Epsilon<T>());
  }

  // Normally A = U U1 S V1 V
  // where U,V are stored in UV, Ubeta, Vbeta
  //
  // If UV2 and Qbeta are set (ie. pointers are not 0), then
  // A = Q U U1 S V1 V
  // where U,V are stored in *UV2, Ubeta, Vbeta
  // and Q is stored in UV, *Qbeta
  //
  // If istrans = true, then all of this is the storage for A.Transpose().
  //
  template <class T> SVDiv<T>::SVDiv(const GenMatrix<T>& A, bool _inplace) :
    istrans(A.colsize() < A.rowsize()), inplace(_inplace),
    Aptr(APTR), UV(UX),
    Ubeta(UV.rowsize()), Vbeta(UV.rowsize()-1),
    UV2(0), Qbeta(0), U1(UV.rowsize(),UV.rowsize()),
    V1(UV.rowsize(),UV.rowsize()), S(UV.rowsize()), det(T(1))
  {
    if (istrans) {
      if (inplace) { TMVAssert(A.Transpose() == UV); }
      else UV = A.Transpose();
    } else {
      if (inplace) { TMVAssert(A == UV); }
      else UV = A;
    }

    SV_Decompose(UV,Ubeta.View(),Vbeta.View(),UV2,Qbeta,
	U1.View(),V1.View(),S.View(),det);

    // Set kmax for actual 0 elements (to within machine precision).
    // Any further cut in the number of singular values to use
    // should be done by the user.
    Thresh(Epsilon<T>());
  }
#undef UX
#undef APTR

  template <class T> SVFDiv<T>::~SVFDiv() 
  { 
    if (!inplace && Aptr) delete [] Aptr; 
    if (U) delete U;
    if (V) delete V;
  }

  template <class T> SVDiv<T>::~SVDiv() 
  { if (!inplace) delete [] Aptr; }

  template <class T> bool SVDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenMatrix<T>* mm = dynamic_cast<const GenMatrix<T>*>(&m);
    TMVAssert(mm);
    if (fout) {
      *fout << "M = "<<(istrans ? mm->Transpose() : mm->View())<<endl;
      *fout << "U = "<<GetU()<<endl;
      *fout << "S = "<<GetS()<<endl;
      *fout << "V = "<<GetV()<<endl;
    }
    Matrix<T> m2 = GetU()*DiagMatrixViewOf(GetS())*GetV();
    RealType(T) nm = Norm(m2-(istrans ? mm->Transpose() : mm->View()));
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    if (fout) {
      *fout << "USV = "<<m2<<endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm<<"  "<<S.size()*Epsilon<T>()<<endl;
    }
    return nm < S.size()*Epsilon<T>();
  }

  template <class T> bool SVFDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenMatrix<T>* mm = dynamic_cast<const GenMatrix<T>*>(&m);
    TMVAssert(mm);
    if (fout) {
      *fout << "M = "<<*mm<<endl;
      *fout << "U = "<<GetU()<<endl;
      *fout << "S = "<<GetS()<<endl;
      *fout << "V = "<<GetV()<<endl;
    }
    Matrix<T> m2 = GetU()*DiagMatrixViewOf(GetS())*GetV();
    RealType(T) nm = Norm(m2-*mm);
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    if (fout) {
      *fout << "USV = "<<m2<<endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm<<"  "<<S.size()*Epsilon<T>()<<endl;
    }
    return nm < S.size()*Epsilon<T>();
  }

  template <class T> void SVFDiv<T>::Inverse(const MatrixView<T>& minv) const
  { 
    TMVAssert(U && V);
    if (istrans) {
      // A^-1 = (Vt S^-1 Ut)T = U* S^-1 V*
      Matrix<T,ColMajor> SinvV = V->Conjugate().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      minv = U->Conjugate().Cols(0,kmax) * SinvV;
    } else {
      // A^-1 = Vt S^-1 Ut
      Matrix<T,ColMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      minv = V->Adjoint().Cols(0,kmax) * SinvUt;
    }
  }

  template <class T> void SVDiv<T>::Inverse(const MatrixView<T>& minv) const
  { 
    if (istrans) {
      TMVAssert(minv.colsize() == UV.colsize());
      TMVAssert(minv.rowsize() == UV.rowsize());
      DoLDiv(Eye<T,ColMajor>(UV.rowsize()),minv);
    } else {
      TMVAssert(minv.colsize() == UV.rowsize());
      TMVAssert(minv.rowsize() == UV.colsize());
      DoRDiv(Eye<T,ColMajor>(UV.rowsize()),minv);
    }
  }

  template <class T> void SVFDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(V);
    // A = U S V
    // At = Vt S Ut
    // AtA = Vt S^2 V
    // (AtA)^-1 = Vt S^-2 V
    //
    // if istrans:
    // AT = U S V
    // At = U* S V*
    // AAt = VT S^2 V*
    // (AAt)^-1 = VT S^-2 V*
    //
    Matrix<T,ColMajor> SinvV = V->Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    if (istrans)
      minv = SinvV.Transpose() * SinvV.Conjugate();
    else
      minv = SinvV.Adjoint() * SinvV;
  }

  template <class T> void SVDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    // A = U U1 S V1 V
    // At = Vt V1t S U1t Ut
    // AtA = Vt V1t S^2 V1 V
    // (AtA)^-1 = Vt V1t S^-2 V1 V
    //
    // if istrans:
    // AT = U S V
    // At = U* S V*
    // AAt = VT S^2 V*
    // (AAt)^-1 = VT S^-2 V*
    
    const size_t N = UV.rowsize();

    Matrix<T,ColMajor> VtSinv = V1.Rows(0,kmax).Adjoint() %
      DiagMatrixViewOf(S.SubVector(0,kmax));
    if (UV2) {
      Q_RDivEq(UV2->SubMatrix(0,N-1,1,N).Transpose(),Vbeta,
	  VtSinv.Rows(1,N).Transpose());
    } else {
      Q_RDivEq(UV.SubMatrix(0,N-1,1,N).Transpose(),Vbeta,
	  VtSinv.Rows(1,N).Transpose());
    }
    if (istrans)
      minv = VtSinv.Conjugate() * VtSinv.Transpose();
    else
      minv = VtSinv * VtSinv.Adjoint();
  }

  template <class T> void SVFDiv<T>::Thresh(RealType(T) toler,
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

  template <class T> void SVDiv<T>::Thresh(RealType(T) toler,
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

  template <class T> void SVFDiv<T>::Top(size_t neigen,
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

  template <class T> void SVDiv<T>::Top(size_t neigen,
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

  template <class T> Matrix<T,ColMajor> SVDiv<T>::GetU() const
  {
    if (UV2) {
      TMVAssert(Qbeta);
      Matrix<T,ColMajor> U=UV;
      GetQFromQR(U.View(),*Qbeta);
      Q_LDivEq(*UV2,Ubeta,U.Adjoint());
      return U*U1;
    } else {
      Matrix<T,ColMajor> U=UV;
      GetQFromQR(U.View(),Ubeta);
      return U*U1;
    }
  }

  template <class T> Matrix<T,ColMajor> SVDiv<T>::GetV() const
  {
    const size_t N = UV.rowsize();
    Matrix<T,ColMajor> V(N,N);
    V.row(0).MakeBasis(0);
    if (UV2) V.Rows(1,N) = UV2->Rows(0,N-1);
    else V.Rows(1,N) = UV.Rows(0,N-1);
    V.col(0,1,N).Zero();
    GetQFromQR(V.SubMatrix(1,N,1,N).Transpose(),Vbeta);
    return V1*V;
  }

  //
  // LDiv
  //

  template <class T1, class T2, class T3> void SV_LDiv(
      const GenMatrix<T1>& U, const GenVector<RealType(T1)>& S, 
      const GenMatrix<T1>& V, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T3>& x)
  {
    // A x = m
    // U S V x = m
    // x = Vt S^-1 Ut m
    TMVAssert(m.colsize() == U.colsize()); // = M
    TMVAssert(x.colsize() == V.rowsize()); // = N
    TMVAssert(x.rowsize() == m.rowsize()); // = R
    TMVAssert(kmax <= V.rowsize()); // = K
    TMVAssert(kmax <= U.colsize());
    Matrix<T3> m2 = U.Adjoint().Rows(0,kmax) * m; // KxR
    m2 /= DiagMatrixViewOf(S.SubVector(0,kmax));
    x = V.Adjoint().Cols(0,kmax) * m2; // NxR
  }

  template <class T1, class T2, class T3> void SV_LDiv(
      const GenMatrix<T1>& UV,
      const GenVector<T1>& Ubeta, const GenVector<T1>& Vbeta, 
      const GenMatrix<T1>* UV2, const GenVector<T1>* Qbeta,
      const GenMatrix<RealType(T1)>& U1, const GenMatrix<RealType(T1)>& V1,
      const GenVector<RealType(T1)>& S, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T3>& x)
  {
    // A x = m
    // U U1 S V1 V x = m
    // x = Vt V1t S^-1 U1t Ut m
    TMVAssert(m.colsize() == UV.colsize()); // = M
    TMVAssert(x.colsize() == UV.rowsize()); // = N
    TMVAssert(x.rowsize() == m.rowsize()); // = R
    TMVAssert(kmax <= UV.rowsize());
    TMVAssert(UV.rowsize() <= UV.colsize());
    
    const size_t N = UV.rowsize();

    Matrix<T3,ColMajor> m2 = m; // MxR
    if (UV2) {
      TMVAssert(Qbeta);
      Q_LDivEq(UV,*Qbeta,m2.View());
      Q_LDivEq(*UV2,Ubeta,m2.Rows(0,N));
    } else {
      Q_LDivEq(UV,Ubeta,m2.View());
    }
    Matrix<T3,ColMajor> m3 = U1.Adjoint().Rows(0,kmax) * m2.Rows(0,N); 
      // = KxR
    m3 /= DiagMatrixViewOf(S.SubVector(0,kmax));
    x = V1.Adjoint().Cols(0,kmax) * m3; // NxR
    if (UV2) {
      Q_RDivEq(UV2->SubMatrix(0,N-1,1,N).Transpose(),Vbeta,
	  x.Rows(1,N).Transpose());
    } else {
      Q_RDivEq(UV.SubMatrix(0,N-1,1,N).Transpose(),Vbeta,
	  x.Rows(1,N).Transpose());
    }
  }

  //
  // RDiv
  //

  template <class T1, class T2, class T3> void SV_RDiv(
      const GenMatrix<T1>& U, const GenVector<RealType(T1)>& S, 
      const GenMatrix<T1>& V, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T3>& x) 
  {
    // x A = m
    // x U S V = m
    // x = m Vt S^-1 Ut
    TMVAssert(m.rowsize() == V.rowsize()); // = N
    TMVAssert(x.rowsize() == U.colsize()); // = M
    TMVAssert(x.colsize() == m.colsize()); // = R
    TMVAssert(kmax <= U.colsize()); // = K
    TMVAssert(kmax <= V.rowsize());

    Matrix<T3,ColMajor> m2 = m * V.Adjoint().Cols(0,kmax); // = RxK
    m2 %= DiagMatrixViewOf(S.SubVector(0,kmax));
    x = m2 * U.Adjoint().Rows(0,kmax); // = RxM
  }

  template <class T1, class T2, class T3> void SV_RDiv(
      const GenMatrix<T1>& UV,
      const GenVector<T1>& Ubeta, const GenVector<T1>& Vbeta, 
      const GenMatrix<T1>* UV2, const GenVector<T1>* Qbeta,
      const GenMatrix<RealType(T1)>& U1, const GenMatrix<RealType(T1)>& V1,
      const GenVector<RealType(T1)>& S, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T3>& x)
  {
    // x A = m
    // x U U1 S V1 V = m
    // x = m Vt V1t S^-1 U1t Ut
    TMVAssert(m.rowsize() == UV.rowsize()); // = N
    TMVAssert(x.rowsize() == UV.colsize()); // = M
    TMVAssert(x.colsize() == m.colsize()); // = R
    TMVAssert(kmax <= UV.rowsize()); // = K
    TMVAssert(UV.rowsize() <= UV.colsize());
    
    const size_t N = UV.rowsize();

    Matrix<T3,ColMajor> m2 = m; // = RxN
    if (UV2) {
      Q_LDivEq(UV2->SubMatrix(0,N-1,1,N).Transpose(),Vbeta,
	  m2.Cols(1,N).Transpose());
    } else {
      Q_LDivEq(UV.SubMatrix(0,N-1,1,N).Transpose(),Vbeta,
	  m2.Cols(1,N).Transpose());
    }
    Matrix<T3,ColMajor> m3 = m2 * V1.Adjoint().Cols(0,kmax); // = RxK
    m3 %= DiagMatrixViewOf(S.SubVector(0,kmax));
    x.Cols(0,N) = m3 * U1.Adjoint().Rows(0,kmax); // = RxN
    x.Cols(N,x.rowsize()).Zero(); // x = RxM
    if (UV2) {
      TMVAssert(Qbeta);
      Q_RDivEq(*UV2,Ubeta,x.Cols(0,N));
      Q_RDivEq(UV,*Qbeta,x);
    } else {
      Q_RDivEq(UV,Ubeta,x);
    }
  }

  //
  // Bidiagonalize
  //
  
  template <class T> void NonBlockBidiagonalize(
      const MatrixView<T>& A, const VectorView<T>& Ubeta,
      const VectorView<T>& Vbeta, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, T& det)
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif
    // Decompose A into U B V
    // The Bidiagonal Matrix B is stored as two vectors: D, E
    // D is the diagonal, E is the super-diagonal
    // A along with Ubeta and Vbeta hold the U and V matrices.
    const size_t M = A.colsize();
    const size_t N = A.rowsize();

    TMVAssert(N <= M);
    TMVAssert(N > 0);
    TMVAssert(Ubeta.size() == N);
    TMVAssert(Vbeta.size() == N-1);
    TMVAssert(D.size() == N);
    TMVAssert(E.size() == N-1);

    // We use Householder reflections to reduce A to the bidiagonal form:
    for(size_t j=0;j<N-1;++j) {
      Ubeta(j) = Householder_Reflect(A.SubMatrix(j,M,j,N),det);
      Vbeta(j) = Householder_Reflect(A.Transpose().SubMatrix(j+1,N,j,M),det);
    }
    Ubeta(N-1) = Householder_Reflect(A.col(N-1,N-1,M),det);

    // The bidiagonal of A is the bidiagonal we want, so copy it to D,E
    if (IsComplex(T())) {
      TMVAssert(NormInf(A.diag().Imag()) == RealType(T)(0));
      TMVAssert(NormInf(A.diag(1).Imag()) == RealType(T)(0));
    }
    D = A.diag().Real();
    E = A.diag(1).Real();

#ifdef XDEBUG
    Matrix<T> U = A;
    GetQFromQR(U.View(),Ubeta);
    Matrix<T> V(N,N);
    V.SetToIdentity();
    V.SubMatrix(1,N,1,N) = A.SubMatrix(0,N-1,1,N);
    GetQFromQR(V.SubMatrix(1,N,1,N).Transpose(),Vbeta);
    Matrix<RealType(T)> B(N,N,RealType(T)(0));
    B.diag() = D;
    B.diag(1) = E;
    Matrix<T> AA = U*B*V;
    if (Norm(A0-AA) > 0.001*Norm(A0)) {
      cerr<<"Bidiagonalize: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"A = "<<A<<endl;
      cerr<<"Ubeta = "<<Ubeta<<endl;
      cerr<<"Vbeta = "<<Vbeta<<endl;
      cerr<<"U = "<<U<<endl;
      cerr<<"B = "<<B<<endl;
      cerr<<"V = "<<V<<endl;
      cerr<<"UBV = "<<AA<<endl;
      abort();
    }
#endif
  }

  template <class T> void BlockBidiagonalize(
      const MatrixView<T>& A, const VectorView<T>& Ubeta,
      const VectorView<T>& Vbeta, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, T& det)
  {
    // Normally we keep the Z matrix for block Householder matrices where 
    // the block Householder is I - YZYt (and Z is upper triangular).
    //
    // However, since the bidiagonalizing process proceeds along both 
    // the rows and columns, we have a problem.  If we just kept the Z
    // matrix for each set, then we would not be able to update the 
    // resulting submatrix at the end of the block.
    // m' = (I-YZYt) m (I-XtWX)
    //
    // For the left multiply, the m needs to have full height for the 
    // Yt m product.  Likewise, for the right multiply, it needs full width.
    // Since we update the first K rows and columns with the Y matrix, 
    // this doesn't work.  So instead of keeping Z,W we are forced to use
    // a bit more temporary storage and store the products ZYtm and mXtW.
    //
    // Furthermore, the m in these products is maintained such that the
    // it already has the appropriate multiplies from the other side.
    // Then, when we are done with the block, the update becomes just:
    //
    // m' = m' - Y (ZYtm) - (mXtW) X
    //
    const size_t M = A.colsize();
    const size_t N = A.rowsize();

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    TMVAssert(N <= M);
    TMVAssert(N > 0);
    TMVAssert(Ubeta.size() == N);
    TMVAssert(Vbeta.size() == N-1);
    TMVAssert(D.size() == N);
    TMVAssert(E.size() == N-1);

    Matrix<T,RowMajor> ZYtm(min(BIDIAG_BLOCKSIZE,N-1),N);
    Matrix<T,ColMajor> mXtW(M,min(BIDIAG_BLOCKSIZE,N-1));
    for(size_t j1=0;j1<N-1;) {
      size_t j2 = min(N-1,j1+BIDIAG_BLOCKSIZE);
      for(size_t j=j1,jmj1=0;j<j2;++j,++jmj1) {

	// Update current column:
	// A(j:M,j) -= Y(j:M,0:j) ZYtm(0:j,j) + mXtW(j:M,0:j) X(0:j,j)
	//
	VectorView<T> u = A.col(j,j,M);
	MatrixView<T> Y0 = A.SubMatrix(j,M,j1,j);
	MatrixView<T> ZYtm0 = ZYtm.SubMatrix(0,jmj1,j,N);
	MatrixView<T> mXtW0 = mXtW.SubMatrix(j,M,0,jmj1);
	MatrixView<T> X0 = A.SubMatrix(j1,j,j,N);
	if (j > j1) {
	  u -= Y0 * ZYtm0.col(0);
	  u -= mXtW0 * X0.col(0);
	}

	// Do the Householder reflection for U
	// Copy the reflection into D(j), and set the top of the 
	// Householder vector to be explicitly 1.  (It makes like easier
	// if it's actually 1 rather than dealing with it implicitly.)
	//
	T bu = Householder_Reflect(u,det);
	Ubeta(j) = bu;
	TMVAssert(IMAG(u(0)) == RealType(T)(0)); 
	D(j) = REAL(u(0));
	u(0) = T(1);

	// Update ZYtm:
	//
	// ZYtm(j,j+1:N) = Z(j,0:j+1) Yt(0:j+1,0:M) m'(0:M,j+1:N)
	// m' is the matrix after taking into account the Householder 
	// multiplies that we have already done:
	//
	// m' = m' - Y ZYtm - mXtW X
	//
	// The new ZYtm(j,j+1:N) = (Z Yt m')(j,j+1:N)
	// = Z(j,0:j+1) Yt(0:j+1,0:M) m'(0:M,j+1:N)
	//
	// Z is Upper Triangular, so Z(j,0:j+1) = 0
	// Also Z(j,j) = bu, so:
	// ZYtm(j,j+1:N) = bu Yt(j,0:M) m'(0:M,j+1:N)
	//
	// Y is Lower Unit Trapezoidal, so Yt(j,0:j) = 0
	// The rest, Yt(j,j:M) is just ut, so:
	// ZYtm(j,j+1:N) = bu ut m'(j:M,j+1:N)
	//
	// Finally, expand out the m':
	// m'(j:M,j+1:N) = A(j:M,j+1:N) 
	//                 - Y(j:M,0:j) ZYtm(0:j,j+1:N)
	//                 - mXtW(j:M,0:j) X(0:j,j+1:N)
	//
	VectorView<T> ZYtmj = ZYtm.row(jmj1,j+1,N);
	VectorView<T> temp = ZYtm.row(jmj1,j1,j);
	VectorView<T> ut = u.Conjugate();
	MatrixView<T> ZYtm0a = ZYtm0.Cols(1,ZYtm0.rowsize());
	MatrixView<T> X0a = X0.Cols(1,X0.rowsize());
	ZYtmj = ut * A.SubMatrix(j,M,j+1,N);
	ZYtmj -= (temp = ut*Y0) * ZYtm0a;
	ZYtmj -= (temp = ut*mXtW0) * X0a;
	ZYtmj *= bu;

	// Update the current row:
	// A(j,j+1:N) -= Y(j,0:j+1) ZYtm(0:j+1,j+1:N) + mXtW(j,0:j) X(0:j:j+1,N)
	//
	MatrixView<T> Y1 = A.SubMatrix(j,M,j1,j+1);
	MatrixView<T> ZYtm1 = ZYtm.SubMatrix(0,jmj1+1,j+1,N);
	VectorView<T> v = A.row(j,j+1,N);
	v -= Y1.row(0) * ZYtm1;
	v -= mXtW0.row(0) * X0a;

	// Do the Householder reflection for V
	//
	T bv = Householder_Reflect(v,det);
	Vbeta(j) = bv;
	TMVAssert(IMAG(v(0)) == RealType(T)(0));
	E(j) = REAL(v(0));
	v(0) = T(1);

	// Update mXtW:
	//
	// mXtW(j+1:M,j) = m'(j+1:M,0:N) Xt(0:N,0:j+1) W(0:j+1,j)
	// = bv m'(j+1:M,j+1:N) vt
	//
	// And m' is:
	//
	// m'(j+1:M,j+1:N) = A(j+1:M,j+1:N) 
	//                   - Y(j+1:M,0:j+1) ZYtm(0:j+1,j+1:N) 
	//                   - mXtW(j+1:M,0:j) X(0:j,j+1:N)
	//
	VectorView<T> mXtWj = mXtW.col(jmj1,j+1,M);
	VectorView<T> temp1 = mXtW.col(jmj1,j1,j+1);
	VectorView<T> temp2 = mXtW.col(jmj1,j1,j);
	VectorView<T> vt = v.Conjugate();
	MatrixView<T> Y1a = Y1.Rows(1,Y1.colsize());
	MatrixView<T> mXtW0a = mXtW0.Rows(1,mXtW0.colsize());
	mXtWj = A.SubMatrix(j+1,M,j+1,N)*vt;
	mXtWj -= Y1a * (temp1 = ZYtm1*vt);
	mXtWj -= mXtW0a * (temp2 = X0a*vt);
	mXtWj *= bv;
      }

      // Update the rest of the matrix:
      A.SubMatrix(j2,M,j2,N) -= A.SubMatrix(j2,M,j1,j2) *
	ZYtm.SubMatrix(0,j2-j1,j2,N);
      A.SubMatrix(j2,M,j2,N) -= mXtW.SubMatrix(j2,M,0,j2-j1) *
	A.SubMatrix(j1,j2,j2,N);
      j1 = j2;
    }

    // Do the last U Householder vector:
    Ubeta(N-1) = Householder_Reflect(A.col(N-1,N-1,M),det);
    TMVAssert(IMAG(A(N-1,N-1)) == RealType(T)(0));
    D(N-1) = REAL(A(N-1,N-1));

#ifdef XDEBUG
    Matrix<T> U = A;
    GetQFromQR(U.View(),Ubeta);
    Matrix<T> V(N,N);
    V.SetToIdentity();
    V.SubMatrix(1,N,1,N) = A.SubMatrix(0,N-1,1,N);
    GetQFromQR(V.SubMatrix(1,N,1,N).Transpose(),Vbeta);
    Matrix<RealType(T)> B(N,N,RealType(T)(0));
    B.diag() = D;
    B.diag(1) = E;
    Matrix<T> AA = U*B*V;
    if (Norm(A0-AA) > 0.001*Norm(A0)) {
      cerr<<"Bidiagonalize: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"A = "<<A<<endl;
      cerr<<"U = "<<U<<endl;
      cerr<<"B = "<<B<<endl;
      cerr<<"V = "<<V<<endl;
      cerr<<"UBV = "<<AA<<endl;
      Matrix<T,ColMajor> A2 = A0;
      Vector<T> Ub2(Ubeta.size());
      Vector<T> Vb2(Vbeta.size());
      Vector<RealType(T)> D2(D.size());
      Vector<RealType(T)> E2(E.size());
      T det2;
      NonBlockBidiagonalize(A2.View(),Ub2.View(),Vb2.View(),D2.View(),
	  E2.View(),det2);
      cerr<<"NonBlock: "<<A2<<endl;
      cerr<<"Ubeta = "<<Ubeta<<endl;
      cerr<<"Nonblock: "<<Ub2<<endl;
      cerr<<"Vbeta = "<<Vbeta<<endl;
      cerr<<"Nonblock: "<<Vb2<<endl;
      cerr<<"D = "<<D<<endl;
      cerr<<"D2 = "<<D2<<endl;
      cerr<<"E = "<<E<<endl;
      cerr<<"E2 = "<<E2<<endl;
      
      abort();
    }
#endif
  }

  template <class T> void NonLapBidiagonalize(
      const MatrixView<T>& A, const VectorView<T>& Ubeta,
      const VectorView<T>& Vbeta, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, T& det)
  {
    TMVAssert(A.rowsize() <= A.colsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(Ubeta.size() == A.rowsize());
    TMVAssert(Vbeta.size() == A.rowsize()-1);
    TMVAssert(D.size() == A.rowsize());
    TMVAssert(E.size() == A.rowsize()-1);

    if (A.rowsize() > BIDIAG_BLOCKSIZE)
      BlockBidiagonalize(A,Ubeta,Vbeta,D,E,det);
    else
      NonBlockBidiagonalize(A,Ubeta,Vbeta,D,E,det);
  }

#ifdef LAP
  template <class T> void LapBidiagonalize(
      const MatrixView<T>& A, const VectorView<T>& Ubeta,
      const VectorView<T>& Vbeta, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, T& det)
  { NonLapBidiagonalize(A,Ubeta,Vbeta,D,E,det); }
  template <> void LapBidiagonalize(
      const MatrixView<double>& A, const VectorView<double>& Ubeta,
      const VectorView<double>& Vbeta, const VectorView<double>& D,
      const VectorView<double>& E, double& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(Ubeta.size() == A.rowsize());
    TMVAssert(Vbeta.size() == A.rowsize()-1);
    TMVAssert(A.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);

    int m = A.colsize();
    int n = A.rowsize();
    int ldu = A.stepj();
    Vector<double> Vbeta2(n);  
    // Stupid LAPack requires an extra element in the Vbeta vector
    // which it sets to 0 (!!!) rather than ignores.
    // So we need to create a temporary Vector which is size n.
    int lwork = (m+n)*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    dgebrd(&m,&n,A.ptr(),&ldu,D.ptr(),E.ptr(),Ubeta.ptr(),Vbeta2.ptr(),
	work,&lwork,&info);
    Vbeta = Vbeta2.SubVector(0,n-1);
    LAP_Results(info,int(work[0]),m,n,lwork,"dgebrd");
    if (det) {
      for(size_t i=0;i<Ubeta.size();++i) if (Ubeta(i) != double(0)) 
	det = -det;
      for(size_t i=0;i<Vbeta.size();++i) if (Vbeta(i) != double(0)) 
	det = -det;
    }
  }
  template <> void LapBidiagonalize(
      const MatrixView<complex<double> >& A, 
      const VectorView<complex<double> >& Ubeta, 
      const VectorView<complex<double> >& Vbeta,
      const VectorView<double>& D, const VectorView<double>& E,
      complex<double>& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(Ubeta.size() == A.rowsize());
    TMVAssert(Vbeta.size() == A.rowsize()-1);
    TMVAssert(A.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int ldu = A.stepj();
    Vector<complex<double> > Vbeta2(n);
    int lwork = (m+n)*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;

    zgebrd(&m,&n,LAP_Complex(A.ptr()),&ldu,D.ptr(),E.ptr(),
	LAP_Complex(Ubeta.ptr()),LAP_Complex(Vbeta2.ptr()),
	LAP_Complex(work),&lwork,&info);
    Ubeta.ConjugateSelf();
    Vbeta = Vbeta2.SubVector(0,n-1);
    LAP_Results(info,int(real(work[0])),m,n,lwork,"zgebrd");
    if (det!=double(0)) {
      for(size_t i=0;i<Ubeta.size();++i) if (Ubeta(i) != double(0)) {
	det *= conj(Ubeta(i)*Ubeta(i))/norm(Ubeta(i));
      }
      for(size_t i=0;i<Vbeta.size();++i) if (Vbeta(i) != double(0)) {
	det *= -conj(Vbeta(i)*Vbeta(i))/norm(Vbeta(i));
      }
    }
  }
#ifndef NOFLOAT
  template <> void LapBidiagonalize(
      const MatrixView<float>& A, const VectorView<float>& Ubeta,
      const VectorView<float>& Vbeta, const VectorView<float>& D,
      const VectorView<float>& E, float& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(Ubeta.size() == A.rowsize());
    TMVAssert(Vbeta.size() == A.rowsize()-1);
    TMVAssert(A.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);

    int m = A.colsize();
    int n = A.rowsize();
    int ldu = A.stepj();
    Vector<float> Vbeta2(n);  
    int lwork = (m+n)*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    sgebrd(&m,&n,A.ptr(),&ldu,D.ptr(),E.ptr(),Ubeta.ptr(),Vbeta2.ptr(),
	work,&lwork,&info);
    Vbeta = Vbeta2.SubVector(0,n-1);
    LAP_Results(info,int(work[0]),m,n,lwork,"dgebrd");
    if (det) {
      for(size_t i=0;i<Ubeta.size();++i) if (Ubeta(i) != float(0)) 
	det = -det;
      for(size_t i=0;i<Vbeta.size();++i) if (Vbeta(i) != float(0)) 
	det = -det;
    }
  }
  template <> void LapBidiagonalize(
      const MatrixView<complex<float> >& A, 
      const VectorView<complex<float> >& Ubeta, 
      const VectorView<complex<float> >& Vbeta,
      const VectorView<float>& D, const VectorView<float>& E,
      complex<float>& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(Ubeta.size() == A.rowsize());
    TMVAssert(Vbeta.size() == A.rowsize()-1);
    TMVAssert(A.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int ldu = A.stepj();
    Vector<complex<float> > Vbeta2(n);
    int lwork = (m+n)*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;

    cgebrd(&m,&n,LAP_Complex(A.ptr()),&ldu,D.ptr(),E.ptr(),
	LAP_Complex(Ubeta.ptr()),LAP_Complex(Vbeta2.ptr()),
	LAP_Complex(work),&lwork,&info);
    Ubeta.ConjugateSelf();
    Vbeta = Vbeta2.SubVector(0,n-1);
    LAP_Results(info,int(real(work[0])),m,n,lwork,"zgebrd");
    if (det!=float(0)) {
      for(size_t i=0;i<Ubeta.size();++i) if (Ubeta(i) != float(0)) {
	det *= conj(Ubeta(i)*Ubeta(i))/norm(Ubeta(i));
      }
      for(size_t i=0;i<Vbeta.size();++i) if (Vbeta(i) != float(0)) {
	det *= -conj(Vbeta(i)*Vbeta(i))/norm(Vbeta(i));
      }
    }
  }
#endif // NOFLOAT
#endif // LAP

  template <class T> void Bidiagonalize(
      const MatrixView<T>& A, const VectorView<T>& Ubeta,
      const VectorView<T>& Vbeta, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == D.size());
    TMVAssert(Ubeta.size() == A.rowsize());
    TMVAssert(Vbeta.size() == A.rowsize()-1);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(Ubeta.step() == 1);
    TMVAssert(Vbeta.step() == 1);
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    
    if (A.rowsize() > 0) {
#ifdef LAP
      if (A.iscm()) 
	LapBidiagonalize(A,Ubeta,Vbeta,D,E,det);
      else 
#endif
	NonLapBidiagonalize(A,Ubeta,Vbeta,D,E,det);
    }
  }

  //
  // Decompose_From_Bidiagonal: QR method
  //
  
  template <class T> void BidiagonalChopSmallElements(
      const VectorView<T>& D, const VectorView<T>& E)
  {
    // This routines sets to 0 any elements in E or D which
    // are essentially 0, given the machine precision:
    // if |E(i)| < Epsilon() * (|D(i)| + |D(i+1)|) then E(i) <- 0
    // if |D(i)| < Epsilon() * |B| then D(i) <- 0
    TMVAssert(E.size() == D.size()-1);
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    RealType(T) eps = Epsilon<T>();
    RealType(T) dthresh = eps * SQRT(NormSq(D) + NormSq(E));

    VIt<T,Unit,NonConj> Di = D.begin();
    VIt<T,Unit,NonConj> Ei = E.begin();
    for(size_t k=E.size();k>0;--k,++Di,++Ei) {
      if (abs(*Di) < dthresh) *Di = T(0);
      if (abs(*Ei) < eps * (abs(*Di)+abs(*(Di+1)))) *Ei = T(0);
    }
    if (abs(*Di) < dthresh) *Di = T(0);
  }

  template <class T, class TU> void BidiagonalZeroFirstRow(
      const MatrixView<TU>* U, const VectorView<T>& D, 
      const VectorView<T>& E)
  {
    // Input D,E form a bidiagonal matrix with the first element of D = 0:
    // (eg. for N = 5)
    //     [ 0 x 0 0 0 ]
    //     [ 0 x x 0 0 ]
    // B = [ 0 0 x x 0 ]
    //     [ 0 0 0 x x ]
    //     [ 0 0 0 0 x ]
    // Zero out the first row maintaining the constancy of U B
    // using Givens transformations.
    const size_t N = D.size();
    if (N <= 1) return; 
    TMVAssert(E.size() == N-1);
    if (U) TMVAssert(U->rowsize() == N);
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);

    VIt<T,Unit,NonConj> Di = D.begin();
    VIt<T,Unit,NonConj> Ei = E.begin();
    TMVAssert(*Di == T(0));

    T x = *Ei;
    if (x != T(0)) {
      *Ei = T(0);
      ++Ei; ++Di;
      // Loop Invariant: x = B(0,i)
      for(size_t i=1; i<N; ++i,++Di,++Ei) {
	Givens<T> G = Givens_Rotate(*Di,x);
	// Make new B = G B
	if (i<N-1) G.Mult(*Ei,x);
	// Make new U = U Gt
	if (U) G.ConjMult(U->ColPair(i,0).Transpose());
      }
    }
  }

  template <class T, class TV> void BidiagonalZeroLastCol(
      const VectorView<T>& D, const VectorView<T>& E,
      const MatrixView<TV>* V)
  {
    // Input D,E form a bidiagonal matrix with the last element of D = 0:
    // (eg. for N = 5)
    //     [ x x 0 0 0 ]
    //     [ 0 x x 0 0 ]
    // B = [ 0 0 x x 0 ]
    //     [ 0 0 0 x x ]
    //     [ 0 0 0 0 0 ]
    // Zero out the last col maintaining the constancy of B V
    // using Givens transformations.
    const size_t N = D.size();
    if (N <= 1) return; 
    TMVAssert(E.size() == N-1);
    if (V) { TMVAssert(V->colsize() == N); }
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(D(N-1) == T(0));

    VIt<T,Unit,NonConj> Di = D.begin()+N-2;
    VIt<T,Unit,NonConj> Ei = E.begin()+N-2;

    T x = *Ei;
    if (x != T(0)) {
      *Ei = T(0);
      // Loop Invariant: x = B(i,N-1)
      for(int i=N-2; i>=0; --i,--Di) {
	Givens<T> G = Givens_Rotate(*Di,x);
	// Make new B = B GT
	if (i>0) G.Mult(*(--Ei),x);
	// Make new V = G* V 
	if (V) G.ConjMult(V->RowPair(i,N-1));
      }
    }
  }

  template <class T> RealType(T) BidiagonalTrailingEigenValue(
      const VectorView<T>& D, const VectorView<T>& E)
  {
    // Return the Wilkinson choice for an eigenvalue of T = BtB, namely
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
    const size_t N = D.size();
    TMVAssert(E.size() == N-1);
    TMVAssert(N > 1);

    RealType(T) a = NORM(D(N-2)) + (N>2 ? NORM(E(N-3)) : RealType(T)(0));
    RealType(T) c = NORM(D(N-1)) + NORM(E(N-2));
    T b = D(N-2)*E(N-2);
    RealType(T) bsq = NORM(b);
    RealType(T) d = (c-a)/2;
    RealType(T) dsq = SQR(d);
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

  template <class TUV, class T> void ReduceUnredBidiagonal_QR(
      const MatrixView<TUV>* U, const VectorView<T>& D,
      const VectorView<T>& E, const MatrixView<TUV>* V)
  {
    // Reduce the superdiagonal elements of unreduced Bidiagonal Matrix B 
    // (given by D,E) while maintaining U B V. 
    // Note: the input B must be unreduced - ie. all entries are non-zero.
    const size_t N = D.size();
    TMVAssert(N>0);
    TMVAssert(E.size() == N-1);
    if (U) { TMVAssert(U->rowsize() == N); }
    if (V) { TMVAssert(V->colsize() == N); }
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);

    if (N == 1) return;

    // The reduction is based on the QR algorithm to diagonalize the
    // unreduced symmetric tridiagonal matrix T = BtB
    // The basic idea is as follows:
    // (see Golub and van Loan, chapter 8 for a derivation)
    //
    // if T is a symmetric tridiagonal matrix
    // and mu is (approximately) an eigenvalue of T
    // and the QR decomposition of T - mu I = V R
    // Then, T' = R V + mu I will be tridiagonal with the last 
    // subdiagonal element small.
    // (Note: T' = Vt (T-muI) V + muI = Vt T V.)
    //
    // Wilkinson (1968) suggested that a good choice for mu is
    // the eigenvalue of the trailing 2x2 block of T that is 
    // closer to the trailing diagonal element.
    //
    // Rather than explicitly forming T = BtB and doing this
    // procedure, Golub and van Load show that it can be done
    // in place.
    // If T' = Vt T V, 
    // then Bt'B' = Vt Bt B V
    // B' = U B V for some U
    //
    // So, start with the first Givens matrix in the QR algorithm for T:
    // G0 [ T00 - mu ] = [ x ]
    //    [   T10    ]   [ 0 ]
    // We apply this to the right of B which has the effect: (for N=5)
    //              [ x x 0 0 0 ]
    //              [ + x x 0 0 ]
    // B <- B G0T = [ 0 0 x x 0 ]
    //              [ 0 0 0 x x ]
    //              [ 0 0 0 0 x ]
    // The + is the element which screws up the bidiagonal structure.
    // The rest of the procedure simply involves chasing this + down
    // the diagonal using Givens rotations.
    // For each Givens rotation we use, we also multiply U or V by the
    // adjoint to maintain the constancy of U B V.
    //
    // At the end of this procedure, E(N-1) should be smaller than it was.
    // Note: This procedure works exactly if N=2.
    VIt<T,Unit,NonConj> Di = D.begin();
    VIt<T,Unit,NonConj> Ei = E.begin();

    T mu = BidiagonalTrailingEigenValue(D,E);
    T y = NORM(*Di) - mu;  // = T00 - mu
    T x = CONJ(*Di)*(*Ei);  // = T10
    Givens<T> G = Givens_Rotate(y,x);
    for(size_t i=1;i<N;++i) {
      G.Mult(*Di,*Ei);
      if (V) G.ConjMult(V->RowPair(i-1,i));
      TMVAssert(x==T(0));
      G.Mult(x,*(++Di)); // x = B(i,i-1)
      G = Givens_Rotate(*(Di-1),x);
      G.Mult(*Ei,*Di);
      if (U) G.ConjMult(U->ColPair(i-1,i).Transpose());
      if (i < N-1) {
	TMVAssert(x==T(0));
	G.Mult(x,*(++Ei)); // x = B(i-1,i+1)
	G = Givens_Rotate(*(Ei-1),x);
      } 
    }
  }

  // MJ: This is still in the works.  Right now it just passes throug
  // to the non Divide and Conquer version.
  template <class TUV, class T> void ReduceUnredBidiagonal(
      const MatrixView<TUV>* U, const VectorView<T>& D,
      const VectorView<T>& E, const MatrixView<TUV>* V)
  {
    // Reduce the superdiagonal elements of unreduced Bidiagonal Matrix B 
    // (given by D,E) while maintaining U B V. 
    // Note: the input B must be unreduced - ie. all entries are non-zero.
    // This routine implements the divide and conquer approach, calling
    // itself for the recursion, and ReduceUnredBidiagonal_QR when the 
    // size gets too small for divide and conquer to be efficient.
    //
    // The basic idea of the divide and conquer algorithm is to split
    // the bidiagonal matrix B into two parts, B1 and B2 and a joining
    // element:
    //     [ D0 E0                  ]   [        |          ]
    //     [    D1 E1               ]   [  B1    |          ] N1
    //     [       .. ..            ]   [     Dx | Ex       ]
    // B = [          Dx Ex         ] = [-------------------]
    //     [             Dx+1 ..    ]   [        |          ]
    //     [                  .. .. ]   [        |    B2    ] N2
    //     [                     DN ]   [        |          ]
    //                                     N1         N2
    //
    //           [           |            ]
    //           [  B1t B1   |            ]
    //           [     + Dx^2|DxEx        ]
    // T = BtB = [------------------------]
    //           [       DxEx|Ex^2 +      ]
    //           [           |    B2t B2  ]
    //           [           |            ]
    // where the Dx^2 and Ex^2 are added to the lower left element of B1t B1
    // and the upper left element of B2t B2 respectively.
    //
    // This symmetric tridiagonal matrix can be written as:
    // [ T1  0 ] + w wt
    // [ 0  T2 ] 
    // where T1 = B1t B1, T2 = B2t B2 and w = (0,0,0...Dx,Ex,....0,0)
    //
    // The smaller bidiagonal matrices are then reduced recursively:
    // B1 = U1 S1 V1
    // B2 = U2 S2 V2
    //
    // The full symmetric matrix, T, is then:
    //
    // T = [ V1t S1^2 V1      0      ] + w wt
    //     [      0      V2t S2^2 V2 ]
    //
    //   = [ V1  0 ]t [ S1^2   0  ] [ V1  0 ] + w wt
    //     [ 0  V2 ]  [  0   S2^2 ] [ 0  V2 ]
    //
    //   = Vt S V + w wt
    //   = Vt (S + z zt) V
    //   where z = V w
    //
    // It turns out that it is relatively easy to find the eigenvalues
    // and eigenvectors of a rank-1 modification of a diagonal matrix.
    // (See the function SV_Rank1Solve for an explanation of this step.)
    // Thus, S^2 + z zt -> Qt S'^2 Q
    //
    // Then we take S_out = S', and V_out = Q V
    //
    // For U_out, we need to form B Bt instead:
    //
    //            [            |          ]
    //            [ B1 B1t     |          ]
    //            [            |          ]
    //            [   Dx^2+Ex^2|ExDx+1    ]
    // T~ = BBt = [-----------------------]
    //            [      ExDx+1|          ]
    //            [            |  B2 B2t  ]
    //            [            |          ]
    //
    // T~ = [ U1  0 ] [ S1^2  0  ] [ U1  0 ]t
    //      [ 0  U2 ] [  0  S2^2 ] [ 0  U2 ]
    // 
    //
    const size_t N = D.size();
    TMVAssert(N>0);
    TMVAssert(E.size() == N-1);
    if (U) { TMVAssert(U->rowsize() == N); }
    if (V) { TMVAssert(V->colsize() == N); }
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);

    //if (N < DC_LIMIT)
      return ReduceUnredBidiagonal_QR(U,D,E,V);
  }

  template <class TUV, class T> void ReduceBidiagonal(
      const MatrixView<TUV>* U, const VectorView<T>& D,
      const VectorView<T>& E, const MatrixView<TUV>* V)
  {
    // Reduce the superdiagonal elements of Bidiagonal Matrix B 
    // (given by D,E) while maintaining U B V. 
    const size_t N = D.size();
    TMVAssert(E.size() == N-1);
    if (U) { TMVAssert(U->rowsize() == N); }
    if (V) { TMVAssert(V->colsize() == N); }

    // The input E(i) are all assumed to be non-zero.
    // If there are any zeros in D, we can zero the corresponding
    // E's (above and right) directly, so look for these first.
    // Loop invariant: all D(i) with p<=i<q are non-zero.
    size_t p=0; 
    for(size_t q=0; q<N; ++q) {
      if (D(q) == T(0)) {
	if (p<q) {
	  if (V) {
	    MatrixView<TUV> V1 = V->Rows(0,q+1);
	    BidiagonalZeroLastCol(D.SubVector(p,q+1), E.SubVector(p,q), &V1);
	    MatrixView<TUV> V2 = V->Rows(p,q);
	    if (U) {
	      MatrixView<TUV> U1 = U->Cols(p,q);
	      ReduceUnredBidiagonal(&U1, D.SubVector(p,q), E.SubVector(p,q-1), 
		  &V2);
	    } else {
	      ReduceUnredBidiagonal(ZMV<TUV>(), D.SubVector(p,q), 
		  E.SubVector(p,q-1), &V2);
	    }
	  } else {
	    BidiagonalZeroLastCol(D.SubVector(p,q+1), E.SubVector(p,q),
		ZMV<TUV>());
	    if (U) {
	      MatrixView<TUV> U1 = U->Cols(p,q);
	      ReduceUnredBidiagonal(&U1, D.SubVector(p,q), E.SubVector(p,q-1), 
		  ZMV<TUV>());
	    } else {
	      ReduceUnredBidiagonal(ZMV<TUV>(), D.SubVector(p,q),
		  E.SubVector(p,q-1), ZMV<TUV>());
	    }
	  }
	}
	if (q<N-1) {
	  if (U) {
	    MatrixView<TUV> U1 = U->Cols(q,N);
	    BidiagonalZeroFirstRow(&U1, D.SubVector(q,N), E.SubVector(q,N-1));
	  } else {
	    BidiagonalZeroFirstRow(ZMV<TUV>(), D.SubVector(q,N), 
		E.SubVector(q,N-1));
	  }
	}
	p=q+1;
      }
    }
    if (p<N) {
      if (U) {
	MatrixView<TUV> U1 = U->Cols(p,N);
	if (V) {
	  MatrixView<TUV> V1 = V->Rows(p,N);
	  ReduceUnredBidiagonal(&U1, D.SubVector(p,N), E.SubVector(p,N-1), &V1);
	} else {
	  ReduceUnredBidiagonal(&U1, D.SubVector(p,N), E.SubVector(p,N-1),
	      ZMV<TUV>());
	}
      } else {
	if (V) {
	  MatrixView<TUV> V1 = V->Rows(p,N);
	  ReduceUnredBidiagonal(ZMV<TUV>(), D.SubVector(p,N), 
	      E.SubVector(p,N-1), &V1);
	} else {
	  ReduceUnredBidiagonal(ZMV<TUV>(), D.SubVector(p,N), 
	      E.SubVector(p,N-1), ZMV<TUV>());
	}
      }
    }
  }

  template <class T> void NonLapSV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V, bool SetUV)
  {
    if (SetUV) {
      TMVAssert(U && V);
      U->SetToIdentity();
      V->SetToIdentity();
    }

#ifdef XDEBUG
    cerr<<"Start Decompose from Bidiag:\n";
    if (U) cerr<<"U = "<<Type(*U)<<"  "<<*U<<endl;
    cerr<<"D = "<<Type(D)<<"  step "<<D.step()<<"  "<<D<<endl;
    cerr<<"E = "<<Type(E)<<"  step "<<E.step()<<"  "<<E<<endl;
    if (V) cerr<<"V = "<<Type(*V)<<"  "<<*V<<endl;
    cerr<<"SetUV = "<<SetUV<<endl;
    Matrix<RealType(T)> B(D.size(),D.size(),RealType(T)(0));
    B.diag() = D;
    B.diag(1) = E;
    Matrix<T> A0(U&&V ? U->colsize() : D.size(),D.size());
    if (U && V && !SetUV) A0 = (*U) * B * (*V);
    else A0 = B;
    cerr<<"A0 = "<<A0<<endl;
#endif

    TMVAssert(D.size() == E.size()+1);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
    }
    if (V) {
      TMVAssert(V->rowsize() == V->colsize());
      TMVAssert(V->rowsize() == D.size());
    }

    const size_t N = D.size();
    const size_t M = U ? U->colsize() : 0;

    // We successively reduce the superdiagonal of B (E) to 0
    // using a sequence of Givens rotations. 
    // The reduction procedure tends to push the values up and left, so it 
    // makes sense to start at the lower right and work back up the matrix.
    // We also set to zero any very small values based on machine precision.
    // Loop invariant: all E(i) with i>=q are 0.
    // Initially q = N-1. (ie. All E(i) are potentially non-zero.)
    // When q = 0, we are done.
    BidiagonalChopSmallElements(D,E);

    for(size_t q = N-1; q>0; ) {
      if (E(q-1) == T(0)) --q;
      else {
	size_t p=q-1;
	while (p > 0 && (E(p-1) != T(0))) --p; 
	// Set p such that E(p-1) = 0 and all E(i) with p<=i<q are non-zero.
	if (U)
	  if (V) {
	    MatrixView<T> U1 = U->Cols(p,q+1);
	    MatrixView<T> V1 = V->Rows(p,q+1);
	    ReduceBidiagonal(&U1,D.SubVector(p,q+1),E.SubVector(p,q),&V1);
	  } else {
	    MatrixView<T> U1 = U->Cols(p,q+1);
	    ReduceBidiagonal(&U1,D.SubVector(p,q+1),
		E.SubVector(p,q),(const MatrixView<T>*)(0));
	  }
	else
	  if (V) {
	    MatrixView<T> V1 = V->Rows(p,q+1);
	    ReduceBidiagonal((const MatrixView<T>*)(0),D.SubVector(p,q+1),
		E.SubVector(p,q),&V1);
	  } else
	    ReduceBidiagonal((const MatrixView<T>*)(0),D.SubVector(p,q+1),
		E.SubVector(p,q),(const MatrixView<T>*)(0));
	BidiagonalChopSmallElements(D,E);
      }
    }

    // Make all of the singular values positive
    VIt<RealType(T),Unit,NonConj> Di = D.begin();
    for(size_t i=0;i<N;++i,++Di) if (*Di < 0) {
      *Di = -(*Di);
      if (V) V->row(i) = -V->row(i);
    }

    // Now A = U * S * V
    // Sort output singular values 
    size_t sortp[D.size()];
    D.Sort(sortp,DESCEND);
    if (U) U->PermuteCols(sortp);
    if (V) V->PermuteRows(sortp);

#ifdef XDEBUG
    if (U && V) {
      Matrix<T> AA = (*U) * DiagMatrixViewOf(D) * (*V);
      if (Norm(AA-A0) > 0.001*Norm(A0)) {
	cerr<<"SV_DecomposeFromBidiagonal: \n";
	cerr<<"input B = "<<B<<endl;
	cerr<<"UBV = "<<A0<<endl;
	cerr<<"USV = "<<AA<<endl;
	cerr<<"U = "<<*U<<endl;
	cerr<<"S = "<<D<<endl;
	cerr<<"V = "<<*V<<endl;
	abort();
      }
    }
#endif
  }

  template <class T> void NonLapSV_Decompose_From_Bidiagonal_DelayedU(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V, bool SetUV)
  {

#ifdef XDEBUG
    cerr<<"Start Decompose from Bidiag:\n";
    if (U) cerr<<"U = "<<Type(*U)<<"  "<<*U<<endl;
    cerr<<"D = "<<Type(D)<<"  step "<<D.step()<<"  "<<D<<endl;
    cerr<<"E = "<<Type(E)<<"  step "<<E.step()<<"  "<<E<<endl;
    if (V) cerr<<"V = "<<Type(*V)<<"  "<<*V<<endl;
    cerr<<"SetUV = "<<SetUV<<endl;
    Matrix<RealType(T)> B(D.size(),D.size(),RealType(T)(0));
    B.diag() = D;
    B.diag(1) = E;
    Matrix<T> A0(U&&V ? U->colsize() : D.size(),D.size());
    if (U && V && !SetUV) A0 = (*U) * B * (*V);
    else A0 = B;
    cerr<<"A0 = "<<A0<<endl;
#endif

    TMVAssert(D.size() == E.size()+1);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
    }
    if (V) {
      TMVAssert(V->rowsize() == V->colsize());
      TMVAssert(V->rowsize() == D.size());
    }

    const size_t N = D.size();
    const size_t M = U ? U->colsize() : 0;

    Vector<RealType(T)>* D0=0;
    Vector<RealType(T)>* E0=0;
    Matrix<RealType(T)>* VV=0;
    if (SetUV) {
      D0 = new Vector<RealType(T)>(D);
      E0 = new Vector<RealType(T)>(E);
      V->SetToIdentity();
    } else if (U || V) {
      D0 = new Vector<RealType(T)>(D);
      E0 = new Vector<RealType(T)>(E);
      VV = new Matrix<RealType(T)>(N,N);
      VV->SetToIdentity();
    }

    // We successively reduce the superdiagonal of B (E) to 0
    // using a sequence of Givens rotations. 
    // The reduction procedure tends to push the values up and left, so it 
    // makes sense to start at the lower right and work back up the matrix.
    // We also set to zero any very small values based on machine precision.
    // Loop invariant: all E(i) with i>=q are 0.
    // Initially q = N-1. (ie. All E(i) are potentially non-zero.)
    // When q = 0, we are done.
    BidiagonalChopSmallElements(D,E);

    for(size_t q = N-1; q>0; ) {
      if (E(q-1) == T(0)) --q;
      else {
	size_t p=q-1;
	while (p > 0 && (E(p-1) != T(0))) --p; 
	// Set p such that E(p-1) = 0 and all E(i) with p<=i<q are non-zero.
	if (SetUV) {
	  MatrixView<T> V1 = V->Rows(p,q+1);
	  ReduceBidiagonal((const MatrixView<T>*)(0),D.SubVector(p,q+1),
	      E.SubVector(p,q),&V1);
	} else if (U || V) {
	  MatrixView<RealType(T)> V1 = VV->Rows(p,q+1);
	  ReduceBidiagonal((const MatrixView<RealType(T)>*)(0),
	      D.SubVector(p,q+1),E.SubVector(p,q),&V1);
	} else
	  ReduceBidiagonal((const MatrixView<RealType(T)>*)(0),
	      D.SubVector(p,q+1), E.SubVector(p,q),
	      (const MatrixView<RealType(T)>*)(0));
	BidiagonalChopSmallElements(D,E);
      }
    }

    // Make all of the singular values positive
    VIt<RealType(T),Unit,NonConj> Di = D.begin();
    for(size_t i=0;i<N;++i,++Di) if (*Di < 0) {
      *Di = -(*Di);
      if (SetUV) V->row(i) *= RealType(T)(-1);
      else if (U || V) VV->row(i) *= RealType(T)(-1);
    }

    // Now A = U * S * V
    // Sort output singular values 
    size_t sortp[D.size()];
    D.Sort(sortp,DESCEND);
    if (SetUV) V->PermuteRows(sortp);
    else if (U || V) VV->PermuteRows(sortp);

    if (!SetUV && V) *V = (*VV)*(*V);

    // Make U
    // US = BVt
    // Start with U = BVt
    // Then QR decompose this into US
    if (SetUV) {
      for(size_t i=0;i<N-1;++i) 
	U->row(i) = (*D0)(i) * V->col(i) + (*E0)(i) * V->col(i+1);
      U->row(N-1) = (*D0)(N-1) * V->col(N-1);
      Vector<T> Ubeta(N);
      QR_Decompose(*U,Ubeta.View());
      bool neg[N];
      VectorView<RealType(T)> Udiag = U->diag().Real();
      for(size_t i=0;i<N;++i) neg[i] = (Udiag[i] < RealType(T)(0));
      GetQFromQR(*U,Ubeta);
      for(size_t i=0;i<N;++i) if (neg[i]) U->col(i) *= RealType(T)(-1);
    } else if (U) {
      // Do everything in transposed version in this case.
      for(size_t i=0;i<N-1;++i) 
	VV->col(i) = (*D0)(i) * VV->col(i) + (*E0)(i) * VV->col(i+1);
      VV->col(N-1) *= (*D0)(N-1);
      Vector<RealType(T)> Ubeta(N);
      QR_Decompose(VV->Transpose(),Ubeta.View());
      bool neg[N];
      VectorView<RealType(T)> Udiag = VV->diag();
      for(size_t i=0;i<N;++i) neg[i] = (Udiag[i] < RealType(T)(0));
      Q_RMultEq(VV->Transpose(),Ubeta,*U);
      for(size_t i=0;i<N;++i) if (neg[i]) U->col(i) *= RealType(T)(-1);
    }

#ifdef XDEBUG
    if (U && V) {
      Matrix<T> AA = (*U) * DiagMatrixViewOf(D) * (*V);
      if (Norm(AA-A0) > 0.001*Norm(A0)) {
	cerr<<"SV_DecomposeFromBidiagonal: \n";
	cerr<<"input B = "<<B<<endl;
	cerr<<"UBV = "<<A0<<endl;
	cerr<<"USV = "<<AA<<endl;
	cerr<<"U = "<<*U<<endl;
	cerr<<"S = "<<D<<endl;
	cerr<<"V = "<<*V<<endl;
	abort();
      }
    }
#endif

    if (D0) delete D0;
    if (E0) delete E0;
    if (VV) delete VV;
  }

#ifdef LAP 
  template <class T> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V,
      bool SetUV)
  { NonLapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV); }
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<double>* U, const VectorView<double>& D, 
      const VectorView<double>& E, const MatrixView<double>* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
    }
    if (SetUV) {
      TMVAssert(U && V);
      TMVAssert(U->stor() == V->stor());
      TMVAssert(U->ct()==NonConj);
      TMVAssert(V->ct()==NonConj);
    }
    char u = 'U';
    int n = D.size();
    int lwork = (3*n+4)*n;
    double* work = LAP_DWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
    int info;
    if (SetUV) {
      char c = 'I';
      TMVAssert(U && V);
      if (U->iscm()) {
	TMVAssert(V->iscm());
	int ldu = U->stepj();
	int ldv = V->stepj();
	dbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U->ptr(),&ldu,V->ptr(),&ldv,0,0,
	    work,iwork,&info);
      } else {
	u = 'L';
	TMVAssert(U->isrm());
	TMVAssert(V->isrm());
	int ldu = U->stepi();
	int ldv = V->stepi();
	dbdsdc(&u,&c,&n,D.ptr(),E.ptr(),V->ptr(),&ldv,U->ptr(),&ldu,0,0,
	    work,iwork,&info);
      }
    } else if (U || V) {
      char c = 'I';
      Matrix<double,ColMajor> U1(n,n);
      Matrix<double,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      dbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,V1.ptr(),&ldv,0,0,
	  work,iwork,&info);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
      dbdsdc(&u,&c,&n,D.ptr(),E.ptr(),0,&ldu,0,&ldv,0,0,
	  work,iwork,&info);
    }
    LAP_Results(info,"dbdsdc");
  }
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<complex<double> >* U, const VectorView<double>& D, 
      const VectorView<double>& E, const MatrixView<complex<double> >* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
    }
    TMVAssert(!SetUV);

    char u = 'U';
    int n = D.size();
    int lwork = (3*n+4)*n;
    double* work = LAP_DWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
    int info;
    if (U || V) {
      char c = 'I';
      Matrix<double,ColMajor> U1(n,n);
      Matrix<double,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      dbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,V1.ptr(),&ldv,0,0,
	  work,iwork,&info);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
      dbdsdc(&u,&c,&n,D.ptr(),E.ptr(),0,&ldu,0,&ldv,0,0,
	  work,iwork,&info);
    }
    LAP_Results(info,"dbdsdc");
  }
#ifndef NOFLOAT
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<float>* U, const VectorView<float>& D, 
      const VectorView<float>& E, const MatrixView<float>* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
    }
    if (SetUV) {
      TMVAssert(U && V);
      TMVAssert(U->stor() == V->stor());
      TMVAssert(U->ct()==NonConj);
      TMVAssert(V->ct()==NonConj);
    }
    char u = 'U';
    int n = D.size();
    int lwork = (3*n+4)*n;
    float* work = LAP_SWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
    int info;
    if (SetUV) {
      char c = 'I';
      TMVAssert(U && V);
      if (U->iscm()) {
	TMVAssert(V->iscm());
	int ldu = U->stepj();
	int ldv = V->stepj();
	sbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U->ptr(),&ldu,V->ptr(),&ldv,0,0,
	    work,iwork,&info);
      } else {
	u = 'L';
	TMVAssert(U->isrm());
	TMVAssert(V->isrm());
	int ldu = U->stepi();
	int ldv = V->stepi();
	sbdsdc(&u,&c,&n,D.ptr(),E.ptr(),V->ptr(),&ldv,U->ptr(),&ldu,0,0,
	    work,iwork,&info);
      }
    } else if (U || V) {
      char c = 'I';
      Matrix<float,ColMajor> U1(n,n);
      Matrix<float,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      sbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,V1.ptr(),&ldv,0,0,
	  work,iwork,&info);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      char c = 'N';
      int ldu = n;
      int ldv = n;
      sbdsdc(&u,&c,&n,D.ptr(),E.ptr(),0,&ldu,0,&ldv,0,0,
	  work,iwork,&info);
    }
    LAP_Results(info,"sbdsdc");
  }
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<complex<float> >* U, const VectorView<float>& D, 
      const VectorView<float>& E, const MatrixView<complex<float> >* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
    }
    TMVAssert(!SetUV);
    char u = 'U';
    int n = D.size();
    int lwork = (3*n+4)*n;
    float* work = LAP_SWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
    int info;
    if (U || V) {
      char c = 'I';
      Matrix<float,ColMajor> U1(n,n);
      Matrix<float,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      sbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,V1.ptr(),&ldv,0,0,
	  work,iwork,&info);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      char c = 'N';
      int ldu = n;
      int ldv = n;
      sbdsdc(&u,&c,&n,D.ptr(),E.ptr(),0,&ldu,0,&ldv,0,0,
	  work,iwork,&info);
    }
    LAP_Results(info,"sbdsdc");
  }
#endif
#endif // LAP

  template <class T> void SV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->ct()==NonConj);
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
      TMVAssert(V->ct()==NonConj);
    }
    
    if (D.size() > 0) {
#ifdef LAP
      if ((!U || U->iscm() || U->isrm()) && 
	  D.step() == 1 && E.step() == 1 && 
	  (!V || V->iscm() || V->isrm()) && 
	  (!U || !V || !SetUV || U->stor()==V->stor())) {
	LapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV);
      }
      else 
#endif
	//NonLapSV_Decompose_From_Bidiagonal_DelayedU(U,D,E,V,SetUV);
	NonLapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV);
    }
  }

  //
  // Main SVD Drivers
  // (With and without explicit formation of U,V)
  //
  
  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>* V, T& det, bool StoreU)
  {
    // Decompose A (input as U) into U S V
    // where S is a diagonal real matrix, and U,V are unitary matrices.
    // A,U are M x N (M >= N)
    // S,V are N x N
    // The determinant is returned in det.
    // (Technically, det is multiplied by the determinant, so det should
    // be set to 1 on entry.)
    const size_t M = U.colsize();
    const size_t N = U.rowsize();
    if (N == 0) return;

    TMVAssert(N <= M);
    if (V) {
      TMVAssert(V->colsize() == N);
      TMVAssert(V->rowsize() == N);
    }
    TMVAssert(S.size() == N);

    // If M is much larger than N (technically M > 5/3 N), then it is quicker
    // to start by doing a QR decomposition and then do SVD on the square
    // R matrix.  Thus, the final U of the SVD is Q (from the QR decomp)
    // times U from R's SVD.
    if (M > 5*N/3) {
      if (StoreU) {
	Matrix<T,ColMajor> R(N,N);
	QR_Decompose(U,R.View(),det);
	SV_Decompose(R.View(),S,V,det,StoreU);
	// Now R is a Unitary Matrix U'.  Need to multiply U by U'
	U = U*R;
      } else {
	Vector<T> Qbeta(N);
	QR_Decompose(U,Qbeta.View(),det);
	LowerTriMatrixViewOf(U.Rows(0,N)).OffDiag().Zero();
	SV_Decompose(U.Rows(0,N),S,V,det,StoreU);
      }
    } else {
      // First we reduce A to bidiagonal form: A = U * B * V
      // using a series of Householder transformations.
      // The diagonal of the Bidiagonal Matrix B is stored in D.
      // The superdiagonal is stored in E.
      Vector<RealType(T)> E(N-1);
      Vector<T> Ubeta(N);
      Vector<T> Vbeta(N-1);
      Bidiagonalize(U,Ubeta.View(),Vbeta.View(),S,E.View(),det);
      // The determinant of B is just the product of the diagonal elements:
      if (det!=T(0)) det *= DiagMatrixViewOf(S).Det();

      // Now UV stores Householder vectors for U in lower diagonal columns 
      // (HLi) and Householder vectors for V in upper diagonal rows (HRi)
      // The Householder matrices for U are actually the adjoints of the 
      // matrices that bidiagonalize A, and for V are the transposes:
      // U = HLn-1t ... HL1t HL0t A HR0T HR1T ... HRn-2T
      // Using the fact that H Ht = I, we get A = U B V with:
      // U = HL0 ... HLn-1 
      if (V) {
	V->row(0).MakeBasis(0);
	V->Rows(1,N) = U.Rows(0,N-1);
	V->col(0,1,N).Zero();
	GetQFromQR(V->SubMatrix(1,N,1,N).Transpose(),Vbeta);
      }
      if (StoreU) {
	GetQFromQR(U,Ubeta);
      }

      if (StoreU) SV_Decompose_From_Bidiagonal(&U,S,E.View(),V);
      else SV_Decompose_From_Bidiagonal((const MatrixView<T>*)(0),S,E.View(),V);

    }
  }

  template <class T> void SV_Decompose(
      const MatrixView<T>& UV,
      const VectorView<T>& Ubeta, const VectorView<T>& Vbeta,
      Matrix<T,ColMajor>*& UV2, Vector<T>*& Qbeta,
      const MatrixView<RealType(T)>& U1, const MatrixView<RealType(T)>& V1,
      const VectorView<RealType(T)>& S, T& det)
  {
    // Decompose A (input as UV) into U U1 S V1 V
    // with U,V stored in UV, Ubeta, Vbeta
    // or Q U U1 S V1 V with Q stored in UV, Qbeta and U,V stored
    // in UV2, Ubeta, Vbeta
    
    const size_t M = UV.colsize();
    const size_t N = UV.rowsize();
    if (N == 0) return;

    TMVAssert(N <= M);
    TMVAssert(Ubeta.size() == N);
    TMVAssert(Vbeta.size() == N-1);
    TMVAssert(U1.colsize() == N);
    TMVAssert(U1.rowsize() == N);
    TMVAssert(V1.colsize() == N);
    TMVAssert(V1.rowsize() == N);
    TMVAssert(S.size() == N);

    if (M > 5*N/3) {
      UV2 = new Matrix<T,ColMajor>(N,N);
      Qbeta = new Vector<T>(N);

      QR_Decompose(UV,Qbeta->View(),det);
      
      *UV2 = UpperTriMatrixViewOf(UV);

      Matrix<T,ColMajor>* temp1=0;
      Vector<T>* temp2=0;
      SV_Decompose(UV2->View(),Ubeta,Vbeta,temp1,temp2,U1,V1,S,det);
      TMVAssert(!temp1 && !temp2);
    } else {
      Vector<RealType(T)> E(N-1);
      Bidiagonalize(UV,Ubeta,Vbeta,S,E.View(),det);
      if (det!=T(0)) det *= DiagMatrixViewOf(S).Det();

      SV_Decompose_From_Bidiagonal(&U1,S,E.View(),&V1,true);
    }
  }

  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>& V, T& det, bool StoreU)
  { SV_Decompose(U,S,&V,det,StoreU); }

  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      T& det, bool StoreU)
  { SV_Decompose(U,S,(const MatrixView<T>*)(0),det,StoreU); }

#define InstFile "TMV_SVDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


