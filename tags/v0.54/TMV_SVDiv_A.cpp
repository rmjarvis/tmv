
#include "TMV.h"
#include "TMV_Householder.h"
#include "TMV_Givens.h"
#include "TMV_Diag.h"
#include "TMV_Tri.h"

//#define XDEBUG
//#define TIME

#ifdef TIME
#include <sys/time.h>
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define BIDIAG_BLOCKSIZE TMV_BLOCKSIZE
#define DC_LIMIT TMV_BLOCKSIZE
#else
#define BIDIAG_BLOCKSIZE 16
#define DC_LIMIT 32
#endif

  template <class T> inline const MatrixView<T>* ZMV() 
  { return (const MatrixView<T>*)(0); }

  //
  // Bidiagonalize
  //
  
  template <class T> inline void NonBlockBidiagonalize(
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
    TMVAssert(A.iscm() || A.isrm());

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

  template <class T> inline void BlockBidiagonalize(
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

    Matrix<T,RowMajor> ZYtm(min(size_t(BIDIAG_BLOCKSIZE),N-1),N);
    Matrix<T,ColMajor> mXtW(M,min(size_t(BIDIAG_BLOCKSIZE),N-1));
    for(size_t j1=0;j1<N-1;) {
      size_t j2 = min(N-1,j1+BIDIAG_BLOCKSIZE);
      for(size_t j=j1,jj=0;j<j2;++j,++jj) { // jj = j-j1

	// Update current column:
	// A(j:M,j) -= Y(j:M,0:j) ZYtm(0:j,j) + mXtW(j:M,0:j) X(0:j,j)
	//
	VectorView<T> u = A.col(j,j,M);
	MatrixView<T> Y0 = A.SubMatrix(j,M,j1,j);
	MatrixView<T> ZYtm0 = ZYtm.SubMatrix(0,jj,j,N);
	MatrixView<T> mXtW0 = mXtW.SubMatrix(j,M,0,jj);
	MatrixView<T> X0 = A.SubMatrix(j1,j,j,N);
	if (jj > 0) {
	  u -= Y0 * ZYtm0.col(0);
	  u -= mXtW0 * X0.col(0);
	}

	// Do the Householder reflection for U
	// Copy the reflection into D(j), and set the top of the 
	// Householder vector to be explicitly 1.  (It makes life easier
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
	VectorView<T> ZYtmj = ZYtm.row(jj,j+1,N);
	VectorView<T> temp = ZYtm.row(jj,j1,j);
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
	MatrixView<T> ZYtm1 = ZYtm.SubMatrix(0,jj+1,j+1,N);
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
	VectorView<T> mXtWj = mXtW.col(jj,j+1,M);
	VectorView<T> temp1 = mXtW.col(jj,j1,j+1);
	VectorView<T> temp2 = mXtW.col(jj,j1,j);
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

  template <class T> inline void NonLapBidiagonalize(
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
  template <class T> inline void LapBidiagonalize(
      const MatrixView<T>& A, const VectorView<T>& Ubeta,
      const VectorView<T>& Vbeta, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, T& det)
  { NonLapBidiagonalize(A,Ubeta,Vbeta,D,E,det); }
  template <> inline void LapBidiagonalize(
      const MatrixView<double>& A, const VectorView<double>& Ubeta,
      const VectorView<double>& Vbeta, const VectorView<double>& D,
      const VectorView<double>& E, double& det)
  {
#ifdef TIME
    cerr<<"LapBidiag: \n";
    timeval tp;
    gettimeofday(&tp,0);
    double t0 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
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
#ifndef LAPNOWORK
    int lwork = (m+n)*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
#endif
#ifdef TIME
    gettimeofday(&tp,0);
    double t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
    LAPNAME(dgebrd) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
	LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
	LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef TIME
    gettimeofday(&tp,0);
    double t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
    Vbeta = Vbeta2.SubVector(0,n-1);
#ifdef LAPNOWORK
    LAP_Results("dgebrd");
#else
    LAP_Results(int(work[0]),m,n,lwork,"dgebrd");
#endif
    if (det) {
      for(size_t i=0;i<Ubeta.size();++i) if (Ubeta(i) != double(0)) 
	det = -det;
      for(size_t i=0;i<Vbeta.size();++i) if (Vbeta(i) != double(0)) 
	det = -det;
    }
#ifdef TIME
    gettimeofday(&tp,0);
    double t3 = tp.tv_sec + tp.tv_usec/1.e6;
    cerr<<"  dbebrd time = "<<t2-t1<<endl;
    cerr<<"  other time = "<<t3-t2+t1-t0<<endl;
#endif
  }
  template <> inline void LapBidiagonalize(
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
#ifndef LAPNOWORK
    int lwork = (m+n)*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
#endif
    LAPNAME(zgebrd) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
	LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
	LAPWK(work) LAPVWK(lwork) LAPINFO);
    Ubeta.ConjugateSelf();
    Vbeta = Vbeta2.SubVector(0,n-1);
#ifdef LAPNOWORK
    LAP_Results("zgebrd");
#else
    LAP_Results(int(real(work[0])),m,n,lwork,"zgebrd");
#endif
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
  template <> inline void LapBidiagonalize(
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
#ifndef LAPNOWORK
    int lwork = (m+n)*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
#endif
    LAPNAME(sgebrd) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
	LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
	LAPWK(work) LAPVWK(lwork) LAPINFO);
    Vbeta = Vbeta2.SubVector(0,n-1);
#ifdef LAPNOWORK
    LAP_Results("sgebrd");
#else
    LAP_Results(int(work[0]),m,n,lwork,"sgebrd");
#endif
    if (det) {
      for(size_t i=0;i<Ubeta.size();++i) if (Ubeta(i) != float(0)) 
	det = -det;
      for(size_t i=0;i<Vbeta.size();++i) if (Vbeta(i) != float(0)) 
	det = -det;
    }
  }
  template <> inline void LapBidiagonalize(
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
#ifndef LAPNOWORK
    int lwork = (m+n)*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
#endif
    LAPNAME(cgebrd) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
	LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
	LAPWK(work) LAPVWK(lwork) LAPINFO);
    Ubeta.ConjugateSelf();
    Vbeta = Vbeta2.SubVector(0,n-1);
#ifdef LAPNOWORK
    LAP_Results("cgebrd");
#else
    LAP_Results(int(real(work[0])),m,n,lwork,"cgebrd");
#endif
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

  template <class T> inline void Bidiagonalize(
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
  
  template <class T> inline void BidiagonalChopSmallElements(
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

    T* Di = D.ptr();
    T* Ei = E.ptr();
    RealType(T) absDim1 = abs(*Di);
    if (absDim1 < dthresh) *Di = T(0);
    ++Di;
    for(size_t k=E.size();k>0;--k,++Di,++Ei) {
      RealType(T) absDi = abs(*Di);
      if (absDi < dthresh) *Di = T(0);
      RealType(T) ethresh = eps * (absDi + absDim1);
      if (ethresh == RealType(T)(0)) ethresh = dthresh;
      if (abs(*Ei) < ethresh) *Ei = T(0);
      absDim1 = absDi;
    }
  }

  template <class T, class TU> inline void BidiagonalZeroFirstRow(
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

    T* Di = D.ptr();
    T* Ei = E.ptr();
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

  template <class T, class TV> inline void BidiagonalZeroLastCol(
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
    if (V) TMVAssert(V->colsize() == N);
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(D(N-1) == T(0));

    T* Di = D.ptr()+N-2;
    T* Ei = E.ptr()+N-2;

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

  template <class T> inline RealType(T) BidiagonalTrailingEigenValue(
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

  template <class TUV, class T> inline void ReduceUnredBidiagonal(
      const MatrixView<TUV>* U, const VectorView<T>& D,
      const VectorView<T>& E, const MatrixView<TUV>* V)
  {
    // Reduce the superdiagonal elements of unreduced Bidiagonal Matrix B 
    // (given by D,E) while maintaining U B V. 
    // Note: the input B must be unreduced - ie. all entries are non-zero.
    const size_t N = D.size();
    TMVAssert(N>0);
    TMVAssert(E.size() == N-1);
    if (U) TMVAssert(U->rowsize() == N);
    if (V) TMVAssert(V->colsize() == N);
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
    T* Di = D.ptr();
    T* Ei = E.ptr();

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

  template <class TUV, class T> inline void ReduceBidiagonal(
      const MatrixView<TUV>* U, const VectorView<T>& D,
      const VectorView<T>& E, const MatrixView<TUV>* V)
  {
    // Reduce the superdiagonal elements of Bidiagonal Matrix B 
    // (given by D,E) while maintaining U B V. 
    const size_t N = D.size();
    TMVAssert(E.size() == N-1);
    if (U) TMVAssert(U->rowsize() == N); 
    if (V) TMVAssert(V->colsize() == N);

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

  template <class T> inline void SV_Decompose_From_Bidiagonal_QR(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V)
  {
    TMVAssert(D.size() == E.size()+1);
    if (U) {
      TMVAssert(U->rowsize() == D.size());
    }
    if (V) {
      TMVAssert(V->colsize() == D.size());
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

    for(size_t q = E.size(); q>0; ) {
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
	    ReduceBidiagonal(&U1,D.SubVector(p,q+1),E.SubVector(p,q),ZMV<T>());
	  }
	else
	  if (V) {
	    MatrixView<T> V1 = V->Rows(p,q+1);
	    ReduceBidiagonal(ZMV<T>(),D.SubVector(p,q+1),E.SubVector(p,q),&V1);
	  } else
	    ReduceBidiagonal(ZMV<T>(),D.SubVector(p,q+1),E.SubVector(p,q),
		ZMV<T>());
	BidiagonalChopSmallElements(D,E);
      }
    }
  }

  // MJ: This is unfinished - write divide and conquer version.
  template <class T> inline void SV_Decompose_From_Bidiagonal_DC(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, const MatrixView<T>* V)
  {
    // Solve the SVD of unreduced Bidiagonal Matrix B (given by D,E).
    // Note: the input B must be unreduced - ie. all entries are non-zero.
    // This routine implements the divide and conquer approach, calling
    // itself for the recursion, and ReduceUnredBidiagonal_QR when the
    // size gets too small for divide and conquer to be efficient.
    //
    // The basic idea of the divide and conquer algorithm is to split
    // the bidiagonal matrix B into two parts, B1 and B2 and a joining
    // element:
    //     [ D0 E0                    ]   [        |          ]
    //     [    D1 E1                 ]   [  B1    |          ] N1+1
    //     [       .. ..              ]   [     Dx | Ex       ]
    // B = [          Dx Ex           ] = [-------------------]
    //     [             D(x+1) ..    ]   [        |          ]
    //     [                    .. .. ]   [        |    B2    ] N2
    //     [                       DN ]   [        |          ]
    //                                       N1+1       N2
    //
    // The smaller bidiagonal matrices are first reduced recursively:
    // B1 = U1 S1 V1
    // B2 = U2 S2 V2
    //
    //
    //     [ B1    |     ]   [ U1 S1 V1    |            ]
    //     [    Dx | Ex  ]   [          Dx | Ex         ]
    // B = [-------------] = [--------------------------]
    //     [       |     ]   [             |            ]
    //     [       |  B2 ]   [             |   U2 S2 V2 ]
    //
    //     [ U1 0  0  ] [ S1 0  0  ] [  V1   0  ]
    //   = [ 0  1  0  ] [ w1 wx w2 ] [ 0  1  0  ]
    //     [ 0  0  U2 ] [ 0  0  S2 ] [ 0  0  V2 ]
    //
    // Note that U1, U2, V2 are all square, but V1 is not, since B1 is not.
    // The vector w = [ w1 wx w2 ] is such that w V = [ 0 Dx Ex 0 ]
    //
    //
    TMVAssert(D.size()>0);
    TMVAssert(E.size()+1 == D.size());
    if (U) TMVAssert(U->rowsize() == D.size()); 
    if (V) TMVAssert(V->colsize() == D.size()); 
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);

  }
   
  template <class T> inline void DoSV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V)
  {
    SV_Decompose_From_Bidiagonal_QR(U,D,E,V);
  }

  /*
  template <class T> inline void NonLapSV_Decompose_From_Bidiagonal_DelayedU(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V, bool SetUV)
  {

#ifdef XDEBUG
    //cerr<<"Start Decompose from Bidiag:\n";
    //if (U) cerr<<"U = "<<Type(*U)<<"  "<<*U<<endl;
    //cerr<<"D = "<<Type(D)<<"  step "<<D.step()<<"  "<<D<<endl;
    //cerr<<"E = "<<Type(E)<<"  step "<<E.step()<<"  "<<E<<endl;
    //if (V) cerr<<"V = "<<Type(*V)<<"  "<<*V<<endl;
    //cerr<<"SetUV = "<<SetUV<<endl;
    Matrix<RealType(T)> B(D.size(),D.size(),RealType(T)(0));
    B.diag() = D;
    B.diag(1) = E;
    Matrix<T> A0(U&&V ? U->colsize() : D.size(),D.size());
    if (U && V && !SetUV) A0 = (*U) * B * (*V);
    else A0 = B;
    //cerr<<"A0 = "<<A0<<endl;
#endif

#ifdef TIME
    timeval tp;
    gettimeofday(&tp,0);
    double t0 = tp.tv_sec + tp.tv_usec/1.e6;
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

    auto_ptr<Vector<RealType(T)> > D0=0;
    auto_ptr<Vector<RealType(T)> > E0=0;
    auto_ptr<Matrix<RealType(T)> > VV=0;
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

#ifdef TIME
    gettimeofday(&tp,0);
    double t1 = tp.tv_sec + tp.tv_usec/1.e6;
    cerr<<"time for first chop: "<<t1-t0<<endl;
    int nloops = 0;
    double t_a=0., t_b=0., t_c=0.;
#endif

    for(size_t q = N-1; q>0; ) {
#ifdef TIME
      ++nloops;
      gettimeofday(&tp,0);
      double xt0 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      if (E(q-1) == T(0)) --q;
      else {
	size_t p=q-1;
	while (p > 0 && (E(p-1) != T(0))) --p; 
#ifdef TIME
	gettimeofday(&tp,0);
	double xt1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
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
#ifdef TIME
	gettimeofday(&tp,0);
	double xt2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
	BidiagonalChopSmallElements(D,E);
#ifdef TIME
	gettimeofday(&tp,0);
	double xt3 = tp.tv_sec + tp.tv_usec/1.e6;
	t_a += xt1-xt0;
	t_b += xt2-xt1;
	t_c += xt3-xt2;
#endif
      }
    }
#ifdef TIME
    cerr<<"nloops = "<<nloops<<endl;
    gettimeofday(&tp,0);
    double t2 = tp.tv_sec + tp.tv_usec/1.e6;
    cerr<<"time : "<<t2-t1<<endl;
    cerr<<"time for --p step = "<<t_a<<endl;
    cerr<<"time for reduce bidiag = "<<t_b<<endl;
    cerr<<"time for chop = "<<t_c<<endl;
#endif

    // Make all of the singular values positive
    RealType(T)* Di = D.ptr();
    for(size_t i=0;i<N;++i,++Di) if (*Di < 0) {
      *Di = -(*Di);
      if (SetUV) V->row(i) *= RealType(T)(-1);
      else if (U || V) VV->row(i) *= RealType(T)(-1);
    }

#ifdef TIME
    gettimeofday(&tp,0);
    double t3 = tp.tv_sec + tp.tv_usec/1.e6;
    cerr<<"time for positive-ify : "<<t3-t2<<endl;
#endif

    // Now A = U * S * V
    // Sort output singular values 
    size_t sortp[D.size()];
    D.Sort(sortp,DESCEND);
    if (SetUV) V->PermuteRows(sortp);
    else if (U || V) VV->PermuteRows(sortp);

#ifdef TIME
    gettimeofday(&tp,0);
    double t4 = tp.tv_sec + tp.tv_usec/1.e6;
    cerr<<"time for sort : "<<t4-t3<<endl;
#endif

    if (!SetUV && V) *V = (*VV)*(*V);

#ifdef TIME
    gettimeofday(&tp,0);
    double t5 = tp.tv_sec + tp.tv_usec/1.e6;
    cerr<<"time for mult V: "<<t5-t4<<endl;
#endif

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
      //GetQFromQR(VV->Transpose(),Ubeta);
      //for(size_t i=0;i<N;++i) if (neg[i]) VV->row(i) *= RealType(T)(-1);
      // *U = *U * VV->Transpose();
    }
#ifdef TIME
    gettimeofday(&tp,0);
    double t6 = tp.tv_sec + tp.tv_usec/1.e6;
    cerr<<"time for make U: "<<t6-t5<<endl;
#endif

#ifdef XDEBUG
    if (U && V) {
      Matrix<T> AA = (*U) * DiagMatrixViewOf(D) * (*V);
      if (Norm(A0-AA) > 0.001*Norm(A0)) {
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

  }*/

  template <class T> inline void NonLapSV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V, bool SetUV)
  {
#ifdef XDEBUG
    //cerr<<"Start Decompose from Bidiag:\n";
    //if (U) cerr<<"U = "<<Type(*U)<<"  "<<*U<<endl;
    //cerr<<"D = "<<Type(D)<<"  step "<<D.step()<<"  "<<D<<endl;
    //cerr<<"E = "<<Type(E)<<"  step "<<E.step()<<"  "<<E<<endl;
    //if (V) cerr<<"V = "<<Type(*V)<<"  "<<*V<<endl;
    //cerr<<"SetUV = "<<SetUV<<endl;
    Matrix<RealType(T)> B(D.size(),D.size(),RealType(T)(0));
    B.diag() = D;
    B.diag(1) = E;
    Matrix<T> A0(U&&V ? U->colsize() : D.size(),D.size());
    if (U && V && !SetUV) A0 = (*U) * B * (*V);
    else A0 = B;
    //cerr<<"A0 = "<<A0<<endl;
#endif

    const size_t N = D.size();

    if (SetUV) {
      TMVAssert(U && V);
      U->SetToIdentity();
      V->SetToIdentity();
    }

    // First chop any small elements in D,E
    BidiagonalChopSmallElements(D,E);

    // Find sub-problems to solve:
    for(size_t q = N-1; q>0; ) {
      if (E(q-1) == T(0)) --q;
      else {
	// Set p such that E(p-1) = 0 and all E(i) with p<=i<q are non-zero.
	size_t p=q-1;
	while (p > 0 && (E(p-1) != T(0))) --p; 
	if (U)
	  if (V) {
	    MatrixView<T> U1 = U->Cols(p,q+1);
	    MatrixView<T> V1 = V->Rows(p,q+1);
	    DoSV_Decompose_From_Bidiagonal(&U1,D.SubVector(p,q+1),
		E.SubVector(p,q),&V1);
	  } else {
	    MatrixView<T> U1 = U->Cols(p,q+1);
	    DoSV_Decompose_From_Bidiagonal(&U1,D.SubVector(p,q+1),
		E.SubVector(p,q),(const MatrixView<T>*)(0));
	  }
	else
	  if (V) {
	    MatrixView<T> V1 = V->Rows(p,q+1);
	    DoSV_Decompose_From_Bidiagonal((const MatrixView<T>*)(0),
		D.SubVector(p,q+1),E.SubVector(p,q),&V1);
	  } else {
	    DoSV_Decompose_From_Bidiagonal((const MatrixView<T>*)(0),
		D.SubVector(p,q+1),E.SubVector(p,q),
		(const MatrixView<T>*)(0));
	  }
	q = p > 0 ? p-1 : 0;
      }
    }

    // Make all of the singular values positive
    RealType(T)* Di = D.ptr();
    for(size_t i=0;i<N;++i,++Di) if (*Di < 0) {
      *Di = -(*Di);
      if (V) V->row(i) = -V->row(i);
    }

    // Now A = U * S * V
    // Sort output singular values 
    auto_array<size_t> sortp(new size_t[N]);
    D.Sort(sortp.get(),DESCEND);
    if (U) U->PermuteCols(sortp.get());
    if (V) V->PermuteRows(sortp.get());

#ifdef XDEBUG
    if (U && V) {
      Matrix<T> AA = (*U) * DiagMatrixViewOf(D) * (*V);
      if (Norm(A0-AA) > 0.001*Norm(A0)) {
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

#ifdef LAP 
  template <class T> inline void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V,
      bool SetUV)
  { NonLapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV); }
  template <> inline void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<double>* U, const VectorView<double>& D, 
      const VectorView<double>& E, const MatrixView<double>* V,
      bool SetUV)
  {
#ifdef TIME
    timeval tp;
    gettimeofday(&tp,0);
    double t0 = tp.tv_sec + tp.tv_usec/1.e6;
    double t1,t2;
    cerr<<"LapSV_Decompose_From_Bidiagonal\n";
    cerr<<"  U,V,SetUV = "<<bool(U)<<','<<bool(V)<<','<<SetUV<<endl;
#endif
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
      if (U && SetUV) TMVAssert(U->stor() == V->stor());
    }
    char u = 'U';
    int n = D.size();
#ifndef LAPNOWORK
    int lwork = (3*n+4)*n;
    double* work = LAP_DWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
#endif
    if (SetUV) {
      char c = 'I';
      TMVAssert(U && V);
      if (U->iscm()) {
	TMVAssert(V->iscm());
	int ldu = U->stepj();
	int ldv = V->stepj();
#ifdef TIME
	gettimeofday(&tp,0);
	t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
	LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	    LAPP(D.ptr()),LAPP(E.ptr()),
	    LAPP(U->ptr()),LAPV(ldu),LAPP(V->ptr()),LAPV(ldv),0,0
	    LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
	gettimeofday(&tp,0);
	t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      } else {
	u = 'L';
	TMVAssert(U->isrm());
	TMVAssert(V->isrm());
	int ldu = U->stepi();
	int ldv = V->stepi();
#ifdef TIME
	gettimeofday(&tp,0);
	t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
	LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	    LAPP(D.ptr()),LAPP(E.ptr()),
	    LAPP(V->ptr()),LAPV(ldv),LAPP(U->ptr()),LAPV(ldu),0,0
	    LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
	gettimeofday(&tp,0);
	t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      }
    } else if (U || V) {
      char c = 'I';
      Matrix<double,ColMajor> U1(n,n);
      Matrix<double,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
#ifdef TIME
      gettimeofday(&tp,0);
      t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
      gettimeofday(&tp,0);
      t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
#ifdef TIME
      gettimeofday(&tp,0);
      t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  0,LAPV(ldu),0,LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
      gettimeofday(&tp,0);
      t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
    }
    LAP_Results("dbdsdc");
#ifdef TIME
    gettimeofday(&tp,0);
    double t3 = tp.tv_sec + tp.tv_usec/1.e6;
    cerr<<"  dbdsdc time = "<<t2-t1<<endl;
    cerr<<"  other time = "<<t3-t2+t1-t0<<endl;
#endif
  }
  template <> inline void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<complex<double> >* U, const VectorView<double>& D, 
      const VectorView<double>& E, const MatrixView<complex<double> >* V, 
      bool
#ifdef TMVDEBUG
      SetUV
#endif
      )
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
    TMVAssert(!SetUV);

    char u = 'U';
    int n = D.size();
#ifndef LAPNOWORK
    int lwork = (3*n+4)*n;
    double* work = LAP_DWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
#endif
    if (U || V) {
      char c = 'I';
      Matrix<double,ColMajor> U1(n,n);
      Matrix<double,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
      LAPNAME(dbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  0,LAPV(ldu),0,LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
    }
    LAP_Results("dbdsdc");
  }
#ifndef NOFLOAT
  template <> inline void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<float>* U, const VectorView<float>& D, 
      const VectorView<float>& E, const MatrixView<float>* V,
      bool SetUV)
  {
#ifdef TIME
    timeval tp;
    gettimeofday(&tp,0);
    float t0 = tp.tv_sec + tp.tv_usec/1.e6;
    float t1,t2;
    cerr<<"LapSV_Decompose_From_Bidiagonal\n";
    cerr<<"  U,V,SetUV = "<<bool(U)<<','<<bool(V)<<','<<SetUV<<endl;
#endif
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
      if (U && SetUV) TMVAssert(U->stor() == V->stor());
    }
    char u = 'U';
    int n = D.size();
#ifndef LAPNOWORK
    int lwork = (3*n+4)*n;
    float* work = LAP_SWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
#endif
    if (SetUV) {
      char c = 'I';
      TMVAssert(U && V);
      if (U->iscm()) {
	TMVAssert(V->iscm());
	int ldu = U->stepj();
	int ldv = V->stepj();
#ifdef TIME
	gettimeofday(&tp,0);
	t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
	LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	    LAPP(D.ptr()),LAPP(E.ptr()),
	    LAPP(U->ptr()),LAPV(ldu),LAPP(V->ptr()),LAPV(ldv),0,0
	    LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
	gettimeofday(&tp,0);
	t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      } else {
	u = 'L';
	TMVAssert(U->isrm());
	TMVAssert(V->isrm());
	int ldu = U->stepi();
	int ldv = V->stepi();
#ifdef TIME
	gettimeofday(&tp,0);
	t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
	LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),LAPP(D.ptr()),LAPP(E.ptr()),
	    LAPP(V->ptr()),LAPV(ldv),LAPP(U->ptr()),LAPV(ldu),0,0
	    LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
	gettimeofday(&tp,0);
	t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      }
    } else if (U || V) {
      char c = 'I';
      Matrix<float,ColMajor> U1(n,n);
      Matrix<float,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
#ifdef TIME
      gettimeofday(&tp,0);
      t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
      gettimeofday(&tp,0);
      t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
#ifdef TIME
      gettimeofday(&tp,0);
      t1 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
      LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  0,LAPV(ldu),0,LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
#ifdef TIME
      gettimeofday(&tp,0);
      t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
    }
    LAP_Results("sbdsdc");
#ifdef TIME
    gettimeofday(&tp,0);
    float t3 = tp.tv_sec + tp.tv_usec/1.e6;
    cerr<<"  sbdsdc time = "<<t2-t1<<endl;
    cerr<<"  other time = "<<t3-t2+t1-t0<<endl;
#endif
  }
  template <> inline void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<complex<float> >* U, const VectorView<float>& D, 
      const VectorView<float>& E, const MatrixView<complex<float> >* V, 
      bool
#ifdef TMVDEBUG
      SetUV
#endif
      )
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
    TMVAssert(!SetUV);

    char u = 'U';
    int n = D.size();
#ifndef LAPNOWORK
    int lwork = (3*n+4)*n;
    float* work = LAP_SWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
#endif
    if (U || V) {
      char c = 'I';
      Matrix<float,ColMajor> U1(n,n);
      Matrix<float,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
      LAPNAME(sbdsdc) (LAPCM LAPV(u),LAPV(c),LAPV(n),
	  LAPP(D.ptr()),LAPP(E.ptr()),
	  0,LAPV(ldu),0,LAPV(ldv),0,0
	  LAPWK(work) LAPWK(iwork) LAPINFO LAP1 LAP1);
    }
    LAP_Results("sbdsdc");
  }
#endif // FLOAT
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
    TMVAssert((!U || U->iscm() || U->isrm()));
    TMVAssert((!V || V->iscm() || V->isrm()));
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(!U || !V || !SetUV || U->stor()==V->stor());

    if (D.size() > 0) {
#ifdef LAP
      LapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV);
#else 
      //NonLapSV_Decompose_From_Bidiagonal_DelayedU(U,D,E,V,SetUV);
      NonLapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV);
#endif
    }
  }

  //
  // Main SVD Drivers
  //
  
  template <class T> inline void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>* V, T& det, bool StoreU)
  {
#ifdef TIME
    cerr<<"Start SV_Decompose #1\n";
    timeval tp;
    gettimeofday(&tp,0);
    double t0 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
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
    TMVAssert(U.iscm() || U.isrm());

    // If M is much larger than N (technically M > 5/3 N), then it is quicker
    // to start by doing a QR decomposition and then do SVD on the square
    // R matrix.  Thus, the final U of the SVD is Q (from the QR decomp)
    // times U from R's SVD.
    if (M > 5*N/3) {
      if (StoreU) {
	Matrix<T,ColMajor> R(N,N);
	LowerTriMatrixViewOf(R).OffDiag().Zero();
	QR_Decompose(U,UpperTriMatrixViewOf(R),det);
#ifdef TIME
	gettimeofday(&tp,0);
	double t1 = tp.tv_sec + tp.tv_usec/1.e6;
	cerr<<"QR_Decompose: "<<t1-t0<<" seconds\n";
#endif
	SV_Decompose(R.View(),S,V,det,StoreU);
#ifdef TIME
	gettimeofday(&tp,0);
	double t2 = tp.tv_sec + tp.tv_usec/1.e6;
#endif
	// Now R is a Unitary Matrix U'.  Need to multiply U by U'
	U = U*R;
#ifdef TIME
	gettimeofday(&tp,0);
	double t3 = tp.tv_sec + tp.tv_usec/1.e6;
	cerr<<"Mult U*R: "<<t3-t2<<" seconds\n";
#endif
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
#ifdef TIME
      gettimeofday(&tp,0);
      double t1 = tp.tv_sec + tp.tv_usec/1.e6;
      cerr<<"Bidiagonalize: "<<t1-t0<<" seconds\n";
#endif
      // The determinant of B is just the product of the diagonal elements:
      if (det!=T(0)) det *= DiagMatrixViewOf(S).Det();
#ifdef TIME
      gettimeofday(&tp,0);
      double t1b = tp.tv_sec + tp.tv_usec/1.e6;
      cerr<<"Calc det: "<<t1b-t1<<" seconds\n";
#endif

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
#ifdef TIME
      gettimeofday(&tp,0);
      double t2 = tp.tv_sec + tp.tv_usec/1.e6;
      cerr<<"Make U,V: "<<t2-t1b<<" seconds\n";
#endif

      if (StoreU) SV_Decompose_From_Bidiagonal(&U,S,E.View(),V);
      else SV_Decompose_From_Bidiagonal((const MatrixView<T>*)(0),S,E.View(),V);

#ifdef TIME
      gettimeofday(&tp,0);
      double t3 = tp.tv_sec + tp.tv_usec/1.e6;
      cerr<<"Decompose From Bidiagonal: "<<t3-t2<<" seconds\n";
#endif
    }
#ifdef TIME
    gettimeofday(&tp,0);
    double tx = tp.tv_sec + tp.tv_usec/1.e6;
    cerr<<"Total time: "<<tx-t0<<" seconds\n";
#endif
  }

  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>& V, T& det, bool StoreU)
  { SV_Decompose(U,S,&V,det,StoreU); }

  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      T& det, bool StoreU)
  { SV_Decompose(U,S,(const MatrixView<T>*)(0),det,StoreU); }

#define InstFile "TMV_SVDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


