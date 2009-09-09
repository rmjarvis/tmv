
#include "TMV.h"
#include "TMV_Householder.h"
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define QR_BLOCKSIZE TMV_BLOCKSIZE
#else
#define QR_BLOCKSIZE 64
#endif

  //
  // QR Decompose
  //

  template <class T> inline void NonBlockQR_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, T& det)
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    // Decompose A into A = Q R 
    // where Q is unitary, and R is upper triangular
    // Q and R are stored in the same matrix (output of A), 
    // with the beta's for the Householder matrices returned in beta.
    const size_t M = A.colsize();
    const size_t N = A.rowsize();

    for(size_t j=0;j<N;++j) {
      // Apply the Householder Reflection for this column
      beta(j) = Householder_Reflect(A.SubMatrix(j,M,j,N),det);
    }
#ifdef XDEBUG
    Matrix<T> R = UpperTriMatrixViewOf(A);
    Matrix<T> Q = A;
    GetQFromQR(Q.View(),beta);
    Matrix<T> AA = Q*R;
    if (Norm(AA-A0) > 0.0001*Norm(Q)*Norm(R)) {
      cerr<<"NonBlockQR_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"Q*R = "<<Q*R<<endl;
      cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
      abort(); 
    }
#endif
  }

  template <class T> inline void RecursiveQR_Decompose(
      const MatrixView<T>& A, const UpperTriMatrixView<T>& Z, T& det,
      bool makeZ)
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Z.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(Z.iscm());
    // This is very similar to the BlockHouseholder_MakeZ function
    // in Householder.cpp.  The difference is the addition of the 
    // Householder_Reflects.
    // The makeZ parameter should be set to true if you want the Z
    // matrix to be correct on output.  If you don't need the Z matrix
    // after this call, setting makeZ to false will speed it up slightly.
    // (In either case, the diagonal of Z is correctly set to be 
    // beta.Conjugate().)
    const size_t M = A.colsize();
    const size_t N = A.rowsize();

    if (N==1) {
      T b = Householder_Reflect(A.col(0),det);
      *Z.ptr() = CONJ(b);
    } else if (N==2) {
      T b0 = Householder_Reflect(A,det);
      T* Zptr = Z.ptr();
      *Zptr = CONJ(b0);
      T b1 = Householder_Reflect(A.col(1,1,M),det);
      const int Zstepj = Z.stepj();
      *(Zptr+Zstepj+1) = CONJ(b1);

      if (makeZ) {
	T temp = A.col(0,2,M).Conjugate()*A.col(1,2,M);
	temp += CONJ(A(1,0));
	*(Zptr+Zstepj) = -CONJ(b0*b1)*temp;
      }
    } else {
      size_t j1 = (N+1)/2;
      MatrixView<T> A1 = A.Cols(0,j1);
      UpperTriMatrixView<T> Z1 = Z.SubTriMatrix(0,j1);
      RecursiveQR_Decompose(A1,Z1,det,true);

      BlockHouseholder_LDiv(A1,Z1,A.Cols(j1,N));

      MatrixView<T> A2 = A.SubMatrix(j1,M,j1,N);
      UpperTriMatrixView<T> Z2 = Z.SubTriMatrix(j1,N);
      RecursiveQR_Decompose(A2,Z2,det,makeZ);

      if (makeZ) {
	MatrixView<T> Z3 = Z.SubMatrix(0,j1,j1,N);
	Z3 = A1.Rows(j1,N).Adjoint() *
	  LowerTriMatrixViewOf(A.SubMatrix(j1,N,j1,N),UnitDiag);
	Z3 += A1.Rows(N,M).Adjoint() * A.SubMatrix(N,M,j1,N);
	Z3 = -Z1*Z3;
	Z3 *= Z2;
      }
    }
#ifdef XDEBUG
    Matrix<T> R = UpperTriMatrixViewOf(A);
    Matrix<T> Q = A;
    GetQFromQR(Q.View(),Z.diag().Conjugate());
    Matrix<T> AA = Q*R;
    if (Norm(AA-A0) > 0.0001*Norm(Q)*Norm(R)) {
      cerr<<"RecursiveQR_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"Z = "<<Z<<endl;
      cerr<<"beta = "<<Z.diag().Conjugate()<<endl;
      cerr<<"Q*R = "<<Q*R<<endl;
      cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
      Matrix<T> A2 = A0;
      Vector<T> beta2(Z.size());
      T det2(0);
      NonBlockQR_Decompose(A2.View(),beta2.View(),det2);
      cerr<<"NonBlock "<<A2<<endl;
      cerr<<"beta = "<<beta2<<endl;
      cerr<<"Norm(beta-beta2) = "<<Norm(beta-beta2)<<endl;
      cerr<<"Norm(A-A2) = "<<Norm(A-A2)<<endl;
      abort(); 
    }
#endif
  }

  template <class T> inline void BlockQR_Decompose(
      const MatrixView<T>& A, const VectorView<T>& beta, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif
    const size_t M = A.colsize();
    const size_t N = A.rowsize();

    UpperTriMatrix<T,NonUnitDiag,ColMajor> BaseZ(min(size_t(QR_BLOCKSIZE),N));
    for(size_t j1=0;j1<N;) {
      size_t j2 = min(N,j1+QR_BLOCKSIZE);
      MatrixView<T> A1 = A.SubMatrix(j1,M,j1,j2);
      UpperTriMatrixView<T> Z = BaseZ.SubTriMatrix(0,j2-j1);

      RecursiveQR_Decompose(A1,Z,det,j2<N);
      beta.SubVector(j1,j2) = Z.diag().Conjugate();

      /*
      for(size_t j=j1;j<j2;++j) {
	T b = Householder_Reflect(A.SubMatrix(j,M,j,j2),det);
	beta(j) = b;
	if (j2 < N)
	  BlockHouseholder_Augment(A.SubMatrix(j1,M,j1,j+1),
	      Z.SubTriMatrix(0,j+1),CONJ(b));
      }
      */

      if (j2 < N) 
	BlockHouseholder_LDiv(A1,Z,A.SubMatrix(j1,M,j2,N));
      j1 = j2;
    }
#ifdef XDEBUG
    Matrix<T> R = UpperTriMatrixViewOf(A);
    Matrix<T> Q = A;
    GetQFromQR(Q.View(),beta);
    Matrix<T> AA = Q*R;
    if (Norm(AA-A0) > 0.0001*Norm(Q)*Norm(R)) {
      cerr<<"BlockQR_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"Q*R = "<<Q*R<<endl;
      cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
      Matrix<T> A2 = A0;
      Vector<T> beta2(beta.size());
      T det2(0);
      NonBlockQR_Decompose(A2.View(),beta2.View(),det2);
      cerr<<"NonBlock "<<A2<<endl;
      cerr<<"beta = "<<beta2<<endl;
      abort(); 
    }
#endif
  }

  template <class T> inline void NonLapQR_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);

    if (A.rowsize() > QR_BLOCKSIZE)
      BlockQR_Decompose(A,beta,det);
    else {

      UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(A.rowsize());
      RecursiveQR_Decompose(A,Z.View(),det,false);
      beta = Z.diag().Conjugate();

      /*
      NonBlockQR_Decompose(A,beta,det);
      */
    }
  }

#ifdef LAP
  template <class T> inline void LapQR_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, T& det)
  { NonLapQR_Decompose(A,beta,det); }
  template <> inline void LapQR_Decompose(const MatrixView<double>& A,
      const VectorView<double>& beta, double& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    int m = A.colsize();
    int n = A.rowsize();
#ifndef LAPNOWORK
    int lwork = 2*n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
#endif
    if (A.isrm()) {
      int lda = A.stepi();
      LAPNAME(dgelqf) (LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("dgelqf");
#else
      LAP_Results(int(work[0]),m,n,lwork,"dgelqf");
#endif
    } else {
      int lda = A.stepj();
      LAPNAME(dgeqrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("dgeqrf");
#else
      LAP_Results(int(work[0]),m,n,lwork,"dgeqrf");
#endif
    }
    if (det) for(size_t i=0;i<beta.size();++i) if (beta(i) != 0.) det = -det;
  }
  template <> inline void LapQR_Decompose(
      const MatrixView<complex<double> >& A,
      const VectorView<complex<double> >& beta, complex<double>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    int m = A.colsize();
    int n = A.rowsize();
#ifndef LAPNOWORK
    int lwork = 2*n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
#endif
    if (A.isrm()) {
      int lda = A.stepi();
      LAPNAME(zgelqf) (LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("zgelqf");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"zgelqf");
#endif
    } else {
      int lda = A.stepj();
      LAPNAME(zgeqrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("zgeqrf");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"zgeqrf");
#endif
      beta.ConjugateSelf();
    }
    if (det!=double(0)) {
      for(size_t i=0;i<beta.size();++i)
	if (imag(beta(i)) != 0.) 
	  det *= -conj(beta(i)*beta(i))/norm(beta(i));
	else if (real(beta(i)) != 0.)
	  det = -det;
    }
  }
#ifndef NOFLOAT
  template <> inline void LapQR_Decompose(const MatrixView<float>& A,
      const VectorView<float>& beta, float& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    int m = A.colsize();
    int n = A.rowsize();
#ifndef LAPNOWORK
    int lwork = 2*n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
#endif
    if (A.isrm()) {
      int lda = A.stepi();
      LAPNAME(sgelqf) (LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("sgelqf");
#else
      LAP_Results(int(work[0]),m,n,lwork,"sgelqf");
#endif
    } else {
      int lda = A.stepj();
      LAPNAME(sgeqrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("sgeqrf");
#else
      LAP_Results(int(work[0]),m,n,lwork,"sgeqrf");
#endif
    }
    if (det) for(size_t i=0;i<beta.size();++i) if (beta(i) != 0.) det = -det;
  }
  template <> inline void LapQR_Decompose(
      const MatrixView<complex<float> >& A,
      const VectorView<complex<float> >& beta, complex<float>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    int m = A.colsize();
    int n = A.rowsize();
#ifndef LAPNOWORK
    int lwork = 2*n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
#endif
    if (A.isrm()) {
      int lda = A.stepi();
      LAPNAME(cgelqf) (LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("cgelqf");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"cgelqf");
#endif
    } else {
      int lda = A.stepj();
      LAPNAME(cgeqrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("cgeqrf");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"cgeqrf");
#endif
      beta.ConjugateSelf();
    }
    if (det!=float(0)) {
      for(size_t i=0;i<beta.size();++i)
	if (imag(beta(i)) != 0.) 
	  det *= -conj(beta(i)*beta(i))/norm(beta(i));
	else if (real(beta(i)) != 0.)
	  det = -det;
    }
  }
#endif
#endif
  template <class T> void QR_Decompose(
      const MatrixView<T>& A, const VectorView<T>& beta, T& det)
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    if (A.rowsize() > 0) {
#ifdef LAP
      LapQR_Decompose(A,beta,det);
#else
      NonLapQR_Decompose(A,beta,det);
#endif
    }
#ifdef XDEBUG
    Matrix<T> R = UpperTriMatrixViewOf(A);
    Matrix<T> Q = A;
    GetQFromQR(Q.View(),beta);
    Matrix<T> AA = Q*R;
    if (Norm(AA-A0) > 0.0001*Norm(Q)*Norm(R)) {
      cerr<<"BlockQR_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"Q*R = "<<Q*R<<endl;
      cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
      Matrix<T> A2 = A0;
      Vector<T> beta2(beta.size());
      T det2(0);
      NonBlockQR_Decompose(A2.View(),beta2.View(),det2);
      cerr<<"NonBlock "<<A2<<endl;
      cerr<<"beta = "<<beta2<<endl;
      abort(); 
    }
#endif
  }


  //
  // QR Decompose - Unpacked
  //

  template <class T> void QR_Decompose(
      const MatrixView<T>& Q, const UpperTriMatrixView<T>& R, T& det)
  {
    // Decompose A (input as Q) into A = Q R 
    // where Q is unitary and R is upper triangular

    const size_t N = Q.rowsize();
    TMVAssert(Q.colsize() >= N);
    TMVAssert(R.colsize() == N);
    TMVAssert(R.rowsize() == N);
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(R.ct() == NonConj);

    Vector<T> beta(N);
    QR_Decompose(Q,beta.View(),det);
    R = UpperTriMatrixViewOf(Q);
    GetQFromQR(Q.View(),beta.View());
  }

#define InstFile "TMV_QRDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


