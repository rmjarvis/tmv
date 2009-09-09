
#include "TMV.h"
#include "TMV_Householder.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

  bool StrictQRP = false;
  bool RecursiveQR = true;

#ifdef TMV_BLOCKSIZE
  const size_t QR_BLOCKSIZE = TMV_BLOCKSIZE;
  const size_t QRP_BLOCKSIZE = TMV_BLOCKSIZE;
#else
  const size_t QR_BLOCKSIZE = 64;
  const size_t QRP_BLOCKSIZE = 32;
#endif

#define APTR inplace ? A.NonConst().ptr() : new T[A.colsize()*A.rowsize()]
#define QRX istrans ? \
  inplace ? A.NonConst().Transpose() : \
  MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor) : \
  inplace ? A.NonConst().View() : \
  MatrixViewOf(Aptr,A.colsize(),A.rowsize(),ColMajor)

  template <class T> QRDiv<T>::QRDiv(const GenMatrix<T>& A,
      bool _inplace) :
    istrans(A.colsize()<A.rowsize()),
    inplace(_inplace), Aptr(APTR), QRx(QRX), beta(QRx.rowsize()),
    det(T(1)), donedet(false)
  {
    if (istrans) {
      if (inplace) { TMVAssert(A.Transpose() == QRx); }
      else QRx = A.Transpose();
    }
    else {
      if (inplace) { TMVAssert(A == QRx); }
      else QRx = A;
    }
    QR_Decompose(QRx,beta.View(),det);
  }

  template <class T> QRPDiv<T>::QRPDiv(const GenMatrix<T>& A, bool _inplace) :
    istrans(A.colsize()<A.rowsize()),
    inplace(_inplace), Aptr(APTR), QRx(QRX), beta(QRx.rowsize()),
    P(new size_t[beta.size()]),
    det(T(1)), donedet(false), N1(beta.size())
  {
    if (istrans) {
      if (inplace) { TMVAssert(A.Transpose() == QRx); }
      else QRx = A.Transpose();
    }
    else {
      if (inplace) { TMVAssert(A == QRx); }
      else QRx = A;
    }
    QRP_Decompose(QRx,beta.View(),P,det);
    while(N1>0 && QRx.diag()(N1-1)==T(0)) --N1;
  }
#undef QRX
#undef APTR

  template <class T> QRDiv<T>::~QRDiv<T>()
  { if (!inplace) delete[] Aptr; }

  template <class T> QRPDiv<T>::~QRPDiv<T>()
  { delete[] P; if (!inplace) delete[] Aptr; }

  template <class T> bool QRDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenMatrix<T>* mm = dynamic_cast<const GenMatrix<T>*>(&m);
    TMVAssert(mm);
#ifdef XDEBUG
    bool printmat = fout && m.colsize() < 100 && m.rowsize() < 100;
    if (printmat) {
      *fout << "M = "<<tmv::Type(*mm)<<"  ";
      *fout << (istrans ? mm->Transpose() : mm->View()) <<endl;
      *fout << "Q = "<<GetQ()<<endl;
      *fout << "R = "<<GetR()<<endl;
    }
#endif
    Matrix<T> qr = GetQ()*GetR();
    RealType(T) nm = Norm(qr- (istrans ? mm->Transpose() : mm->View()) );
    nm /= Norm(GetQ())*Norm(GetR());
#ifdef XDEBUG
    if (printmat) {
      *fout << "QR = "<<qr<<endl;
    }
#endif
    if (fout) {
      *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<"  "<<QRx.rowsize()*Epsilon<T>()<<endl;
    }
    Matrix<T> m2 = *mm;
    m2.DivideUsing(SVS);
    m2.SetDiv();
    return nm < m2.Condition()*m2.colsize()*Epsilon<T>();
  }

  template <class T> bool QRPDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenMatrix<T>* mm = dynamic_cast<const GenMatrix<T>*>(&m);
    TMVAssert(mm);
#ifdef XDEBUG
    bool printmat = fout && m.colsize() < 100 && m.rowsize() < 100;
    if (printmat) {
      *fout << "M = "<<tmv::Type(*mm)<<"  ";
      *fout << (istrans ? mm->Transpose() : mm->View()) <<endl;
      *fout << "Q = "<<GetQ()<<endl;
      *fout << "R = "<<GetR()<<endl;
      *fout << "P = ";
      for(size_t i=0;i<QRx.rowsize();i++) *fout<<P[i]<<" ";
      *fout<<endl;
    }
#endif
    Matrix<T> qr = GetQ()*GetR();
    qr.ReversePermuteCols(GetP());
    RealType(T) nm = Norm(qr- (istrans ? mm->Transpose() : mm->View()) );
    nm /= Norm(GetQ())*Norm(GetR());
#ifdef XDEBUG
    if (printmat) {
      *fout << "QRP = "<<qr<<endl;
    }
#endif
    if (fout) {
      *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<"  "<<QRx.rowsize()*Epsilon<T>()<<endl;
    }
    Matrix<T> m2 = *mm;
    m2.DivideUsing(SVS);
    m2.SetDiv();
    return nm < m2.Condition()*Epsilon<T>();
  }

  // 
  // Get Q from QR
  // 
 
  template <class T> void NonBlockGetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    UpperTriMatrixViewOf(Q).Zero();
    for(int i=N-1;i>=0;i--) {
      Householder_Unpack(Q.SubMatrix(i,M,i,N),beta(i));
    }
  }

  template <class T> void BlockGetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
#ifdef XDEBUG
    Matrix<T> Q0 = Q;
    Matrix<T> Q2 = Q;
    NonBlockGetQFromQR(Q2.View(),beta);
#endif
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    UpperTriMatrixViewOf(Q).Zero();
    UpperTriMatrix<T,NonUnitDiag,ColMajor> BaseZ(min(QR_BLOCKSIZE,N));
    for(size_t j2=N;j2>0;) {
      size_t j1 = j2 > QR_BLOCKSIZE ? j2-QR_BLOCKSIZE : 0;
      MatrixView<T> Y = Q.SubMatrix(j1,M,j1,j2);
      UpperTriMatrixView<T> Z = BaseZ.SubTriMatrix(0,Y.rowsize());
      BlockHouseholder_MakeZ(Y,Z,beta.SubVector(j1,j2));
      BlockHouseholder_Unpack(Y,Z,Q.SubMatrix(j1,M,j2,N));
      j2 = j1;
    }
#ifdef XDEBUG
    if (Norm(Q-Q2) > 0.001*max(RealType(T)(1),Norm(Q2))) {
      cerr<<"BlockGetQ: Q = "<<Type(Q)<<"  "<<Q0<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"-> Q = "<<Q<<endl;
      cerr<<"NonBlock Q = "<<Q2<<endl;
      abort(); 
    }
#endif
  }

  template <class T> inline void NonLapGetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    if (Q.rowsize() >= QR_BLOCKSIZE)
      BlockGetQFromQR(Q,beta);
    else
      NonBlockGetQFromQR(Q,beta);
  }
#ifdef LAP
  template <class T> inline void LapGetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& beta)
  { NonLapGetQFromQR(Q,beta); }
  template <> void LapGetQFromQR(const MatrixView<double>& Q,
      const GenVector<double>& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
    int lwork = n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    if (Q.isrm()) {
      int lda = Q.stepi();
      dorglq(&n,&m,&n,Q.ptr(),&lda,const_cast<double*>(beta.cptr()),
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dorglq");
    } else {
      int lda = Q.stepj();
      dorgqr(&m,&n,&n,Q.ptr(),&lda,const_cast<double*>(beta.cptr()),
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dorgqr");
    }
  }
  template <> void LapGetQFromQR(
      const MatrixView<complex<double> >& Q,
      const GenVector<complex<double> >& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
    int lwork = n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    if (Q.isrm()) {
      int lda = Q.stepi();
      zunglq(&n,&m,&n,LAP_Complex(Q.ptr()),&lda,LAP_Complex(beta.cptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"zunglq");
    } else {
      int lda = Q.stepj();
      Vector<complex<double> > conjbeta = beta.Conjugate();
      zungqr(&m,&n,&n,LAP_Complex(Q.ptr()),&lda,LAP_Complex(conjbeta.cptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"zungqr");
    }
  }
#ifndef NOFLOAT
  template <> void LapGetQFromQR(const MatrixView<float>& Q,
      const GenVector<float>& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
    int lwork = n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    if (Q.isrm()) {
      int lda = Q.stepi();
      sorglq(&n,&m,&n,Q.ptr(),&lda,const_cast<float*>(beta.cptr()),
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sorglq");
    } else {
      int lda = Q.stepj();
      sorgqr(&m,&n,&n,Q.ptr(),&lda,const_cast<float*>(beta.cptr()),
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sorgqr");
    }
  }
  template <> void LapGetQFromQR(
      const MatrixView<complex<float> >& Q,
      const GenVector<complex<float> >& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
    int lwork = n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    if (Q.isrm()) {
      int lda = Q.stepi();
      cunglq(&n,&m,&n,LAP_Complex(Q.ptr()),&lda,LAP_Complex(beta.cptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"cunglq");
    } else {
      int lda = Q.stepj();
      Vector<complex<float> > conjbeta = beta.Conjugate();
      cungqr(&m,&n,&n,LAP_Complex(Q.ptr()),&lda,LAP_Complex(conjbeta.cptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"cungqr");
    }
  }
#endif
#endif
  template <class T> void GetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& beta)
  {
    // Extract the Q matrix from a combined QRx matrix
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(Q.isrm() || Q.iscm());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
#ifdef LAP
    LapGetQFromQR(Q,beta);
#else
    NonLapGetQFromQR(Q,beta);
#endif
  }

  //
  // QR Decompose
  //

  template <class T> void NonBlockQR_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
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
  }

  template <class T> void RecursiveQR_Decompose(
      const MatrixView<T>& A, const UpperTriMatrixView<T>& Z, T& det,
      bool makeZ)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
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
      Z(0,0) = CONJ(b);
    } else if (N==2) {
      T b0 = Householder_Reflect(A,det);
      Z(0,0) = CONJ(b0);
      T b1 = Householder_Reflect(A.col(1,1,M),det);
      Z(1,1) = CONJ(b1);

      if (makeZ) {
	T temp = A.col(0,2,M).Conjugate()*A.col(1,2,M);
	temp += CONJ(A(1,0));
	Z(0,1) = -Z(0,0)*Z(1,1)*temp;
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
  }

  template <class T> void BlockQR_Decompose(
      const MatrixView<T>& A, const VectorView<T>& beta, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> A2 = A;
    Vector<T> beta2(A.rowsize());
    T det2 = det;
    NonBlockQR_Decompose(A2.View(),beta2.View(),det2);
#endif
    const size_t M = A.colsize();
    const size_t N = A.rowsize();

    UpperTriMatrix<T,NonUnitDiag,ColMajor> BaseZ(min(QR_BLOCKSIZE,N));
    for(size_t j1=0;j1<N;) {
      size_t j2 = min(N,j1+QR_BLOCKSIZE);
      MatrixView<T> A1 = A.SubMatrix(j1,M,j1,j2);
      UpperTriMatrixView<T> Z = BaseZ.SubTriMatrix(0,j2-j1);
      if (RecursiveQR) {
	RecursiveQR_Decompose(A1,Z,det,j2<N);
	beta.SubVector(j1,j2) = Z.diag().Conjugate();
      } else {
	for(size_t j=j1;j<j2;++j) {
	  T b = Householder_Reflect(A.SubMatrix(j,M,j,j2),det);
	  beta(j) = b;
	  if (j2 < N)
	    BlockHouseholder_Augment(A.SubMatrix(j1,M,j1,j+1),
		Z.SubTriMatrix(0,j+1),CONJ(b));
	}
      }
      if (j2 < N) 
	BlockHouseholder_LDiv(A1,Z,A.SubMatrix(j1,M,j2,N));
      j1 = j2;
    }
#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*max(RealType(T)(1),Norm(A2))) {
      cerr<<"BlockQR_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"beta = "<<beta<<endl;
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
    TMVAssert(beta.ct() == NonConj);

    if (A.rowsize() > QR_BLOCKSIZE)
      BlockQR_Decompose(A,beta,det);
    else if (RecursiveQR) {
      UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(A.rowsize());
      RecursiveQR_Decompose(A,Z.View(),det,false);
      beta = Z.diag().Conjugate();
    } else NonBlockQR_Decompose(A,beta,det);
  }

#ifdef LAP
  template <class T> inline void LapQR_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, T& det)
  { NonLapQR_Decompose(A,beta,det); }
  template <> void LapQR_Decompose(const MatrixView<double>& A,
      const VectorView<double>& beta, double& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lwork = 2*n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    if (A.isrm()) {
      int lda = A.stepi();
      dgelqf(&n,&m,A.ptr(),&lda,beta.ptr(),work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dgelqf");
    } else {
      int lda = A.stepj();
      dgeqrf(&m,&n,A.ptr(),&lda,beta.ptr(),work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dgeqrf");
    }
    if (det) for(size_t i=0;i<beta.size();++i) if (beta(i) != 0.) det = -det;
  }
  template <> void LapQR_Decompose(
      const MatrixView<complex<double> >& A,
      const VectorView<complex<double> >& beta, complex<double>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lwork = 2*n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    if (A.isrm()) {
      int lda = A.stepi();
      zgelqf(&n,&m,LAP_Complex(A.ptr()),&lda,LAP_Complex(beta.ptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"zgelqf");
    } else {
      int lda = A.stepj();
      zgeqrf(&m,&n,LAP_Complex(A.ptr()),&lda,LAP_Complex(beta.ptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"zgeqrf");
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
  template <> void LapQR_Decompose(const MatrixView<float>& A,
      const VectorView<float>& beta, float& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lwork = 2*n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    if (A.isrm()) {
      int lda = A.stepi();
      sgelqf(&n,&m,A.ptr(),&lda,beta.ptr(),work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sgelqf");
    } else {
      int lda = A.stepj();
      sgeqrf(&m,&n,A.ptr(),&lda,beta.ptr(),work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sgeqrf");
    }
    if (det) for(size_t i=0;i<beta.size();++i) if (beta(i) != 0.) det = -det;
  }
  template <> void LapQR_Decompose(
      const MatrixView<complex<float> >& A,
      const VectorView<complex<float> >& beta, complex<float>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lwork = 2*n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    if (A.isrm()) {
      int lda = A.stepi();
      cgelqf(&n,&m,LAP_Complex(A.ptr()),&lda,LAP_Complex(beta.ptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"cgelqf");
    } else {
      int lda = A.stepj();
      cgeqrf(&m,&n,LAP_Complex(A.ptr()),&lda,LAP_Complex(beta.ptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"cgeqrf");
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
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    if (A.rowsize() > 0)
#ifdef LAP
      if (beta.step()==1 && (A.isrm() || A.iscm()))
	LapQR_Decompose(A,beta,det);
      else 
#endif
	NonLapQR_Decompose(A,beta,det);
  }

  //
  // QR Update
  //

  template <class T> void NonBlockQR_Update(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);
    // Given that A0 = Q0 R0
    // Find R1, so that [ A0 ] = Q1 R1
    //                  [ A  ] 
    // Input R is R0, output is R1

    const size_t N = A.rowsize();

    T* Rdiag = R.ptr();
    const size_t ds = R.stepi()+R.stepj();
    T det(0);

    for(size_t j=0;j<N;++j,Rdiag+=ds) {
      // Apply the Householder Reflection for this column
      const VectorView<T> v = A.col(j);
      T beta = Householder_Reflect(*Rdiag,v,det);
      if (beta != T(0))
	Householder_LMult(v,beta,R.row(j,j+1,N),A.Cols(j+1,N));
    }
  }

  template <class T> void RecursiveQR_Update(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A,
      const UpperTriMatrixView<T>& Z, bool makeZ)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

    const size_t N = A.rowsize();
    T det(0);

    if (N==1) {
      T b = Householder_Reflect(R(0,0),A.col(0),det);
      Z(0,0) = CONJ(b);
    } else if (N==2) {
      T b0 = Householder_Reflect(R(0,0),A.col(0),det);
      if (b0 != T(0)) {
	T temp = b0*(A.col(0).Conjugate()*A.col(1) + R(0,1));
	R(0,1) -= temp;
	A.col(1) -= temp * A.col(0);
      }
      Z(0,0) = CONJ(b0);
      T b1 = Householder_Reflect(R(1,1),A.col(1),det);
      Z(1,1) = CONJ(b1);

      if (makeZ) {
	T temp = A.col(0).Conjugate()*A.col(1);
	Z(0,1) = -Z(0,0)*Z(1,1)*temp;
      }
    } else {
      size_t j1 = N/2;

      UpperTriMatrixView<T> R1 = R.SubTriMatrix(0,j1);
      MatrixView<T> Rx = R.SubMatrix(0,j1,j1,N);
      UpperTriMatrixView<T> R2 = R.SubTriMatrix(j1,N);

      MatrixView<T> A1 = A.Cols(0,j1);
      MatrixView<T> A2 = A.Cols(j1,N);

      UpperTriMatrixView<T> Z1 = Z.SubTriMatrix(0,j1);
      MatrixView<T> Zx = Z.SubMatrix(0,j1,j1,N);
      UpperTriMatrixView<T> Z2 = Z.SubTriMatrix(j1,N);

      RecursiveQR_Update(R1,A1,Z1,true);

      // Zx is a temporary here - it happens to be the right shape.
      Zx = A1.Adjoint() * A2; 
      Zx += Rx;
      Zx = Z1.Adjoint()*Zx;
      Rx -= Zx;
      A2 -= A1 * Zx;

      RecursiveQR_Update(R2,A2,Z2,makeZ);

      if (makeZ) {
	Zx = A1.Adjoint() * A2; 
	Zx = -Z1*Zx;
	Zx *= Z2;
      }
    }
  }

  template <class T> void BlockQR_Update(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

    const size_t N = A.rowsize();
    const size_t ds = R.stepi() + R.stepj();
    T* Rdiag = R.ptr();

    UpperTriMatrix<T,NonUnitDiag,ColMajor> BaseZ(min(QR_BLOCKSIZE,N));
    for(size_t j1=0;j1<N;) {
      size_t j2 = min(N,j1+QR_BLOCKSIZE);
      MatrixView<T> A1 = A.Cols(j1,j2);
      UpperTriMatrixView<T> R1 = R.SubTriMatrix(j1,j2);
      UpperTriMatrixView<T> Z = BaseZ.SubTriMatrix(0,j2-j1);
      if (RecursiveQR) {
	RecursiveQR_Update(R1,A1,Z,j2<N);
      } else {
	T det(0);
	for(size_t j=j1,jj=0;j<j2;++j,++jj,Rdiag+=ds) {
	  const VectorView<T> v = A.col(j);
	  T b = Householder_Reflect(*Rdiag,v,det);
	  if (b != T(0)) {
	    Householder_LMult(v,b,R.row(j,j+1,j2),A.Cols(j+1,j2));
	    if (j2 < N) {
	      if (jj > 0) {
		VectorView<T> z = Z.col(jj,0,jj);
		z = A.Cols(j1,j).Adjoint() * v;
		z = -b * Z.SubTriMatrix(0,jj) * z;
	      }
	      Z(jj,jj) = CONJ(b);
	    }
	  } else {
	    Z.col(jj,0,jj+1).Zero();
	  }
	}
      }
      if (j2 < N) {
	Matrix<T,ColMajor> ZtYtm = A.Cols(j1,j2).Adjoint() * A.Cols(j2,N);
	ZtYtm += R.SubMatrix(j1,j2,j2,N);
	ZtYtm = Z.Adjoint() * ZtYtm;
	R.SubMatrix(j1,j2,j2,N) -= ZtYtm;
        A.Cols(j2,N) -= A.Cols(j1,j2) * ZtYtm;
      }
      j1 = j2;
    }
  }

  template <class T> void QR_Update(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

#ifdef XDEBUG
    UpperTriMatrix<T> R0 = R;
    UpperTriMatrix<T> R2 = R;
    Matrix<T> A0 = A;
    Matrix<T> A2 = A;
    NonBlockQR_Update(R2.View(),A2.View());
#endif
    if (A.rowsize() > 0)
      if (A.rowsize() > QR_BLOCKSIZE)
	BlockQR_Update(R,A);
      else if (RecursiveQR) {
	UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(A.rowsize());
	RecursiveQR_Update(R,A,Z.View(),false);
      }
      else 
	NonBlockQR_Update(R,A);
#ifdef XDEBUG
    if (Norm(R2-R) > 1.e-5*Norm(R)) {
      cerr<<"QR_Update\n";
      cerr<<"R0 = "<<Type(R)<<"  "<<R0<<endl;
      cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"R -> "<<R<<endl;
      cerr<<"NonBlock R -> "<<R2<<endl;
      abort();
    }
#endif
  }

  //
  // QR Downdate
  //

  template <class T> bool NonBlockQR_Downdate(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);
    // Given that A0 = Q0 R0
    // Given that [ A0 ] = Q1 R1
    //            [ A  ] 
    // Find R0 so that A0 = Q0 R0
    // Input R is R1, output is R0

    const size_t N = A.rowsize();

    T* Rdiag = R.ptr();
    const size_t ds = R.stepi()+R.stepj();
    T det(0);

    for(size_t j=0;j<N;++j,Rdiag+=ds) {
      // Apply the Householder Reflection for this column
      const VectorView<T> v = A.col(j);
      T beta;
      if (!(Householder_UnReflect(*Rdiag,v,beta,det))) return false;
      TMVAssert(beta != T(1));
      VectorView<T> m0 = R.row(j,j+1,N);
      MatrixView<T> mx = A.Cols(j+1,N);

      // m0' = m0 - beta m0 - beta vtmx
      // m0 = (m0' + beta btmx)/(1-beta)
      Vector<T> bvtmx = beta*v.Conjugate()*mx;
      m0 += bvtmx;
      m0 /= T(1)-beta;

      // mx' = mx - beta v (m0 + vtmx)
      bvtmx += beta*m0;
      mx -= v ^ bvtmx;
    }
    return true;
  }

  template <class T> bool RecursiveQR_Downdate(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A,
      const UpperTriMatrixView<T>& Z, bool makeZ)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

    const size_t N = A.rowsize();
    T det(0);

    if (N==1) {
      T b;
      if (!(Householder_UnReflect(R(0,0),A.col(0),b,det))) return false;
      Z(0,0) = CONJ(b);
    } else if (N==2) {
      T b0;
      if (!(Householder_UnReflect(R(0,0),A.col(0),b0,det))) return false;
      Z(0,0) = CONJ(b0);
      if (b0 != T(0)) {
	TMVAssert(b0 != T(1));
	T vtmx = A.col(0).Conjugate() * A.col(1);
	R(0,1) = (R(0,1) + b0*vtmx)/(T(1)-b0);
	A.col(1) -= b0*(vtmx + R(0,1)) * A.col(0);
      }

      T b1;
      if (!(Householder_UnReflect(R(1,1),A.col(1),b1,det))) return false;
      Z(1,1) = CONJ(b1);

      if (makeZ) {
	T vtmx = A.col(0).Conjugate() * A.col(1);
	Z(0,1) = -CONJ(b0*b1)*vtmx;
      }
    } else {
      size_t j1 = N/2;

      UpperTriMatrixView<T> R1 = R.SubTriMatrix(0,j1);
      MatrixView<T> Rx = R.SubMatrix(0,j1,j1,N);
      UpperTriMatrixView<T> R2 = R.SubTriMatrix(j1,N);

      MatrixView<T> A1 = A.Cols(0,j1);
      MatrixView<T> A2 = A.Cols(j1,N);

      UpperTriMatrixView<T> Z1 = Z.SubTriMatrix(0,j1);
      MatrixView<T> Zx = Z.SubMatrix(0,j1,j1,N);
      UpperTriMatrixView<T> Z2 = Z.SubTriMatrix(j1,N);

      if (!(RecursiveQR_Downdate(R1,A1,Z1,true))) return false;

      Zx = A1.Adjoint() * A2;
      Zx = Z1.Adjoint() * Zx;

      Rx += Zx;
      LowerTriMatrix<T> ImZt = T(1)-Z1.Adjoint();
      Rx /= ImZt;

      Zx += Z1.Adjoint() * Rx;
      A2 -= A1 * Zx;

      if (!(RecursiveQR_Downdate(R2,A2,Z2,makeZ))) return false;

      if (makeZ) {
	Zx = A1.Adjoint() * A2; 
	Zx = -Z1*Zx;
	Zx *= Z2;
      }
    }
    return true;
  }

  template <class T> bool BlockQR_Downdate(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

    const size_t N = A.rowsize();
    const size_t ds = R.stepi() + R.stepj();
    T* Rdiag = R.ptr();

    UpperTriMatrix<T,NonUnitDiag,ColMajor> BaseZ(min(QR_BLOCKSIZE,N));
    for(size_t j1=0;j1<N;) {
      size_t j2 = min(N,j1+QR_BLOCKSIZE);
      MatrixView<T> A1 = A.Cols(j1,j2);
      UpperTriMatrixView<T> R1 = R.SubTriMatrix(j1,j2);
      UpperTriMatrixView<T> Z = BaseZ.SubTriMatrix(0,j2-j1);
      if (RecursiveQR) {
	if (!(RecursiveQR_Downdate(R1,A1,Z,j2<N))) return false;
      } else {
	T det(0);
	for(size_t j=j1,jj=0;j<j2;++j,++jj,Rdiag+=ds) {
	  const VectorView<T> v = A.col(j);
	  T b;
	  if (!(Householder_UnReflect(*Rdiag,v,b,det))) return false;
	  if (b != T(0)) {
	    TMVAssert(b != T(1));
	    VectorView<T> m0 = R.row(j,j+1,j2);
	    MatrixView<T> mx = A.Cols(j+1,j2);

	    // m0' = m0 - beta m0 - beta vtmx
	    // m0 = (m0' + beta btmx)/(1-beta)
	    Vector<T> bvtmx = b*v.Conjugate()*mx;
	    m0 += bvtmx;
	    m0 /= T(1)-b;

	    // mx' = mx - beta v (m0 + vtmx)
	    bvtmx += b*m0;
	    mx -= v ^ bvtmx;

	    if (j2 < N) {
	      if (jj > 0) {
		VectorView<T> z = Z.col(jj,0,jj);
		z = A.Cols(j1,j).Adjoint() * v;
		z = -b * Z.SubTriMatrix(0,jj) * z;
	      }
	      Z(jj,jj) = CONJ(b);
	    }
	  } else {
	    Z.col(jj,0,jj+1).Zero();
	  }
	}
      }
      if (j2 < N) {
	// m0' = m0 - Zt(Ytmx+m0)
	// m0' + ZtYtmx = (I-Zt) m0;
	MatrixView<T> m0 = R.SubMatrix(j1,j2,j2,N);
	MatrixView<T> Y = A.Cols(j1,j2);
	MatrixView<T> mx = A.Cols(j2,N);

	Matrix<T,ColMajor> ZtYtm = Y.Adjoint() * mx;
	ZtYtm = Z.Adjoint() * ZtYtm;

	m0 += ZtYtm;
	LowerTriMatrix<T> ImZt = T(1)-Z.Adjoint();
	m0 /= ImZt;

	ZtYtm += Z.Adjoint() * m0;
        mx -= Y * ZtYtm;
      }
      j1 = j2;
    }
    return true;
  }

  template <class T> bool QR_Downdate(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A)
  {
    TMVAssert(A.rowsize() == R.size());
    TMVAssert(R.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(R.dt() == NonUnitDiag);

#ifdef XDEBUG
    UpperTriMatrix<T> R0 = R;
    UpperTriMatrix<T> R2 = R;
    Matrix<T> A0 = A;
    Matrix<T> A2 = A;
    bool ret2 = NonBlockQR_Downdate(R2.View(),A2.View());
#endif

    bool ret;
    if (A.rowsize() > 0) {
      if (A.rowsize() > QR_BLOCKSIZE)
	return BlockQR_Downdate(R,A);
      else if (RecursiveQR) {
	UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(A.rowsize());
	return RecursiveQR_Downdate(R,A,Z.View(),false);
      }
      else 
        ret = NonBlockQR_Downdate(R,A);
    }
    else ret = true;

#ifdef XDEBUG
    if (ret && (!ret2 || Norm(R2-R) > 1.e-5*Norm(R)) ) {
      cerr<<"QR_Downdate\n";
      cerr<<"Succeed? = "<<ret<<endl;
      cerr<<"R0 = "<<Type(R)<<"  "<<R0<<endl;
      cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"R -> "<<R<<endl;
      cerr<<"NonBlock R -> "<<R2<<endl;
      abort();
    }
#endif

    return ret;
  }

  //
  // QRP Decompose
  //

  template <class T> void NonBlockQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, size_t* P, T& det)
  {
    // Decompose A into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (output of A), 
    // with the beta's for the Householder matrices returned in beta.
    
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(beta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    const size_t M = A.colsize();
    const size_t N = A.rowsize();
    const int Astepj = A.stepj();

    // Keep track of the norm of each column
    // When considering column j, these are actually just the norm
    // of each column from j:M, not 0:M.
    Vector<RealType(T)> colnormsq(N);
    for(size_t j=0;j<N;++j) colnormsq(j) = NormSq(A.col(j));
    RealType(T) anormsq = colnormsq.SumElements();
    RealType(T) thresh = RealType(T)(N) * SQR(Epsilon<T>()) * anormsq;
    //cerr<<"thresh = "<<thresh<<endl;
    // Set to 0 any diag element whose norm is < epsilon * |A|
    RealType(T) recalcthresh(0);
    // recalcthresh is the threshold for recalculating the norm to account for
    // rounding errors in the subtractions which keep track of it..
    // The is set to sqrt(Epsilon) * the largest normsq whenever we 
    // recalculate the norms. 

    for(size_t j=0;j<N;++j) {
      //cerr<<"j = "<<j<<" colnormsq = "<<colnormsq(j)<<endl;

      if (StrictQRP || j==0 || colnormsq(j) < recalcthresh) {
	// Find the column with the largest norm
	size_t jpiv;
	RealType(T) maxnormsq = colnormsq.SubVector(j,N).MaxElement(&jpiv);
	if (j==0) recalcthresh = 4*SqrtEpsilon<T>() * maxnormsq;
	// Note: jpiv is relative to the SubVector(j,N)

	//cerr<<"jpiv = "<<jpiv+j<<", maxnormsq = "<<maxnormsq<<endl;

	// If the largest colnormsq is lower than the recalulation threshold,
	// then recalc all colnormsq's, and redetermine max.
	if (maxnormsq < recalcthresh) {
	  for(size_t k=j;k<N;++k) colnormsq(k) = NormSq(A.col(k,j,M));
	  maxnormsq = colnormsq.SubVector(j,N).MaxElement(&jpiv);
	  recalcthresh = 4*SqrtEpsilon<T>() * maxnormsq;
	  if (recalcthresh < thresh) recalcthresh = thresh;
	}

	// If maxnormsq = 0 (technically < thresh to account for rounding)
	// then the rest of the R matrix is 0, and the Householder matrices 
	// are identities (indicated by 0's in the Q part of the matrix).
	if (maxnormsq < thresh) {
	  //cerr<<"Do Zero\n";
	  A.SubMatrix(j,M,j,N).Zero();
	  // Already essentially zero - make it exact
	  beta.SubVector(j,N).Zero();
	  // Set the Householder matrices for these to identities
	  for(;j<N;j++) P[j] = j;
	  break;
	} else {

	  // Swap the column with the largest norm into the current column
	  if (jpiv != 0) {
	    // Add j to get real index
	    jpiv += j;
	    colnormsq.Swap(j,jpiv);
	    A.SwapCols(j,jpiv);
	    det = -det;
	    P[j] = jpiv;
	  } else {
	    P[j] = j;
	  }
	}
      } else P[j] = j;

      // Apply the Householder Reflection for this column
      beta(j) = Householder_Reflect(A.SubMatrix(j,M,j,N),det);

      // And update the norms for use with the next column
      const T* Ajk = A.row(j,j+1,N).cptr();
      for(size_t k=j+1;k<N;++k,Ajk+=Astepj) {
	colnormsq(k) -= NORM(*Ajk);
      }
    }
  }

  template <class T> void StrictBlockQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, size_t* P, T& det)
  {
    // Decompose A (input as A) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (A), with the beta's for
    // the Householder matrices returned in beta.
    
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(beta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    const size_t M = A.colsize();
    const size_t N = A.rowsize();
    const int Astepj = A.stepj();

    Vector<RealType(T)> colnormsq(N);
    for(size_t j=0;j<N;++j) colnormsq(j) = NormSq(A.col(j));
    RealType(T) anormsq = colnormsq.SumElements();
    RealType(T) thresh = RealType(T)(N) * SQR(Epsilon<T>()) * anormsq;
    RealType(T) recalcthresh(0);

    Matrix<T,RowMajor> ZYtA(min(QRP_BLOCKSIZE,M),N);
    // We keep track of the product ZYtA(j1:M,j1:N) [ stored as ZYtA ]
    // since this is the product that we need.  We update this one 
    // row at a time.
    
    for(size_t j1=0;j1<N;) {
      size_t j2 = min(N,j1+QRP_BLOCKSIZE);

      for(size_t j=j1,jmj1=0; j<j2; ++j,++jmj1) {
	size_t jpiv;
	RealType(T) maxnormsq = colnormsq.SubVector(j,N).MaxElement(&jpiv);
	if (recalcthresh == RealType(T)(0)) 
	  recalcthresh = 4*SqrtEpsilon<T>() * maxnormsq;

	if (maxnormsq < recalcthresh) {
	  for(size_t k=j;k<N;++k) colnormsq(k) = NormSq(A.col(k,j,M));
	  recalcthresh = RealType(T)(0);
	  j2 = j;
	} else if (maxnormsq < thresh) {
	  if (j==j1) {
	    // If first in set, then just zero the rest out and 
	    // indicate that we are done.
	    A.SubMatrix(j,M,j,N).Zero();
	    beta.SubVector(j,N).Zero();
	    for(;j<N;j++) P[j] = j;
	    j2 = N; 
	  } else {
	    // Otherwise, do the block Householder transforms on the 
	    // previous columns first.  (The next time through the block
	    // loop should still result in maxnormsq < thresh.)
	    j2 = j; 
	  }
	} else {

	  // Pivot
	  if (jpiv != 0) {
	    jpiv += j;
	    colnormsq.Swap(j,jpiv);
	    A.SwapCols(j,jpiv);
	    ZYtA.Rows(0,jmj1).SwapCols(j,jpiv);
	    det = -det;
	    P[j] = jpiv;
	  } else {
	    P[j] = j;
	  }

	  // Update the pivot column with Block Householder so far:
	  // A(j1:M,j) -= Y Z Yt A(j1:M,j)
	  // A(j1:j,j) has already been updated, so we only need to do A(j:M,j)
	  // A(j:M,j) -= Y(j:M,0:j) (ZYtA)(0:j,j)
	  A.col(j,j,M) -= A.SubMatrix(j,M,j1,j) * ZYtA.col(j,0,jmj1);

	  // Find Householder matrix for this column
	  beta(j) = Householder_Reflect(A.col(j,j,M),det);

	  // Update ZYtA:
	  if (beta(j) != T(0)) {
	    // (I-beta v vt)(I-Y Z Yt) A
	    // = I - Y (ZYtA) - v (beta vt A) + v (beta vt Y (ZYtA))
	    // The augmented Y now includes v in the j column, 
	    // so the augmented ZYtA now has to include in the j row:
	    // beta (vt A - vt Y ZYtA)
	    VectorView<T> vt = A.col(j,j+1,M).Conjugate();
	    // Remember, this doesn't include the implicit 1 at the top of v.
	    ZYtA.row(jmj1,j1,j+1).Zero();
	    ZYtA.row(jmj1,j+1,N) = vt * A.SubMatrix(j+1,M,j+1,N);
	    ZYtA.row(jmj1,j+1,N) += A.row(j,j+1,N);
	    Vector<T> vtY = vt * A.SubMatrix(j+1,M,j1,j);
	    vtY += A.row(j,j1,j);
	    ZYtA.row(jmj1,j1,N) -= vtY * ZYtA.SubMatrix(0,jmj1,j1,N);
	    ZYtA.row(jmj1,j1,N) *= beta(j);
	  } else ZYtA.row(jmj1,j1,N).Zero();

	  // Update row j of the rest of the matrix:
	  // A(j,j+1:N) -= (Y ZYtA)(j,j+1:N) = Y(j,j1:j+1) ZYtA(j1:j+1,j+1:N)
	  VectorView<T> Arowj = A.row(j,j+1,N);
	  Arowj -= A.row(j,j1,j)*ZYtA.SubMatrix(0,jmj1,j+1,N);
	  Arowj -= ZYtA.row(jmj1,j+1,N);

	  // Update the colnormsq values
	  const T* Ajk = Arowj.cptr();
	  for(size_t k=j+1;k<N;++k,Ajk+=Astepj) colnormsq(k) -= tmv::NORM(*Ajk);
	}
      }
      // Do the Block Householder update of the rest of the matrix:
      // A(j2:M,j2:N) -= Y(j2:M,j1:j2) ZYtA(j1:j2,j1:N)
      A.SubMatrix(j2,M,j2,N) -= A.SubMatrix(j2,M,j1,j2) * 
	ZYtA.SubMatrix(0,j2-j1,j2,N);
      j1 = j2;
    }

#ifdef XDEBUG
    Matrix<T> Q = A;
    GetQFromQR(Q.View(),beta);
    Matrix<T> AA = Q*UpperTriMatrixViewOf(A);
    AA.ReversePermuteCols(P);
    if (Norm(AA-A0) > 0.001*max(RealType(T)(1),Norm(A0))) {
      cerr<<"StrictBlockQRP_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"P = ";
      for(size_t i=0;i<A.rowsize();i++) cerr<<P[i]<<" ";
      cerr<<endl;
      cerr<<"QRP = "<<AA<<endl;
      abort(); 
    }
#endif
  }

  template <class T> void MoveLowColsToEnd(
      Vector<RealType(T)>& colnormsq, RealType(T) thresh,
      size_t j1, size_t& j2, size_t& j3, const MatrixView<T>& A, size_t* P)
  {
    // Move all columns of A (whose norms are in colnormsq) with norms
    // less than thresh to the end.  j1 is the first column we need to 
    // look at.  j2 is one past the last column we need to look at.
    // j3 is one past the last good column overall.
    // On output j2 and j3 are updated lower if necessary.
    
    --j3; // temporarily, j3 is the last good column.
    while (j3 > j1 && colnormsq(j3) < thresh) --j3;
    if (j3==j1) {
      if (colnormsq(j1) < thresh) j2 = j1;
      else { j2 = j3 = j1+1; P[j1] = j1; }
      return;
    }
    if (j3 < j2) j2 = j3+1;
    for(size_t i=j1;i<j2;++i) {
      if (colnormsq(i) < thresh) {
	//cerr<<"  "<<i<<" is low ("<<colnormsq(i)<<") - swap with "<<j3<<endl;
	A.SwapCols(i,j3);
	colnormsq.Swap(i,j3);
	P[i] = j3;
	while (colnormsq(--j3) < thresh);
	if (j3 < j2) j2 = j3+1;
      } else P[i] = i;
    }
    ++j3; // j3 back to being first bad column.
  }

#ifdef XDEBUG
  void CheckIndex(const Vector<double>& index, const size_t* P, size_t j1)
  {
    const size_t N = index.size();
    Vector<double> index2(N);
    for(size_t k=0;k<N;k++) index2(k) = double(k);
    for(size_t k=0;k<j1;k++) index2.Swap(k,P[k]);
    if (Norm(index-index2) > 0.01) {
      cerr<<"index = "<<index<<endl;
      cerr<<"index2 = "<<index2<<endl;
      cerr<<"norm(diff) = "<<Norm(index-index2)<<endl;
      abort();
    }
  }
#endif

  template <class T> void LooseBlockQRP_Decompose(
      const MatrixView<T>& A, const VectorView<T>& beta, size_t* P, T& det)
  {
    // Decompose A (input as A) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (A), with the beta's for
    // the Householder matrices returned in beta.
    //
    // This loose version doesn't sort the diagonal of R exactly.
    // It only sorts them enough to make sure the 0's fall at the end.
    
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(beta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    const size_t M = A.colsize();
    const size_t N = A.rowsize();
    const int Astepj = A.stepj();

    Vector<RealType(T)> colnormsq(N);
    for(size_t j=0;j<N;++j) colnormsq(j) = NormSq(A.col(j));
    //cerr<<"colnormsq = "<<colnormsq<<endl;
    RealType(T) anormsq = colnormsq.SumElements();
    //cerr<<"anormsq = "<<anormsq<<endl;
    //cerr<<"eps = "<<Epsilon<T>()<<endl;
    RealType(T) thresh = RealType(T)(N) * SQR(Epsilon<T>()) * anormsq;
    //cerr<<"thresh = "<<thresh<<endl;

#ifdef XDEBUG
    Vector<double> index(N);
    for(size_t k=0;k<N;k++) index(k)=double(k);
#endif

    for (size_t j1 = 0; j1 < N;) {
      //cerr<<"Start outer loop: j1 = "<<j1<<endl;
      // Do as many columns as possible such that none have to have 
      // their norms recalculated.
      size_t j3=N; // j3 will be how many we have done in this loop
      // Invariant: all columns from j3..N are known to have norms that
      // need to be recalculated.
      // The recalculation is done at the end of the loop.

      size_t jpiv0;
      RealType(T) maxnormsq = colnormsq.SubVector(j1,N).MaxElement(&jpiv0);
      //cerr<<"maxnormsq = "<<maxnormsq<<endl;

      if (maxnormsq < thresh) {
	//cerr<<"Do Zero\n";
	// Zero the rest out and we are done.
	A.SubMatrix(j1,M,j1,N).Zero();
	beta.SubVector(j1,N).Zero();
	for(;j1<N;j1++) P[j1] = j1;
	break;
      } 

      // Move max column to the front:
      if (jpiv0 != 0) {
	jpiv0 += j1;
	//cerr<<"swap j1 with jpiv0 = "<<jpiv0<<endl;
	A.SwapCols(j1,jpiv0);
	colnormsq.Swap(j1,jpiv0);
	P[j1] = jpiv0;
#ifdef XDEBUG
	index.Swap(j1,jpiv0);
#endif
      } else P[j1] = j1;
#ifdef XDEBUG
      CheckIndex(index,P,j1+1);
#endif

      RealType(T) recalcthresh = RealType(T)(N)*SqrtEpsilon<T>()*maxnormsq;
      if (recalcthresh < thresh) recalcthresh = thresh;
      //cerr<<"recalcthresh = "<<recalcthresh<<endl;

      TMVAssert(j1<j3);
      size_t j1x = j1+1; 
      // The first pass through, we don't want to include j1 in the 
      // MoveLowColsToEnd call.

      // Work on this one block at a time:
      while (j1 < j3) {
	//cerr<<"Start inner loop: j1 = "<<j1<<", j3 = "<<j3<<endl;
	size_t j2 = min(j3,j1+QRP_BLOCKSIZE);
	//cerr<<"j2 = "<<j2<<endl;
	TMVAssert(j1 < j2);
	MoveLowColsToEnd(colnormsq,recalcthresh,j1x,j2,j3,A,P);
	//cerr<<"After MoveLowColsToEnd: \n";
	//cerr<<"j2, j3 = "<<j2<<','<<j3<<endl;
#ifdef XDEBUG
	for(size_t k=j1x;k<j2;k++) index.Swap(k,P[k]);
	CheckIndex(index,P,j2);
#endif

	size_t origj2 = j2;
	UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(j2-j1);
	//cerr<<"Block "<<j1<<".."<<j2<<endl;

	//cerr<<"colnormsq = "<<colnormsq.SubVector(j1,j2)<<endl;

	for(size_t j=j1; j<j2; ++j) {

	  if (colnormsq(j) < recalcthresh) {
#ifdef XDEBUG
	    CheckIndex(index,P,origj2);
#endif
	    --j2;
	    //cerr<<" "<<j<<" is low ("<<colnormsq(j)<<") - swap with "<<j2<<endl;
	    //cerr<<" Current Ps are "<<P[j]<<"  "<<P[j2]<<endl;
	    if (j==j2) break;
	    A.SwapCols(j,j2);
	    colnormsq.Swap(j,j2);
#ifdef XDEBUG
	    index.Swap(j,j2);
#endif
	    if (P[j2] > j2) {
	      if (P[j] > j2) {
		swap(P[j],P[j2]);
		A.SwapCols(P[j],P[j2]);
		colnormsq.Swap(P[j],P[j2]);
#ifdef XDEBUG
		index.Swap(P[j],P[j2]);
#endif
	      } else {
		P[j] = P[j2];
	      }
	    } else {
	      if (P[j] > j2) P[j2] = P[j];
	      P[j] = j2;
	    }
#ifdef XDEBUG
	    CheckIndex(index,P,origj2);
#endif
	  }

	  // Find Householder matrix for this column
	  // This multiplies through to the end of the original block.
	  // This way, when we are done, the whole block has had the
	  // same Householder reflections applied to it.
	  beta(j) = Householder_Reflect(A.SubMatrix(j,M,j,origj2),det);

	  // Update Z:
	  BlockHouseholder_Augment(A.SubMatrix(j1,M,j1,j+1),
	      Z.SubTriMatrix(0,j-j1+1),CONJ(beta(j)));

	  // Update the colnormsq values within this block
	  // (No need to go all the way to origj2, since the j2..origj2
	  // columns are those with low norm already - we don't need
	  // those values until we recalculate them from scratch anyway.)
	  const T* Ajk = A.row(j,j+1,j2).cptr();
	  for(size_t k=j+1;k<j2;++k,Ajk+=Astepj) 
	    colnormsq(k) -= tmv::NORM(*Ajk);
	}
	//cerr<<"done: j2 = "<<j2<<endl;

	if (j1 < j2) {

	  // Do the Block Householder update of the rest of the matrix:
	  BlockHouseholder_LDiv(A.SubMatrix(j1,M,j1,j2),
	      Z.SubTriMatrix(0,j2-j1),A.SubMatrix(j1,M,origj2,N));

	  // Update the colnormsq values for the rest of the matrix:
	  if (M-j2 > j2-j1)
	    for(size_t k=origj2;k<N;++k) colnormsq(k) -= 
	      A.col(k,j1,j2).NormSq();
	  else 
	    for(size_t k=origj2;k<N;++k) colnormsq(k) = 
	      A.col(k,j2,M).NormSq();
	}

	if (j2 < origj2) {
	  //cerr<<"j2 = "<<j2<<", origj2 = "<<origj2<<endl;
	  //cerr<<"P["<<j2<<".."<<origj2<<"] = ";
	  //for(size_t j=j2; j<origj2; ++j) cerr<<P[j]<<" ";
	  
#ifdef XDEBUG
	  CheckIndex(index,P,origj2);
#endif
	  // Put the bad columns back where they started before this loop:
	  for(size_t j=j2; j<origj2; ++j) if (P[j] > j2) {
	    //cerr<<"Sending "<<j<<" back to "<<P[j]<<endl;
	    A.SwapCols(j,P[j]);
	    colnormsq.Swap(j,P[j]);
#ifdef XDEBUG
	    index.Swap(j,P[j]);
#endif
	  }
#ifdef XDEBUG
	  CheckIndex(index,P,j2);
#endif
	}

	j1 = j1x = j2;
	//cerr<<"For next loop: j1 = "<<j1<<", j2 = "<<j2<<endl;
      }
      //cerr<<"Done main loop pass\n";

      if (j3 < N) {
	// Then need to recalculate some of the colnorms:
	for(size_t k=j3;k<N;++k) colnormsq(k) = NormSq(A.col(k,j3,M));
	//cerr<<"New colnorms are "<<colnormsq.SubVector(j3,N)<<endl;
      }
    }

    if (det != T(0)) {
      for(size_t i=0;i<N;++i) if (P[i] != i) det = -det;
    }

#ifdef XDEBUG
    CheckIndex(index,P,N);
    Matrix<T> Q = A;
    GetQFromQR(Q.View(),beta);
    Matrix<T> AA = Q*UpperTriMatrixViewOf(A);
    AA.ReversePermuteCols(P);
    if (Norm(AA-A0) > 0.001*max(RealType(T)(1),Norm(A0))) {
      cerr<<"LooseBlockQRP_Decompose: \n";
      cerr<<"A = "<<Type(A)<<endl;
      if (N < 100) {
	cerr<<"  "<<A0<<endl;
	cerr<<"-> "<<A<<endl;
	cerr<<"beta = "<<beta<<endl;
	cerr<<"P = ";
	for(size_t i=0;i<N;i++) cerr<<P[i]<<" ";
	cerr<<endl;
	cerr<<"QRP = "<<AA<<endl;
	Matrix<T> diff = AA-A0;
	diff.Clip(0.0001);
	cerr<<"diff = "<<diff<<endl;
      }
      cerr<<"Rdiag = "<<A.diag()<<endl;
      cerr<<"Norm(A-QRP) = "<<Norm(AA-A0)<<endl;
      abort(); 
    }
#endif
  }

  template <class T> void NonLapQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, size_t* P, T& det)
  {
    // Decompose A (input as A) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (A), with the beta's for
    // the Householder matrices returned in beta.

    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(beta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    if (A.rowsize() > QRP_BLOCKSIZE)
      if (StrictQRP)
	StrictBlockQRP_Decompose(A,beta,P,det);
      else
	LooseBlockQRP_Decompose(A,beta,P,det);
    else
      NonBlockQRP_Decompose(A,beta,P,det);
  }

#ifdef XLAP
  template <class T> inline void LapQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, size_t* P, T& det)
  { NonLapQR_Decompose(A,beta,det); }
  template <> void LapQRP_Decompose(const MatrixView<double>& A,
      const VectorView<double>& beta, size_t* P, double& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lap_p[n];
    int lwork = 3*n;
    double* work = LAP_DWork(lwork);
    int info;
    int lda = A.stepj();
    char cc = 'F';
    double thresh = Epsilon<double>()*dlange(&cc,&m,&n,A.ptr(),&lda,0);
    dgeqpf(&m,&n,A.ptr(),&lda,lap_p,beta.ptr(),work,&info);
    LAP_Results(info,"dgeqpf");
    for(size_t i=0;i<A.rowsize();++i) {
      if (abs(A(i,i)) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	beta.SubVector(i,beta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<beta.size();++i) {
      if (det) if (beta(i) != 0.) det = -det;
      P[i] = lap_p[i]-1;
      if (det) if (P[i] != i) det = -det;
    }
  }
  template <> void LapQRP_Decompose(
      const MatrixView<complex<double> >& A,
      const VectorView<complex<double> >& beta, size_t* P,
      complex<double>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lap_p[n];
    int lwork = 2*n;
    double* rwork = LAP_DWork(lwork);
    lwork = 3*n;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    int lda = A.stepj();
    char cc = 'F';
    double thresh = Epsilon<double>()*
      zlange(&cc,&m,&n,LAP_Complex(A.ptr()),&lda,0);
    zgeqpf(&m,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(beta.ptr()),
	LAP_Complex(work),rwork,&info);
    LAP_Results(info,"zgeqpf");
    beta.ConjugateSelf();
    for(size_t i=0;i<A.rowsize();++i) {
      TMVAssert(abs(imag(A(i,i))) < Epsilon<double>());
      if (abs(real(A(i,i))) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	beta.SubVector(i,beta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<beta.size();++i) {
      if (det!=double(0)) {
	if (imag(beta(i)) != 0.) 
	  det *= -conj(beta(i)*beta(i))/norm(beta(i));
	else if (real(beta(i)) != 0.)
	  det = -det;
      }
      P[i] = lap_p[i]-1;
      if (det!=double(0)) if (P[i] != i) det = -det;
    }
  }
#ifndef NOFLOAT
  template <> void LapQRP_Decompose(const MatrixView<float>& A,
      const VectorView<float>& beta, size_t* P, float& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lap_p[n];
    int lwork = 3*n;
    float* work = LAP_SWork(lwork);
    int info;
    int lda = A.stepj();
    char cc = 'F';
    float thresh = Epsilon<float>()*slange(&cc,&m,&n,A.ptr(),&lda,0);
    sgeqpf(&m,&n,A.ptr(),&lda,lap_p,beta.ptr(),work,&info);
    LAP_Results(info,"sgeqpf");
    for(size_t i=0;i<A.rowsize();++i) {
      if (abs(A(i,i)) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	beta.SubVector(i,beta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<beta.size();++i) {
      if (det) if (beta(i) != 0.) det = -det;
      P[i] = lap_p[i]-1;
      if (det) if (P[i] != i) det = -det;
    }
  }
  template <> void LapQRP_Decompose(
      const MatrixView<complex<float> >& A,
      const VectorView<complex<float> >& beta, size_t* P,
      complex<float>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lap_p[n];
    int lwork = 2*n;
    float* rwork = LAP_SWork(lwork);
    lwork = 3*n;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    int lda = A.stepj();
    char cc = 'F';
    float thresh = Epsilon<float>()*
      clange(&cc,&m,&n,LAP_Complex(A.ptr()),&lda,0);
    cgeqpf(&m,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(beta.ptr()),
	LAP_Complex(work),rwork,&info);
    LAP_Results(info,"cgeqpf");
    for(size_t i=0;i<A.rowsize();++i) {
      TMVAssert(abs(imag(A(i,i))) < Epsilon<double>());
      if (abs(real(A(i,i))) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	beta.SubVector(i,beta.size()).Zero();
	break;
      }
    }
    beta.ConjugateSelf();
    for(size_t i=0;i<beta.size();++i) {
      if (det!=float(0)) {
	if (imag(beta(i)) != 0.) 
	  det *= -conj(beta(i)*beta(i))/norm(beta(i));
	else if (real(beta(i)) != 0.)
	  det = -det;
      }
      P[i] = lap_p[i]-1;
      if (det!=float(0)) if (P[i] != i) det = -det;
    }
  }
#endif
  template <class T> inline void NewLapQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, size_t* P, T& det)
  { NonLapQR_Decompose(A,beta,det); }
  template <> void NewLapQRP_Decompose(const MatrixView<double>& A,
      const VectorView<double>& beta, size_t* P, double& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lap_p[n];
    int lwork = 3*n+1;
    double* work = LAP_DWork(lwork);
    int info;
    int lda = A.stepj();
    char cc = 'F';
    double thresh = Epsilon<double>()*dlange(&cc,&m,&n,A.ptr(),&lda,0);
    dgeqp3(&m,&n,A.ptr(),&lda,lap_p,beta.ptr(),work,&lwork,&info);
    LAP_Results(info,"dgeqp3");
    for(size_t i=0;i<A.rowsize();++i) {
      if (abs(A(i,i)) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	beta.SubVector(i,beta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<beta.size();++i) {
      if (det) if (beta(i) != 0.) det = -det;
      P[i] = lap_p[i]-1;
      if (det) if (P[i] != i) det = -det;
    }
  }
  template <> void NewLapQRP_Decompose(
      const MatrixView<complex<double> >& A,
      const VectorView<complex<double> >& beta, size_t* P,
      complex<double>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lap_p[n];
    int lwork = 2*n;
    double* rwork = LAP_DWork(lwork);
    lwork = n+1;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    int lda = A.stepj();
    char cc = 'F';
    double thresh = Epsilon<double>()*
      zlange(&cc,&m,&n,LAP_Complex(A.ptr()),&lda,0);
    zgeqp3(&m,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(beta.ptr()),
	LAP_Complex(work),&lwork,rwork,&info);
    LAP_Results(info,"zgeqp3");
    beta.ConjugateSelf();
    for(size_t i=0;i<A.rowsize();++i) {
      TMVAssert(abs(imag(A(i,i))) < Epsilon<double>());
      if (abs(real(A(i,i))) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	beta.SubVector(i,beta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<beta.size();++i) {
      if (det!=double(0)) {
	if (imag(beta(i)) != 0.) 
	  det *= -conj(beta(i)*beta(i))/norm(beta(i));
	else if (real(beta(i)) != 0.)
	  det = -det;
      }
      P[i] = lap_p[i]-1;
      if (det!=double(0)) if (P[i] != i) det = -det;
    }
  }
#ifndef NOFLOAT
  template <> void NewLapQRP_Decompose(const MatrixView<float>& A,
      const VectorView<float>& beta, size_t* P, float& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lap_p[n];
    int lwork = 3*n+1;
    float* work = LAP_SWork(lwork);
    int info;
    int lda = A.stepj();
    char cc = 'F';
    float thresh = Epsilon<float>()*slange(&cc,&m,&n,A.ptr(),&lda,0);
    sgeqp3(&m,&n,A.ptr(),&lda,lap_p,beta.ptr(),work,&lwork,&info);
    LAP_Results(info,"sgeqp3");
    for(size_t i=0;i<A.rowsize();++i) {
      if (abs(A(i,i)) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	beta.SubVector(i,beta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<beta.size();++i) {
      if (det) if (beta(i) != 0.) det = -det;
      P[i] = lap_p[i]-1;
      if (det) if (P[i] != i) det = -det;
    }
  }
  template <> void NewLapQRP_Decompose(
      const MatrixView<complex<float> >& A,
      const VectorView<complex<float> >& beta, size_t* P,
      complex<float>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lap_p[n];
    int lwork = 2*n;
    float* rwork = LAP_SWork(lwork);
    lwork = n+1;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    int lda = A.stepj();
    char cc = 'F';
    float thresh = Epsilon<float>()*
      clange(&cc,&m,&n,LAP_Complex(A.ptr()),&lda,0);
    cgeqp3(&m,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(beta.ptr()),
	LAP_Complex(work),&lwork,rwork,&info);
    LAP_Results(info,"cgeqp3");
    for(size_t i=0;i<A.rowsize();++i) {
      TMVAssert(abs(imag(A(i,i))) < Epsilon<double>());
      if (abs(real(A(i,i))) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	beta.SubVector(i,beta.size()).Zero();
	break;
      }
    }
    beta.ConjugateSelf();
    for(size_t i=0;i<beta.size();++i) {
      if (det!=float(0)) {
	if (imag(beta(i)) != 0.) 
	  det *= -conj(beta(i)*beta(i))/norm(beta(i));
	else if (real(beta(i)) != 0.)
	  det = -det;
      }
      P[i] = lap_p[i]-1;
      if (det!=float(0)) if (P[i] != i) det = -det;
    }
  }
#endif
#endif
  template <class T> void QRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, size_t* P, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    if (A.rowsize() > 0) {
      // Neither LAP version of QRP seems to do any actual pivoting.
      // At least for the Intel MKL version.
      // So until I figure out why and how to get them to work, leave
      // this is XLAP (ie. ignore the LAP routines.
#ifdef XLAP
      //cerr<<"A.stor = "<<Text(A.stor())<<endl;
      if (A.iscm())
	//LapQRP_Decompose(A,beta,P,det);
	NewLapQRP_Decompose(A,beta,P,det);
      else {
#ifdef TMVDEBUG
	cout<<"Lap QRDecomp: A, beta are wrong step:\n";
	cout<<"A isrm = "<<A.isrm()<<", stepi = "<<A.stepi()<<endl;
	cout<<"beta step = "<<beta.step()<<endl;
#endif
	NonLapQRP_Decompose(A,beta,P,det);
      }
#else
      NonLapQRP_Decompose(A,beta,P,det);
#endif
    }
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
    GetQFromQR(Q,beta.View());
  }

  template <class T> void QRP_Decompose(
      const MatrixView<T>& Q, const UpperTriMatrixView<T>& R, size_t* P, T& det)
  {
    // Decompose A (input as Q) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation

    const size_t N = Q.rowsize();
    TMVAssert(Q.colsize() >= N);
    TMVAssert(R.colsize() == N);
    TMVAssert(R.rowsize() == N);
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(R.ct() == NonConj);

    Vector<T> beta(N);
    QRP_Decompose(Q,beta.View(),P,det);
    R = UpperTriMatrixViewOf(Q);
    GetQFromQR(Q,beta.View());
  }

  //
  // Packed Q - LDivEq
  //

  template <class T, class T1> void NonBlockQ_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == m.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    // Solve Q x = m in place 
    // where Q is stored as Householder vectors along with beta
    //
    // Q is H0t H1t ... H_N-1t
    // So x = H_N-1 .. H1 H0 m
    //
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    for(size_t j=0;j<N;++j) if (beta(j) != T1(0)) {
      Householder_LMult(Q.col(j,j+1,M),beta(j),m.Rows(j,M));
    }
  }

  template <class T, class T1> void BlockQ_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == m.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T> m0 = m;
    Matrix<T> m2 = m;
    NonBlockQ_LDivEq(Q,beta,m2.View());
#endif

    // x = H_N-1 .. H1 H0 m
    // In first block step:
    // m = (Hr .. H0) m
    //   = (H0t .. Hrt)t m
    // So form Y,Z from Ht's, rather than H's, and then call LDiv
    // Ht just means use CONJ(beta) rather than beta.
    
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    UpperTriMatrix<T1,NonUnitDiag,ColMajor> BaseZ(min(QR_BLOCKSIZE,N));
    for(size_t j1=0;j1<N;) {
      size_t j2 = min(N,j1+QR_BLOCKSIZE);
      ConstMatrixView<T1> Y = Q.SubMatrix(j1,M,j1,j2);
      UpperTriMatrixView<T1> Z = BaseZ.SubTriMatrix(0,Y.rowsize());
      BlockHouseholder_MakeZ(Y,Z,beta.SubVector(j1,j2));
      BlockHouseholder_LDiv(Y,Z,m.Rows(j1,M));
      j1 = j2;
    }
#ifdef XDEBUG
    if (Norm(m-m2) > 0.001*max(RealType(T)(1),Norm(m2))) {
      cerr<<"BlockQ_LDivEq: Q = "<<Type(Q)<<"  "<<Q<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"NonBlock m = "<<m2<<endl;
      abort(); 
    }
#endif
  }

  template <class T, class T1> inline void NonLapQ_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == m.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    if (Q.rowsize() > QR_BLOCKSIZE && m.rowsize() > QR_BLOCKSIZE)
      BlockQ_LDivEq(Q,beta,m);
    else
      NonBlockQ_LDivEq(Q,beta,m);
  }
#ifdef LAP
  template <class T, class T1> inline void LapQ_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  { NonLapQ_LDivEq(Q,beta,m); }
  template <> void LapQ_LDivEq(const GenMatrix<double>& Q,
      const GenVector<double>& beta, const MatrixView<double>& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == x.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    char side = x.isrm() ? 'R' : 'L';
    int ldx = x.isrm() ? x.stepi() : x.stepj();
    int m = x.isrm() ? x.rowsize() : x.colsize();
    int n = x.isrm() ? x.colsize() : x.rowsize();
    int k = Q.rowsize();
    int lwork = x.rowsize()*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    if (Q.isrm()) {
      char trans = x.isrm() ? 'T' : 'N';
      int ldq = Q.stepi();
      dormlq(&side,&trans,&m,&n,&k,const_cast<double*>(Q.cptr()),&ldq,
	  const_cast<double*>(beta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dormlq");
    } else {
      char trans = x.isrm() ? 'N' : 'T';
      int ldq = Q.stepj();
      dormqr(&side,&trans,&m,&n,&k,const_cast<double*>(Q.cptr()),&ldq,
	  const_cast<double*>(beta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dormqr");
    }
  }
  template <> void LapQ_LDivEq(
      const GenMatrix<complex<double> >& Q,
      const GenVector<complex<double> >& beta,
      const MatrixView<complex<double> >& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == x.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int k = Q.rowsize();
    int lwork = x.rowsize()*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    if (Q.isrm()) {
      if (x.isrm()) {
	char side = 'R';
	char trans = 'C';
	int ldx = x.stepi();
	int m = x.rowsize();
	int n = x.colsize();
	int ldq = Q.stepi();
	if (x.isconj()) x.ConjugateSelf();
	zunmlq(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	    LAP_Complex(beta.cptr()),LAP_Complex(x.ptr()),&ldx,
	    LAP_Complex(work),&lwork,&info);
	if (x.isconj()) x.ConjugateSelf();
        LAP_Results(info,int(real(work[0])),m,n,lwork,"zunmlq");
      } else {
	char side = 'L';
	char trans = 'N';
	int ldx = x.stepj();
	int m = x.colsize();
	int n = x.rowsize();
	int ldq = Q.stepi();
	if (!x.isconj()) x.ConjugateSelf();
	zunmlq(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	    LAP_Complex(beta.cptr()),LAP_Complex(x.ptr()),&ldx,
	    LAP_Complex(work),&lwork,&info);
	if (!x.isconj()) x.ConjugateSelf();
        LAP_Results(info,int(real(work[0])),m,n,lwork,"zunmlq");
      }
    } else {
      if (x.isrm()) {
	char side = 'R';
	char trans = 'N';
	int ldx = x.stepi();
	int m = x.rowsize();
	int n = x.colsize();
	int ldq = Q.stepj();
	if (!x.isconj()) x.ConjugateSelf();
	Vector<complex<double> > conjbeta = beta.Conjugate();
	zunmqr(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	    LAP_Complex(conjbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
	    LAP_Complex(work),&lwork,&info);
	if (!x.isconj()) x.ConjugateSelf();
        LAP_Results(info,int(real(work[0])),m,n,lwork,"zunmqr");
      } else {
	char side = 'L';
	char trans = 'C';
	int ldx = x.stepj();
	int m = x.colsize();
	int n = x.rowsize();
	int ldq = Q.stepj();
	if (x.isconj()) x.ConjugateSelf();
	Vector<complex<double> > conjbeta = beta.Conjugate();
	zunmqr(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	    LAP_Complex(conjbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
	    LAP_Complex(work),&lwork,&info);
	if (x.isconj()) x.ConjugateSelf();
        LAP_Results(info,int(real(work[0])),m,n,lwork,"zunmqr");
      }
    }
  }
#ifndef NOFLOAT
  template <> void LapQ_LDivEq(const GenMatrix<float>& Q,
      const GenVector<float>& beta, const MatrixView<float>& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == x.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    char side = x.isrm() ? 'R' : 'L';
    int ldx = x.isrm() ? x.stepi() : x.stepj();
    int m = x.isrm() ? x.rowsize() : x.colsize();
    int n = x.isrm() ? x.colsize() : x.rowsize();
    int k = Q.rowsize();
    int lwork = x.rowsize()*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    if (Q.isrm()) {
      char trans = x.isrm() ? 'T' : 'N';
      int ldq = Q.stepi();
      sormlq(&side,&trans,&m,&n,&k,const_cast<float*>(Q.cptr()),&ldq,
	  const_cast<float*>(beta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sormlq");
    } else {
      char trans = x.isrm() ? 'N' : 'T';
      int ldq = Q.stepj();
      sormqr(&side,&trans,&m,&n,&k,const_cast<float*>(Q.cptr()),&ldq,
	  const_cast<float*>(beta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sormqr");
    }
  }
  template <> void LapQ_LDivEq(
      const GenMatrix<complex<float> >& Q,
      const GenVector<complex<float> >& beta,
      const MatrixView<complex<float> >& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == x.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int k = Q.rowsize();
    int lwork = x.rowsize()*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    if (Q.isrm()) {
      if (x.isrm()) {
	char side = 'R';
	char trans = 'C';
	int ldx = x.stepi();
	int m = x.rowsize();
	int n = x.colsize();
	int ldq = Q.stepi();
	if (x.isconj()) x.ConjugateSelf();
	cunmlq(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	    LAP_Complex(beta.cptr()),LAP_Complex(x.ptr()),&ldx,
	    LAP_Complex(work),&lwork,&info);
	if (x.isconj()) x.ConjugateSelf();
	LAP_Results(info,int(real(work[0])),m,n,lwork,"cunmlq");
      } else {
	char side = 'L';
	char trans = 'N';
	int ldx = x.stepj();
	int m = x.colsize();
	int n = x.rowsize();
	int ldq = Q.stepi();
	if (!x.isconj()) x.ConjugateSelf();
	cunmlq(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	    LAP_Complex(beta.cptr()),LAP_Complex(x.ptr()),&ldx,
	    LAP_Complex(work),&lwork,&info);
	if (!x.isconj()) x.ConjugateSelf();
	LAP_Results(info,int(real(work[0])),m,n,lwork,"cunmlq");
      }
    } else {
      if (x.isrm()) {
	char side = 'R';
	char trans = 'N';
	int ldx = x.stepi();
	int m = x.rowsize();
	int n = x.colsize();
	int ldq = Q.stepj();
	if (!x.isconj()) x.ConjugateSelf();
	Vector<complex<float> > conjbeta = beta.Conjugate();
	cunmqr(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	    LAP_Complex(conjbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
	    LAP_Complex(work),&lwork,&info);
	if (!x.isconj()) x.ConjugateSelf();
	LAP_Results(info,int(real(work[0])),m,n,lwork,"cunmqr");
      } else {
	char side = 'L';
	char trans = 'C';
	int ldx = x.stepj();
	int m = x.colsize();
	int n = x.rowsize();
	int ldq = Q.stepj();
	if (x.isconj()) x.ConjugateSelf();
	Vector<complex<float> > conjbeta = beta.Conjugate();
	cunmqr(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	    LAP_Complex(conjbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
	    LAP_Complex(work),&lwork,&info);
	if (x.isconj()) x.ConjugateSelf();
	LAP_Results(info,int(real(work[0])),m,n,lwork,"cunmqr");
      }
    }
  }
#endif
#endif

  template <class T, class T1> void Q_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(m.colsize() == Q.colsize());
    TMVAssert(Q.isrm() || Q.iscm());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    if (m.colsize() > 0 && m.rowsize() > 0) {
#ifdef LAP
      if ( m.isrm() || m.iscm() )
	LapQ_LDivEq(Q,beta,m);
      else
#endif
	NonLapQ_LDivEq(Q,beta,m);
    }
  }

  //
  // Packed Q - RDivEq
  //

  template <class T, class T1> void NonBlockQ_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(m.rowsize() == Q.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    // Solve x Q = m in place
    // where Q is stored as Householder vectors along with beta
    //
    // x = m Qt 
    // Qt is H_N-1 H_N-2 ... H1 H0
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    for(int j=N-1;j>=0;--j) if (beta(j) != T1(0)) {
      Householder_LMult(Q.col(j,j+1,M).Conjugate(),beta(j),
	  m.Cols(j,M).Transpose());
    }
  }

  template <class T, class T1> void BlockQ_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == m.rowsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T> m0 = m;
    Matrix<T> m2 = m;
    NonBlockQ_RDivEq(Q,beta,m2.View());
#endif

    // x = m Qt 
    // x = m H_N-1 H_N-2 ... H1 H0
    // Again form Y,Z from Ht's, rather than H's, and then call RDiv
    
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    UpperTriMatrix<T1,NonUnitDiag,ColMajor> BaseZ(min(QR_BLOCKSIZE,N));
    for(size_t j2=N;j2>0;) {
      size_t j1 = j2 > QR_BLOCKSIZE ? j2-QR_BLOCKSIZE : 0;
      ConstMatrixView<T1> Y = Q.SubMatrix(j1,M,j1,j2);
      UpperTriMatrixView<T1> Z = BaseZ.SubTriMatrix(0,Y.rowsize());
      BlockHouseholder_MakeZ(Y,Z,beta.SubVector(j1,j2));
      BlockHouseholder_LMult(Y,Z,m.Cols(j1,M).Adjoint());
      j2 = j1;
    }
#ifdef XDEBUG
    if (Norm(m-m2) > 0.001*max(RealType(T)(1),Norm(m2))) {
      cerr<<"BlockQ_RDivEq: Q = "<<Type(Q)<<"  "<<Q<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"NonBlock m = "<<m2<<endl;
      abort(); 
    }
#endif
  }

  template <class T, class T1> inline void NonLapQ_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == m.rowsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    if (Q.rowsize() > QR_BLOCKSIZE && m.colsize() > QR_BLOCKSIZE)
      BlockQ_RDivEq(Q,beta,m);
    else
      NonBlockQ_RDivEq(Q,beta,m);
  }

#ifdef LAP
  template <class T, class T1> inline void LapQ_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  { NonLapQ_RDivEq(Q,beta,m); }
  template <> void LapQ_RDivEq(const GenMatrix<double>& Q,
      const GenVector<double>& beta, const MatrixView<double>& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(x.rowsize() == Q.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    char side = x.isrm() ? 'L' : 'R';
    int ldx = x.isrm() ? x.stepi() : x.stepj();
    int m = x.isrm() ? x.rowsize() : x.colsize();
    int n = x.isrm() ? x.colsize() : x.rowsize();
    int k = Q.rowsize();
    int lwork = x.colsize()*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    if (Q.isrm()) {
      char trans = x.isrm() ? 'T' : 'N';
      int ldq = Q.stepi();
      dormlq(&side,&trans,&m,&n,&k,const_cast<double*>(Q.cptr()),&ldq,
	  const_cast<double*>(beta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dormlq");
    } else {
      char trans = x.isrm() ? 'N' : 'T';
      int ldq = Q.stepj();
      dormqr(&side,&trans,&m,&n,&k,const_cast<double*>(Q.cptr()),&ldq,
	  const_cast<double*>(beta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dormqr");
    }
  }
  template <> void LapQ_RDivEq(const GenMatrix<complex<double> >& Q,
      const GenVector<complex<double> >& beta,
      const MatrixView<complex<double> >& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(x.rowsize() == Q.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int k = Q.rowsize();
    int lwork = x.colsize()*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    if (Q.isrm()) {
      if (x.isrm()) {
	char side = 'L';
	char trans = 'C';
	int ldx = x.stepi();
	int m = x.rowsize();
	int n = x.colsize();
	int ldq = Q.stepi();
	if (x.isconj()) x.ConjugateSelf();
	zunmlq(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	    LAP_Complex(beta.cptr()),LAP_Complex(x.ptr()),&ldx,
	    LAP_Complex(work),&lwork,&info);
	if (x.isconj()) x.ConjugateSelf();
	LAP_Results(info,int(real(work[0])),m,n,lwork,"zunmlq");
      } else {
	char side = 'R';
	char trans = 'N';
	int ldx = x.stepj();
	int m = x.colsize();
	int n = x.rowsize();
	int ldq = Q.stepi();
	if (!x.isconj()) x.ConjugateSelf();
	zunmlq(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	    LAP_Complex(beta.cptr()),LAP_Complex(x.ptr()),&ldx,
	    LAP_Complex(work),&lwork,&info);
	if (!x.isconj()) x.ConjugateSelf();
	LAP_Results(info,int(real(work[0])),m,n,lwork,"zunmlq");
      }
    } else {
      if (x.isrm()) {
	char side = 'L';
	char trans = 'N';
	int ldx = x.stepi();
	int m = x.rowsize();
	int n = x.colsize();
	int ldq = Q.stepj();
	if (!x.isconj()) x.ConjugateSelf();
	Vector<complex<double> > conjbeta = beta.Conjugate();
	zunmqr(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	    LAP_Complex(conjbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
	    LAP_Complex(work),&lwork,&info);
	if (!x.isconj()) x.ConjugateSelf();
	LAP_Results(info,int(real(work[0])),m,n,lwork,"zunmqr");
      } else {
	char side = 'R';
	char trans = 'C';
	int ldx = x.stepj();
	int m = x.colsize();
	int n = x.rowsize();
	int ldq = Q.stepj();
	if (x.isconj()) x.ConjugateSelf();
	Vector<complex<double> > conjbeta = beta.Conjugate();
	zunmqr(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	    LAP_Complex(conjbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
	    LAP_Complex(work),&lwork,&info);
	if (x.isconj()) x.ConjugateSelf();
	LAP_Results(info,int(real(work[0])),m,n,lwork,"zunmqr");
      }
    }
  }
#ifndef NOFLOAT
  template <> void LapQ_RDivEq(const GenMatrix<float>& Q,
      const GenVector<float>& beta, const MatrixView<float>& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(x.rowsize() == Q.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    char side = x.isrm() ? 'L' : 'R';
    int ldx = x.isrm() ? x.stepi() : x.stepj();
    int m = x.isrm() ? x.rowsize() : x.colsize();
    int n = x.isrm() ? x.colsize() : x.rowsize();
    int k = Q.rowsize();
    int lwork = x.colsize()*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    if (Q.isrm()) {
      char trans = x.isrm() ? 'T' : 'N';
      int ldq = Q.stepi();
      sormlq(&side,&trans,&m,&n,&k,const_cast<float*>(Q.cptr()),&ldq,
	  const_cast<float*>(beta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sormlq");
    } else {
      char trans = x.isrm() ? 'N' : 'T';
      int ldq = Q.stepj();
      sormqr(&side,&trans,&m,&n,&k,const_cast<float*>(Q.cptr()),&ldq,
	  const_cast<float*>(beta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sormqr");
    }
  }
  template <> void LapQ_RDivEq(const GenMatrix<complex<float> >& Q,
      const GenVector<complex<float> >& beta,
      const MatrixView<complex<float> >& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(x.rowsize() == Q.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int k = Q.rowsize();
    int lwork = x.colsize()*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    if (Q.isrm()) {
      if (x.isrm()) {
	char side = 'L';
	char trans = 'C';
	int ldx = x.stepi();
	int m = x.rowsize();
	int n = x.colsize();
	int ldq = Q.stepi();
	if (x.isconj()) x.ConjugateSelf();
	  cunmlq(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	      LAP_Complex(beta.cptr()),LAP_Complex(x.ptr()),&ldx,
	      LAP_Complex(work),&lwork,&info);
	if (x.isconj()) x.ConjugateSelf();
        LAP_Results(info,int(real(work[0])),m,n,lwork,"cunmlq");
      } else {
	char side = 'R';
	char trans = 'N';
	int ldx = x.stepj();
	int m = x.colsize();
	int n = x.rowsize();
	int ldq = Q.stepi();
	if (!x.isconj()) x.ConjugateSelf();
	  cunmlq(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	      LAP_Complex(beta.cptr()),LAP_Complex(x.ptr()),&ldx,
	      LAP_Complex(work),&lwork,&info);
	if (!x.isconj()) x.ConjugateSelf();
        LAP_Results(info,int(real(work[0])),m,n,lwork,"cunmlq");
      }
    } else {
      if (x.isrm()) {
	char side = 'L';
	char trans = 'N';
	int ldx = x.stepi();
	int m = x.rowsize();
	int n = x.colsize();
	int ldq = Q.stepj();
	if (!x.isconj()) x.ConjugateSelf();
	  Vector<complex<float> > conjbeta = beta.Conjugate();
	  cunmqr(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	      LAP_Complex(conjbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
	      LAP_Complex(work),&lwork,&info);
	if (!x.isconj()) x.ConjugateSelf();
        LAP_Results(info,int(real(work[0])),m,n,lwork,"cunmqr");
      } else {
	char side = 'R';
	char trans = 'C';
	int ldx = x.stepj();
	int m = x.colsize();
	int n = x.rowsize();
	int ldq = Q.stepj();
	if (x.isconj()) x.ConjugateSelf();
	  Vector<complex<float> > conjbeta = beta.Conjugate();
	  cunmqr(&side,&trans,&m,&n,&k,LAP_Complex(Q.cptr()),&ldq,
	      LAP_Complex(conjbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
	      LAP_Complex(work),&lwork,&info);
	if (x.isconj()) x.ConjugateSelf();
        LAP_Results(info,int(real(work[0])),m,n,lwork,"cunmqr");
      }
    }
  }
#endif
#endif

  template <class T, class T1> void Q_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(m.rowsize() == Q.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    if (m.colsize() > 0 && m.rowsize() > 0) {
#ifdef LAP
      if (m.isrm() || m.iscm())
	LapQ_RDivEq(Q,beta,m);
      else
#endif
	NonLapQ_RDivEq(Q,beta,m);
    }
  }

  //
  // LDiv
  //

  template <class T, class T1, class T2> void QR_LDiv(
      const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const size_t* P,
      const GenMatrix<T2>& m, const MatrixView<T>& x, size_t N1)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(beta.size() == QRx.rowsize());
    TMVAssert(m.colsize() == QRx.colsize());
    TMVAssert(x.colsize() == QRx.rowsize());
    TMVAssert(x.rowsize() == m.rowsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    // Solve Q R P x = m
    // where Q and R are stored in QRx, and beta are the beta

    // First Solve Q y = m
    if (QRx.IsSquare()) {
      x = m;
      Q_LDivEq(QRx,beta,x);
    } else {
      if (m.isrm()) {
	Matrix<T,RowMajor> m1 = m;
	// Q is Q1 [ I ]
	//         [ 0 ]
	// where Q1 is the part of Q that is stored in QRx and beta
	// m1 = Q^-1 m1
	Q_LDivEq(QRx,beta,m1.View());
	// y = [ I 0 ] m1
	x = m1.Rows(0,x.colsize()); // x = y here
      } else {
	Matrix<T,ColMajor> m1 = m;
	Q_LDivEq(QRx,beta,m1.View());
	x = m1.Rows(0,x.colsize()); // x = y here
      }
    }

    // Now solve R z = y
    x.Rows(N1,x.colsize()).Zero();
    x.Rows(0,N1) /= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);

    // Finally P x = z
    if (P) x.ReversePermuteRows(P);
  }

  //
  // LDivEq
  //

  template <class T, class T1> void QR_LDivEq(
      const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const size_t* P,
      const MatrixView<T>& m, size_t N1)
  {
    TMVAssert(QRx.colsize() == QRx.rowsize());
    TMVAssert(beta.size() == QRx.rowsize());
    TMVAssert(m.colsize() == QRx.colsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    // Solves Q R P x = m in place (m <- x)
    Q_LDivEq(QRx,beta,m);
    m.Rows(N1,m.colsize()).Zero();
    m.Rows(0,N1) /= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);
    if (P) m.ReversePermuteRows(P);
  }


  //
  // RDiv
  //

  template <class T, class T1, class T2> void QR_RDiv(
      const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const size_t* P,
      const GenMatrix<T2>& m, const MatrixView<T>& x, size_t N1)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(beta.size() == QRx.rowsize());
    TMVAssert(x.rowsize() == QRx.colsize());
    TMVAssert(m.rowsize() == QRx.rowsize());
    TMVAssert(x.colsize() == m.colsize());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    // Solve x Q R P = m
    // where Q and R are stored in QRx, and beta are the beta

    // First solve y P = m
    x.Cols(0,m.rowsize()) = m;
    if (P) x.Cols(0,m.rowsize()).PermuteCols(P);

    // Next solve z R = y by forward substitution
    x.Cols(N1,x.rowsize()).Zero();
    x.Cols(0,N1) %= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);

    // Finally solve x Q = z
    // Q = Q1 [ I ]
    //        [ 0 ]
    // where Q1 is the part of Q that is stored in QRx and beta
    // We've already dealt with the first part by zeroing out the 
    // right columns of x.
    Q_RDivEq(QRx,beta,x);
  }

  //
  // RDivEq
  //

  template <class T, class T1> void QR_RDivEq(
      const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const size_t* P,
      const MatrixView<T>& m, size_t N1)
  {
    TMVAssert(QRx.colsize() == QRx.rowsize());
    TMVAssert(beta.size() == QRx.rowsize());
    TMVAssert(m.rowsize() == QRx.colsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    // Solve x Q R P = m in place (m <- x)
 
    if (P) m.PermuteCols(P);
    m.Cols(N1,m.rowsize()).Zero();
    m.Cols(0,N1) %= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);
    Q_RDivEq(QRx,beta,m);
  }

  template <class T> T QRDiv<T>::Det() const
  {
    if (!donedet) {
      det *= DiagMatrixViewOf(QRx.diag()).Det();
      donedet = true;
    }
    return det;
  }

  template <class T> T QRPDiv<T>::Det() const
  {
    if (!donedet) {
      det *= DiagMatrixViewOf(QRx.diag()).Det();
      donedet = true;
    }
    return det;
  }

  template <class T> void QRDiv<T>::Inverse(const MatrixView<T>& minv) const
  {
    if (istrans) {
      TMVAssert(minv.colsize() == QRx.colsize());
      TMVAssert(minv.rowsize() == QRx.rowsize());
      LDiv(Eye<RealType(T),ColMajor>(QRx.rowsize()),minv);
    } else {
      TMVAssert(minv.colsize() == QRx.rowsize());
      TMVAssert(minv.rowsize() == QRx.colsize());
      RDiv(Eye<RealType(T),ColMajor>(QRx.rowsize()),minv);
    }
  }

  template <class T> void QRPDiv<T>::Inverse(const MatrixView<T>& minv) const
  {
    if (istrans) {
      TMVAssert(minv.colsize() == QRx.colsize());
      TMVAssert(minv.rowsize() == QRx.rowsize());
      LDiv(Eye<RealType(T),ColMajor>(QRx.rowsize()),minv);
    } else {
      TMVAssert(minv.colsize() == QRx.rowsize());
      TMVAssert(minv.rowsize() == QRx.colsize());
      RDiv(Eye<RealType(T),ColMajor>(QRx.rowsize()),minv);
    }
  }

  template <class T> void QRDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    // At A = Rt Qt Q R = Rt R
    // (At A)^-1 = (Rt R)^-1 = R^-1 * Rt^-1
    const size_t N = QRx.rowsize();

    UpperTriMatrix<T,NonUnitDiag,ColMajor> rinv(N);
    rinv.SetToIdentity();
    rinv /= UpperTriMatrixViewOf(QRx);
    minv = rinv * rinv.Adjoint();
  }

  template <class T> void QRPDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    // At A = Rt R
    // (At A)^-1 = (Rt R)^-1 = R^-1 * Rt^-1
    // AT A = PT RT R P
    // (At A)^-1 = (Pt Rt R P)^-1 = P^-1 R^-1 * (P^-1 R^-1)t
    const size_t N = QRx.rowsize();

    Matrix<T,ColMajor> temp(N,N);
    temp.SetToIdentity();
    UpperTriMatrixViewOf(temp).SubTriMatrix(0,N1) /= 
      UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);
    temp.ReversePermuteRows(P);
    minv = temp * temp.Adjoint();
  }

  template <class T> Matrix<T> QRDiv<T>::GetQ() const
  {
    Matrix<T> Q = QRx;
    GetQFromQR(Q.View(),beta);
    return Q;
  }

  template <class T> Matrix<T> QRPDiv<T>::GetQ() const
  {
    Matrix<T> Q = QRx;
    GetQFromQR(Q.View(),beta);
    return Q;
  }

#define InstFile "TMV_QRDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


