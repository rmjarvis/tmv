
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
    inplace(_inplace), Aptr(APTR), QRx(QRX), Qbeta(QRx.rowsize()),
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
    QR_Decompose(QRx,Qbeta.View(),det);
  }

  template <class T> QRPDiv<T>::QRPDiv(const GenMatrix<T>& A, bool _inplace) :
    istrans(A.colsize()<A.rowsize()),
    inplace(_inplace), Aptr(APTR), QRx(QRX), Qbeta(QRx.rowsize()),
    P(new size_t[Qbeta.size()]),
    det(T(1)), donedet(false), N1(Qbeta.size())
  {
    if (istrans) {
      if (inplace) { TMVAssert(A.Transpose() == QRx); }
      else QRx = A.Transpose();
    }
    else {
      if (inplace) { TMVAssert(A == QRx); }
      else QRx = A;
    }
    QRP_Decompose(QRx,Qbeta.View(),P,det);
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
      *fout << "M = "<< (istrans ? mm->Transpose() : mm->View()) <<endl;
      *fout << "Q = "<<GetQ()<<endl;
      *fout << "R = "<<GetR()<<endl;
    }
#endif
    Matrix<T> m2 = GetQ()*GetR();
    RealType(T) nm = Norm(m2- (istrans ? mm->Transpose() : mm->View()) );
    nm /= Norm(GetQ())*Norm(GetR());
#ifdef XDEBUG
    if (printmat) {
      *fout << "QR = "<<m2<<endl;
    }
#endif
    if (fout) {
      *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<"  "<<QRx.rowsize()*Epsilon<T>()<<endl;
    }
    return nm < QRx.rowsize()*Epsilon<T>();
  }

  template <class T> bool QRPDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenMatrix<T>* mm = dynamic_cast<const GenMatrix<T>*>(&m);
    TMVAssert(mm);
#ifdef XDEBUG
    bool printmat = fout && m.colsize() < 100 && m.rowsize() < 100;
    if (printmat) {
      *fout << "M = "<< (istrans ? mm->Transpose() : mm->View()) <<endl;
      *fout << "Q = "<<GetQ()<<endl;
      *fout << "R = "<<GetR()<<endl;
      *fout << "P = ";
      for(size_t i=0;i<QRx.rowsize();i++) *fout<<P[i]<<" ";
      *fout<<endl;
    }
#endif
    Matrix<T> m2 = GetQ()*GetR();
    m2.ReversePermuteCols(GetP());
    RealType(T) nm = Norm(m2- (istrans ? mm->Transpose() : mm->View()) );
    nm /= Norm(GetQ())*Norm(GetR());
#ifdef XDEBUG
    if (printmat) {
      *fout << "QRP = "<<m2<<endl;
    }
#endif
    if (fout) {
      *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<"  "<<QRx.rowsize()*Epsilon<T>()<<endl;
    }
    return nm < QRx.rowsize()*Epsilon<T>();
  }

  // 
  // Get Q from QR
  // 
 
  template <class T> void NonBlockGetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& Qbeta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    UpperTriMatrixViewOf(Q).Zero();
    for(int i=N-1;i>=0;i--) {
      Householder_Unpack(Q.SubMatrix(i,M,i,N),Qbeta(i));
    }
  }

  template <class T> void BlockGetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& Qbeta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
#ifdef XDEBUG
    Matrix<T> Q0 = Q;
    Matrix<T> Q2 = Q;
    NonBlockGetQFromQR(Q2.View(),Qbeta);
#endif
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    UpperTriMatrixViewOf(Q).Zero();
    UpperTriMatrix<T,NonUnitDiag,ColMajor> BaseZ(min(QR_BLOCKSIZE,N));
    for(size_t j2=N;j2>0;) {
      size_t j1 = j2 > QR_BLOCKSIZE ? j2-QR_BLOCKSIZE : 0;
      MatrixView<T> Y = Q.SubMatrix(j1,M,j1,j2);
      UpperTriMatrixView<T> Z = BaseZ.SubTriMatrix(0,Y.rowsize());
      BlockHouseholder_MakeZ(Y,Z,Qbeta.SubVector(j1,j2));
      BlockHouseholder_Unpack(Y,Z,Q.SubMatrix(j1,M,j2,N));
      j2 = j1;
    }
#ifdef XDEBUG
    if (Norm(Q-Q2) > 0.001*Norm(Q2)) {
      cerr<<"BlockGetQ: Q = "<<Type(Q)<<"  "<<Q0<<endl;
      cerr<<"Qbeta = "<<Qbeta<<endl;
      cerr<<"-> Q = "<<Q<<endl;
      cerr<<"NonBlock Q = "<<Q2<<endl;
      abort(); 
    }
#endif
  }

  template <class T> inline void NonLapGetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& Qbeta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    if (Q.rowsize() >= QR_BLOCKSIZE)
      BlockGetQFromQR(Q,Qbeta);
    else
      NonBlockGetQFromQR(Q,Qbeta);
  }
#ifdef LAP
  template <class T> inline void LapGetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& Qbeta)
  { NonLapGetQFromQR(Q,Qbeta); }
  template <> void LapGetQFromQR(const MatrixView<double>& Q,
      const GenVector<double>& Qbeta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
    int lwork = n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    if (Q.isrm()) {
      int lda = Q.stepi();
      dorglq(&n,&m,&n,Q.ptr(),&lda,const_cast<double*>(Qbeta.cptr()),
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dorglq");
    } else {
      int lda = Q.stepj();
      dorgqr(&m,&n,&n,Q.ptr(),&lda,const_cast<double*>(Qbeta.cptr()),
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dorgqr");
    }
  }
  template <> void LapGetQFromQR(
      const MatrixView<complex<double> >& Q,
      const GenVector<complex<double> >& Qbeta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
    int lwork = n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    if (Q.isrm()) {
      int lda = Q.stepi();
      zunglq(&n,&m,&n,LAP_Complex(Q.ptr()),&lda,LAP_Complex(Qbeta.cptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"zunglq");
    } else {
      int lda = Q.stepj();
      Vector<complex<double> > conjbeta = Qbeta.Conjugate();
      zungqr(&m,&n,&n,LAP_Complex(Q.ptr()),&lda,LAP_Complex(conjbeta.cptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"zungqr");
    }
  }
#ifndef NOFLOAT
  template <> void LapGetQFromQR(const MatrixView<float>& Q,
      const GenVector<float>& Qbeta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
    int lwork = n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    if (Q.isrm()) {
      int lda = Q.stepi();
      sorglq(&n,&m,&n,Q.ptr(),&lda,const_cast<float*>(Qbeta.cptr()),
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sorglq");
    } else {
      int lda = Q.stepj();
      sorgqr(&m,&n,&n,Q.ptr(),&lda,const_cast<float*>(Qbeta.cptr()),
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sorgqr");
    }
  }
  template <> void LapGetQFromQR(
      const MatrixView<complex<float> >& Q,
      const GenVector<complex<float> >& Qbeta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
    int lwork = n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    if (Q.isrm()) {
      int lda = Q.stepi();
      cunglq(&n,&m,&n,LAP_Complex(Q.ptr()),&lda,LAP_Complex(Qbeta.cptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"cunglq");
    } else {
      int lda = Q.stepj();
      Vector<complex<float> > conjbeta = Qbeta.Conjugate();
      cungqr(&m,&n,&n,LAP_Complex(Q.ptr()),&lda,LAP_Complex(conjbeta.cptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"cungqr");
    }
  }
#endif
#endif
  template <class T> void GetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& Qbeta)
  {
    // Extract the Q matrix from a combined QRx matrix
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Qbeta.size() == Q.rowsize());
    TMVAssert(Q.isrm() || Q.iscm());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
#ifdef LAP
    LapGetQFromQR(Q,Qbeta);
#else
    NonLapGetQFromQR(Q,Qbeta);
#endif
  }

  /*
  template <class T> void GetQFromQR(
      const MatrixView<T>& Q, const vector<GenUpperTriMatrix<T>*>& ZZ)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
#ifdef XDEBUG
    Matrix<T> Q0 = Q;
    Matrix<T> Q2 = Q;
    Vector<T> Qbeta(Q.rowsize());
    for(size_t k=0,i=0;k<ZZ.size();i+=ZZ[k]->size(),++k)
      Qbeta.SubVector(i,i+ZZ[k]->size()) = ZZ[k]->diag().Conjugate();
    NonBlockGetQFromQR(Q2.View(),Qbeta);
#endif
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    UpperTriMatrixViewOf(Q).Zero();
    for(size_t k=ZZ.size()-1,j2=N;j2>0;--k) {
      const GenUpperTriMatrix<T>& Z = *(ZZ[k]);
      size_t j1 = j2-Z.size();
      MatrixView<T> Y = Q.SubMatrix(j1,M,j1,j2);
      BlockHouseholder_Unpack(Y,Z,Q.SubMatrix(j1,M,j2,N));
      if (k==0) { TMVAssert(j1==0); }
      j2 = j1;
    }
#ifdef XDEBUG
    if (Norm(Q-Q2) > 0.001*Norm(Q2)) {
      cerr<<"BlockGetQ: Q = "<<Type(Q)<<"  "<<Q0<<endl;
      cerr<<"Qbeta = "<<Qbeta<<endl;
      cerr<<"-> Q = "<<Q<<endl;
      cerr<<"NonBlock Q = "<<Q2<<endl;
      abort(); 
    }
#endif
  }*/

  //
  // QR Decompose
  //

  template <class T> void NonBlockQR_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& Qbeta, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    // Decompose A into A = Q R 
    // where Q is unitary, and R is upper triangular
    // Q and R are stored in the same matrix (output of A), 
    // with the beta's for the Householder matrices returned in Qbeta.
    const size_t M = A.colsize();
    const size_t N = A.rowsize();

    for(size_t j=0;j<N;++j) {
      // Apply the Householder Reflection for this column
      Qbeta(j) = Householder_Reflect(A.SubMatrix(j,M,j,N),det);
    }
  }

  template <class T> void RecursiveQR_Decompose(
      const MatrixView<T>& A, const UpperTriMatrixView<T>& Z, T& det,
      bool makeZ)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.ct() == NonConj);
    // This is very similar to the BlockHouseholder_MakeZ function
    // in Householder.cpp.  The difference is the addition of the 
    // Householder_Reflects.
    // The makeZ parameter should be set to true if you want the Z
    // matrix to be correct on output.  If you don't need the Z matrix
    // after this call, setting makeZ to false will speed it up slightly.
    // (In either case, the diagonal of Z is correctly set to be 
    // Qbeta.Conjugate().)
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
      const MatrixView<T>& A, const VectorView<T>& Qbeta, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> A2 = A;
    Vector<T> Qbeta2(A.rowsize());
    T det2 = det;
    NonBlockQR_Decompose(A2.View(),Qbeta2.View(),det2);
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
	Qbeta.SubVector(j1,j2) = Z.diag().Conjugate();
      } else {
	for(size_t j=j1;j<j2;++j) {
	  T b = Householder_Reflect(A.SubMatrix(j,M,j,j2),det);
	  Qbeta(j) = b;
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
    if (Norm(A-A2) > 0.001*Norm(A2)) {
      cerr<<"BlockQR_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"Qbeta = "<<Qbeta<<endl;
      Matrix<T> A2 = A0;
      Vector<T> Qbeta2(Qbeta.size());
      T det2(0);
      NonBlockQR_Decompose(A2.View(),Qbeta2.View(),det2);
      cerr<<"NonBlock "<<A2<<endl;
      cerr<<"Qbeta = "<<Qbeta2<<endl;
      abort(); 
    }
#endif
  }

  template <class T> inline void NonLapQR_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& Qbeta, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

    if (A.rowsize() > QR_BLOCKSIZE)
      BlockQR_Decompose(A,Qbeta,det);
    else if (RecursiveQR) {
      UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(A.rowsize());
      RecursiveQR_Decompose(A,Z.View(),det,false);
      Qbeta = Z.diag().Conjugate();
    }
    else
      NonBlockQR_Decompose(A,Qbeta,det);
  }

#ifdef LAP
  template <class T> inline void LapQR_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& Qbeta, T& det)
  { NonLapQR_Decompose(A,Qbeta,det); }
  template <> void LapQR_Decompose(const MatrixView<double>& A,
      const VectorView<double>& Qbeta, double& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lwork = 2*n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    if (A.isrm()) {
      int lda = A.stepi();
      dgelqf(&n,&m,A.ptr(),&lda,Qbeta.ptr(),work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dgelqf");
    } else {
      int lda = A.stepj();
      dgeqrf(&m,&n,A.ptr(),&lda,Qbeta.ptr(),work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dgeqrf");
    }
    if (det) for(size_t i=0;i<Qbeta.size();++i) if (Qbeta(i) != 0.) det = -det;
  }
  template <> void LapQR_Decompose(
      const MatrixView<complex<double> >& A,
      const VectorView<complex<double> >& Qbeta, complex<double>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lwork = 2*n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    if (A.isrm()) {
      int lda = A.stepi();
      zgelqf(&n,&m,LAP_Complex(A.ptr()),&lda,LAP_Complex(Qbeta.ptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"zgelqf");
    } else {
      int lda = A.stepj();
      zgeqrf(&m,&n,LAP_Complex(A.ptr()),&lda,LAP_Complex(Qbeta.ptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"zgeqrf");
      Qbeta.ConjugateSelf();
    }
    if (det!=double(0)) {
      for(size_t i=0;i<Qbeta.size();++i)
	if (imag(Qbeta(i)) != 0.) 
	  det *= -conj(Qbeta(i)*Qbeta(i))/norm(Qbeta(i));
	else if (real(Qbeta(i)) != 0.)
	  det = -det;
    }
  }
#ifndef NOFLOAT
  template <> void LapQR_Decompose(const MatrixView<float>& A,
      const VectorView<float>& Qbeta, float& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lwork = 2*n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    if (A.isrm()) {
      int lda = A.stepi();
      sgelqf(&n,&m,A.ptr(),&lda,Qbeta.ptr(),work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sgelqf");
    } else {
      int lda = A.stepj();
      sgeqrf(&m,&n,A.ptr(),&lda,Qbeta.ptr(),work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sgeqrf");
    }
    if (det) for(size_t i=0;i<Qbeta.size();++i) if (Qbeta(i) != 0.) det = -det;
  }
  template <> void LapQR_Decompose(
      const MatrixView<complex<float> >& A,
      const VectorView<complex<float> >& Qbeta, complex<float>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lwork = 2*n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    if (A.isrm()) {
      int lda = A.stepi();
      cgelqf(&n,&m,LAP_Complex(A.ptr()),&lda,LAP_Complex(Qbeta.ptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"cgelqf");
    } else {
      int lda = A.stepj();
      cgeqrf(&m,&n,LAP_Complex(A.ptr()),&lda,LAP_Complex(Qbeta.ptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"cgeqrf");
      Qbeta.ConjugateSelf();
    }
    if (det!=float(0)) {
      for(size_t i=0;i<Qbeta.size();++i)
	if (imag(Qbeta(i)) != 0.) 
	  det *= -conj(Qbeta(i)*Qbeta(i))/norm(Qbeta(i));
	else if (real(Qbeta(i)) != 0.)
	  det = -det;
    }
  }
#endif
#endif
  template <class T> void QR_Decompose(
      const MatrixView<T>& A, const VectorView<T>& Qbeta, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    if (A.rowsize() > 0) {
#ifdef LAP
      if (Qbeta.step()==1 && (A.isrm() || A.iscm()))
	LapQR_Decompose(A,Qbeta,det);
      else 
#endif
	NonLapQR_Decompose(A,Qbeta,det);
    }
  }

  /*
  template <class T> void QR_Decompose(
      const MatrixView<T>& A, vector<GenUpperTriMatrix<T>*>& ZZ, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> A2 = A;
    Vector<T> Qbeta2(A.rowsize());
    T det2 = det;
    NonBlockQR_Decompose(A2.View(),Qbeta2.View(),det2);
#endif
    const size_t M = A.colsize();
    const size_t N = A.rowsize();

    for(size_t j1=0;j1<N;) {
      size_t j2 = min(N,j1+QR_BLOCKSIZE);
      MatrixView<T> A1 = A.SubMatrix(j1,M,j1,j2);
      UpperTriMatrix<T,NonUnitDiag,ColMajor>* Z = 
	new UpperTriMatrix<T,NonUnitDiag,ColMajor>(A1.rowsize());
      if (RecursiveQR) {
	RecursiveQR_Decompose(A1,Z->View(),det,true);
      } else {
	for(size_t j=j1;j<j2;++j) {
	  T b = Householder_Reflect(A.SubMatrix(j,M,j,j2),det);
	  BlockHouseholder_Augment(A.SubMatrix(j1,M,j1,j+1),
	      Z->SubTriMatrix(0,j+1),CONJ(b));
	}
      }
      BlockHouseholder_LDiv(A1,*Z,A.SubMatrix(j1,M,j2,N));
      ZZ.push_back(Z);
      j1 = j2;
    }
#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*Norm(A2)) {
      cerr<<"BlockQR_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      UpperTriMatrix<T> Z(A.rowsize());
      for(size_t k=0,j=0;k<ZZ.size();j+=ZZ[k]->size(),k++)
	Z.SubTriMatrix(j,j+ZZ[k]->size()) = *ZZ[k];
      cerr<<"Z = "<<Z<<endl;
      cerr<<"NonBlock "<<A2<<endl;
      cerr<<"Qbeta = "<<Qbeta2<<endl;
      abort(); 
    }
#endif
  }*/

  //
  // QRP Decompose
  //

  template <class T> void NonBlockQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& Qbeta, size_t* P, T& det)
  {
    // Decompose A into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (output of A), 
    // with the beta's for the Householder matrices returned in Qbeta.
    
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(Qbeta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    const size_t M = A.colsize();
    const size_t N = A.rowsize();
    const RealType(T) RT0(0);

    // Keep track of the norm of each column
    // When considering column j, these are actually just the norm
    // of each column from j:M, not 0:M.
    Vector<RealType(T)> colnormsq(N);
    for(size_t j=0;j<N;++j) colnormsq(j) = NormSq(A.col(j));
    RealType(T) anormsq = colnormsq.SumElements();
    RealType(T) thresh = RealType(T)(N) * SQR(Epsilon<T>()) * anormsq;
    //cerr<<"thresh = "<<thresh<<endl;
    // Set to 0 any diag element whose norm is < epsilon * |A|
    RealType(T) recalcthresh=RT0;
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
	  Qbeta.SubVector(j,N).Zero();
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
      Qbeta(j) = Householder_Reflect(A.SubMatrix(j,M,j,N),det);

      // And update the norms for use with the next column
      CVIt<T,Step,NonConj> Ajk = A.row(j,j+1,N).begin();
      for(size_t k=j+1;k<N;++k,++Ajk) {
	colnormsq(k) -= NORM(*Ajk);
      }
    }
  }

  template <class T> void StrictBlockQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& Qbeta, size_t* P, T& det)
  {
    // Decompose A (input as A) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (A), with the beta's for
    // the Householder matrices returned in Qbeta.
    
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(Qbeta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    const size_t M = A.colsize();
    const size_t N = A.rowsize();
    const RealType(T) RT0(0);

    Vector<RealType(T)> colnormsq(N);
    for(size_t j=0;j<N;++j) colnormsq(j) = NormSq(A.col(j));
    RealType(T) anormsq = colnormsq.SumElements();
    RealType(T) thresh = RealType(T)(N) * SQR(Epsilon<T>()) * anormsq;
    RealType(T) recalcthresh=RT0;

    Matrix<T,RowMajor> ZYtA(min(QRP_BLOCKSIZE,M),N);
    // We keep track of the product ZYtA(j1:M,j1:N) [ stored as ZYtA ]
    // since this is the product that we need.  We update this one 
    // row at a time.
    
    for(size_t j1=0;j1<N;) {
      size_t j2 = min(N,j1+QRP_BLOCKSIZE);

      for(size_t j=j1,jmj1=0; j<j2; ++j,++jmj1) {
	size_t jpiv;
	RealType(T) maxnormsq = colnormsq.SubVector(j,N).MaxElement(&jpiv);
	if (recalcthresh==RT0) recalcthresh = 4*SqrtEpsilon<T>() * maxnormsq;

	if (maxnormsq < recalcthresh) {
	  for(size_t k=j;k<N;++k) colnormsq(k) = NormSq(A.col(k,j,M));
	  recalcthresh = RT0;
	  j2 = j;
	} else if (maxnormsq < thresh) {
	  if (j==j1) {
	    // If first in set, then just zero the rest out and 
	    // indicate that we are done.
	    A.SubMatrix(j,M,j,N).Zero();
	    Qbeta.SubVector(j,N).Zero();
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
	  Qbeta(j) = Householder_Reflect(A.col(j,j,M),det);

	  // Update ZYtA:
	  if (Qbeta(j) != T(0)) {
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
	    ZYtA.row(jmj1,j1,N) *= Qbeta(j);
	  } else ZYtA.row(jmj1,j1,N).Zero();

	  // Update row j of the rest of the matrix:
	  // A(j,j+1:N) -= (Y ZYtA)(j,j+1:N) = Y(j,j1:j+1) ZYtA(j1:j+1,j+1:N)
	  VectorView<T> Arowj = A.row(j,j+1,N);
	  Arowj -= A.row(j,j1,j)*ZYtA.SubMatrix(0,jmj1,j+1,N);
	  Arowj -= ZYtA.row(jmj1,j+1,N);

	  // Update the colnormsq values
	  CVIt<T,Step,NonConj> Ajk = Arowj.begin();
	  //for(size_t k=j+1;k<N;++k,++Ajk) colnormsq(k) -= tmv::NORM(*Ajk);
	  VIt<RealType(T),Unit,NonConj> cnk = colnormsq.begin()+j+1;
	  for(size_t k=N-j-1;k>0;--k,++cnk,++Ajk) *cnk -= tmv::NORM(*Ajk);
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
    GetQFromQR(Q.View(),Qbeta);
    Matrix<T> AA = Q*UpperTriMatrixViewOf(A);
    AA.ReversePermuteCols(P);
    if (Norm(AA-A0) > 1.e-12*Norm(A0)) {
      cerr<<"StrictBlockQRP_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"Qbeta = "<<Qbeta<<endl;
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
      const MatrixView<T>& A, const VectorView<T>& Qbeta, size_t* P, T& det)
  {
    // Decompose A (input as A) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (A), with the beta's for
    // the Householder matrices returned in Qbeta.
    //
    // This loose version doesn't sort the diagonal of R exactly.
    // It only sorts them enough to make sure the 0's fall at the end.
    
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(Qbeta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    const size_t M = A.colsize();
    const size_t N = A.rowsize();
    const RealType(T) RT0(0);

    Vector<RealType(T)> colnormsq(N);
    for(size_t j=0;j<N;++j) colnormsq(j) = NormSq(A.col(j));
    //cerr<<"colnormsq = "<<colnormsq<<endl;
    RealType(T) anormsq = colnormsq.SumElements();
    //cerr<<"anormsq = "<<anormsq<<endl;
    //cerr<<"eps = "<<Epsilon<T>()<<endl;
    RealType(T) thresh = RealType(T)(N) * SQR(Epsilon<T>()) * anormsq;
    //cerr<<"thresh = "<<thresh<<endl;
    RealType(T) recalcthresh(0);

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
	Qbeta.SubVector(j1,N).Zero();
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

#ifdef XDEBUG
      RealType(T) recalcthresh = 0.01*maxnormsq;
#else
      RealType(T) recalcthresh = RealType(T)(N)*SqrtEpsilon<T>()*maxnormsq;
#endif
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
	  Qbeta(j) = Householder_Reflect(A.SubMatrix(j,M,j,origj2),det);

	  // Update Z:
	  BlockHouseholder_Augment(A.SubMatrix(j1,M,j1,j+1),
	      Z.SubTriMatrix(0,j-j1+1),CONJ(Qbeta(j)));

	  // Update the colnormsq values within this block
	  // (No need to go all the way to origj2, since the j2..origj2
	  // columns are those with low norm already - we don't need
	  // those values until we recalculate them from scratch anyway.)
	  VectorView<T> Arowj = A.row(j,j+1,j2);
	  CVIt<T,Step,NonConj> Ajk = Arowj.begin();
	  for(size_t k=j+1;k<j2;++k,++Ajk) colnormsq(k) -= tmv::NORM(*Ajk);
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
    GetQFromQR(Q.View(),Qbeta);
    Matrix<T> AA = Q*UpperTriMatrixViewOf(A);
    AA.ReversePermuteCols(P);
    if (Norm(AA-A0) > 1.e-12*Norm(A0)) {
      cerr<<"LooseBlockQRP_Decompose: \n";
      cerr<<"A = "<<Type(A)<<endl;
      if (N < 100) {
	cerr<<"  "<<A0<<endl;
	cerr<<"-> "<<A<<endl;
	cerr<<"Qbeta = "<<Qbeta<<endl;
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
      const VectorView<T>& Qbeta, size_t* P, T& det)
  {
    // Decompose A (input as A) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (A), with the beta's for
    // the Householder matrices returned in Qbeta.

    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(Qbeta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

    if (A.rowsize() > QRP_BLOCKSIZE)
      if (StrictQRP)
	StrictBlockQRP_Decompose(A,Qbeta,P,det);
      else
	LooseBlockQRP_Decompose(A,Qbeta,P,det);
    else
      NonBlockQRP_Decompose(A,Qbeta,P,det);
  }

#ifdef XLAP
  template <class T> inline void LapQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& Qbeta, size_t* P, T& det)
  { NonLapQR_Decompose(A,Qbeta,det); }
  template <> void LapQRP_Decompose(const MatrixView<double>& A,
      const VectorView<double>& Qbeta, size_t* P, double& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lap_p[n];
    int lwork = 3*n;
    double* work = LAP_DWork(lwork);
    int info;
    int lda = A.stepj();
    char cc = 'F';
    double thresh = Epsilon<double>()*dlange(&cc,&m,&n,A.ptr(),&lda,0);
    dgeqpf(&m,&n,A.ptr(),&lda,lap_p,Qbeta.ptr(),work,&info);
    LAP_Results(info,"dgeqpf");
    for(size_t i=0;i<A.rowsize();++i) {
      if (abs(A(i,i)) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	Qbeta.SubVector(i,Qbeta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<Qbeta.size();++i) {
      if (det) if (Qbeta(i) != 0.) det = -det;
      P[i] = lap_p[i]-1;
      if (det) if (P[i] != i) det = -det;
    }
  }
  template <> void LapQRP_Decompose(
      const MatrixView<complex<double> >& A,
      const VectorView<complex<double> >& Qbeta, size_t* P,
      complex<double>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
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
    zgeqpf(&m,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(Qbeta.ptr()),
	LAP_Complex(work),rwork,&info);
    LAP_Results(info,"zgeqpf");
    Qbeta.ConjugateSelf();
    for(size_t i=0;i<A.rowsize();++i) {
      TMVAssert(abs(imag(A(i,i))) < Epsilon<double>());
      if (abs(real(A(i,i))) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	Qbeta.SubVector(i,Qbeta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<Qbeta.size();++i) {
      if (det!=double(0)) {
	if (imag(Qbeta(i)) != 0.) 
	  det *= -conj(Qbeta(i)*Qbeta(i))/norm(Qbeta(i));
	else if (real(Qbeta(i)) != 0.)
	  det = -det;
      }
      P[i] = lap_p[i]-1;
      if (det!=double(0)) if (P[i] != i) det = -det;
    }
  }
#ifndef NOFLOAT
  template <> void LapQRP_Decompose(const MatrixView<float>& A,
      const VectorView<float>& Qbeta, size_t* P, float& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lap_p[n];
    int lwork = 3*n;
    float* work = LAP_SWork(lwork);
    int info;
    int lda = A.stepj();
    char cc = 'F';
    float thresh = Epsilon<float>()*slange(&cc,&m,&n,A.ptr(),&lda,0);
    sgeqpf(&m,&n,A.ptr(),&lda,lap_p,Qbeta.ptr(),work,&info);
    LAP_Results(info,"sgeqpf");
    for(size_t i=0;i<A.rowsize();++i) {
      if (abs(A(i,i)) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	Qbeta.SubVector(i,Qbeta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<Qbeta.size();++i) {
      if (det) if (Qbeta(i) != 0.) det = -det;
      P[i] = lap_p[i]-1;
      if (det) if (P[i] != i) det = -det;
    }
  }
  template <> void LapQRP_Decompose(
      const MatrixView<complex<float> >& A,
      const VectorView<complex<float> >& Qbeta, size_t* P,
      complex<float>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
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
    cgeqpf(&m,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(Qbeta.ptr()),
	LAP_Complex(work),rwork,&info);
    LAP_Results(info,"cgeqpf");
    for(size_t i=0;i<A.rowsize();++i) {
      TMVAssert(abs(imag(A(i,i))) < Epsilon<double>());
      if (abs(real(A(i,i))) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	Qbeta.SubVector(i,Qbeta.size()).Zero();
	break;
      }
    }
    Qbeta.ConjugateSelf();
    for(size_t i=0;i<Qbeta.size();++i) {
      if (det!=float(0)) {
	if (imag(Qbeta(i)) != 0.) 
	  det *= -conj(Qbeta(i)*Qbeta(i))/norm(Qbeta(i));
	else if (real(Qbeta(i)) != 0.)
	  det = -det;
      }
      P[i] = lap_p[i]-1;
      if (det!=float(0)) if (P[i] != i) det = -det;
    }
  }
#endif
  template <class T> inline void NewLapQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& Qbeta, size_t* P, T& det)
  { NonLapQR_Decompose(A,Qbeta,det); }
  template <> void NewLapQRP_Decompose(const MatrixView<double>& A,
      const VectorView<double>& Qbeta, size_t* P, double& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lap_p[n];
    int lwork = 3*n+1;
    double* work = LAP_DWork(lwork);
    int info;
    int lda = A.stepj();
    char cc = 'F';
    double thresh = Epsilon<double>()*dlange(&cc,&m,&n,A.ptr(),&lda,0);
    dgeqp3(&m,&n,A.ptr(),&lda,lap_p,Qbeta.ptr(),work,&lwork,&info);
    LAP_Results(info,"dgeqp3");
    for(size_t i=0;i<A.rowsize();++i) {
      if (abs(A(i,i)) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	Qbeta.SubVector(i,Qbeta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<Qbeta.size();++i) {
      if (det) if (Qbeta(i) != 0.) det = -det;
      P[i] = lap_p[i]-1;
      if (det) if (P[i] != i) det = -det;
    }
  }
  template <> void NewLapQRP_Decompose(
      const MatrixView<complex<double> >& A,
      const VectorView<complex<double> >& Qbeta, size_t* P,
      complex<double>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
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
    zgeqp3(&m,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(Qbeta.ptr()),
	LAP_Complex(work),&lwork,rwork,&info);
    LAP_Results(info,"zgeqp3");
    Qbeta.ConjugateSelf();
    for(size_t i=0;i<A.rowsize();++i) {
      TMVAssert(abs(imag(A(i,i))) < Epsilon<double>());
      if (abs(real(A(i,i))) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	Qbeta.SubVector(i,Qbeta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<Qbeta.size();++i) {
      if (det!=double(0)) {
	if (imag(Qbeta(i)) != 0.) 
	  det *= -conj(Qbeta(i)*Qbeta(i))/norm(Qbeta(i));
	else if (real(Qbeta(i)) != 0.)
	  det = -det;
      }
      P[i] = lap_p[i]-1;
      if (det!=double(0)) if (P[i] != i) det = -det;
    }
  }
#ifndef NOFLOAT
  template <> void NewLapQRP_Decompose(const MatrixView<float>& A,
      const VectorView<float>& Qbeta, size_t* P, float& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int lap_p[n];
    int lwork = 3*n+1;
    float* work = LAP_SWork(lwork);
    int info;
    int lda = A.stepj();
    char cc = 'F';
    float thresh = Epsilon<float>()*slange(&cc,&m,&n,A.ptr(),&lda,0);
    sgeqp3(&m,&n,A.ptr(),&lda,lap_p,Qbeta.ptr(),work,&lwork,&info);
    LAP_Results(info,"sgeqp3");
    for(size_t i=0;i<A.rowsize();++i) {
      if (abs(A(i,i)) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	Qbeta.SubVector(i,Qbeta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<Qbeta.size();++i) {
      if (det) if (Qbeta(i) != 0.) det = -det;
      P[i] = lap_p[i]-1;
      if (det) if (P[i] != i) det = -det;
    }
  }
  template <> void NewLapQRP_Decompose(
      const MatrixView<complex<float> >& A,
      const VectorView<complex<float> >& Qbeta, size_t* P,
      complex<float>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
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
    cgeqp3(&m,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(Qbeta.ptr()),
	LAP_Complex(work),&lwork,rwork,&info);
    LAP_Results(info,"cgeqp3");
    for(size_t i=0;i<A.rowsize();++i) {
      TMVAssert(abs(imag(A(i,i))) < Epsilon<double>());
      if (abs(real(A(i,i))) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	Qbeta.SubVector(i,Qbeta.size()).Zero();
	break;
      }
    }
    Qbeta.ConjugateSelf();
    for(size_t i=0;i<Qbeta.size();++i) {
      if (det!=float(0)) {
	if (imag(Qbeta(i)) != 0.) 
	  det *= -conj(Qbeta(i)*Qbeta(i))/norm(Qbeta(i));
	else if (real(Qbeta(i)) != 0.)
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
      const VectorView<T>& Qbeta, size_t* P, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == Qbeta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

    if (A.rowsize() > 0) {
      // Neither LAP version of QRP seems to do any actual pivoting.
      // At least for the Intel MKL version.
      // So until I figure out why and how to get them to work, leave
      // this is XLAP (ie. ignore the LAP routines.
#ifdef XLAP
      //cerr<<"A.stor = "<<Text(A.stor())<<endl;
      if (A.iscm())
	//LapQRP_Decompose(A,Qbeta,P,det);
	NewLapQRP_Decompose(A,Qbeta,P,det);
      else {
#ifdef TMVDEBUG
	cout<<"Lap QRDecomp: A, Qbeta are wrong step:\n";
	cout<<"A isrm = "<<A.isrm()<<", stepi = "<<A.stepi()<<endl;
	cout<<"Qbeta step = "<<Qbeta.step()<<endl;
#endif
	NonLapQRP_Decompose(A,Qbeta,P,det);
      }
#else
      NonLapQRP_Decompose(A,Qbeta,P,det);
#endif
    }
  }

  //
  // QR Decompose - Unpacked
  //

  template <class T> void QR_Decompose(
      const MatrixView<T>& Q, const MatrixView<T>& R, T& det)
  {
    // Decompose A (input as Q) into A = Q R 
    // where Q is unitary and R is upper triangular

    const size_t N = Q.rowsize();
    TMVAssert(Q.colsize() >= N);
    TMVAssert(R.colsize() == N);
    TMVAssert(R.rowsize() == N);
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(R.ct() == NonConj);

    Vector<T> Qbeta(N);
    QR_Decompose(Q,Qbeta.View(),det);
    R = UpperTriMatrixViewOf(Q);
    GetQFromQR(Q,Qbeta.View());
  }

  template <class T> void QRP_Decompose(
      const MatrixView<T>& Q, const MatrixView<T>& R, size_t* P, T& det)
  {
    // Decompose A (input as Q) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation

    const size_t N = Q.rowsize();
    TMVAssert(Q.colsize() >= N);
    TMVAssert(R.colsize() == N);
    TMVAssert(R.rowsize() == N);
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(R.ct() == NonConj);

    Vector<T> Qbeta(N);
    QRP_Decompose(Q,Qbeta.View(),P,det);
    R = UpperTriMatrixViewOf(Q);
    GetQFromQR(Q,Qbeta.View());
  }

  //
  // Packed Q - LDivEq
  //

  template <class T1, class T2> void NonBlockQ_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.colsize() == m.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    // Solve Q x = m in place 
    // where Q is stored as Householder vectors along with Qbeta
    //
    // Q is H0t H1t ... H_N-1t
    // So x = H_N-1 .. H1 H0 m
    //
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    for(size_t j=0;j<N;++j) if (Qbeta(j) != T1(0)) {
      Householder_LMult(Q.col(j,j+1,M),Qbeta(j),m.Rows(j,M));
    }
  }

  template <class T1, class T2> void BlockQ_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.colsize() == m.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T2> m0 = m;
    Matrix<T2> m2 = m;
    NonBlockQ_LDivEq(Q,Qbeta,m2.View());
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
      BlockHouseholder_MakeZ(Y,Z,Qbeta.SubVector(j1,j2));
      BlockHouseholder_LDiv(Y,Z,m.Rows(j1,M));
      j1 = j2;
    }
#ifdef XDEBUG
    if (Norm(m-m2) > 0.001*Norm(m2)) {
      cerr<<"BlockQ_LDivEq: Q = "<<Type(Q)<<"  "<<Q<<endl;
      cerr<<"Qbeta = "<<Qbeta<<endl;
      cerr<<"m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"NonBlock m = "<<m2<<endl;
      abort(); 
    }
#endif
  }

  template <class T1, class T2> inline void NonLapQ_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.colsize() == m.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

    if (Q.rowsize() > QR_BLOCKSIZE && m.rowsize() > QR_BLOCKSIZE)
      BlockQ_LDivEq(Q,Qbeta,m);
    else
      NonBlockQ_LDivEq(Q,Qbeta,m);
  }
#ifdef LAP
  template <class T1, class T2> inline void LapQ_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  { NonLapQ_LDivEq(Q,Qbeta,m); }
  template <> void LapQ_LDivEq(const GenMatrix<double>& Q,
      const GenVector<double>& Qbeta, const MatrixView<double>& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.colsize() == x.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
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
	  const_cast<double*>(Qbeta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dormlq");
    } else {
      char trans = x.isrm() ? 'N' : 'T';
      int ldq = Q.stepj();
      dormqr(&side,&trans,&m,&n,&k,const_cast<double*>(Q.cptr()),&ldq,
	  const_cast<double*>(Qbeta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dormqr");
    }
  }
  template <> void LapQ_LDivEq(
      const GenMatrix<complex<double> >& Q,
      const GenVector<complex<double> >& Qbeta,
      const MatrixView<complex<double> >& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.colsize() == x.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
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
	    LAP_Complex(Qbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
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
	    LAP_Complex(Qbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
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
	Vector<complex<double> > conjbeta = Qbeta.Conjugate();
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
	Vector<complex<double> > conjbeta = Qbeta.Conjugate();
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
      const GenVector<float>& Qbeta, const MatrixView<float>& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.colsize() == x.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
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
	  const_cast<float*>(Qbeta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sormlq");
    } else {
      char trans = x.isrm() ? 'N' : 'T';
      int ldq = Q.stepj();
      sormqr(&side,&trans,&m,&n,&k,const_cast<float*>(Q.cptr()),&ldq,
	  const_cast<float*>(Qbeta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sormqr");
    }
  }
  template <> void LapQ_LDivEq(
      const GenMatrix<complex<float> >& Q,
      const GenVector<complex<float> >& Qbeta,
      const MatrixView<complex<float> >& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.colsize() == x.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
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
	    LAP_Complex(Qbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
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
	    LAP_Complex(Qbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
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
	Vector<complex<float> > conjbeta = Qbeta.Conjugate();
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
	Vector<complex<float> > conjbeta = Qbeta.Conjugate();
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

  template <class T1, class T2> void Q_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Qbeta.size() == Q.rowsize());
    TMVAssert(m.colsize() == Q.colsize());
    TMVAssert(Q.isrm() || Q.iscm());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    if (m.colsize() > 0 && m.rowsize() > 0) {
#ifdef LAP
      if ( m.isrm() || m.iscm() )
	LapQ_LDivEq(Q,Qbeta,m);
      else
#endif
	NonLapQ_LDivEq(Q,Qbeta,m);
    }
  }

  /*
  template <class T1, class T2> void Q_LDivEq(
      const GenMatrix<T1>& Q, const vector<GenUpperTriMatrix<T1>*>& ZZ,
      const MatrixView<T2>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.colsize() == m.colsize());
    TMVAssert(Q.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T2> m0 = m;
    Matrix<T2> m2 = m;
    Vector<T1> Qbeta(Q.rowsize());
    for(size_t k=0,i=0;k<ZZ.size();i+=ZZ[k]->size(),++k)
      Qbeta.SubVector(i,i+ZZ[k]->size()) = ZZ[k]->diag().Conjugate();
    NonBlockQ_LDivEq(Q,Qbeta,m2.View());
#endif

    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    for(size_t k=0,j1=0;k<ZZ.size();++k) {
      const GenUpperTriMatrix<T1>& Z = *(ZZ[k]);
      size_t j2 = j1+Z.size();
      ConstMatrixView<T1> Y = Q.SubMatrix(j1,M,j1,j2);
      BlockHouseholder_LDiv(Y,Z,m.Rows(j1,M));
      if (k==ZZ.size()-1) { TMVAssert(j2==N); }
      j1 = j2;
    }
#ifdef XDEBUG
    if (Norm(m-m2) > 0.001*Norm(m2)) {
      cerr<<"BlockQ_LDivEq: Q = "<<Type(Q)<<"  "<<Q<<endl;
      cerr<<"Qbeta = "<<Qbeta<<endl;
      cerr<<"m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"NonBlock m = "<<m2<<endl;
      abort(); 
    }
#endif
  }*/

  //
  // Packed Q - RDivEq
  //

  template <class T1, class T2> void NonBlockQ_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Qbeta.size() == Q.rowsize());
    TMVAssert(m.rowsize() == Q.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    // Solve x Q = m in place
    // where Q is stored as Householder vectors along with Qbeta
    //
    // x = m Qt 
    // Qt is H_N-1 H_N-2 ... H1 H0
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    for(int j=N-1;j>=0;--j) if (Qbeta(j) != T1(0)) {
      Householder_RMult(Q.col(j,j+1,M),Qbeta(j),m.Cols(j,M));
    }
  }

  template <class T1, class T2> void BlockQ_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.colsize() == m.rowsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T2> m0 = m;
    Matrix<T2> m2 = m;
    NonBlockQ_RDivEq(Q,Qbeta,m2.View());
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
      BlockHouseholder_MakeZ(Y,Z,Qbeta.SubVector(j1,j2));
      BlockHouseholder_RDiv(Y,Z,m.Cols(j1,M));
      j2 = j1;
    }
#ifdef XDEBUG
    if (Norm(m-m2) > 0.001*Norm(m2)) {
      cerr<<"BlockQ_RDivEq: Q = "<<Type(Q)<<"  "<<Q<<endl;
      cerr<<"Qbeta = "<<Qbeta<<endl;
      cerr<<"m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"NonBlock m = "<<m2<<endl;
      abort(); 
    }
#endif
  }

  template <class T1, class T2> inline void NonLapQ_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.colsize() == m.rowsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

    if (Q.rowsize() > QR_BLOCKSIZE && m.colsize() > QR_BLOCKSIZE)
      BlockQ_RDivEq(Q,Qbeta,m);
    else
      NonBlockQ_RDivEq(Q,Qbeta,m);
  }

#ifdef LAP
  template <class T1, class T2> inline void LapQ_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  { NonLapQ_RDivEq(Q,Qbeta,m); }
  template <> void LapQ_RDivEq(const GenMatrix<double>& Q,
      const GenVector<double>& Qbeta, const MatrixView<double>& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Qbeta.size() == Q.rowsize());
    TMVAssert(x.rowsize() == Q.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
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
	  const_cast<double*>(Qbeta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dormlq");
    } else {
      char trans = x.isrm() ? 'N' : 'T';
      int ldq = Q.stepj();
      dormqr(&side,&trans,&m,&n,&k,const_cast<double*>(Q.cptr()),&ldq,
	  const_cast<double*>(Qbeta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dormqr");
    }
  }
  template <> void LapQ_RDivEq(const GenMatrix<complex<double> >& Q,
      const GenVector<complex<double> >& Qbeta,
      const MatrixView<complex<double> >& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Qbeta.size() == Q.rowsize());
    TMVAssert(x.rowsize() == Q.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
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
	    LAP_Complex(Qbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
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
	    LAP_Complex(Qbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
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
	Vector<complex<double> > conjbeta = Qbeta.Conjugate();
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
	Vector<complex<double> > conjbeta = Qbeta.Conjugate();
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
      const GenVector<float>& Qbeta, const MatrixView<float>& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Qbeta.size() == Q.rowsize());
    TMVAssert(x.rowsize() == Q.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
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
	  const_cast<float*>(Qbeta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sormlq");
    } else {
      char trans = x.isrm() ? 'N' : 'T';
      int ldq = Q.stepj();
      sormqr(&side,&trans,&m,&n,&k,const_cast<float*>(Q.cptr()),&ldq,
	  const_cast<float*>(Qbeta.cptr()),x.ptr(),&ldx,
	  work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sormqr");
    }
  }
  template <> void LapQ_RDivEq(const GenMatrix<complex<float> >& Q,
      const GenVector<complex<float> >& Qbeta,
      const MatrixView<complex<float> >& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Qbeta.size() == Q.rowsize());
    TMVAssert(x.rowsize() == Q.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
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
	      LAP_Complex(Qbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
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
	      LAP_Complex(Qbeta.cptr()),LAP_Complex(x.ptr()),&ldx,
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
	  Vector<complex<float> > conjbeta = Qbeta.Conjugate();
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
	  Vector<complex<float> > conjbeta = Qbeta.Conjugate();
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

  template <class T1, class T2> void Q_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Qbeta.size() == Q.rowsize());
    TMVAssert(m.rowsize() == Q.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    if (m.colsize() > 0 && m.rowsize() > 0) {
#ifdef LAP
      if (m.isrm() || m.iscm())
	LapQ_RDivEq(Q,Qbeta,m);
      else
#endif
	NonLapQ_RDivEq(Q,Qbeta,m);
    }
  }

  /*
  template <class T1, class T2> void Q_RDivEq(
      const GenMatrix<T1>& Q, const vector<GenUpperTriMatrix<T1>*>& ZZ,
      const MatrixView<T2>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.colsize() == m.rowsize());
    TMVAssert(Q.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T2> m0 = m;
    Matrix<T2> m2 = m;
    Vector<T1> Qbeta(Q.rowsize());
    for(size_t k=0,i=0;k<ZZ.size();i+=ZZ[k]->size(),++k)
      Qbeta.SubVector(i,i+ZZ[k]->size()) = ZZ[k]->diag().Conjugate();
    NonBlockQ_RDivEq(Q,Qbeta,m2.View());
#endif

    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    for(size_t k=ZZ.size()-1,j2=N;j2>0;--k) {
      const GenUpperTriMatrix<T1>& Z = *(ZZ[k]);
      size_t j1 = j2-Z.size();
      ConstMatrixView<T1> Y = Q.SubMatrix(j1,M,j1,j2);
      BlockHouseholder_RDiv(Y,Z,m.Cols(j1,M));
      if (k==0) { TMVAssert(j1==0); }
      j2 = j1;
    }
#ifdef XDEBUG
    if (Norm(m-m2) > 0.001*Norm(m2)) {
      cerr<<"BlockQ_LDivEq: Q = "<<Type(Q)<<"  "<<Q<<endl;
      cerr<<"Qbeta = "<<Qbeta<<endl;
      cerr<<"m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"NonBlock m = "<<m2<<endl;
      abort(); 
    }
#endif
  }*/

  //
  // LDiv
  //

  template <class T1, class T2, class T3> void QR_LDiv(
      const GenMatrix<T1>& QRx, const GenVector<T1>& Qbeta, const size_t* P,
      const GenMatrix<T2>& m, const MatrixView<T3>& x, size_t N1)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(Qbeta.size() == QRx.rowsize());
    TMVAssert(m.colsize() == QRx.colsize());
    TMVAssert(x.colsize() == QRx.rowsize());
    TMVAssert(x.rowsize() == m.rowsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

    // Solve Q R P x = m
    // where Q and R are stored in QRx, and Qbeta are the beta

    // First Solve Q y = m
    if (QRx.IsSquare()) {
      x = m;
      Q_LDivEq(QRx,Qbeta,x);
    } else {
      if (m.isrm()) {
	Matrix<T3,RowMajor> m1 = m;
	// Q is Q1 [ I ]
	//         [ 0 ]
	// where Q1 is the part of Q that is stored in QRx and Qbeta
	// m1 = Q^-1 m1
	Q_LDivEq(QRx,Qbeta,m1.View());
	// y = [ I 0 ] m1
	x = m1.Rows(0,x.colsize()); // x = y here
      } else {
	Matrix<T3,ColMajor> m1 = m;
	Q_LDivEq(QRx,Qbeta,m1.View());
	x = m1.Rows(0,x.colsize()); // x = y here
      }
    }

    // Now solve R z = y
    x.Rows(N1,x.colsize()).Zero();
    x.Rows(0,N1) /= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);

    // Finally P x = z
    if (P) x.ReversePermuteRows(P);
  }

  /*
  template <class T1, class T2, class T3> void QR_LDiv(
      const GenMatrix<T1>& QRx, const vector<GenUpperTriMatrix<T1>*>& ZZ,
      const size_t* P, const GenMatrix<T2>& m, const MatrixView<T3>& x,
      size_t N1)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(m.colsize() == QRx.colsize());
    TMVAssert(x.colsize() == QRx.rowsize());
    TMVAssert(x.rowsize() == m.rowsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);

    if (QRx.IsSquare()) {
      x = m;
      Q_LDivEq(QRx,ZZ,x);
    } else {
      if (m.isrm()) {
	Matrix<T3,RowMajor> m1 = m;
	Q_LDivEq(QRx,ZZ,m1.View());
	x = m1.Rows(0,x.colsize());
      } else {
	Matrix<T3,ColMajor> m1 = m;
	Q_LDivEq(QRx,ZZ,m1.View());
	x = m1.Rows(0,x.colsize());
      }
    }

    x.Rows(N1,x.colsize()).Zero();
    x.Rows(0,N1) /= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);

    if (P) x.ReversePermuteRows(P);
  }*/

  //
  // LDivEq
  //

  template <class T1, class T2> void QR_LDivEq(
      const GenMatrix<T1>& QRx, const GenVector<T1>& Qbeta, const size_t* P,
      const MatrixView<T2>& m, size_t N1)
  {
    TMVAssert(QRx.colsize() == QRx.rowsize());
    TMVAssert(Qbeta.size() == QRx.rowsize());
    TMVAssert(m.colsize() == QRx.colsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

    // Solves Q R P x = m in place (m <- x)
    Q_LDivEq(QRx,Qbeta,m);
    m.Rows(N1,m.colsize()).Zero();
    m.Rows(0,N1) /= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);
    if (P) m.ReversePermuteRows(P);
  }

  /*
  template <class T1, class T2> void QR_LDivEq(
      const GenMatrix<T1>& QRx, const vector<GenUpperTriMatrix<T1>*>& ZZ,
      const size_t* P, const MatrixView<T2>& m, size_t N1)
  {
    TMVAssert(QRx.colsize() == QRx.rowsize());
    TMVAssert(m.colsize() == QRx.colsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);

    Q_LDivEq(QRx,ZZ,m);
    m.Rows(N1,m.colsize()).Zero();
    m.Rows(0,N1) /= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);
    if (P) m.ReversePermuteRows(P);
  }*/

  //
  // RDiv
  //

  template <class T1, class T2, class T3> void QR_RDiv(
      const GenMatrix<T1>& QRx, const GenVector<T1>& Qbeta, const size_t* P,
      const GenMatrix<T2>& m, const MatrixView<T3>& x, size_t N1)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(Qbeta.size() == QRx.rowsize());
    TMVAssert(x.rowsize() == QRx.colsize());
    TMVAssert(m.rowsize() == QRx.rowsize());
    TMVAssert(x.colsize() == m.colsize());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

    // Solve x Q R P = m
    // where Q and R are stored in QRx, and Qbeta are the beta

    // First solve y P = m
    x.Cols(0,m.rowsize()) = m;
    if (P) x.Cols(0,m.rowsize()).PermuteCols(P);

    // Next solve z R = y by forward substitution
    x.Cols(N1,x.rowsize()).Zero();
    x.Cols(0,N1) %= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);

    // Finally solve x Q = z
    // Q = Q1 [ I ]
    //        [ 0 ]
    // where Q1 is the part of Q that is stored in QRx and Qbeta
    // We've already dealt with the first part by zeroing out the 
    // right columns of x.
    Q_RDivEq(QRx,Qbeta,x);
  }

  /*
  template <class T1, class T2, class T3> void QR_RDiv(
      const GenMatrix<T1>& QRx, const vector<GenUpperTriMatrix<T1>*>& ZZ,
      const size_t* P, const GenMatrix<T2>& m, const MatrixView<T3>& x,
      size_t N1)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(x.rowsize() == QRx.colsize());
    TMVAssert(m.rowsize() == QRx.rowsize());
    TMVAssert(x.colsize() == m.colsize());
    TMVAssert(QRx.ct() == NonConj);

    x.Cols(0,m.rowsize()) = m;
    if (P) x.Cols(0,m.rowsize()).PermuteCols(P);

    x.Cols(N1,x.rowsize()).Zero();
    x.Cols(0,N1) %= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);

    Q_RDivEq(QRx,ZZ,x);
  }*/

  //
  // RDivEq
  //

  template <class T1, class T2> void QR_RDivEq(
      const GenMatrix<T1>& QRx, const GenVector<T1>& Qbeta, const size_t* P,
      const MatrixView<T2>& m, size_t N1)
  {
    TMVAssert(QRx.colsize() == QRx.rowsize());
    TMVAssert(Qbeta.size() == QRx.rowsize());
    TMVAssert(m.rowsize() == QRx.colsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

    // Solve x Q R P = m in place (m <- x)
 
    if (P) m.PermuteCols(P);
    m.Cols(N1,m.rowsize()).Zero();
    m.Cols(0,N1) %= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);
    Q_RDivEq(QRx,Qbeta,m);
  }

  /*
  template <class T1, class T2> void QR_RDivEq(
      const GenMatrix<T1>& QRx, const vector<GenUpperTriMatrix<T1>*>& ZZ,
      const size_t* P, const MatrixView<T2>& m, size_t N1)
  {
    TMVAssert(QRx.colsize() == QRx.rowsize());
    TMVAssert(m.rowsize() == QRx.colsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);

    if (P) m.PermuteCols(P);
    m.Cols(N1,m.rowsize()).Zero();
    m.Cols(0,N1) %= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);
    Q_RDivEq(QRx,ZZ,m);
  }*/

  // 
  // DownDate QR
  //
 
  template <class T> bool QR_DownDate(
      const MatrixView<T>& R, const GenVector<T>& zz)
  {
    // Given: A = P [ A1 ] = Q R  
    //              [ z  ]
    //
    // where P is a permutation, so z is any row of A.
    // 
    // Update R->R1 where A1 = Q1 R1
    //
    // The return value is whether the downdate was completed successfully.
    //
    // Note: P is an arbitrary permutation, which means you can remove
    // any row of A.
    // 
    // The input and output Q are ignored, so this doesn't fully
    // update a QR decomposition.  Just the R part.
    // This is useful when using QR for solving normal equations:
    // A b = c -> At A b = At c
    //
    // Initially solve using A = Q R -> b = (Rt R)^-1 At c
    // Don't need to keep Q for this.
    //
    // To remove a row i (eg for outlier rejection), use:
    // QR_DownDate(R,A.row(i)) and
    // At c -= A.row(i).Conjugate() * c(i)
    //
    // For how to do the downdate, see Golub and Van Loan, section 12.5.4
    // Here is the gist of it:
    //
    // A1t A1 = At A - zt z = Rt R - zt z 
    //        = [ Rt zt ] [ I 0  ] [ R ] 
    //                    [ 0 -1 ] [ z ]
    //        = [ Rt zt ] S [ R ]
    //                      [ z ]
    //
    // We find Hyperbolic rotation matrix H such that:
    // Ht S H = S and
    // H [ R ]  = [ R1 ]
    //   [ z ]    [ 0  ]
    //
    // Then, A1t A1 = R1t R1, so R1 is as desired.
    //
    // We zero out each element of z by multiplying the rows
    // i,n+1 by a matrix Hi where 
    // Hi = [ c  s  ] 
    //      [ s* c* ] 
    // with |c|^2 - |s|^2 = +-1, 
    // c Rii + s zi = Rii', and
    // s* Rii = -c* zi
    //
    // Each Hi zeros the ith element of z, and it preserves Ht S H = S.
    //
    // As with Givens, the phase of either c or s is unconstrained.
    // The solution with real c is:
    //
    // c = |Rii|/sqrt(|Rii|^2 - |zi|^2)
    // s = -(zi)*(Rii/|Rii|)/sqrt(|Rii|^2 - |zi|^2)
    //
    // or:
    //
    // c = 1/sqrt(1-|zi|^2/|Rii|^2)
    // s = -(zi/Rii)*/sqrt(1-|zi|^2/|Rii|^2)
    //
    // This produces a new Rii' = Rii sqrt(1-|zi|^2/|Rii|^2)
    //

    const size_t N = zz.size();
    TMVAssert(R.rowsize() == zz.size());
    TMVAssert(R.colsize() >= zz.size());

    Vector<T> z = zz;
    VectorView<T> Rdiag = R.diag();
    Vector<T> dz(N);
    for(size_t i=0; i < N; ++i) {
      if (z(i) == T(0)) continue;

      const T zoverr = z(i)/Rdiag(i);
      const RealType(T) normzoverr = tmv::NORM(zoverr);
      if (normzoverr >= 1) {
#ifdef TMVDEBUG
	cout<<"z = "<<z<<endl;
	cout<<"Rdiag = "<<Rdiag<<endl;
	cout<<"z.step = "<<z.step()<<endl;
	cout<<"zz.step = "<<z.step()<<endl;
	cout<<"Norm(zz-z) = "<<Norm(zz-z)<<endl;
	cout<<"|z("<<i<<")| = "<<abs(z(i))<<" >= |Rdiag("<<i<<")| = "<<abs(Rdiag(i))<<endl;
#endif
	return false;
      }

      const RealType(T) sqrtfactor = tmv::SQRT(1-normzoverr);
      // c = 1/sqrt(1-n)
      // c-1 = (1-sqrt(1-n))/sqrt(1-n)
      //     = n/sqrt(1-n)/(1+sqrt(1-n))
      //     = n/(1-n+sqrt(1-n))
      const RealType(T) cm1 = normzoverr/(1-normzoverr+sqrtfactor);
      const T s = -tmv::CONJ(zoverr)/sqrtfactor;
      // Rii' = Rii*sqrt(1-n) = Rii - Rii*(1-sqrt(1-n)) 
      //      = Rii - Rii*n/(1+sqrt(1-n))
      Rdiag(i) -= Rdiag(i)*normzoverr/(1+sqrtfactor);
      z(i) = T(0);

      dz.SubVector(i+1,N) = CONJ(s)*R.row(i,i+1,N) + cm1*z.SubVector(i+1,N);
      R.row(i,i+1,N) += cm1*R.row(i,i+1,N) + s*z.SubVector(i+1,N);
      z.SubVector(i+1,N) += dz.SubVector(i+1,N);
    }
    return true;
  }

  template <class T> bool QR_DownDate(const MatrixView<T>& Q,
      const MatrixView<T>& R, const GenVector<T>& zz)
  {
    // Given: A = [ A1 ] = Q R  
    //            [ z  ]
    //
    // Update R->R1 and Q->[ Q1 ], where A1 = Q1 R1
    //                     [ 0  ]
    //
    // The return value is whether the downdate was completed successfully.
    //
    // Note: Unlike the above version of DownDate, here z must be the 
    // last row of A, so appropriate permutations should be applied to 
    // Q if necessary before calling this function.
    // 
    // The method for doing this is very similar to above.
    //
    // Write Q = [ Q' ]
    //           [ q  ]
    //
    // [ A1 ] = [ Q' 0 ] [ I  0 ] [ R ] = [ Q' 0 ] S [ R ]
    // [ 0  ]   [ q  1 ] [ 0 -1 ] [ z ]   [ q  1 ]   [ z ]
    //
    // We find the same Hyperbolic rotation matrix H as above such that:
    // Ht S H = S and
    // H [ R ]  = [ R1 ]
    //   [ z ]    [ 0  ]
    //
    // Then we calculate [ Q' 0 ] Ht = [ Q1 ? ]  
    //                   [ q  1 ]      [ 0  ? ]
    // (The final column is garbage for our purposes.)
    //
    // Then, A1 = Q1 R1
    //

    TMVAssert(Q.rowsize() == zz.size());
    TMVAssert(Q.colsize() >= zz.size());
    TMVAssert(R.rowsize() == zz.size());
    TMVAssert(R.IsSquare());

    const size_t N = zz.size();

    Vector<T> z = zz;
    vector<RealType(T)> cm1(N);
    vector<T> s(N);

    VectorView<T> Rdiag = R.diag();
    Vector<T> dz(N);
    for(size_t i=0; i < N; ++i) {
      if (z(i) == T(0)) {
	s[i] = cm1[i] = RealType(T)(0);
	continue;
      }

      const T zoverr = z(i)/Rdiag(i);
      const RealType(T) normzoverr = tmv::NORM(zoverr);
      if (normzoverr >= 1) {
#ifdef TMVDEBUG
	cout<<"z = "<<z<<endl;
	cout<<"Rdiag = "<<Rdiag<<endl;
	cout<<"|z("<<i<<")| = "<<abs(z(i))<<" >= |Rdiag("<<i<<")| = "<<abs(Rdiag(i))<<endl;
#endif
	return false;
      }

      const RealType(T) sqrtfactor = tmv::SQRT(1-normzoverr);
      cm1[i] = normzoverr/(1-normzoverr+sqrtfactor);
      s[i] = -tmv::CONJ(zoverr)/sqrtfactor;
      Rdiag(i) -= Rdiag(i)*normzoverr/(1+sqrtfactor);
      z(i) = T(0);
      dz.SubVector(i+1,N) = 
	CONJ(s[i])*R.row(i,i+1,N) + cm1[i]*z.SubVector(i+1,N);
      R.row(i,i+1,N) += 
	cm1[i]*R.row(i,i+1,N) + s[i]*z.SubVector(i+1,N);
      z.SubVector(i+1,N) += dz.SubVector(i+1,N);
    }

    for(size_t k=0;k<Q.colsize()-1;++k) {
      // t = the element in the final column of the augmented matrix
      T t = T(0);
      for(size_t i=0; i<N; ++i) {
	T dqi = cm1[i]*Q(k,i) + tmv::CONJ(s[i])*t;
	t += s[i]*Q(k,i) + cm1[i]*t;
	Q(k,i) += dqi;
      }
    }

    return true;
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
    GetQFromQR(Q.View(),Qbeta);
    return Q;
  }

  template <class T> Matrix<T> QRPDiv<T>::GetQ() const
  {
    Matrix<T> Q = QRx;
    GetQFromQR(Q.View(),Qbeta);
    return Q;
  }

#define InstFile "TMV_QRDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


