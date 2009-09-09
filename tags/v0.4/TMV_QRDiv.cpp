
#include "TMV.h"
#include "TMV_Householder.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

namespace tmv {

  // 
  // Get Q from QR
  // 
 
  // MJ: Write level 3 version of this:
  template <class T> inline void NonLapGetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& Qbeta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    for(int i=N-1;i>=0;i--) {
      Q.row(i,i,N).Zero();
      Householder_Unpack(Q.SubMatrix(i,M,i,N),Qbeta(i));
    }
  }
#ifdef LAP
  template <class T> inline void LapGetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& Qbeta)
  { NonLapGetQFromQR(Q,Qbeta); }
  template <> inline void LapGetQFromQR(const MatrixView<double>& Q,
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
  template <> inline void LapGetQFromQR(
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
  template <> inline void LapGetQFromQR(const MatrixView<float>& Q,
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
  template <> inline void LapGetQFromQR(
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

  //
  // QR Decompose
  //

  // MJ: Write Level 3 version (with blocked Householder matrices)
  template <class T> inline void NonLapQR_Decompose(
      const MatrixView<T>& QRx,
      const VectorView<T>& Qbeta, T& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    // Decompose A (input as QRx) into A = Q R 
    // where Q is unitary, and R is upper triangular
    // Q and R are stored in the same matrix (QRx), with the beta's for 
    // the Householder matrices returned in Qbeta.
    const size_t M = QRx.colsize();
    const size_t N = QRx.rowsize();

    for(size_t j=0;j<N;++j) {
      // Apply the Householder Reflection for this column
      Qbeta(j) = Householder_Reflect(QRx.SubMatrix(j,M,j,N),det);
    }
  }
#ifdef LAP
  template <class T> inline void LapQR_Decompose(
      const MatrixView<T>& QRx,
      const VectorView<T>& Qbeta, T& det)
  { NonLapQR_Decompose(QRx,Qbeta,det); }
  template <> inline void LapQR_Decompose(const MatrixView<double>& QRx,
      const VectorView<double>& Qbeta, double& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = QRx.colsize();
    int n = QRx.rowsize();
    int lwork = 2*n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    if (QRx.isrm()) {
      int lda = QRx.stepi();
      dgelqf(&n,&m,QRx.ptr(),&lda,Qbeta.ptr(),work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dgelqf");
    } else {
      int lda = QRx.stepj();
      dgeqrf(&m,&n,QRx.ptr(),&lda,Qbeta.ptr(),work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"dgeqrf");
    }
    if (det) for(size_t i=0;i<Qbeta.size();++i) if (Qbeta(i) != 0.) det = -det;
  }
  template <> inline void LapQR_Decompose(
      const MatrixView<complex<double> >& QRx,
      const VectorView<complex<double> >& Qbeta, complex<double>& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = QRx.colsize();
    int n = QRx.rowsize();
    int lwork = 2*n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    if (QRx.isrm()) {
      int lda = QRx.stepi();
      zgelqf(&n,&m,LAP_Complex(QRx.ptr()),&lda,LAP_Complex(Qbeta.ptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"zgelqf");
    } else {
      int lda = QRx.stepj();
      zgeqrf(&m,&n,LAP_Complex(QRx.ptr()),&lda,LAP_Complex(Qbeta.ptr()),
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
  template <> inline void LapQR_Decompose(const MatrixView<float>& QRx,
      const VectorView<float>& Qbeta, float& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = QRx.colsize();
    int n = QRx.rowsize();
    int lwork = 2*n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    if (QRx.isrm()) {
      int lda = QRx.stepi();
      sgelqf(&n,&m,QRx.ptr(),&lda,Qbeta.ptr(),work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sgelqf");
    } else {
      int lda = QRx.stepj();
      sgeqrf(&m,&n,QRx.ptr(),&lda,Qbeta.ptr(),work,&lwork,&info);
      LAP_Results(info,int(work[0]),m,n,lwork,"sgeqrf");
    }
    if (det) for(size_t i=0;i<Qbeta.size();++i) if (Qbeta(i) != 0.) det = -det;
  }
  template <> inline void LapQR_Decompose(
      const MatrixView<complex<float> >& QRx,
      const VectorView<complex<float> >& Qbeta, complex<float>& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = QRx.colsize();
    int n = QRx.rowsize();
    int lwork = 2*n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    if (QRx.isrm()) {
      int lda = QRx.stepi();
      cgelqf(&n,&m,LAP_Complex(QRx.ptr()),&lda,LAP_Complex(Qbeta.ptr()),
	  LAP_Complex(work),&lwork,&info);
      LAP_Results(info,int(real(work[0])),m,n,lwork,"cgelqf");
    } else {
      int lda = QRx.stepj();
      cgeqrf(&m,&n,LAP_Complex(QRx.ptr()),&lda,LAP_Complex(Qbeta.ptr()),
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
      const MatrixView<T>& QRx,
      const VectorView<T>& Qbeta, T& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
#ifdef LAP
    if (Qbeta.step()==1 && (QRx.isrm() || QRx.iscm()))
      LapQR_Decompose(QRx,Qbeta,det);
    else 
#endif
      NonLapQR_Decompose(QRx,Qbeta,det);
  }

  //
  // QRP Decompose
  //

  // MJ: Try to figure out level 3 version, supposedly with look-aheads
  // to determine which pivots will be necessary somehow.
  template <class T> void NonLapQRP_Decompose(
      const MatrixView<T>& QRx,
      const VectorView<T>& Qbeta, Permutation& P, T& det)
  {
    // Decompose A (input as QRx) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (QRx), with the beta's for
    // the Householder matrices returned in Qbeta.
    const size_t M = QRx.colsize();
    const size_t N = QRx.rowsize();
    TMVAssert(M >= N);
    TMVAssert(Qbeta.size() == N);
    TMVAssert(P.size() == N);
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);

    P.SetToIdentity();

    // Keep track of the norm of each column
    // When considering column j, these are actually just the norm
    // of each column from j:M, not 0:M.
    Vector<RealType(T)> colnormsq(N);
    for(size_t j=0;j<N;++j) colnormsq(j) = NormSq(QRx.col(j));
    RealType(T) anormsq = colnormsq.SumElements();
    RealType(T) thresh2 = tmv::Epsilon<T>() * anormsq;
    RealType(T) thresh = tmv::Epsilon<T>() * thresh2;
    // Set to 0 any diag element whose norm is < epsilon * |A|
    // thresh2 is the threshold for recalculating the norm to account for
    // rounding errors in the subtractions which keep track of it..

    for(size_t j=0;j<N;++j) {
      // Find the column with the largest norm
      size_t k;
      RealType(T) maxnormsq = MaxElement(colnormsq.SubVector(j,N),&k);
      // Note: k is relative to the SubVector(j,N), so add j to get real index
      k += j;
      maxnormsq = NormSq(QRx.col(k,j,M));

      // If largest norm is 0 (technically < thresh to account for rounding)
      // then the rest of the R matrix is 0, and the Householder matrices 
      // are identities (indicated by 0's in the Q part of the matrix).
      if (maxnormsq < thresh) {
	QRx.SubMatrix(j,M,j,N).Zero();
	// Already essentially zero - make it exact
	Qbeta.SubVector(j,N).Zero();
	// Set the Householder matrices for these to identities
	break;
      }

      // Swap the column with the largest norm into the current column
      if (k != j) {
	colnormsq.Swap(j,k);
	P.SwapRows(j,k);
	QRx.SwapCols(j,k);
      }
      // Apply the Householder Reflection for this column
      Qbeta(j) = Householder_Reflect(QRx.SubMatrix(j,M,j,N),det);

      // And update the norms for use with the next column
      for(size_t k=j+1;k<N;++k) {
	colnormsq(k) -= tmv::NORM(QRx(j,k));
      }
    }
    if (det!=T(0)) det *= P.Det();
  }
#ifdef XLAP
  template <class T> inline void LapQRP_Decompose(
      const MatrixView<T>& QRx,
      const VectorView<T>& Qbeta, Permutation& P, T& det)
  { NonLapQR_Decompose(QRx,Qbeta,det); }
  template <> inline void LapQRP_Decompose(const MatrixView<double>& QRx,
      const VectorView<double>& Qbeta, Permutation& P, double& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = QRx.colsize();
    int n = QRx.rowsize();
    int lap_p[P.size()];
    int lwork = 3*n;
    double* work = LAP_DWork(lwork);
    int info;
    int lda = QRx.stepj();
    char cc = 'F';
    double thresh = Epsilon<double>()*dlange(&cc,&m,&n,QRx.ptr(),&lda,0);
    dgeqpf(&m,&n,QRx.ptr(),&lda,lap_p,Qbeta.ptr(),work,&info);
    LAP_Results(info,"dgeqpf");
    for(size_t i=0;i<QRx.rowsize();++i) {
      if (abs(QRx(i,i)) < thresh) {
	TMVAssert(NormInf(QRx.diag(0,i+1,QRx.rowsize())) < thresh);
	QRx.SubMatrix(i,QRx.colsize(),i,QRx.rowsize()).Zero();
	Qbeta.SubVector(i,Qbeta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<Qbeta.size();++i) {
      if (det) if (Qbeta(i) != 0.) det = -det;
      lap_p[i]--;
    }
    P = Permutation(P.size(),lap_p,false);
    if (det) det *= P.Det();
  }
  template <> inline void LapQRP_Decompose(
      const MatrixView<complex<double> >& QRx,
      const VectorView<complex<double> >& Qbeta, Permutation& P,
      complex<double>& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = QRx.colsize();
    int n = QRx.rowsize();
    int lap_p[P.size()];
    int lwork = 2*n;
    double* rwork = LAP_DWork(lwork);
    lwork = 3*n;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    int lda = QRx.stepj();
    char cc = 'F';
    double thresh = Epsilon<double>()*
      zlange(&cc,&m,&n,LAP_Complex(QRx.ptr()),&lda,0);
    zgeqpf(&m,&n,LAP_Complex(QRx.ptr()),&lda,lap_p,LAP_Complex(Qbeta.ptr()),
	LAP_Complex(work),rwork,&info);
    LAP_Results(info,"zgeqpf");
    Qbeta.ConjugateSelf();
    for(size_t i=0;i<QRx.rowsize();++i) {
      TMVAssert(abs(imag(QRx(i,i))) < Epsilon<double>());
      if (abs(real(QRx(i,i))) < thresh) {
	TMVAssert(NormInf(QRx.diag(0,i+1,QRx.rowsize())) < thresh);
	QRx.SubMatrix(i,QRx.colsize(),i,QRx.rowsize()).Zero();
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
      lap_p[i]--;
    }
    P = Permutation(P.size(),lap_p,false);
    if (det!=double(0)) det *= P.Det();
  }
#ifndef NOFLOAT
  template <> inline void LapQRP_Decompose(const MatrixView<float>& QRx,
      const VectorView<float>& Qbeta, Permutation& P, float& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = QRx.colsize();
    int n = QRx.rowsize();
    int lap_p[P.size()];
    int lwork = 3*n;
    float* work = LAP_SWork(lwork);
    int info;
    int lda = QRx.stepj();
    char cc = 'F';
    float thresh = Epsilon<float>()*slange(&cc,&m,&n,QRx.ptr(),&lda,0);
    sgeqpf(&m,&n,QRx.ptr(),&lda,lap_p,Qbeta.ptr(),work,&info);
    LAP_Results(info,"sgeqpf");
    for(size_t i=0;i<QRx.rowsize();++i) {
      if (abs(QRx(i,i)) < thresh) {
	TMVAssert(NormInf(QRx.diag(0,i+1,QRx.rowsize())) < thresh);
	QRx.SubMatrix(i,QRx.colsize(),i,QRx.rowsize()).Zero();
	Qbeta.SubVector(i,Qbeta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<Qbeta.size();++i) {
      if (det) if (Qbeta(i) != 0.) det = -det;
      lap_p[i]--;
    }
    P = Permutation(P.size(),lap_p,false);
    if (det) det *= P.Det();
  }
  template <> inline void LapQRP_Decompose(
      const MatrixView<complex<float> >& QRx,
      const VectorView<complex<float> >& Qbeta, Permutation& P,
      complex<float>& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = QRx.colsize();
    int n = QRx.rowsize();
    int lap_p[P.size()];
    int lwork = 2*n;
    float* rwork = LAP_SWork(lwork);
    lwork = 3*n;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    int lda = QRx.stepj();
    char cc = 'F';
    float thresh = Epsilon<float>()*
      clange(&cc,&m,&n,LAP_Complex(QRx.ptr()),&lda,0);
    cgeqpf(&m,&n,LAP_Complex(QRx.ptr()),&lda,lap_p,LAP_Complex(Qbeta.ptr()),
	LAP_Complex(work),rwork,&info);
    LAP_Results(info,"cgeqpf");
    for(size_t i=0;i<QRx.rowsize();++i) {
      TMVAssert(abs(imag(QRx(i,i))) < Epsilon<double>());
      if (abs(real(QRx(i,i))) < thresh) {
	TMVAssert(NormInf(QRx.diag(0,i+1,QRx.rowsize())) < thresh);
	QRx.SubMatrix(i,QRx.colsize(),i,QRx.rowsize()).Zero();
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
      lap_p[i]--;
    }
    P = Permutation(P.size(),lap_p,false);
    if (det!=float(0)) det *= P.Det();
  }
#endif
  template <class T> inline void NewLapQRP_Decompose(
      const MatrixView<T>& QRx,
      const VectorView<T>& Qbeta, Permutation& P, T& det)
  { NonLapQR_Decompose(QRx,Qbeta,det); }
  template <> inline void NewLapQRP_Decompose(const MatrixView<double>& QRx,
      const VectorView<double>& Qbeta, Permutation& P, double& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = QRx.colsize();
    int n = QRx.rowsize();
    int lap_p[P.size()];
    int lwork = 3*n+1;
    double* work = LAP_DWork(lwork);
    int info;
    int lda = QRx.stepj();
    char cc = 'F';
    double thresh = Epsilon<double>()*dlange(&cc,&m,&n,QRx.ptr(),&lda,0);
    dgeqp3(&m,&n,QRx.ptr(),&lda,lap_p,Qbeta.ptr(),work,&lwork,&info);
    LAP_Results(info,"dgeqp3");
    for(size_t i=0;i<QRx.rowsize();++i) {
      if (abs(QRx(i,i)) < thresh) {
	TMVAssert(NormInf(QRx.diag(0,i+1,QRx.rowsize())) < thresh);
	QRx.SubMatrix(i,QRx.colsize(),i,QRx.rowsize()).Zero();
	Qbeta.SubVector(i,Qbeta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<Qbeta.size();++i) {
      if (det) if (Qbeta(i) != 0.) det = -det;
      lap_p[i]--;
    }
    P = Permutation(P.size(),lap_p,false);
    if (det) det *= P.Det();
  }
  template <> inline void NewLapQRP_Decompose(
      const MatrixView<complex<double> >& QRx,
      const VectorView<complex<double> >& Qbeta, Permutation& P,
      complex<double>& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = QRx.colsize();
    int n = QRx.rowsize();
    int lap_p[P.size()];
    int lwork = 2*n;
    double* rwork = LAP_DWork(lwork);
    lwork = n+1;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    int lda = QRx.stepj();
    char cc = 'F';
    double thresh = Epsilon<double>()*
      zlange(&cc,&m,&n,LAP_Complex(QRx.ptr()),&lda,0);
    zgeqp3(&m,&n,LAP_Complex(QRx.ptr()),&lda,lap_p,LAP_Complex(Qbeta.ptr()),
	LAP_Complex(work),&lwork,rwork,&info);
    LAP_Results(info,"zgeqp3");
    Qbeta.ConjugateSelf();
    for(size_t i=0;i<QRx.rowsize();++i) {
      TMVAssert(abs(imag(QRx(i,i))) < Epsilon<double>());
      if (abs(real(QRx(i,i))) < thresh) {
	TMVAssert(NormInf(QRx.diag(0,i+1,QRx.rowsize())) < thresh);
	QRx.SubMatrix(i,QRx.colsize(),i,QRx.rowsize()).Zero();
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
      lap_p[i]--;
    }
    P = Permutation(P.size(),lap_p,false);
    if (det!=double(0)) det *= P.Det();
  }
#ifndef NOFLOAT
  template <> inline void NewLapQRP_Decompose(const MatrixView<float>& QRx,
      const VectorView<float>& Qbeta, Permutation& P, float& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = QRx.colsize();
    int n = QRx.rowsize();
    int lap_p[P.size()];
    int lwork = 3*n+1;
    float* work = LAP_SWork(lwork);
    int info;
    int lda = QRx.stepj();
    char cc = 'F';
    float thresh = Epsilon<float>()*slange(&cc,&m,&n,QRx.ptr(),&lda,0);
    sgeqp3(&m,&n,QRx.ptr(),&lda,lap_p,Qbeta.ptr(),work,&lwork,&info);
    LAP_Results(info,"sgeqp3");
    for(size_t i=0;i<QRx.rowsize();++i) {
      if (abs(QRx(i,i)) < thresh) {
	TMVAssert(NormInf(QRx.diag(0,i+1,QRx.rowsize())) < thresh);
	QRx.SubMatrix(i,QRx.colsize(),i,QRx.rowsize()).Zero();
	Qbeta.SubVector(i,Qbeta.size()).Zero();
	break;
      }
    }
    for(size_t i=0;i<Qbeta.size();++i) {
      if (det) if (Qbeta(i) != 0.) det = -det;
      lap_p[i]--;
    }
    P = Permutation(P.size(),lap_p,false);
    if (det) det *= P.Det();
  }
  template <> inline void NewLapQRP_Decompose(
      const MatrixView<complex<float> >& QRx,
      const VectorView<complex<float> >& Qbeta, Permutation& P,
      complex<float>& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    int m = QRx.colsize();
    int n = QRx.rowsize();
    int lap_p[P.size()];
    int lwork = 2*n;
    float* rwork = LAP_SWork(lwork);
    lwork = n+1;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    int lda = QRx.stepj();
    char cc = 'F';
    float thresh = Epsilon<float>()*
      clange(&cc,&m,&n,LAP_Complex(QRx.ptr()),&lda,0);
    cgeqp3(&m,&n,LAP_Complex(QRx.ptr()),&lda,lap_p,LAP_Complex(Qbeta.ptr()),
	LAP_Complex(work),&lwork,rwork,&info);
    LAP_Results(info,"cgeqp3");
    for(size_t i=0;i<QRx.rowsize();++i) {
      TMVAssert(abs(imag(QRx(i,i))) < Epsilon<double>());
      if (abs(real(QRx(i,i))) < thresh) {
	TMVAssert(NormInf(QRx.diag(0,i+1,QRx.rowsize())) < thresh);
	QRx.SubMatrix(i,QRx.colsize(),i,QRx.rowsize()).Zero();
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
      lap_p[i]--;
    }
    P = Permutation(P.size(),lap_p,false);
    if (det!=float(0)) det *= P.Det();
  }
#endif
#endif
  template <class T> void QRP_Decompose(
      const MatrixView<T>& QRx,
      const VectorView<T>& Qbeta, Permutation& P, T& det)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(QRx.rowsize() == Qbeta.size());
    TMVAssert(P.size() == QRx.rowsize());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    // MJ: Neither LAP version of QRP seems to do any actual pivoting.
    // At least for the Intel MKL version.
    // So until I figure out why and how to get them to work, leave
    // this is XLAP (ie. ignore the LAP routines.
#ifdef XLAP
    if (QRx.iscm())
      LapQRP_Decompose(QRx,Qbeta,P,det);
    else {
#ifdef TMVDEBUG
      cout<<"Lap QRDecomp: QRx, Qbeta are wrong step:\n";
      cout<<"QRx isrm = "<<QRx.isrm()<<", stepi = "<<QRx.stepi()<<endl;
      cout<<"Qbeta step = "<<Qbeta.step()<<endl;
#endif
      NonLapQRP_Decompose(QRx,Qbeta,P,det);
    }
#else
    NonLapQRP_Decompose(QRx,Qbeta,P,det);
#endif
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
      const MatrixView<T>& Q, const MatrixView<T>& R, Permutation& P, T& det)
  {
    // Decompose A (input as Q) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation

    const size_t N = Q.rowsize();
    TMVAssert(Q.colsize() >= N);
    TMVAssert(R.colsize() == N);
    TMVAssert(R.rowsize() == N);
    TMVAssert(P.size() == N);
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

  template <class T1, class T2> inline void NonLapQ_LDivEq(
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
      Householder_Mult(Q.col(j,j+1,M),Qbeta(j),m.Rows(j,M));
    }
  }
#ifdef LAP
  template <class T1, class T2> inline void LapQ_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  { NonLapQ_LDivEq(Q,Qbeta,m); }
  template <> inline void LapQ_LDivEq(const GenMatrix<double>& Q,
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
  template <> inline void LapQ_LDivEq(
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
  template <> inline void LapQ_LDivEq(const GenMatrix<float>& Q,
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
  template <> inline void LapQ_LDivEq(
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

  //
  // Packed Q - RDivEq
  //

  template <class T1, class T2> inline void NonLapQ_RDivEq(
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
    // xT = Q* mT (t = Adjoint, T = Transpose, * = Conjugate, and Q Qt = I)
    // Q* is H0T H1T ... H_N-1T
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    for(int j=N-1;j>=0;j--) if (Qbeta(j) != T1(0)) {
      Householder_ConjMult(Q.col(j,j+1,M),CONJ(Qbeta(j)),
	  m.QuickTranspose().Rows(j,M));
    }
  }
#ifdef LAP
  template <class T1, class T2> inline void LapQ_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  { NonLapQ_RDivEq(Q,Qbeta,m); }
  template <> inline void LapQ_RDivEq(const GenMatrix<double>& Q,
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
  template <> inline void LapQ_RDivEq(const GenMatrix<complex<double> >& Q,
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
  template <> inline void LapQ_RDivEq(const GenMatrix<float>& Q,
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
  template <> inline void LapQ_RDivEq(const GenMatrix<complex<float> >& Q,
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

  //
  // LDiv
  //

  template <class T1, class T2, class T3> void QR_LDiv(
      const GenMatrix<T1>& QRx, const GenVector<T1>& Qbeta,
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
    if (N1 == 0) N1 = Qbeta.size();

    // Solve Q R x = m
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
	// MJ: Look into writing version of Q_LDivEq that only does N rows
	Q_LDivEq(QRx,Qbeta,m1.QuickView());
	// y = [ I 0 ] m1
	x = m1.Rows(0,x.colsize()); // x = y here
      } else {
	Matrix<T3,ColMajor> m1 = m;
	Q_LDivEq(QRx,Qbeta,m1.QuickView());
	x = m1.Rows(0,x.colsize()); // x = y here
      }
    }

    // Now solve R x = y
    x.Rows(N1,x.colsize()).Zero();
    x.Rows(0,N1) /= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);
  }

  //
  // LDivEq
  //

  template <class T1, class T2> void QR_LDivEq(
      const GenMatrix<T1>& QRx, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m, size_t N1)
  {
    TMVAssert(QRx.colsize() == QRx.rowsize());
    TMVAssert(Qbeta.size() == QRx.rowsize());
    TMVAssert(m.colsize() == QRx.colsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    if (N1 == 0) N1 = Qbeta.size();

    // Solves Q R x = m in place (m <- x)
    Q_LDivEq(QRx,Qbeta,m);
    m.Rows(N1,m.colsize()).Zero();
    m.Rows(0,N1) /= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);
  }

  //
  // RDiv
  //

  template <class T1, class T2, class T3> void QR_RDiv(
      const GenMatrix<T1>& QRx, const GenVector<T1>& Qbeta, 
      const GenMatrix<T2>& m, const MatrixView<T3>& x, size_t N1)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(Qbeta.size() == QRx.rowsize());
    TMVAssert(x.rowsize() == QRx.colsize());
    TMVAssert(m.rowsize() == QRx.rowsize());
    TMVAssert(x.colsize() == m.colsize());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    if (N1 == 0) N1 = Qbeta.size();

    // Solve x Q R = m
    // where Q and R are stored in QRx, and Qbeta are the beta

    // First solve y R = m by forward substitution
    x.Cols(N1,x.rowsize()).Zero();
    x.Cols(0,N1) = m.Cols(0,N1) % UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);

    // Now solve x Q = y
    // Q = Q1 [ I ]
    //        [ 0 ]
    // where Q1 is the part of Q that is stored in QRx and Qbeta
    // We've already dealt with the first part by zeroing out the 
    // bottom of x and using only the top of x above.
    Q_RDivEq(QRx,Qbeta,x);
  }

  //
  // RDivEq
  //

  template <class T1, class T2> void QR_RDivEq(
      const GenMatrix<T1>& QRx, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m, size_t N1)
  {
    TMVAssert(QRx.colsize() == QRx.rowsize());
    TMVAssert(Qbeta.size() == QRx.rowsize());
    TMVAssert(m.rowsize() == QRx.colsize());
    TMVAssert(QRx.isrm() || QRx.iscm());
    TMVAssert(QRx.ct() == NonConj);
    TMVAssert(Qbeta.ct() == NonConj);
    if (N1 == 0) N1 = Qbeta.size();

    // Solve x Q R = m in place (m <- x)
 
    m.Cols(N1,m.rowsize()).Zero();
    m.Cols(0,N1) %= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);
    Q_RDivEq(QRx,Qbeta,m);
  }

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
    // with |c|^2 - |s|^2 = 1, 
    // c Rii + s zi = Rii', and
    // s* Rii + c* zi = 0
    //
    // Each Hi zeros the ith element of z, and it preserves Ht S H = S.
    //
    // As with Givens, the phase of either c or s is unconstrained.
    // The solution with real c is:
    //
    // c = +- 1/sqrt(1 - |zi/Rii|^2)
    // s = -c (zi/Rii)*
    //
    // or:
    //
    // This produces a new Rii' = +- Rii sqrt(1-|zi|^2/|Rii|^2)
    //
    // I don't see any reason to use the - sign, so I always choose + here.
    //

    const size_t N = zz.size();
    TMVAssert(R.rowsize() == zz.size());
    TMVAssert(R.colsize() >= zz.size());
    RealType(T) RT1(1);

    Vector<T> z = zz;
    VectorView<T> Rdiag = R.diag();
    Vector<T> newz(N);
    for(size_t i=0; i < N; ++i) {
      if (z(i) == T(0)) continue;

      const T zoverr = z(i)/Rdiag(i);
      const RealType(T) normzoverr = tmv::NORM(zoverr);
      if (normzoverr >= RT1) {
#ifdef TMVDEBUG
	cout<<"(orig) z = "<<zz<<endl;
	cout<<"Rdiag = "<<Rdiag<<endl;
	cout<<"|z("<<i<<")| = "<<abs(z(i))<<" >= |Rdiag("<<i<<")| = "<<abs(Rdiag(i))<<endl;
#endif
	return false;
      }

      const RealType(T) sqrtfactor = tmv::SQRT(RT1-normzoverr);
      // c = 1/sqrt(1-n)
      //const RealType(T) c = RT1/sqrtfactor;
      const RealType(T) c = -RT1/sqrtfactor;
      const T s = -CONJ(zoverr)*c;
      //Rdiag(i) *= sqrtfactor;
      Rdiag(i) *= -sqrtfactor;
      if (REAL(zz(0)) > 19.559 && REAL(zz(0)) < 19.560 && N == 2310 && i <= 2012) {
	cout<<"z(2012) = "<<z(2012)<<"  s,c = "<<s<<", "<<c;
      }
      z(i) = T(0);

      newz.SubVector(i+1,N) = CONJ(s)*R.row(i,i+1,N) + c*z.SubVector(i+1,N);
      R.row(i,i+1,N) *= c;
      R.row(i,i+1,N) += s*z.SubVector(i+1,N);
      z.SubVector(i+1,N) = newz.SubVector(i+1,N);
      if (REAL(zz(0)) > 19.559 && REAL(zz(0)) < 19.560 && N == 2310 && i <= 2012) {
	cout<<"  --> z(2012) = "<<z(2012)<<endl;
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

  template <class T> Matrix<T,ColMajor> QRDiv<T>::Inverse() const
  {
    if (istrans) {
      Matrix<T,ColMajor> inv(QRx.colsize(),QRx.rowsize());
      LDiv(Eye<RealType(T),ColMajor>(QRx.rowsize()),inv.QuickView());
      return inv;
    } else {
      Matrix<T,ColMajor> inv(QRx.rowsize(),QRx.colsize());
      RDiv(Eye<RealType(T),ColMajor>(QRx.rowsize()),inv.QuickView());
      return inv;
    }
  }

  template <class T> Matrix<T,ColMajor> QRDiv<T>::DoInverseATA() const
  {
    // At A = Rt Qt Q R = Rt R
    // (At A)^-1 = (Rt R)^-1 = R^-1 * Rt^-1
    const size_t N = istrans ? QRx.colsize() : QRx.rowsize();

    UpperTriMatrix<T,NonUnitDiag,ColMajor> rinv(N);
    rinv.SetToIdentity();
    if (istrans) {
      rinv.SubTriMatrix(0,N1) /= 
	UpperTriMatrixViewOf(QRx.QuickTranspose()).SubTriMatrix(0,N1);
    } else {
      rinv.SubTriMatrix(0,N1) /= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);
    }
    if (isshort)
      return Matrix<T,ColMajor>(rinv.QuickConjugate() * rinv.QuickTranspose());
    else
      return Matrix<T,ColMajor>(rinv * rinv.QuickAdjoint());
  }

  template <class T> Matrix<T,ColMajor> QRDiv<T>::QR_GetQ() const
  {
    Matrix<T,ColMajor> Q = (istrans != isshort ? QRx.QuickTranspose() :
	QRx.QuickView());
    GetQFromQR(Q.QuickView(),Qbeta);
    return Q;
  }

  template <class T> Matrix<T,ColMajor> QRDiv<T>::QR_GetR() const
  { 
    if (istrans != isshort)
      return Matrix<T,ColMajor>(UpperTriMatrixViewOf(QRx.QuickTranspose())); 
    else
      return Matrix<T,ColMajor>(UpperTriMatrixViewOf(QRx)); 
  }

  template <class T> T QRPDiv<T>::Det() const
  {
    if (!donedet) {
      det *= DiagMatrixViewOf(QRx.diag()).Det();
      donedet = true;
    }
    return det;
  }

  template <class T> Matrix<T,ColMajor> QRPDiv<T>::Inverse() const
  {
    if (istrans) {
      Matrix<T,ColMajor> inv(QRx.colsize(),QRx.rowsize());
      LDiv(Eye<RealType(T),ColMajor>(QRx.rowsize()),inv.QuickView());
      return inv;
    } else {
      Matrix<T,ColMajor> inv(QRx.rowsize(),QRx.colsize());
      RDiv(Eye<RealType(T),ColMajor>(QRx.rowsize()),inv.QuickView());
      return inv;
    }
  }

  template <class T> Matrix<T,ColMajor> QRPDiv<T>::DoInverseATA() const
  {
    // At A = Rt R
    // (At A)^-1 = (Rt R)^-1 = R^-1 * Rt^-1
    // AT A = PT RT R P
    // (At A)^-1 = (Pt Rt R P)^-1 = P^-1 R^-1 * (P^-1 R^-1)t
    const size_t N = QRx.rowsize();

    Matrix<T,ColMajor> temp(N,N);
    temp.SetToIdentity();
    temp.SubMatrix(0,N1,0,N1) /= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);
    temp /= P;
    if (istrans) 
      return temp.QuickConjugate() * temp.QuickTranspose();
    else
      return temp * temp.QuickAdjoint();
  }

  template <class T> Matrix<T,ColMajor> QRPDiv<T>::QR_GetQ() const
  {
    Matrix<T,ColMajor> Q = QRx;
    GetQFromQR(Q.QuickView(),Qbeta);
    return Q;
  }

  template <class T> Matrix<T,ColMajor> QRPDiv<T>::QR_GetR() const
  { return Matrix<T,ColMajor>(UpperTriMatrixViewOf(QRx)); }

#define InstFile "TMV_QRDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


