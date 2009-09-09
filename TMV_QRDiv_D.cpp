
#include "TMV.h"
#include "TMV_Householder.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

#define RecursiveQR

#ifdef TMV_BLOCKSIZE
  const size_t QR_BLOCKSIZE = TMV_BLOCKSIZE;
#else
  const size_t QR_BLOCKSIZE = 64;
#endif

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

#define InstFile "TMV_QRDiv_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


