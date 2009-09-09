
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
  // Get Q from QR
  // 
 
  template <class T> inline void NonBlockGetQFromQR(
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

  template <class T> inline void BlockGetQFromQR(
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
    UpperTriMatrix<T,NonUnitDiag,ColMajor> BaseZ(min(size_t(QR_BLOCKSIZE),N));
    for(size_t j2=N;j2>0;) {
      size_t j1 = j2 > QR_BLOCKSIZE ? j2-QR_BLOCKSIZE : 0;
      MatrixView<T> Y = Q.SubMatrix(j1,M,j1,j2);
      UpperTriMatrixView<T> Z = BaseZ.SubTriMatrix(0,Y.rowsize());
      BlockHouseholder_MakeZ(Y,Z,beta.SubVector(j1,j2));
      BlockHouseholder_Unpack(Y,Z,Q.SubMatrix(j1,M,j2,N));
      j2 = j1;
    }
#ifdef XDEBUG
    if (Norm(Q-Q2) > 0.001*Norm(Q)*Norm(beta)) {
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
  template <> inline void LapGetQFromQR(const MatrixView<double>& Q,
      const GenVector<double>& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
#endif
    if (Q.isrm()) {
      int ldq = Q.stepi();
      LAPNAME(dorglq) (LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
	  LAPP(beta.cptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("dorglq");
#else
      LAP_Results(int(work[0]),m,n,lwork,"dorglq");
#endif
    } else {
      int ldq = Q.stepj();
      LAPNAME(dorgqr) (LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
	  LAPP(beta.cptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("dorgqr");
#else
      LAP_Results(int(work[0]),m,n,lwork,"dorgqr");
#endif
    }
  }
  template <> inline void LapGetQFromQR(
      const MatrixView<complex<double> >& Q,
      const GenVector<complex<double> >& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
#endif
    if (Q.isrm()) {
      int ldq = Q.stepi();
      LAPNAME(zunglq) (LAPCM LAPV(n),LAPV(m),LAPV(n),
	  LAPP(Q.ptr()),LAPV(ldq),LAPP(beta.cptr())
	  LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("zunglq");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"zunglq");
#endif
    } else {
      int ldq = Q.stepj();
      Vector<complex<double> > conjbeta = beta.Conjugate();
      LAPNAME(zungqr) (LAPCM LAPV(m),LAPV(n),LAPV(n),
	  LAPP(Q.ptr()),LAPV(ldq),LAPP(conjbeta.cptr())
	  LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("zungqr");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"zungqr");
#endif
    }
  }
#ifndef NOFLOAT
  template <> inline void LapGetQFromQR(const MatrixView<float>& Q,
      const GenVector<float>& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
#endif
    if (Q.isrm()) {
      int ldq = Q.stepi();
      LAPNAME(sorglq) (LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
	  LAPP(beta.cptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("sorglq");
#else
      LAP_Results(int(work[0]),m,n,lwork,"sorglq");
#endif
    } else {
      int ldq = Q.stepj();
      LAPNAME(sorgqr) (LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
	  LAPP(beta.cptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("sorgqr");
#else
      LAP_Results(int(work[0]),m,n,lwork,"sorgqr");
#endif
    }
  }
  template <> inline void LapGetQFromQR(
      const MatrixView<complex<float> >& Q,
      const GenVector<complex<float> >& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
#endif
    if (Q.isrm()) {
      int ldq = Q.stepi();
      LAPNAME(cunglq) (LAPCM LAPV(n),LAPV(m),LAPV(n),
	  LAPP(Q.ptr()),LAPV(ldq),LAPP(beta.cptr())
	  LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("cunglq");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"cunglq");
#endif
    } else {
      int ldq = Q.stepj();
      Vector<complex<float> > conjbeta = beta.Conjugate();
      LAPNAME(cungqr) (LAPCM LAPV(m),LAPV(n),LAPV(n),
	  LAPP(Q.ptr()),LAPV(ldq),LAPP(conjbeta.cptr())
	  LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("cungqr");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"cungqr");
#endif
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
#ifdef LAP
    if (beta.isconj()) {
      Vector<T> b2 = beta;
      LapGetQFromQR(Q,b2.View());
    } else {
      LapGetQFromQR(Q,beta);
    }
#else
    NonLapGetQFromQR(Q,beta);
#endif
  }

#define InstFile "TMV_QRDiv_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


