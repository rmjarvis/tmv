
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
  // Packed Q - LDivEq
  //

  template <class T, class T1> inline void NonBlockQ_LDivEq(
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

  template <class T, class T1> inline void BlockQ_LDivEq(
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
    UpperTriMatrix<T1,NonUnitDiag,ColMajor> BaseZ(min(size_t(QR_BLOCKSIZE),N));
    for(size_t j1=0;j1<N;) {
      size_t j2 = min(N,j1+QR_BLOCKSIZE);
      ConstMatrixView<T1> Y = Q.SubMatrix(j1,M,j1,j2);
      UpperTriMatrixView<T1> Z = BaseZ.SubTriMatrix(0,Y.rowsize());
      BlockHouseholder_MakeZ(Y,Z,beta.SubVector(j1,j2));
      BlockHouseholder_LDiv(Y,Z,m.Rows(j1,M));
      j1 = j2;
    }
#ifdef XDEBUG
    if (Norm(m-m2) > 0.001*Norm(Q)*Norm(m0)) {
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
  template <> inline void LapQ_LDivEq(const GenMatrix<double>& Q,
      const GenVector<double>& beta, const MatrixView<double>& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == x.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int ldx = x.isrm() ? x.stepi() : x.stepj();
    int m = x.isrm() ? x.rowsize() : x.colsize();
    int n = x.isrm() ? x.colsize() : x.rowsize();
    int k = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = x.rowsize()*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
#endif
    if (Q.isrm()) {
      int ldq = Q.stepi();
      LAPNAME(dormlq) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
	  x.isrm()?LAPCH_T:LAPCH_NT,
	  LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
      LAP_Results("dormlq");
#else
      LAP_Results(int(work[0]),m,n,lwork,"dormlq");
#endif
    } else {
      int ldq = Q.stepj();
      LAPNAME(dormqr) (LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
	  x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
	  LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
      LAP_Results("dormqr");
#else
      LAP_Results(int(work[0]),m,n,lwork,"dormqr");
#endif
    }
  }
  template <> inline void LapQ_LDivEq(
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
#ifndef LAPNOWORK
    int lwork = x.rowsize()*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
#endif
    if (Q.isrm()) {
      int ldx = x.isrm() ? x.stepi() : x.stepj();
      int m = x.isrm() ? x.rowsize() : x.colsize();
      int n = x.isrm() ? x.colsize() : x.rowsize();
      int ldq = Q.stepi();
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
      LAPNAME(zunmlq) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
	  x.isrm()?LAPCH_CT:LAPCH_NT,
	  LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("zunmlq");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"zunmlq");
#endif
    } else {
      int ldx = x.isrm() ? x.stepi() : x.stepj();
      int m = x.isrm() ? x.rowsize() : x.colsize();
      int n = x.isrm() ? x.colsize() : x.rowsize();
      int ldq = Q.stepj();
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
      Vector<complex<double> > conjbeta = beta.Conjugate();
      LAPNAME(zunmqr) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
	  x.isrm()?LAPCH_NT:LAPCH_CT,
	  LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("zunmqr");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"zunmqr");
#endif
    }
  }
#ifndef NOFLOAT
  template <> inline void LapQ_LDivEq(const GenMatrix<float>& Q,
      const GenVector<float>& beta, const MatrixView<float>& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == x.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int ldx = x.isrm() ? x.stepi() : x.stepj();
    int m = x.isrm() ? x.rowsize() : x.colsize();
    int n = x.isrm() ? x.colsize() : x.rowsize();
    int k = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = x.rowsize()*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
#endif
    if (Q.isrm()) {
      int ldq = Q.stepi();
      LAPNAME(sormlq) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
	  x.isrm()?LAPCH_T:LAPCH_NT,
	  LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
      LAP_Results("sormlq");
#else
      LAP_Results(int(work[0]),m,n,lwork,"sormlq");
#endif
    } else {
      int ldq = Q.stepj();
      LAPNAME(sormqr) (LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
	  x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
	  LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
      LAP_Results("sormqr");
#else
      LAP_Results(int(work[0]),m,n,lwork,"sormqr");
#endif
    }
  }
 
  template <> inline void LapQ_LDivEq(
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
#ifndef LAPNOWORK
    int lwork = x.rowsize()*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
#endif
    if (Q.isrm()) {
      int ldx = x.isrm() ? x.stepi() : x.stepj();
      int m = x.isrm() ? x.rowsize() : x.colsize();
      int n = x.isrm() ? x.colsize() : x.rowsize();
      int ldq = Q.stepi();
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
      LAPNAME(cunmlq) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
	  x.isrm()?LAPCH_CT:LAPCH_NT,
	  LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("cunmlq");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"cunmlq");
#endif
    } else {
      int ldx = x.isrm() ? x.stepi() : x.stepj();
      int m = x.isrm() ? x.rowsize() : x.colsize();
      int n = x.isrm() ? x.colsize() : x.rowsize();
      int ldq = Q.stepj();
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
      Vector<complex<float> > conjbeta = beta.Conjugate();
      LAPNAME(cunmqr) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
	  x.isrm()?LAPCH_NT:LAPCH_CT,
	  LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("cunmqr");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"cunmqr");
#endif
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

  template <class T, class T1> inline void NonBlockQ_RDivEq(
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

  template <class T, class T1> inline void BlockQ_RDivEq(
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
    UpperTriMatrix<T1,NonUnitDiag,ColMajor> BaseZ(min(size_t(QR_BLOCKSIZE),N));
    for(size_t j2=N;j2>0;) {
      size_t j1 = j2 > QR_BLOCKSIZE ? j2-QR_BLOCKSIZE : 0;
      ConstMatrixView<T1> Y = Q.SubMatrix(j1,M,j1,j2);
      UpperTriMatrixView<T1> Z = BaseZ.SubTriMatrix(0,Y.rowsize());
      BlockHouseholder_MakeZ(Y,Z,beta.SubVector(j1,j2));
      BlockHouseholder_LMult(Y,Z,m.Cols(j1,M).Adjoint());
      j2 = j1;
    }
#ifdef XDEBUG
    if (Norm(m-m2) > 0.001*Norm(Q)*Norm(m0)) {
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
  template <> inline void LapQ_RDivEq(const GenMatrix<double>& Q,
      const GenVector<double>& beta, const MatrixView<double>& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == x.rowsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    int ldx = x.isrm() ? x.stepi() : x.stepj();
    int m = x.isrm() ? x.rowsize() : x.colsize();
    int n = x.isrm() ? x.colsize() : x.rowsize();
    int k = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = x.colsize()*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
#endif
    if (Q.isrm()) {
      int ldq = Q.stepi();
      LAPNAME(dormlq) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
	  x.isrm()?LAPCH_T:LAPCH_NT,
	  LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
      LAP_Results("dormlq");
#else
      LAP_Results(int(work[0]),m,n,lwork,"dormlq");
#endif
    } else {
      int ldq = Q.stepj();
      LAPNAME(dormqr) (LAPCM x.isrm()?LAPCH_L:LAPCH_R, 
	  x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
	  LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
      LAP_Results("dormqr");
#else
      LAP_Results(int(work[0]),m,n,lwork,"dormqr");
#endif
    }
  }
  template <> inline void LapQ_RDivEq(const GenMatrix<complex<double> >& Q,
      const GenVector<complex<double> >& beta,
      const MatrixView<complex<double> >& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == x.rowsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int k = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = x.colsize()*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
#endif
    if (Q.isrm()) {
      int ldx = x.isrm() ? x.stepi() : x.stepj();
      int m = x.isrm() ? x.rowsize() : x.colsize();
      int n = x.isrm() ? x.colsize() : x.rowsize();
      int ldq = Q.stepi();
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
      LAPNAME(zunmlq) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
	  x.isrm()?LAPCH_CT:LAPCH_NT,
	  LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("zunmlq");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"zunmlq");
#endif
    } else {
      int ldx = x.isrm() ? x.stepi() : x.stepj();
      int m = x.isrm() ? x.rowsize() : x.colsize();
      int n = x.isrm() ? x.colsize() : x.rowsize();
      int ldq = Q.stepj();
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
      Vector<complex<double> > conjbeta = beta.Conjugate();
      LAPNAME(zunmqr) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
	  x.isrm()?LAPCH_NT:LAPCH_CT,
	  LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("zunmqr");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"zunmqr");
#endif
    }
  }
#ifndef NOFLOAT
  template <> inline void LapQ_RDivEq(const GenMatrix<float>& Q,
      const GenVector<float>& beta, const MatrixView<float>& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == x.rowsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    int ldx = x.isrm() ? x.stepi() : x.stepj();
    int m = x.isrm() ? x.rowsize() : x.colsize();
    int n = x.isrm() ? x.colsize() : x.rowsize();
    int k = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = x.colsize()*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
#endif
    if (Q.isrm()) {
      int ldq = Q.stepi();
      LAPNAME(sormlq) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
	  x.isrm()?LAPCH_T:LAPCH_NT,
	  LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
      LAP_Results("sormlq");
#else
      LAP_Results(int(work[0]),m,n,lwork,"sormlq");
#endif
    } else {
      int ldq = Q.stepj();
      LAPNAME(sormqr) (LAPCM x.isrm()?LAPCH_L:LAPCH_R, 
	  x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
	  LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
      LAP_Results("sormqr");
#else
      LAP_Results(int(work[0]),m,n,lwork,"sormqr");
#endif
    }
  }
  template <> inline void LapQ_RDivEq(const GenMatrix<complex<float> >& Q,
      const GenVector<complex<float> >& beta,
      const MatrixView<complex<float> >& x)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == x.rowsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int k = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = x.colsize()*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
#endif
    if (Q.isrm()) {
      int ldx = x.isrm() ? x.stepi() : x.stepj();
      int m = x.isrm() ? x.rowsize() : x.colsize();
      int n = x.isrm() ? x.colsize() : x.rowsize();
      int ldq = Q.stepi();
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
      LAPNAME(cunmlq) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
	  x.isrm()?LAPCH_CT:LAPCH_NT,
	  LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("cunmlq");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"cunmlq");
#endif
    } else {
      int ldx = x.isrm() ? x.stepi() : x.stepj();
      int m = x.isrm() ? x.rowsize() : x.colsize();
      int n = x.isrm() ? x.colsize() : x.rowsize();
      int ldq = Q.stepj();
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
      Vector<complex<float> > conjbeta = beta.Conjugate();
      LAPNAME(cunmqr) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
	  x.isrm()?LAPCH_NT:LAPCH_CT,
	  LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
	  LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
	  LAPWK(work) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("cunmqr");
#else
      LAP_Results(int(real(work[0])),m,n,lwork,"cunmqr");
#endif
    }
  }
#endif // FLOAT
#endif // LAP

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

#define InstFile "TMV_QRDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


