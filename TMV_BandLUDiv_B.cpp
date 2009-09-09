
#include "TMV_Band.h"

//#define XDEBUG

namespace tmv {

  //
  // LDivEq
  //

  template <class T, class T1> inline void NonLapBandLU_LDivEq(
      const GenBandMatrix<T1>& LUx,
      const size_t* p, const MatrixView<T>& m) 
  { 
    // Solve A x = m given that A = L U
    // L U x = m
    TMVAssert(m.colsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    const size_t N = LUx.colsize();
    const size_t nlo = LUx.nlo();

    // Solve L y = m by forward substitution
    // Remember L is really:
    //
    // L = (  1   0 ) ( 1  0   0 ) ( I  0   0 )
    //     ( P0L0 I ) ( 0  1   0 ) ( 0  1   0 ) ...
    //                ( 0 P1L1 I ) ( 0 P2L2 I )
    //
    // where the Li are columns of nlo length which are
    // stored in the lower band of LUx,
    // and each Pi is a row swap of i with p[i]
    //
    if (nlo > 0) {
      size_t jn=nlo+1;  // jn = j+nlo+1
      const size_t* pj = p;
      for(size_t j=0; j+1<N; ++j,++pj) {
	TMVAssert(*pj<m.colsize());
	m.SwapRows(j,*pj);
	m.Rows(j+1,jn) -= LUx.col(j,j+1,jn) ^ m.row(j);
	if (jn<N) ++jn;
      }
    }

    // Next solve U x = y by back substitution
    BandTriLDivEq(BandMatrixViewOf(LUx,0,LUx.nhi()),m,NonUnitDiag);
  }

#ifdef LAP
  template <class T, class T1> inline void LapBandLU_LDivEq(
      const GenBandMatrix<T1>& LUx, const size_t* P, const MatrixView<T>& m)
  { NonLapBandLU_LDivEq(LUx,P,m); }
  template <> inline void LapBandLU_LDivEq(
      const GenBandMatrix<double>& LUx,
      const size_t* P, const MatrixView<double>& m) 
  {
    TMVAssert(m.colsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() > LUx.nlo());
    TMVAssert(LUx.nlo() > 0);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(LUx.ct() == NonConj);
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.iscm());
    TMVAssert(m.iscm());

    int n = LUx.colsize();
    int kl = LUx.nlo();
    int ku = LUx.nhi()-LUx.nlo();
    int nrhs = m.rowsize();
    int lda = LUx.stepj()+1;
    int ldm = m.stepj();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
    LAPNAME(dgbtrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
	LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
	LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("dgbtrs");
  }
  template <> inline void LapBandLU_LDivEq(
      const GenBandMatrix<complex<double> >& LUx,
      const size_t* P, const MatrixView<complex<double> >& m)
  {
    TMVAssert(m.colsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() > LUx.nlo());
    TMVAssert(LUx.nlo() > 0);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(LUx.ct() == NonConj);
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.iscm());
    TMVAssert(m.iscm());

    int n = LUx.colsize();
    int kl = LUx.nlo();
    int ku = LUx.nhi()-LUx.nlo();
    int nrhs = m.rowsize();
    int lda = LUx.diagstep();
    int ldm = m.stepj();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
    LAPNAME(zgbtrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
	LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
	LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("zgbtrs");
  }
  template <class T, class T1> inline void LapTriDiagLU_LDivEq(
      const GenBandMatrix<T1>& LUx, const size_t* P, const MatrixView<T>& m)
  { NonLapBandLU_LDivEq(LUx,P,m); }
  template <> inline void LapTriDiagLU_LDivEq(
      const GenBandMatrix<double>& LUx,
      const size_t* P, const MatrixView<double>& m) 
  {
    TMVAssert(m.colsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() == 2);
    TMVAssert(LUx.nlo() == 1);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(LUx.ct() == NonConj);
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.isdm());
    TMVAssert(m.iscm());

    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int ldm = m.stepj();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
    LAPNAME(dgttrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
	LAPP(LUx.diag(-1).cptr()),LAPP(LUx.diag().cptr()),
	LAPP(LUx.diag(1).cptr()),LAPP(LUx.diag(2).cptr()),
	LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("dgttrs");
  }
  template <> inline void LapTriDiagLU_LDivEq(
      const GenBandMatrix<complex<double> >& LUx,
      const size_t* P, const MatrixView<complex<double> >& m)
  {
    TMVAssert(m.colsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() == 2);
    TMVAssert(LUx.nlo() == 1);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(LUx.ct() == NonConj);
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.isdm());
    TMVAssert(m.iscm());

    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int ldm = m.stepj();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
    LAPNAME(zgttrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
	LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
	LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()), 
	LAPP(lap_p.get()), LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("zgttrs");
  }
#ifndef NOFLOAT
  template <> inline void LapBandLU_LDivEq(
      const GenBandMatrix<float>& LUx,
      const size_t* P, const MatrixView<float>& m) 
  {
    TMVAssert(m.colsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() > LUx.nlo());
    TMVAssert(LUx.nlo() > 0);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(LUx.ct() == NonConj);
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.iscm());
    TMVAssert(m.iscm());

    int n = LUx.colsize();
    int kl = LUx.nlo();
    int ku = LUx.nhi()-LUx.nlo();
    int nrhs = m.rowsize();
    int lda = LUx.stepj()+1;
    int ldm = m.stepj();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
    LAPNAME(sgbtrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
	LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
	LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("sgbtrs");
  }
  template <> inline void LapBandLU_LDivEq(
      const GenBandMatrix<complex<float> >& LUx,
      const size_t* P, const MatrixView<complex<float> >& m)
  {
    TMVAssert(m.colsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() > LUx.nlo());
    TMVAssert(LUx.nlo() > 0);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(LUx.ct() == NonConj);
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.iscm());
    TMVAssert(m.iscm());

    int n = LUx.colsize();
    int kl = LUx.nlo();
    int ku = LUx.nhi()-LUx.nlo();
    int nrhs = m.rowsize();
    int lda = LUx.stepj()+1;
    int ldm = m.stepj();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
    LAPNAME(cgbtrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
	LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
	LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("cgbtrs");
  }
  template <> inline void LapTriDiagLU_LDivEq(
      const GenBandMatrix<float>& LUx,
      const size_t* P, const MatrixView<float>& m) 
  {
    TMVAssert(m.colsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() == 2);
    TMVAssert(LUx.nlo() == 1);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(LUx.ct() == NonConj);
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.isdm());
    TMVAssert(m.iscm());

    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int ldm = m.stepj();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
    LAPNAME(sgttrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
	LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
	LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()),
	LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("sgttrs");
  }
  template <> inline void LapTriDiagLU_LDivEq(
      const GenBandMatrix<complex<float> >& LUx,
      const size_t* P, const MatrixView<complex<float> >& m)
  {
    TMVAssert(m.colsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() == 2);
    TMVAssert(LUx.nlo() == 1);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(LUx.ct() == NonConj);
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.isdm());
    TMVAssert(m.iscm());

    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int ldm = m.stepj();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
    LAPNAME(cgttrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
	LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
	LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()),
	LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("cgttrs");
  }
#endif // FLOAT
#endif // LAP

  template <class T, class T1> void BandLU_LDivEq(
      const GenBandMatrix<T1>& LUx,
      const size_t* P, const MatrixView<T>& m) 
  {
    TMVAssert(m.colsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());

    if (m.colsize() > 0 && m.rowsize() > 0) {
#ifdef LAP
      if (m.iscm() && LUx.iscm() && !LUx.isconj() && LUx.nlo() > 0)
	LapBandLU_LDivEq(LUx,P,m);
      else if (m.iscm() && LUx.isdm() && LUx.nlo() == 1 && LUx.nhi() == 2 && 
	  !LUx.isconj())
	LapTriDiagLU_LDivEq(LUx,P,m);
      else
#endif
	NonLapBandLU_LDivEq(LUx,P,m);
    }
  }

  //
  // RDivEq
  //

  template <class T, class T1> inline void NonLapBandLU_RDivEq(
      const GenBandMatrix<T1>& LUx,
      const size_t* p, const MatrixView<T>& m) 
  { 
    // Solve x A = m given that A = L U
    // x L U = m
    TMVAssert(m.rowsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());

    const size_t N = LUx.colsize();
    const size_t nlo = LUx.nlo();

    // First solve y U = m by forward substitution
    // Or: UT yT = mT
    BandTriLDivEq(Transpose(BandMatrixViewOf(LUx,0,LUx.nhi())),
	Transpose(m),NonUnitDiag);

    // Next solve z L = y by back substitution with L = :
    //
    // L = (  1   0 ) ( 1  0   0 ) ( I  0   0 )     ( I    0     0 ) ( I 0 )
    //     ( P0L0 I ) ( 0  1   0 ) ( 0  1   0 ) ... ( 0    1     0 ) ( 0 1 )
    //                ( 0 P1L1 I ) ( 0 P2L2 I )     ( 0 Pn-1Ln-1 1 )
    //
    if (nlo > 0) {
      size_t jn=N;
      size_t k=nlo-1;
      const size_t* pj = p+N-1;
      for(size_t j=N-1;j>0;) {
	--j; --pj;
	m.col(j) -= m.Cols(j+1,jn) * LUx.col(j,j+1,jn);
	TMVAssert(*pj<m.rowsize());
	m.SwapCols(j,*pj);
	if (k>0) --k; else --jn;
      }
    }
  }

#ifdef LAP
  template <class T, class T1> inline void LapBandLU_RDivEq(
      const GenBandMatrix<T1>& LUx, const size_t* P, const MatrixView<T>& m)
  { NonLapBandLU_RDivEq(LUx,P,m); }
  template <> inline void LapBandLU_RDivEq(
      const GenBandMatrix<double>& LUx,
      const size_t* P, const MatrixView<double>& m) 
  {
    TMVAssert(m.rowsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() > LUx.nlo());
    TMVAssert(LUx.nlo() > 0);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(LUx.ct() == NonConj);
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.iscm());
    TMVAssert(m.isrm());

    int n = LUx.colsize();
    int kl = LUx.nlo();
    int ku = LUx.nhi()-LUx.nlo();
    int nrhs = m.colsize();
    int lda = LUx.stepj()+1;
    int ldm = m.stepi();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
    LAPNAME(dgbtrs) (LAPCM LAPCH_T,LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
	LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
	LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("dgbtrs");
  }
  template <> inline void LapBandLU_RDivEq(
      const GenBandMatrix<complex<double> >& LUx,
      const size_t* P, const MatrixView<complex<double> >& m)
  {
    TMVAssert(m.rowsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() > LUx.nlo());
    TMVAssert(LUx.nlo() > 0);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.iscm());
    TMVAssert(m.isrm());

    int n = LUx.colsize();
    int kl = LUx.nlo();
    int ku = LUx.nhi()-LUx.nlo();
    int nrhs = m.colsize();
    int lda = LUx.stepj()+1;
    int ldm = m.stepi();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
    LAPNAME(zgbtrs) (LAPCM LUx.isconj()?LAPCH_CT:LAPCH_T,
	LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
	LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
	LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("zgbtrs");
  }
  template <class T, class T1> inline void LapTriDiagLU_RDivEq(
      const GenBandMatrix<T1>& LUx, const size_t* P, const MatrixView<T>& m)
  { NonLapBandLU_RDivEq(LUx,P,m); }
  template <> inline void LapTriDiagLU_RDivEq(
      const GenBandMatrix<double>& LUx,
      const size_t* P, const MatrixView<double>& m) 
  {
    TMVAssert(m.rowsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() == 2);
    TMVAssert(LUx.nlo() == 1);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(LUx.ct() == NonConj);
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.isdm());
    TMVAssert(m.isrm());

    int n = LUx.colsize();
    int nrhs = m.colsize();
    int ldm = m.stepi();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
    LAPNAME(dgttrs) (LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
	LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()),
	LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("dgttrs");
  }
  template <> inline void LapTriDiagLU_RDivEq(
      const GenBandMatrix<complex<double> >& LUx,
      const size_t* P, const MatrixView<complex<double> >& m)
  {
    TMVAssert(m.rowsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() == 2);
    TMVAssert(LUx.nlo() == 1);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.isdm());
    TMVAssert(m.isrm());

    int n = LUx.colsize();
    int nrhs = m.colsize();
    int ldm = m.stepi();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
    LAPNAME(zgttrs) (LAPCM LUx.isconj()?LAPCH_CT:LAPCH_T, LAPV(n),LAPV(nrhs),
	LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
	LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()),
	LAPP(lap_p.get()), LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("zgttrs");
  }
#ifndef NOFLOAT
  template <> inline void LapBandLU_RDivEq(
      const GenBandMatrix<float>& LUx,
      const size_t* P, const MatrixView<float>& m) 
  {
    TMVAssert(m.rowsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() > LUx.nlo());
    TMVAssert(LUx.nlo() > 0);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(LUx.ct() == NonConj);
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.iscm());
    TMVAssert(m.isrm());

    int n = LUx.colsize();
    int kl = LUx.nlo();
    int ku = LUx.nhi()-LUx.nlo();
    int nrhs = m.colsize();
    int lda = LUx.stepj()+1;
    int ldm = m.stepi();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
    LAPNAME(sgbtrs) (LAPCM LAPCH_T,LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
	LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
	LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("sgbtrs");
  }
  template <> inline void LapBandLU_RDivEq(
      const GenBandMatrix<complex<float> >& LUx,
      const size_t* P, const MatrixView<complex<float> >& m)
  {
    TMVAssert(m.rowsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() > LUx.nlo());
    TMVAssert(LUx.nlo() > 0);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.iscm());
    TMVAssert(m.isrm());

    int n = LUx.colsize();
    int kl = LUx.nlo();
    int ku = LUx.nhi()-LUx.nlo();
    int nrhs = m.colsize();
    int lda = LUx.stepj()+1;
    int ldm = m.stepi();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
    LAPNAME(cgbtrs) (LAPCM LUx.isconj()?LAPCH_CT:LAPCH_T,
	LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
	LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
	LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("cgbtrs");
  }
  template <> inline void LapTriDiagLU_RDivEq(
      const GenBandMatrix<float>& LUx,
      const size_t* P, const MatrixView<float>& m) 
  {
    TMVAssert(m.rowsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() == 2);
    TMVAssert(LUx.nlo() == 1);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(LUx.ct() == NonConj);
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.isdm());
    TMVAssert(m.isrm());

    int n = LUx.colsize();
    int nrhs = m.colsize();
    int ldm = m.stepi();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
    LAPNAME(sgttrs) (LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
	LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()),
	LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("sgttrs");
  }
  template <> inline void LapTriDiagLU_RDivEq(
      const GenBandMatrix<complex<float> >& LUx,
      const size_t* P, const MatrixView<complex<float> >& m)
  {
    TMVAssert(m.rowsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(m.rowsize() > 0);
    TMVAssert(m.colsize() > 0);
    TMVAssert(LUx.nhi() == 2);
    TMVAssert(LUx.nlo() == 1);
    TMVAssert(LUx.nhi() < int(LUx.colsize()));
    TMVAssert(LUx.nlo() < int(LUx.colsize()));
    TMVAssert(m.ct() == NonConj);
    TMVAssert(LUx.isdm());
    TMVAssert(m.isrm());

    int n = LUx.colsize();
    int nrhs = m.colsize();
    int ldm = m.stepi();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
    LAPNAME(cgttrs) (LAPCM LUx.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
	LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()),
	LAPP(lap_p.get()), LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
    LAP_Results("cgttrs");
  }
#endif
#endif // LAP

  template <class T, class T1> void BandLU_RDivEq(
      const GenBandMatrix<T1>& LUx,
      const size_t* P, const MatrixView<T>& m) 
  {
    TMVAssert(m.rowsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());

    if (m.colsize() > 0 && m.rowsize() > 0) {
#ifdef LAP
      if (m.isrm() && LUx.iscm() && LUx.nlo()>0)
	LapBandLU_RDivEq(LUx,P,m);
      else if (m.isrm() && LUx.isdm() && LUx.nlo() == 1 && LUx.nhi() == 2)
	LapTriDiagLU_RDivEq(LUx,P,m);
      else
#endif
	NonLapBandLU_RDivEq(LUx,P,m);
    }
  }

  //
  // BandTriLDivEq
  //

  template <bool rm, bool ca, bool ua, class T, class Ta> 
    inline void DoRowUpperBandTriLDivEq(
	const GenBandMatrix<Ta>& A, const VectorView<T>& b)
    {
      TMVAssert(b.step()==1);
      TMVAssert(A.IsSquare());
      TMVAssert(b.size() == A.colsize());
      TMVAssert(b.size() > 0);
      TMVAssert(b.ct() == NonConj);
      TMVAssert(A.nlo() == 0);
      TMVAssert(A.nhi() > 0);
      TMVAssert(rm == A.isrm());
      TMVAssert(ca == A.isconj());

      const size_t N = b.size();

      const int sj = (rm?1:A.stepj());
      const int ds = A.diagstep();
      const Ta* Aii = A.cptr() + (ua ? N-2 : N-1)*ds;
      T* bi = b.ptr() + (ua ? N-2 : N-1);

      if (!ua) {
	if (*Aii==Ta(0)) 
	  throw SingularBandLU<Ta>(A);
	*bi /= (ca ? CONJ(*Aii) : *Aii);
	Aii -= ds;
	--bi;
      }
      if (N==1) return;

      size_t k = A.nhi()-1;
      for(size_t i=N-1,len=1;i>0;--i,Aii-=ds,--bi) {
	// Actual row being done is i-1, not i

	// *bi -= A.row(i,i+1,j2) * b.SubVector(i+1,j2);
	const T* bj = bi+1;
	const Ta* Aij = Aii+sj;
	for(size_t j=len;j>0;--j,++bj,(rm?++Aij:Aij+=sj))
	  *bi -= *bj * (ca ? CONJ(*Aij) : *Aij);

	if (!ua) {
	  if (*Aii==Ta(0)) 
	    throw SingularBandLU<Ta>(A);
	  *bi /= (ca ? CONJ(*Aii) : *Aii);
	}
	if (k > 0) { --k; ++len; }
      } 
    }

  template <bool rm, bool ua, class T, class Ta> 
    inline void RowUpperBandTriLDivEq(
	const GenBandMatrix<Ta>& A, const VectorView<T>& b)
    {
      if (A.isconj())
	DoRowUpperBandTriLDivEq<rm,true,ua>(A,b);
      else
	DoRowUpperBandTriLDivEq<rm,false,ua>(A,b);
    }
                  
  template <bool cm, bool ca, bool ua, class T, class Ta> 
    inline void DoColUpperBandTriLDivEq(
	const GenBandMatrix<Ta>& A, const VectorView<T>& b)
    {
      TMVAssert(b.step()==1);
      TMVAssert(A.IsSquare());
      TMVAssert(b.size() == A.colsize());
      TMVAssert(b.size() > 0);
      TMVAssert(b.ct() == NonConj);
      TMVAssert(A.nlo() == 0);
      TMVAssert(A.nhi() > 0);
      TMVAssert(cm == A.iscm());
      TMVAssert(ca == A.isconj());

      const size_t N = b.size();

      const int si = (cm ? 1 : A.stepi());
      const int sj = A.stepj();
      const int ds = A.diagstep();
      const size_t hi = A.nhi();

      size_t i1 = N-1;
      if (i1 > hi) i1 -= hi;
      else i1 = 0;

      const Ta* Ai1j = A.cptr()+(N-1)*sj+i1*si;
      const Ta* Ajj = (ua ? 0 : A.cptr()+(N-1)*ds); // if unit, this isn't used
      T* bi1 = b.ptr()+i1;
      T* bj = b.ptr()+N-1;

      for(size_t len=N-1-i1;len>0;--bj) {
	if (*bj != T(0)) {
	  if (!ua) {
	    if (*Ajj==Ta(0)) 
	      throw SingularBandLU<Ta>(A);
	    *bj /= (ca ? CONJ(*Ajj) : *Ajj);
	    Ajj -= ds;
	  }

	  // b.SubVector(i1,j) -= (*bj) * A.col(j,i1,j);
	  T* bi = bi1;
	  const Ta* Aij = Ai1j;
	  for(size_t i=len;i>0;--i,++bi,(cm?++Aij:Aij+=si))
	    *bi -= *bj * (ca ? CONJ(*Aij) : *Aij);
	}
	else if (!ua) Ajj -= ds;

	if (i1 > 0) { --i1; --bi1; Ai1j-=ds; } 
	else { --len; Ai1j-=sj; }
      } 
      if (!ua && *bj != T(0)) {
	if (*Ajj==Ta(0)) 
	  throw SingularBandLU<Ta>(A);
	*bj /= (ca ? CONJ(*Ajj) : *Ajj);
      }
    }

  template <bool cm, bool ua, class T, class Ta> 
    inline void ColUpperBandTriLDivEq(
	const GenBandMatrix<Ta>& A, const VectorView<T>& b)
    {
      if (A.isconj())
	DoColUpperBandTriLDivEq<cm,true,ua>(A,b);
      else
	DoColUpperBandTriLDivEq<cm,false,ua>(A,b);
    }
                  
  template <bool rm, bool ca, bool ua, class T, class Ta> 
    inline void DoRowLowerBandTriLDivEq(
	const GenBandMatrix<Ta>& A, const VectorView<T>& b)
    {
      TMVAssert(b.step()==1);
      TMVAssert(A.IsSquare());
      TMVAssert(b.size() == A.colsize());
      TMVAssert(b.size() > 0);
      TMVAssert(b.ct() == NonConj);
      TMVAssert(A.nhi() == 0);
      TMVAssert(A.nlo() > 0);
      TMVAssert(rm == A.isrm());
      TMVAssert(ca == A.isconj());

      const size_t N = A.colsize();

      const int sj = (rm ? 1 : A.stepj());
      const int si = A.stepi();
      const int ds = A.diagstep();

      const Ta* Aij1 = A.cptr();
      const T* bj1 = b.cptr();
      T* bi = b.ptr();

      if (!ua) {
	if (*Aij1==Ta(0)) 
	  throw SingularBandLU<Ta>(A);
	*bi /= (ca ? CONJ(*Aij1) : *Aij1);
      }

      ++bi;
      Aij1 += si;
      size_t k=A.nlo()-1;

      for(size_t i=1,len=1;i<N;++i,++bi) {
	// *bi -= A.row(i,j1,i) * b.SubVector(j1,i);
	const Ta* Aij = Aij1;
	const T* bj = bj1;
	for(size_t j=len;j>0;--j,++bj,(rm?++Aij:Aij+=sj))
	  *bi -= *bj * (ca ? CONJ(*Aij) : *Aij);
	if (!ua) {
	  // Aij is Aii after the above for loop
	  if (*Aij == Ta(0))
	    throw SingularBandLU<Ta>(A);
	  *bi /= (ca ? CONJ(*Aij) : *Aij);
	}
	if (k>0) { --k; ++len; Aij1+=A.stepi(); } 
	else { ++bj1; Aij1+=ds; }
      }
    }

  template <bool rm, bool ua, class T, class Ta> 
    inline void RowLowerBandTriLDivEq(
	const GenBandMatrix<Ta>& A, const VectorView<T>& b)
    {
      if (A.isconj())
	DoRowLowerBandTriLDivEq<rm,true,ua>(A,b);
      else
	DoRowLowerBandTriLDivEq<rm,false,ua>(A,b);
    }
                  
  template <bool cm, bool ca, bool ua, class T, class Ta> 
    inline void DoColLowerBandTriLDivEq(
	const GenBandMatrix<Ta>& A, const VectorView<T>& b)
    {
      TMVAssert(b.step()==1);
      TMVAssert(A.IsSquare());
      TMVAssert(b.size() == A.colsize());
      TMVAssert(b.size() > 0);
      TMVAssert(b.ct() == NonConj);
      TMVAssert(A.nhi() == 0);
      TMVAssert(A.nlo() > 0);
      TMVAssert(cm == A.iscm());
      TMVAssert(ca == A.isconj());

      const size_t N = A.colsize();

      const int si = (cm ? 1 : A.stepi());
      const int ds = A.diagstep();

      const Ta* Ajj= A.cptr();
      T* bj = b.ptr();

      size_t i2=min(size_t(A.nlo())+1,A.colsize());

      for(size_t len=i2-1;len>0;++bj,Ajj+=ds) {
	if (*bj != T(0)) {
	  if (!ua) {
	    if (*Ajj==Ta(0)) 
	      throw SingularBandLU<Ta>(A);
	    *bj /= (ca ? CONJ(*Ajj) : *Ajj);
	  }
	  // b.SubVector(j+1,i2) -= *bj * A.col(j,j+1,i2)
	  T* bi = bj+1;
	  const Ta* Aij = Ajj+si;
	  for(size_t i=len;i>0;--i,++bi,(cm?++Aij:Aij+=si))
	    *bi -= *bj * (ca ? CONJ(*Aij) : *Aij);
	}
	if (i2 < N) ++i2; else --len;
      }
      if (!ua && *bj != T(0)) {
	if (*Ajj==Ta(0)) 
	  throw SingularBandLU<Ta>(A);
	*bj /= (ca ? CONJ(*Ajj) : *Ajj);
      }
    }

  template <bool rm, bool ua, class T, class Ta> 
    inline void ColLowerBandTriLDivEq(
	const GenBandMatrix<Ta>& A, const VectorView<T>& b)
    {
      if (A.isconj())
	DoColLowerBandTriLDivEq<rm,true,ua>(A,b);
      else
	DoColLowerBandTriLDivEq<rm,false,ua>(A,b);
    }
                  
  template <class T, class Ta> inline void UpperBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const VectorView<T>& b, DiagType dt)
  {
    if (dt == UnitDiag)
      if (A.isrm()) RowUpperBandTriLDivEq<true,true>(A,b);
      else if (A.iscm()) ColUpperBandTriLDivEq<true,true>(A,b);
      else RowUpperBandTriLDivEq<false,true>(A,b);
    else
      if (A.isrm()) RowUpperBandTriLDivEq<true,false>(A,b);
      else if (A.iscm()) ColUpperBandTriLDivEq<true,false>(A,b);
      else RowUpperBandTriLDivEq<false,false>(A,b);
  }

  template <class T, class Ta> inline void LowerBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const VectorView<T>& b, DiagType dt)
  {
    if (dt == UnitDiag)
      if (A.isrm()) RowLowerBandTriLDivEq<true,true>(A,b);
      else if (A.iscm()) ColLowerBandTriLDivEq<true,true>(A,b);
      else RowLowerBandTriLDivEq<false,true>(A,b);
    else
      if (A.isrm()) RowLowerBandTriLDivEq<true,false>(A,b);
      else if (A.iscm()) ColLowerBandTriLDivEq<true,false>(A,b);
      else RowLowerBandTriLDivEq<false,false>(A,b);
  }

  template <class T, class Ta> inline void NonBlasBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const VectorView<T>& b, DiagType dt)
  {
    // Solve A x = y  where A is an upper or lower band triangle matrix
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(b.size() > 0);
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(b.ct() == NonConj);

    const size_t N = b.size();
    if (b.step() == 1) {
      if (A.nlo() == 0) {
	size_t i2 = N;
	for(const T* b2 = b.cptr()+i2-1; i2>0 && *b2==T(0); --i2,--b2);
	if (i2==0) return;
	else if (i2 == 1) {
	  if (dt == NonUnitDiag) b(0) /= A(0,0);
	}
	else if (i2 == N) 
	  UpperBandTriLDivEq(A,b,dt);
	else if (size_t(A.nhi()) < i2)
	  UpperBandTriLDivEq(A.SubBandMatrix(0,i2,0,i2,0,A.nhi()),
	      b.SubVector(0,i2),dt);
	else
	  UpperBandTriLDivEq(A.SubBandMatrix(0,i2,0,i2,0,i2-1),
	      b.SubVector(0,i2),dt);
      } else {
	size_t i1 = 0;
	for(const T* b1 = b.cptr(); i1<N && *b1==T(0); ++i1,++b1);
	if (i1==N) return;
	else if (i1 == 0) 
	  LowerBandTriLDivEq(A,b,dt);
	else if (i1 == N-1) {
	  if (dt == NonUnitDiag) b(N-1) /= A(N-1,N-1);
	}
	else if (size_t(A.nlo()) < (N-i1))
	  LowerBandTriLDivEq(A.SubBandMatrix(i1,N,i1,N,A.nlo(),0),
	      b.SubVector(i1,N),dt);
	else
	  LowerBandTriLDivEq(A.SubBandMatrix(i1,N,i1,N,N-i1-1,0),
	      b.SubVector(i1,N),dt);
      }
    } else {
      Vector<T> bb = b;
      NonBlasBandTriLDivEq(A,bb.View(),dt);
      b = bb;
    }
  }

#ifdef BLAS
  template <class T, class Ta> inline void BlasBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const VectorView<T>& b, DiagType dt)
  { NonBlasBandTriLDivEq(A,b,dt); }
  template <> inline void BlasBandTriLDivEq(
      const GenBandMatrix<double>& A, const VectorView<double>& b, DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    int n = A.colsize();
    int kd = A.nlo()==0 ? A.nhi() : A.nlo();
    int aoffset = A.isrm() ? A.nlo() : A.nhi();
    int ds = A.diagstep();
    int s = b.step();
    BLASNAME(dtbsv) (BLASCM 
	(A.nlo()==0 == A.isrm()) ?BLASCH_LO:BLASCH_UP, 
	A.isrm() ? BLASCH_T : BLASCH_NT, 
	dt==UnitDiag ? BLASCH_U : BLASCH_NU,
	BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
	BLASP(b.ptr()), BLASV(s) BLAS1 BLAS1 BLAS1);
  }
  template <> inline void BlasBandTriLDivEq(
      const GenBandMatrix<complex<double> >& A,
      const VectorView<complex<double> >& b, DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(b.ct() == NonConj);

    int n = A.colsize();
    int kd = A.nlo()==0 ? A.nhi() : A.nlo();
    int aoffset = A.isrm() ? A.nlo() : A.nhi();
    int ds = A.diagstep();
    int s = b.step();
    if (A.iscm() && A.isconj()) {
#ifdef CBLAS
      BLASNAME(ztbsv) (BLASRM
	  A.nlo()==0 ? BLASCH_LO : BLASCH_UP, BLASCH_CT,
	  dt==UnitDiag ? BLASCH_U : BLASCH_NU, 
	  BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
	  BLASP(b.ptr()), BLASV(s) BLAS1 BLAS1 BLAS1);
#else
      b.ConjugateSelf();
      BLASNAME(ztbsv) (BLASCM A.nlo()==0?BLASCH_UP:BLASCH_LO, BLASCH_NT, 
	  dt==UnitDiag ? BLASCH_U : BLASCH_NU,
	  BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
	  BLASP(b.ptr()), BLASV(s) BLAS1 BLAS1 BLAS1);
      b.ConjugateSelf();
#endif
    } else {
      BLASNAME(ztbsv) (BLASCM 
	  (A.nlo()==0 == A.isrm()) ?BLASCH_LO:BLASCH_UP, 
	  A.isrm() ? A.isconj() ? BLASCH_CT : BLASCH_T : BLASCH_NT, 
	  dt==UnitDiag ? BLASCH_U : BLASCH_NU,
	  BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
	  BLASP(b.ptr()), BLASV(s) BLAS1 BLAS1 BLAS1);
    }
  }
#ifndef NOFLOAT
  template <> inline void BlasBandTriLDivEq(
      const GenBandMatrix<float>& A, const VectorView<float>& b, DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    int n = A.colsize();
    int kd = A.nlo()==0 ? A.nhi() : A.nlo();
    int aoffset = A.isrm() ? A.nlo() : A.nhi();
    int ds = A.diagstep();
    int s = b.step();
    BLASNAME(stbsv) (BLASCM 
	(A.nlo()==0 == A.isrm()) ?BLASCH_LO:BLASCH_UP, 
	A.isrm() ? BLASCH_T : BLASCH_NT, 
	dt==UnitDiag ? BLASCH_U : BLASCH_NU,
	BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
	BLASP(b.ptr()), BLASV(s) BLAS1 BLAS1 BLAS1);
  }
  template <> inline void BlasBandTriLDivEq(
      const GenBandMatrix<complex<float> >& A,
      const VectorView<complex<float> >& b, DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(b.ct() == NonConj);

    int n = A.colsize();
    int kd = A.nlo()==0 ? A.nhi() : A.nlo();
    int aoffset = A.isrm() ? A.nlo() : A.nhi();
    int ds = A.diagstep();
    int s = b.step();
    if (A.iscm() && A.isconj()) {
#ifdef CBLAS
      BLASNAME(ctbsv) (BLASRM
	  A.nlo()==0 ? BLASCH_LO : BLASCH_UP, BLASCH_CT,
	  dt==UnitDiag ? BLASCH_U : BLASCH_NU, 
	  BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
	  BLASP(b.ptr()), BLASV(s) BLAS1 BLAS1 BLAS1);
#else
      b.ConjugateSelf();
      BLASNAME(ctbsv) (BLASCM A.nlo()==0?BLASCH_UP:BLASCH_LO, BLASCH_NT, 
	  dt==UnitDiag ? BLASCH_U : BLASCH_NU,
	  BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
	  BLASP(b.ptr()), BLASV(s) BLAS1 BLAS1 BLAS1);
      b.ConjugateSelf();
#endif
    } else {
      BLASNAME(ctbsv) (BLASCM 
	  (A.nlo()==0 == A.isrm()) ?BLASCH_LO:BLASCH_UP, 
	  A.isrm() ? A.isconj() ? BLASCH_CT : BLASCH_T : BLASCH_NT, 
	  dt==UnitDiag ? BLASCH_U : BLASCH_NU,
	  BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
	  BLASP(b.ptr()), BLASV(s) BLAS1 BLAS1 BLAS1);
    }
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Ta> void BandTriLDivEq(
      const GenBandMatrix<Ta>& A, const VectorView<T>& b, DiagType dt)
  {
#ifdef XDEBUG
    Vector<T> b0 = b;
#endif
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.nhi() >= 0);
    TMVAssert(A.nlo() >= 0);
    TMVAssert(A.nhi() < int(A.colsize()));
    TMVAssert(A.nlo() < int(A.colsize()));
    if (b.isconj())
      BandTriLDivEq(A.Conjugate(),b.Conjugate(),dt);
    else
#ifdef BLAS
      if (A.isrm() || A.iscm())
	BlasBandTriLDivEq(A,b,dt);
      else 
#endif
	NonBlasBandTriLDivEq(A,b,dt);

#ifdef XDEBUG
    Vector<T> bb = A*b;
    if (Norm(bb-b0) > 0.001*Norm(b0)) {
      cerr<<"BandTriLDivEq Vector:\n";
      cerr<<"A = "<<Type(A)<<A<<endl;
      cerr<<"b = "<<Type(b)<<b0<<endl;
      cerr<<"--> b = "<<b<<endl;
      cerr<<"A*b = "<<bb<<endl;
      abort();
    }
#endif
  }

  template <bool ua, class T, class Ta> inline void RowUpperBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const MatrixView<T>& B)
  {
    TMVAssert(B.isrm());
    TMVAssert(A.IsSquare());
    TMVAssert(B.colsize() == A.colsize());
    TMVAssert(A.colsize() > 0);
    TMVAssert(B.rowsize() > 0);
    TMVAssert(A.nlo() == 0);
    TMVAssert(B.ct() == NonConj);

    size_t N = B.colsize();

    size_t k = A.nhi();
    if (ua) {
      for(int i=N-1; i>=0; --i) {
	B.row(i) -= A.row(i,i+1,N) * B.Rows(i+1,N);
	if (k > 0) --k; else --N;
      }
    } else {
      const int ds = A.diagstep();
      const Ta* Aii = A.cptr() + (N-1)*ds;
      for(int i=N-1; i>=0; --i,Aii-=ds) {
	B.row(i) -= A.row(i,i+1,N) * B.Rows(i+1,N);
	if (*Aii==Ta(0)) 
	  throw SingularBandLU<Ta>(A);
	B.row(i) /= (A.isconj() ? CONJ(*Aii) : *Aii);
	if (k > 0) --k; else --N;
      }
    } 
  }

  template <bool ua, class T, class Ta> inline void ColUpperBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const MatrixView<T>& B)
  {
    TMVAssert(B.isrm());
    TMVAssert(A.IsSquare());
    TMVAssert(B.colsize() == A.colsize());
    TMVAssert(A.colsize() > 0);
    TMVAssert(B.rowsize() > 0);
    TMVAssert(A.nlo() == 0);
    TMVAssert(int(A.colsize())>A.nhi());
    TMVAssert(B.ct() == NonConj);

    size_t N = A.colsize();

    size_t i1 = N-1-A.nhi();
    if (ua) {
      for(int j=N-1; j>0; --j) {
	B.Rows(i1,j) -= A.col(j,i1,j) ^ B.row(j);
	if (i1 > 0) --i1;
      }
    } else {
      const int ds = A.diagstep();
      const Ta* Ajj = A.cptr() + (N-1)*ds;
      for(int j=N-1; j>=0; --j,Ajj-=ds) {
	if (*Ajj==Ta(0)) 
	  throw SingularBandLU<Ta>(A);
	B.row(j) /= (A.isconj() ? CONJ(*Ajj) : *Ajj);
	B.Rows(i1,j) -= A.col(j,i1,j) ^ B.row(j);
	if (i1 > 0) --i1;
      }
    } 
  }

  template <bool ua, class T, class Ta> inline void RowLowerBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const MatrixView<T>& B)
  {
    TMVAssert(B.isrm());
    TMVAssert(A.IsSquare());
    TMVAssert(B.colsize() == A.colsize());
    TMVAssert(A.colsize() > 0);
    TMVAssert(B.rowsize() > 0);
    TMVAssert(A.nhi() == 0);
    TMVAssert(B.ct() == NonConj);

    const size_t N = B.colsize();

    size_t i1=0;
    size_t k=A.nlo();
    if (ua) {
      for(size_t i=0; i<N; ++i) {
	B.row(i) -= A.row(i,i1,i) * B.Rows(i1,i);
	if (k>0) --k; else ++i1;
      }
    } else {
      const int ds = A.diagstep();
      const Ta* Aii = A.cptr();
      for(size_t i=0; i<N; ++i,Aii+=ds) {
	B.row(i) -= A.row(i,i1,i) * B.Rows(i1,i);
	if (*Aii==Ta(0)) 
	  throw SingularBandLU<Ta>(A);
	B.row(i) /= (A.isconj() ? CONJ(*Aii) : *Aii);
	if (k>0) --k; else ++i1;
      }
    }
  }

  template <bool ua, class T, class Ta> inline void ColLowerBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const MatrixView<T>& B)
  {
    TMVAssert(B.isrm());
    TMVAssert(A.IsSquare());
    TMVAssert(B.colsize() == A.colsize());
    TMVAssert(A.colsize() > 0);
    TMVAssert(B.rowsize() > 0);
    TMVAssert(A.nhi() == 0);
    TMVAssert(B.ct() == NonConj);

    const size_t N = B.colsize();

    size_t i2=A.nlo()+1;
    if (ua) {
      for(size_t j=0; j<N; ++j) {
	B.Rows(j+1,i2) -= A.col(j,j+1,i2) ^ B.row(j);
	if (i2 < N) ++i2;
      }
    } else {
      const int ds = A.diagstep();
      const Ta* Ajj = A.cptr();
      for(size_t j=0; j<N; ++j,Ajj+=ds) {
	if (*Ajj==Ta(0)) 
	  throw SingularBandLU<Ta>(A);
	B.row(j) /= (A.isconj() ? CONJ(*Ajj) : *Ajj);
	B.Rows(j+1,i2) -= A.col(j,j+1,i2) ^ B.row(j);
	if (i2 < N) ++i2;
      }
    }
  }

  template <bool ua, class T, class Ta> inline void NonLapBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const MatrixView<T>& B)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(B.colsize() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.colsize() > 0);
    TMVAssert(B.rowsize() > 1);
    TMVAssert(A.nhi() >= 0);
    TMVAssert(A.nlo() >= 0);
    TMVAssert(A.nhi() < int(A.colsize()));
    TMVAssert(A.nlo() < int(A.colsize()));
    TMVAssert(B.ct() == NonConj);

    if (B.isrm()) {
      if (A.nlo()==0)
	if (A.isrm()) RowUpperBandTriLDivEq<ua>(A,B);
	else if (A.iscm()) ColUpperBandTriLDivEq<ua>(A,B);
	else RowUpperBandTriLDivEq<ua>(A,B);
      else
	if (A.isrm()) RowLowerBandTriLDivEq<ua>(A,B);
	else if (A.iscm()) ColLowerBandTriLDivEq<ua>(A,B);
	else RowLowerBandTriLDivEq<ua>(A,B);
    } else {
      for(size_t j=0;j<B.rowsize();++j) 
	BandTriLDivEq(A,B.col(j),ua?UnitDiag:NonUnitDiag);
    }
  }

#ifdef LAP
  template <class T, class Ta> inline void LapBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const MatrixView<T>& B, DiagType dt)
  { 
    if (dt == UnitDiag) NonLapBandTriLDivEq<true>(A,B);
    else NonLapBandTriLDivEq<false>(A,B);
  }
  template <> inline void LapBandTriLDivEq(
      const GenBandMatrix<double>& A, const MatrixView<double>& B, DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(B.colsize() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.colsize() > 0);
    TMVAssert(B.rowsize() > 1);
    TMVAssert(A.nhi() >= 0);
    TMVAssert(A.nlo() >= 0);
    TMVAssert(A.nhi() < int(A.colsize()));
    TMVAssert(A.nlo() < int(A.colsize()));
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(B.iscm());
    int n = A.colsize();
    int kd = A.nlo()==0 ? A.nhi() : A.nlo();
    int nrhs = B.rowsize();
    int aoffset = 
      (A.nlo()==0 && A.iscm()) ? A.nhi() :
      (A.nhi()==0 && A.isrm()) ? A.nlo() : 0;
    int lda = A.diagstep();
    int ldb = B.stepj();
    LAPNAME(dtbtrs) (LAPCM (A.iscm()?A.nlo():A.nhi())==0?LAPCH_UP:LAPCH_LO,
	A.iscm()?LAPCH_NT:LAPCH_T, dt==UnitDiag?LAPCH_U:LAPCH_NU,
	LAPV(n),LAPV(kd),LAPV(nrhs), LAPP(A.cptr()-aoffset),LAPV(lda),
	LAPP(B.ptr()),LAPV(ldb) LAPINFO LAP1 LAP1 LAP1);
    LAP_Results("dtbtrs");
  }
  template <> inline void LapBandTriLDivEq(
      const GenBandMatrix<complex<double> >& A,
      const MatrixView<complex<double> >& B, DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(B.colsize() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.colsize() > 0);
    TMVAssert(B.rowsize() > 1);
    TMVAssert(A.nhi() >= 0);
    TMVAssert(A.nlo() >= 0);
    TMVAssert(A.nhi() < int(A.colsize()));
    TMVAssert(A.nlo() < int(A.colsize()));
    TMVAssert(B.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(B.iscm());
    int n = A.colsize();
    int kd = A.nlo()==0 ? A.nhi() : A.nlo();
    int nrhs = B.rowsize();
    int aoffset = 
      (A.nlo()==0 && A.iscm()) ? A.nhi() :
      (A.nhi()==0 && A.isrm()) ? A.nlo() : 0;
    int lda = A.diagstep();
    int ldb = B.stepj();
    if (A.iscm() && A.isconj()) {
      B.ConjugateSelf();
      LAPNAME(ztbtrs) (LAPCM A.nlo()==0?LAPCH_UP:LAPCH_LO,
	  LAPCH_NT, dt==UnitDiag?LAPCH_U:LAPCH_NU,
	  LAPV(n),LAPV(kd),LAPV(nrhs), LAPP(A.cptr()-aoffset),LAPV(lda),
	  LAPP(B.ptr()),LAPV(ldb) LAPINFO LAP1 LAP1 LAP1);
      B.ConjugateSelf();
    } else {
      LAPNAME(ztbtrs) (LAPCM (A.iscm()?A.nlo():A.nhi())==0?LAPCH_UP:LAPCH_LO,
	  A.iscm()?LAPCH_NT:A.isconj()?LAPCH_CT:LAPCH_T, 
	  dt==UnitDiag?LAPCH_U:LAPCH_NU,
	  LAPV(n),LAPV(kd),LAPV(nrhs), LAPP(A.cptr()-aoffset),LAPV(lda),
	  LAPP(B.ptr()),LAPV(ldb) LAPINFO LAP1 LAP1 LAP1);
    }
    LAP_Results("ztbtrs");
  }
#ifndef NOFLOAT
  template <> inline void LapBandTriLDivEq(
      const GenBandMatrix<float>& A, const MatrixView<float>& B, DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(B.colsize() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.colsize() > 0);
    TMVAssert(B.rowsize() > 1);
    TMVAssert(A.nhi() >= 0);
    TMVAssert(A.nlo() >= 0);
    TMVAssert(A.nhi() < int(A.colsize()));
    TMVAssert(A.nlo() < int(A.colsize()));
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(B.iscm());
    int n = A.colsize();
    int kd = A.nlo()==0 ? A.nhi() : A.nlo();
    int nrhs = B.rowsize();
    int aoffset = 
      (A.nlo()==0 && A.iscm()) ? A.nhi() :
      (A.nhi()==0 && A.isrm()) ? A.nlo() : 0;
    int lda = A.diagstep();
    int ldb = B.stepj();
    LAPNAME(stbtrs) (LAPCM (A.iscm()?A.nlo():A.nhi())==0?LAPCH_UP:LAPCH_LO,
	A.iscm()?LAPCH_NT:LAPCH_T, dt==UnitDiag?LAPCH_U:LAPCH_NU,
	LAPV(n),LAPV(kd),LAPV(nrhs), LAPP(A.cptr()-aoffset),LAPV(lda),
	LAPP(B.ptr()),LAPV(ldb) LAPINFO LAP1 LAP1 LAP1);
    LAP_Results("stbtrs");
  }
  template <> inline void LapBandTriLDivEq(
      const GenBandMatrix<complex<float> >& A,
      const MatrixView<complex<float> >& B, DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(B.colsize() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.colsize() > 0);
    TMVAssert(B.rowsize() > 1);
    TMVAssert(A.nhi() >= 0);
    TMVAssert(A.nlo() >= 0);
    TMVAssert(A.nhi() < int(A.colsize()));
    TMVAssert(A.nlo() < int(A.colsize()));
    TMVAssert(B.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(B.iscm());
    int n = A.colsize();
    int kd = A.nlo()==0 ? A.nhi() : A.nlo();
    int nrhs = B.rowsize();
    int aoffset = 
      (A.nlo()==0 && A.iscm()) ? A.nhi() :
      (A.nhi()==0 && A.isrm()) ? A.nlo() : 0;
    int lda = A.diagstep();
    int ldb = B.stepj();
    if (A.iscm() && A.isconj()) {
      B.ConjugateSelf();
      LAPNAME(ctbtrs) (LAPCM A.nlo()==0?LAPCH_UP:LAPCH_LO,
	  LAPCH_NT, dt==UnitDiag?LAPCH_U:LAPCH_NU,
	  LAPV(n),LAPV(kd),LAPV(nrhs), LAPP(A.cptr()-aoffset),LAPV(lda),
	  LAPP(B.ptr()),LAPV(ldb) LAPINFO LAP1 LAP1 LAP1);
      B.ConjugateSelf();
    } else {
      LAPNAME(ctbtrs) (LAPCM (A.iscm()?A.nlo():A.nhi())==0?LAPCH_UP:LAPCH_LO,
	  A.iscm()?LAPCH_NT:A.isconj()?LAPCH_CT:LAPCH_T, 
	  dt==UnitDiag?LAPCH_U:LAPCH_NU,
	  LAPV(n),LAPV(kd),LAPV(nrhs), LAPP(A.cptr()-aoffset),LAPV(lda),
	  LAPP(B.ptr()),LAPV(ldb) LAPINFO LAP1 LAP1 LAP1);
    }
    LAP_Results("ctbtrs");
  }
#endif
#endif // LAP

  template <class T, class Ta> void BandTriLDivEq(
      const GenBandMatrix<Ta>& A, const MatrixView<T>& B, DiagType dt)
  {
#ifdef XDEBUG
    Matrix<T> B0 = B;
#endif
    TMVAssert(A.IsSquare());
    TMVAssert(B.colsize() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.colsize() > 0);
    if (B.rowsize() == 0) return;
    TMVAssert(A.nhi() >= 0);
    TMVAssert(A.nlo() >= 0);
    TMVAssert(A.nhi() < int(A.colsize()));
    TMVAssert(A.nlo() < int(A.colsize()));

    if (B.rowsize() == 1) BandTriLDivEq(A,B.col(0),dt);
    else if (B.isconj()) 
      BandTriLDivEq(A.Conjugate(),B.Conjugate(),dt);
    else
#ifdef LAP
      if (!(B.iscm() && B.stepj()>0)) {
	Matrix<T,ColMajor> BB=B;
	BandTriLDivEq(A,BB.View(),dt);
	B = BB;
      } else if ( !((A.iscm() && A.stepj()>0) || (A.isrm() && A.stepi()>0))) {
	BandMatrix<T,ColMajor> AA = A;
	LapBandTriLDivEq(AA,B,dt);
      } else 
	LapBandTriLDivEq(A,B,dt);
#else
	if (dt == UnitDiag) NonLapBandTriLDivEq<true>(A,B);
	else NonLapBandTriLDivEq<false>(A,B);
#endif

#ifdef XDEBUG
    Matrix<T> BB = A*B;
    if (Norm(BB-B0) > 0.001*Norm(B0)) {
      cerr<<"BandTriLDivEq Matrix:\n";
      cerr<<"A = "<<Type(A)<<A<<endl;
      cerr<<"B = "<<Type(B)<<B0<<endl;
      cerr<<"--> B = "<<B<<endl;
      cerr<<"A*B = "<<BB<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_BandLUDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


