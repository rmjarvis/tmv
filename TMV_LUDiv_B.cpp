
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

  //
  // LDivEq
  //

  template <class T, class T1> inline void NonLapLU_LDivEq(
      const GenMatrix<T1>& LUx, const MatrixView<T>& m)
  {
    // Solve L U x = m:
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.colsize());
    m /= LowerTriMatrixViewOf(LUx,UnitDiag);
    m /= UpperTriMatrixViewOf(LUx,NonUnitDiag);
  }


#ifdef ALAP
  // ALAP, not LAP, since ATLAS has these routines
  template <class T, class T1> inline void LapLU_LDivEq(
      const GenMatrix<T1>& LUx, const MatrixView<T>& m)
  { NonLapLU_LDivEq(LUx,m); }
  template <> inline void LapLU_LDivEq(
      const GenMatrix<double>& LUx, const MatrixView<double>& m)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.colsize());
    TMVAssert(LUx.iscm());
    TMVAssert(m.iscm());
    TMVAssert(LUx.ct()==NonConj);
    TMVAssert(m.ct()==NonConj);

    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int lda = LUx.stepj();
    int ldb = m.stepj();
    LAPNAME(dgetrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
	LAPP(LUx.cptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
    LAP_Results("dgetrs");
  }
  template <> inline void LapLU_LDivEq(
      const GenMatrix<complex<double> >& LUx,
      const MatrixView<complex<double> >& m)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.colsize());
    TMVAssert(LUx.iscm());
    TMVAssert(m.iscm());
    TMVAssert(LUx.ct()==NonConj);
    TMVAssert(m.ct()==NonConj);

    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int lda = LUx.stepj();
    int ldb = m.stepj();
    LAPNAME(zgetrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
	LAPP(LUx.cptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
    LAP_Results("zgetrs");
  }
#ifndef NOFLOAT
#ifndef MKL
  template <> inline void LapLU_LDivEq(
      const GenMatrix<float>& LUx, const MatrixView<float>& m)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.colsize());
    TMVAssert(LUx.iscm());
    TMVAssert(m.iscm());
    TMVAssert(LUx.ct()==NonConj);
    TMVAssert(m.ct()==NonConj);

    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int lda = LUx.stepj();
    int ldb = m.stepj();
    LAPNAME(sgetrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
	LAPP(LUx.cptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
    LAP_Results("sgetrs");
  }
  template <> inline void LapLU_LDivEq(
      const GenMatrix<complex<float> >& LUx,
      const MatrixView<complex<float> >& m)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.colsize());
    TMVAssert(LUx.iscm());
    TMVAssert(m.iscm());
    TMVAssert(LUx.ct()==NonConj);
    TMVAssert(m.ct()==NonConj);

    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int lda = LUx.stepj();
    int ldb = m.stepj();
    LAPNAME(cgetrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
	LAPP(LUx.cptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
    LAP_Results("cgetrs");
  }
#endif // MKL
#endif // FLOAT
#endif // ALAP

  template <class T, class T1> void LU_LDivEq(
      const GenMatrix<T1>& LUx, const size_t* P, 
      const MatrixView<T>& m)
  {
    TMVAssert(m.colsize() == LUx.rowsize()); 
    TMVAssert(LUx.rowsize() == LUx.colsize());

#ifdef XDEBUG
    Matrix<T> m0 = m;
    Matrix<T1> L = LowerTriMatrixViewOf(LUx,UnitDiag);
    Matrix<T1> U = UpperTriMatrixViewOf(LUx,NonUnitDiag);
    Matrix<T1> LU = L*U;
    LU.ReversePermuteRows(P);
    //cerr<<"Start LU_LDivEq\n";
    //cerr<<"LU = "<<LU<<endl;
    //cerr<<"m0 = "<<m0<<endl;
#endif

    m.PermuteRows(P); 

#ifdef ALAP
    if (m.iscm() && !m.isconj() && LUx.iscm() && !LUx.isconj())
      LapLU_LDivEq(LUx,m);
    else 
#endif
      NonLapLU_LDivEq(LUx,m); 

#ifdef XDEBUG
    //cerr<<"m-> "<<m<<endl;
    Matrix<T> mm = LU*m;
    //cerr<<"mm = "<<mm<<endl;
    if (Norm(mm-m0) > 0.001*Norm(L)*Norm(U)*Norm(m0)) {
      cerr<<"LU_LDivEq: m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"L = "<<L<<endl;
      cerr<<"U = "<<U<<endl;
      cerr<<"LU = "<<LU<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"LU*m = "<<mm<<endl;
      abort();
    }
#endif
  }

  //
  // RDivEq Matrix
  //

  template <class T, class T1> inline void NonLapLU_RDivEq(
      const GenMatrix<T1>& LUx, const MatrixView<T>& m)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.rowsize());
    // m = m (LU)^-1 
    //   = m U^-1 L^-1
    m %= UpperTriMatrixViewOf(LUx,NonUnitDiag);
    m %= LowerTriMatrixViewOf(LUx,UnitDiag);

  }

#ifdef ALAP
  template <class T, class T1> inline void LapLU_RDivEq(
      const GenMatrix<T1>& LUx, const MatrixView<T>& m)
  { NonLapLU_RDivEq(LUx,m); }
  template <> inline void LapLU_RDivEq(
      const GenMatrix<double>& LUx, const MatrixView<double>& m)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.rowsize());
    TMVAssert(LUx.iscm());
    TMVAssert(m.isrm());
    TMVAssert(m.ct()==NonConj);
    TMVAssert(LUx.ct()==NonConj);

    int n = LUx.colsize();
    int nrhs = m.colsize();
    int lda = LUx.stepj();
    int ldb = m.stepi();
    LAPNAME(dgetrs) (LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPP(LUx.cptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
    LAP_Results("dgetrs");
  }
  template <> inline void LapLU_RDivEq(
      const GenMatrix<complex<double> >& LUx,
      const MatrixView<complex<double> >& m)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.rowsize());
    TMVAssert(LUx.iscm());
    TMVAssert(m.isrm());
    TMVAssert(LUx.ct()==NonConj);

    int n = LUx.colsize();
    int nrhs = m.colsize();
    int lda = LUx.stepj();
    int ldb = m.stepi();
    LAPNAME(zgetrs) (LAPCM m.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPP(LUx.cptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
    LAP_Results("zgetrs");
  }
#ifndef NOFLOAT
#ifndef MKL
  template <> inline void LapLU_RDivEq(
      const GenMatrix<float>& LUx, const MatrixView<float>& m)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.rowsize());
    TMVAssert(LUx.iscm());
    TMVAssert(m.isrm());
    TMVAssert(m.ct()==NonConj);
    TMVAssert(LUx.ct()==NonConj);

    int n = LUx.colsize();
    int nrhs = m.colsize();
    int lda = LUx.stepj();
    int ldb = m.stepi();
    LAPNAME(sgetrs) (LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPP(LUx.cptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
    LAP_Results("sgetrs");
  }
  template <> inline void LapLU_RDivEq(
      const GenMatrix<complex<float> >& LUx,
      const MatrixView<complex<float> >& m)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.rowsize());
    TMVAssert(LUx.iscm());
    TMVAssert(m.isrm());
    TMVAssert(LUx.ct()==NonConj);

    int n = LUx.colsize();
    int nrhs = m.colsize();
    int lda = LUx.stepj();
    int ldb = m.stepi();
    LAPNAME(cgetrs) (LAPCM m.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPP(LUx.cptr()),LAPV(lda),
	LAPP(LAP_IPiv(n)),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
    LAP_Results("cgetrs");
  }
#endif // MKL
#endif // FLOAT
#endif // ALAP

  template <class T, class T1> void LU_RDivEq(
      const GenMatrix<T1>& LUx, const size_t* P, 
      const MatrixView<T>& m)
    // Solve x P L U = m:
  {
#ifdef XDEBUG
    Matrix<T> m0 = m;
    Matrix<T1> L = LowerTriMatrixViewOf(LUx,UnitDiag);
    Matrix<T1> U = UpperTriMatrixViewOf(LUx,NonUnitDiag);
    Matrix<T1> LU = L*U;
    LU.ReversePermuteRows(P);
    //cerr<<"Start LU_RDivEq\n";
    //cerr<<"LU = "<<LU<<endl;
    //cerr<<"m0 = "<<m0<<endl;
#endif

    TMVAssert(m.rowsize() == LUx.rowsize()); 
    TMVAssert(LUx.rowsize() == LUx.colsize());

#ifdef ALAP
    if (LUx.iscm() && !LUx.isconj() && m.isrm())
      LapLU_RDivEq(LUx,m);
    else 
#endif
      NonLapLU_RDivEq(LUx,m); 

    m.ReversePermuteCols(P); 

#ifdef XDEBUG
    //cerr<<"m-> "<<m<<endl;
    Matrix<T> mm = m*LU;
    //cerr<<"mm = "<<mm<<endl;
    if (Norm(mm-m0) > 0.001*Norm(L)*Norm(U)*Norm(m0)) {
      cerr<<"LU_RDivEq: m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"L = "<<L<<endl;
      cerr<<"U = "<<U<<endl;
      cerr<<"LU = "<<LU<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"m*LU = "<<mm<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_LUDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


