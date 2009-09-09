
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

namespace tmv {

#define RecursiveLU

#ifdef TMV_BLOCKSIZE
  const size_t LU_BLOCKSIZE = TMV_BLOCKSIZE;
  const size_t LU_BLOCKSIZE2 = TMV_BLOCKSIZE/2;
#else
  const size_t LU_BLOCKSIZE = 64;
  const size_t LU_BLOCKSIZE2 = 2;
#endif

  //
  // LDivEq
  //

  template <class T, class T1> void NonLapLU_LDivEq(
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
#ifdef CLAP
    clapack_dgetrs(CblasColMajor,CblasNoTrans,LUx.colsize(),m.rowsize(),
	LUx.cptr(),LUx.stepj(),LAP_IPiv(n),m.ptr(),m.stepj());
#else
    char c = 'N';
    int nrhs = m.rowsize();
    int lda = LUx.stepj();
    int ldb = m.stepj();
    int info;
    dgetrs(&c,&n,&nrhs,const_cast<double*>(LUx.cptr()),&lda,
	LAP_IPiv(n),m.ptr(),&ldb,&info);
    if (info < 0) tmv_error("dgetrs returned info < 0");
#endif
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
#ifdef CLAP
    clapack_zgetrs(CblasColMajor,CblasNoTrans,LUx.colsize(),m.rowsize(),
	LUx.cptr(),LUx.stepj(),LAP_IPiv(n),m.ptr(),m.stepj());
#else
    char c = 'N';
    int nrhs = m.rowsize();
    int lda = LUx.stepj();
    int ldb = m.stepj();
    int info;
    zgetrs(&c,&n,&nrhs,LAP_Complex(LUx.cptr()),&lda,
	LAP_IPiv(n),LAP_Complex(m.ptr()),&ldb,&info);
    if (info < 0) tmv_error("zgetrs returned info < 0");
#endif
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
#ifdef CLAP
    clapack_sgetrs(CblasColMajor,CblasNoTrans,LUx.colsize(),m.rowsize(),
	LUx.cptr(),LUx.stepj(),LAP_IPiv(n),m.ptr(),m.stepj());
#else
    char c = 'N';
    int nrhs = m.rowsize();
    int lda = LUx.stepj();
    int ldb = m.stepj();
    int info;
    sgetrs(&c,&n,&nrhs,const_cast<float*>(LUx.cptr()),&lda,
	LAP_IPiv(n),m.ptr(),&ldb,&info);
    if (info < 0) tmv_error("sgetrs returned info < 0");
#endif
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
#ifdef CLAP
    clapack_cgetrs(CblasColMajor,CblasNoTrans,LUx.colsize(),m.rowsize(),
	LUx.cptr(),LUx.stepj(),LAP_IPiv(n),m.ptr(),m.stepj());
#else
    char c = 'N';
    int nrhs = m.rowsize();
    int lda = LUx.stepj();
    int ldb = m.stepj();
    int info;
    cgetrs(&c,&n,&nrhs,LAP_Complex(LUx.cptr()),&lda,
	LAP_IPiv(n),LAP_Complex(m.ptr()),&ldb,&info);
    if (info < 0) tmv_error("cgetrs returned info < 0");
#endif
  }
#endif // MKL
#endif // FLOAT
#endif // ALAP

  template <class T, class T1> inline void LU_LDivEq(
      const GenMatrix<T1>& LUx, const size_t* P, 
      const MatrixView<T>& m)
  {
    TMVAssert(m.colsize() == LUx.rowsize()); 
    TMVAssert(LUx.rowsize() == LUx.colsize());

    m.PermuteRows(P); 

#ifdef ALAP
    if (m.iscm() && !m.isconj() && LUx.iscm() && !LUx.isconj())
      LapLU_LDivEq(LUx,m);
    else 
#endif
      NonLapLU_LDivEq(LUx,m); 
  }

  //
  // RDivEq Matrix
  //

  template <class T, class T1> void NonLapLU_RDivEq(
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
#ifdef CLAP
    clapack_dgetrs(CblasColMajor,CblasTrans,LUx.colsize(),m.colsize(),
	LUx.cptr(),LUx.stepj(),LAP_IPiv(n),m.ptr(),m.stepi());
#else
    char c = 'T';
    int nrhs = m.colsize();
    int lda = LUx.stepj();
    int ldb = m.stepi();
    int info;
    dgetrs(&c,&n,&nrhs,const_cast<double*>(LUx.cptr()),&lda,
	LAP_IPiv(n),m.ptr(),&ldb,&info);
    if (info < 0) tmv_error("dgetrs returned info < 0");
#endif
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
#ifdef CLAP
    clapack_zgetrs(CblasColMajor,CblasTrans,LUx.colsize(),m.colsize(),
	LUx.cptr(),LUx.stepj(),LAP_IPiv(n),m.ptr(),m.stepi());
#else
    char c = (m.isconj() ? 'C' : 'T');
    int nrhs = m.colsize();
    int lda = LUx.stepj();
    int ldb = m.stepi();
    int info;
    zgetrs(&c,&n,&nrhs,LAP_Complex(LUx.cptr()),&lda,
	LAP_IPiv(n),LAP_Complex(m.ptr()),&ldb,&info);
    if (info < 0) tmv_error("zgetrs returned info < 0");
#endif
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
#ifdef CLAP
    clapack_sgetrs(CblasColMajor,CblasTrans,LUx.colsize(),m.colsize(),
	LUx.cptr(),LUx.stepj(),LAP_IPiv(n),m.ptr(),m.stepi());
#else
    char c = 'T';
    int nrhs = m.colsize();
    int lda = LUx.stepj();
    int ldb = m.stepi();
    int info;
    sgetrs(&c,&n,&nrhs,const_cast<float*>(LUx.cptr()),&lda,
	LAP_IPiv(n),m.ptr(),&ldb,&info);
    if (info < 0) tmv_error("sgetrs returned info < 0");
#endif
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
#ifdef CLAP
    clapack_cgetrs(CblasColMajor,CblasTrans,LUx.colsize(),m.colsize(),
	LUx.cptr(),LUx.stepj(),LAP_IPiv(n),m.ptr(),m.stepi());
#else
    char c = (m.isconj() ? 'C' : 'T');
    int nrhs = m.colsize();
    int lda = LUx.stepj();
    int ldb = m.stepi();
    int info;
    cgetrs(&c,&n,&nrhs,LAP_Complex(LUx.cptr()),&lda,
	LAP_IPiv(n),LAP_Complex(m.ptr()),&ldb,&info);
    if (info < 0) tmv_error("cgetrs returned info < 0");
#endif
  }
#endif // MKL
#endif // FLOAT
#endif // ALAP

  template <class T, class T1> inline void LU_RDivEq(
      const GenMatrix<T1>& LUx, const size_t* P, 
      const MatrixView<T>& m)
    // Solve x P L U = m:
  {
    TMVAssert(m.rowsize() == LUx.rowsize()); 
    TMVAssert(LUx.rowsize() == LUx.colsize());

#ifdef ALAP
    if (LUx.iscm() && !LUx.isconj() && m.isrm())
      LapLU_RDivEq(LUx,m);
    else 
#endif
      NonLapLU_RDivEq(LUx,m); 

    m.ReversePermuteCols(P); 
  }

#define InstFile "TMV_LUDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


