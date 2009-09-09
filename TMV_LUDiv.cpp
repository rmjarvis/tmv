
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

namespace tmv {

  //
  // LDivEq
  //

  template <class T1, class T2> void NonLapLU_LDivEq(
      const GenMatrix<T1>& LUx, const MatrixView<T2>& m)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.colsize());
    // m = (LU)^-1 m
    //   = U^-1 L^-1 m
    m /= LowerTriMatrixViewOf(LUx,UnitDiag);
    m /= UpperTriMatrixViewOf(LUx,NonUnitDiag);
  }

#ifdef LAP
  template <class T1, class T2> inline void LapLU_LDivEq(
      const GenMatrix<T1>& LUx, const MatrixView<T2>& m)
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
    char c = 'N';
    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int lda = LUx.stepj();
    int ldb = m.stepj();
    int info;
    dgetrs(&c,&n,&nrhs,const_cast<double*>(LUx.cptr()),&lda,
	LAP_IPiv(n),m.ptr(),&ldb,&info);
    if (info < 0) tmv_error("dgetrs returned info < 0");
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
    char c = 'N';
    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int lda = LUx.stepj();
    int ldb = m.stepj();
    int info;
    zgetrs(&c,&n,&nrhs,LAP_Complex(LUx.cptr()),&lda,
	LAP_IPiv(n),LAP_Complex(m.ptr()),&ldb,&info);
    if (info < 0) tmv_error("zgetrs returned info < 0");
  }
#ifndef NOFLOAT
  template <> inline void LapLU_LDivEq(
      const GenMatrix<float>& LUx, const MatrixView<float>& m)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.colsize());
    TMVAssert(LUx.iscm());
    TMVAssert(m.iscm());
    TMVAssert(LUx.ct()==NonConj);
    TMVAssert(m.ct()==NonConj);
    char c = 'N';
    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int lda = LUx.stepj();
    int ldb = m.stepj();
    int info;
    sgetrs(&c,&n,&nrhs,const_cast<float*>(LUx.cptr()),&lda,
	LAP_IPiv(n),m.ptr(),&ldb,&info);
    if (info < 0) tmv_error("sgetrs returned info < 0");
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
    char c = 'N';
    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int lda = LUx.stepj();
    int ldb = m.stepj();
    int info;
    cgetrs(&c,&n,&nrhs,LAP_Complex(LUx.cptr()),&lda,
	LAP_IPiv(n),LAP_Complex(m.ptr()),&ldb,&info);
    if (info < 0) tmv_error("cgetrs returned info < 0");
  }
#endif
#endif // LAP

  template <class T1, class T2> inline void LU_LDivEq(
      const GenMatrix<T1>& LUx, const Permutation& P, 
      const MatrixView<T2>& m)
    // Solve P L U x = m:
    // y = m / P
    // Solve L z = y
    // Solve U x = z
  {
    TMVAssert(m.colsize() == LUx.rowsize()); 
    TMVAssert(m.colsize() == P.size()); 
    TMVAssert(LUx.rowsize() == LUx.colsize());
    TMVAssert(LUx.iscm());
    TMVAssert(LUx.ct() == NonConj);

    m /= P; 

#ifdef LAP
    if (m.iscm() && !m.isconj())
      LapLU_LDivEq(LUx,m);
    else 
#endif
      NonLapLU_LDivEq(LUx,m); 
  }

  //
  // RDivEq Matrix
  //

  template <class T1, class T2> void NonLapLU_RDivEq(
      const GenMatrix<T1>& LUx, const MatrixView<T2>& m)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.rowsize());
    // m = m (LU)^-1 
    //   = m U^-1 L^-1
    m %= UpperTriMatrixViewOf(LUx,NonUnitDiag);
    m %= LowerTriMatrixViewOf(LUx,UnitDiag);
  }

#ifdef LAP
  template <class T1, class T2> 
    inline void LapLU_RDivEq(
	const GenMatrix<T1>& LUx, const MatrixView<T2>& m)
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
    char c = 'T';
    int n = LUx.colsize();
    int nrhs = m.colsize();
    int lda = LUx.stepj();
    int ldb = m.stepi();
    int info;
    dgetrs(&c,&n,&nrhs,const_cast<double*>(LUx.cptr()),&lda,
	LAP_IPiv(n),m.ptr(),&ldb,&info);
    if (info < 0) tmv_error("dgetrs returned info < 0");
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
    char c = (m.isconj() ? 'C' : 'T');
    int n = LUx.colsize();
    int nrhs = m.colsize();
    int lda = LUx.stepj();
    int ldb = m.stepi();
    int info;
    zgetrs(&c,&n,&nrhs,LAP_Complex(LUx.cptr()),&lda,
	LAP_IPiv(n),LAP_Complex(m.ptr()),&ldb,&info);
    if (info < 0) tmv_error("zgetrs returned info < 0");
  }
#ifndef NOFLOAT
  template <> inline void LapLU_RDivEq(
      const GenMatrix<float>& LUx, const MatrixView<float>& m)
  {
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.rowsize());
    TMVAssert(LUx.iscm());
    TMVAssert(m.isrm());
    TMVAssert(m.ct()==NonConj);
    TMVAssert(LUx.ct()==NonConj);
    char c = 'T';
    int n = LUx.colsize();
    int nrhs = m.colsize();
    int lda = LUx.stepj();
    int ldb = m.stepi();
    int info;
    sgetrs(&c,&n,&nrhs,const_cast<float*>(LUx.cptr()),&lda,
	LAP_IPiv(n),m.ptr(),&ldb,&info);
    if (info < 0) tmv_error("sgetrs returned info < 0");
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
    char c = (m.isconj() ? 'C' : 'T');
    int n = LUx.colsize();
    int nrhs = m.colsize();
    int lda = LUx.stepj();
    int ldb = m.stepi();
    int info;
    cgetrs(&c,&n,&nrhs,LAP_Complex(LUx.cptr()),&lda,
	LAP_IPiv(n),LAP_Complex(m.ptr()),&ldb,&info);
    if (info < 0) tmv_error("cgetrs returned info < 0");
  }
#endif
#endif // LAP

  template <class T1, class T2> inline void LU_RDivEq(
      const GenMatrix<T1>& LUx, const Permutation& P, 
      const MatrixView<T2>& m)
    // Solve x P L U = m:
    // Solve y U = m
    // Solve z L = y
    // x = z % P
  {
    TMVAssert(m.rowsize() == LUx.rowsize()); 
    TMVAssert(LUx.rowsize() == LUx.colsize());
    TMVAssert(LUx.iscm());
    TMVAssert(LUx.stepi() == 1);
    TMVAssert(LUx.ct() == NonConj);
    TMVAssert(m.rowsize() == P.size()); 

#ifdef LAP
    if (m.isrm())
      LapLU_RDivEq(LUx,m);
    else 
#endif
      NonLapLU_RDivEq(LUx,m); 

    m %= P; 
  }

  template <class T> T LUDiv<T>::Det() const
  {
    if (!donedet) {
      det = P.Det();
      det *= DiagMatrixViewOf(LUx.diag()).Det();
      donedet = true;
    }
    return det;
  }

  template <class T> Matrix<T,ColMajor> LUDiv<T>::Inverse() const
  {
    Matrix<T,ColMajor> temp = Eye<T,ColMajor>(LUx.colsize());
    LDivEq(temp.QuickView());
    return temp;
  }

  template <class T> Matrix<T,ColMajor> LUDiv<T>::InverseATA() const
  {
    Matrix<T,ColMajor> temp = Eye<T,ColMajor>(LUx.colsize());
    LDivEq(temp.QuickView());
    return temp*temp.QuickAdjoint();
  }

  template <class T> Matrix<T,ColMajor> LUDiv<T>::LU_GetL() const 
  { return Matrix<T,ColMajor>(LowerTriMatrixViewOf(LUx,UnitDiag)); }

  template <class T> Matrix<T,ColMajor> LUDiv<T>::LU_GetU() const 
  { return Matrix<T,ColMajor>(UpperTriMatrixViewOf(LUx,NonUnitDiag)); }

  //
  // Decompose
  //

  template <class T> void NonLapLU_Decompose(
      const MatrixView<T>& A, Permutation& P)
  {
    // LU Decompostion with partial pivoting and row scaling.  
    //
    // We want to decompose the matrix (input as A) into P * L * U
    // where P is a permutation, L is a lower triangle matrix with 1's for 
    // the diagonal, and U is an upper triangle matrix.  
    //
    // We do this calculation one column at a time.  There are other versions
    // of this algorithm which use Rank 1 updates.  This version uses mostly
    // martix-vector products and forward substitutions.  The different
    // algorithms all require the same number of operation, but this
    // one requires fewer vector touches which often means it will be 
    // slightly faster.
    //
    // After doing j cols of the calculation, we have calculated
    // U(0:j,0:j) and L(0:M,0:j) 
    // (where my a:b notation does not include the index b)
    //
    // The equation A = LU gives for the j col:
    //
    // A(0:M,j) = L(0:M,0:R) U(0:R,j)
    // where R = min(M,N)
    //
    // which breaks up into:
    //
    // (1) A(0:j,j) = L(0:j,0:R) U(0:R,j)
    // (2) A(j:M,j) = L(j:M,0:R) U(0:R,j)
    //
    // The first of these (1) simplifies to:
    // 
    // (1*) A(0:j,j) = L(0:j,0:j) U(0:j,j)
    //
    // since L is lower triangular, so L(0:j,j:R) = 0.
    // L(0:j,0:j) is already known, so this equation can be solved for
    // U(0:j,j) by forward substitution.
    // 
    // The second equation (2) simplifies to:
    //
    //      A(j:M,j) = L(j:M,0:j+1) U(0:j+1,j)
    // (2*)          = L(j:M,0:j) U(0:j,j) + L(j:M,j) U(j,j)
    //
    // since U is upper triangular so U(j+1:R,j) = 0.
    // Since we now know U(0:j,j) from (1*) above, this equation can
    // be solved for the product L(j:M,j) U(j,j)
    // 
    // This means we have some leeway on the values for L(j,j) and U(j,j),
    // as only their product is specified.
    //
    // If we take U to have unit diagonal, then L(j,j) is set here along
    // with the rest of the L(j:M,j) column.  However, this will mean that the
    // forward substutions in the (1*) steps will require divisions by the 
    // non-unit-diagonal elements of L.  It is faster to take L to have
    // unit-diagonal elements, and do the division by U(j,j) here, since then
    // we can calculate 1/U(j,j) and multiply.  So 1 division and M-j
    // which is generally faster than M-j divisions.
    //
    // However, another potential problem is that U(j,j) could be 0, or close
    // to 0.  This would lead to either an error or inaccurate results.
    // Thus we add a step in the middle of the (2*) calculation:
    //
    // Define v(j:M) = A(j:M,j) - L(j:M,0:j) U(0:j,j)
    // 
    // We search v for the element with the largest absolute value and apply
    // a permutation to swap it into the j spot.  This element then becomes 
    // U(j,j), which is then the divisor for the rest of the vector.  This 
    // will minimize the possibility of roundoff errors due to small U(j,j)'s.
    // 
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == P.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    const size_t M = A.colsize();
    const size_t N = A.rowsize();
    const size_t R = min(M,N);

    // We implicitly rescale each row by the largest (absolute) value.
    // If a particular row (and the vector on the right hand side)
    // is multiplied by a million, the essential equations would be
    // the same, but this row would almost certainly be the first 
    // pivot element.  Rescaling each row by the largest element
    // for the purposes of determining the pivot avoids this type of problem.
    Vector<RealType(T)> scale(M);
    for (size_t i=0; i<M; ++i) scale(i) = MaxAbsElement(A.row(i));

    for (size_t j=0; j<R; ++j)
    {
      if (j > 0) {
	// Solve for U(0:j,j))
	A.col(j,0,j) /= LowerTriMatrixViewOf(A.SubMatrix(0,j,0,j),UnitDiag);

	// Solve for v = L(j:M,j) U(j,j)
	A.col(j,j,M) -= A.SubMatrix(j,M,0,j) * A.col(j,0,j);
      }

      // Find the pivot element
      size_t ip=j;
      RealType(T) pivot = abs(A(j,j));
      if (scale(j)) pivot /= scale(j);
      for(size_t i=j+1; i<M; ++i) {
	RealType(T) temp = abs(A(i,j));
	if (scale(i)) temp /= scale(i);
	if (temp > pivot) { pivot = temp; ip = i; }
      }

      // Swap the pivot row with j if necessary
      if (ip != j) {
	A.SwapRows(ip,j);  // This does both Lkb and A'
	P.SwapCols(ip,j);
	scale.Swap(ip,j);
      }

      // Solve for L(j+1:M,j)
      // If Ujj is 0, then all of the L's are 0.
      // ie. Ujj Lij = 0 for all i>j
      // Any value for Lij is valid, so leave them 0.
      T Ujj = A(j,j);
      if (Ujj != T(0)) A.col(j,j+1,M) /= Ujj;
    }
  }

#ifdef LAP
  template <class T> inline void LapLU_Decompose(
      const MatrixView<T>& A, Permutation& P)
  { NonLapLU_Decompose(A,P); }
  template <> inline void LapLU_Decompose(
      const MatrixView<double>& A, Permutation& P)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == P.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int lap_p[P.size()];
    int m = A.colsize();
    int n = A.rowsize();
    int lda = A.stepj();
    int info;
    dgetrf(&m,&n,A.ptr(),&lda,lap_p,&info);
    if (info < 0) tmv_error("dgetrf returned info < 0");
    for(size_t i=0;i<P.size();++i) {
      if (lap_p[i]-1 != int(i)) P.SwapCols(i,lap_p[i]-1);
    }
  }
  template <> inline void LapLU_Decompose(
      const MatrixView<complex<double> >& A, Permutation& P)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == P.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int lap_p[P.size()];
    int m = A.colsize();
    int n = A.rowsize();
    int lda = A.stepj();
    int info;
    zgetrf(&m,&n,LAP_Complex(A.ptr()),&lda,lap_p,&info);
    if (info < 0) tmv_error("zgetrf returned info < 0");
    for(size_t i=0;i<P.size();++i) 
      if (lap_p[i]-1 != int(i)) P.SwapCols(i,lap_p[i]-1);
  }
#ifndef NOFLOAT
  template <> inline void LapLU_Decompose(
      const MatrixView<float>& A, Permutation& P)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == P.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int lap_p[P.size()];
    int m = A.colsize();
    int n = A.rowsize();
    int lda = A.stepj();
    int info;
    sgetrf(&m,&n,A.ptr(),&lda,lap_p,&info);
    if (info < 0) tmv_error("sgetrf returned info < 0");
    for(size_t i=0;i<P.size();++i) {
      if (lap_p[i]-1 != int(i)) P.SwapCols(i,lap_p[i]-1);
    }
  }
  template <> inline void LapLU_Decompose(
      const MatrixView<complex<float> >& A, Permutation& P)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == P.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int lap_p[P.size()];
    int m = A.colsize();
    int n = A.rowsize();
    int lda = A.stepj();
    int info;
    cgetrf(&m,&n,LAP_Complex(A.ptr()),&lda,lap_p,&info);
    if (info < 0) tmv_error("cgetrf returned info < 0");
    for(size_t i=0;i<P.size();++i) 
      if (lap_p[i]-1 != int(i)) P.SwapCols(i,lap_p[i]-1);
  }
#endif // NOFLOAT
#endif // LAP
  template <class T> inline void LU_Decompose(
      const MatrixView<T>& A, Permutation& P)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == P.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    if (A.colsize() > 0 && A.rowsize() > 0) {
#ifdef LAP
      LapLU_Decompose(A,P);
#else
      NonLapLU_Decompose(A,P);
#endif
    }
  }

#define InstFile "TMV_LUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


