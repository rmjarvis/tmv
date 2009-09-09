
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"
#include "TMV_VectorArith_Inline.h"

//#define XDEBUG

namespace tmv {

  bool RecursiveLU = true;

#ifdef TMV_BLOCKSIZE
  const size_t LU_BLOCKSIZE = TMV_BLOCKSIZE;
#else
  const size_t LU_BLOCKSIZE = 64;
#endif

#define APTR inplace ? A.NonConst().ptr() : new T[A.colsize()*A.rowsize()]
#define LUX istrans ? \
  inplace ? A.NonConst().Transpose() : \
  MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor) : \
  inplace ? A.NonConst().View() : \
  MatrixViewOf(Aptr,A.colsize(),A.rowsize(),ColMajor)

  template <class T> LUDiv<T>::LUDiv(const GenMatrix<T>& A, bool _inplace) :
    istrans(A.isrm()), inplace(_inplace),
    Aptr(APTR), LUx(LUX), P(new size_t[A.colsize()]),
    det(T(1)), donedet(false)
  {
    TMVAssert(A.IsSquare());
    if (istrans) {
      if (inplace) { TMVAssert(A.Transpose() == LUx); }
      else LUx = A.Transpose();
    }
    else {
      if (inplace) { TMVAssert(A == LUx); }
      else LUx = A;
    }
    LU_Decompose(LUx,P,det);
  }
#undef LUX
#undef APTR

  template <class T> LUDiv<T>::~LUDiv()
  { if (!inplace) delete[] Aptr; delete[] P; }

  template <class T> bool LUDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenMatrix<T>* mm = dynamic_cast<const GenMatrix<T>*>(&m);
    TMVAssert(mm);
    if (fout) {
      *fout << "M = "<< (istrans ? mm->Transpose() : mm->View()) <<endl;
      *fout << "L = "<<GetL()<<endl;
      *fout << "U = "<<GetU()<<endl;
    }
    Matrix<T> m2 = GetL()*GetU();
    m2.ReversePermuteRows(GetP());
    RealType(T) nm = Norm(m2- (istrans ? mm->Transpose() : mm->View()) );
    nm /= Norm(GetL())*Norm(GetU());
    if (fout) {
      *fout << "PLU = "<<m2<<endl;
      *fout << "Norm(M-PLU)/Norm(PLU) = "<<nm<<endl;
    }
    return nm < Epsilon<T>();
  }

  //
  // LDivEq
  //

  template <class T1, class T2> void NonLapLU_LDivEq(
      const GenMatrix<T1>& LUx, const MatrixView<T2>& m)
  {
    // Solve L U x = m:
    TMVAssert(LUx.IsSquare());
    TMVAssert(LUx.rowsize() == m.colsize());
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
      const GenMatrix<T1>& LUx, const size_t* P, 
      const MatrixView<T2>& m)
  {
    TMVAssert(m.colsize() == LUx.rowsize()); 
    TMVAssert(LUx.rowsize() == LUx.colsize());

    m.PermuteRows(P); 

#ifdef LAP
    if (m.iscm() && !m.isconj() && LUx.iscm() && !LUx.isconj())
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
  template <class T1, class T2> inline void LapLU_RDivEq(
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
      const GenMatrix<T1>& LUx, const size_t* P, 
      const MatrixView<T2>& m)
    // Solve x P L U = m:
  {
    TMVAssert(m.rowsize() == LUx.rowsize()); 
    TMVAssert(LUx.rowsize() == LUx.colsize());

#ifdef LAP
    if (LUx.iscm() && !LUx.isconj() && m.isrm())
      LapLU_RDivEq(LUx,m);
    else 
#endif
      NonLapLU_RDivEq(LUx,m); 

    m.ReversePermuteCols(P); 
  }

  template <class T> T LUDiv<T>::Det() const
  {
    if (!donedet) {
      det *= DiagMatrixViewOf(LUx.diag()).Det();
      donedet = true;
    }
    return det;
  }

  template <class T> void LUDiv<T>::Inverse(const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LUx.colsize());
    TMVAssert(minv.rowsize() == LUx.colsize());
    minv.SetToIdentity();
    LDivEq(minv);
  }

  template <class T> void LUDiv<T>::InverseATA(const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LUx.colsize());
    TMVAssert(minv.rowsize() == LUx.colsize());
    Matrix<T,ColMajor> temp(LUx.colsize(),LUx.colsize());
    Inverse(temp.View());
    minv = temp*temp.Adjoint();
  }

  //
  // Decompose
  //

  template <class T> void NonBlockLU_Decompose(
      const MatrixView<T>& A, size_t* P, T& det)
  {
    // LU Decompostion with partial pivoting.
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
    // U(0:j,0:j) and L(0:N,0:j) 
    // (where my a:b notation does not include the index b)
    //
    // The equation A = LU gives for the j col:
    //
    // A(0:N,j) = L(0:N,0:N) U(0:N,j)
    //
    // which breaks up into:
    //
    // (1) A(0:j,j) = L(0:j,0:N) U(0:N,j)
    // (2) A(j:N,j) = L(j:N,0:N) U(0:N,j)
    //
    // The first of these (1) simplifies to:
    // 
    // (1*) A(0:j,j) = L(0:j,0:j) U(0:j,j)
    //
    // since L is lower triangular, so L(0:j,j:N) = 0.
    // L(0:j,0:j) is already known, so this equation can be solved for
    // U(0:j,j) by forward substitution.
    // 
    // The second equation (2) simplifies to:
    //
    //      A(j:N,j) = L(j:N,0:j+1) U(0:j+1,j)
    // (2*)          = L(j:N,0:j) U(0:j,j) + L(j:N,j) U(j,j)
    //
    // since U is upper triangular so U(j+1:N,j) = 0.
    // Since we now know U(0:j,j) from (1*) above, this equation can
    // be solved for the product L(j:N,j) U(j,j)
    // 
    // This means we have some leeway on the values for L(j,j) and U(j,j),
    // as only their product is specified.
    //
    // If we take U to have unit diagonal, then L(j,j) is set here along
    // with the rest of the L(j:N,j) column.  However, this will mean that the
    // forward substutions in the (1*) steps will require divisions by the 
    // non-unit-diagonal elements of L.  It is faster to take L to have
    // unit-diagonal elements, and do the division by U(j,j) here, since then
    // we can calculate 1/U(j,j) and multiply.  So 1 division and N-j
    // multiplies which is generally faster than N-j divisions.
    //
    // However, another potential problem is that U(j,j) could be 0, or close
    // to 0.  This would lead to either an error or inaccurate results.
    // Thus we add a step in the middle of the (2*) calculation:
    //
    // Define v(j:N) = A(j:N,j) - L(j:N,0:j) U(0:j,j)
    // 
    // We search v for the element with the largest absolute value and apply
    // a permutation to swap it into the j spot.  This element then becomes 
    // U(j,j), which is then the divisor for the rest of the vector.  This 
    // will minimize the possibility of roundoff errors due to small U(j,j)'s.
    // 
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    const size_t N = A.rowsize();
    const size_t M = A.colsize();
    const size_t R = min(N,M);
#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    VIt<T,Step,NonConj> Ujj = A.diag().begin();
    size_t* Pj = P;

    for (size_t j=0; j<R; ++j,++Ujj,++Pj)
    {
      if (j > 0) {
	// Solve for U(0:j,j))
	A.col(j,0,j) /= LowerTriMatrixViewOf(A.SubMatrix(0,j,0,j),UnitDiag);

	// Solve for v = L(j:M,j) U(j,j)
	A.col(j,j,M) -= A.SubMatrix(j,M,0,j) * A.col(j,0,j);
      }

      // Find the pivot element
      size_t ip;
      A.col(j,j,M).MaxAbsElement(&ip);
      // ip is relative to j index, not absolute.

      // Swap the pivot row with j if necessary
      if (ip != 0) {
	ip += j;
	A.SwapRows(ip,j);  // This does both Lkb and A'
	*Pj = ip;
	det = -det;
      } else *Pj = j;

      // Solve for L(j+1:M,j)
      // If Ujj is 0, then all of the L's are 0.
      // ie. Ujj Lij = 0 for all i>j
      // Any value for Lij is valid, so leave them 0.
      if (*Ujj != T(0)) 
	A.col(j,j+1,M) /= *Ujj;
    }
    if (N > M) {
      // Solve for U(0:M,M:N))
      A.Cols(M,N) /= LowerTriMatrixViewOf(A.Cols(0,M),UnitDiag);
    }
#ifdef XDEBUG
    //cerr<<"A = "<<A<<endl;
    Matrix<T> L = A;
    UpperTriMatrixViewOf(L).OffDiag().Zero();
    L.diag().SetAllTo(T(1));
    //cerr<<"L = "<<L<<endl;
    Matrix<T> U = UpperTriMatrixViewOf(A,NonUnitDiag);
    //cerr<<"U = "<<U<<endl;
    Matrix<T> AA = L * U;
    AA.ReversePermuteRows(P,0,R);
    if (Norm(AA-A0) > 0.001*Norm(A0)) {
      cerr<<"Done NonBlock LU: \n";
      cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"LU = "<<A<<endl;
      cerr<<"P = (";
      for(size_t i=0;i<R;i++) cerr<<P[i]<<" ";
      cerr<<")\n";
      cerr<<"AA = "<<AA<<endl;
      cerr<<"Norm(A0-AA) = "<<Norm(AA-A0)<<endl;
      abort();
    }
#endif
  }

  template <class T> void BlockLU_Decompose(
      const MatrixView<T>& A, size_t* P, T& det)
  {
    // If A is large, we can take advantage of Blas Level 3 speed
    // by partitioning the matrix by columns.
    //
    // We do this calculation one block at a time down the diagonal. 
    // Each block contains k columns (expect for possibly the last 
    // block which may be fewer).
    //
    // For the first block, we decompose A into:
    //
    // ( A00 A01 ) = ( L00  0  ) ( U00 U01 )
    // ( A10 A11 )   ( L10  A~ ) (  0   1  )
    //
    // From this we obtain:
    //
    // (1) A00 = L00 U00
    // (2) A10 = L10 U00
    // (3) A01 = L00 U01
    // (4) A11 = L10 U01 + A~
    //
    // For (1) we decompose A00 in place using the non-blocked algorithm.
    // (2) and (3) then give us L10 and U01.
    // Finally, (4) lets us solve for A~.
    // Repeat until done.
    //
    // With pivoting, the only real change is to combine equations (1),(2)
    // and solve both together with the non-blocked algorithm.
    //
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T,ColMajor> A2 = A;
    size_t P2[A.colsize()];
    T det2=1;
    NonBlockLU_Decompose(A2.View(),P2,det2);
#endif

    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    const size_t N = A.rowsize();
    const size_t M = A.colsize();
    const size_t R = min(N,M);

    for (size_t jk=0; jk<R; jk+=LU_BLOCKSIZE)
    {
      size_t jkpk = min(jk+LU_BLOCKSIZE,R);
      // Solve for L00, U00

      NonBlockLU_Decompose(A.SubMatrix(jk,M,jk,jkpk),P+jk,det);

      // Apply the permutation to the rest of the matrix
      if (jk > 0) {
	A.SubMatrix(jk,M,0,jk).PermuteRows(P+jk,0,jkpk-jk);
      }
      if (jkpk < N) {
	A.SubMatrix(jk,M,jkpk,N).PermuteRows(P+jk,0,jkpk-jk);

	// Solve for U01
	A.SubMatrix(jk,jkpk,jkpk,N) /= 
	  LowerTriMatrixViewOf(A.SubMatrix(jk,jkpk,jk,jkpk),UnitDiag);

	// Solve for A~
	if (jkpk < M)
	  A.SubMatrix(jkpk,M,jkpk,N) -= A.SubMatrix(jkpk,M,jk,jkpk) *
	    A.SubMatrix(jk,jkpk,jkpk,N);
      }
      for(size_t i=jk;i<jkpk;++i) P[i]+=jk;
    }
#ifdef XDEBUG
    Matrix<T> L = A;
    UpperTriMatrixViewOf(L).OffDiag().Zero();
    L.diag().SetAllTo(T(1));
    Matrix<T> U = UpperTriMatrixViewOf(A,NonUnitDiag);
    Matrix<T> AA = L * U;
    AA.ReversePermuteRows(P,0,R);
    if (Norm(AA-A0) > 0.001*Norm(A0)) {
      cerr<<"Done Block LU: \n";
      cerr<<"A0 = "<<A0<<endl;
      cerr<<"LU = "<<A<<endl;
      cerr<<"P = (";
      for(size_t i=0;i<R;i++) cerr<<P[i]<<" ";
      cerr<<")\n";
      cerr<<"A2 = "<<A2<<endl;
      cerr<<"Norm(A-A2) = "<<Norm(A-A2)<<endl;
      cerr<<"Norm(AA-A0) = "<<Norm(AA-A0)<<endl;
      cerr<<"correct P = (";
      for(size_t i=0;i<R;i++) cerr<<P2[i]<<" ";
      cerr<<")\n";
      abort();
    }
#endif
  }
  template <class T> void RecursiveLU_Decompose(
      const MatrixView<T>& A, size_t* P, T& det)
  {
    // The recursive LU algorithm is similar to the block algorithm, except 
    // that the block is roughly half the size of the whole matrix.
    // We keep dividing the matrix in half (column-wise) until we get down
    // to an Mx2 or Mx1 matrix.

#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T,ColMajor> A2 = A;
    size_t P2[A.colsize()];
    T det2=1;
    NonBlockLU_Decompose(A2.View(),P2,det2);
#endif

    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    const size_t N = A.rowsize();
    const size_t M = A.colsize();
    const size_t R = min(N,M);

    if (R > 2) {
      // Split N in half, with N1 being rounded to multiple of BLOCKSIZE
      // if appropriate.
      size_t N1 = R/2;
      if (N1 > LU_BLOCKSIZE) N1 = (N1/LU_BLOCKSIZE)*LU_BLOCKSIZE;
      // Decompose left half into PLU
      RecursiveLU_Decompose(A.Cols(0,N1),P,det);

      // Apply the permutation to the right half of the matrix
      A.Cols(N1,N).PermuteRows(P,0,N1);

      // Solve for U01
      A.SubMatrix(0,N1,N1,N) /= 
	LowerTriMatrixViewOf(A.SubMatrix(0,N1,0,N1),UnitDiag);

      // Solve for A~
      A.SubMatrix(N1,M,N1,N) -= A.SubMatrix(N1,M,0,N1) *
	A.SubMatrix(0,N1,N1,N);

      // Decompose A~ into PLU
      RecursiveLU_Decompose(A.SubMatrix(N1,M,N1,N),P+N1,det);
      for(size_t i=N1;i<R;++i) P[i]+=N1;

      // Apply the new permutations to the left half
      A.Cols(0,N1).PermuteRows(P,N1,R);
    } else if (R == 2) {
      // Same as NonBlock version, but with R==2 hard coded
      VectorView<T> A0 = A.col(0);
      VectorView<T> A1 = A.col(1);

      // Find the pivot element
      size_t ip;
      RealType(T) piv = A0.MaxAbsElement(&ip);
      if (ip != 0) {
	A0.Swap(ip,0);
	A1.Swap(ip,0);
	*P = ip;
	det = -det;
      } else *P = 0;

      // Solve for L(1:M,0)
      if (piv != RealType(T)(0)) {
	A0.SubVector(1,M) /= A0(0);
	A1.SubVector(1,M) -= A0.SubVector(1,M) * A1(0);
      }

      piv = A1.SubVector(1,M).MaxAbsElement(&ip); 
      if (ip != 0) {
	++ip;
	A1.Swap(ip,1);
	A0.Swap(ip,1);
	P[1] = ip;
	det = -det;
      } else P[1] = 1;

      if (piv != RealType(T)(0)) 
	A1.SubVector(2,M) /= A1(1);

      if (N > 2) {
	// M=2, N>2, so solve for U(0:2,2:N))
	// A.Cols(2,N).PermuteRows(P);
	if (P[0] == 1) A.Cols(2,N).SwapRows(0,1);
	// A.Cols(2,N) /= LowerTriMatrixViewOf(A.Cols(0,2),UnitDiag);
	A.row(1,2,N) -= A(0,1) * A.row(0,2,N);
      }
    } else if (R == 1) {
      // Same as NonBlock version, but with R==1 hard coded
      VectorView<T> A0 = A.col(0);

      // Find the pivot element
      size_t ip;
      RealType(T) piv = A0.MaxAbsElement(&ip);
      if (ip != 0) {
	A0.Swap(ip,0);
	*P = ip;
	det = -det;
      } else *P = 0;

      // Solve for L(1:M,0)
      if (piv != RealType(T)(0)) {
	A0.SubVector(1,M) /= A0(0);
      }
    }
#ifdef XDEBUG
    Matrix<T> L = A;
    UpperTriMatrixViewOf(L).OffDiag().Zero();
    L.diag().SetAllTo(T(1));
    Matrix<T> U = UpperTriMatrixViewOf(A,NonUnitDiag);
    Matrix<T> AA = L * U;
    AA.ReversePermuteRows(P,0,R);
    if (Norm(AA-A0) > 0.001*Norm(A0)) {
      cerr<<"Done Recursive LU: \n";
      cerr<<"A0 = "<<A0<<endl;
      cerr<<"LU = "<<A<<endl;
      cerr<<"P = (";
      for(size_t i=0;i<R;i++) cerr<<P[i]<<" ";
      cerr<<")\n";
      cerr<<"A2 = "<<A2<<endl;
      cerr<<"Norm(A-A2) = "<<Norm(A-A2)<<endl;
      cerr<<"Norm(AA-A0) = "<<Norm(AA-A0)<<endl;
      cerr<<"correct P = (";
      for(size_t i=0;i<R;i++) cerr<<P2[i]<<" ";
      cerr<<")\n";
      abort();
    }
#endif
  }

  template <class T> void NonLapLU_Decompose(
      const MatrixView<T>& A, size_t* P, T& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);

    if (RecursiveLU)
      RecursiveLU_Decompose(A,P,det);
    else 
      if (A.rowsize() >= 2*LU_BLOCKSIZE) 
	BlockLU_Decompose(A,P,det);
      else
	NonBlockLU_Decompose(A,P,det);
  }

#ifdef LAP
  template <class T> inline void LapLU_Decompose(
      const MatrixView<T>& A, size_t* P, T& det)
  { NonLapLU_Decompose(A,P); }
  template <> inline void LapLU_Decompose(
      const MatrixView<double>& A, size_t* P, double& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int lap_p[A.colsize()];
    int m = A.colsize();
    int n = A.rowsize();
    int lda = A.stepj();
    int info;
    dgetrf(&m,&n,A.ptr(),&lda,lap_p,&info);
    if (info < 0) tmv_error("dgetrf returned info < 0");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = lap_p[i]-1;
      if (P[i]!=i) det = -det;
    }
  }
  template <> inline void LapLU_Decompose(
      const MatrixView<complex<double> >& A, size_t* P,
      complex<double>& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int lap_p[A.colsize()];
    int m = A.colsize();
    int n = A.rowsize();
    int lda = A.stepj();
    int info;
    zgetrf(&m,&n,LAP_Complex(A.ptr()),&lda,lap_p,&info);
    if (info < 0) tmv_error("zgetrf returned info < 0");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = lap_p[i]-1;
      if (P[i]!=i) det = -det;
    }
  }
#ifndef NOFLOAT
  template <> inline void LapLU_Decompose(
      const MatrixView<float>& A, size_t* P, float& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int lap_p[A.colsize()];
    int m = A.colsize();
    int n = A.rowsize();
    int lda = A.stepj();
    int info;
    sgetrf(&m,&n,A.ptr(),&lda,lap_p,&info);
    if (info < 0) tmv_error("sgetrf returned info < 0");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = lap_p[i]-1;
      if (P[i]!=i) det = -det;
    }
  }
  template <> inline void LapLU_Decompose(
      const MatrixView<complex<float> >& A, size_t* P, 
      complex<float>& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int lap_p[A.colsize()];
    int m = A.colsize();
    int n = A.rowsize();
    int lda = A.stepj();
    int info;
    cgetrf(&m,&n,LAP_Complex(A.ptr()),&lda,lap_p,&info);
    if (info < 0) tmv_error("cgetrf returned info < 0");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = lap_p[i]-1;
      if (P[i]!=i) det = -det;
    }
  }
#endif // NOFLOAT
#endif // LAP
  template <class T> inline void LU_Decompose(
      const MatrixView<T>& A, size_t* P, T& det)
  {
    TMVAssert(!A.isrm());
    TMVAssert(A.ct()==NonConj);
    if (A.colsize() > 0 && A.rowsize() > 0) {
#ifdef LAP
      LapLU_Decompose(A,P,det);
#else
      NonLapLU_Decompose(A,P,det);
#endif
    }
  }

#define InstFile "TMV_LUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


