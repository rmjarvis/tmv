
#include "TMV_Band.h"
#include "TMV_VectorArith_Inline.h"
#include "TMV_Diag.h"

namespace tmv {

  template <class T> T BandLUDiv<T>::Det() const
  {
    if (!donedet) {
      const CVIt<T,Step,NonConj> _end = LUx.diag().end();
      for(CVIt<T,Step,NonConj> it = LUx.diag().begin(); it!=_end; ++it)
	det *= *it;
      donedet = true;
    }         
    return det;  
  }                  

  template <class T> Matrix<T,ColMajor> BandLUDiv<T>::Inverse() const
  {
    Matrix<T,ColMajor> temp = Eye<T,ColMajor>(LUx.colsize());
    LDivEq(temp.QuickView());
    return temp;
  }

  template <class T> Matrix<T,ColMajor> BandLUDiv<T>::InverseATA() const
  {
    Matrix<T,ColMajor> temp = Eye<T,ColMajor>(LUx.colsize());
    LDivEq(temp.QuickView());
    return temp*Transpose(temp);
  }

  template <class T> Matrix<T,ColMajor> GetLFromBandLU(
      const GenBandMatrix<T>& LUx, const vector<size_t>& p)
  {
    size_t N = LUx.colsize();
    int nlo = LUx.nlo();
    Matrix<T,ColMajor> L(N,N,T(0));
    for(size_t i=0;i<N;++i) {
      if (p[i] != i) L.SubMatrix(0,N,0,i).SwapRows(i,p[i]);
      size_t end = min(i+nlo+1,N);
      L.col(i).SubVector(i+1,end) = LUx.col(i,i+1,end);
      L(i,i) = T(1);
    }
    return L;
  }

  template <class T> Matrix<T,ColMajor> BandLUDiv<T>::LU_GetL() const 
  {
    if (istrans) 
      return Matrix<T,ColMajor>(BandMatrixViewOf(LUx,0,LUx.nhi()).Transpose()); 
    else return GetLFromBandLU(LUx,p);
  }

  template <class T> Matrix<T,ColMajor> BandLUDiv<T>::LU_GetU() const 
  {
    if (istrans) return Matrix<T,ColMajor>(GetLFromBandLU(LUx,p).Transpose());
    else return Matrix<T,ColMajor>(BandMatrixViewOf(LUx,0,LUx.nhi())); 
  }

  template <class T> Permutation BandLUDiv<T>::LU_GetP() const 
  {
    if (istrans) return Permutation(LUx.colsize());
    else {
      Permutation P(LUx.colsize());
      for(size_t i=0;i<p.size();++i) if (p[i] != i) P.SwapCols(i,p[i]);
      return P;
    }
  }

  template <class T> Permutation BandLUDiv<T>::LU_GetQ() const 
  {
    if (istrans) {
      Permutation Q(LUx.colsize());
      for(size_t i=0;i<p.size();++i) if (p[i] != i) Q.SwapRows(i,p[i]);
      return Q;
    }
    else return Permutation(LUx.colsize());
  }

  //
  // Decompose
  //

  template <class T> void BandLU_Decompose(
      const BandMatrixView<T>& A, vector<size_t>& p, T& det)
  {
    // LU Decompostion with partial pivoting and row scaling.  
    //
    // For band matrices, we use a somewhat different algorithm than for
    // regular matrices.  With regular matrices, the main operations were 
    // Matrix * Vector.  With the band matrix, these become difficult
    // to implement, since the submatrix that is doing the multiplying
    // goes off the edge of the bands.  So we choose to implement an
    // algorithm based on outerproduct updates which works better, since
    // the matrix being updated is completely within the band.
    //
    // On input, A contains the original band matrix with nlo subdiagonals
    // and nhi superdiagonals.  A must be created with nlo subdiagonals
    // and nlo+nhi superdiagonals.
    //
    // On each step, we calculate the j column of L, the j row of U
    // and the diagonal Ujj.  here is the first step:
    // 
    // A = L0 U0
    // ( A00 A0x ) = (  1  0 ) ( U00 U0x )
    // ( Ax0 Axx )   ( Lx0 I ) (  0  A'  )
    //
    // In Ax0, only A_10..A_nlo,0 are nonzero.
    // In A0x, only A_01..A_0,nhi are nonzero.
    // Axx is also band diagonal with the same nlo,nhi
    //
    // The formulae for L,U components are:
    //
    // U00 = A00
    // Lx0 = Ax0/U00
    // U0x = A0x
    // Uxx = Axx - Lx0 U0x
    //
    // It is apparent that Lx0 and U0x will have the same nonzero structure 
    // as Ax0 and A0x respectively.  This continues down the recursion,
    // so when we are done, L is lower banded with nlo subdiagonals, and
    // U is upper banded with nhi superdiagonals.
    //  
    // Unfortunately, this gets messed up a bit with pivoting.  
    // If the pivot element is as low as it can be, A_nlo,0, then 
    // swapping rows nlo and 0 will put A_nlo,nlo+nhi into A_0,nlo+nhi.
    // So we need to expand the upper band storage to nlo+nhi.
    //
    // The other problem is a bit more subtle.  Swapping rows j+nlo and j
    // also moves data from the j row down to the j+nlo in some of the 
    // previous columns which store data for L.  This would also screw up
    // the band structure, but we don't actually need to do this swap.
    // If we just keep track of what the swaps would be, we can just swap
    // rows in the remaining parts of A without swapping rows for L.
    // We simply keep track of the permutation swap for each step in a 
    // vector of swap indices.
    //
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == p.size());
    const size_t N = A.colsize();
    TMVAssert(A.nhi() > A.nlo());

    Vector<RealType(T)> scale(N);
    size_t j1=0;
    size_t j2=A.nhi()-A.nlo()+1;
    size_t k=A.nlo();
    for (size_t i=0; i<N; ++i) {
      scale(i) = MaxAbsElement(A.row(i,j1,j2));
      if (k>0) --k; else ++j1;
      if (j2<A.rowsize()) ++j2;
    }

    size_t endcol = A.nlo()+1;
    size_t endrow = A.nhi()+1;
    for (size_t j=0; j<N; ++j)
    {
      // Find the pivot element
      size_t i_pivot=j;
      RealType(T) pivot = abs(A(j,j));
      if (scale(j)) pivot /= scale(j);
      for(size_t i=j+1; i<endcol; ++i) {
	RealType(T) temp = abs(A(i,j));
	if (scale(i)) temp /= scale(i);
	if (temp > pivot) { pivot = temp; i_pivot = i; }
      }

      p[j] = i_pivot;
      // Swap the pivot row with j if necessary
      if (i_pivot != j) {
	Swap(A.row(i_pivot,j,endrow),A.row(j,j,endrow));
	scale.Swap(i_pivot,j);
	det = -det;
      }

      // If Ujj is 0, then all of the L's are 0.
      // ie. Ujj Lij = 0 for all i>j
      // Any value for Lij is valid, so leave them 0.
      const T& Ujj = A(j,j);
      if (Ujj != T(0)) {
	if (A.iscm()) A.col(j,j+1,endcol) /= Ujj;
	else A.col(j,j+1,endcol) /= Ujj;
      }

      A.SubMatrix(j+1,endcol,j+1,endrow) -= 
        (A.col(j,j+1,endcol) ^ A.row(j,j+1,endrow));
      if (endcol < N) ++endcol;
      if (endrow < N) ++endrow;
    }
  }

  //
  // LDivEq
  //

  template <class T1, class T2> void BandLU_LDivEq(
      const GenBandMatrix<T1>& LUx,
      const vector<size_t>& p, const MatrixView<T2>& m) 
  { 
    // Solve A x = m given that A = L U
    // L U x = m
    TMVAssert(m.colsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(p.size() == LUx.colsize());
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
    size_t j=0;
    size_t j1=1;      // j1 = j+1
    size_t jn=nlo+1;  // jn = j+nlo+1
    for(; j1<N; ++j,++j1) {
      if (p[j] != j) m.SwapRows(j,p[j]);
      m.Rows(j1,jn) -= LUx.col(j,j1,jn) ^ m.row(j);
      if (jn<N) ++jn;
    }

    // Next solve U x = y by back substitution
    BandTriLDivEq(BandMatrixViewOf(LUx,0,LUx.nhi()),m,NonUnitDiag);
  }

  //
  // RDivEq
  //

  template <class T1, class T2> void BandLU_RDivEq(
      const GenBandMatrix<T1>& LUx,
      const vector<size_t>& p, const MatrixView<T2>& m) 
  { 
    // Solve x A = m given that A = L U
    // x L U = m
    TMVAssert(m.rowsize() == LUx.colsize());
    TMVAssert(LUx.IsSquare());
    TMVAssert(p.size() == LUx.colsize());

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
      size_t j=N-2;
      size_t j1=N-1;
      size_t jn=N;
      size_t k=nlo-1;
      for(;j1>0;--j,--j1) {
	m.col(j) -= m.Cols(j1,jn) * LUx.col(j,j1,jn);
	if (p[j] != j) m.SwapCols(j,p[j]);
	if (k>0) --k; else --jn;
      }
    }
  }

  //
  // BandTriLDivEq
  //

  // MJ: convert to Aptr version
  template <class T1, class T2> void RowMajorUpperBandTriLDivEq(
      const GenBandMatrix<T1>& A, const VectorView<T2>& b, DiagType dt)
  {
    TMVAssert(A.isrm());
    TMVAssert(b.step()==1);
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.nhi() > 0);
    TMVAssert(b.ct() == NonConj);

    size_t N = b.size();
    VIt<T2,Unit,NonConj> bi = b.begin()+N-1;
    for(;N>0 && *bi==T2(0);--N,--bi);
    if (N==0) return;

    size_t k = A.nhi()-1;
    if (dt == UnitDiag) {
      if (N == 1) return;
      for(int i=N-2,ip1=N-1,len=1;i>=0;--i,--ip1) {
	T2 temp;
	if (A.isconj())
	  DoMultVV(temp,CVIt<T1,Unit,Conj>(A.row(i,ip1,N).begin()),
	      CVIt<T2,Unit,NonConj>(bi),len); // here, bi = b(i+1)
	else
	  DoMultVV(temp,CVIt<T1,Unit,NonConj>(A.row(i,ip1,N).begin()),
	      CVIt<T2,Unit,NonConj>(bi),len); // here, bi = b(i+1)
	*(--bi) -= temp; 
	if (k > 0) { --k; ++len; } else --N;
      }
    } else {
      CVIter<T1> Aii = A.diag().begin()+N-1;
      if (*Aii==T1(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
      *bi /= *Aii;
      for(int i=N-2,len=1;i>=0;--i) {
	T2 temp;
	if (A.isconj())
	  DoMultVV(temp,CVIt<T1,Unit,Conj>(A.row(i,i+1,N).begin()),
	      CVIt<T2,Unit,NonConj>(bi),len); // here, bi = b(i+1)
	else
	  DoMultVV(temp,CVIt<T1,Unit,NonConj>(A.row(i,i+1,N).begin()),
	      CVIt<T2,Unit,NonConj>(bi),len); // here, bi = b(i+1)
	*(--bi) -= temp;
	if (*(--Aii)==T1(0)) 
	  tmv_error("Singular Matrix found in BandTriLDivEq");
	*bi /= *Aii;
	if (k > 0) { --k; ++len; } else --N;
      }
    } 
  }

  template <class T1, class T2> void RowUpperBandTriLDivEq(
      const GenBandMatrix<T1>& A, const VectorView<T2>& b, DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.nhi() > 0);
    TMVAssert(b.ct() == NonConj);

    size_t N = b.size();
    VIt<T2,Step,NonConj> bi = b.begin()+N-1;
    for(;N>0 && *bi==T2(0);--N,--bi);
    if (N==0) return;

    size_t k = A.nhi();
    if (dt == UnitDiag) {
      for(int i=N-1;i>=0;--i,--bi) {
	*bi -= A.row(i,i+1,N) * b.SubVector(i+1,N);
	if (k > 0) --k; else --N;
      }
    } else {
      CVIter<T1> Aii = A.diag().begin()+N-1;
      for(int i=N-1;i>=0;--i,--bi,--Aii) {
	*bi -= A.row(i,i+1,N) * b.SubVector(i+1,N);
	if (*Aii==T1(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	*bi /= *Aii;
	if (k > 0) --k; else --N;
      }
    } 
  }

  // MJ: convert to Aptr version
  template <class T1, class T2> void ColMajorUpperBandTriLDivEq(
      const GenBandMatrix<T1>& A, const VectorView<T2>& b, DiagType dt)
  {
    TMVAssert(A.iscm());
    TMVAssert(b.step()==1);
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.nhi() > 0);
    TMVAssert(b.ct() == NonConj);

    int N = b.size();
    VIt<T2,Unit,NonConj> bj = b.begin()+N-1;
    for(;N>0 && *bj==T2(0);--N,--bj);
    if (N==0) return;

    int i1 = max(N-1-A.nhi(),int(0));
    VIt<T2,Unit,NonConj> bi1 = b.begin()+i1;
    if (dt == UnitDiag) {
      if (N<=i1+1) return;
      for(int j=N-1,len=j-i1;j>0&&len>0;--j,--bj) {
	if (*bj != T2(0)) {
	  if (A.isconj())
	    DoAddVV(-*bj,CVIt<T1,Unit,Conj>(A.col(j,i1,j).begin()),bi1,len);
	  else
	    DoAddVV(-*bj,CVIt<T1,Unit,NonConj>(A.col(j,i1,j).begin()),bi1,len);
	}
	if (i1 > 0) { --i1; --bi1; } else --len;
      }
    } else {
      CVIter<T1> Ajj = A.diag().begin()+N-1;
      for(int j=N-1,len=j-i1;j>=0;--j,--bj,--Ajj) {
	if (*bj != T2(0)) {
	  if (*Ajj==T1(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	  *bj /= *Ajj;
	  if (len>0) {
	    if (A.isconj())
	      DoAddVV(-*bj,CVIt<T1,Unit,Conj>(A.col(j,i1,j).begin()),bi1,len);
	    else
	      DoAddVV(-*bj,CVIt<T1,Unit,NonConj>(A.col(j,i1,j).begin()),bi1,
		  len);
	  }
	}
	if (i1 > 0) { --i1; --bi1; } else --len;
      }
    } 
  }

  // MJ: convert to Aptr version
  template <class T1, class T2> void RowMajorLowerBandTriLDivEq(
      const GenBandMatrix<T1>& A, const VectorView<T2>& b, DiagType dt)
  {
    TMVAssert(A.isrm());
    TMVAssert(b.step()==1);
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nhi() == 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.nlo() > 0);
    TMVAssert(b.ct() == NonConj);

    const size_t N = A.colsize();
    size_t i1=0; // all b(i<i1) == 0
    VIt<T2,Unit,NonConj> bi1 = b.begin();
    for(;i1<N && *bi1==T2(0);++i1,++bi1);
    if (i1==N) return;

    VIt<T2,Unit,NonConj> bi = bi1+1;
    size_t k=A.nlo()-1;
    if (dt == UnitDiag) {
      for(size_t i=i1+1,len=1;i<N;++i,++bi) {
	T2 temp;
	if (A.isconj())
	  DoMultVV(temp,CVIt<T1,Unit,Conj>(A.row(i,i1,i).begin()),
	      CVIt<T2,Unit,NonConj>(bi1),len);
	else
	  DoMultVV(temp,CVIt<T1,Unit,NonConj>(A.row(i,i1,i).begin()),
	      CVIt<T2,Unit,NonConj>(bi1),len);
	*bi -= temp;
	if (k>0) { --k; ++len; } else { ++i1; ++bi1; }
      }
    } else {
      CVIter<T1> Aii = A.diag().begin()+i1;
      if (*Aii==T1(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
      *bi1 /= *Aii;
      for(size_t i=i1+1,len=1;i<N;++i,++bi) {
	T2 temp;
	if (A.isconj())
	  DoMultVV(temp,CVIt<T1,Unit,Conj>(A.row(i,i1,i).begin()),
	      CVIt<T2,Unit,NonConj>(bi1),len);
	else
	  DoMultVV(temp,CVIt<T1,Unit,NonConj>(A.row(i,i1,i).begin()),
	      CVIt<T2,Unit,NonConj>(bi1),len);
	*bi -= temp;
	if (*(++Aii)==T1(0)) 
	  tmv_error("Singular Matrix found in BandTriLDivEq");
	*bi /= *Aii;
	if (k>0) { --k; ++len; } else { ++i1; ++bi1; }
      }
    }
  }

  template <class T1, class T2> void RowLowerBandTriLDivEq(
      const GenBandMatrix<T1>& A, const VectorView<T2>& b, DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nhi() == 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.nlo() > 0);
    TMVAssert(b.ct() == NonConj);

    const size_t N = A.colsize();
    size_t i1=0; // all b(i<i1) == 0
    VIt<T2,Step,NonConj> bi = b.begin()+i1;
    for(;i1<N && *bi==T2(0);++i1,++bi);
    if (i1==N) return;

    size_t k=A.nlo();
    if (dt == UnitDiag) {
      for(size_t i=i1;i<N;++i,++bi) {
	*bi -= A.row(i,i1,i) * b.SubVector(i1,i);
	if (k>0) --k; else ++i1;
      }
    } else {
      CVIter<T1> Aii = A.diag().begin()+i1;
      for(size_t i=i1;i<N;++i,++bi,++Aii) {
	*bi -= A.row(i,i1,i) * b.SubVector(i1,i);
	if (*Aii==T1(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	*bi /= *Aii;
	if (k>0) --k; else ++i1;
      }
    }
  }

  // MJ: convert to Aptr version
  template <class T1, class T2> void ColMajorLowerBandTriLDivEq(
      const GenBandMatrix<T1>& A, const VectorView<T2>& b, DiagType dt)
  {
    TMVAssert(A.iscm());
    TMVAssert(b.step()==1);
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nhi() == 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.nlo() > 0);
    TMVAssert(b.ct() == NonConj);

    const size_t N = A.colsize();
    size_t i1=0; // all b(i<i1) == 0
    VIt<T2,Unit,NonConj> bj = b.begin()+i1;
    for(;i1<N && *bj==T2(0);++i1,++bj);
    if (i1==N) return;
    VIt<T2,Unit,NonConj> bjp1 = bj+1;

    size_t i2=min(i1+A.nlo()+1,A.colsize());
    if (dt == UnitDiag) {
      for(size_t j=i1,jp1=j+1,len=i2-jp1;j<N&&len>0;++j,++jp1,++bj,++bjp1) {
	if (*bj != T2(0)) {
	  if (A.isconj())
	    DoAddVV(-*bj,CVIt<T1,Unit,Conj>(A.col(j,jp1,i2).begin()),bjp1,len);
	  else
	    DoAddVV(-*bj,CVIt<T1,Unit,NonConj>(A.col(j,jp1,i2).begin()),bjp1,len);
	}
	if (i2 < N) ++i2; else --len;
      }
    } else {
      CVIter<T1> Ajj = A.diag().begin()+i1;
      size_t j=i1;
      for(size_t jp1=j+1,len=i2-jp1;j<N&&len>0;++j,++jp1,++bj,++bjp1,++Ajj) {
	if (*bj != T2(0)) {
	  if (*Ajj==T1(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	  *bj /= *Ajj;
	  if (len > 0) {
	    if (A.isconj())
	      DoAddVV(-*bj,CVIt<T1,Unit,NonConj>(A.col(j,jp1,i2).begin()),bjp1,
		  len);
	    else
	      DoAddVV(-*bj,CVIt<T1,Unit,Conj>(A.col(j,jp1,i2).begin()),bjp1,
		  len);
	  }
	}
	if (i2 < N) ++i2; else --len;
      }
      for(;j<N;++j,++bj,++Ajj) {
	if (*bj != T2(0)) {
	  if (*Ajj==T1(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	  *bj /= *Ajj;
	}
      }
    }
  }

  template <class T1, class T2> void NonBlasBandTriLDivEq(
      const GenBandMatrix<T1>& A, const VectorView<T2>& b, DiagType dt)
  {
    // Solve A x = y  where A is an upper or lower band triangle matrix
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.nlo() > 0 || A.nhi() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(b.ct() == NonConj);

    if (A.nlo() == 0) {
      if (b.step() == 1) {
	if (A.iscm()) ColMajorUpperBandTriLDivEq(A,b,dt);
	else if (A.isrm()) RowMajorUpperBandTriLDivEq(A,b,dt);
	else RowUpperBandTriLDivEq(A,b,dt);
      }
      else RowUpperBandTriLDivEq(A,b,dt);
    }
    else {
      if (b.step() == 1) {
	if (A.iscm()) ColMajorLowerBandTriLDivEq(A,b,dt);
	else if (A.isrm()) RowMajorLowerBandTriLDivEq(A,b,dt);
	else RowLowerBandTriLDivEq(A,b,dt);
      }
      else RowLowerBandTriLDivEq(A,b,dt);
    }
  }

#ifdef BLAS
  template <class T1, class T2> inline void BlasBandTriLDivEq(
      const GenBandMatrix<T1>& A, const VectorView<T2>& b, DiagType dt)
  { NonBlasBandTriLDivEq(A,b,dt); }
  template <> inline void BlasBandTriLDivEq(
      const GenBandMatrix<double>& A, const VectorView<double>& b,
      DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.nlo() > 0 || A.nhi() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    cblas_dtbsv(
	A.isrm() ? CblasRowMajor : CblasColMajor,
	A.nlo()==0 ? CblasUpper : CblasLower, CblasNoTrans,
	dt==UnitDiag ? CblasUnit : CblasNonUnit, A.colsize(),
	A.nlo()==0 ? A.nhi() : A.nlo(),
	A.isrm() ? A.cptr()-A.nlo() : A.cptr()-A.nhi(),
	A.diagstep(), b.ptr(), b.step());
  }
  template <> inline void BlasBandTriLDivEq(
      const GenBandMatrix<complex<double> >& A,
      const VectorView<complex<double> >& b, DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.nlo() > 0 || A.nhi() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(b.ct() == NonConj);
    if (A.isconj())
      cblas_ztbsv(
	  A.isrm() ? CblasColMajor : CblasRowMajor,
	  A.nlo()==0 ? CblasLower : CblasUpper, CblasConjTrans,
	  dt==UnitDiag ? CblasUnit : CblasNonUnit, A.colsize(),
	  A.nlo()==0 ? A.nhi() : A.nlo(),
	  A.isrm() ? A.cptr()-A.nlo() : A.cptr()-A.nhi(),
	  A.diagstep(), b.ptr(), b.step());
    else
      cblas_ztbsv(
	  A.isrm() ? CblasRowMajor : CblasColMajor,
	  A.nlo()==0 ? CblasUpper : CblasLower, CblasNoTrans,
	  dt==UnitDiag ? CblasUnit : CblasNonUnit, A.colsize(),
	  A.nlo()==0 ? A.nhi() : A.nlo(),
	  A.isrm() ? A.cptr()-A.nlo() : A.cptr()-A.nhi(),
	  A.diagstep(), b.ptr(), b.step());
  }
  template <> inline void BlasBandTriLDivEq(
      const GenBandMatrix<double>& A,
      const VectorView<complex<double> >& b, DiagType dt)
  {
    TMVAssert(b.ct() == NonConj);
    BlasBandTriLDivEq(A,b.Real(),dt);
    BlasBandTriLDivEq(A,b.Imag(),dt);
  }
#ifndef NOFLOAT
  template <> inline void BlasBandTriLDivEq(
      const GenBandMatrix<float>& A, const VectorView<float>& b,
      DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.nlo() > 0 || A.nhi() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    cblas_stbsv(
	A.isrm() ? CblasRowMajor : CblasColMajor,
	A.nlo()==0 ? CblasUpper : CblasLower, CblasNoTrans,
	dt==UnitDiag ? CblasUnit : CblasNonUnit, A.colsize(),
	A.nlo()==0 ? A.nhi() : A.nlo(),
	A.isrm() ? A.cptr()-A.nlo() : A.cptr()-A.nhi(),
	A.diagstep(), b.ptr(), b.step());
  }
  template <> inline void BlasBandTriLDivEq(
      const GenBandMatrix<complex<float> >& A,
      const VectorView<complex<float> >& b, 
      DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.nlo() > 0 || A.nhi() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(b.ct() == NonConj);
    if (A.isconj())
      cblas_ctbsv(
	  A.isrm() ? CblasColMajor : CblasRowMajor,
	  A.nlo()==0 ? CblasLower : CblasUpper, CblasConjTrans,
	  dt==UnitDiag ? CblasUnit : CblasNonUnit, A.colsize(),
	  A.nlo()==0 ? A.nhi() : A.nlo(),
	  A.isrm() ? A.cptr()-A.nlo() : A.cptr()-A.nhi(),
	  A.diagstep(), b.ptr(), b.step());
    else
      cblas_ctbsv(
	  A.isrm() ? CblasRowMajor : CblasColMajor,
	  A.nlo()==0 ? CblasUpper : CblasLower, CblasNoTrans,
	  dt==UnitDiag ? CblasUnit : CblasNonUnit, A.colsize(),
	  A.nlo()==0 ? A.nhi() : A.nlo(),
	  A.isrm() ? A.cptr()-A.nlo() : A.cptr()-A.nhi(),
	  A.diagstep(), b.ptr(), b.step());
  }
  template <> inline void BlasBandTriLDivEq(
      const GenBandMatrix<float>& A,
      const VectorView<complex<float> >& b, DiagType dt)
  {
    TMVAssert(b.ct() == NonConj);
    BlasBandTriLDivEq(A,b.Real(),dt);
    BlasBandTriLDivEq(A,b.Imag(),dt);
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T1, class T2> void BandTriLDivEq(
      const GenBandMatrix<T1>& A, const VectorView<T2>& b, DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.nhi() >= 0);
    TMVAssert(A.nlo() >= 0);
    TMVAssert(A.nhi() < int(A.colsize()));
    TMVAssert(A.nlo() < int(A.colsize()));
    if (A.nlo() == 0 && A.nhi() == 0) {
      if (dt == NonUnitDiag) b /= DiagMatrixViewOf(A.diag());
    } else {
      if (b.isconj())
	BandTriLDivEq(A.QuickConjugate(),b.Conjugate(),dt);
      else
#ifdef BLAS
	if (A.isrm() || A.iscm()) BlasBandTriLDivEq(A,b,dt);
	else 
#endif
	  NonBlasBandTriLDivEq(A,b,dt);
    }
  }

  template <class T1, class T2> void RowUpperBandTriLDivEq(
      const GenBandMatrix<T1>& A, const MatrixView<T2>& B, DiagType dt)
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
    if (dt == UnitDiag) {
      for(int i=N-1; i>=0; --i) {
	B.row(i) -= A.row(i,i+1,N) * B.Rows(i+1,N);
	if (k > 0) --k; else --N;
      }
    } else {
      CVIter<T1> Aii = A.diag().begin()+N-1;
      for(int i=N-1; i>=0; --i,--Aii) {
	B.row(i) -= A.row(i,i+1,N) * B.Rows(i+1,N);
	if (*Aii==T1(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	B.row(i) /= *Aii;
	if (k > 0) --k; else --N;
      }
    } 
  }

  template <class T1, class T2> void ColUpperBandTriLDivEq(
      const GenBandMatrix<T1>& A, const MatrixView<T2>& B, DiagType dt)
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
    if (dt == UnitDiag) {
      for(int j=N-1; j>0; --j) {
	B.Rows(i1,j) -= A.col(j,i1,j) ^ B.row(j);
	if (i1 > 0) --i1;
      }
    } else {
      CVIter<T1> Ajj = A.diag().begin()+N-1;
      for(int j=N-1; j>=0; --j,--Ajj) {
	if (*Ajj==T1(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	B.row(j) /= *Ajj;
	B.Rows(i1,j) -= A.col(j,i1,j) ^ B.row(j);
	if (i1 > 0) --i1;
      }
    } 
  }

  template <class T1, class T2> void RowLowerBandTriLDivEq(
      const GenBandMatrix<T1>& A, const MatrixView<T2>& B, DiagType dt)
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
    if (dt == UnitDiag) {
      for(size_t i=0; i<N; ++i) {
	B.row(i) -= A.row(i,i1,i) * B.Rows(i1,i);
	if (k>0) --k; else ++i1;
      }
    } else {
      CVIter<T1> Aii = A.diag().begin();
      for(size_t i=0; i<N; ++i,++Aii) {
	B.row(i) -= A.row(i,i1,i) * B.Rows(i1,i);
	if (*Aii==T1(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	B.row(i) /= *Aii;
	if (k>0) --k; else ++i1;
      }
    }
  }

  template <class T1, class T2> void ColLowerBandTriLDivEq(
      const GenBandMatrix<T1>& A, const MatrixView<T2>& B, DiagType dt)
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
    if (dt == UnitDiag) {
      for(size_t j=0; j<N; ++j) {
	B.Rows(j+1,i2) -= A.col(j,j+1,i2) ^ B.row(j);
	if (i2 < N) ++i2;
      }
    } else {
      CVIter<T1> Ajj = A.diag().begin();
      for(size_t j=0; j<N; ++j,++Ajj) {
	if (*Ajj==T1(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	B.row(j) /= *Ajj;
	B.Rows(j+1,i2) -= A.col(j,j+1,i2) ^ B.row(j);
	if (i2 < N) ++i2;
      }
    }
  }

  // MJ: There is an LAPack function for this: tbtrs
  template <class T1, class T2> void BandTriLDivEq(
      const GenBandMatrix<T1>& A, const MatrixView<T2>& B, DiagType dt)
  {
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
      BandTriLDivEq(A.QuickConjugate(),B.QuickConjugate(),dt);
    else if (B.isrm()) {
      if (A.nlo()==0)
	if (A.isrm()) RowUpperBandTriLDivEq(A,B,dt);
	else if (A.iscm()) ColUpperBandTriLDivEq(A,B,dt);
	else RowUpperBandTriLDivEq(A,B,dt);
      else
	if (A.isrm()) RowLowerBandTriLDivEq(A,B,dt);
	else if (A.iscm()) ColLowerBandTriLDivEq(A,B,dt);
	else RowLowerBandTriLDivEq(A,B,dt);
    } else {
      for(size_t j=0;j<B.rowsize();++j) 
	BandTriLDivEq(A,B.col(j),dt);
    }
  }

#define InstFile "TMV_BandLUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


