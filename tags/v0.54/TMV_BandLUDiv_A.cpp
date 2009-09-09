
#include "TMV_Band.h"

//#define XDEBUG

namespace tmv {

  //
  // Decompose
  //

  template <class T> inline void NonLapBandLU_Decompose(
      const BandMatrixView<T>& A, size_t* P, T& det)
  {
    // LU Decompostion with partial pivoting.
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
    // This makes the LDivEq, RDivEq functions a bit more complicated.
    //
    TMVAssert(A.nhi() > A.nlo());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.IsSquare());
    const size_t N = A.rowsize();

    const int ds = A.diagstep();
    T* Ujj = A.ptr();
    size_t* Pj=P;

    size_t endcol = A.nlo()+1;
    size_t endrow = A.nhi()+1;
    for (size_t j=0; j<N-1; ++j,Ujj+=ds,++Pj)
    {
      // Find the pivot element
      size_t ip;
      A.col(j,j,endcol).MaxAbsElement(&ip);
      if (ip != 0) {
	ip += j;
	Swap(A.row(ip,j,endrow),A.row(j,j,endrow));
	*Pj = ip;
	det = -det;
      } else *Pj = j;

      // If Ujj is 0, then all of the L's are 0.
      // ie. Ujj Lij = 0 for all i>j
      // Any value for Lij is valid, so leave them 0.
      if (*Ujj != T(0)) A.col(j,j+1,endcol) /= *Ujj;

      A.SubMatrix(j+1,endcol,j+1,endrow) -= 
        (A.col(j,j+1,endcol) ^ A.row(j,j+1,endrow));
      if (endcol < N) ++endcol;
      if (endrow < N) ++endrow;
    }
    // j == N-1
    *Pj = N-1;
  }

  template <class T> inline void NonLapTriDiagLU_Decompose(
      const BandMatrixView<T>& A, size_t* P, T& det)
  {
    // LU Decompostion for TriDiagonal BandMatrix
    //
    // For TriDiagonal BandMatrices, there are only two choices for 
    // each pivot, so we can specialize some of the calculations.
    // Otherwise, this is the same algorithm as above.
    //
    TMVAssert(A.isdm());
    TMVAssert(A.nlo() == 1);
    TMVAssert(A.nhi() == 2);
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.IsSquare());
    const size_t N = A.rowsize();

    T* Ujj = A.ptr(); // = U(j,j)
    T* Lj = Ujj+A.stepi();  // = L(j+1,j)
    T* Uj1 = Ujj+A.stepj(); // = U(j,j+1)
    T* Uj2 = Uj1+A.stepj(); // = U(j,j+2)
    size_t* Pj=P;

    for (size_t j=0; j<N-1; ++j,++Ujj,++Lj,++Uj1,++Uj2,++Pj)
    {
      bool pivot = abs(*Lj) > abs(*Ujj);

      if (pivot) {
	//Swap(A.row(j+1,j,endrow),A.row(j,j,endrow));
	swap(*Lj,*Ujj);
	if (j+1<N) swap(*(Ujj+1),*Uj1);
	if (j+2<N) swap(*(Uj1+1),*Uj2);
	*Pj = j+1;
	det = -det;
      } else *Pj = j;

      if (*Ujj != T(0)) *Lj /= *Ujj;

      //A.row(j+1,j+1,endrow) -= *Lj * A.row(j,j+1,endrow);
      *(Ujj+1) -= *Lj * (*Uj1);
      if (pivot && j+2<N) *(Uj1+1) -= *Lj * (*Uj2); // (if !pivot, *Uj2 == 0)
    }
    // j == N-1
    *Pj = N-1;
  }

#ifdef LAP
  template <class T> inline void LapBandLU_Decompose(
      const BandMatrixView<T>& A, size_t* P, T& det)
  { NonLapBandLU_Decompose(A,P,det); }
  template <> inline void LapBandLU_Decompose(
      const BandMatrixView<double>& A, size_t* P, double& det)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int n = A.rowsize();
    int kl = A.nlo();
    int ku = A.nhi()-kl;
    int lda = A.stepj()+1;
    auto_array<int> lap_p(new int[n]);
    LAPNAME(dgbtrf) (LAPCM LAPV(n),LAPV(n),LAPV(kl),LAPV(ku),
	LAPP(A.ptr()-A.nhi()),LAPV(lda),LAPP(lap_p.get()) LAPINFO );
    LAP_Results("dgbtrf");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = (lap_p.get())[i] LAPMINUS1;
      if (P[i]!=i) det = -det;
    }
  }
  template <> inline void LapBandLU_Decompose(
      const BandMatrixView<complex<double> >& A, size_t* P,
      complex<double>& det)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int n = A.rowsize();
    int kl = A.nlo();
    int ku = A.nhi()-kl;
    int lda = A.stepj()+1;
    auto_array<int> lap_p(new int[n]);
    LAPNAME(zgbtrf) (LAPCM LAPV(n),LAPV(n),LAPV(kl),LAPV(ku),
	LAPP(A.ptr()-A.nhi()),LAPV(lda),LAPP(lap_p.get()) LAPINFO );
    LAP_Results("zgbtrf");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = (lap_p.get())[i] LAPMINUS1;
      if (P[i]!=i) det = -det;
    }
  }
#ifndef NOFLOAT
  template <> inline void LapBandLU_Decompose(
      const BandMatrixView<float>& A, size_t* P, float& det)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int n = A.rowsize();
    int kl = A.nlo();
    int ku = A.nhi()-kl;
    int lda = A.stepj()+1;
    auto_array<int> lap_p(new int[n]);
    LAPNAME(sgbtrf) (LAPCM LAPV(n),LAPV(n),LAPV(kl),LAPV(ku),
	LAPP(A.ptr()-A.nhi()),LAPV(lda),LAPP(lap_p.get()) LAPINFO );
    LAP_Results("sgbtrf");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = (lap_p.get())[i] LAPMINUS1;
      if (P[i]!=i) det = -det;
    }
  }
  template <> inline void LapBandLU_Decompose(
      const BandMatrixView<complex<float> >& A, size_t* P,
      complex<float>& det)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int n = A.rowsize();
    int kl = A.nlo();
    int ku = A.nhi()-kl;
    int lda = A.stepj()+1;
    auto_array<int> lap_p(new int[n]);
    LAPNAME(cgbtrf) (LAPCM LAPV(n),LAPV(n),LAPV(kl),LAPV(ku),
	LAPP(A.ptr()-A.nhi()),LAPV(lda),LAPP(lap_p.get()) LAPINFO );
    LAP_Results("cgbtrf");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = (lap_p.get())[i] LAPMINUS1;
      if (P[i]!=i) det = -det;
    }
  }
#endif // NOFLOAT
  // Now Lap version of DiagMajor Tridiagonal:
  template <class T> inline void LapTriDiagLU_Decompose(
      const BandMatrixView<T>& A, size_t* P, T& det)
  { NonLapTriDiagLU_Decompose(A,P,det); }
  template <> inline void LapTriDiagLU_Decompose(
      const BandMatrixView<double>& A, size_t* P, double& det)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.isdm());
    TMVAssert(A.nlo()==1);
    TMVAssert(A.nhi()==2);
    TMVAssert(A.ct()==NonConj);
    int n = A.colsize();
    auto_array<int> lap_p(new int[n]);
    LAPNAME(dgttrf) (LAPCM LAPV(n),LAPP(A.diag(-1).ptr()),LAPP(A.diag().ptr()),
	LAPP(A.diag(1).ptr()),LAPP(A.diag(2).ptr()),LAPP(lap_p.get()) LAPINFO );
    LAP_Results("dgttrf");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = (lap_p.get())[i] LAPMINUS1;
      if (P[i]!=i) det = -det;
    }
  }
  template <> inline void LapTriDiagLU_Decompose(
      const BandMatrixView<complex<double> >& A, size_t* P,
      complex<double>& det)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.isdm());
    TMVAssert(A.nlo()==1);
    TMVAssert(A.nhi()==2);
    TMVAssert(A.ct()==NonConj);
    int n = A.colsize();
    auto_array<int> lap_p(new int[n]);
    LAPNAME(zgttrf) (LAPCM LAPV(n),LAPP(A.diag(-1).ptr()),LAPP(A.diag().ptr()), 
	LAPP(A.diag(1).ptr()),LAPP(A.diag(2).ptr()),LAPP(lap_p.get()) LAPINFO);
    LAP_Results("zgttrf");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = (lap_p.get())[i] LAPMINUS1;
      if (P[i]!=i) det = -det;
    }
  }
#ifndef NOFLOAT
  template <> inline void LapTriDiagLU_Decompose(
      const BandMatrixView<float>& A, size_t* P, float& det)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.isdm());
    TMVAssert(A.nlo()==1);
    TMVAssert(A.nhi()==2);
    TMVAssert(A.ct()==NonConj);
    int n = A.colsize();
    auto_array<int> lap_p(new int[n]);
    LAPNAME(sgttrf) (LAPCM LAPV(n),LAPP(A.diag(-1).ptr()),LAPP(A.diag().ptr()),
	LAPP(A.diag(1).ptr()),LAPP(A.diag(2).ptr()),LAPP(lap_p.get()) LAPINFO );
    LAP_Results("sgttrf");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = (lap_p.get())[i] LAPMINUS1;
      if (P[i]!=i) det = -det;
    }
  }
  template <> inline void LapTriDiagLU_Decompose(
      const BandMatrixView<complex<float> >& A, size_t* P,
      complex<float>& det)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.isdm());
    TMVAssert(A.nlo()==1);
    TMVAssert(A.nhi()==2);
    TMVAssert(A.ct()==NonConj);
    int n = A.colsize();
    auto_array<int> lap_p(new int[n]);
    LAPNAME(cgttrf) (LAPCM LAPV(n),LAPP(A.diag(-1).ptr()),LAPP(A.diag().ptr()), 
	LAPP(A.diag(1).ptr()),LAPP(A.diag(2).ptr()),LAPP(lap_p.get()) LAPINFO);
    LAP_Results("cgttrf");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = (lap_p.get())[i] LAPMINUS1;
      if (P[i]!=i) det = -det;
    }
  }
#endif // NOFLOAT
#endif // LAP
  
  template <class T> void BandLU_Decompose(
      const BandMatrixView<T>& A, size_t* P, T& det, int 
#ifdef LAP
      Anhi
#endif
      )
  {
#ifdef XDEBUG
    //cerr<<"BandLU_Decompose: A = "<<Type(A)<<"  "<<A<<endl;
    BandMatrix<T> A0 = A;
#endif

    TMVAssert(A.IsSquare());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.nlo()>0);
    TMVAssert(A.iscm() || (A.isdm() && A.nlo()==1 && A.nhi()==2));
    if (A.colsize() > 0 && A.rowsize() > 0) {
      if (A.iscm()) {
#ifdef LAP
	if (A.nlo()+Anhi+1 > int(A.rowsize())) {
	  TMVAssert(size_t(A.nhi()+1) == A.rowsize());
	  NonLapBandLU_Decompose(A,P,det);
	} else {
	  LapBandLU_Decompose(A,P,det);
	}
#else
	NonLapBandLU_Decompose(A,P,det);
#endif
      } else {
#ifdef LAP
	LapTriDiagLU_Decompose(A,P,det);
#else
	NonLapTriDiagLU_Decompose(A,P,det);
#endif
      }
    }

#ifdef XDEBUG
    size_t N = A.colsize();
    int nlo = A.nlo();
    Matrix<T> L(N,N,T(0));
    L.diag().SetAllTo(T(1));
    for(size_t i=0;i<N;++i) {
      Swap(L.row(i,0,i),L.row(P[i],0,i));
      size_t end = min(i+nlo+1,N);
      L.col(i,i+1,end) = A.col(i,i+1,end);
    }
    Matrix<T> U = BandMatrixViewOf(A,0,A.nhi());
    Matrix<T> AA = L*U;
    AA.ReversePermuteRows(P);
    if (!(Norm(AA-A0) < 0.001 * Norm(A0))) {
      cerr<<"BandLU_Decompose: A = "<<Type(A)<<A0<<endl;
      cerr<<"AA = "<<AA<<endl;
      cerr<<"Norm(diff) = "<<Norm(A0-AA)<<endl;
      cerr<<"BandLU = "<<A<<endl;
      cerr<<"L = "<<L<<endl;
      cerr<<"U = "<<U<<endl;
      cerr<<"P = ";
      for(size_t i=0;i<N;i++) cerr<<P[i]<<" ";
      cerr<<endl;
#ifdef LAP
      BandMatrix<T,ColMajor> A2 = A0;
      auto_array<size_t> P2(new size_t[A.colsize()]);
      T det2(0);
      NonLapBandLU_Decompose(A2.View(),P2.get(),det2);
      cerr<<"NonLap version = "<<A2<<endl;
      cerr<<"P2 = ";
      for(size_t i=0;i<N;i++) cerr<<(P2.get())[i]<<" ";
      cerr<<endl;
#endif
      abort();
    }
#endif
  }

#define InstFile "TMV_BandLUDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


