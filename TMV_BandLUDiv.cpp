
#include "TMV_Band.h"
#include "TMV_Diag.h"
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

#define NEWLO min(A.nlo(),A.nhi())
#define NEWHI min(A.nlo()+A.nhi(),int(A.colsize())-1)
#define APTR inplace ? A.NonConst().ptr() : \
  new T[BandStorageLength(ColMajor,A.colsize(),A.colsize(),NEWLO,NEWHI)]

#define LUX istrans ? \
  inplace ? A.NonConst().Transpose() : \
  BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),A.nhi(), NEWHI, RowMajor) : \
  inplace ? A.NonConst().View() : \
  BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),A.nlo(), NEWHI, ColMajor)

  template <class T> BandLUDiv<T>::BandLUDiv(const GenBandMatrix<T>& A,
      bool _inplace) :
    istrans(A.nhi()<A.nlo() || (A.nhi()==A.nlo() && A.isrm())),
    inplace(_inplace || NEWLO == 0),
    Aptr(APTR), LUx(LUX),
    P(new size_t[A.colsize()]), det(T(1)), donedet(false)
  {
    TMVAssert(A.IsSquare());
    if (inplace) {
      // For inplace decomposition, make sure the original band matrix
      // has room for the extra upper diagonals...
      // if isrm stepi >= (2*A.nlo()+A.nhi())
      // if iscm stepj >= (2*A.nlo()+A.nhi())
      // if isdm extra diags appear at end, so can't really check
      TMVAssert(!LUx.isrm() || LUx.stepi()>=2*LUx.nlo()+LUx.nhi());
      TMVAssert(!LUx.iscm() || LUx.stepj()>=2*LUx.nlo()+LUx.nhi());
      TMVAssert(LUx == (istrans ? A.Transpose() : A.View()));
    } else {
      if (istrans) BandMatrixViewOf(LUx,A.nhi(),A.nlo()) = A.Transpose();
      else BandMatrixViewOf(LUx,A.nlo(),A.nhi()) = A;
    }
    if (LUx.nlo() > 0) {
      for(int i=LUx.nhi()-LUx.nlo()+1;i<=LUx.nhi();++i) LUx.diag(i).Zero();
      BandLU_Decompose(LUx,P,det);
    } else {
      P[0] = 1; // A tag to indicate that P is not set yet.
    }
  }

#undef LUX
#undef APTR
#undef NEWLO
#undef NEWHI

  template <class T> BandLUDiv<T>::~BandLUDiv()
  { delete [] P; if (!inplace) delete[] Aptr; }

  template <class T> bool BandLUDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenBandMatrix<T>* bm = dynamic_cast<const GenBandMatrix<T>*>(&m);
    TMVAssert(bm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*bm)<<"  ";
      *fout <<(istrans ? bm->Transpose() : bm->View()) <<endl;
      *fout << "L = "<<GetL()<<endl;
      *fout << "U = "<<GetU()<<endl;
    }
    Matrix<T> lu = GetL() * Matrix<T>(GetU());
    lu.ReversePermuteRows(GetP());
    RealType(T) nm = Norm(lu- (istrans ? bm->Transpose() : bm->View()) );
    nm /= Norm(GetL())*Norm(GetU());
    if (fout) {
      *fout << "PLU = "<<lu<<endl;
      *fout << "Norm(M-PLU)/Norm(PLU) = "<<nm<<endl;
    }
    BandMatrix<T> bm2 = *bm;
    bm2.DivideUsing(SVS);
    bm2.SetDiv();
    return nm < bm2.Condition()*bm2.colsize()*Epsilon<T>();
  }

  template <class T> T BandLUDiv<T>::Det() const
  {
    if (!donedet) {
      det *= DiagMatrixViewOf(LUx.diag()).Det();
      donedet = true;
    }         
    return det;  
  }                  

  template <class T> void BandLUDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LUx.colsize());
    TMVAssert(minv.rowsize() == LUx.colsize());
    minv.SetToIdentity();
    LDivEq(minv);
  }

  template <class T> void BandLUDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    Matrix<T,ColMajor> temp = Eye<T,ColMajor>(LUx.colsize());
    Inverse(temp.View());
    minv = temp*Transpose(temp);
  }

  template <class T> LowerTriMatrix<T,UnitDiag> BandLUDiv<T>::GetL() const
  {
    size_t N = LUx.colsize();
    int nlo = LUx.nlo();
    LowerTriMatrix<T,UnitDiag> L(N,T(0));
    if (nlo == 0) 
      L.SetToIdentity();
    else {
      for(size_t i=0;i<N;++i) {
	Swap(L.row(i,0,i),L.row(P[i],0,i));
	size_t end = min(i+nlo+1,N);
	L.col(i,i+1,end) = LUx.col(i,i+1,end);
      }
    }
    return L;
  }

  template <class T> const size_t* BandLUDiv<T>::GetP() const
  {
    if (LUx.nlo() == 0 && LUx.colsize() > 0 && P[0] == 1) {
      for(size_t i=0;i<LUx.colsize();++i) P[i] = i;
    }
    return P;
  }

  //
  // Decompose
  //

  // MJ: Is there a better way to do RowMajor and DiagMajor?
  template <class T> void NonLapBandLU_Decompose(
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
    const size_t N = A.colsize();
    const size_t M = A.rowsize();
    const size_t R = min(N,M);
    TMVAssert(A.nhi() > A.nlo());

    const int ds = A.diagstep();
    T* Ujj = A.ptr();
    size_t* Pj=P;

    size_t endcol = A.nlo()+1;
    size_t endrow = A.nhi()+1;
    for (size_t j=0; j<R; ++j,Ujj+=ds,++Pj)
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
      if (*Ujj != T(0)) 
	A.col(j,j+1,endcol) /= *Ujj;

      A.SubMatrix(j+1,endcol,j+1,endrow) -= 
        (A.col(j,j+1,endcol) ^ A.row(j,j+1,endrow));
      if (endcol < M) ++endcol;
      if (endrow < N) ++endrow;
    }
  }

#ifdef LAP
  template <class T> inline void LapBandLU_Decompose(
      const BandMatrixView<T>& A, size_t* P, T& det)
  { NonLapBandLU_Decompose(A,P,det); }
  template <> inline void LapBandLU_Decompose(
      const BandMatrixView<double>& A, size_t* P, double& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int kl = A.nlo();
    int ku = A.nhi()-kl;
    int lda = A.stepj()+1;
    int lap_p[A.colsize()];
    int info;
    dgbtrf(&m,&n,&kl,&ku,A.ptr()-A.nhi(),&lda,lap_p,&info);
    if (info < 0) tmv_error("dgbtrf returned info < 0");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = lap_p[i]-1;
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
    int m = A.colsize();
    int n = A.rowsize();
    int kl = A.nlo();
    int ku = A.nhi()-kl;
    int lda = A.stepj()+1;
    int lap_p[A.colsize()];
    int info;
    zgbtrf(&m,&n,&kl,&ku,LAP_Complex(A.ptr()-A.nhi()),&lda,lap_p,&info);
    if (info < 0) tmv_error("zgbtrf returned info < 0");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = lap_p[i]-1;
      if (P[i]!=i) det = -det;
    }
  }
#ifndef NOFLOAT
  template <> inline void LapBandLU_Decompose(
      const BandMatrixView<float>& A, size_t* P, float& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);
    int m = A.colsize();
    int n = A.rowsize();
    int kl = A.nlo();
    int ku = A.nhi()-kl;
    int lda = A.stepj()+1;
    int lap_p[A.colsize()];
    int info;
    sgbtrf(&m,&n,&kl,&ku,A.ptr()-A.nhi(),&lda,lap_p,&info);
    if (info < 0) tmv_error("sgbtrf returned info < 0");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = lap_p[i]-1;
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
    int m = A.colsize();
    int n = A.rowsize();
    int kl = A.nlo();
    int ku = A.nhi()-kl;
    int lda = A.stepj()+1;
    int lap_p[A.colsize()];
    int info;
    cgbtrf(&m,&n,&kl,&ku,LAP_Complex(A.ptr()-A.nhi()),&lda,lap_p,&info);
    if (info < 0) tmv_error("cgbtrf returned info < 0");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = lap_p[i]-1;
      if (P[i]!=i) det = -det;
    }
  }
#endif // NOFLOAT
  // Now Lap version of DiagMajor Tridiagonal:
  template <class T> inline void LapTriDiagLU_Decompose(
      const BandMatrixView<T>& A, size_t* P, T& det)
  { NonLapBandLU_Decompose(A,P,det); }
  template <> inline void LapTriDiagLU_Decompose(
      const BandMatrixView<double>& A, size_t* P, double& det)
  {
    TMVAssert(A.isdm());
    TMVAssert(A.IsSquare());
    TMVAssert(A.nlo()==1);
    TMVAssert(A.nhi()==2);
    TMVAssert(A.ct()==NonConj);
    int n = A.colsize();
    int lap_p[A.colsize()];
    int info;
    dgttrf(&n,A.diag(-1).ptr(),A.diag().ptr(),A.diag(1).ptr(),A.diag(2).ptr(),
	lap_p,&info);
    if (info < 0) tmv_error("dgttrf returned info < 0");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = lap_p[i]-1;
      if (P[i]!=i) det = -det;
    }
  }
  template <> inline void LapTriDiagLU_Decompose(
      const BandMatrixView<complex<double> >& A, size_t* P,
      complex<double>& det)
  {
    TMVAssert(A.isdm());
    TMVAssert(A.IsSquare());
    TMVAssert(A.nlo()==1);
    TMVAssert(A.nhi()==2);
    TMVAssert(A.ct()==NonConj);
    int n = A.colsize();
    int lap_p[A.colsize()];
    int info;
    zgttrf(&n,LAP_Complex(A.diag(-1).ptr()),LAP_Complex(A.diag().ptr()),
	LAP_Complex(A.diag(1).ptr()),LAP_Complex(A.diag(2).ptr()),
	lap_p,&info);
    if (info < 0) tmv_error("zgttrf returned info < 0");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = lap_p[i]-1;
      if (P[i]!=i) det = -det;
    }
  }
#ifndef NOFLOAT
  template <> inline void LapTriDiagLU_Decompose(
      const BandMatrixView<float>& A, size_t* P, float& det)
  {
    TMVAssert(A.isdm());
    TMVAssert(A.IsSquare());
    TMVAssert(A.nlo()==1);
    TMVAssert(A.nhi()==2);
    TMVAssert(A.ct()==NonConj);
    int n = A.colsize();
    int lap_p[A.colsize()];
    int info;
    sgttrf(&n,A.diag(-1).ptr(),A.diag().ptr(),A.diag(1).ptr(),A.diag(2).ptr(),
	lap_p,&info);
    if (info < 0) tmv_error("sgttrf returned info < 0");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = lap_p[i]-1;
      if (P[i]!=i) det = -det;
    }
  }
  template <> inline void LapTriDiagLU_Decompose(
      const BandMatrixView<complex<float> >& A, size_t* P,
      complex<float>& det)
  {
    TMVAssert(A.isdm());
    TMVAssert(A.IsSquare());
    TMVAssert(A.nlo()==1);
    TMVAssert(A.nhi()==2);
    TMVAssert(A.ct()==NonConj);
    int n = A.colsize();
    int lap_p[A.colsize()];
    int info;
    cgttrf(&n,LAP_Complex(A.diag(-1).ptr()),LAP_Complex(A.diag().ptr()),
	LAP_Complex(A.diag(1).ptr()),LAP_Complex(A.diag(2).ptr()),
	lap_p,&info);
    if (info < 0) tmv_error("cgttrf returned info < 0");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = lap_p[i]-1;
      if (P[i]!=i) det = -det;
    }
  }
#endif // NOFLOAT
#endif // LAP
  
  template <class T> inline void BandLU_Decompose(
      const BandMatrixView<T>& A, size_t* P, T& det)
  {
#ifdef XDEBUG
    //cerr<<"BandLU_Decompose: A = "<<Type(A)<<A<<endl;
    BandMatrix<T> A0 = A;
#ifdef LAP
    BandMatrix<T> A2 = A;
    size_t P2[A.colsize()];
    T det2 = det;
    NonLapBandLU_Decompose(A2.View(),P2,det2);
#endif
#endif

    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.nlo()>0);
    if (A.colsize() > 0 && A.rowsize() > 0) {
#ifdef LAP
      if (A.iscm())
	LapBandLU_Decompose(A,P,det);
      else if (A.isdm() && A.IsSquare() && A.nlo() == 1 && A.nhi() == 2)
	LapTriDiagLU_Decompose(A,P,det);
      else
#endif
	NonLapBandLU_Decompose(A,P,det);
    }

#ifdef XDEBUG
    size_t N = A.colsize();
    int nlo = A.nlo();
    LowerTriMatrix<T,UnitDiag> L(N,T(0));
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
      cerr<<"NonLap version = "<<A2<<endl;
      cerr<<"P2 = ";
      for(size_t i=0;i<N;i++) cerr<<P2[i]<<" ";
      cerr<<endl;
#endif
      abort();
    }
#endif
  }

  //
  // LDivEq
  //

  template <class T, class T1> void NonLapBandLU_LDivEq(
      const GenBandMatrix<T1>& LUx,
      const size_t* P, const MatrixView<T>& m) 
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
      const size_t* Pj = P;
      for(size_t j=0; j+1<N; ++j,++Pj) {
	m.SwapRows(j,*Pj);
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

    char c = 'N';
    int n = LUx.colsize();
    int kl = LUx.nlo();
    int ku = LUx.nhi()-LUx.nlo();
    int nrhs = m.rowsize();
    int lda = LUx.stepj()+1;
    int ldm = m.stepj();
    int info;
    int lap_p[n];
    for(int i=0;i<n;i++) lap_p[i] = P[i]+1;
    dgbtrs(&c,&n,&kl,&ku,&nrhs,const_cast<double*>(LUx.cptr()-LUx.nhi()),&lda,
	lap_p,m.ptr(),&ldm,&info);
    if (info < 0) tmv_error("dgbtrs returned info < 0");
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

    char c = 'N';
    int n = LUx.colsize();
    int kl = LUx.nlo();
    int ku = LUx.nhi()-LUx.nlo();
    int nrhs = m.rowsize();
    int lda = LUx.stepj()+1;
    int ldm = m.stepj();
    int info;
    int lap_p[n];
    for(int i=0;i<n;i++) lap_p[i] = P[i]+1;
    zgbtrs(&c,&n,&kl,&ku,&nrhs,LAP_Complex(LUx.cptr()-LUx.nhi()),&lda,
	lap_p,LAP_Complex(m.ptr()),&ldm,&info);
    if (info < 0) tmv_error("zgbtrs returned info < 0");
  }
#ifndef NOFLOAT
#endif
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

    char c = 'N';
    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int ldm = m.stepj();
    int info;
    int lap_p[n];
    for(int i=0;i<n;i++) lap_p[i] = P[i]+1;
    dgttrs(&c,&n,&nrhs,const_cast<double*>(LUx.diag(-1).cptr()),
	const_cast<double*>(LUx.diag().cptr()),
	const_cast<double*>(LUx.diag(1).cptr()),
	const_cast<double*>(LUx.diag(2).cptr()),lap_p,m.ptr(),&ldm,&info);
    if (info < 0) tmv_error("dgttrs returned info < 0");
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

    char c = 'N';
    int n = LUx.colsize();
    int nrhs = m.rowsize();
    int ldm = m.stepj();
    int info;
    int lap_p[n];
    for(int i=0;i<n;i++) lap_p[i] = P[i]+1;
    zgttrs(&c,&n,&nrhs,LAP_Complex(LUx.diag(-1).cptr()),
	LAP_Complex(LUx.diag().cptr()),LAP_Complex(LUx.diag(1).cptr()),
	LAP_Complex(LUx.diag(2).cptr()),lap_p,
	LAP_Complex(m.ptr()),&ldm,&info);
    if (info < 0) tmv_error("zgttrs returned info < 0");
  }
#ifndef NOFLOAT
#endif
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

  template <class T, class T1> void NonLapBandLU_RDivEq(
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

    char c = 'T';
    int n = LUx.colsize();
    int kl = LUx.nlo();
    int ku = LUx.nhi()-LUx.nlo();
    int nrhs = m.colsize();
    int lda = LUx.stepj()+1;
    int ldm = m.stepi();
    int info;
    int lap_p[n];
    for(int i=0;i<n;i++) lap_p[i] = P[i]+1;
    dgbtrs(&c,&n,&kl,&ku,&nrhs,const_cast<double*>(LUx.cptr()-LUx.nhi()),&lda,
	lap_p,m.ptr(),&ldm,&info);
    if (info < 0) tmv_error("dgbtrs returned info < 0");
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

    char c = LUx.isconj() ? 'C' : 'T';
    int n = LUx.colsize();
    int kl = LUx.nlo();
    int ku = LUx.nhi()-LUx.nlo();
    int nrhs = m.colsize();
    int lda = LUx.stepj()+1;
    int ldm = m.stepi();
    int info;
    int lap_p[n];
    for(int i=0;i<n;i++) lap_p[i] = P[i]+1;
    zgbtrs(&c,&n,&kl,&ku,&nrhs,LAP_Complex(LUx.cptr()-LUx.nhi()),&lda,
	lap_p,LAP_Complex(m.ptr()),&ldm,&info);
    if (info < 0) tmv_error("zgbtrs returned info < 0");
  }
#ifndef NOFLOAT
#endif
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

    char c = 'T';
    int n = LUx.colsize();
    int nrhs = m.colsize();
    int ldm = m.stepi();
    int info;
    int lap_p[n];
    for(int i=0;i<n;i++) lap_p[i] = P[i]+1;
    dgttrs(&c,&n,&nrhs,const_cast<double*>(LUx.diag(-1).cptr()),
	const_cast<double*>(LUx.diag().cptr()),
	const_cast<double*>(LUx.diag(1).cptr()),
	const_cast<double*>(LUx.diag(2).cptr()),lap_p,m.ptr(),&ldm,&info);
    if (info < 0) tmv_error("dgttrs returned info < 0");
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

    char c = LUx.isconj() ? 'C' : 'T';
    int n = LUx.colsize();
    int nrhs = m.colsize();
    int ldm = m.stepi();
    int info;
    int lap_p[n];
    for(int i=0;i<n;i++) lap_p[i] = P[i]+1;
    zgttrs(&c,&n,&nrhs,LAP_Complex(LUx.diag(-1).cptr()),
	LAP_Complex(LUx.diag().cptr()),LAP_Complex(LUx.diag(1).cptr()),
	LAP_Complex(LUx.diag(2).cptr()),lap_p,
	LAP_Complex(m.ptr()),&ldm,&info);
    if (info < 0) tmv_error("zgttrs returned info < 0");
  }
#ifndef NOFLOAT
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
    void DoRowUpperBandTriLDivEq(
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
	if (*Aii==Ta(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
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
	    tmv_error("Singular Matrix found in BandTriLDivEq");
	  *bi /= (ca ? CONJ(*Aii) : *Aii);
	}
	if (k > 0) { --k; ++len; }
      } 
    }

  template <bool rm, bool ua, class T, class Ta> void RowUpperBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isconj())
      DoRowUpperBandTriLDivEq<rm,true,ua>(A,b);
    else
      DoRowUpperBandTriLDivEq<rm,false,ua>(A,b);
  }
                  
  template <bool cm, bool ca, bool ua, class T, class Ta> 
    void DoColUpperBandTriLDivEq(
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
	      tmv_error("Singular Matrix found in BandTriLDivEq");
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
	if (*Ajj==Ta(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	*bj /= (ca ? CONJ(*Ajj) : *Ajj);
      }
    }

  template <bool cm, bool ua, class T, class Ta> void ColUpperBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isconj())
      DoColUpperBandTriLDivEq<cm,true,ua>(A,b);
    else
      DoColUpperBandTriLDivEq<cm,false,ua>(A,b);
  }
                  
  template <bool rm, bool ca, bool ua, class T, class Ta> 
    void DoRowLowerBandTriLDivEq(
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
	if (*Aij1==Ta(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
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
	    tmv_error("Singular Matrix found in BandTriLDivEq");
	  *bi /= (ca ? CONJ(*Aij) : *Aij);
	}
	if (k>0) { --k; ++len; Aij1+=A.stepi(); } 
	else { ++bj1; Aij1+=ds; }
      }
    }

  template <bool rm, bool ua, class T, class Ta> void RowLowerBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isconj())
      DoRowLowerBandTriLDivEq<rm,true,ua>(A,b);
    else
      DoRowLowerBandTriLDivEq<rm,false,ua>(A,b);
  }
                  
  template <bool cm, bool ca, bool ua, class T, class Ta> 
    void DoColLowerBandTriLDivEq(
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
	      tmv_error("Singular Matrix found in BandTriLDivEq");
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
	if (*Ajj==Ta(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	*bj /= (ca ? CONJ(*Ajj) : *Ajj);
      }
    }

  template <bool rm, bool ua, class T, class Ta> void ColLowerBandTriLDivEq(
      const GenBandMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isconj())
      DoColLowerBandTriLDivEq<rm,true,ua>(A,b);
    else
      DoColLowerBandTriLDivEq<rm,false,ua>(A,b);
  }
                  
  template <class T, class Ta> void UpperBandTriLDivEq(
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

  template <class T, class Ta> void LowerBandTriLDivEq(
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

  template <class T, class Ta> void NonBlasBandTriLDivEq(
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
      const VectorView<complex<float> >& b, DiagType dt)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(b.size() == A.colsize());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
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

  template <bool ua, class T, class Ta> void RowUpperBandTriLDivEq(
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
	if (*Aii==Ta(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	B.row(i) /= (A.isconj() ? CONJ(*Aii) : *Aii);
	if (k > 0) --k; else --N;
      }
    } 
  }

  template <bool ua, class T, class Ta> void ColUpperBandTriLDivEq(
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
	if (*Ajj==Ta(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	B.row(j) /= (A.isconj() ? CONJ(*Ajj) : *Ajj);
	B.Rows(i1,j) -= A.col(j,i1,j) ^ B.row(j);
	if (i1 > 0) --i1;
      }
    } 
  }

  template <bool ua, class T, class Ta> void RowLowerBandTriLDivEq(
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
	if (*Aii==Ta(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	B.row(i) /= (A.isconj() ? CONJ(*Aii) : *Aii);
	if (k>0) --k; else ++i1;
      }
    }
  }

  template <bool ua, class T, class Ta> void ColLowerBandTriLDivEq(
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
	if (*Ajj==Ta(0)) tmv_error("Singular Matrix found in BandTriLDivEq");
	B.row(j) /= (A.isconj() ? CONJ(*Ajj) : *Ajj);
	B.Rows(j+1,i2) -= A.col(j,j+1,i2) ^ B.row(j);
	if (i2 < N) ++i2;
      }
    }
  }

  template <bool ua, class T, class Ta> void NonLapBandTriLDivEq(
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
    char u,c;
    if (A.iscm()) {
      c = 'N';
      u = A.nlo() == 0 ? 'U' : 'L';
    } else {
      c = 'T';
      u = A.nhi() == 0 ? 'U' : 'L';
    }
    char d = dt==UnitDiag ? 'U' : 'N';
    int n = A.colsize();
    int kd = A.nlo()==0 ? A.nhi() : A.nlo();
    int nrhs = B.rowsize();
    double* a = const_cast<double*>(A.cptr());
    if (A.nlo()==0 && A.iscm()) a -= A.nhi();
    else if (A.nhi()==0 && A.isrm()) a -= A.nlo();
    int lda = A.iscm() ? A.stepj()+1 : A.stepi()+1;
    int ldb = B.stepj();
    int info;
    dtbtrs(&u,&c,&d,&n,&kd,&nrhs,a,&lda,B.ptr(),&ldb,&info);
    if (info < 0) tmv_error("dtbtrs returned info < 0");
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
    TMVAssert(A.isrm() || A.ct()==NonConj);

    char u,c;
    if (A.iscm()) {
      c = 'N';
      u = A.nlo() == 0 ? 'U' : 'L';
    } else {
      c = A.isconj() ? 'C' : 'T';
      u = A.nhi() == 0 ? 'U' : 'L';
    }
    char d = dt==UnitDiag ? 'U' : 'N';
    int n = A.colsize();
    int kd = A.nlo()==0 ? A.nhi() : A.nlo();
    int nrhs = B.rowsize();
    MKL_Complex16* a = LAP_Complex(A.cptr());
    if (A.nlo()==0 && A.iscm()) a -= A.nhi();
    else if (A.nhi()==0 && A.isrm()) a -= A.nlo();
    int lda = A.iscm() ? A.stepj()+1 : A.stepi()+1;
    int ldb = B.stepj();
    int info;
    ztbtrs(&u,&c,&d,&n,&kd,&nrhs,a,&lda,LAP_Complex(B.ptr()),&ldb,&info);
    if (info < 0) tmv_error("ztbtrs returned info < 0");
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
    char u,c;
    if (A.iscm()) {
      c = 'N';
      u = A.nlo() == 0 ? 'U' : 'L';
    } else {
      c = 'T';
      u = A.nhi() == 0 ? 'U' : 'L';
    }
    char d = dt==UnitDiag ? 'U' : 'N';
    int n = A.colsize();
    int kd = A.nlo()==0 ? A.nhi() : A.nlo();
    int nrhs = B.rowsize();
    float* a = const_cast<float*>(A.cptr());
    if (A.nlo()==0 && A.iscm()) a -= A.nhi();
    else if (A.nhi()==0 && A.isrm()) a -= A.nlo();
    int lda = A.iscm() ? A.stepj()+1 : A.stepi()+1;
    int ldb = B.stepj();
    int info;
    stbtrs(&u,&c,&d,&n,&kd,&nrhs,a,&lda,B.ptr(),&ldb,&info);
    if (info < 0) tmv_error("stbtrs returned info < 0");
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
    TMVAssert(A.isrm() || A.ct()==NonConj);

    char u,c;
    if (A.iscm()) {
      c = 'N';
      u = A.nlo() == 0 ? 'U' : 'L';
    } else {
      c = A.isconj() ? 'C' : 'T';
      u = A.nhi() == 0 ? 'U' : 'L';
    }
    char d = dt==UnitDiag ? 'U' : 'N';
    int n = A.colsize();
    int kd = A.nlo()==0 ? A.nhi() : A.nlo();
    int nrhs = B.rowsize();
    MKL_Complex8* a = LAP_Complex(A.cptr());
    if (A.nlo()==0 && A.iscm()) a -= A.nhi();
    else if (A.nhi()==0 && A.isrm()) a -= A.nlo();
    int lda = A.iscm() ? A.stepj()+1 : A.stepi()+1;
    int ldb = B.stepj();
    int info;
    ctbtrs(&u,&c,&d,&n,&kd,&nrhs,a,&lda,LAP_Complex(B.ptr()),&ldb,&info);
    if (info < 0) tmv_error("ctbtrs returned info < 0");
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
      if (((A.iscm() && A.stepj()>0) || (A.isrm() && A.stepi()>0)) &&
	  (B.iscm() && B.stepj()>0) && (A.isrm() || !A.isconj()))
	LapBandTriLDivEq(A,B,dt);
      else
#endif
	if (dt == UnitDiag) NonLapBandTriLDivEq<true>(A,B);
	else NonLapBandTriLDivEq<false>(A,B);

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

#define InstFile "TMV_BandLUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


