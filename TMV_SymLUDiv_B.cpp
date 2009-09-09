
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"
#include "TMV_Diag.h"
#include "TMV_SymLUDiv_A.h"

//#define XDEBUG

namespace tmv {

  //
  // LDivEq
  //

  template <class T, class T1ab, class T1c> inline void LMultEq_2x2(
      T1ab a, T1ab b, T1c c, T1c cc, T& x, T& y)
  {
    // Solve [ x ] <- [ a  cc ] [ x ]
    //       [ y ]    [ c  b  ] [ y ]
    T tempx = x;
    x = (a*x+cc*y);
    y = (b*y+c*tempx);
  }

  template <class T, class T1ab, class T1c> inline void LMultEq_2x2(
      T1ab a, T1ab b, T1c c, T1c cc, const MatrixView<T>& m)
  {
    // Solve m <- [ a  cc ] m 
    //            [ c  b  ]
    TMVAssert(m.colsize() == 2);
    TMVAssert(m.ct() == NonConj);
    if (m.isrm()) {
      T* m0 = m.ptr();
      T* m1 = m0 + m.stepi();
      for (size_t k=m.rowsize();k>0;--k,++m0,++m1) 
	LMultEq_2x2(a,b,c,cc,*m0,*m1);
    } else {
      const int sj = m.stepj();
      T* m0 = m.ptr();
      T* m1 = m0 + m.stepi();
      for (size_t k=m.rowsize();k>0;--k,m0+=sj,m1+=sj) 
	LMultEq_2x2(a,b,c,cc,*m0,*m1);
    }
  }

  template <bool herm, class T, class T1> inline void SymLMultEq_2x2(
      T1 a, T1 b, T1 c, const MatrixView<T>& m)
  {
    TMVAssert(m.colsize() == 2);
    if (herm)
      if (m.isconj())
	LMultEq_2x2(REAL(a),REAL(b),CONJ(c),c,m.Conjugate());
      else
	LMultEq_2x2(REAL(a),REAL(b),c,CONJ(c),m);
    else
      if (m.isconj())
	LMultEq_2x2(CONJ(a),CONJ(b),CONJ(c),CONJ(c),m.Conjugate());
      else
	LMultEq_2x2(a,b,c,c,m);
  }

  template <bool herm, class T, class T1> inline void SymRMultEq_2x2(
      T1 a, T1 b, T1 c, const MatrixView<T>& m)
  {
    TMVAssert(m.rowsize() == 2);
    if (herm)
      if (m.isconj())
	LMultEq_2x2(REAL(a),REAL(b),c,CONJ(c),m.Adjoint());
      else
	LMultEq_2x2(REAL(a),REAL(b),CONJ(c),c,m.Transpose());
    else
      if (m.isconj())
	LMultEq_2x2(CONJ(a),CONJ(b),CONJ(c),CONJ(c),m.Adjoint());
      else
	LMultEq_2x2(a,b,c,c,m.Transpose());
  }

  template <bool herm, class T, class T1> void PseudoDiag_LDivEq(
      const GenVector<T1>& D, const GenVector<T1>& xD,
      const MatrixView<T>& m)
  {
    TMVAssert(D.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());
    TMVAssert(D.ct() == NonConj);
    TMVAssert(xD.ct() == NonConj);

#ifdef XDEBUG
    //cerr<<"Start PseudoDiag_LDivEq\n";
    //cerr<<"xD = "<<xD<<endl;
    Matrix<T> m0 = m;
    Matrix<T1> DD(D.size(),D.size(),T1(0));
    DD.diag() = D;
    DD.diag(-1) = xD;
    DD.diag(1) = herm ? xD.Conjugate() : xD.View();
    //cerr<<"xD = "<<xD<<endl;
#endif

    const T1* Di = D.cptr();
    const T1* xDi = xD.cptr();

    const size_t N = D.size();
    const int sd = D.step();
    const int sx = xD.step();

    for(size_t i=0;i<N;) {
      if (i==N-1 || *xDi == T1(0)) {
	if (herm)
	  m.row(i) /= REAL(*Di);
	else
	  m.row(i) /= *Di;
	Di+=sd; xDi+=sx; ++i;
      } else {
	T1 x = *Di;
	T1 y = *(Di+=sd);
	T1 z = *xDi;
	SymInvert_2x2<herm>(x,y,z);
	SymLMultEq_2x2<herm>(x,y,z,m.Rows(i,i+2));
	Di+=sd,xDi+=2*sx,i+=2;
      }
    }
    //cerr<<"xD = "<<xD<<endl;
#ifdef XDEBUG
    Matrix<T> m2 = DD * m;
    if (Norm(m2-m0) > 0.001*Norm(DD)*Norm(m0)) {
      cerr<<"PseudoDiag_LDivEq: m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"D = "<<D<<endl;
      cerr<<"xD = "<<xD<<endl;
      cerr<<"DD = "<<DD<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"m2 = DD*m "<<m2<<endl;
      abort();
    }
    //cerr<<"Done PseudoDiag_LDivEq: xD = "<<xD<<endl;
#endif
  }

  template <bool herm, class T, class T1> void PseudoDiag_LMultEq(
      const GenVector<T1>& D, const GenVector<T1>& xD,
      const MatrixView<T>& m)
  {
    TMVAssert(D.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());
    TMVAssert(D.ct() == NonConj);
    TMVAssert(xD.ct() == NonConj);

    //cerr<<"Start PseudoDiag_LMultEq: xD = "<<xD<<endl;
#ifdef XDEBUG
    Matrix<T> m0 = m;
    Matrix<T1> DD(D.size(),D.size(),T1(0));
    DD.diag() = D;
    DD.diag(-1) = xD;
    DD.diag(1) = herm ? xD.Conjugate() : xD.View();
    Matrix<T> m2 = DD * m;
    //cerr<<"In PseudoLMultEq:\n";
    //cerr<<"D = "<<D<<endl;
    //cerr<<"xD = "<<xD<<endl;
    //cerr<<"DD = "<<DD<<endl;
    //cerr<<"m = "<<m<<endl;
    //cerr<<"m0 = "<<m0<<endl;
    //cerr<<"DD*m = "<<m2<<endl;
#endif

    const T1* Di = D.cptr();
    const T1* xDi = xD.cptr();

    const size_t N = D.size();
    const int sd = D.step();
    const int sx = xD.step();

    for(size_t i=0;i<N;) {
      if (i==N-1 || *xDi == T1(0)) {
	if (herm) 
	  m.row(i) *= REAL(*Di);
	else
	  m.row(i) *= *Di;
	Di+=sd; xDi+=sx; ++i;
      } else {
	T1 x = *Di;
	T1 y = *(Di+=sd);
	T1 z = *xDi;
	SymLMultEq_2x2<herm>(x,y,z,m.Rows(i,i+2));
	Di+=sd,xDi+=2*sx,i+=2;
      }
    }
    //cerr<<"xD = "<<xD<<endl;
    //cerr<<"Done: m = "<<m<<endl;
#ifdef XDEBUG
    if (Norm(m2-m) > 0.00001*Norm(DD)*Norm(m0)) {
      cerr<<"PseudoDiag_LMultEq: m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"D = "<<D<<endl;
      cerr<<"xD = "<<xD<<endl;
      cerr<<"DD = "<<DD<<endl;
      cerr<<"m2 = "<<m2<<endl;
      cerr<<"-> m = "<<m<<endl;
      abort();
    }
#endif
    //cerr<<"Done PseudoDiag_LMultEq: xD = "<<xD<<endl;
  }

  // Note that we do not use the LAPACK division routines here.
  // LAPACK stores L and D in a complicated way which seems to be
  // intended to avoid the extra O(N) storage of xD.  
  // They store the xD vector in the subdiagonal of L, since the 
  // L matrix has 0's in these locations.
  //
  // However, this storage method makes the division _much_ slower,
  // because they don't permute all of the L matrix along the way.
  // This means that the permutations need to be done during the 
  // division routine, and they are mixed in with the multiplications,
  // which is a lot slower than doing all the permutations at the 
  // beginning and then dividing by a regular triangle matrix.
  //
  // Here are some timing measurements on my computer for 
  // 2000 x 2000 SymMatrix<double,Lower,ColMajor> decomposition and division,
  // both using my TMV code and directly calling the LAOPACK routines
  // (not through the TMV library).
  // All times are in seconds - I took the fastest of 3 runs.
  //
  //                      TMV Non-Lap code     Direct LAPACK calls
  //
  // Decomposition              1.4	              1.3
  // (dsytrf)
  //
  // Divide into
  // 2000 x 2000 matrix         4.9                  54.1
  // (dsytrs)
  //
  // Divide into
  // 2000 x 200 matrix          0.54                  5.4
  // (dsytrs)
  //
  // Divide into
  // 2000 x 20 matrix           0.09                  0.17
  // (dsytrs)
  //
  // Divide into
  // 2000 x 1 vector            0.018                 0.022
  // (dsytrs)
  //
  // Find Inverse               2.9                   6.3
  // (dsytri)
  //
  // Clearly, unless you are only dividing into a vector or a very thin 
  // right hand side matrix, the slight savings in the time for 
  // the LAPACK decomposition is more than lost in performing the division 
  // or calculate the inverse.  For division, a 1 to 100 size ratio on the 
  // right seems to be about the break even point.
  //
  // Hence our decision to forego the LAPACK routines here and
  // in calculating the inverse.
  //
  template <class T, class T1> void SymLU_LDivEq(
      const GenSymMatrix<T1>& LL, const GenVector<T1>& xD, const size_t* P, 
      const MatrixView<T>& m)
  {
    // Solve P L D Lt Pt x = m:
    TMVAssert(LL.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());
    TMVAssert(LL.ct() == NonConj);
    TMVAssert(xD.ct() == NonConj);

    //cerr<<"Start SymLU_LDivEq:\n";
    //cerr<<"xD = "<<xD<<endl;
#ifdef XDEBUG
    Matrix<T> m0 = m;
    Matrix<T1> DD(LL.size(),LL.size(),T1(0));
    DD.diag() = LL.diag();
    DD.diag(-1) = xD;
    DD.diag(1) = LL.isherm() ? xD.Conjugate() : xD.View();
    LowerTriMatrix<T1> L = LL.LowerTri(UnitDiag);
    Matrix<T1> LDL = L*DD*(LL.isherm() ? L.Adjoint() : L.Transpose());
    LDL.ReversePermuteRows(P);
    LDL.ReversePermuteCols(P);
#endif
    //cerr<<"xD = "<<xD<<endl;

    m.PermuteRows(P);
    m /= LL.LowerTri(UnitDiag);
    if (LL.isherm())
      PseudoDiag_LDivEq<true>(LL.diag(),xD,m);
    else
      PseudoDiag_LDivEq<false>(LL.diag(),xD,m);
    m /= LL.UpperTri(UnitDiag);
    m.ReversePermuteRows(P);

    //cerr<<"xD = "<<xD<<endl;
#ifdef XDEBUG
    Matrix<T> mm = LDL*m;
    if (Norm(mm-m0) > 0.001*Norm(LDL)*Norm(m0)) {
      cerr<<"SymLU_LDivEq: m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"L = "<<L<<endl;
      cerr<<"D = "<<LL.diag()<<endl;
      cerr<<"xD = "<<xD<<endl;
      cerr<<"DD = "<<DD<<endl;
      cerr<<"LDL = "<<LDL<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"LDL*m = "<<mm<<endl;
      abort();
    }
#endif
    //cerr<<"Done SymLU_LDivEq: xD = "<<xD<<endl;
  }

  //
  // RDivEq Matrix
  //

  template <class T, class T1> void SymLU_RDivEq(
      const GenSymMatrix<T1>& LL, const GenVector<T1>& xD, const size_t* P, 
      const MatrixView<T>& m)
  {
    // Solve x P L D Lt Pt = m:
    // P L Dt Lt Pt xt = mt
    TMVAssert(LL.size() == m.rowsize());
    TMVAssert(xD.size()+1 == m.rowsize());
    TMVAssert(LL.ct() == NonConj);
    TMVAssert(xD.ct() == NonConj);
    SymLU_LDivEq(LL,xD,P,LL.isherm()?m.Adjoint():m.Transpose());
  }

#define InstFile "TMV_SymLUDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


