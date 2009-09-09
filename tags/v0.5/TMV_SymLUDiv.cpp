
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

#define APTR (inplace ? A.NonConst().ptr() : new T[A.size()*A.size()])
#define LLX \
  (inplace ? A.uplo()==Upper ? A.NonConst().Adjoint() : A.NonConst() : \
  HermMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.LowerTri())))
#define XDPTR (inplace ? new T[A.size()-1] : Aptr+LLx.stepj())
#define XD VectorViewOf(xDptr,A.size()-1,inplace ? 1 : int(A.size()+1))

  template <class T> HermLUDiv<T>::HermLUDiv(const GenSymMatrix<T>& A,
      bool _inplace) :
    inplace(_inplace), Aptr(APTR), LLx(LLX), xDptr(XDPTR), xD(XD),
    P(new size_t[A.colsize()]), det(T(1))
  {
    TMVAssert(IsReal(T()) || A.isherm());
    if (inplace) { TMVAssert(A == LLx); }
    else LLx = A;

    SymLU_Decompose(LLx,xD,P,det);
  }

#undef LLX
#define LLX \
  inplace ? A.uplo()==Upper ? A.NonConst().Transpose() : A.NonConst() : \
  SymMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.LowerTri()))

  template <class T> SymLUDiv<T>::SymLUDiv(const GenSymMatrix<T>& A,
      bool _inplace) :
    inplace(_inplace), Aptr(APTR), LLx(LLX), xDptr(XDPTR), xD(XD),
    P(new size_t[A.colsize()]), det(T(1))
  {
    TMVAssert(IsComplex(T()) && !A.isherm());
    if (inplace) { TMVAssert(A == LLx); }
    else LLx = A;

    SymLU_Decompose(LLx,xD,P,det);
  }

#undef APTR
#undef LLX
#undef XDPTR
#undef XD

  template <class T> HermLUDiv<T>::~HermLUDiv() 
  { delete[] P; if (!inplace) delete[] Aptr; else delete[] xDptr; }

  template <class T> SymLUDiv<T>::~SymLUDiv() 
  { delete[] P; if (!inplace) delete[] Aptr; else delete[] xDptr; }

  template <class T> bool SymLUDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenSymMatrix<T>* sm = dynamic_cast<const GenSymMatrix<T>*>(&m);
    TMVAssert(sm);
    if (fout) {
      *fout << "M = "<<*sm<<endl;
      *fout << "L = "<<GetL()<<endl;
      *fout << "D = "<<GetD()<<endl;
    }
    Matrix<T> m2 = GetL()*GetD()*GetL().Transpose();
    m2.ReversePermuteRows(P);
    m2.ReversePermuteCols(P);
    RealType(T) nm = Norm(m2-*sm);
    nm /= SQR(Norm(GetL()))*Norm(GetD());
    if (fout) {
      *fout << "LDLT = "<<m2<<endl;
      *fout << "Norm(M-LDLT)/Norm(LDLT) = "<<nm<<endl;
    }
    return nm < Epsilon<T>();
  }

  template <class T> bool HermLUDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenSymMatrix<T>* sm = dynamic_cast<const GenSymMatrix<T>*>(&m);
    TMVAssert(sm);
    if (fout) {
      *fout << "M = "<<*sm<<endl;
      *fout << "L = "<<GetL()<<endl;
      *fout << "D = "<<GetD()<<endl;
    }
    Matrix<T> m2 = GetL()*GetD()*GetL().Adjoint();
    m2.ReversePermuteRows(P);
    m2.ReversePermuteCols(P);
    RealType(T) nm = Norm(m2-*sm);
    nm /= SQR(Norm(GetL()))*Norm(GetD());
    if (fout) {
      *fout << "LDLt = "<<m2<<endl;
      *fout << "Norm(M-LDLt)/Norm(LDLt) = "<<nm<<endl;
    }
    return nm < Epsilon<T>();
  }

  //
  // LDivEq
  //

  template <class T> void SymInvert_2x2(T& a, T& b, T& c, T* dd=0)
  {
    // Invert matrix [ a  c ] -->  1/(ab-c^2) [ b -c ]
    //               [ c  b ]                 [ -c a ]
    T d = a*b-c*c;
    swap(a,b);
    a/=d;
    b/=d;
    c = -c/d;
    if (dd) *dd = d;
  }

  template <class T> void HermInvert_2x2(
      RealType(T)& a, RealType(T)& b, T& c, RealType(T)* dd=0)
  {
    // Invert matrix [ a  c* ] -->  1/(ab-|c|^2) [ b -c* ]
    //               [ c  b  ]                   [ -c a  ]
    RealType(T) d = a*b-NORM(c);
    swap(a,b);
    a/=d;
    b/=d;
    c = -c/d;
    if (dd) *dd = d;
  }

  template <class T1ab, class T1c, class T2> inline void LMultEq_2x2(
      T1ab a, T1ab b, T1c c, T1c cc, T2& x, T2& y)
  {
    // Solve [ x ] <- [ a  cc ] [ x ]
    //       [ y ]    [ c  b  ] [ y ]
    T2 tempx = x;
    x = (a*x+cc*y);
    y = (b*y+c*tempx);
  }

  template <class T1ab, class T1c, class T2> void LMultEq_2x2(
      T1ab a, T1ab b, T1c c, T1c cc, const MatrixView<T2>& m)
  {
    // Solve m <- [ a  cc ] m 
    //            [ c  b  ]
    TMVAssert(m.colsize() == 2);
    TMVAssert(m.ct() == NonConj);
    if (m.isrm()) {
      VIt<T2,Unit,NonConj> m0 = m.row(0).begin();
      VIt<T2,Unit,NonConj> m1 = m.row(1).begin();
      for (size_t k=m.rowsize();k>0;--k,++m0,++m1) 
	LMultEq_2x2(a,b,c,cc,*m0,*m1);
    } else {
      VIt<T2,Step,NonConj> m0 = m.row(0).begin();
      VIt<T2,Step,NonConj> m1 = m.row(1).begin();
      for (size_t k=m.rowsize();k>0;--k,++m0,++m1) 
	LMultEq_2x2(a,b,c,cc,*m0,*m1);
    }
  }

  template <class T1, class T2> inline void SymLMultEq_2x2(
      T1 a, T1 b, T1 c, const MatrixView<T2>& m)
  {
    TMVAssert(m.colsize() == 2);
    if (m.isconj())
      LMultEq_2x2(CONJ(a),CONJ(b),CONJ(c),CONJ(c),m.Conjugate());
    else
      LMultEq_2x2(a,b,c,c,m);
  }

  template <class T1, class T2> inline void SymRMultEq_2x2(
      T1 a, T1 b, T1 c, const MatrixView<T2>& m)
  {
    TMVAssert(m.rowsize() == 2);
    if (m.isconj())
      LMultEq_2x2(CONJ(a),CONJ(b),CONJ(c),CONJ(c),m.Adjoint());
    else
      LMultEq_2x2(a,b,c,c,m.Transpose());
  }

  template <class T1, class T2> inline void HermLMultEq_2x2(
      RealType(T1) a, RealType(T1) b, T1 c, const MatrixView<T2>& m)
  {
    TMVAssert(m.colsize() == 2);
    if (m.isconj())
      LMultEq_2x2(a,b,CONJ(c),c,m.Conjugate());
    else
      LMultEq_2x2(a,b,c,CONJ(c),m);
  }

  template <class T1, class T2> inline void HermRMultEq_2x2(
      RealType(T1) a, RealType(T1) b, T1 c, const MatrixView<T2>& m)
  {
    TMVAssert(m.rowsize() == 2);
    if (m.isconj())
      LMultEq_2x2(a,b,c,CONJ(c),m.Adjoint());
    else
      LMultEq_2x2(a,b,CONJ(c),c,m.Transpose());
  }

  template <class T1, class T2> void SymPseudoDiag_LDivEq(
      const GenVector<T1>& D, const GenVector<T1>& xD,
      const MatrixView<T2>& m)
  {
    TMVAssert(D.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());

    const size_t N = D.size();
    CVIter<T1> Di = D.begin();
    CVIter<T1> xDi = xD.begin();
    for(size_t i=0;i<N;) {
      if (i==N-1 || *xDi == T1(0)) {
	m.row(i) /= *Di;
	++Di; ++xDi; ++i;
      } else {
	T1 x = *Di;
	T1 y = *(++Di);
	T1 z = *xDi;
	SymInvert_2x2(x,y,z);
	SymLMultEq_2x2(x,y,z,m.Rows(i,i+2));
	++Di,xDi+=2,i+=2;
      }
    }

  }

  template <class T1, class T2> void HermPseudoDiag_LDivEq(
      const GenVector<RealType(T1)>& D, const GenVector<T1>& xD,
      const MatrixView<T2>& m)
  {
    TMVAssert(D.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());

    const size_t N = D.size();
    CVIter<RealType(T1)> Di = D.begin();
    CVIter<T1> xDi = xD.begin();
    for(size_t i=0;i<N;) {
      if (i==N-1 || *xDi == T1(0)) {
	m.row(i) /= *Di;
	++Di; ++xDi; ++i;
      } else {
	RealType(T1) x = *Di;
	RealType(T1) y = *(++Di);
	T1 z = *xDi;
	HermInvert_2x2(x,y,z);
	HermLMultEq_2x2(x,y,z,m.Rows(i,i+2));
	++Di,xDi+=2,i+=2;
      }
    }
  }

  template <class T1, class T2> void SymPseudoDiag_LMultEq(
      const GenVector<T1>& D, const GenVector<T1>& xD,
      const MatrixView<T2>& m)
  {
    TMVAssert(D.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());

    const size_t N = D.size();
    CVIter<T1> Di = D.begin();
    CVIter<T1> xDi = xD.begin();
    for(size_t i=0;i<N;) {
      if (i==N-1 || *xDi == T1(0)) {
	m.row(i) *= *Di;
	++Di; ++xDi; ++i;
      } else {
	T1 x = *Di;
	T1 y = *(++Di);
	T1 z = *xDi;
	SymLMultEq_2x2(x,y,z,m.Rows(i,i+2));
	++Di,xDi+=2,i+=2;
      }
    }
  }

  template <class T1, class T2> void HermPseudoDiag_LMultEq(
      const GenVector<RealType(T1)>& D, const GenVector<T1>& xD,
      const MatrixView<T2>& m)
  {
    TMVAssert(D.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());

    const size_t N = D.size();
    CVIter<RealType(T1)> Di = D.begin();
    CVIter<T1> xDi = xD.begin();
    for(size_t i=0;i<N;) {
      if (i==N-1 || *xDi == T1(0)) {
	m.row(i) *= *Di;
	++Di; ++xDi; ++i;
      } else {
	RealType(T1) x = *Di;
	RealType(T1) y = *(++Di);
	T1 z = *xDi;
	HermLMultEq_2x2(x,y,z,m.Rows(i,i+2));
	++Di,xDi+=2,i+=2;
      }
    }
  }

  template <class T1, class T2> void SymLU_LDivEq(
      const GenSymMatrix<T1>& LL, const GenVector<T1>& xD, const size_t* P, 
      const MatrixView<T2>& m)
  {
    // Solve P L D LT Pt x = m:
    TMVAssert(IsComplex(T1()) && !LL.isherm());
    TMVAssert(LL.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());

#ifdef XDEBUG
    //cerr<<"SymLU_LDivEq\n";
    Matrix<T1> DD(LL.size(),LL.size());
    DD.SetToIdentity();
    SymPseudoDiag_LMultEq(LL.diag(),xD,DD.View());
    LowerTriMatrix<T1> L = LL.LowerTri().MakeUnitDiag();
    Matrix<T1> LDL = L*DD*L.Transpose();
    LDL.ReversePermuteRows(P);
    LDL.ReversePermuteCols(P);
    Matrix<T2> m0 = m;
#endif

    m.PermuteRows(P);
    m /= LL.LowerTri().MakeUnitDiag();
    SymPseudoDiag_LDivEq(LL.diag(),xD,m);
    m /= LL.UpperTri().MakeUnitDiag();
    m.ReversePermuteRows(P);

#ifdef XDEBUG
    Matrix<T2> mm = LDL*m;
    if (Norm(mm-m0) > 0.001*Norm(m0)) {
      cerr<<"SymPseudoDiag_LDivEq: m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"LDL = "<<LDL<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"LDL*m = "<<mm<<endl;
      abort();
    }
#endif
  }

  template <class T1, class T2> void HermLU_LDivEq(
      const GenSymMatrix<T1>& LL, const GenVector<T1>& xD, const size_t* P, 
      const MatrixView<T2>& m)
  {
    // Solve P L D Lt Pt x = m:
    TMVAssert(IsReal(T1()) || LL.isherm());
    TMVAssert(LL.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());

#ifdef XDEBUG
    //cerr<<"HermLU_LDivEq\n";
    //cerr<<"LL = "<<Type(LL)<<"  "<<LL<<endl;
    //cerr<<"LL.LowerTri = "<<LL.LowerTri()<<endl;
    //cerr<<"L = "<<LL.LowerTri().MakeUnitDiag()<<endl;
    Matrix<T2> m0 = m;
    //cerr<<"m = "<<Type(m)<<"  "<<m0<<endl;
    Matrix<T1> DD(LL.size(),LL.size());
    DD.SetToIdentity();
    HermPseudoDiag_LMultEq(LL.diag().Real(),xD,DD.View());
    LowerTriMatrix<T1> L = LL.LowerTri().MakeUnitDiag();
    Matrix<T1> LDL = L*DD*L.Adjoint();
    LDL.ReversePermuteRows(P);
    LDL.ReversePermuteCols(P);
    //cerr<<"LDL = "<<LDL<<endl;
#endif

    m.PermuteRows(P);
    m /= LL.LowerTri().MakeUnitDiag();
    HermPseudoDiag_LDivEq(LL.diag().Real(),xD,m);
    m /= LL.UpperTri().MakeUnitDiag();
    m.ReversePermuteRows(P);

#ifdef XDEBUG
    Matrix<T2> mm = LDL*m;
    //cerr<<"-> m = "<<m<<endl;
    //cerr<<"LDL*m = "<<mm<<endl;
    if (Norm(mm-m0) > 0.001*Norm(m0)) {
      cerr<<"HermPseudoDiag_LDivEq: m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"LDL = "<<LDL<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"LDL*m = "<<mm<<endl;
      abort();
    }
#endif
  }

  //
  // RDivEq Matrix
  //

  template <class T1, class T2> inline void SymLU_RDivEq(
      const GenSymMatrix<T1>& LL, const GenVector<T1>& xD, const size_t* P, 
      const MatrixView<T2>& m)
  {
    // Solve x P L D LT Pt = m:
    TMVAssert(LL.size() == m.rowsize());
    TMVAssert(xD.size()+1 == m.rowsize());
    SymLU_LDivEq(LL,xD,P,m.Transpose());
  }

  template <class T1, class T2> inline void HermLU_RDivEq(
      const GenSymMatrix<T1>& LL, const GenVector<T1>& xD, const size_t* P, 
      const MatrixView<T2>& m)
  {
    // Solve x P L D Lt Pt = m:
    TMVAssert(LL.size() == m.rowsize());
    TMVAssert(xD.size()+1 == m.rowsize());
    HermLU_LDivEq(LL,xD,P,m.Adjoint());
  }

  template <class T> void HermLUDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LLx.size());
    TMVAssert(minv.rowsize() == LLx.size());

    minv.SetToIdentity();
    LDivEq(minv.View());
  }

  template <class T> void HermLUDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LLx.size());
    TMVAssert(minv.rowsize() == LLx.size());

    Matrix<T,ColMajor> temp = Eye<T,ColMajor>(LLx.size());
    LDivEq(temp.View());
    minv =  temp*temp.Adjoint();
  }

  template <class T> void SymLUDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LLx.size());
    TMVAssert(minv.rowsize() == LLx.size());

    minv.SetToIdentity();
    LDivEq(minv.View());
  }

  template <class T> void SymLUDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LLx.size());
    TMVAssert(minv.rowsize() == LLx.size());

    Matrix<T,ColMajor> temp = Eye<T,ColMajor>(LLx.size());
    LDivEq(temp.View());
    minv =  temp*temp.Adjoint();
  }

  template <class T> const Matrix<T> HermLUDiv<T>::GetD() const 
  { 
    Matrix<T> temp(LLx.size(),LLx.size());
    temp.SetToIdentity();
    HermPseudoDiag_LMultEq(LLx.diag().Real(),xD,temp.View());
    return temp;
  }

  template <class T> const Matrix<T> SymLUDiv<T>::GetD() const 
  { 
    Matrix<T> temp(LLx.size(),LLx.size());
    temp.SetToIdentity();
    SymPseudoDiag_LMultEq(LLx.diag(),xD,temp.View());
    return temp;
  }

  //
  // Decompose
  //

  // MJ: Write Level 3 version
  template <class T> void NonLapHermLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, T& det)
  {
    // Bunch-Kauffman algorithm for LDL Decomposition of a symmetric matrix.
    //
    // We want to decompose the matrix (input as A) into P * L * D * Lt * Pt,
    // where P is a permutation, L is a lower triangle matrix with 1's for 
    // the diagonal, and D is a pseudo-diagonal matrix, meaning that 
    // there may be some 2x2 blocks along the diagonal.
    //
    // The reason LU decomposition is hard for symmetric matrices lies in
    // the pivoting.  If we pivot an off-diagonal element onto the diagonal,
    // then it destroys the symmetric structure of the matrix.
    // We can symmetrically pivot to move diagonal elements up or down the 
    // diagonal, but if a diagonal element is 0 (or very small compared to
    // an off-diagonal element) dividing by it causes problems.
    //
    // Bunch and Parlett showed that we can pivot off-diagonal elements 
    // to the first off diagonal and leave the 2x2 block unreduced in the 
    // D matrix.
    // When we are done, dividing by the 2x2 blocks is trivial to do
    // with correct pivoting.
    // 
    // (Following the explanation of Golub and van Loan, section 4.4.4:)
    // Consider the first step in the decomposition where we find
    // P1, E, C, B such that
    //
    // P1 A P1t = [ E Ct ]
    //            [ C B  ]
    // where E is sxs, C is (n-s)xs, B is (n-s)x(n-s), and s is either 1 or 2.
    // 
    // If A is non-zero, then E can be chosed to be non-singular, so:
    //
    // [ E Ct ] = [   I   0 ] [ E     0     ] [ I E^-1Ct ]
    // [ C B  ]   [ CE^-1 I ] [ 0 B-CE^-1Ct ] [ 0   I    ]
    //
    // CE^-1 is the first portion of the L matrix will will be making.
    // E is the first part of D.
    // B-CE^-1Ct = A~ is the submatrix on which to continue the process.
    //
    // Bunch and Parlett showed how to find the appropriate P,E to 
    // minimize the growth of A~.
    // Define mu0 = max_i,j |A(i,j)|
    //        mu1 = max_i |A(i,i)|
    //        
    // When the largest element is on the diagonal (mu0 = mu1), then 
    // it seems obvious that we would want to use it for E.
    // When the largest element is off-diagonal and it is much larger than
    // any on-diagonal element mu0 >> mu1, we want E to be 2x2,
    // with mu0 being the off-diagonal of E.
    //
    // When mu0 is only slightly larger than mu1, it is less clear which
    // strategy is better.
    //
    // The Bunch-Parlett strategy is to take:
    // s=1, E=mu1     when mu1 > alpha * mu0.
    // s=2, E_10=mu0  when mu1 < alpha * mu0.
    // where alpha is some parameter in [0,1].
    //
    // To determine what is a good value for alpha, find the bound on
    // the values of A~ in each case:
    // s=1: |a~_ij| < max(|B-CE^-1Ct|) < max(B) + max(CCt/mu1) 
    //              < mu0 + max(CCt)/(alpha*mu0) 
    //              < mu0 + mu0^2/(alpha*mu0)
    //              = mu0 ( 1 + 1/alpha )
    // s=2: |a~_ij| < max(|B-CE^-1Ct|) < max(B) + max(CE^-1Ct)
    //              < mu0 + max(|[ mu0 mu0 ] [ x  mu0 ]^-1 [ mu0 ]|)
    //                                       [ mu0 y  ]    [ mu0 ]
    //              [ max over (x,y) with -alpha*mu0 < x,y < alpha*mu0 ]
    //              = mu0 + max(|[ mu0 mu0 ] [ y  -mu0 ] [ mu0 ] / (xy-mu0^2)|)
    //                                       [ -mu0  x ] [ mu0 ] 
    //              = mu0 + max( |x+y-2mu0| mu0^2 / |xy-mu0^2| )
    //              [ Numerator is largest when x+y = -2alpha*mu0.        ]
    //              [ Denominator is smallest when xy = alpha^2 mu0^2.    ]
    //              [ These are consistent with each other, so that's it. ]
    //              = mu0 + 2(1+alpha)mu0/(1-alpha^2)
    //              = mu0 ( 1-alpha^2+2+2*alpha ) / (1-alpha^2)
    //              = mu0 ( 3 + 2alpha - alpha^2) / (1-alpha^2)
    //              = mu0 ( 3 - alpha ) / ( 1 - alpha )
    //
    // The optimal alpha is where the growth from 2 s=1 steps equals the 
    // growth from a single s=2 step:
    //
    // (1+1/a)^2 = (3-a)/(1-a)
    // 1+a-a^2-a^3 = 3a^2 - a^3
    // 1+a-4a^2 = 0
    // alpha = (1+sqrt(17))/8
    //
    // The only trouble with this algorithm is that finding mu0 requires
    // O(n^2) comparisons, which needs to be done O(n) times, so these become 
    // a significant fraction of the computation time.
    //
    // Bunch and Kauffman modified the algorithm slightly to require only
    // O(n) comparisons for each step.  The values to calculate are:
    // a00 = |A(0,0)|
    // ap0 = |A(p,0)| = max_i |A(i,0)| 
    // apq = |A(p,q)| = max_j!=p |A(p,j)|
    // app = |A(p,p)|
    //
    // Then their tests are:
    //
    // if a00 > alpha * ap0
    //   s = 1, E = a00  (No need to calculate arp.)
    // else if a00*apq > alpha * ap0^2
    //   s = 1, E = a00
    // else if app > alpha * apq
    //   s = 1, E = app
    // else 
    //   s = 2, E = (a00,app,ap0)
    // [Note: Golub and van Loan wrongly say to put apq in E, rather than ap0.]
    // 
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    const size_t N = A.size();
    static const RealType(T) alpha = (SQRT(RealType(T)(17))+1)/8;

    xD.Zero();
    VectorView<RealType(T)> D = A.diag().Real();
#ifdef XDEBUG
    Matrix<T> L(A.size(),A.size());
    L.SetToIdentity();
    Matrix<T> DD = A;
    Matrix<T> A0 = A;
    Matrix<T> A2 = L*DD*L.Adjoint();
    for(size_t i=0;i<A.size();i++) P[i]=i;
    A2.ReversePermuteRows(P);
    A2.ReversePermuteCols(P);
    TMVAssert(Norm(A0-A2) <= 0.001*Norm(A0));
#endif

    bool seq1;
    for (size_t j=0; j<N;) // ++j or j+=2 done below
    {
      seq1 = true;
      RealType(T) ajj = abs(D(j));
      size_t p; // p is relative to j index, not absolute.
      RealType(T) apj = A.col(j,j,N).MaxAbsElement(&p);

      if (p == 0 || ajj >= alpha * apj) {
	// No permutation
	P[j] = j;
      } else {
	p+=j;
	RealType(T) app = abs(D(p));
	RealType(T) apq = A.row(p,j,p).MaxAbsElement();
	if (p+1 < N) {
	  RealType(T) apq2 = A.col(p,p+1,N).MaxAbsElement();
	  apq = max(apq,apq2);
	}
	if (ajj*apq >= alpha * apj * apj) {
	  // No permutation
	  P[j] = j;
	} else if (app >= alpha * apq) {
	  // Permute p diagonal into j spot
	  A.SwapRowsCols(j,p);
	  P[j] = p;
	  ajj = app;
#ifdef XDEBUG
	  L.SubMatrix(0,N,0,j).SwapRows(j,p);
#endif
	} else {
	  // Permute pj element into j+1,j spot
	  // This also permutes pp element into j+1,j+1
	  seq1 = false;
	  P[j] = j;
	  P[j+1] = p;
	  if (p != j+1) {
	    A.SwapRowsCols(j+1,p);
#ifdef XDEBUG
	    L.SubMatrix(0,N,0,j).SwapRows(j+1,p);
#endif
	  }
	}
      }

      // Now the LU solving:
      if (seq1) {
	if (ajj == T(0))
	  tmv_error("Zero pivot found in HermLU_Decompose");
	det *= D(j);
	A.col(j,j+1,N) /= D(j);
	// A.SubSymMatrix(j+1,N) -= 
	//     D(j)*A.col(j,j+1,N)^A.col(j,j+1,N).Conjugate()
	Rank1Update(T(-D(j)),A.col(j,j+1,N),A.SubSymMatrix(j+1,N));
#ifdef XDEBUG
	DD.col(j,j+1,N).Zero();
	DD.row(j,j+1,N).Zero();
	DD(j,j) = D(j);
	L.col(j,j+1,N) = A.col(j,j+1,N);
	DD.SubMatrix(j+1,N,j+1,N) = A.SubSymMatrix(j+1,N);
        A2 = L*DD*L.Adjoint();
	A2.ReversePermuteRows(P);
	A2.ReversePermuteCols(P);
	if (Norm(A0-A2) > 0.001*Norm(A0)) {
	  cerr<<"Herm: s==1\n";
	  cerr<<"j = "<<j<<endl;
	  cerr<<"A = "<<A<<endl;
	  cerr<<"D = "<<D<<endl;
	  cerr<<"xD = "<<xD<<endl;
	  cerr<<"L = "<<L<<endl;
	  cerr<<"DD = "<<DD<<endl;
	  cerr<<"A2 = "<<A2<<endl;
	  cerr<<"A0 = "<<A0<<endl;
	  cerr<<"Norm(A0-A2) = "<<Norm(A0-A2)<<endl;
	  abort();
	}
#endif
	++j;
      } else {
	// Invert E:  E^-1 = [ x z* ]
	//                   [ z y  ]
	RealType(T) x = D(j);
	RealType(T) y = D(j+1);
	T z = xD(j) = A(j+1,j);
	A(j+1,j)=T(0);
	RealType(T) d;
	HermInvert_2x2(x,y,z,&d);
	det *= d;

	if (A.isrm()) {
	  // A(j+2:N,j:j+2) = CE^-1
	  // A(j+2:N,j+2:N) -= CE^-1Ct
	  // A(p,q) -= CEinv(p,0)*(C(q,0)*) + CEinv(p,1)*(C(q,1)*)
	  //   p >= q, so loop from bottom up.
	  size_t len = N-j-2;
	  VIt<T,Step,NonConj> Cp0=A.col(j,j+2,N).begin()+len-1;
	  VIt<T,Step,NonConj> Cp1=A.col(j+1,j+2,N).begin()+len-1;
	  for(size_t p = N-1; p>=j+2; --p,--Cp0,--Cp1,--len) {
	    T CEinvp0 = *Cp0 * x + *Cp1 * z;
	    T CEinvp1 = *Cp0 * CONJ(z) + *Cp1 * y;
	    A.row(p,j+2,p+1) -= CEinvp0 * A.col(j,j+2,p+1).Conjugate();
	    A.row(p,j+2,p+1) -= CEinvp1 * A.col(j+1,j+2,p+1).Conjugate();
	    *Cp0 = CEinvp0;
	    *Cp1 = CEinvp1;
	  }
	} else {
	  // A(j+2:N,j:j+2) = CE^-1
	  // A(j+2:N,j+2:N) -= CE^-1Ct 
	  //                -= C (CE^-1)t (since E^-1 is Hermitian)
	  // A(p,q) -= C(p,0)*(CEinv(q,0)*) + C(p,1)*(CEinv(q,1)*)
	  //   p >= q, so loop from top down
	  size_t len = 0;
	  VIt<T,Unit,NonConj> Cq0=A.col(j,j+2,N).begin();
	  VIt<T,Unit,NonConj> Cq1=A.col(j+1,j+2,N).begin();
	  for(size_t q = j+2; q<N; ++q,++Cq0,++Cq1,++len) {
	    T CEinvq0 = *Cq0 * x + *Cq1 * z;
	    T CEinvq1 = *Cq0 * CONJ(z) + *Cq1 * y;
	    A.col(q,q,N) -= CONJ(CEinvq0) * A.col(j,q,N);
	    A.col(q,q,N) -= CONJ(CEinvq1) * A.col(j+1,q,N);
	    *Cq0 = CEinvq0;
	    *Cq1 = CEinvq1;
	  }
	}
#ifdef XDEBUG
	DD.col(j,j+2,N).Zero();
	DD.col(j+1,j+2,N).Zero();
	DD.row(j,j+2,N).Zero();
	DD.row(j+1,j+2,N).Zero();
	DD(j,j) = D(j);
	DD(j+1,j+1) = D(j+1);
	DD(j+1,j) = xD(j);
	DD(j,j+1) = CONJ(xD(j));
	L.col(j,j+2,N) = A.col(j,j+2,N);
	L.col(j+1,j+2,N) = A.col(j+1,j+2,N);
	DD.SubMatrix(j+2,N,j+2,N) = A.SubSymMatrix(j+2,N);
        A2 = L*DD*L.Adjoint();
	A2.ReversePermuteRows(P);
	A2.ReversePermuteCols(P);
	if (Norm(A0-A2) > 0.001*Norm(A0)) {
	  cerr<<"Herm: s==2\n";
	  cerr<<"j = "<<j<<endl;
	  cerr<<"A = "<<A<<endl;
	  cerr<<"D = "<<D<<endl;
	  cerr<<"xD = "<<xD<<endl;
	  cerr<<"L = "<<L<<endl;
	  cerr<<"DD = "<<DD<<endl;
	  cerr<<"A2 = "<<A2<<endl;
	  cerr<<"A0 = "<<A0<<endl;
	  cerr<<"Norm(A0-A2) = "<<Norm(A0-A2)<<endl;
	  abort();
	}
#endif
	j+=2;
      }
    }
  }

  template <class T> void NonLapSymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, T& det)
  {
    // Same thing, but for Symmetric (not Hermitian) Matrix
    TMVAssert(IsComplex(T()));
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);

    //cerr<<"Start Symm LDL Decompose: A = "<<Type(A)<<"  "<<A<<endl;

    const size_t N = A.size();
    static const RealType(T) alpha = (SQRT(RealType(T)(17))+1)/8;

    xD.Zero();
    VectorView<T> D = A.diag();
#ifdef XDEBUG
    Matrix<T> L(A.size(),A.size());
    L.SetToIdentity();
    Matrix<T> DD = A;
    Matrix<T> A0 = A;
    Matrix<T> A2 = L*DD*L.Transpose();
    for(size_t i=0;i<A.size();i++) P[i]=i;
    A2.ReversePermuteRows(P);
    A2.ReversePermuteCols(P);
    TMVAssert(Norm(A0-A2) <= 0.001*Norm(A0));
#endif

    bool seq1;
    for (size_t j=0; j<N;) // ++j or j+=2 done below
    {
      seq1 = true;
      RealType(T) ajj = abs(D(j));
      size_t p; // p is relative to j index, not absolute.
      RealType(T) apj = A.col(j,j,N).MaxAbsElement(&p);

      if (p == 0 || ajj >= alpha * apj) {
	// No permutation
	P[j] = j;
      } else {
	p+=j;
	RealType(T) app = abs(D(p));
        RealType(T) apq = A.row(p,j,p).MaxAbsElement();
	if (p+1 < N) {
	  RealType(T) apq2 = A.col(p,p+1,N).MaxAbsElement();
	  apq = max(apq,apq2);
	}
	if (ajj*apq >= alpha * apj * apj) {
	  // No permutation
	  P[j] = j;
	} else if (app >= alpha * apq) {
	  // Permute p diagonal into j spot
	  A.SwapRowsCols(j,p);
	  P[j] = p;
	  ajj = app;
#ifdef XDEBUG
	  L.SubMatrix(0,N,0,j).SwapRows(j,p);
#endif
	} else {
	  // Permute pj element into j+1,j spot
	  // This also permutes pp element into j+1,j+1
	  seq1 = false;
	  P[j] = j;
	  P[j+1] = p;
	  if (p != j+1) {
	    A.SwapRowsCols(j+1,p);
#ifdef XDEBUG
	    L.SubMatrix(0,N,0,j).SwapRows(j+1,p);
#endif
	  }
	}
      }

      // Now the LU solving:
      if (seq1) {
	if (ajj == T(0)) 
	  tmv_error("Zero pivot found in SymLU_Decompose");
	det *= D(j);
	A.col(j,j+1,N) /= D(j);
	// A.SubSymMatrix(j+1,N) -= D(j)*A.col(j,j+1,N)^A.col(j,j+1,N)
	Rank1Update(-D(j),A.col(j,j+1,N),A.SubSymMatrix(j+1,N));
#ifdef XDEBUG
	DD.col(j,j+1,N).Zero();
	DD.row(j,j+1,N).Zero();
	DD(j,j) = D(j);
	L.col(j,j+1,N) = A.col(j,j+1,N);
	DD.SubMatrix(j+1,N,j+1,N) = A.SubSymMatrix(j+1,N);
        A2 = L*DD*L.Transpose();
	A2.ReversePermuteRows(P);
	A2.ReversePermuteCols(P);
	if (Norm(A0-A2) > 0.001*Norm(A0)) {
	  cerr<<"Sym: s==1\n";
	  cerr<<"j = "<<j<<endl;
	  cerr<<"A = "<<A<<endl;
	  cerr<<"D = "<<D<<endl;
	  cerr<<"xD = "<<xD<<endl;
	  cerr<<"L = "<<L<<endl;
	  cerr<<"DD = "<<DD<<endl;
	  cerr<<"A2 = "<<A2<<endl;
	  cerr<<"A0 = "<<A0<<endl;
	  cerr<<"Norm(A0-A2) = "<<Norm(A0-A2)<<endl;
	  abort();
	}
#endif
	++j;
      } else {
	// Invert E:  E^-1 = [ x z ]
	//                   [ z y ]
	T x = D(j);
	T y = D(j+1);
	T z = xD(j) = A(j+1,j);
	A(j+1,j) = T(0);
	T d;
	SymInvert_2x2(x,y,z,&d);
	det *= d;

	if (A.isrm()) {
	  // A(j+2:N,j:j+2) = CE^-1
	  // A(j+2:N,j+2:N) -= CE^-1CT
	  // A(p,q) -= CEinv(p,0)*C(q,0) + CEinv(p,1)*C(q,1)
	  //   p >= q, so loop from bottom up.
	  size_t len = N-j-2;
	  VIt<T,Step,NonConj> Cp0=A.col(j,j+2,N).begin()+len-1;
	  VIt<T,Step,NonConj> Cp1=A.col(j+1,j+2,N).begin()+len-1;
	  for(size_t p = N-1; p>=j+2; --p,--Cp0,--Cp1,--len) {
	    T CEinvp0 = *Cp0 * x + *Cp1 * z;
	    T CEinvp1 = *Cp0 * z + *Cp1 * y;
	    A.row(p,j+2,p+1) -= CEinvp0 * A.col(j,j+2,p+1);
	    A.row(p,j+2,p+1) -= CEinvp1 * A.col(j+1,j+2,p+1);
	    *Cp0 = CEinvp0;
	    *Cp1 = CEinvp1;
	  }
	} else {
	  // A(j+2:N,j:j+2) = CE^-1
	  // A(j+2:N,j+2:N) -= CE^-1CT 
	  //                -= C (CE^-1)T (since E^-1 is Symmetric)
	  // A(p,q) -= C(p,0)*CEinv(q,0) + C(p,1)*CEinv(q,1)
	  //   p >= q, so loop from top down
	  size_t len = 0;
	  VIt<T,Unit,NonConj> Cq0=A.col(j,j+2,N).begin();
	  VIt<T,Unit,NonConj> Cq1=A.col(j+1,j+2,N).begin();
	  for(size_t q = j+2; q<N; ++q,++Cq0,++Cq1,++len) {
	    T CEinvq0 = *Cq0 * x + *Cq1 * z;
	    T CEinvq1 = *Cq0 * z + *Cq1 * y;
	    A.col(q,q,N) -= CEinvq0 * A.col(j,q,N);
	    A.col(q,q,N) -= CEinvq1 * A.col(j+1,q,N);
	    *Cq0 = CEinvq0;
	    *Cq1 = CEinvq1;
	  }
	}
#ifdef XDEBUG
	DD.col(j,j+2,N).Zero();
	DD.col(j+1,j+2,N).Zero();
	DD.row(j,j+2,N).Zero();
	DD.row(j+1,j+2,N).Zero();
	DD(j,j) = D(j);
	DD(j+1,j+1) = D(j+1);
	DD(j+1,j) = xD(j);
	DD(j,j+1) = xD(j);
	L.col(j,j+2,N) = A.col(j,j+2,N);
	L.col(j+1,j+2,N) = A.col(j+1,j+2,N);
	DD.SubMatrix(j+2,N,j+2,N) = A.SubSymMatrix(j+2,N);
        A2 = L*DD*L.Transpose();
	A2.ReversePermuteRows(P);
	A2.ReversePermuteCols(P);
	if (Norm(A0-A2) > 0.001*Norm(A0)) {
	  cerr<<"Sym: s==2\n";
	  cerr<<"j = "<<j<<endl;
	  cerr<<"A = "<<A<<endl;
	  cerr<<"D = "<<D<<endl;
	  cerr<<"xD = "<<xD<<endl;
	  cerr<<"L = "<<L<<endl;
	  cerr<<"DD = "<<DD<<endl;
	  cerr<<"A2 = "<<A2<<endl;
	  cerr<<"A0 = "<<A0<<endl;
	  cerr<<"Norm(A0-A2) = "<<Norm(A0-A2)<<endl;
	  abort();
	}
#endif
	j+=2;
      }
    }
  }

#ifdef LAP
  template <class T> inline void LapSymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD,
      size_t* P, T& det)
  { NonLapSymLU_Decompose(A,xD,P,det); }
  template <class T> inline void LapHermLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD,
      size_t* P, T& det)
  { NonLapHermLU_Decompose(A,xD,P,det); }
  template <> inline void LapHermLU_Decompose(
      const SymMatrixView<double>& A, const VectorView<double>& xD,
      size_t* P, double& det)
  {
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    char u = 'L';
    int n = A.size();
    int lda = A.stepj();
    int lap_p[A.size()];
    int lwork = n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    dsytrf(&u,&n,A.ptr(),&lda,lap_p,work,&lwork,&info);
    if (info < 0) tmv_error("dsytrf returned info < 0");
    xD.Zero();
    int* pi = lap_p;
    for(size_t i=0;i<A.size();++i,++pi) {
      if (*pi-1 != int(i)) {
	if (*pi > 0) {
	  Swap(A.row(i,0,i),A.row(*pi-1,0,i));
	  P[i] = *pi-1;
	  det *= A(i,i);
	} else {
	  TMVAssert(*pi < 0);
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  P[i] = i;
	  P[i+1] = -*pi-1;
	  if (int(i)+1 != -*pi-1) {
	    Swap(A.row(i+1,0,i),A.row(-*pi-1,0,i));
	  }
	  double& xDi = A(i+1,i);
	  xD(i) = xDi;
	  det *= A(i,i)*A(i+1,i+1)-xDi*xDi;
	  xDi = double(0);
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	det *= A(i,i);
	P[i] = i;
      }
    }
  }
  template <> inline void LapHermLU_Decompose(
      const SymMatrixView<complex<double> >& A,
      const VectorView<complex<double> >& xD,
      size_t* P, complex<double>& det)
  {
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    TMVAssert(A.isherm());
    char u = 'L';
    int n = A.size();
    int lda = A.stepj();
    int lap_p[A.size()];
    int lwork = n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    zhetrf(&u,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(work),&lwork,
	&info);
    if (info < 0) tmv_error("zhetrf returned info < 0");
    xD.Zero();
    int* pi = lap_p;
    for(size_t i=0;i<A.size();++i,++pi) {
      if (*pi-1 != int(i)) {
	if (*pi > 0) {
	  Swap(A.row(i,0,i),A.row(*pi-1,0,i));
	  P[i] = *pi-1;
	  det *= REAL(A(i,i));
	} else {
	  TMVAssert(*pi < 0);
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  P[i] = i;
	  P[i+1] = -*pi-1;
	  if (int(i)+1 != -*pi-1) {
	    Swap(A.row(i+1,0,i),A.row(-*pi-1,0,i));
	  }
	  complex<double>& xDi = A(i+1,i).GetRef();
	  xD(i) = xDi;
	  det *= REAL(A(i,i))*REAL(A(i+1,i+1))-NORM(xDi);
	  xDi = double(0);
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	det *= REAL(A(i,i));
	P[i] = i;
      }
    }
  }
  template <> inline void LapSymLU_Decompose(
      const SymMatrixView<complex<double> >& A,
      const VectorView<complex<double> >& xD,
      size_t* P, complex<double>& det)
  {
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    TMVAssert(!A.isherm());
    char u = 'L';
    int n = A.size();
    int lda = A.stepj();
    int lap_p[A.size()];
    int lwork = n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    zsytrf(&u,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(work),&lwork,
	&info);
    if (info < 0) tmv_error("zsytrf returned info < 0");
    xD.Zero();
    int* pi = lap_p;
    for(size_t i=0;i<A.size();++i,++pi) {
      if (*pi-1 != int(i)) {
	if (*pi > 0) {
	  Swap(A.row(i,0,i),A.row(*pi-1,0,i));
	  P[i] = *pi-1;
	  det *= A(i,i);
	} else {
	  TMVAssert(*pi < 0);
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  P[i] = i;
	  P[i+1] = -*pi-1;
	  if (int(i)+1 != -*pi-1) {
	    Swap(A.row(i+1,0,i),A.row(-*pi-1,0,i));
	  }
	  complex<double>& xDi = A(i+1,i).GetRef();
	  xD(i) = xDi;
	  det *= A(i,i)*A(i+1,i+1)-xDi*xDi;
	  xDi = double(0);
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	det *= A(i,i);
	P[i] = i;
      }
    }
  }
#ifndef NOFLOAT
  template <> inline void LapHermLU_Decompose(
      const SymMatrixView<float>& A, const VectorView<float>& xD,
      size_t* P, float& det)
  {
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    char u = 'L';
    int n = A.size();
    int lda = A.stepj();
    int lap_p[A.size()];
    int lwork = n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    ssytrf(&u,&n,A.ptr(),&lda,lap_p,work,&lwork,&info);
    if (info < 0) tmv_error("ssytrf returned info < 0");
    xD.Zero();
    int* pi = lap_p;
    for(size_t i=0;i<A.size();++i,++pi) {
      if (*pi-1 != int(i)) {
	if (*pi > 0) {
	  Swap(A.row(i,0,i),A.row(*pi-1,0,i));
	  P[i] = *pi-1;
	  det *= A(i,i);
	} else {
	  TMVAssert(*pi < 0);
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  P[i] = i;
	  P[i+1] = -*pi-1;
	  if (int(i)+1 != -*pi-1) {
	    Swap(A.row(i+1,0,i),A.row(-*pi-1,0,i));
	  }
	  float& xDi = A(i+1,i);
	  xD(i) = xDi;
	  det *= A(i,i)*A(i+1,i+1)-xDi*xDi;
	  xDi = float(0);
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	det *= A(i,i);
	P[i] = i;
      }
    }
  }
  template <> inline void LapHermLU_Decompose(
      const SymMatrixView<complex<float> >& A,
      const VectorView<complex<float> >& xD,
      size_t* P, complex<float>& det)
  {
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    TMVAssert(A.isherm());
    char u = 'L';
    int n = A.size();
    int lda = A.stepj();
    int lap_p[A.size()];
    int lwork = n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    chetrf(&u,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(work),&lwork,
	&info);
    if (info < 0) tmv_error("chetrf returned info < 0");
    xD.Zero();
    int* pi = lap_p;
    for(size_t i=0;i<A.size();++i,++pi) {
      if (*pi-1 != int(i)) {
	if (*pi > 0) {
	  Swap(A.row(i,0,i),A.row(*pi-1,0,i));
	  P[i] = *pi-1;
	  det *= REAL(A(i,i));
	} else {
	  TMVAssert(*pi < 0);
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  P[i] = i;
	  P[i+1] = -*pi-1;
	  if (int(i)+1 != -*pi-1) {
	    Swap(A.row(i+1,0,i),A.row(-*pi-1,0,i));
	  }
	  complex<float>& xDi = A(i+1,i).GetRef();
	  xD(i) = xDi;
	  det *= REAL(A(i,i))*REAL(A(i+1,i+1))-NORM(xDi);
	  xDi = float(0);
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	det *= REAL(A(i,i));
	P[i] = i;
      }
    }
  }
  template <> inline void LapSymLU_Decompose(
      const SymMatrixView<complex<float> >& A,
      const VectorView<complex<float> >& xD,
      size_t* P, complex<float>& det)
  {
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    TMVAssert(!A.isherm());
    char u = 'L';
    int n = A.size();
    int lda = A.stepj();
    int lap_p[A.size()];
    int lwork = n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    csytrf(&u,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(work),&lwork,
	&info);
    if (info < 0) tmv_error("csytrf returned info < 0");
    xD.Zero();
    int* pi = lap_p;
    for(size_t i=0;i<A.size();++i,++pi) {
      if (*pi-1 != int(i)) {
	if (*pi > 0) {
	  Swap(A.row(i,0,i),A.row(*pi-1,0,i));
	  P[i] = *pi-1;
	  det *= A(i,i);
	} else {
	  TMVAssert(*pi < 0);
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  P[i] = i;
	  P[i+1] = -*pi-1;
	  if (int(i)+1 != -*pi-1) {
	    Swap(A.row(i+1,0,i),A.row(-*pi-1,0,i));
	  }
	  complex<float>& xDi = A(i+1,i).GetRef();
	  xD(i) = xDi;
	  det *= A(i,i)*A(i+1,i+1)-xDi*xDi;
	  xDi = float(0);
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	det *= A(i,i);
	P[i] = i;
      }
    }
  }
#endif // NOFLOAT
#endif // LAP
  template <class T> inline void SymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, T& det)
  {
    TMVAssert(A.size() > 0);
    TMVAssert(xD.size()+1 == A.size());
    if (A.isconj()) SymLU_Decompose(A.Conjugate(),xD,P,det);
    else if (A.uplo() == Upper) {
      if (A.isherm()) SymLU_Decompose(A.Adjoint(),xD,P,det);
      else SymLU_Decompose(A.Transpose(),xD,P,det);
    }
    else {
#ifdef LAP
      if (A.iscm())
	if (IsReal(T()) || A.isherm())
	  LapHermLU_Decompose(A,xD,P,det);
	else
	  LapSymLU_Decompose(A,xD,P,det);
      else
#endif
	if (IsReal(T()) || A.isherm())
	  NonLapHermLU_Decompose(A,xD,P,det);
	else
	  NonLapSymLU_Decompose(A,xD,P,det);
    }
  }

#define InstFile "TMV_SymLUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


