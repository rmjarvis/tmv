
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t SYM_LU_BLOCKSIZE = TMV_BLOCKSIZE;
  const size_t SYM_LU_BLOCKSIZE1 = TMV_BLOCKSIZE;
  const size_t SYM_LU_BLOCKSIZE2 = TMV_BLOCKSIZE/2;
#else
  const size_t SYM_LU_BLOCKSIZE = 48;
  const size_t SYM_LU_BLOCKSIZE1 = 64;
  const size_t SYM_LU_BLOCKSIZE2 = 32;
#endif

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

  template <class T, class T1ab, class T1c> void LMultEq_2x2(
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

  template <class T, class T1> inline void SymLMultEq_2x2(
      T1 a, T1 b, T1 c, const MatrixView<T>& m)
  {
    TMVAssert(m.colsize() == 2);
    if (m.isconj())
      LMultEq_2x2(CONJ(a),CONJ(b),CONJ(c),CONJ(c),m.Conjugate());
    else
      LMultEq_2x2(a,b,c,c,m);
  }

  template <class T, class T1> inline void SymRMultEq_2x2(
      T1 a, T1 b, T1 c, const MatrixView<T>& m)
  {
    TMVAssert(m.rowsize() == 2);
    if (m.isconj())
      LMultEq_2x2(CONJ(a),CONJ(b),CONJ(c),CONJ(c),m.Adjoint());
    else
      LMultEq_2x2(a,b,c,c,m.Transpose());
  }

  template <class T, class T1> inline void HermLMultEq_2x2(
      RealType(T1) a, RealType(T1) b, T1 c, const MatrixView<T>& m)
  {
    TMVAssert(m.colsize() == 2);
    if (m.isconj())
      LMultEq_2x2(a,b,CONJ(c),c,m.Conjugate());
    else
      LMultEq_2x2(a,b,c,CONJ(c),m);
  }

  template <class T, class T1> inline void HermRMultEq_2x2(
      RealType(T1) a, RealType(T1) b, T1 c, const MatrixView<T>& m)
  {
    TMVAssert(m.rowsize() == 2);
    if (m.isconj())
      LMultEq_2x2(a,b,c,CONJ(c),m.Adjoint());
    else
      LMultEq_2x2(a,b,CONJ(c),c,m.Transpose());
  }

  template <class T> void SymInvert_2x2(T& a, T& b, T& c)
  {
    // Invert matrix [ a  c ] -->  1/(ab-c^2) [ b -c ]
    //               [ c  b ]                 [ -c a ]
    T d = a*b-c*c;
    swap(a,b);
    a/=d;
    b/=d;
    c = -c/d;
  }

  template <class T> void HermInvert_2x2(RealType(T)& a, RealType(T)& b, T& c)
  {
    // Invert matrix [ a  c* ] -->  1/(ab-|c|^2) [ b -c* ]
    //               [ c  b  ]                   [ -c a  ]
    RealType(T) d = a*b-NORM(c);
    swap(a,b);
    a/=d;
    b/=d;
    c = -c/d;
  }

  template <class T, class T1> void SymPseudoDiag_LDivEq(
      const GenVector<T1>& D, const GenVector<T1>& xD,
      const MatrixView<T>& m)
  {
    TMVAssert(D.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());

    const T1* Di = D.cptr();
    const T1* xDi = xD.cptr();

    const size_t N = D.size();
    const int sd = D.step();
    const int sx = xD.step();

    for(size_t i=0;i<N;) {
      if (i==N-1 || *xDi == T1(0)) {
	m.row(i) /= (D.isconj() ? CONJ(*Di) : *Di);
	Di+=sd; xDi+=sx; ++i;
      } else {
	T1 x = (D.isconj() ? CONJ(*Di) : *Di);
	T1 y = (D.isconj() ? CONJ(*(Di+=sd)) : *(Di+=sd));
	T1 z = *xDi;
	SymInvert_2x2(x,y,z);
	SymLMultEq_2x2(x,y,z,m.Rows(i,i+2));
	Di+=sd,xDi+=2*sx,i+=2;
      }
    }
  }

  template <class T, class T1> void HermPseudoDiag_LDivEq(
      const GenVector<RealType(T1)>& D, const GenVector<T1>& xD,
      const MatrixView<T>& m)
  {
    TMVAssert(D.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());

    const RealType(T1)* Di = D.cptr();
    const T1* xDi = xD.cptr();

    const size_t N = D.size();
    const int sd = D.step();
    const int sx = xD.step();

    for(size_t i=0;i<N;) {
      if (i==N-1 || *xDi == T1(0)) {
	m.row(i) /= *Di;
	Di+=sd; xDi+=sx; ++i;
      } else {
	RealType(T1) x = *Di;
	RealType(T1) y = *(Di+=sd);
	T1 z = (xD.isconj() ? CONJ(*xDi) : *xDi);
	HermInvert_2x2(x,y,z);
	HermLMultEq_2x2(x,y,z,m.Rows(i,i+2));
	Di+=sd,xDi+=2*sx,i+=2;
      }
    }
  }

  template <class T, class T1> void SymPseudoDiag_LMultEq(
      const GenVector<T1>& D, const GenVector<T1>& xD,
      const MatrixView<T>& m)
  {
    TMVAssert(D.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());

    const T1* Di = D.cptr();
    const T1* xDi = xD.cptr();

    const size_t N = D.size();
    const int sd = D.step();
    const int sx = xD.step();

    for(size_t i=0;i<N;) {
      if (i==N-1 || *xDi == T1(0)) {
	m.row(i) *= (D.isconj() ? CONJ(*Di) : *Di);
	Di+=sd; xDi+=sx; ++i;
      } else {
	T1 x = (D.isconj() ? CONJ(*Di) : *Di);
	T1 y = (D.isconj() ? CONJ(*(Di+=sd)) : *(Di+=sd));
	T1 z = (xD.isconj() ? CONJ(*xDi) : *xDi);
	SymLMultEq_2x2(x,y,z,m.Rows(i,i+2));
	Di+=sd,xDi+=2*sx,i+=2;
      }
    }
  }

  template <class T, class T1> void HermPseudoDiag_LMultEq(
      const GenVector<RealType(T1)>& D, const GenVector<T1>& xD,
      const MatrixView<T>& m)
  {
    TMVAssert(D.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());

    const RealType(T1)* Di = D.cptr();
    const T1* xDi = xD.cptr();

    const size_t N = D.size();
    const int sd = D.step();
    const int sx = xD.step();

    for(size_t i=0;i<N;) {
      if (i==N-1 || *xDi == T1(0)) {
	m.row(i) *= *Di;
	Di+=sd; xDi+=sx; ++i;
      } else {
	RealType(T1) x = *Di;
	RealType(T1) y = *(Di+=sd);
	T1 z = (xD.isconj() ? CONJ(*xDi) : *xDi);
	HermLMultEq_2x2(x,y,z,m.Rows(i,i+2));
	Di+=sd,xDi+=2*sx,i+=2;
      }
    }
  }

  template <class T, class T1> void SymLU_LDivEq(
      const GenSymMatrix<T1>& LL, const GenVector<T1>& xD, const size_t* P, 
      const MatrixView<T>& m)
  {
    // Solve P L D LT Pt x = m:
    TMVAssert(IsComplex(T1()) && !LL.isherm());
    TMVAssert(LL.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());

#ifdef XDEBUG
    Matrix<T1> DD(LL.size(),LL.size(),T1(0));
    DD.diag() = LL.diag();
    DD.diag(-1) = xD;
    DD.diag(1) = xD;
    LowerTriMatrix<T1> L = LL.LowerTri().MakeUnitDiag();
    Matrix<T1> LDL = L*DD*L.Transpose();
    LDL.ReversePermuteRows(P);
    LDL.ReversePermuteCols(P);
    Matrix<T> m0 = m;
#endif

    m.PermuteRows(P);
    m /= LL.LowerTri().MakeUnitDiag();
    SymPseudoDiag_LDivEq(LL.diag(),xD,m);
    m /= LL.UpperTri().MakeUnitDiag();
    m.ReversePermuteRows(P);

#ifdef XDEBUG
    Matrix<T> mm = LDL*m;
    if (Norm(mm-m0) > 0.001*max(RealType(T)(1),Norm(m0))) {
      cerr<<"SymLU_LDivEq: m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"LDL = "<<LDL<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"LDL*m = "<<mm<<endl;
      abort();
    }
#endif
  }

  template <class T, class T1> void HermLU_LDivEq(
      const GenSymMatrix<T1>& LL, const GenVector<T1>& xD, const size_t* P, 
      const MatrixView<T>& m)
  {
    // Solve P L D Lt Pt x = m:
    TMVAssert(IsReal(T1()) || LL.isherm());
    TMVAssert(LL.size() == m.colsize());
    TMVAssert(xD.size()+1 == m.colsize());

#ifdef XDEBUG
    Matrix<T> m0 = m;
    Matrix<T1> DD(LL.size(),LL.size(),T1(0));
    DD.diag() = LL.diag().Real();
    DD.diag(-1) = xD;
    DD.diag(1) = xD.Conjugate();
    LowerTriMatrix<T1> L = LL.LowerTri().MakeUnitDiag();
    Matrix<T1> LDL = L*DD*L.Adjoint();
    LDL.ReversePermuteRows(P);
    LDL.ReversePermuteCols(P);
#endif

    m.PermuteRows(P);
    m /= LL.LowerTri().MakeUnitDiag();
    HermPseudoDiag_LDivEq(LL.diag().Real(),xD,m);
    m /= LL.UpperTri().MakeUnitDiag();
    m.ReversePermuteRows(P);

#ifdef XDEBUG
    Matrix<T> mm = LDL*m;
    if (Norm(mm-m0) > 0.001*max(RealType(T)(1),Norm(m0))) {
      cerr<<"HermPseudoDiag_LDivEq: m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"L = "<<L<<endl;
      cerr<<"DD = "<<DD<<endl;
      cerr<<"D = "<<LL.diag().Real()<<endl;
      cerr<<"xD = "<<xD<<endl;
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

  template <class T, class T1> inline void SymLU_RDivEq(
      const GenSymMatrix<T1>& LL, const GenVector<T1>& xD, const size_t* P, 
      const MatrixView<T>& m)
  {
    // Solve x P L D LT Pt = m:
    TMVAssert(LL.size() == m.rowsize());
    TMVAssert(xD.size()+1 == m.rowsize());
    SymLU_LDivEq(LL,xD,P,m.Transpose());
  }

  template <class T, class T1> inline void HermLU_RDivEq(
      const GenSymMatrix<T1>& LL, const GenVector<T1>& xD, const size_t* P, 
      const MatrixView<T>& m)
  {
    // Solve x P L D Lt Pt = m:
    TMVAssert(LL.size() == m.rowsize());
    TMVAssert(xD.size()+1 == m.rowsize());
    HermLU_LDivEq(LL,xD,P,m.Adjoint());
  }

#define InstFile "TMV_SymLUDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


