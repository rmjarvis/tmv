#ifndef TMV_BandMatrixArith_H
#define TMV_BandMatrixArith_H

#include "TMV_BandMatrix.h"

#define CT complex<T>
#define CCT ConjRef<complex<T> >
#define VCT VarConjRef<complex<T> >

namespace tmv {

  // y = alpha * A * x + beta * y
  template <class T, class Ta, class Tx> void MultMV(const T alpha,
      const GenBandMatrix<Ta>& A, const GenVector<Tx>& x, 
      const T beta, const VectorView<T>& y);

  // A = alpha * A
  template <class T> void MultXM(const T alpha, const BandMatrixView<T>& A);

  // B = alpha * A + beta * B
  template <class T, class Ta> void AddMM(const T alpha,
      const GenBandMatrix<Ta>& A, const T beta, const BandMatrixView<T>& B);
  template <class T, class Ta> void AddMM(const T alpha,
      const GenBandMatrix<Ta>& A, const T beta, const MatrixView<T>& B);

  // C = alpha * A * B + beta * C
  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B, 
      const T beta, const MatrixView<T>& C);
  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B, 
      const T beta, const BandMatrixView<T>& C);
  template <class T, class Ta, class Tb> inline void MultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenBandMatrix<Tb>& B, 
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    MultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose()); 
  }

  template <class T, class Ta> void ElementProd(
      const T alpha, const GenBandMatrix<Ta>& A, const BandMatrixView<T>& B);
  template <class T, class Ta, class Tb> void AddElementProd(
      const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const BandMatrixView<T>& C);

  template <class T> class BandMatrixComposite : 
    public GenBandMatrix<T>
  {
    public:

      BandMatrixComposite(StorageType s) : 
	GenBandMatrix<T>(s,NonConj,s==DiagMajor?0:1), 
        itsm(0), itsm1(0), si(0), sj(0), ds(0) {}
      BandMatrixComposite(const BandMatrixComposite<T>& rhs) :
	GenBandMatrix<T>(rhs),
        itsm(rhs.itsm), itsm1(0),
	si(rhs.si), sj(rhs.sj), ds(rhs.ds) {}
      virtual ~BandMatrixComposite() {}
      virtual void AssignTo(const BandMatrixView<T>&) const = 0;
      inline const T* cptr() const
      { 
	if (!itsm1.get()) {
	  size_t len = BandStorageLength(stor(),colsize(),rowsize(),
	      nlo(),nhi());
	  itsm1.reset(new T[len]);
	  itsm = isdm() ? itsm1.get()-nlo()*stepi() : itsm1.get();
	  BandMatrixView<T>(itsm,colsize(),rowsize(),nlo(),nhi(),
	      stepi(),stepj(),diagstep(),stor(),NonConj,isdm()?0:len 
	      FIRSTLAST1(itsm1.get(),itsm1.get()+len) ) = *this;
	}
	return itsm;
      }
      inline int stepi() const
      {
	if (si == 0) si = iscm() ? 1 : isrm() ? nlo()+nhi() :
	  rowsize()>=colsize() ? -int(colsize())+1 : -int(rowsize());
	return si;
      }
      inline int stepj() const
      {
	if (sj == 0) sj = isrm() ? 1 : iscm() ? nlo()+nhi() :
	  rowsize()>=colsize() ? int(colsize()) : int(rowsize())+1;
	return sj;
      }
      inline int diagstep() const
      {
	if (ds == 0) ds = isdm() ? 1 : nlo()+nhi()+1;
	return ds;
      }
      using GenBandMatrix<T>::colsize;
      using GenBandMatrix<T>::rowsize;
      using GenBandMatrix<T>::nlo;
      using GenBandMatrix<T>::nhi;
      using GenBandMatrix<T>::stor;
      using GenBandMatrix<T>::isrm;
      using GenBandMatrix<T>::iscm;
      using GenBandMatrix<T>::isdm;

    private:
      mutable T* itsm;
      mutable auto_array<T> itsm1;
      mutable int si;
      mutable int sj;
      mutable int ds;
  };

  template <class T, StorageType S, class Tx> 
    inline BandMatrix<T,S>& operator+=(BandMatrix<T,S>& m, const Tx& x) 
    { m.View() += x; return m; }

  template <class T, StorageType S, class Tx> 
    inline BandMatrix<T,S>& operator-=(BandMatrix<T,S>& m, const Tx& x) 
    { m.View() -= x; return m; }

  template <class T, StorageType S, class Tx> 
    inline BandMatrix<T,S>& operator*=(BandMatrix<T,S>& m, const Tx& x) 
    { m.View() *= x; return m; }

  template <class T, StorageType S, class Tx> 
    inline BandMatrix<T,S>& operator/=(BandMatrix<T,S>& m, const Tx& x) 
    { m.View() /= x; return m; }

  template <class T, StorageType S, class Tx>
    inline BandMatrix<T,S>& operator%=(BandMatrix<T,S>& m, const Tx& x) 
    { m.View() %= x; return m; }

  template <class T, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m1, const GenBandMatrix<T2>& m2)
  { BandMatrixViewOf(m1,m2.nlo(),m2.nhi()) += m2; return m1; }

  template <class T, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m1, const GenBandMatrix<T2>& m2)
  { BandMatrixViewOf(m1,m2.nlo(),m2.nhi()) -= m2; return m1; }

  //
  // Scalar * BandMatrix
  //

  template <class T, class Tm> class ProdXB : 
    public BandMatrixComposite<T> 
  {
    public:
      ProdXB(const T _x, const GenBandMatrix<Tm>& _m) :
	BandMatrixComposite<T>(_m.stor()!=NoMajor ? _m.stor() : DiagMajor),
	x(_x), m(&_m) {}
      inline size_t colsize() const { return m->colsize(); }
      inline size_t rowsize() const { return m->rowsize(); }
      inline int nlo() const { return m->nlo(); }
      inline int nhi() const { return m->nhi(); }
      inline T GetX() const { return x; }
      inline const GenBandMatrix<Tm>& GetM() const { return *m; }
      inline void AssignTo(const BandMatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() == nlo());
	TMVAssert(m0.nhi() == nhi());
	m0 = *m;
	MultXM(x,m0);
      }
    private:
      const T x;
      const GenBandMatrix<Tm>*const m;
  };

  // m*=x
  template <class T> inline const BandMatrixView<T>& operator*=(
      const BandMatrixView<T>& m, T x) 
  { MultXM(x,m); return m; }

  template <class T> inline const BandMatrixView<CT>& operator*=(
      const BandMatrixView<CT>& m, T x) 
  { MultXM(CT(x),m); return m; }

  template <class T> inline const BandMatrixView<CT>& operator*=(
      const BandMatrixView<CT>& m, CCT x) 
  { MultXM(CT(x),m); return m; }

  template <class T> inline const BandMatrixView<CT>& operator*=(
      const BandMatrixView<CT>& m, VCT x) 
  { MultXM(CT(x),m); return m; }

  // m/=x
  template <class T> inline const BandMatrixView<T>& operator/=(
      const BandMatrixView<T>& m, T x) 
  { MultXM(RealType(T)(1)/x,m); return m; }

  template <class T> inline const BandMatrixView<CT>& operator/=(
      const BandMatrixView<CT>& m, T x) 
  { MultXM(CT(T(1)/x),m); return m; }

  template <class T> inline const BandMatrixView<CT>& operator/=(
      const BandMatrixView<CT>& m, CCT x) 
  { MultXM(T(1)/CT(x),m); return m; }

  template <class T> inline const BandMatrixView<CT>& operator/=(
      const BandMatrixView<CT>& m, VCT x) 
  { MultXM(T(1)/CT(x),m); return m; }

#define GENMATRIX GenBandMatrix
#define PRODXM ProdXB
#include "TMV_AuxProdXM.h"


  //
  // BandMatrix + Scalar
  //

  template <class T, class Tm> class SumBX : 
    public BandMatrixComposite<T>
  {
    public:
      SumBX(T _x1, const GenBandMatrix<Tm>& _m, T _x2) :
	BandMatrixComposite<T>(_m.stor()!=NoMajor ? _m.stor() : DiagMajor),
	x1(_x1), m(&_m), x2(_x2)
	{ TMVAssert(m->IsSquare()); }
      inline size_t colsize() const { return m->colsize(); }
      inline size_t rowsize() const { return m->rowsize(); }
      inline int nlo() const { return m->nlo(); }
      inline int nhi() const { return m->nhi(); }
      inline T GetX1() const { return *x1; }
      inline const GenBandMatrix<Tm>& GetM() const { return *m; }
      inline T GetX2() const { return *x2; }
      inline void AssignTo(const BandMatrixView<T>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() == nlo());
	TMVAssert(m0.nhi() == nhi());
	if (x1 == T(1)) m0 = *m;
	else m0 = x1*(*m);
	m0.diag().AddToAll(x2);
      }
    private:
      const T x1;
      const GenBandMatrix<Tm>* m;
      const T x2;
  };

  // m+=x
  template <class T> inline const BandMatrixView<T>& operator+=(
      const BandMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(x); return m; 
  }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m, T x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(x); return m; 
  }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(x)); return m; 
  }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(x)); return m; 
  }

  // m-=x
  template <class T> inline const BandMatrixView<T>& operator-=(
      const BandMatrixView<T>& m, T x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(-x); return m; 
  }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m, T x) 
  { 
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(-x); return m; 
  }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(-CT(x)); return m; 
  }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(-CT(x)); return m; 
  }

#define SUMMX SumBX
#include "TMV_AuxSumMX.h"
#undef SUMMX


  //
  // BandMatrix + BandMatrix
  //
  
  template <class T, class T1, class T2> class SumBB : 
    public BandMatrixComposite<T> 
  {
    public:
      SumBB(T _x1, const GenBandMatrix<T1>& _m1, 
	  T _x2, const GenBandMatrix<T2>& _m2) :
	BandMatrixComposite<T>(_m1.stor()!=NoMajor ? _m1.stor() :
	    _m2.stor()!=NoMajor ? _m2.stor() : DiagMajor),
	x1(_x1),m1(&_m1),x2(_x2),m2(&_m2),
	itsnlo(max(m1->nlo(),m2->nlo())), itsnhi(max(m1->nhi(),m2->nhi()))
	{ 
	  TMVAssert(m1->rowsize() == m2->rowsize()); 
	  TMVAssert(m1->colsize() == m2->colsize()); 
	}
      inline size_t colsize() const { return m1->colsize(); }
      inline size_t rowsize() const { return m1->rowsize(); }
      inline int nlo() const { return itsnlo; }
      inline int nhi() const { return itsnhi; }
      inline T GetX1() const { return x1; }
      inline const GenBandMatrix<T1>& GetM1() const { return *m1; }
      inline T GetX2() const { return x2; }
      inline const GenBandMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const BandMatrixView<T>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	if (m0.nlo() > nlo() || m0.nhi() > nhi()) {
	  AssignTo(m0.Diags(-nlo(),nhi()+1));
	  m0.Diags(-m0.nlo(),-nlo()).Zero();
	  m0.Diags(nhi()+1,m0.nhi()+1).Zero();
	} else if (m0.SameStorageAs(*m1)) {
	  m0 = *m1;
	  AddMM(x2,*m2,x1,m0);
	} else if (m0.SameStorageAs(*m2) || x1 != T(1)) {
	  m0 = *m2;
	  AddMM(x1,*m1,x2,m0);
	} else {
	  m0 = *m1;
	  AddMM(x2,*m2,x1,m0);
        }
      }
    private:
      T x1;
      const GenBandMatrix<T1>* m1;
      T x2;
      const GenBandMatrix<T2>* m2;
      int itsnlo,itsnhi;
  };

  template <class T> inline const BandMatrixView<T>& operator+=(
      const BandMatrixView<T>& m1, const GenBandMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVAssert(m1.nlo() >= m2.nlo());
    TMVAssert(m1.nhi() >= m2.nhi());
    AddMM(T(1),m2,T(1),m1); return m1; 
  }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m1, const GenBandMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVAssert(m1.nlo() >= m2.nlo());
    TMVAssert(m1.nhi() >= m2.nhi());
    AddMM(CT(1),m2,CT(1),m1); return m1; 
  }

  template <class T> inline const BandMatrixView<T>& operator-=(
      const BandMatrixView<T>& m1, const GenBandMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVAssert(m1.nlo() >= m2.nlo());
    TMVAssert(m1.nhi() >= m2.nhi());
    AddMM(T(-1),m2,T(1),m1); return m1; 
  }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m1, const GenBandMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVAssert(m1.nlo() >= m2.nlo());
    TMVAssert(m1.nhi() >= m2.nhi());
    AddMM(CT(-1),m2,CT(1),m1); return m1; 
  }

  template <class T, class T2> inline const BandMatrixView<T>& operator+=(
      const BandMatrixView<T>& m, const ProdXB<T,T2>& pxm)
  {
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    TMVAssert(m.nlo() >= pxm.nlo());
    TMVAssert(m.nhi() >= pxm.nhi());
    AddMM(pxm.GetX(),pxm.GetM(),T(1),m);
    return m;
  }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m, const ProdXB<T,T>& pxm)
  {
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    TMVAssert(m.nlo() >= pxm.nlo());
    TMVAssert(m.nhi() >= pxm.nhi());
    AddMM(CT(pxm.GetX()),pxm.GetM(),CT(1),m);
    return m;
  }

  template <class T, class T2> inline const BandMatrixView<T>& operator-=(
      const BandMatrixView<T>& m, const ProdXB<T,T2>& pxm)
  {
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    TMVAssert(m.nlo() >= pxm.nlo());
    TMVAssert(m.nhi() >= pxm.nhi());
    AddMM(-pxm.GetX(),pxm.GetM(),T(1),m);
    return m;
  }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m, const ProdXB<T,T>& pxm)
  {
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    TMVAssert(m.nlo() >= pxm.nlo());
    TMVAssert(m.nhi() >= pxm.nhi());
    AddMM(CT(-pxm.GetX()),pxm.GetM(),CT(1),m);
    return m;
  }

#define SUMMM SumBB
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXB
#include "TMV_AuxSumMM.h"
#undef SUMMM


  //
  // BandMatrix * BandMatrix
  //

  template <class T, class T1, class T2> class ProdBB : 
    public BandMatrixComposite<T>
  {
    public:
      ProdBB(T _x, const GenBandMatrix<T1>& _m1,
	  const GenBandMatrix<T2>& _m2) :
	BandMatrixComposite<T>(_m1.isrm() ? RowMajor :
	    _m2.iscm() ? ColMajor :
	    (_m1.iscm() && _m2.isrm()) ? ColMajor : DiagMajor),
	x(_x), m1(&_m1), m2(&_m2),
	itsnlo(min(int(colsize()-1),m1->nlo()+m2->nlo())), 
	itsnhi(min(int(rowsize()-1),m1->nhi()+m2->nhi()))
      { TMVAssert( m1->rowsize() == m2->colsize()); }
      inline size_t colsize() const { return m1->colsize(); }
      inline size_t rowsize() const { return m2->rowsize(); }
      inline int nlo() const { return itsnlo; }
      inline int nhi() const { return itsnhi; }
      inline T GetX() const { return x; }
      inline const GenBandMatrix<T1>& GetM1() const { return *m1; }
      inline const GenBandMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const BandMatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	MultMM(x,*m1,*m2,T(0),m0);
      }
    private:
      T x;
      const GenBandMatrix<T1>*const m1;
      const GenBandMatrix<T2>*const m2;
      int itsnlo, itsnhi;
  };

  template <class T> inline const BandMatrixView<T>& operator*=(
      const BandMatrixView<T>& m1, const GenBandMatrix<T>& m2)
  {
    TMVAssert(m2.nlo() == 0 || m1.nlo() == int(m1.colsize())-1);
    TMVAssert(m2.nhi() == 0 || m1.nhi() == int(m1.rowsize())-1);
    MultMM(T(1),m1,m2,T(0),m1); return m1; 
  }

  template <class T> inline const BandMatrixView<CT>& operator*=(
      const BandMatrixView<CT>& m1, const GenBandMatrix<T>& m2)
  {
    TMVAssert(m2.nlo() == 0 || m1.nlo() == int(m1.colsize())-1);
    TMVAssert(m2.nhi() == 0 || m1.nhi() == int(m1.rowsize())-1);
    MultMM(CT(1),m1,m2,CT(0),m1); return m1; 
  }

  template <class T, class T1, class T2>
    inline const BandMatrixView<T>& operator+=(
	const BandMatrixView<T>& m, const ProdBB<T,T1,T2>& pmm)
    {
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      TMVAssert(m.nlo() >= pmm.nlo());
      TMVAssert(m.nhi() >= pmm.nhi());
      MultMM(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); 
      return m; 
    }

  template <class T, class T1, class T2>
    inline const BandMatrixView<T>& operator-=(
	const BandMatrixView<T>& m, const ProdBB<T,T1,T2>& pmm)
    { 
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      TMVAssert(m.nlo() >= pmm.nlo());
      TMVAssert(m.nhi() >= pmm.nhi());
      MultMM(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); 
      return m; 
    }

#define PRODMM ProdBB
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM



#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenMatrix
#define PRODXM1 ProdXM
#define SUMMM SumMB
#define PRODMM ProdMB
#define QUOTXM QuotXB
#define QUOTMM QuotMB
#define RQUOTMM RQuotMB
#define TQUOTMM TransientQuotMB
#define TRQUOTMM TransientRQuotMB

#include "TMV_AuxMatComposite1.h"
#include "TMV_AuxMatComposite3.h"

#define PRODMV ProdBV
#define QUOTVM QuotVB
#define RQUOTVM RQuotVB
#include "TMV_AuxVecComposite.h"
#undef PRODMV
#undef QUOTVM
#undef RQUOTVM
#undef PRODXM
#undef GENMATRIX

  // B/B -> Matrix.  Use TransientQuot
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenBandMatrix
#define PRODXM1 ProdXB
#include "TMV_AuxTQuotMM.h"
  

#undef GENMATRIX2
#undef PRODXM2
#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM
#define GENMATRIX2 GenMatrix
#define PRODXM2 ProdXM
#define SUMMMa SumMB
#define SUMMM SumBM
#define PRODMM ProdBM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM
#define QUOTXM QuotXM

#include "TMV_AuxMatComposite2.h"
#include "TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#undef SUMMM
#undef SUMMMa
#undef PRODMM
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif
