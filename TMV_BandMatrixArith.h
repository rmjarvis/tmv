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
  template <class T, class Ta> inline void AddMM(const T alpha,
      const GenBandMatrix<Ta>& A, const T beta, const MatrixView<T>& B)
  { 
    AddMM(alpha,A,beta,BandMatrixViewOf(B,A.nlo(),A.nhi())); 
    if (beta != T(1)) {
      for(int i=-int(B.colsize())+1;i<-A.nlo();i++) B.diag(i) *= beta;
      for(int i=A.nhi()+1;i<int(B.rowsize());i++) B.diag(i) *= beta;
    }
  }

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
  { MultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose()); }

  template <class T> class BandMatrixComposite : 
    public GenBandMatrix<T>
  {
    public:

      BandMatrixComposite(StorageType s) : 
	GenBandMatrix<T>(s,NonConj), inst(0) {}
      virtual ~BandMatrixComposite() { if (inst) delete inst; };
      virtual void AssignTo(const BandMatrixView<T>&) const = 0;
      inline const T* cptr() const
      { 
	if (!inst) {
	  if (stor() == RowMajor) inst = new BandMatrix<T,RowMajor>(*this);
	  else if (stor() == ColMajor) inst = new BandMatrix<T,ColMajor>(*this);
	  else inst = new BandMatrix<T,DiagMajor>(*this);
	}
	return inst->cptr(); 
      }
      inline int stepi() const 
      { 
	if (inst) return inst->stepi();
	else if (this->isrm()) return this->nlo()+this->nhi();
	else if (this->iscm()) return 1;
	else if (this->rowsize() >= this->colsize()) 
	  return -int(this->colsize())+1;
	else return -int(this->rowsize());
      }
      inline int stepj() const 
      { 
	if (inst) return inst->stepj();
	else if (this->isrm()) return 1;
	else if (this->iscm()) return this->nlo()+this->nhi();
	else if (this->rowsize() >= this->colsize()) 
	  return int(this->colsize());
	else return int(this->rowsize())+1;
      }
      inline int diagstep() const
      { 
	if (inst) return inst->diagstep();
	else if (this->isdm()) return 1;
	else return this->nlo()+this->nhi()+1;
      }
      using GenBandMatrix<T>::stor;

    private:
      mutable GenBandMatrix<T>* inst;
  };

  template <class T, StorageType S, class T2> 
    inline BandMatrix<T,S>& operator+=(BandMatrix<T,S>& m1, const T2& x2) 
    { m1.View() += x2; return m1; }

  template <class T, StorageType S, class T2> 
    inline BandMatrix<T,S>& operator-=(BandMatrix<T,S>& m1, const T2& x2) 
    { m1.View() -= x2; return m1; }

  template <class T, StorageType S, class T2> 
    inline BandMatrix<T,S>& operator*=(BandMatrix<T,S>& m1, const T2& x2) 
    { m1.View() *= x2; return m1; }

  template <class T, StorageType S, class T2> 
    inline BandMatrix<T,S>& operator/=(BandMatrix<T,S>& m1, const T2& x2) 
    { m1.View() /= x2; return m1; }

  template <class T, StorageType S, class T2>
    inline BandMatrix<T,S>& operator%=(BandMatrix<T,S>& m1, const T2& x2) 
    { m1.View() %= x2; return m1; }

  template <class T, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m1, const GenBandMatrix<T2>& m2)
  { BandMatrixViewOf(m1,m2.nlo(),m2.nhi()) += m2; return m1; }

  template <class T, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m1, const GenBandMatrix<T2>& m2)
  { BandMatrixViewOf(m1,m2.nlo(),m2.nhi()) -= m2; return m1; }

  //
  // BandMatrix + Scalar
  //

  template <class T, class T2> class SumBX : 
    public BandMatrixComposite<T>
  {
    public:
      SumBX(T _x1, const GenBandMatrix<T2>& _m2, T _x3) :
	BandMatrixComposite<T>(_m2.stor()!=NoMajor ? _m2.stor() : DiagMajor),
	x1(_x1), m2(&_m2), x3(_x3)
	{ TMVAssert(m2->IsSquare()); }
      inline size_t colsize() const { return m2->colsize(); }
      inline size_t rowsize() const { return m2->rowsize(); }
      inline int nlo() const { return m2->nlo(); }
      inline int nhi() const { return m2->nhi(); }
      inline T GetX1() const { return *x3; }
      inline const GenBandMatrix<T2>& GetM2() const { return *m2; }
      inline T GetX3() const { return *x3; }
      inline void AssignTo(const BandMatrixView<T>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() == nlo());
	TMVAssert(m0.nhi() == nhi());
	if (x1 == T(1)) m0 = *m2;
	else m0 = x1*(*m2);
	m0.diag().AddToAll(x3);
      }
    private:
      const T x1;
      const GenBandMatrix<T2>* m2;
      const T x3;
  };

  // m+=x
  template <class T> inline const BandMatrixView<T>& operator+=(
      const BandMatrixView<T>& m1, T x2) 
  { m1.diag().AddToAll(x2); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m1, T x2) 
  { m1.diag().AddToAll(x2); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m1, CCT x2) 
  { m1.diag().AddToAll(CT(x2)); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m1, VCT x2) 
  { m1.diag().AddToAll(CT(x2)); return m1; }

  // m-=x
  template <class T> inline const BandMatrixView<T>& operator-=(
      const BandMatrixView<T>& m1, T x2) 
  { m1.diag().AddToAll(-x2); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m1, T x2) 
  { m1.diag().AddToAll(-x2); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m1, CCT x2) 
  { m1.diag().AddToAll(-CT(x2)); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m1, VCT x2) 
  { m1.diag().AddToAll(-CT(x2)); return m1; }


  //
  // Scalar * BandMatrix
  //

  template <class T, class T2> class ProdXB : 
    public BandMatrixComposite<T> 
  {
    public:
      ProdXB(const T _x1, const GenBandMatrix<T2>& _m2) :
	BandMatrixComposite<T>(_m2.stor()!=NoMajor ? _m2.stor() : DiagMajor),
	x1(_x1), m2(&_m2) {}
      inline size_t colsize() const { return m2->colsize(); }
      inline size_t rowsize() const { return m2->rowsize(); }
      inline int nlo() const { return m2->nlo(); }
      inline int nhi() const { return m2->nhi(); }
      inline T GetX1() const { return x1; }
      inline const GenBandMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const BandMatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() == nlo());
	TMVAssert(m0.nhi() == nhi());
	m0 = *m2;
	MultXM(x1,m0);
      }
    private:
      const T x1;
      const GenBandMatrix<T2>*const m2;
  };

  // m*=x
  template <class T> inline const BandMatrixView<T>& operator*=(
      const BandMatrixView<T>& m1, T x2) 
  { MultXM(x2,m1); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator*=(
      const BandMatrixView<CT>& m1, T x2) 
  { MultXM(CT(x2),m1); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator*=(
      const BandMatrixView<CT>& m1, CCT x2) 
  { MultXM(CT(x2),m1); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator*=(
      const BandMatrixView<CT>& m1, VCT x2) 
  { MultXM(CT(x2),m1); return m1; }

  // m/=x
  template <class T> inline const BandMatrixView<T>& operator/=(
      const BandMatrixView<T>& m1, T x2) 
  { MultXM(RealType(T)(1)/x2,m1); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator/=(
      const BandMatrixView<CT>& m1, T x2) 
  { MultXM(CT(T(1)/x2),m1); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator/=(
      const BandMatrixView<CT>& m1, CCT x2) 
  { MultXM(T(1)/CT(x2),m1); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator/=(
      const BandMatrixView<CT>& m1, VCT x2) 
  { MultXM(T(1)/CT(x2),m1); return m1; }

  //
  // BandMatrix + BandMatrix
  //
  
  template <class T, class T2, class T4> class SumBB : 
    public BandMatrixComposite<T> 
  {
    public:
      SumBB(T _x1, const GenBandMatrix<T2>& _m2, 
	  T _x3, const GenBandMatrix<T4>& _m4) :
	BandMatrixComposite<T>(_m2.stor()!=NoMajor ? _m2.stor() :
	    _m4.stor()!=NoMajor ? _m4.stor() : DiagMajor),
	x1(_x1),m2(&_m2),x3(_x3),m4(&_m4),
	itsnlo(max(m2->nlo(),m4->nlo())), itsnhi(max(m2->nhi(),m4->nhi()))
	{ 
	  TMVAssert(m2->rowsize() == m4->rowsize()); 
	  TMVAssert(m2->colsize() == m4->colsize()); 
	}
      inline size_t colsize() const { return m2->colsize(); }
      inline size_t rowsize() const { return m2->rowsize(); }
      inline int nlo() const { return itsnlo; }
      inline int nhi() const { return itsnhi; }
      inline T GetX1() const { return x1; }
      inline const GenBandMatrix<T2>& GetM2() const { return *m2; }
      inline T GetX3() const { return x3; }
      inline const GenBandMatrix<T4>& GetM4() const { return *m4; }
      inline void AssignTo(const BandMatrixView<T>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() == nlo());
	TMVAssert(m0.nhi() == nhi());
	if (m0.SameStorageAs(*m4)) {
	  m0 = *m4;
	  AddMM(x1,*m2,x3,m0);
	} else {
	  m0 = *m2;
	  AddMM(x3,*m4,x1,m0);
        }
      }
    private:
      T x1;
      const GenBandMatrix<T2>* m2;
      T x3;
      const GenBandMatrix<T4>* m4;
      int itsnlo,itsnhi;
  };

  template <class T> inline const BandMatrixView<T>& operator+=(
      const BandMatrixView<T>& m1, const GenBandMatrix<T>& m2) 
  { AddMM(T(1),m2,T(1),m1); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m1, const GenBandMatrix<T>& m2) 
  { AddMM(CT(1),m2,CT(1),m1); return m1; }

  template <class T> inline const BandMatrixView<T>& operator-=(
      const BandMatrixView<T>& m1, const GenBandMatrix<T>& m2) 
  { AddMM(T(-1),m2,T(1),m1); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m1, const GenBandMatrix<T>& m2) 
  { AddMM(CT(-1),m2,CT(1),m1); return m1; }

  //
  // BandMatrix * BandMatrix
  //

  template <class T, class T2, class T3, bool trans> class ProdBB : 
    public BandMatrixComposite<T>
  {
    public:
      ProdBB(T _x1, const GenBandMatrix<T2>& _m2,
	  const GenBandMatrix<T3>& _m3) :
	BandMatrixComposite<T>(_m2.isrm() ? RowMajor :
	    _m3.iscm() ? ColMajor :
	    (_m2.iscm() && _m3.isrm()) ? ColMajor : DiagMajor),
	x1(_x1), m2(&_m2), m3(&_m3),
	itsnlo(m2->nlo()+m3->nlo()), itsnhi(m2->nhi()+m3->nhi())
      { TMVAssert( m2->rowsize() == m3->colsize()); }
      inline size_t colsize() const { return m2->colsize(); }
      inline size_t rowsize() const { return m3->rowsize(); }
      inline int nlo() const { return itsnlo; }
      inline int nhi() const { return itsnhi; }
      inline T GetX1() const { return x1; }
      inline const GenBandMatrix<T2>& GetM2() const { return *m2; }
      inline const GenBandMatrix<T3>& GetM3() const { return *m3; }
      inline void AssignTo(const BandMatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() == nlo());
	TMVAssert(m0.nhi() == nhi());
	if (trans)
	  MultMM(x1,*m2,m3->Transpose(),T(0),m0);
	else
	  MultMM(x1,*m2,*m3,T(0),m0);
      }
    private:
      T x1;
      const GenBandMatrix<T2>*const m2;
      const GenBandMatrix<T3>*const m3;
      int itsnlo, itsnhi;
  };

  template <class T> inline const BandMatrixView<T>& operator*=(
      const BandMatrixView<T>& m1, const GenBandMatrix<T>& m2)
  { MultMM(T(1),m1,m2,T(0),m1); return m1; }

  template <class T> inline const BandMatrixView<CT>& operator*=(
      const BandMatrixView<CT>& m1, const GenBandMatrix<T>& m2)
  { MultMM(CT(1),m1,m2,CT(0),m1); return m1; }

  template <class T, class T2, class T3, bool trans>
    inline const BandMatrixView<T>& operator+=(
	const BandMatrixView<T>& m, const ProdBB<T,T2,T3,trans>& pmm)
    {
      MultMM(pmm.GetX1(),pmm.GetM2(), 
	trans ? pmm.GetM3().Transpose() : pmm.GetM3(),T(1),m); 
      return m; 
    }

  template <class T, class T2, class T3, bool trans>
    inline const BandMatrixView<T>& operator-=(
	const BandMatrixView<T>& m, const ProdBB<T,T2,T3,trans>& pmm)
    { 
      MultMM(-pmm.GetX1(),pmm.GetM2(),
	  trans ? pmm.GetM3().Transpose() : pmm.GetM3(),T(1),m); 
      return m; 
    }

  // First all the (M) op (B) possibilities

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXB

#define GENVECTOR1 GenVector
#define GENVECTOR2 GenVector
#define PRODXV1 ProdXV
#define PRODXV2 ProdXV

#define SUMMM SumMB
#define SUMMX SumBX
#define PRODXM ProdXB
#define PRODMV ProdBV
#define PRODMM ProdMB
#define QUOTVM QuotVB
#define RQUOTVM RQuotVB
#define QUOTMM QuotMB
#define RQUOTMM RQuotMB
#define TQUOTMM TransientQuotMB
#define TRQUOTMM TransientRQuotMB
#define QUOTXM QuotXB

#include "TMV_AuxVecComposite.h"
#include "TMV_AuxMatComposite1.h"
#include "TMV_AuxMatComposite3.h"

#include "TMV_AuxSumMM.h"
#include "TMV_AuxSumXM.h"
#include "TMV_AuxProdXM.h"
#include "TMV_AuxProdVM.h"
#include "TMV_AuxProdMM.h"
#include "TMV_AuxQuotVM.h"
#include "TMV_AuxQuotMM.h"
#include "TMV_AuxQuotXM.h"

  // Next (B) op (B)
  
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenBandMatrix
#define PRODXM1 ProdXB

#undef SUMMM
#undef PRODMM
#define SUMMM SumBB
#define PRODMM ProdBB

#include "TMV_AuxSumMM.h"
#include "TMV_AuxProdMM.h"
#include "TMV_AuxTQuotMM.h"

  // Finally (B) op (M)
  
#undef GENMATRIX2
#undef PRODXM2
#define GENMATRIX2 GenMatrix
#define PRODXM2 ProdXM

#undef SUMMM
#undef PRODMM
#undef TQUOTMM
#undef TRQUOTMM
#define SUMMM1 SumMB
#define SUMMM SumBM
#define PRODMM ProdBM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM

#include "TMV_AuxMatComposite2.h"

#include "TMV_AuxSumMM.h"
#include "TMV_AuxProdMM.h"
#include "TMV_AuxTQuotMM.h"

#undef GENVECTOR1
#undef GENVECTOR2
#undef PRODXV1
#undef PRODXV2

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#undef SUMMM
#undef SUMMM1
#undef SUMMX
#undef PRODXM
#undef PRODMV
#undef PRODMM
#undef QUOTVM
#undef RQUOTVM
#undef QUOTMM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM

}; // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif
