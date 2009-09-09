#ifndef TMV_SymMatrixArith_H
#define TMV_SymMatrixArith_H

#include "TMV_SymMatrix.h"

#define CT complex<T>
#define CCT ConjRef<complex<T> >
#define VCT VarConjRef<complex<T> >

namespace tmv {

  // y = alpha * A * x + beta * y
  template <class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y);

  // A = alpha * A
  template <class T> void MultXM(const T alpha, const SymMatrixView<T>& A);

  // B = alpha * A + beta * B
  template <class T, class Ta> void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A,
      const T beta, const SymMatrixView<T>& B);
  template <class T, class Ta> void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A,
      const T beta, const MatrixView<T>& B);

  // C = alpha * A * B + beta * C
  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const GenMatrix<Tb>& B, const T beta, const MatrixView<T>& C);
  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const GenSymMatrix<Tb>& B, const T beta, const MatrixView<T>& C)
  { MultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose()); }
  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const GenSymMatrix<Tb>& B, 
      const T beta, const MatrixView<T>& C);

  // A += alpha * (x^xT) (or x* if A is Herm)
  template <class T, class Tx> void Rank1Update(const T alpha,
      const GenVector<Tx>& x, const SymMatrixView<T>& A);
  // B += alpha * (A * AT) (or At if B is Herm)
  template <class T, class Ta> void RankKUpdate(const T alpha,
      const GenMatrix<Ta>& A, const SymMatrixView<T>& B);

  // These two don't have += forms: they must called explicitly
  // A += alpha * (x^y + y^x) (or x^y* + y^x* is A is Herm)
  template <class T, class Tx, class Ty> void Rank2Update(const T alpha,
      const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A);
  // C += alpha * (A * BT + B*AT) (or A*Bt + B*At if C is Herm)
  template <class T, class Ta, class Tb> void Rank2KUpdate(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const SymMatrixView<T>& C);

  template <class T> class SymMatrixComposite : 
    public GenSymMatrix<T>
  {
    public:

      SymMatrixComposite(SymType _sym, UpLoType _uplo, StorageType _stor) : 
	GenSymMatrix<T>(_sym,_uplo,_stor,NonConj), itsm(0) {}
      virtual ~SymMatrixComposite() { if (itsm) delete [] itsm; };
      virtual void AssignTo(const SymMatrixView<T>&) const = 0;
      inline const T* cptr() const
      { 
	if (!itsm) {
	  itsm = new T[size()*size()];
	  if (isherm())
	    HermMatrixViewOf(itsm,size(),uplo(),stor()) = *this;
	  else
	    SymMatrixViewOf(itsm,size(),uplo(),stor()) = *this;
	}
	return itsm;
      }
      inline int stepi() const { return isrm() ? size() : 1; }
      inline int stepj() const { return isrm() ? 1 : size(); }
      inline int diagstep() const { return size()+1; }
      using GenSymMatrix<T>::isrm;
      using GenSymMatrix<T>::isherm;
      using GenSymMatrix<T>::size;
      using GenSymMatrix<T>::sym;
      using GenSymMatrix<T>::uplo;
      using GenSymMatrix<T>::stor;
    private:
      mutable T* itsm;
  }; 

  template <class T, UpLoType U, StorageType S, class T2>
    inline SymMatrix<T,U,S>& operator+=(
	SymMatrix<T,U,S>& m1, const T2& x2) 
    { m1.View() += x2; return m1; }

  template <class T, UpLoType U, StorageType S, class T2>
    inline SymMatrix<T,U,S>& operator-=(
	SymMatrix<T,U,S>& m1, const T2& x2) 
    { m1.View() -= x2; return m1; }

  template <class T, UpLoType U, StorageType S, class T2>
    inline SymMatrix<T,U,S>& operator*=(
	SymMatrix<T,U,S>& m1, const T2& x2) 
    { m1.View() *= x2; return m1; }

  template <class T, UpLoType U, StorageType S, class T2>
    inline SymMatrix<T,U,S>& operator/=(
	SymMatrix<T,U,S>& m1, const T2& x2) 
    { m1.View() /= x2; return m1; }

  template <class T, UpLoType U, StorageType S, class T2>
    inline SymMatrix<T,U,S>& operator%=(
	SymMatrix<T,U,S>& m1, const T2& x2) 
    { m1.View() %= x2; return m1; }

  template <class T, UpLoType U, StorageType S, class T2>
    inline HermMatrix<T,U,S>& operator+=(
	HermMatrix<T,U,S>& m1, const T2& x2) 
    { m1.View() += x2; return m1; }

  template <class T, UpLoType U, StorageType S, class T2>
    inline HermMatrix<T,U,S>& operator-=(
	HermMatrix<T,U,S>& m1, const T2& x2) 
    { m1.View() -= x2; return m1; }

  template <class T, UpLoType U, StorageType S, class T2>
    inline HermMatrix<T,U,S>& operator*=(
	HermMatrix<T,U,S>& m1, const T2& x2) 
    { m1.View() *= x2; return m1; }

  template <class T, UpLoType U, StorageType S, class T2>
    inline HermMatrix<T,U,S>& operator/=(
	HermMatrix<T,U,S>& m1, const T2& x2) 
    { m1.View() /= x2; return m1; }

  template <class T, UpLoType U, StorageType S, class T2>
    inline HermMatrix<T,U,S>& operator%=(
	HermMatrix<T,U,S>& m1, const T2& x2) 
    { m1.View() %= x2; return m1; }


  //
  // SymMatrix + Scalar
  //

  template <class T, class T2> class SumSX : 
    public SymMatrixComposite<T>
  {
    public:
      SumSX(T _x1, const GenSymMatrix<T2>& _m2, T _x3) :
	SymMatrixComposite<T>(_m2.sym(),_m2.uplo(),BaseStorOf(_m2)),
	x1(_x1), m2(&_m2), x3(_x3)
      {
	TMVAssert(m2->IsSquare()); 
	TMVAssert(!m2->isherm() || 
	    (IMAG(x1)==RealType(T)(0) && IMAG(x3)==RealType(T)(0)) ) ;
      }
      inline size_t size() const { return m2->size(); }
      inline T GetX1() const { return *x3; }
      inline const GenSymMatrix<T2>& GetM2() const { return *m2; }
      inline T GetX3() const { return *x3; }
      inline void AssignTo(const SymMatrixView<T>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(m0.isherm() == this->isherm());
	if (x1 == T(1)) m0 = *m2;
	else m0 = x1*(*m2);
	m0.diag().AddToAll(x3);
      }
    private:
      const T x1;
      const GenSymMatrix<T2>* m2;
      const T x3;
  };

  // m+=x
  template <class T> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m1, T x2) 
  { 
    TMVAssert(!m1.isherm() || IMAG(x2) == RealType(T)(0));
    m1.diag().AddToAll(x2);
    return m1; 
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m1, T x2) 
  { 
    m1.diag().AddToAll(x2);
    return m1; 
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m1, CCT x2) 
  { 
    TMVAssert(!m1.isherm() || IMAG(x2) == T(0));
    m1.diag().AddToAll(CT(x2));
    return m1; 
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m1, VCT x2) 
  { 
    TMVAssert(!m1.isherm() || IMAG(x2) == T(0));
    m1.diag().AddToAll(CT(x2));
    return m1; 
  }

  // m-=x
  template <class T> inline const SymMatrixView<T>& operator-=(
      const SymMatrixView<T>& m1, T x2) 
  { 
    TMVAssert(!m1.isherm() || IMAG(x2) == RealType(T)(0));
    m1.diag().AddToAll(-x2);
    return m1; 
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m1, T x2) 
  { 
    m1.diag().AddToAll(-x2);
    return m1; 
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m1, CCT x2) 
  { 
    TMVAssert(!m1.isherm() || IMAG(x2) == T(0));
    m1.diag().AddToAll(-CT(x2));
    return m1; 
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m1, VCT x2) 
  { 
    TMVAssert(!m1.isherm() || IMAG(x2) == T(0));
    m1.diag().AddToAll(-CT(x2));
    return m1; 
  }


  //
  // Scalar * SymMatrix
  //

  template <class T, class T2> class ProdXS : 
    public SymMatrixComposite<T> 
  {
    public:
      ProdXS(const T _x1, const GenSymMatrix<T2>& _m2) :
	SymMatrixComposite<T>(_m2.sym(),_m2.uplo(), BaseStorOf(_m2)),
	x1(_x1), m2(&_m2) 
	{ TMVAssert(!m2->isherm() || IMAG(x1)==RealType(T)(0) ) ; }
      inline size_t size() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenSymMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const SymMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(m0.isherm() == this->isherm());
	TMVAssert(!m2->isherm() || IMAG(x1) == RealType(T)(0));
	m0 = *m2;
	MultXM(x1,m0);
      }
    private:
      const T x1;
      const GenSymMatrix<T2>*const m2;
  };

  // m*=x
  template <class T> inline const SymMatrixView<T>& operator*=(
      const SymMatrixView<T>& m1, T x2) 
  { 
    TMVAssert(!m1.isherm() || IMAG(x2) == RealType(T)(0));
    MultXM(x2,m1);
    return m1;
  }

  template <class T> inline const SymMatrixView<CT>& operator*=(
      const SymMatrixView<CT>& m1, T x2) 
  { 
    MultXM(CT(x2),m1);
    return m1;
  }

  template <class T> inline const SymMatrixView<CT>& operator*=(
      const SymMatrixView<CT>& m1, CCT x2) 
  { 
    TMVAssert(!m1.isherm() || IMAG(x2) == T(0));
    MultXM(CT(x2),m1);
    return m1;
  }

  template <class T> inline const SymMatrixView<CT>& operator*=(
      const SymMatrixView<CT>& m1, VCT x2) 
  { 
    TMVAssert(!m1.isherm() || IMAG(x2) == T(0));
    MultXM(CT(x2),m1);
    return m1;
  }

  // m/=x
  template <class T> inline const SymMatrixView<T>& operator/=(
      const SymMatrixView<T>& m1, T x2) 
  { 
    TMVAssert(!m1.isherm() || IMAG(x2) == RealType(T)(0));
    MultXM(RealType(T)(1)/x2,m1);
    return m1;
  }

  template <class T> inline const SymMatrixView<CT>& operator/=(
      const SymMatrixView<CT>& m1, T x2) 
  { 
    MultXM(CT(T(1)/x2),m1);
    return m1;
  }

  template <class T> inline const SymMatrixView<CT>& operator/=(
      const SymMatrixView<CT>& m1, CCT x2) 
  { 
    TMVAssert(!m1.isherm() || IMAG(x2) == T(0));
    MultXM(T(1)/CT(x2),m1);
    return m1;
  }

  template <class T> inline const SymMatrixView<CT>& operator/=(
      const SymMatrixView<CT>& m1, VCT x2) 
  { 
    TMVAssert(!m1.isherm() || IMAG(x2) == T(0));
    MultXM(T(1)/CT(x2),m1);
    return m1;
  }

  //
  // SymMatrix + SymMatrix
  //

  template <class T, class T2, class T4> class SumSS : 
    public SymMatrixComposite<T> 
  {
    public:
      SumSS(T _x1, const GenSymMatrix<T2>& _m2, 
	  T _x3, const GenSymMatrix<T4>& _m4) :
	SymMatrixComposite<T>(IsReal(T2()) ? _m4.sym() : _m2.sym(),
	    _m2.uplo(),BaseStorOf(_m2)),
	x1(_x1),m2(&_m2),x3(_x3),m4(&_m4)
	{ 
	  TMVAssert(m2->size() == m4->size()); 
	  // Technically, I should allow one of m2,m4 to be Sym
	  // and the other Herm, but this would require separating
	  // the GenSymMatrix class into a Sym and Herm version separately.
	  // This might be a good idea, but I think this is the only place
	  // the math would require it, and I suspect this is a rare
	  // thing to want to do, so we'll put it off for now.
	  // The workaround is to cast one of m2,m4 as a Matrix and
	  // then += the other.
	  TMVAssert(IsReal(T2()) || IsReal(T4()) || m2->sym() == m4->sym());
	}
      inline size_t size() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenSymMatrix<T2>& GetM2() const { return *m2; }
      inline T GetX3() const { return x3; }
      inline const GenSymMatrix<T4>& GetM4() const { return *m4; }
      inline void AssignTo(const SymMatrixView<T>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(m0.isherm() == this->isherm());
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
      const GenSymMatrix<T2>* m2;
      T x3;
      const GenSymMatrix<T4>* m4;
  };

  template <class T> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m1, const GenSymMatrix<T>& m2) 
  { AddMM(T(1),m2,T(1),m1); return m1; }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m1, const GenSymMatrix<T>& m2) 
  { AddMM(CT(1),m2,CT(1),m1); return m1; }

  template <class T> inline const SymMatrixView<T>& operator-=(
      const SymMatrixView<T>& m1, const GenSymMatrix<T>& m2) 
  { AddMM(T(-1),m2,T(1),m1); return m1; }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m1, const GenSymMatrix<T>& m2) 
  { AddMM(CT(-1),m2,CT(1),m1); return m1; }


  //
  // SymMatrix * SymMatrix
  //

  template <class T, class T2, class T3> class ProdSS : 
    public MatrixComposite<T>
  {
    public:
      ProdSS(T _x1, const GenSymMatrix<T2>& _m2,
	  const GenSymMatrix<T3>& _m3) :
	MatrixComposite<T>(BaseStorOf(_m2)), x1(_x1), m2(&_m2), m3(&_m3)
      { TMVAssert(m2->size() == m3->size()); }
      inline size_t colsize() const { return m2->size(); }
      inline size_t rowsize() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenSymMatrix<T2>& GetM2() const { return *m2; }
      inline const GenSymMatrix<T3>& GetM3() const { return *m3; }
      inline void AssignTo(const MatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	MultMM(x1,*m2,*m3,T(0),m0);
      }
    private:
      T x1;
      const GenSymMatrix<T2>*const m2;
      const GenSymMatrix<T3>*const m3;
  };

  template <class T, class T2, class T3> 
    inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdSS<T,T2,T3>& pmm)
  { 
    MultMM(pmm.GetX1(),pmm.GetM2(),pmm.GetM3(),T(1),m); 
    return m; 
  }

  template <class T, class T2, class T3> 
    inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdSS<T,T2,T3>& pmm)
  { 
    MultMM(-pmm.GetX1(),pmm.GetM2(),pmm.GetM3(),T(1),m); 
    return m; 
  }

  //
  // Vector ^ Vector
  //

  // m += (x*v^v)
  template <class T, class T2> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m0, const OProdVV<T,T2,T2>& opvv)
  {
    TMVAssert(m0.size() == opvv.colsize());
    TMVAssert(m0.size() == opvv.rowsize());
    if (m0.isherm()) 
    { TMVAssert(opvv.GetV2().SameAs(opvv.GetV3().Conjugate())); }
    else { TMVAssert(opvv.GetV2().SameAs(opvv.GetV3())); }
    Rank1Update(opvv.GetX1(), opvv.GetV2(), m0);
    return m0;
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m0, const OProdVV<T,T,T>& opvv)
  {
    TMVAssert(m0.size() == opvv.colsize());
    TMVAssert(m0.size() == opvv.rowsize());
    TMVAssert(opvv.GetV2().SameAs(opvv.GetV3()));
    Rank1Update(CT(opvv.GetX1()), opvv.GetV2(), m0);
    return m0;
  }

  // m -= (x*v^v)
  template <class T, class T2> inline const SymMatrixView<T>& operator-=(
      const SymMatrixView<T>& m0, const OProdVV<T,T2,T2>& opvv)
  {
    TMVAssert(m0.size() == opvv.colsize());
    TMVAssert(m0.size() == opvv.rowsize());
    if (m0.isherm()) 
    { TMVAssert(opvv.GetV2().SameAs(opvv.GetV3().Conjugate())); }
    else { TMVAssert(opvv.GetV2().SameAs(opvv.GetV3())); }
    Rank1Update(-opvv.GetX1(), opvv.GetV2(), m0);
    return m0;
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m0, const OProdVV<T,T,T>& opvv)
  {
    TMVAssert(m0.size() == opvv.colsize());
    TMVAssert(m0.size() == opvv.rowsize());
    TMVAssert(opvv.GetV2().SameAs(opvv.GetV3()));
    Rank1Update(CT(-opvv.GetX1()), opvv.GetV2(), m0);
    return m0;
  }

  //
  // Matrix ^ Matrix
  //

  // m += (x*m^m)
  template <class T, class T2> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m0, const ProdMM<T,T2,T2>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM2().SameAs(
	  m0.isherm() ? opmm.GetM3().Adjoint() : opmm.GetM3().Transpose()));
    RankKUpdate(opmm.GetX1(), opmm.GetM2(), m0);
    return m0;
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m0, const ProdMM<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM2().SameAs(opmm.GetM3().Transpose()));
    RankKUpdate(CT(opmm.GetX1()), opmm.GetM2(), m0);
    return m0;
  }

  // m -= (x*m^m)
  template <class T, class T2> inline const SymMatrixView<T>& operator-=(
	const SymMatrixView<T>& m0, const ProdMM<T,T2,T2>& opmm)
    {
      TMVAssert(m0.size() == opmm.colsize());
      TMVAssert(m0.size() == opmm.rowsize());
      TMVAssert(opmm.GetM2().SameAs(
	    m0.isherm() ? opmm.GetM3().Adjoint() : opmm.GetM3().Transpose()));
      RankKUpdate(-opmm.GetX1(), opmm.GetM2(), m0);
      return m0;
    }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m0, const ProdMM<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM2().SameAs(opmm.GetM3().Transpose()));
    RankKUpdate(CT(-opmm.GetX1()), opmm.GetM2(), m0);
    return m0;
  }

  // First all the (M) op (S) possibilities

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXS

#define GENVECTOR1 GenVector
#define GENVECTOR2 GenVector
#define PRODXV1 ProdXV
#define PRODXV2 ProdXV

#define SUMMM SumMS
#define SUMMX SumSX
#define PRODXM ProdXS
#define PRODMV ProdSV
#define PRODMM ProdMS
#define QUOTVM QuotVS
#define RQUOTVM RQuotVS
#define QUOTMM QuotMS
#define RQUOTMM RQuotMS
#define TQUOTMM TransientQuotMS
#define TRQUOTMM TransientRQuotMS
#define QUOTXM QuotXS

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

  // Next (S) op (S)
  
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenSymMatrix
#define PRODXM1 ProdXS

#undef SUMMM
#undef PRODMM
#define SUMMM SumSS
#define PRODMM ProdSS

#include "TMV_AuxSumMM.h"
#include "TMV_AuxProdMM.h"
#include "TMV_AuxTQuotMM.h"

  // Finally (S) op (M)
  
#undef GENMATRIX2
#undef PRODXM2
#define GENMATRIX2 GenMatrix
#define PRODXM2 ProdXM

#undef SUMMM
#undef PRODMM
#undef TQUOTMM
#undef TRQUOTMM
#define SUMMM1 SumMS
#define SUMMM SumSM
#define PRODMM ProdSM
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
