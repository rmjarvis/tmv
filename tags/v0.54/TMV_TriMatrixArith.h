#ifndef TMV_TriMatrixArith_H
#define TMV_TriMatrixArith_H

#include "TMV_TriMatrix.h"

#define CT complex<T>
#define CCT ConjRef<complex<T> >
#define VCT VarConjRef<complex<T> >

namespace tmv {

  // y = alpha * A * x + beta * y
  template <class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenUpperTriMatrix<Ta>& A, 
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y);
  template <class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenLowerTriMatrix<Ta>& A, 
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y);

  // A = alpha * A
  template <class T> void MultXM(const T alpha, const UpperTriMatrixView<T>& A);
  template <class T> inline void MultXM(const T alpha, const LowerTriMatrixView<T>& A)
  { MultXM(alpha,A.Transpose()); }

  // B = alpha * A + beta * B
  template <class T, class Ta> void AddMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const T beta, const UpperTriMatrixView<T>& B);
  template <class T, class Ta> inline void AddMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const T beta, const LowerTriMatrixView<T>& B)
  { AddMM(alpha,A.Transpose(),beta,B.Transpose()); }
  template <class T, class Ta> inline void AddMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const T beta, const MatrixView<T>& B)
  { 
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == B.rowsize());
    AddMM(alpha,A,beta,UpperTriMatrixViewOf(B)); 
    if (A.size() > 0) LowerTriMatrixViewOf(B).OffDiag() *= beta;
  }
  template <class T, class Ta> inline void AddMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const T beta, const MatrixView<T>& B)
  { AddMM(alpha,A.Transpose(),beta,B.Transpose()); }
  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const T beta, const GenUpperTriMatrix<Tb>& B, const MatrixView<T>& C);

  // C = alpha * A * B + beta * C
  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A, 
      const GenMatrix<Tb>& B, const T beta, const MatrixView<T>& C);
  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A, 
      const GenMatrix<Tb>& B, const T beta, const MatrixView<T>& C);
  template <class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const GenUpperTriMatrix<Tb>& B, const T beta, const MatrixView<T>& C)
  { MultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose()); }
  template <class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const GenLowerTriMatrix<Tb>& B, const T beta, const MatrixView<T>& C)
  { MultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose()); }

  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A, 
      const GenUpperTriMatrix<Tb>& B, 
      const T beta, const UpperTriMatrixView<T>& C);
  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A, 
      const GenLowerTriMatrix<Tb>& B, const T beta, const MatrixView<T>& C);
  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A, 
      const GenUpperTriMatrix<Tb>& B, const T beta, const MatrixView<T>& C);
  template <class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A, 
      const GenLowerTriMatrix<Tb>& B, 
      const T beta, const LowerTriMatrixView<T>& C)
  { MultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose()); }

  template <class T, class Ta> void ElementProd(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const UpperTriMatrixView<T>& B);
  template <class T, class Ta, class Tb> void AddElementProd(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenUpperTriMatrix<Tb>& B, const UpperTriMatrixView<T>& C);
  template <class T, class Ta> inline void ElementProd(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const LowerTriMatrixView<T>& B)
  { ElementProd(alpha,A.Transpose(),B.Transpose()); }
  template <class T, class Ta, class Tb> inline void AddElementProd(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenLowerTriMatrix<Tb>& B, const LowerTriMatrixView<T>& C)
  { AddElementProd(alpha,A.Transpose(),B.Transpose(),C.Transpose()); }
 

  //
  // First do everything which returns an UpperTriMatrixComposite
  // or only deals with an UpperTriMatrix, and not a LowerTriMatrix
  //
  
  template <class T> class UpperTriMatrixComposite : 
    public GenUpperTriMatrix<T>
  {
    public:

      UpperTriMatrixComposite(DiagType d, StorageType s) : 
	GenUpperTriMatrix<T>(d,s,NonConj), itsm(0) {}
      UpperTriMatrixComposite(const UpperTriMatrixComposite<T>& rhs) :
	GenUpperTriMatrix<T>(rhs), itsm(0) {}
      virtual ~UpperTriMatrixComposite() {}
      virtual void AssignTo(const UpperTriMatrixView<T>&) const = 0;
      inline const T* cptr() const
      { 
	if (!itsm.get()) {
	  itsm.reset(new T[size()*size()]);
	  UpperTriMatrixViewOf(itsm.get(),size(),dt(),stor()) = *this;
	}
	return itsm.get();
      }
      inline int stepi() const { return isrm() ? size() : 1; }
      inline int stepj() const { return isrm() ? 1 : size(); }
      inline int diagstep() const { return size()+1; }
      using GenUpperTriMatrix<T>::isrm;
      using GenUpperTriMatrix<T>::size;
      using GenUpperTriMatrix<T>::dt;
      using GenUpperTriMatrix<T>::stor;
    private:
      mutable auto_array<T> itsm;
  }; 

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline UpperTriMatrix<T,D,S,I>& operator+=(
	UpperTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() += x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline UpperTriMatrix<T,D,S,I>& operator-=(
	UpperTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() -= x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline UpperTriMatrix<T,D,S,I>& operator*=(
	UpperTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() *= x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline UpperTriMatrix<T,D,S,I>& operator/=(
	UpperTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() /= x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline UpperTriMatrix<T,D,S,I>& operator%=(
	UpperTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() %= x; return m; }

  template <class T, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m1, const GenUpperTriMatrix<T2>& m2) 
  {
    TMVAssert(m1.IsSquare());
    TMVAssert(m1.colsize() == m2.size());
    UpperTriMatrixViewOf(m1,NonUnitDiag) += m2; 
    return m1; 
  }

  template <class T, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m1, const GenUpperTriMatrix<T2>& m2) 
  {
    TMVAssert(m1.IsSquare());
    TMVAssert(m1.colsize() == m2.size());
    UpperTriMatrixViewOf(m1,NonUnitDiag) -= m2; 
    return m1; 
  }

  //
  // Scalar * TriMatrix
  //

  template <class T, class Tm> class ProdXU : 
    public UpperTriMatrixComposite<T> 
  {
    public:
      ProdXU(const T _x, const GenUpperTriMatrix<Tm>& _m) :
	UpperTriMatrixComposite<T>(NonUnitDiag, BaseStorOf(_m)),
	x(_x), m(&_m) {}
      inline size_t size() const { return m->size(); }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<Tm>& GetM() const { return *m; }
      inline void AssignTo(const UpperTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	m0 = *m;
	if (x != T(1)) MultXM(x,m0);
      }
    private:
      const T x;
      const GenUpperTriMatrix<Tm>*const m;
  };

  // m*=x
  template <class T> inline const UpperTriMatrixView<T>& operator*=(
      const UpperTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(x,m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator*=(
      const UpperTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(CT(x),m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator*=(
      const UpperTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(CT(x),m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator*=(
      const UpperTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(CT(x),m);
    return m;
  }

  // m/=x
  template <class T> inline const UpperTriMatrixView<T>& operator/=(
      const UpperTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(RealType(T)(1)/x,m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator/=(
      const UpperTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(CT(T(1)/x),m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator/=(
      const UpperTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(T(1)/CT(x),m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator/=(
      const UpperTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(T(1)/CT(x),m);
    return m;
  }

#define GENMATRIX GenUpperTriMatrix
#define PRODXM ProdXU
#include "TMV_AuxProdXM.h"

  //
  // TriMatrix + Scalar
  //

  template <class T, class Tm> class SumUX : 
    public UpperTriMatrixComposite<T>
  {
    // x1*m + x2
    public:
      SumUX(T _x1, const GenUpperTriMatrix<Tm>& _m, T _x2) :
	UpperTriMatrixComposite<T>(NonUnitDiag,BaseStorOf(_m)),
	x1(_x1), m(&_m), x2(_x2)
      { TMVAssert(m->IsSquare()); }
      inline size_t size() const { return m->size(); }
      inline T GetX1() const { return *x1; }
      inline const GenUpperTriMatrix<Tm>& GetM() const { return *m; }
      inline T GetX2() const { return *x2; }
      inline void AssignTo(const UpperTriMatrixView<T>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	if (x1 == T(1)) m0 = *m;
	else m0 = x1*(*m);
	m0.diag().AddToAll(x2);
      }
    private:
      const T x1;
      const GenUpperTriMatrix<Tm>* m;
      const T x2;
  };

  // m+=x
  template <class T> inline const UpperTriMatrixView<T>& operator+=(
      const UpperTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(x);
    return m; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(x);
    return m; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(CT(x));
    return m; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(CT(x));
    return m; 
  }

  // m-=x
  template <class T> inline const UpperTriMatrixView<T>& operator-=(
      const UpperTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-x);
    return m; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-x);
    return m; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-CT(x));
    return m; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-CT(x));
    return m; 
  }

#define SUMMX SumUX
#include "TMV_AuxSumMX.h"
#undef SUMMX


  //
  // TriMatrix + TriMatrix
  //

  template <class T, class T1, class T2> class SumUU : 
    public UpperTriMatrixComposite<T> 
  {
    public:
      SumUU(T _x1, const GenUpperTriMatrix<T1>& _m1, 
	  T _x2, const GenUpperTriMatrix<T2>& _m2) :
	UpperTriMatrixComposite<T>(NonUnitDiag, BaseStorOf(_m1)),
	x1(_x1),m1(&_m1),x2(_x2),m2(&_m2)
      { TMVAssert(m1->size() == m2->size()); }
      inline size_t size() const { return m1->size(); }
      inline T GetX1() const { return x1; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return *m1; }
      inline T GetX2() const { return x2; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const UpperTriMatrixView<T>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	if (m0.SameStorageAs(*m1)) {
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
      const T x1;
      const GenUpperTriMatrix<T1>* m1;
      const T x2;
      const GenUpperTriMatrix<T2>* m2;
  };

  template <class T> inline const UpperTriMatrixView<T>& operator+=(
      const UpperTriMatrixView<T>& m1, const GenUpperTriMatrix<T>& m2) 
  {
    TMVAssert(m1.size() == m2.size());
    AddMM(T(1),m2,T(1),m1); 
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m1, const GenUpperTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    AddMM(CT(1),m2,CT(1),m1); 
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<T>& operator-=(
      const UpperTriMatrixView<T>& m1, const GenUpperTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    AddMM(T(-1),m2,T(1),m1); 
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m1, const GenUpperTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    AddMM(CT(-1),m2,CT(1),m1); 
    return m1; 
  }

  template <class T, class T2> inline const UpperTriMatrixView<T>& operator+=(
      const UpperTriMatrixView<T>& m, const ProdXU<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(pxm.GetX(),pxm.GetM(),T(1),m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m, const ProdXU<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(CT(pxm.GetX()),pxm.GetM(),CT(1),m);
    return m;
  }

  template <class T, class T2> inline const UpperTriMatrixView<T>& operator-=(
      const UpperTriMatrixView<T>& m, const ProdXU<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(-pxm.GetX(),pxm.GetM(),T(1),m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m, const ProdXU<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(CT(-pxm.GetX()),pxm.GetM(),CT(1),m);
    return m;
  }

#define SUMMM SumUU
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXU
#include "TMV_AuxSumMM.h"
#undef SUMMM

  //
  // TriMatrix * TriMatrix
  //

  template <class T, class T1, class T2> class ProdUU : 
    public UpperTriMatrixComposite<T>
  {
    public:
      ProdUU(T _x, const GenUpperTriMatrix<T1>& _m1,
	  const GenUpperTriMatrix<T2>& _m2) :
	UpperTriMatrixComposite<T>(
	    (_m1.isunit()&&_m2.isunit()&&_x==T(1)) ?  UnitDiag : NonUnitDiag,
	    BaseStorOf(_m1)), x(_x), m1(&_m1), m2(&_m2)
      { TMVAssert(m1->size() == m2->size()); }
      inline size_t size() const { return m1->size(); }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return *m1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const UpperTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	MultMM(x,*m1,*m2,T(0),m0);
      }
    protected:
      T x;
      const GenUpperTriMatrix<T1>*const m1;
      const GenUpperTriMatrix<T2>*const m2;
  };

  template <class T> inline const UpperTriMatrixView<T>& operator*=(
      const UpperTriMatrixView<T>& m1, const GenUpperTriMatrix<T>& m2)
  { MultMM(T(1),m1,m2,T(0),m1); return m1; }

  template <class T> inline const UpperTriMatrixView<CT>& operator*=(
      const UpperTriMatrixView<CT>& m1, const GenUpperTriMatrix<T>& m2)
  { MultMM(CT(1),m1,m2,CT(0),m1); return m1; }

  template <class T, class T1, class T2> 
    inline const UpperTriMatrixView<T>& operator+=(
	const UpperTriMatrixView<T>& m, const ProdUU<T,T1,T2>& pmm)
    { MultMM(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); return m; }

  template <class T, class T1, class T2> 
    inline const UpperTriMatrixView<T>& operator-=(
	const UpperTriMatrixView<T>& m, const ProdUU<T,T1,T2>& pmm)
    { MultMM(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); return m; }

#define PRODMM ProdUU
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
  
  //
  // Scalar / TriMatrix
  //

  template <class T, class Tm> class QuotXU : 
    public UpperTriMatrixComposite<T> 
  {
    public:
      QuotXU(const T _x, const GenUpperTriMatrix<Tm>& _m) :
	UpperTriMatrixComposite<T>(_x==T(1) ? _m.dt() : NonUnitDiag,
	    BaseStorOf(_m)), x(_x), m(&_m)  {}
      inline size_t size() const { return m->size(); }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<Tm>& GetM() const { return *m; }
      inline void AssignTo(const UpperTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	if (x == T(0)) m0.Zero();
	else {
	  // Need temporary, since T, Tm are not the same.
	  // T = Tm overridden below
	  if (m0.isrm()) {
	    if (this->dt() == UnitDiag) {
	      TMVAssert(x == T(1));
	      UpperTriMatrix<Tm,UnitDiag,RowMajor> temp(
		  m0.colsize(),m0.rowsize());
	      m->Inverse(temp.View());
	      m0 = temp;
	    } else {
	      UpperTriMatrix<Tm,NonUnitDiag,RowMajor> temp(
		  m0.colsize(),m0.rowsize());
	      m->Inverse(temp.View());
	      if (x != T(1)) m0 = temp*x;
	      else m0 = temp;
	    }
	  } else {
	    if (this->dt() == UnitDiag) {
	      TMVAssert(x == T(1));
	      UpperTriMatrix<Tm,UnitDiag,ColMajor> temp(
		  m0.colsize(),m0.rowsize());
	      m->Inverse(temp.View());
	      m0 = temp;
	    } else {
	      UpperTriMatrix<Tm,NonUnitDiag,ColMajor> temp(
		  m0.colsize(),m0.rowsize());
	      m->Inverse(temp.View());
	      if (x != T(1)) m0 = temp*x;
	      else m0 = temp;
	    }
	  }
	}
      }
    private:
      const T x;
      const GenUpperTriMatrix<Tm>*const m;
  };

  template <class T> class QuotXU<T,T> : 
    public UpperTriMatrixComposite<T> 
    {
      public:
	QuotXU(const T _x, const GenUpperTriMatrix<T>& _m) :
	  UpperTriMatrixComposite<T>(_x==T(1) ? _m.dt() : NonUnitDiag,
	      BaseStorOf(_m)), x(_x), m(&_m)  {}
	inline size_t size() const { return m->size(); }
	inline T GetX() const { return x; }
	inline const GenUpperTriMatrix<T>& GetM() const { return *m; }
	inline void AssignTo(const UpperTriMatrixView<T>& m0) const
	{
	  TMVAssert(m0.size() == size());
	  TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	  if (x == T(0)) m0.Zero();
	  else {
	    m->Inverse(m0);
	    if (x != T(1)) m0 *= x;
	  }
	}
      private:
	const T x;
	const GenUpperTriMatrix<T>*const m;
    };

#define QUOTXM QuotXU
#include "TMV_AuxQuotXM.h"

  //
  // TriMatrix / % TriMatrix
  //

  template <class T, class T1, class T2> class QuotUU : 
    public UpperTriMatrixComposite<T>
  {
    public:
      QuotUU(const T _x, const GenUpperTriMatrix<T1>& _m1,
	  const GenUpperTriMatrix<T2>& _m2) :
	UpperTriMatrixComposite<T>(
	    _m1.dt()==_m2.dt() ? _m1.dt() : NonUnitDiag,
	    BaseStorOf(_m2)), x(_x), m1(&_m1), m2(&_m2)
      { TMVAssert( m1->size() == m2->size() ); }
      inline size_t size() const { return m1->size(); }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return *m1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const UpperTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || (m0.dt() == this->dt() && x==T(1)));
	if (m0.SameAs(*m1)) m2->LDivEq(m0);
	else if (m0.SameStorageAs(*m1) || m0.SameStorageAs(*m2)) {
	  if (m0.isrm()) {
	    if (m0.isunit()) {
	      UpperTriMatrix<T,UnitDiag,RowMajor> temp(m0.size());
	      m2->LDiv(*m1,temp.View());
	      m0 = temp;
	    } else {
	      UpperTriMatrix<T,NonUnitDiag,RowMajor> temp(m0.size());
	      m2->LDiv(*m1,temp.View());
	      m0 = temp;
	    }
	  } else {
	    if (m0.isunit()) {
	      UpperTriMatrix<T,UnitDiag,ColMajor> temp(m0.size());
	      m2->LDiv(*m1,temp.View());
	      m0 = temp;
	    } else {
	      UpperTriMatrix<T,NonUnitDiag,ColMajor> temp(m0.size());
	      m2->LDiv(*m1,temp.View());
	      m0 = temp;
	    }
	  } 
	} else {
	  m2->LDiv(*m1,m0);
	}
	if (x != T(1)) m0 *= x;
      }
    protected:
      const T x;
      const GenUpperTriMatrix<T1>*const m1;
      const GenUpperTriMatrix<T2>*const m2;
  };

  template <class T, class T1, class T2> class RQuotUU : 
    public UpperTriMatrixComposite<T>
  {
    public:
      RQuotUU(const T _x, const GenUpperTriMatrix<T1>& _m1,
	  const GenUpperTriMatrix<T2>& _m2) :
	UpperTriMatrixComposite<T>(
	    _m1.dt()==_m2.dt() ? _m1.dt() : NonUnitDiag, 
	    BaseStorOf(_m2)), x(_x), m1(&_m1), m2(&_m2)
      { TMVAssert( m1->size() == m2->size() ); }
      inline size_t size() const { return m1->size(); }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return *m1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const UpperTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || (m0.dt() == this->dt() && x==T(1)));
	if (m0.SameAs(*m1)) m2->RDivEq(m0);
	else if (m0.SameStorageAs(*m1) || m0.SameStorageAs(*m2)) {
	  if (m0.isrm()) {
	    if (m0.isunit()) {
	      UpperTriMatrix<T,UnitDiag,RowMajor> temp(m0.size());
	      m2->RDiv(*m1,temp.View());
	      m0 = temp;
	    } else {
	      UpperTriMatrix<T,NonUnitDiag,RowMajor> temp(m0.size());
	      m2->RDiv(*m1,temp.View());
	      m0 = temp;
	    }
	  } else {
	    if (m0.isunit()) {
	      UpperTriMatrix<T,UnitDiag,ColMajor> temp(m0.size());
	      m2->RDiv(*m1,temp.View());
	      m0 = temp;
	    } else {
	      UpperTriMatrix<T,NonUnitDiag,ColMajor> temp(m0.size());
	      m2->RDiv(*m1,temp.View());
	      m0 = temp;
	    }
	  } 
	} else {
	  m2->RDiv(*m1,m0);
	}
	if (x != T(1)) m0 *= x;
      }
    protected:
      const T x;
      const GenUpperTriMatrix<T1>*const m1;
      const GenUpperTriMatrix<T2>*const m2;
  };

  template <class T> inline const UpperTriMatrixView<T>& operator/=(
      const UpperTriMatrixView<T>& m1, const GenUpperTriMatrix<T>& m2)
  { m2.LDivEq(m1); return m1; }

  template <class T> inline const UpperTriMatrixView<CT>& operator/=(
      const UpperTriMatrixView<CT>& m1, const GenUpperTriMatrix<T>& m2)
  { m2.LDivEq(m1); return m1; }

  template <class T> inline const UpperTriMatrixView<T>& operator%=(
      const UpperTriMatrixView<T>& m1, const GenUpperTriMatrix<T>& m2)
  { m2.RDivEq(m1); return m1; }

  template <class T> inline const UpperTriMatrixView<CT>& operator%=(
      const UpperTriMatrixView<CT>& m1, const GenUpperTriMatrix<T>& m2)
  { m2.RDivEq(m1); return m1; }

#define QUOTMM QuotUU
#define RQUOTMM RQuotUU
#include "TMV_AuxQuotMM.h"
#include "TMV_AuxQuotMMa.h"
#undef QUOTMM
#undef RQUOTMM
    
#define PRODMV ProdUV
#define QUOTVM QuotVU
#define RQUOTVM RQuotVU

#include "TMV_AuxVecComposite.h"

#undef PRODMV
#undef QUOTVM
#undef RQUOTVM
#undef PRODXM
#undef GENMATRIX

#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenMatrix
#define PRODXM1 ProdXM
#define SUMMM SumMU
#define PRODMM ProdMU
#define QUOTMM QuotMU
#define RQUOTMM RQuotMU
#define TQUOTMM TransientQuotMU
#define TRQUOTMM TransientRQuotMU

#include "TMV_AuxMatComposite1.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef PRODMM
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM

#define GENMATRIX1 GenUpperTriMatrix
#define PRODXM1 ProdXU
#define GENMATRIX2 GenMatrix
#define PRODXM2 ProdXM
#define SUMMMa SumMU
#define SUMMM SumUM
#define PRODMM ProdUM
#define QUOTXM QuotXM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM

#include "TMV_AuxMatComposite2.h"
#include "TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#undef SUMMM
#undef SUMMMa
#undef PRODMM
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM

  //
  // Now the things which return a LowerTriMatrixComposite
  // or only deal with a LowerTriMatrix
  //

  template <class T> class LowerTriMatrixComposite : 
    public GenLowerTriMatrix<T>
  {
    public:

      LowerTriMatrixComposite(DiagType d, StorageType s) : 
	GenLowerTriMatrix<T>(d,s,NonConj), itsm(0) {}
      LowerTriMatrixComposite(const LowerTriMatrixComposite<T>& rhs) :
	GenLowerTriMatrix<T>(rhs), itsm(0) {}
      virtual ~LowerTriMatrixComposite() {}
      virtual void AssignTo(const LowerTriMatrixView<T>&) const = 0;
      inline const T* cptr() const
      { 
	if (!itsm.get()) {
	  itsm.reset(new T[size()*size()]);
	  LowerTriMatrixViewOf(itsm.get(),size(),dt(),stor()) = *this;
	}
	return itsm.get();
      }
      inline int stepi() const { return isrm() ? size() : 1; }
      inline int stepj() const { return isrm() ? 1 : size(); }
      inline int diagstep() const { return size()+1; }
      using GenLowerTriMatrix<T>::isrm;
      using GenLowerTriMatrix<T>::size;
      using GenLowerTriMatrix<T>::dt;
      using GenLowerTriMatrix<T>::stor;
    private:
      mutable auto_array<T> itsm;
  };

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline LowerTriMatrix<T,D,S,I>& operator+=(
	LowerTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() += x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline LowerTriMatrix<T,D,S,I>& operator-=(
	LowerTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() -= x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline LowerTriMatrix<T,D,S,I>& operator*=(
	LowerTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() *= x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline LowerTriMatrix<T,D,S,I>& operator/=(
	LowerTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() /= x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline LowerTriMatrix<T,D,S,I>& operator%=(
	LowerTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() %= x; return m; }

  template <class T, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m1, const GenLowerTriMatrix<T2>& m2) 
  {
    TMVAssert(m1.IsSquare());
    TMVAssert(m1.colsize() == m2.size());
    LowerTriMatrixViewOf(m1,NonUnitDiag) += m2; 
    return m1; 
  }

  template <class T, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m1, const GenLowerTriMatrix<T2>& m2) 
  { 
    TMVAssert(m1.IsSquare());
    TMVAssert(m1.colsize() == m2.size());
    LowerTriMatrixViewOf(m1,NonUnitDiag) -= m2; 
    return m1; 
  }

  //
  // Scalar * TriMatrix
  //

  template <class T, class Tm> class ProdXL : 
    public LowerTriMatrixComposite<T> 
  {
    public:
      ProdXL(const T _x, const GenLowerTriMatrix<Tm>& _m) :
	LowerTriMatrixComposite<T>(NonUnitDiag, BaseStorOf(_m)),
	x(_x), m(&_m) {}
      inline size_t size() const { return m->size(); }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<Tm>& GetM() const { return *m; }
      inline void AssignTo(const LowerTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	m0 = *m;
	if (x != T(1)) MultXM(x,m0);
      }
    private:
      const T x;
      const GenLowerTriMatrix<Tm>*const m;
  };

  // m*=x
  template <class T> inline const LowerTriMatrixView<T>& operator*=(
      const LowerTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(x,m);
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator*=(
      const LowerTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(CT(x),m);
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator*=(
      const LowerTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(CT(x),m);
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator*=(
      const LowerTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(CT(x),m);
    return m;
  }

  // m/=x
  template <class T> inline const LowerTriMatrixView<T>& operator/=(
      const LowerTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(RealType(T)(1)/x,m);
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator/=(
      const LowerTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(CT(T(1)/x),m);
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator/=(
      const LowerTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(T(1)/CT(x),m);
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator/=(
      const LowerTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(T(1)/CT(x),m);
    return m;
  }

#define GENMATRIX GenLowerTriMatrix
#define PRODXM ProdXL
#include "TMV_AuxProdXM.h"


  //
  // TriMatrix + Scalar
  //

  template <class T, class Tm> class SumLX : 
    public LowerTriMatrixComposite<T>
  {
    // x1*m + x2
    public:
      SumLX(T _x1, const GenLowerTriMatrix<Tm>& _m, T _x2) :
	LowerTriMatrixComposite<T>(NonUnitDiag,BaseStorOf(_m)),
	x1(_x1), m(&_m), x2(_x2)
      { TMVAssert(m->IsSquare()); }
      inline size_t size() const { return m->size(); }
      inline T GetX1() const { return *x1; }
      inline const GenLowerTriMatrix<Tm>& GetM() const { return *m; }
      inline T GetX2() const { return *x2; }
      inline void AssignTo(const LowerTriMatrixView<T>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	if (x1 == T(1)) m0 = *m;
	else m0 = x1*(*m);
	m0.diag().AddToAll(x2);
      }
    private:
      const T x1;
      const GenLowerTriMatrix<Tm>* m;
      const T x2;
  };

  // m+=x
  template <class T> inline const LowerTriMatrixView<T>& operator+=(
      const LowerTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(x);
    return m; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(x);
    return m; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(CT(x));
    return m; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(CT(x));
    return m; 
  }

  // m-=x
  template <class T> inline const LowerTriMatrixView<T>& operator-=(
      const LowerTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-x);
    return m; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-x);
    return m; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-CT(x));
    return m; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-CT(x));
    return m; 
  }

#define SUMMX SumLX
#include "TMV_AuxSumMX.h"
#undef SUMMX


  //
  // TriMatrix + TriMatrix
  // 

  template <class T, class T1, class T2> class SumLL : 
    public LowerTriMatrixComposite<T> 
  {
    public:
      SumLL(T _x1, const GenLowerTriMatrix<T1>& _m1, 
	  T _x2, const GenLowerTriMatrix<T2>& _m2) :
	LowerTriMatrixComposite<T>(NonUnitDiag, BaseStorOf(_m1)),
	x1(_x1),m1(&_m1),x2(_x2),m2(&_m2)
      { TMVAssert(m1->size() == m2->size()); }
      inline size_t size() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return *m1; }
      inline T GetX2() const { return x2; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const LowerTriMatrixView<T>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	if (m0.SameStorageAs(*m1)) {
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
      const GenLowerTriMatrix<T1>* m1;
      T x2;
      const GenLowerTriMatrix<T2>* m2;
  };

  template <class T> inline const LowerTriMatrixView<T>& operator+=(
      const LowerTriMatrixView<T>& m1, const GenLowerTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    AddMM(T(1),m2,T(1),m1); 
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m1, const GenLowerTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    AddMM(CT(1),m2,CT(1),m1); 
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<T>& operator-=(
      const LowerTriMatrixView<T>& m1, const GenLowerTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    AddMM(T(-1),m2,T(1),m1); 
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m1, const GenLowerTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    AddMM(CT(-1),m2,CT(1),m1); 
    return m1; 
  }

  template <class T, class T2> inline const LowerTriMatrixView<T>& operator+=(
      const LowerTriMatrixView<T>& m, const ProdXL<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(pxm.GetX(),pxm.GetM(),T(1),m);
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m, const ProdXL<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(CT(pxm.GetX()),pxm.GetM(),CT(1),m);
    return m;
  }

  template <class T, class T2> inline const LowerTriMatrixView<T>& operator-=(
      const LowerTriMatrixView<T>& m, const ProdXL<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(-pxm.GetX(),pxm.GetM(),T(1),m);
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m, const ProdXL<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(CT(-pxm.GetX()),pxm.GetM(),CT(1),m);
    return m;
  }

#define SUMMM SumLL
#define GENMATRIX1 GenLowerTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXL
#define PRODXM2 ProdXL
#include "TMV_AuxSumMM.h"
#undef SUMMM


  //
  // TriMatrix * TriMatrix
  //

  template <class T, class T1, class T2> class ProdLL : 
    public LowerTriMatrixComposite<T>
  {
    public:
      ProdLL(T _x, const GenLowerTriMatrix<T1>& _m1,
	  const GenLowerTriMatrix<T2>& _m2) :
	LowerTriMatrixComposite<T>(
	    (_m1.isunit()&&_m2.isunit()&&_x==T(1)) ?  UnitDiag : NonUnitDiag,
	    BaseStorOf(_m1)), x(_x), m1(&_m1), m2(&_m2)
      { TMVAssert(m2->size() == m2->size()); }
      inline size_t size() const { return m1->size(); }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return *m1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const LowerTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	MultMM(x,*m1,*m2,T(0),m0);
      }
    protected:
      T x;
      const GenLowerTriMatrix<T1>*const m1;
      const GenLowerTriMatrix<T2>*const m2;
  };

  template <class T> inline const LowerTriMatrixView<T>& operator*=(
      const LowerTriMatrixView<T>& m1, const GenLowerTriMatrix<T>& m2)
  { MultMM(T(1),m1,m2,T(0),m1); return m1; }

  template <class T> inline const LowerTriMatrixView<CT>& operator*=(
      const LowerTriMatrixView<CT>& m1, const GenLowerTriMatrix<T>& m2)
  { MultMM(CT(1),m1,m2,CT(0),m1); return m1; }

  template <class T, class T1, class T2> 
    inline const LowerTriMatrixView<T>& operator+=(
	const LowerTriMatrixView<T>& m, const ProdLL<T,T1,T2>& pmm)
    { MultMM(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); return m; }

  template <class T, class T1, class T2> 
    inline const LowerTriMatrixView<T>& operator-=(
	const LowerTriMatrixView<T>& m, const ProdLL<T,T1,T2>& pmm)
    { MultMM(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); return m; }

#define PRODMM ProdLL
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM

  //
  // Scalar / TriMatrix
  //

  template <class T, class Tm> class QuotXL : 
    public LowerTriMatrixComposite<T> 
  {
    public:
      QuotXL(const T _x, const GenLowerTriMatrix<Tm>& _m) :
	LowerTriMatrixComposite<T>(_x==T(1) ? _m.dt() : NonUnitDiag,
	    BaseStorOf(_m)), x(_x), m(&_m)  {}
      inline size_t size() const { return m->size(); }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<Tm>& GetM() const { return *m; }
      inline void AssignTo(const LowerTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	if (x == T(0)) m0.Zero();
	else {
	  // Need temporary, since T, Tm are not the same.
	  // T = Tm overridden below
	  if (m0.isrm()) {
	    if (this->dt() == UnitDiag) {
	      TMVAssert(x == T(1));
	      LowerTriMatrix<Tm,UnitDiag,RowMajor> temp(
		  m0.colsize(),m0.rowsize());
	      m->Inverse(temp.View());
	      m0 = temp;
	    } else {
	      LowerTriMatrix<Tm,NonUnitDiag,RowMajor> temp(
		  m0.colsize(),m0.rowsize());
	      m->Inverse(temp.View());
	      if (x != T(1)) m0 = temp*x;
	      else m0 = temp;
	    }
	  } else {
	    if (this->dt() == UnitDiag) {
	      TMVAssert(x == T(1));
	      LowerTriMatrix<Tm,UnitDiag,ColMajor> temp(
		  m0.colsize(),m0.rowsize());
	      m->Inverse(temp.View());
	      m0 = temp;
	    } else {
	      LowerTriMatrix<Tm,NonUnitDiag,ColMajor> temp(
		  m0.colsize(),m0.rowsize());
	      m->Inverse(temp.View());
	      if (x != T(1)) m0 = temp*x;
	      else m0 = temp;
	    }
	  }
	}
      }
    private:
      const T x;
      const GenLowerTriMatrix<Tm>*const m;
  };

  template <class T> class QuotXL<T,T> : 
    public LowerTriMatrixComposite<T> 
    {
      public:
	QuotXL(const T _x, const GenLowerTriMatrix<T>& _m) :
	  LowerTriMatrixComposite<T>(_x==T(1) ? _m.dt() : NonUnitDiag,
	      BaseStorOf(_m)), x(_x), m(&_m)  {}
	inline size_t size() const { return m->size(); }
	inline T GetX() const { return x; }
	inline const GenLowerTriMatrix<T>& GetM() const { return *m; }
	inline void AssignTo(const LowerTriMatrixView<T>& m0) const
	{
	  TMVAssert(m0.size() == size());
	  TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	  if (x == T(0)) m0.Zero();
	  else {
	    m->Inverse(m0);
	    if (x != T(1)) m0 *= x;
	  }
	}
      private:
	const T x;
	const GenLowerTriMatrix<T>*const m;
    };

#define QUOTXM QuotXL
#include "TMV_AuxQuotXM.h"

  // 
  // TriMatrix / TriMatrix
  //

  template <class T, class T1, class T2> class QuotLL : 
    public LowerTriMatrixComposite<T>
  {
    public:
      QuotLL(const T _x, const GenLowerTriMatrix<T1>& _m1,
	  const GenLowerTriMatrix<T2>& _m2) :
	LowerTriMatrixComposite<T>(
	    _m1.dt()==_m2.dt() ? _m1.dt() : NonUnitDiag,
	    BaseStorOf(_m2)), x(_x), m1(&_m1), m2(&_m2)
      { TMVAssert( m1->size() == m2->size() ); }
      inline size_t size() const { return m1->size(); }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return *m1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const LowerTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || (m0.dt() == this->dt() && x==T(1)));
	if (m0.SameAs(*m1)) m2->LDivEq(m0);
	else if (m0.SameStorageAs(*m1) || m0.SameStorageAs(*m2)) {
	  if (m0.isrm()) {
	    if (m0.isunit()) {
	      LowerTriMatrix<T,UnitDiag,RowMajor> temp(m0.size());
	      m2->LDiv(*m1,temp.View());
	      m0 = temp;
	    } else {
	      LowerTriMatrix<T,NonUnitDiag,RowMajor> temp(m0.size());
	      m2->LDiv(*m1,temp.View());
	      m0 = temp;
	    }
	  } else {
	    if (m0.isunit()) {
	      LowerTriMatrix<T,UnitDiag,ColMajor> temp(m0.size());
	      m2->LDiv(*m1,temp.View());
	      m0 = temp;
	    } else {
	      LowerTriMatrix<T,NonUnitDiag,ColMajor> temp(m0.size());
	      m2->LDiv(*m1,temp.View());
	      m0 = temp;
	    }
	  } 
	} else {
	  m2->LDiv(*m1,m0);
	}
	if (x != T(1)) m0 *= x;
      }
    protected:
      const T x;
      const GenLowerTriMatrix<T1>*const m1;
      const GenLowerTriMatrix<T2>*const m2;
  };

  template <class T, class T1, class T2> class RQuotLL : 
    public LowerTriMatrixComposite<T>
  {
    public:
      RQuotLL(const T _x, const GenLowerTriMatrix<T1>& _m1,
	  const GenLowerTriMatrix<T2>& _m2) :
	LowerTriMatrixComposite<T>(
	    _m1.dt()==_m2.dt() ? _m1.dt() : NonUnitDiag, 
	    BaseStorOf(_m2)), x(_x), m1(&_m1), m2(&_m2)
      { TMVAssert( m1->size() == m2->size() ); }
      inline size_t size() const { return m1->size(); }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return *m1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const LowerTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || (m0.dt() == this->dt() && x==T(1)));
	if (m0.SameAs(*m1)) m2->RDivEq(m0);
	else if (m0.SameStorageAs(*m1) || m0.SameStorageAs(*m2)) {
	  if (m0.isrm()) {
	    if (m0.isunit()) {
	      LowerTriMatrix<T,UnitDiag,RowMajor> temp(m0.size());
	      m2->RDiv(*m1,temp.View());
	      m0 = temp;
	    } else {
	      LowerTriMatrix<T,NonUnitDiag,RowMajor> temp(m0.size());
	      m2->RDiv(*m1,temp.View());
	      m0 = temp;
	    }
	  } else {
	    if (m0.isunit()) {
	      LowerTriMatrix<T,UnitDiag,ColMajor> temp(m0.size());
	      m2->RDiv(*m1,temp.View());
	      m0 = temp;
	    } else {
	      LowerTriMatrix<T,NonUnitDiag,ColMajor> temp(m0.size());
	      m2->RDiv(*m1,temp.View());
	      m0 = temp;
	    }
	  } 
	} else {
	  m2->RDiv(*m1,m0);
	}
	if (x != T(1)) m0 *= x;
      }
    protected:
      const T x;
      const GenLowerTriMatrix<T1>*const m1;
      const GenLowerTriMatrix<T2>*const m2;
  };

  template <class T> inline const LowerTriMatrixView<T>& operator/=(
      const LowerTriMatrixView<T>& m1, const GenLowerTriMatrix<T>& m2)
  { m2.LDivEq(m1); return m1; }

  template <class T> inline const LowerTriMatrixView<CT>& operator/=(
      const LowerTriMatrixView<CT>& m1, const GenLowerTriMatrix<T>& m2)
  { m2.LDivEq(m1); return m1; }

  template <class T> inline const LowerTriMatrixView<T>& operator%=(
      const LowerTriMatrixView<T>& m1, const GenLowerTriMatrix<T>& m2)
  { m2.RDivEq(m1); return m1; }

  template <class T> inline const LowerTriMatrixView<CT>& operator%=(
      const LowerTriMatrixView<CT>& m1, const GenLowerTriMatrix<T>& m2)
  { m2.RDivEq(m1); return m1; }

#define QUOTMM QuotLL
#define RQUOTMM RQuotLL
#include "TMV_AuxQuotMM.h"
#include "TMV_AuxQuotMMa.h"
#undef QUOTMM
#undef RQUOTMM

#define PRODMV ProdLV
#define QUOTVM QuotVL
#define RQUOTVM RQuotVL

#include "TMV_AuxVecComposite.h"

#undef PRODMV
#undef QUOTVM
#undef RQUOTVM
#undef PRODXM
#undef GENMATRIX

#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenMatrix
#define PRODXM1 ProdXM
#define SUMMM SumML
#define PRODMM ProdML
#define QUOTMM QuotML
#define RQUOTMM RQuotML
#define TQUOTMM TransientQuotML
#define TRQUOTMM TransientRQuotML

#include "TMV_AuxMatComposite1.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef PRODMM
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM

#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL
#define GENMATRIX2 GenMatrix
#define PRODXM2 ProdXM
#define SUMMMa SumML
#define SUMMM SumLM
#define PRODMM ProdLM
#define QUOTXM QuotXM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM

#include "TMV_AuxMatComposite2.h"
#include "TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#undef SUMMM
#undef SUMMMa
#undef PRODMM
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM

  //
  // Finally, the things which mix Upper and Lower
  //
 
  template <class T, class T1, class T2> class SumUL : 
    public MatrixComposite<T> 
  {
    public:
      SumUL(T _x1, const GenUpperTriMatrix<T1>& _m1, 
	  T _x2, const GenLowerTriMatrix<T2>& _m2) :
	MatrixComposite<T>(BaseStorOf(_m1)),
	x1(_x1),m1(&_m1),x2(_x2),m2(&_m2)
      {  TMVAssert(m1->size() == m2->size()); }
      inline size_t colsize() const { return m1->size(); }
      inline size_t rowsize() const { return m1->size(); }
      inline T GetX1() const { return x1; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return *m1; }
      inline T GetX2() const { return x2; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const MatrixView<T>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());

	AddMM(x2,*m2,x1,*m1,m0);
      }
    private:
      T x1;
      const GenUpperTriMatrix<T1>* m1;
      T x2;
      const GenLowerTriMatrix<T2>* m2;
  };

#define SUMMM SumUL
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXL
#include "TMV_AuxSumMM.h"
#undef SUMMM

  template <class T, class T1, class T2> class ProdUL : 
    public MatrixComposite<T>
  {
    public:
      ProdUL(T _x, const GenUpperTriMatrix<T1>& _m1,
	  const GenLowerTriMatrix<T2>& _m2) :
	MatrixComposite<T>(BaseStorOf(_m1)), x(_x), m1(&_m1), m2(&_m2)
      { TMVAssert( m1->size() == m2->size()); }
      inline size_t colsize() const { return m1->size(); }
      inline size_t rowsize() const { return m1->size(); }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return *m1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const MatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	MultMM(x,*m1,*m2,T(0),m0);
      }
    protected:
      T x;
      const GenUpperTriMatrix<T1>*const m1;
      const GenLowerTriMatrix<T2>*const m2;
  };

  template <class T, class T1, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdUL<T,T1,T2>& pmm)
  { MultMM(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); return m; }

  template <class T, class T1, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdUL<T,T1,T2>& pmm)
  { MultMM(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); return m; }

#define PRODMM ProdUL
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
  
#define TQUOTMM TransientQuotML
#define TRQUOTMM TransientRQuotML
#define QUOTXM QuotXL
#include "TMV_AuxTQuotMM.h"
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

  template <class T, class T1, class T2> class SumLU : 
    public MatrixComposite<T> 
  {
    public:
      SumLU(T _x1, const GenLowerTriMatrix<T1>& _m1, 
	  T _x2, const GenUpperTriMatrix<T2>& _m2) :
	MatrixComposite<T>(BaseStorOf(_m1)),
	x1(_x1),m1(&_m1),x2(_x2),m2(&_m2)
      {  TMVAssert(m1->size() == m2->size()); }
      inline size_t colsize() const { return m1->size(); }
      inline size_t rowsize() const { return m1->size(); }
      inline T GetX1() const { return x1; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return *m1; }
      inline T GetX2() const { return x2; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const MatrixView<T>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());

	AddMM(x1,*m1,x2,*m2,m0);
      }
    private:
      T x1;
      const GenLowerTriMatrix<T1>* m1;
      T x2;
      const GenUpperTriMatrix<T2>* m2;
  };

#define SUMMM SumLU
#define GENMATRIX1 GenLowerTriMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXL
#define PRODXM2 ProdXU
#include "TMV_AuxSumMM.h"
#undef SUMMM

  template <class T, class T1, class T2> class ProdLU : 
    public MatrixComposite<T>
  {
    public:
      ProdLU(T _x, const GenLowerTriMatrix<T1>& _m1,
	  const GenUpperTriMatrix<T2>& _m2) :
	MatrixComposite<T>(BaseStorOf(_m1)), x(_x), m1(&_m1), m2(&_m2)
      { TMVAssert( m1->size() == m2->size()); }
      inline size_t colsize() const { return m1->size(); }
      inline size_t rowsize() const { return m1->size(); }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return *m1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const MatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	MultMM(x,*m1,*m2,T(0),m0);
      }
    protected:
      T x;
      const GenLowerTriMatrix<T1>*const m1;
      const GenUpperTriMatrix<T2>*const m2;
  };

  template <class T, class T1, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdLU<T,T1,T2>& pmm)
  { MultMM(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); return m; }

  template <class T, class T1, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdLU<T,T1,T2>& pmm)
  { MultMM(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); return m; }

#define PRODMM ProdLU
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM

#define TQUOTMM TransientQuotMU
#define TRQUOTMM TransientRQuotMU
#define QUOTXM QuotXU
#include "TMV_AuxTQuotMM.h"
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif
