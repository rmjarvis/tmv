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
  template <class T> inline void MultXM(
      const T alpha, const SymMatrixView<T>& A)
  {
    TMVAssert(A.issym() || IMAG(alpha) == RealType(T)(0));
    MultXM(alpha,A.UpperTri());
  }

  // B = alpha * A + beta * B
  template <class T, class Ta> inline void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A,
      const T beta, const SymMatrixView<T>& B)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(IsReal(Ta()) || A.sym() == B.sym());
    TMVAssert(B.issym() || IMAG(alpha) == RealType(T)(0));
    TMVAssert(B.issym() || IMAG(beta) == RealType(T)(0));
    AddMM(alpha,A.UpperTri(),beta,B.UpperTri());
  }
  template <class T, class Ta> inline void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A,
      const T beta, const MatrixView<T>& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == B.rowsize());
    AddMM(alpha,A.UpperTri(),beta,UpperTriMatrixViewOf(B));
    if (A.size() > 0)
      AddMM(alpha,A.LowerTri().OffDiag(),beta,
	  LowerTriMatrixViewOf(B).OffDiag());
  }
 

  // C = alpha * A * B + beta * C
  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const GenMatrix<Tb>& B, const T beta, const MatrixView<T>& C);
  template <class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const GenSymMatrix<Tb>& B, const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.rowsize() == B.size());
    TMVAssert(C.rowsize() == B.size());
    TMVAssert(A.colsize() == C.colsize());
    MultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose()); 
  }
  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const GenSymMatrix<Tb>& B, 
      const T beta, const MatrixView<T>& C);

  // A = beta*A + alpha * (x^xT) (or x* if A is Herm)
  // beta = 0 or 1
  template <class T, class Tx> void Rank1Update(const T alpha,
      const GenVector<Tx>& x, const int beta, const SymMatrixView<T>& A);
  // B = beta*B + alpha * (A * AT) (or At if B is Herm)
  // beta = 0 or 1
  template <class T, class Ta> void RankKUpdate(const T alpha,
      const GenMatrix<Ta>& A, const int beta, const SymMatrixView<T>& B);
  template <class T, class Ta> void RankKUpdate(const T alpha,
      const GenLowerTriMatrix<Ta>& A, const int beta, 
      const SymMatrixView<T>& B);
  template <class T, class Ta> void RankKUpdate(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const int beta,
      const SymMatrixView<T>& B);

  // These two don't have += forms: they must called explicitly
  // A += beta*A + alpha * (x^y + y^x) (or x^y* + y^x* is A is Herm)
  // beta = 0 or 1
  template <class T, class Tx, class Ty> void Rank2Update(const T alpha,
      const GenVector<Tx>& x, const GenVector<Ty>& y, const int beta,
      const SymMatrixView<T>& A);
  // C = beta*C + alpha * (A * BT + B*AT) (or A*Bt + B*At if C is Herm)
  // beta = 0 or 1
  template <class T, class Ta, class Tb> void Rank2KUpdate(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B, const int beta,
      const SymMatrixView<T>& C);

  // C = beta*C + alpha * A * B
  // beta = 0 or 1
  // This also needs to be called explicitly.
  // This is to prevent the programmer from doing this accidentally
  // when alpha * A * B is not necessarily symmetrix/hermitian.
  template <class T, class Ta, class Tb> void SymMultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B, const int beta,
      const SymMatrixView<T>& C);

  template <class T, class Ta> inline void ElementProd(
      const T alpha, const GenSymMatrix<Ta>& A, const SymMatrixView<T>& B)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(IsReal(T()) || A.sym() == B.sym());
    TMVAssert(B.issym() || IMAG(alpha) == RealType(T)(0));
    ElementProd(alpha,A.UpperTri(),B.UpperTri());
  }
  template <class T, class Ta, class Tb> inline void AddElementProd(
      const T alpha, const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
      const SymMatrixView<T>& C)
  {
    TMVAssert(A.size() == C.size());
    TMVAssert(B.size() == C.size());
    TMVAssert(IsReal(T()) || A.sym() == C.sym());
    TMVAssert(IsReal(T()) || B.sym() == C.sym());
    TMVAssert(C.issym() || IMAG(alpha) == RealType(T)(0));
    AddElementProd(alpha,A.UpperTri(),B.UpperTri(),C.UpperTri());
  }
    
  template <class T> class SymMatrixComposite : 
    public GenSymMatrix<T>
  {
    public:

      SymMatrixComposite(SymType _sym, UpLoType _uplo, StorageType _stor) : 
	GenSymMatrix<T>(_sym,_uplo,_stor,NonConj), itsm(0) {}
      SymMatrixComposite(const SymMatrixComposite<T>& rhs) :
	GenSymMatrix<T>(rhs), itsm(0) {}
      virtual ~SymMatrixComposite() {}
      virtual void AssignTo(const SymMatrixView<T>&) const = 0;
      inline const T* cptr() const
      { 
	if (!itsm.get()) {
	  itsm.reset(new T[size()*size()]);
	  if (isherm())
	    HermMatrixViewOf(itsm.get(),size(),uplo(),stor()) = *this;
	  else
	    SymMatrixViewOf(itsm.get(),size(),uplo(),stor()) = *this;
	}
	return itsm.get();
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
      mutable auto_array<T> itsm;
  }; 

  template <class T, UpLoType U, StorageType S, class Tx>
    inline SymMatrix<T,U,S>& operator+=(
	SymMatrix<T,U,S>& m, const Tx& x) 
    { m.View() += x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline SymMatrix<T,U,S>& operator-=(
	SymMatrix<T,U,S>& m, const Tx& x) 
    { m.View() -= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline SymMatrix<T,U,S>& operator*=(
	SymMatrix<T,U,S>& m, const Tx& x) 
    { m.View() *= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline SymMatrix<T,U,S>& operator/=(
	SymMatrix<T,U,S>& m, const Tx& x) 
    { m.View() /= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline SymMatrix<T,U,S>& operator%=(
	SymMatrix<T,U,S>& m, const Tx& x) 
    { m.View() %= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline HermMatrix<T,U,S>& operator+=(
	HermMatrix<T,U,S>& m, const Tx& x) 
    { m.View() += x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline HermMatrix<T,U,S>& operator-=(
	HermMatrix<T,U,S>& m, const Tx& x) 
    { m.View() -= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline HermMatrix<T,U,S>& operator*=(
	HermMatrix<T,U,S>& m, const Tx& x) 
    { m.View() *= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline HermMatrix<T,U,S>& operator/=(
	HermMatrix<T,U,S>& m, const Tx& x) 
    { m.View() /= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline HermMatrix<T,U,S>& operator%=(
	HermMatrix<T,U,S>& m, const Tx& x) 
    { m.View() %= x; return m; }


  //
  // Scalar * SymMatrix
  //

  template <class T, class Tm> class ProdXS : 
    public MatrixComposite<T> 
  {
    public:
      ProdXS(const T _x, const GenSymMatrix<Tm>& _m) :
	MatrixComposite<T>(BaseStorOf(_m)), x(_x), m(&_m) {}
      inline size_t colsize() const { return m->size(); }
      inline size_t rowsize() const { return m->size(); }
      inline T GetX() const { return x; }
      inline const GenSymMatrix<Tm>& GetM() const { return *m; }
      inline void AssignTo(const MatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == m->size());
	TMVAssert(m0.rowsize() == m->size());
	m0 = *m;
	if (x != T(1)) MultXM(x,m0);
      }
      inline void AssignTo(const SymMatrixView<T>& m0) const
      {
	TMVAssert(m->issym() || IMAG(x)==RealType(T)(0) ); 
	TMVAssert(m0.size() == m->size());
	TMVAssert(IsReal(Tm()) || m0.issym() == m->issym());
	m0 = *m;
	if (x != T(1)) MultXM(x,m0);
      }
    private:
      const T x;
      const GenSymMatrix<Tm>*const m;
  };

  // m*=x
  template <class T> inline const SymMatrixView<T>& operator*=(
      const SymMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(x,m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator*=(
      const SymMatrixView<CT>& m, T x) 
  { 
    MultXM(CT(x),m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator*=(
      const SymMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(CT(x),m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator*=(
      const SymMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(CT(x),m);
    return m;
  }

  // m/=x
  template <class T> inline const SymMatrixView<T>& operator/=(
      const SymMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(RealType(T)(1)/x,m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator/=(
      const SymMatrixView<CT>& m, T x) 
  { 
    MultXM(CT(T(1)/x),m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator/=(
      const SymMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(T(1)/CT(x),m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator/=(
      const SymMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(T(1)/CT(x),m);
    return m;
  }

#define GENMATRIX GenSymMatrix
#define PRODXM ProdXS
#include "TMV_AuxProdXM.h"

  //
  // SymMatrix + Scalar
  //

  template <class T, class Tm> class SumSX : 
    public MatrixComposite<T>
  {
    public:
      SumSX(T _x1, const GenSymMatrix<Tm>& _m, T _x2) :
	MatrixComposite<T>(BaseStorOf(_m)),
	x1(_x1), m(&_m), x2(_x2) {}
      inline size_t colsize() const { return m->size(); }
      inline size_t rowsize() const { return m->size(); }
      inline T GetX1() const { return x1; }
      inline const GenSymMatrix<Tm>& GetM() const { return *m; }
      inline T GetX2() const { return x2; }
      inline void AssignTo(const MatrixView<T>& m0) const
      { 
	TMVAssert(m0.colsize() == m->size());
	TMVAssert(m0.rowsize() == m->size());
	if (x1 == T(1)) m0 = *m;
	else m0 = x1*(*m);
	m0.diag().AddToAll(x2);
      } 
      inline void AssignTo(const SymMatrixView<T>& m0) const
      { 
	TMVAssert(IsReal(Tm()) || m0.issym() == m->issym());
	TMVAssert(m->issym() || IMAG(x1)==RealType(T)(0));
	TMVAssert(m->issym() || IMAG(x2)==RealType(T)(0));
	TMVAssert(m0.size() == m->size());
	if (x1 == T(1)) m0 = *m;
	else m0 = x1*(*m);
	m0.diag().AddToAll(x2);
      }
    private:
      const T x1;
      const GenSymMatrix<Tm>* m;
      const T x2;
  };

  // m+=x
  template <class T> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    m.diag().AddToAll(x);
    return m; 
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m, T x) 
  { 
    m.diag().AddToAll(x);
    return m; 
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == T(0));
    m.diag().AddToAll(CT(x));
    return m; 
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == T(0));
    m.diag().AddToAll(CT(x));
    return m; 
  }

  // m-=x
  template <class T> inline const SymMatrixView<T>& operator-=(
      const SymMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    m.diag().AddToAll(-x);
    return m; 
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m, T x) 
  { 
    m.diag().AddToAll(-x);
    return m; 
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == T(0));
    m.diag().AddToAll(-CT(x));
    return m; 
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == T(0));
    m.diag().AddToAll(-CT(x));
    return m; 
  }

#define SUMMX SumSX
#include "TMV_AuxSumMX.h"
#undef SUMMX 

  //
  // SymMatrix + SymMatrix
  //

  template <class T, class T1, class T2> class SumSS : 
    public MatrixComposite<T> 
  {
    public:
      SumSS(T _x1, const GenSymMatrix<T1>& _m1, 
	  T _x2, const GenSymMatrix<T2>& _m2) :
	MatrixComposite<T>(BaseStorOf(_m1)),
	x1(_x1),m1(&_m1),x2(_x2),m2(&_m2)
      { TMVAssert(m1->size() == m2->size()); }
      inline size_t colsize() const { return m1->size(); }
      inline size_t rowsize() const { return m1->size(); }
      inline T GetX1() const { return x1; }
      inline const GenSymMatrix<T1>& GetM1() const { return *m1; }
      inline T GetX2() const { return x2; }
      inline const GenSymMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const MatrixView<T>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
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
      inline void AssignTo(const SymMatrixView<T>& m0) const
      { 
	TMVAssert(IsReal(T1()) || IsReal(T2()) || m1->sym() == m2->sym());
	TMVAssert(m0.size() == m1->size());
	TMVAssert(IsReal(T1()) || m0.issym() == m1->issym());
	TMVAssert(IsReal(T2()) || m0.issym() == m2->issym());
	TMVAssert(m0.issym() || IMAG(x1) == RealType(T)(0));
	TMVAssert(m0.issym() || IMAG(x2) == RealType(T)(0));
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
      const GenSymMatrix<T1>* m1;
      T x2;
      const GenSymMatrix<T2>* m2;
  };

  template <class T> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m1, const GenSymMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(IsReal(T()) || m1.sym() == m2.sym());
    AddMM(T(1),m2,T(1),m1); return m1; 
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m1, const GenSymMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    AddMM(CT(1),m2,CT(1),m1); return m1; 
  }

  template <class T> inline const SymMatrixView<T>& operator-=(
      const SymMatrixView<T>& m1, const GenSymMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(IsReal(T()) || m1.sym() == m2.sym());
    AddMM(T(-1),m2,T(1),m1); return m1; 
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m1, const GenSymMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    AddMM(CT(-1),m2,CT(1),m1); return m1; 
  }

  template <class T, class T2> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m, const ProdXS<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(IsReal(T2()) || m.issym() == pxm.GetM().issym());
    TMVAssert(m.issym() || IMAG(pxm.GetX()) == RealType(T)(0));
    AddMM(pxm.GetX(),pxm.GetM(),T(1),m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m, const ProdXS<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(CT(pxm.GetX()),pxm.GetM(),CT(1),m);
    return m;
  }

  template <class T, class T2> inline const SymMatrixView<T>& operator-=(
      const SymMatrixView<T>& m, const ProdXS<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(IsReal(T2()) || m.issym() == pxm.GetM().issym());
    TMVAssert(m.issym() || IMAG(pxm.GetX()) == RealType(T)(0));
    AddMM(-pxm.GetX(),pxm.GetM(),T(1),m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m, const ProdXS<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(CT(-pxm.GetX()),pxm.GetM(),CT(1),m);
    return m;
  }

#define SUMMM SumSS
#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXS
#define PRODXM2 ProdXS
#include "TMV_AuxSumMM.h"
#undef SUMMM


  //
  // SymMatrix * SymMatrix
  //

  template <class T, class T1, class T2> class ProdSS : 
    public MatrixComposite<T>
  {
    public:
      ProdSS(T _x, const GenSymMatrix<T1>& _m1,
	  const GenSymMatrix<T2>& _m2) :
	MatrixComposite<T>(BaseStorOf(_m1)), x(_x), m1(&_m1), m2(&_m2)
      { TMVAssert(m1->size() == m2->size()); }
      inline size_t colsize() const { return m1->size(); }
      inline size_t rowsize() const { return m1->size(); }
      inline T GetX() const { return x; }
      inline const GenSymMatrix<T1>& GetM1() const { return *m1; }
      inline const GenSymMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const MatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	MultMM(x,*m1,*m2,T(0),m0);
      }
    private:
      T x;
      const GenSymMatrix<T1>*const m1;
      const GenSymMatrix<T2>*const m2;
  };

  template <class T, class T1, class T2> 
    inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdSS<T,T1,T2>& pmm)
    {
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      MultMM(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdSS<T,T1,T2>& pmm)
    { 
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      MultMM(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); 
      return m; 
    }

#define PRODMM ProdSS
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM


  //
  // Vector ^ Vector
  //

  // m += (x*v^v)
  template <class T, class Tv> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m0, const OProdVV<T,Tv,Tv>& opvv)
  {
    TMVAssert(m0.size() == opvv.colsize());
    TMVAssert(m0.size() == opvv.rowsize());
    TMVAssert(opvv.GetV1().SameAs(
	  m0.isherm() ? opvv.GetV2().Conjugate() : opvv.GetV2().View()));
    Rank1Update(opvv.GetX(), opvv.GetV1(), 1, m0);
    return m0;
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m0, const OProdVV<T,T,T>& opvv)
  {
    TMVAssert(m0.size() == opvv.colsize());
    TMVAssert(m0.size() == opvv.rowsize());
    TMVAssert(opvv.GetV1().SameAs(opvv.GetV2()));
    Rank1Update(CT(opvv.GetX()), opvv.GetV1(), 1, m0);
    return m0;
  }

  // m -= (x*v^v)
  template <class T, class Tv> inline const SymMatrixView<T>& operator-=(
      const SymMatrixView<T>& m0, const OProdVV<T,Tv,Tv>& opvv)
  {
    TMVAssert(m0.size() == opvv.colsize());
    TMVAssert(m0.size() == opvv.rowsize());
    TMVAssert(opvv.GetV1().SameAs(
	  m0.isherm() ? opvv.GetV2().Conjugate() : opvv.GetV2().View()));
    Rank1Update(-opvv.GetX(), opvv.GetV1(), 1, m0);
    return m0;
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m0, const OProdVV<T,T,T>& opvv)
  {
    TMVAssert(m0.size() == opvv.colsize());
    TMVAssert(m0.size() == opvv.rowsize());
    TMVAssert(opvv.GetV1().SameAs(opvv.GetV2()));
    Rank1Update(CT(-opvv.GetX()), opvv.GetV1(), 1, m0);
    return m0;
  }

  //
  // Matrix * Matrix.Transpose()
  //

  // m += (x*m*mt)
  template <class T, class Tm> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m0, const ProdMM<T,Tm,Tm>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(
	  m0.isherm() ? opmm.GetM2().Adjoint() : opmm.GetM2().Transpose()));
    RankKUpdate(opmm.GetX(), opmm.GetM1(), 1, m0);
    return m0;
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m0, const ProdMM<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(opmm.GetM2().Transpose()));
    RankKUpdate(CT(opmm.GetX()), opmm.GetM1(), 1, m0);
    return m0;
  }

  template <class T, class Tm> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m0, const ProdLU<T,Tm,Tm>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(
	  m0.isherm() ? opmm.GetM2().Adjoint() : opmm.GetM2().Transpose()));
    RankKUpdate(opmm.GetX(), opmm.GetM1(), 1, m0);
    return m0;
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m0, const ProdLU<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(opmm.GetM2().Transpose()));
    RankKUpdate(CT(opmm.GetX()), opmm.GetM1(), 1, m0);
    return m0;
  }

  template <class T, class Tm> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m0, const ProdUL<T,Tm,Tm>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(
	  m0.isherm() ? opmm.GetM2().Adjoint() : opmm.GetM2().Transpose()));
    RankKUpdate(opmm.GetX(), opmm.GetM1(), 1, m0);
    return m0;
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m0, const ProdUL<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(opmm.GetM2().Transpose()));
    RankKUpdate(CT(opmm.GetX()), opmm.GetM1(), 1, m0);
    return m0;
  }

  // m -= (x*m*mt)
  template <class T, class Tm> inline const SymMatrixView<T>& operator-=(
	const SymMatrixView<T>& m0, const ProdMM<T,Tm,Tm>& opmm)
    {
      TMVAssert(m0.size() == opmm.colsize());
      TMVAssert(m0.size() == opmm.rowsize());
      TMVAssert(opmm.GetM1().SameAs(
	    m0.isherm() ? opmm.GetM2().Adjoint() : opmm.GetM2().Transpose()));
      RankKUpdate(-opmm.GetX(), opmm.GetM1(), 1, m0);
      return m0;
    }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m0, const ProdMM<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(opmm.GetM2().Transpose()));
    RankKUpdate(CT(-opmm.GetX()), opmm.GetM1(), 1, m0);
    return m0;
  }

  template <class T, class Tm> inline const SymMatrixView<T>& operator-=(
	const SymMatrixView<T>& m0, const ProdUL<T,Tm,Tm>& opmm)
    {
      TMVAssert(m0.size() == opmm.colsize());
      TMVAssert(m0.size() == opmm.rowsize());
      TMVAssert(opmm.GetM1().SameAs(
	    m0.isherm() ? opmm.GetM2().Adjoint() : opmm.GetM2().Transpose()));
      RankKUpdate(-opmm.GetX(), opmm.GetM1(), 1, m0);
      return m0;
    }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m0, const ProdUL<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(opmm.GetM2().Transpose()));
    RankKUpdate(CT(-opmm.GetX()), opmm.GetM1(), 1, m0);
    return m0;
  }

  template <class T, class Tm> inline const SymMatrixView<T>& operator-=(
	const SymMatrixView<T>& m0, const ProdLU<T,Tm,Tm>& opmm)
    {
      TMVAssert(m0.size() == opmm.colsize());
      TMVAssert(m0.size() == opmm.rowsize());
      TMVAssert(opmm.GetM1().SameAs(
	    m0.isherm() ? opmm.GetM1().Adjoint() : opmm.GetM2().Transpose()));
      RankKUpdate(-opmm.GetX(), opmm.GetM1(), 1, m0);
      return m0;
    }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m0, const ProdLU<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(opmm.GetM2().Transpose()));
    RankKUpdate(CT(-opmm.GetX()), opmm.GetM1(), 1, m0);
    return m0;
  }



#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenMatrix
#define PRODXM1 ProdXM
#define SUMMM SumMS
#define PRODMM ProdMS
#define QUOTXM QuotXS
#define QUOTMM QuotMS
#define RQUOTMM RQuotMS
#define TQUOTMM TransientQuotMS
#define TRQUOTMM TransientRQuotMS

#include "TMV_AuxMatComposite1.h"
#include "TMV_AuxMatComposite3.h"

#define PRODMV ProdSV
#define QUOTVM QuotVS
#define RQUOTVM RQuotVS
#include "TMV_AuxVecComposite.h"
#undef PRODMV
#undef QUOTVM
#undef RQUOTVM
#undef PRODXM
#undef GENMATRIX

  // S/S -> Matrix.  Use TransientQuot
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenSymMatrix
#define PRODXM1 ProdXS
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
#define SUMMMa SumMS
#define SUMMM SumSM
#define PRODMM ProdSM
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
