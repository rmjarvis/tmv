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
  template <class T> void MultXM(const T alpha, const LowerTriMatrixView<T>& A)
  { MultXM(alpha,A.Transpose()); }

  // B = alpha * A + beta * B
  template <class T, class Ta> void AddMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const T beta, const UpperTriMatrixView<T>& B);
  template <class T, class Ta> void AddMM(
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

  template <class T> class UpperTriMatrixComposite : 
    public GenUpperTriMatrix<T>
  {
    public:

      UpperTriMatrixComposite(DiagType d, StorageType s) : 
	GenUpperTriMatrix<T>(d,s,NonConj), itsm(0) {}
      virtual ~UpperTriMatrixComposite() { if (itsm) delete [] itsm; };
      virtual void AssignTo(const UpperTriMatrixView<T>&) const = 0;
      inline const T* cptr() const
      { 
	if (!itsm) {
	  itsm = new T[size()*size()];
	  UpperTriMatrixViewOf(itsm,size(),dt(),stor()) = *this;
	}
	return itsm;
      }
      inline int stepi() const { return isrm() ? size() : 1; }
      inline int stepj() const { return isrm() ? 1 : size(); }
      inline int diagstep() const { return size()+1; }
      using GenUpperTriMatrix<T>::isrm;
      using GenUpperTriMatrix<T>::size;
      using GenUpperTriMatrix<T>::dt;
      using GenUpperTriMatrix<T>::stor;
    private:
      mutable T* itsm;
  }; 

  template <class T> class LowerTriMatrixComposite : 
    public GenLowerTriMatrix<T>
  {
    public:

      LowerTriMatrixComposite(DiagType d, StorageType s) : 
	GenLowerTriMatrix<T>(d,s,NonConj), itsm(0) {}
      virtual ~LowerTriMatrixComposite() { if (itsm) delete itsm; };
      virtual void AssignTo(const LowerTriMatrixView<T>&) const = 0;
      inline const T* cptr() const
      { 
	if (!itsm) {
	  itsm = new T[size()*size()];
	  LowerTriMatrixViewOf(itsm,size(),dt(),stor()) = *this;
	}
	return itsm;
      }
      inline int stepi() const { return isrm() ? size() : 1; }
      inline int stepj() const { return isrm() ? 1 : size(); }
      inline int diagstep() const { return size()+1; }
      using GenLowerTriMatrix<T>::isrm;
      using GenLowerTriMatrix<T>::size;
      using GenLowerTriMatrix<T>::dt;
      using GenLowerTriMatrix<T>::stor;
    private:
      mutable T* itsm;
  };

  template <class T, DiagType D, StorageType S, class T2>
    inline UpperTriMatrix<T,D,S>& operator+=(
	UpperTriMatrix<T,D,S>& m1, const T2& x2) 
    { m1.View() += x2; return m1; }

  template <class T, DiagType D, StorageType S, class T2>
    inline UpperTriMatrix<T,D,S>& operator-=(
	UpperTriMatrix<T,D,S>& m1, const T2& x2) 
    { m1.View() -= x2; return m1; }

  template <class T, DiagType D, StorageType S, class T2>
    inline UpperTriMatrix<T,D,S>& operator*=(
	UpperTriMatrix<T,D,S>& m1, const T2& x2) 
    { m1.View() *= x2; return m1; }

  template <class T, DiagType D, StorageType S, class T2>
    inline UpperTriMatrix<T,D,S>& operator/=(
	UpperTriMatrix<T,D,S>& m1, const T2& x2) 
    { m1.View() /= x2; return m1; }

  template <class T, DiagType D, StorageType S, class T2>
    inline UpperTriMatrix<T,D,S>& operator%=(
	UpperTriMatrix<T,D,S>& m1, const T2& x2) 
    { m1.View() %= x2; return m1; }

  template <class T, DiagType D, StorageType S, class T2>
    inline LowerTriMatrix<T,D,S>& operator+=(
	LowerTriMatrix<T,D,S>& m1, const T2& x2) 
    { m1.View() += x2; return m1; }

  template <class T, DiagType D, StorageType S, class T2>
    inline LowerTriMatrix<T,D,S>& operator-=(
	LowerTriMatrix<T,D,S>& m1, const T2& x2) 
    { m1.View() -= x2; return m1; }

  template <class T, DiagType D, StorageType S, class T2>
    inline LowerTriMatrix<T,D,S>& operator*=(
	LowerTriMatrix<T,D,S>& m1, const T2& x2) 
    { m1.View() *= x2; return m1; }

  template <class T, DiagType D, StorageType S, class T2>
    inline LowerTriMatrix<T,D,S>& operator/=(
	LowerTriMatrix<T,D,S>& m1, const T2& x2) 
    { m1.View() /= x2; return m1; }

  template <class T, DiagType D, StorageType S, class T2>
    inline LowerTriMatrix<T,D,S>& operator%=(
	LowerTriMatrix<T,D,S>& m1, const T2& x2) 
    { m1.View() %= x2; return m1; }

  template <class T, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m1, const GenUpperTriMatrix<T2>& m2) 
  { UpperTriMatrixViewOf(m1,NonUnitDiag) += m2; return m1; }

  template <class T, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m1, const GenUpperTriMatrix<T2>& m2) 
  { UpperTriMatrixViewOf(m1,NonUnitDiag) -= m2; return m1; }

  template <class T, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m1, const GenLowerTriMatrix<T2>& m2) 
  { LowerTriMatrixViewOf(m1,NonUnitDiag) += m2; return m1; }

  template <class T, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m1, const GenLowerTriMatrix<T2>& m2) 
  { LowerTriMatrixViewOf(m1,NonUnitDiag) -= m2; return m1; }

  //
  // TriMatrix + Scalar
  //

  template <class T, class T2> class SumUX : 
    public UpperTriMatrixComposite<T>
  {
    public:
      SumUX(T _x1, const GenUpperTriMatrix<T2>& _m2, T _x3) :
	UpperTriMatrixComposite<T>(NonUnitDiag,BaseStorOf(_m2)),
	x1(_x1), m2(&_m2), x3(_x3)
      { TMVAssert(m2->IsSquare()); }
      inline size_t size() const { return m2->size(); }
      inline T GetX1() const { return *x3; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline T GetX3() const { return *x3; }
      inline void AssignTo(const UpperTriMatrixView<T>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	if (x1 == T(1)) m0 = *m2;
	else m0 = x1*(*m2);
	m0.diag().AddToAll(x3);
      }
    private:
      const T x1;
      const GenUpperTriMatrix<T2>* m2;
      const T x3;
  };

  // m+=x
  template <class T> inline const UpperTriMatrixView<T>& operator+=(
      const UpperTriMatrixView<T>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(x2);
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(x2);
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m1, CCT x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(CT(x2));
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m1, VCT x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(CT(x2));
    return m1; 
  }

  // m-=x
  template <class T> inline const UpperTriMatrixView<T>& operator-=(
      const UpperTriMatrixView<T>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(-x2);
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(-x2);
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m1, CCT x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(-CT(x2));
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m1, VCT x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(-CT(x2));
    return m1; 
  }

  template <class T, class T2> class SumLX : 
    public LowerTriMatrixComposite<T>
  {
    public:
      SumLX(T _x1, const GenLowerTriMatrix<T2>& _m2, T _x3) :
	LowerTriMatrixComposite<T>(NonUnitDiag,BaseStorOf(_m2)),
	x1(_x1), m2(&_m2), x3(_x3)
      { TMVAssert(m2->IsSquare()); }
      inline size_t size() const { return m2->size(); }
      inline T GetX1() const { return *x3; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline T GetX3() const { return *x3; }
      inline void AssignTo(const LowerTriMatrixView<T>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	if (x1 == T(1)) m0 = *m2;
	else m0 = x1*(*m2);
	m0.diag().AddToAll(x3);
      }
    private:
      const T x1;
      const GenLowerTriMatrix<T2>* m2;
      const T x3;
  };

  // m+=x
  template <class T> inline const LowerTriMatrixView<T>& operator+=(
      const LowerTriMatrixView<T>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(x2);
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(x2);
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m1, CCT x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(CT(x2));
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m1, VCT x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(CT(x2));
    return m1; 
  }

  // m-=x
  template <class T> inline const LowerTriMatrixView<T>& operator-=(
      const LowerTriMatrixView<T>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(-x2);
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(-x2);
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m1, CCT x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(-CT(x2));
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m1, VCT x2) 
  { 
    TMVAssert(!m1.isunit());
    m1.diag().AddToAll(-CT(x2));
    return m1; 
  }


  //
  // Scalar * TriMatrix
  //

  template <class T, class T2> class ProdXU : 
    public UpperTriMatrixComposite<T> 
  {
    public:
      ProdXU(const T _x1, const GenUpperTriMatrix<T2>& _m2) :
	UpperTriMatrixComposite<T>(NonUnitDiag, BaseStorOf(_m2)),
	x1(_x1), m2(&_m2) {}
      inline size_t size() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const UpperTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	m0 = *m2;
	MultXM(x1,m0);
      }
    private:
      const T x1;
      const GenUpperTriMatrix<T2>*const m2;
  };

  // m*=x
  template <class T> inline const UpperTriMatrixView<T>& operator*=(
      const UpperTriMatrixView<T>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(x2,m1);
    return m1;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator*=(
      const UpperTriMatrixView<CT>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(CT(x2),m1);
    return m1;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator*=(
      const UpperTriMatrixView<CT>& m1, CCT x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(CT(x2),m1);
    return m1;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator*=(
      const UpperTriMatrixView<CT>& m1, VCT x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(CT(x2),m1);
    return m1;
  }

  // m/=x
  template <class T> inline const UpperTriMatrixView<T>& operator/=(
      const UpperTriMatrixView<T>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(RealType(T)(1)/x2,m1);
    return m1;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator/=(
      const UpperTriMatrixView<CT>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(CT(T(1)/x2),m1);
    return m1;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator/=(
      const UpperTriMatrixView<CT>& m1, CCT x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(T(1)/CT(x2),m1);
    return m1;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator/=(
      const UpperTriMatrixView<CT>& m1, VCT x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(T(1)/CT(x2),m1);
    return m1;
  }

  template <class T, class T2> class ProdXL : 
    public LowerTriMatrixComposite<T> 
  {
    public:
      ProdXL(const T _x1, const GenLowerTriMatrix<T2>& _m2) :
	LowerTriMatrixComposite<T>(NonUnitDiag, BaseStorOf(_m2)),
	x1(_x1), m2(&_m2) {}
      inline size_t size() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const LowerTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	m0 = *m2;
	MultXM(x1,m0);
      }
    private:
      const T x1;
      const GenLowerTriMatrix<T2>*const m2;
  };

  // m*=x
  template <class T> inline const LowerTriMatrixView<T>& operator*=(
      const LowerTriMatrixView<T>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(x2,m1);
    return m1;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator*=(
      const LowerTriMatrixView<CT>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(CT(x2),m1);
    return m1;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator*=(
      const LowerTriMatrixView<CT>& m1, CCT x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(CT(x2),m1);
    return m1;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator*=(
      const LowerTriMatrixView<CT>& m1, VCT x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(CT(x2),m1);
    return m1;
  }

  // m/=x
  template <class T> inline const LowerTriMatrixView<T>& operator/=(
      const LowerTriMatrixView<T>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(RealType(T)(1)/x2,m1);
    return m1;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator/=(
      const LowerTriMatrixView<CT>& m1, T x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(CT(T(1)/x2),m1);
    return m1;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator/=(
      const LowerTriMatrixView<CT>& m1, CCT x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(T(1)/CT(x2),m1);
    return m1;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator/=(
      const LowerTriMatrixView<CT>& m1, VCT x2) 
  { 
    TMVAssert(!m1.isunit());
    MultXM(T(1)/CT(x2),m1);
    return m1;
  }

  //
  // TriMatrix + TriMatrix
  //

  template <class T, class T2, class T4> class SumUU : 
    public UpperTriMatrixComposite<T> 
  {
    public:
      SumUU(T _x1, const GenUpperTriMatrix<T2>& _m2, 
	  T _x3, const GenUpperTriMatrix<T4>& _m4) :
	UpperTriMatrixComposite<T>(NonUnitDiag, BaseStorOf(_m2)),
	x1(_x1),m2(&_m2),x3(_x3),m4(&_m4)
      { TMVAssert(m2->size() == m4->size()); }
      inline size_t size() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline T GetX3() const { return x3; }
      inline const GenUpperTriMatrix<T4>& GetM4() const { return *m4; }
      inline void AssignTo(const UpperTriMatrixView<T>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
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
      const GenUpperTriMatrix<T2>* m2;
      T x3;
      const GenUpperTriMatrix<T4>* m4;
  };

  template <class T> inline const UpperTriMatrixView<T>& operator+=(
      const UpperTriMatrixView<T>& m1, const GenUpperTriMatrix<T>& m2) 
  { AddMM(T(1),m2,T(1),m1); return m1; }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m1, const GenUpperTriMatrix<T>& m2) 
  { AddMM(CT(1),m2,CT(1),m1); return m1; }

  template <class T> inline const UpperTriMatrixView<T>& operator-=(
      const UpperTriMatrixView<T>& m1, const GenUpperTriMatrix<T>& m2) 
  { AddMM(T(-1),m2,T(1),m1); return m1; }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m1, const GenUpperTriMatrix<T>& m2) 
  { AddMM(CT(-1),m2,CT(1),m1); return m1; }

  template <class T, class T2, class T4> class SumLL : 
    public LowerTriMatrixComposite<T> 
  {
    public:
      SumLL(T _x1, const GenLowerTriMatrix<T2>& _m2, 
	  T _x3, const GenLowerTriMatrix<T4>& _m4) :
	LowerTriMatrixComposite<T>(NonUnitDiag, BaseStorOf(_m2)),
	x1(_x1),m2(&_m2),x3(_x3),m4(&_m4)
      { TMVAssert(m2->size() == m4->size()); }
      inline size_t size() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline T GetX3() const { return x3; }
      inline const GenLowerTriMatrix<T4>& GetM4() const { return *m4; }
      inline void AssignTo(const LowerTriMatrixView<T>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
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
      const GenLowerTriMatrix<T2>* m2;
      T x3;
      const GenLowerTriMatrix<T4>* m4;
  };

  template <class T> inline const LowerTriMatrixView<T>& operator+=(
      const LowerTriMatrixView<T>& m1, const GenLowerTriMatrix<T>& m2) 
  { AddMM(T(1),m2,T(1),m1); return m1; }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m1, const GenLowerTriMatrix<T>& m2) 
  { AddMM(CT(1),m2,CT(1),m1); return m1; }

  template <class T> inline const LowerTriMatrixView<T>& operator-=(
      const LowerTriMatrixView<T>& m1, const GenLowerTriMatrix<T>& m2) 
  { AddMM(T(-1),m2,T(1),m1); return m1; }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m1, const GenLowerTriMatrix<T>& m2) 
  { AddMM(CT(-1),m2,CT(1),m1); return m1; }

  template <class T, class T2, class T4> class SumUL : 
    public MatrixComposite<T> 
  {
    public:
      SumUL(T _x1, const GenUpperTriMatrix<T2>& _m2, 
	  T _x3, const GenLowerTriMatrix<T4>& _m4) :
	MatrixComposite<T>(BaseStorOf(_m2)),
	x1(_x1),m2(&_m2),x3(_x3),m4(&_m4)
      {  TMVAssert(m2->size() == m4->size()); }
      inline size_t colsize() const { return m2->size(); }
      inline size_t rowsize() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline T GetX3() const { return x3; }
      inline const GenLowerTriMatrix<T4>& GetM4() const { return *m4; }
      inline void AssignTo(const MatrixView<T>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());

	AddMM(x3,*m4,x1,*m2,m0);
      }
    private:
      T x1;
      const GenUpperTriMatrix<T2>* m2;
      T x3;
      const GenLowerTriMatrix<T4>* m4;
  };

  template <class T, class T2, class T4> class SumLU : 
    public MatrixComposite<T> 
  {
    public:
      SumLU(T _x1, const GenLowerTriMatrix<T2>& _m2, 
	  T _x3, const GenUpperTriMatrix<T4>& _m4) :
	MatrixComposite<T>(BaseStorOf(_m2)),
	x1(_x1),m2(&_m2),x3(_x3),m4(&_m4)
      {  TMVAssert(m2->size() == m4->size()); }
      inline size_t colsize() const { return m2->size(); }
      inline size_t rowsize() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline T GetX3() const { return x3; }
      inline const GenUpperTriMatrix<T4>& GetM4() const { return *m4; }
      inline void AssignTo(const MatrixView<T>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());

	AddMM(x1,*m2,x3,*m4,m0);
      }
    private:
      T x1;
      const GenLowerTriMatrix<T2>* m2;
      T x3;
      const GenUpperTriMatrix<T4>* m4;
  };

  //
  // TriMatrix * TriMatrix
  //

  template <class T, class T2, class T3> class ProdUU;
  template <class T, class T2, class T3> class ProdLL;
  template <class T, class T2, class T3> class ProdUL;
  template <class T, class T2, class T3> class ProdLU;

  template <class T, class T2, class T3> class ProdUU : 
    public UpperTriMatrixComposite<T>
  {
    public:
      ProdUU(T _x1, const GenUpperTriMatrix<T2>& _m2,
	  const GenUpperTriMatrix<T3>& _m3) :
	UpperTriMatrixComposite<T>((_m2.isunit()&&_m3.isunit()&&_x1==T(1)) ?
	    UnitDiag : NonUnitDiag,
	    BaseStorOf(_m2)), x1(_x1), m2(&_m2), m3(&_m3)
      { TMVAssert(m2->size() == m3->size()); }
      inline size_t size() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline const GenUpperTriMatrix<T3>& GetM3() const { return *m3; }
      inline void AssignTo(const UpperTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	MultMM(x1,*m2,*m3,T(0),m0);
      }
    protected:
      T x1;
      const GenUpperTriMatrix<T2>*const m2;
      const GenUpperTriMatrix<T3>*const m3;
  };

  template <class T> inline const UpperTriMatrixView<T>& operator*=(
      const UpperTriMatrixView<T>& m1, const GenUpperTriMatrix<T>& m2)
  { MultMM(T(1),m1,m2,T(0),m1); return m1; }

  template <class T> inline const UpperTriMatrixView<CT>& operator*=(
      const UpperTriMatrixView<CT>& m1, const GenUpperTriMatrix<T>& m2)
  { MultMM(CT(1),m1,m2,CT(0),m1); return m1; }

  template <class T, class T2, class T3> 
    inline const UpperTriMatrixView<T>& operator+=(
	const UpperTriMatrixView<T>& m, const ProdUU<T,T2,T3>& pmm)
    { MultMM(pmm.GetX1(),pmm.GetM2(),pmm.GetM3(),T(1),m); return m; }

  template <class T, class T2, class T3> 
    inline const UpperTriMatrixView<T>& operator-=(
	const UpperTriMatrixView<T>& m, const ProdUU<T,T2,T3>& pmm)
    { MultMM(-pmm.GetX1(),pmm.GetM2(),pmm.GetM3(),T(1),m); return m; }

  template <class T, class T2, class T3> class ProdLL : 
    public LowerTriMatrixComposite<T>
  {
    public:
      ProdLL(T _x1, const GenLowerTriMatrix<T2>& _m2,
	  const GenLowerTriMatrix<T3>& _m3) :
	LowerTriMatrixComposite<T>((_m2.isunit()&&_m3.isunit()&&_x1==T(1)) ?
	    UnitDiag : NonUnitDiag,
	    BaseStorOf(_m2)), x1(_x1), m2(&_m2), m3(&_m3)
      { TMVAssert(m2->size() == m3->size()); }
      inline size_t size() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline const GenLowerTriMatrix<T3>& GetM3() const { return *m3; }
      inline void AssignTo(const LowerTriMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	MultMM(x1,*m2,*m3,T(0),m0);
      }
    protected:
      T x1;
      const GenLowerTriMatrix<T2>*const m2;
      const GenLowerTriMatrix<T3>*const m3;
  };

  template <class T> inline const LowerTriMatrixView<T>& operator*=(
      const LowerTriMatrixView<T>& m1, const GenLowerTriMatrix<T>& m2)
  { MultMM(T(1),m1,m2,T(0),m1); return m1; }

  template <class T> inline const LowerTriMatrixView<CT>& operator*=(
      const LowerTriMatrixView<CT>& m1, const GenLowerTriMatrix<T>& m2)
  { MultMM(CT(1),m1,m2,CT(0),m1); return m1; }

  template <class T, class T2, class T3> 
    inline const LowerTriMatrixView<T>& operator+=(
	const LowerTriMatrixView<T>& m, const ProdLL<T,T2,T3>& pmm)
    { MultMM(pmm.GetX1(),pmm.GetM2(),pmm.GetM3(),T(1),m); return m; }

  template <class T, class T2, class T3> 
    inline const LowerTriMatrixView<T>& operator-=(
	const LowerTriMatrixView<T>& m, const ProdLL<T,T2,T3>& pmm)
    { MultMM(-pmm.GetX1(),pmm.GetM2(),pmm.GetM3(),T(1),m); return m; }

  template <class T, class T2, class T3> class ProdUL : 
    public MatrixComposite<T>
  {
    public:
      ProdUL(T _x1, const GenUpperTriMatrix<T2>& _m2,
	  const GenLowerTriMatrix<T3>& _m3) :
	MatrixComposite<T>(BaseStorOf(_m2)), x1(_x1), m2(&_m2), m3(&_m3)
      { TMVAssert( m2->size() == m3->size()); }
      inline size_t colsize() const { return m2->size(); }
      inline size_t rowsize() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline const GenLowerTriMatrix<T3>& GetM3() const { return *m3; }
      inline void AssignTo(const MatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	MultMM(x1,*m2,*m3,T(0),m0);
      }
    protected:
      T x1;
      const GenUpperTriMatrix<T2>*const m2;
      const GenLowerTriMatrix<T3>*const m3;
  };

  template <class T, class T2, class T3> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdUL<T,T2,T3>& pmm)
  { MultMM(pmm.GetX1(),pmm.GetM2(),pmm.GetM3(),T(1),m); return m; }

  template <class T, class T2, class T3> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdUL<T,T2,T3>& pmm)
  { MultMM(-pmm.GetX1(),pmm.GetM2(),pmm.GetM3(),T(1),m); return m; }

  template <class T, class T2, class T3> class ProdLU : 
    public MatrixComposite<T>
  {
    public:
      ProdLU(T _x1, const GenLowerTriMatrix<T2>& _m2,
	  const GenUpperTriMatrix<T3>& _m3) :
	MatrixComposite<T>(BaseStorOf(_m2)), x1(_x1), m2(&_m2), m3(&_m3)
      { TMVAssert( m2->size() == m3->size()); }
      inline size_t colsize() const { return m2->size(); }
      inline size_t rowsize() const { return m2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline const GenUpperTriMatrix<T3>& GetM3() const { return *m3; }
      inline void AssignTo(const MatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	MultMM(x1,*m2,*m3,T(0),m0);
      }
    protected:
      T x1;
      const GenLowerTriMatrix<T2>*const m2;
      const GenUpperTriMatrix<T3>*const m3;
  };

  template <class T, class T2, class T3> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdLU<T,T2,T3>& pmm)
  { MultMM(pmm.GetX1(),pmm.GetM2(),pmm.GetM3(),T(1),m); return m; }

  template <class T, class T2, class T3> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdLU<T,T2,T3>& pmm)
  { MultMM(-pmm.GetX1(),pmm.GetM2(),pmm.GetM3(),T(1),m); return m; }

  //
  // TriMatrix / % TriMatrix
  //

  template <class T0, class T1, class T2> class QuotUU : 
    public UpperTriMatrixComposite<T0>
  {
    public:
      QuotUU(const GenUpperTriMatrix<T1>& _m1,
	  const GenUpperTriMatrix<T2>& _m2) :
	UpperTriMatrixComposite<T0>(
	    _m1.dt()==_m2.dt() ? _m1.dt() : NonUnitDiag,
	    BaseStorOf(_m2)), m1(&_m1), m2(&_m2)
      { TMVAssert( m1->size() == m2->size() ); }
      inline size_t size() const { return m1->size(); }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return *m1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const UpperTriMatrixView<T0>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	m2->LDiv(*m1,m0);
      }
    protected:
      const GenUpperTriMatrix<T1>*const m1;
      const GenUpperTriMatrix<T2>*const m2;
  };

  template <class T0, class T1, class T2> class RQuotUU : 
    public UpperTriMatrixComposite<T0>
  {
    public:
      RQuotUU(const GenUpperTriMatrix<T1>& _m1,
	  const GenUpperTriMatrix<T2>& _m2) :
	UpperTriMatrixComposite<T0>(
	    _m1.dt()==_m2.dt() ? _m1.dt() : NonUnitDiag, 
	    BaseStorOf(_m2)), m1(&_m1), m2(&_m2)
      { TMVAssert( m1->size() == m2->size() ); }
      inline size_t size() const { return m1->size(); }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return *m1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const UpperTriMatrixView<T0>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	m2->RDiv(*m1,m0);
      }
    protected:
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

  template <class T0, class T1, class T2> class QuotLL : 
    public LowerTriMatrixComposite<T0>
  {
    public:
      QuotLL(const GenLowerTriMatrix<T1>& _m1,
	  const GenLowerTriMatrix<T2>& _m2) :
	LowerTriMatrixComposite<T0>(
	    _m1.dt()==_m2.dt() ? _m1.dt() : NonUnitDiag,
	    BaseStorOf(_m2)), m1(&_m1), m2(&_m2)
      { TMVAssert( m1->size() == m2->size() ); }
      inline size_t size() const { return m1->size(); }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return *m1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const LowerTriMatrixView<T0>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	m2->LDiv(*m1,m0);
      }
    protected:
      const GenLowerTriMatrix<T1>*const m1;
      const GenLowerTriMatrix<T2>*const m2;
  };

  template <class T0, class T1, class T2> class RQuotLL : 
    public LowerTriMatrixComposite<T0>
  {
    public:
      RQuotLL(const GenLowerTriMatrix<T1>& _m1,
	  const GenLowerTriMatrix<T2>& _m2) :
	LowerTriMatrixComposite<T0>(
	    _m1.dt()==_m2.dt() ? _m1.dt() : NonUnitDiag, 
	    BaseStorOf(_m2)), m1(&_m1), m2(&_m2)
      { TMVAssert( m1->size() == m2->size() ); }
      inline size_t size() const { return m1->size(); }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return *m1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const LowerTriMatrixView<T0>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	m2->RDiv(*m1,m0);
      }
    protected:
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

  //
  // Scalar / TriMatrix
  //

  template <class T0, class T1, class T2> class QuotXU : 
    public UpperTriMatrixComposite<T0> 
  {
    public:
      QuotXU(const T1 _x1, const GenUpperTriMatrix<T2>& _m2) :
	UpperTriMatrixComposite<T0>(_x1==T1(1) ? _m2.dt() : NonUnitDiag,
	    BaseStorOf(_m2)), x1(_x1), m2(&_m2)  {}
      inline size_t size() const { return m2->size(); }
      inline T1 GetX1() const { return x1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const UpperTriMatrixView<T0>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	m0 = UpperTriMatrixViewOf(m2->TInverse(),m0.dt());
	if (x1 != T1(1)) m0 *= x1;
      }
    private:
      const T1 x1;
      const GenUpperTriMatrix<T2>*const m2;
  };

  template <class T0, class T1, class T2> class QuotXL : 
    public LowerTriMatrixComposite<T0> 
  {
    public:
      QuotXL(const T1 _x1, const GenLowerTriMatrix<T2>& _m2) :
	LowerTriMatrixComposite<T0>(_x1==T1(1) ? _m2.dt() : NonUnitDiag,
	    BaseStorOf(_m2)), x1(_x1), m2(&_m2)  {}
      inline size_t size() const { return m2->size(); }
      inline T1 GetX1() const { return x1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const LowerTriMatrixView<T0>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == this->dt());
	m0 = LowerTriMatrixViewOf(m2->TInverse(),m0.dt());
	if (x1 != T1(1)) m0 *= x1;
      }
    private:
      const T1 x1;
      const GenLowerTriMatrix<T2>*const m2;
  };


  // First all the (M) op (T) possibilities

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXU

#define GENVECTOR1 GenVector
#define GENVECTOR2 GenVector
#define PRODXV1 ProdXV
#define PRODXV2 ProdXV

#define SUMMM SumMU
#define SUMMX SumUX
#define PRODXM ProdXU
#define PRODMV ProdUV
#define PRODMM ProdMU
#define QUOTVM QuotVU
#define RQUOTVM RQuotVU
#define QUOTMM QuotMU
#define RQUOTMM RQuotMU
#define TQUOTMM TransientQuotMU
#define TRQUOTMM TransientRQuotMU
#define QUOTXM QuotXU

#include "TMV_AuxVecComposite.h"
#include "TMV_AuxMatComposite1.h"

#include "TMV_AuxSumMM.h"
#include "TMV_AuxSumXM.h"
#include "TMV_AuxProdXM.h"
#include "TMV_AuxProdVM.h"
#include "TMV_AuxProdMM.h"
#include "TMV_AuxQuotVM.h"
#include "TMV_AuxQuotMM.h"
#include "TMV_AuxQuotXM.h"

#undef GENMATRIX2
#undef PRODXM2
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM2 ProdXL

#undef SUMMM
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
#define SUMMM SumML
#define SUMMX SumLX
#define PRODXM ProdXL
#define PRODMV ProdLV
#define PRODMM ProdML
#define QUOTVM QuotVL
#define RQUOTVM RQuotVL
#define QUOTMM QuotML
#define RQUOTMM RQuotML
#define TQUOTMM TransientQuotML
#define TRQUOTMM TransientRQuotML
#define QUOTXM QuotXL

#include "TMV_AuxVecComposite.h"
#include "TMV_AuxMatComposite1.h"

#include "TMV_AuxSumMM.h"
#include "TMV_AuxSumXM.h"
#include "TMV_AuxProdXM.h"
#include "TMV_AuxProdVM.h"
#include "TMV_AuxProdMM.h"
#include "TMV_AuxQuotVM.h"
#include "TMV_AuxQuotMM.h"
#include "TMV_AuxQuotXM.h"

  // Next (T) op (T)
  
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenUpperTriMatrix
#define PRODXM1 ProdXU

#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef RQUOTMM
#define SUMMM SumUL
#define PRODMM ProdUL
#define QUOTMM QuotUL
#define RQUOTMM RQuotUL

#include "TMV_AuxSumMM.h"
#include "TMV_AuxProdMM.h"
#include "TMV_AuxTQuotMM.h"

#undef GENMATRIX2
#undef PRODXM2
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM2 ProdXU

#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef RQUOTMM
#define SUMMM SumUU
#define PRODMM ProdUU
#define QUOTMM QuotUU
#define RQUOTMM RQuotUU

#include "TMV_AuxSumMM.h"
#include "TMV_AuxProdMM.h"
#include "TMV_AuxQuotMM.h"

#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL

#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM
#define SUMMM SumLU
#define PRODMM ProdLU
#define QUOTMM QuotLU
#define RQUOTMM RQuotLU
#define TQUOTMM TransientQuotMU
#define TRQUOTMM TransientRQuotMU

#include "TMV_AuxSumMM.h"
#include "TMV_AuxProdMM.h"
#include "TMV_AuxTQuotMM.h"

#undef GENMATRIX2
#undef PRODXM2
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM2 ProdXL

#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef RQUOTMM
#define SUMMM SumLL
#define PRODMM ProdLL
#define QUOTMM QuotLL
#define RQUOTMM RQuotLL

#include "TMV_AuxSumMM.h"
#include "TMV_AuxProdMM.h"
#include "TMV_AuxQuotMM.h"

  // Finally (T) op (M)
  
#undef GENMATRIX2
#undef PRODXM2
#define GENMATRIX2 GenMatrix
#define PRODXM2 ProdXM

#undef SUMMM
#undef PRODMM
#undef TQUOTMM
#undef TRQUOTMM
#define SUMMM1 SumML
#define SUMMM SumLM
#define PRODMM ProdLM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM

#include "TMV_AuxMatComposite2.h"

#include "TMV_AuxSumMM.h"
#include "TMV_AuxProdMM.h"
#include "TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenUpperTriMatrix
#define PRODXM1 ProdXU

#undef SUMMM1
#undef SUMMM
#undef PRODMM
#define SUMMM1 SumMU
#define SUMMM SumUM
#define PRODMM ProdUM

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
