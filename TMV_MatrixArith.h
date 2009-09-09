//---------------------------------------------------------------------------
#ifndef TMV_MatrixArith_H
#define TMV_MatrixArith_H

#include "TMV_Matrix.h"
#include "TMV_VectorArith.h"

#define CT complex<T>
#define CCT ConjRef<complex<T> >
#define VCT VarConjRef<complex<T> >

namespace tmv {

  // y = alpha * A * x + beta * y
  template <class T, class Ta, class Tx> void MultMV(const T alpha,
      const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y);

  // A = alpha * A
  template <class T> void MultXM(const T alpha, const MatrixView<T>& A);

  // B = alpha * A + beta * B
  template <class T, class Ta> void AddMM(const T alpha,
      const GenMatrix<Ta>& A, const T beta, const MatrixView<T>& B);

  // C = alpha * A * B + beta * C
  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C);

  // A = A + alpha * x * yT
  template <class T, class Tx, class Ty> void Rank1Update(const T alpha,
      const GenVector<Tx>& x, const GenVector<Ty>& y, const MatrixView<T>& A);
       
  template <class T> class MatrixComposite : public GenMatrix<T>
  {
    public:

      MatrixComposite(StorageType s) : GenMatrix<T>(s,NonConj), itsm(0) {}
      virtual ~MatrixComposite() { if (itsm) delete [] itsm; };
      virtual void AssignTo(const MatrixView<T>&) const = 0;
      inline const T* cptr() const
      { 
	if (!itsm) {
	  itsm = new T[colsize()*rowsize()];
	  MatrixViewOf(itsm,colsize(),rowsize(),stor()) = *this;
	}
	return itsm;
      }
      using GenMatrix<T>::rowsize;
      using GenMatrix<T>::colsize;
      using GenMatrix<T>::isrm;
      using GenMatrix<T>::stor;
      inline int stepi() const { return isrm() ? rowsize() : 1; }
      inline int stepj() const { return isrm() ? 1 : colsize(); }
    private:
      mutable T* itsm;
  };

  template <class T, StorageType S, class T2> inline Matrix<T,S>& operator+=(
      Matrix<T,S>& m1, const T2& x2) 
  { m1.View() += x2; return m1; }

  template <class T, StorageType S, class T2> inline Matrix<T,S>& operator-=(
      Matrix<T,S>& m1, const T2& x2) 
  { m1.View() -= x2; return m1; }

  template <class T, StorageType S, class T2> inline Matrix<T,S>& operator*=(
      Matrix<T,S>& m1, const T2& x2) 
  { m1.View() *= x2; return m1; }

  template <class T, StorageType S, class T2> inline Matrix<T,S>& operator/=(
      Matrix<T,S>& m1, const T2& x2) 
  { m1.View() /= x2; return m1; }

  template <class T, StorageType S, class T2> inline Matrix<T,S>& operator%=(
      Matrix<T,S>& m1, const T2& x2) 
  { m1.View() %= x2; return m1; }

  //
  // Matrix + Scalar
  //

  template <class T, class T2> class SumMX : 
    public MatrixComposite<T>
  {
    public:
      SumMX(T _x1, const GenMatrix<T2>& _m2, T _x3) :
	MatrixComposite<T>(BaseStorOf(_m2)),
	x1(_x1), m2(&_m2), x3(_x3)
	{ TMVAssert(m2->IsSquare()); }
      inline size_t colsize() const { return m2->colsize(); }
      inline size_t rowsize() const { return m2->rowsize(); }
      inline T GetX1() const { return *x3; }
      inline const GenMatrix<T2>& GetM2() const { return *m2; }
      inline T GetX3() const { return *x3; }
      inline void AssignTo(const MatrixView<T>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	if (x1 == T(1)) m0 = *m2;
	else m0 = x1*(*m2);
	m0.diag().AddToAll(x3);
      }
    private:
      const T x1;
      const GenMatrix<T2>* m2;
      const T x3;
  };

  // m+=x
  template <class T> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m1, T x2) 
  {
    TMVAssert(m1.IsSquare());
    m1.diag().AddToAll(x2); 
    return m1; 
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m1, T x2) 
  {
    TMVAssert(m1.IsSquare());
    m1.diag().AddToAll(x2); 
    return m1; 
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m1, CCT x2) 
  {
    TMVAssert(m1.IsSquare());
    m1.diag().AddToAll(CT(x2)); 
    return m1; 
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m1, VCT x2) 
  {
    TMVAssert(m1.IsSquare());
    m1.diag().AddToAll(CT(x2)); 
    return m1; 
  }

  // m-=x
  template <class T> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m1, T x2) 
  {
    TMVAssert(m1.IsSquare());
    m1.diag().AddToAll(-x2); 
    return m1; 
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m1, T x2) 
  {
    TMVAssert(m1.IsSquare());
    m1.diag().AddToAll(-x2); 
    return m1; 
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m1, CCT x2) 
  {
    TMVAssert(m1.IsSquare());
    m1.diag().AddToAll(CT(-x2)); 
    return m1; 
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m1, VCT x2) 
  {
    TMVAssert(m1.IsSquare());
    m1.diag().AddToAll(CT(-x2)); 
    return m1; 
  }


  //
  // Scalar * Matrix
  //

  template <class T, class T2> class ProdXM : 
    public MatrixComposite<T> 
  {
    public:
      ProdXM(const T _x1, const GenMatrix<T2>& _m2) :
	MatrixComposite<T>(BaseStorOf(_m2)),
	x1(_x1), m2(&_m2) {}
      inline size_t colsize() const { return m2->colsize(); }
      inline size_t rowsize() const { return m2->rowsize(); }
      inline T GetX1() const { return x1; }
      inline const GenMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const MatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	m0 = *m2;
	MultXM(x1,m0);
      }
    private:
      const T x1;
      const GenMatrix<T2>*const m2;
  };

  // m*=x
  template <class T> inline const MatrixView<T>& operator*=(
      const MatrixView<T>& m1, T x2) 
  { MultXM(x2,m1); return m1; }

  template <class T> inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m1, T x2) 
  { MultXM(CT(x2),m1); return m1; }

  template <class T> inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m1, CCT x2) 
  { MultXM(CT(x2),m1); return m1; }

  template <class T> inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m1, VCT x2) 
  { MultXM(CT(x2),m1); return m1; }

  // m/=x
  template <class T> inline const MatrixView<T>& operator/=(
      const MatrixView<T>& m1, T x2) 
  { MultXM(T(1)/x2,m1); return m1; }

  template <class T> inline const MatrixView<CT>& operator/=(
      const MatrixView<CT>& m1, T x2) 
  { MultXM(CT(T(1)/x2),m1); return m1; }

  template <class T> inline const MatrixView<CT>& operator/=(
      const MatrixView<CT>& m1, CCT x2) 
  { MultXM(T(1)/CT(x2),m1); return m1; }

  template <class T> inline const MatrixView<CT>& operator/=(
      const MatrixView<CT>& m1, VCT x2) 
  { MultXM(T(1)/CT(x2),m1); return m1; }

  //
  // Vector ^ Vector (OuterProduct)
  //

  template <class T, class T2, class T3> class OProdVV : 
    public MatrixComposite<T>
  {
    public:
      OProdVV(const T _x1, const GenVector<T2>& _v2,
	  const GenVector<T3>& _v3) :
	MatrixComposite<T>(RowMajor), // arbitrary
	x1(_x1), v2(&_v2), v3(&_v3) {}
      inline size_t colsize() const { return v2->size(); }
      inline size_t rowsize() const { return v3->size(); }
      inline T GetX1() const { return x1; }
      inline const GenVector<T2>& GetV2() const { return *v2; }
      inline const GenVector<T3>& GetV3() const { return *v3; }
      inline void AssignTo(const MatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	m0.Zero();
	Rank1Update(x1, *v2, *v3, m0);
      }
    private:
      T x1;
      const GenVector<T2>*const v2;
      const GenVector<T3>*const v3;
  };

  // m+=(x*v^v)
  template <class T, class T2, class T3> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m0, const OProdVV<T,T2,T3>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update(opvv.GetX1(), opvv.GetV2(), opvv.GetV3(), m0);
    return m0;
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m0, const OProdVV<T,T,T>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update(CT(opvv.GetX1()), opvv.GetV2(), opvv.GetV3(), m0);
    return m0;
  }

  // m-=(x*v^v)
  template <class T, class T2, class T3> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m0, const OProdVV<T,T2,T3>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update(-opvv.GetX1(), opvv.GetV2(), opvv.GetV3(), m0);
    return m0;
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m0, const OProdVV<T,T,T>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update(CT(-opvv.GetX1()), opvv.GetV2(), opvv.GetV3(), m0);
    return m0;
  }

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenMatrix

#define SUMMM SumMM
#define SUMMX SumMX
#define PRODXM ProdXM
#define PRODXM1 ProdXM
#define PRODXM2 ProdXM
#define PRODMV ProdMV
#define PRODMM ProdMM
#define QUOTVM QuotVM
#define QUOTMM QuotMM
#define QUOTXM QuotXM
#define RQUOTVM RQuotVM
#define RQUOTMM RQuotMM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM
#define OPRODVV OProdVV

#define GENVECTOR1 GenVector
#define GENVECTOR2 GenVector
#define PRODXV1 ProdXV
#define PRODXV2 ProdXV

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
#include "TMV_AuxOProdVV.h"

#undef GENVECTOR1
#undef GENVECTOR2
#undef PRODXV1
#undef PRODXV2

#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef SUMMX
#undef PRODXM
#undef PRODXM1
#undef PRODXM2
#undef PRODMV
#undef PRODMM
#undef QUOTVM
#undef QUOTMM
#undef QUOTXM
#undef RQUOTVM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM
#undef OPRODVV

}; // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif
