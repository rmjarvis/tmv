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

  // A = beta*A + alpha * x * yT
  // beta must be 0 or 1
  template <class T, class Tx, class Ty> void Rank1Update(const T alpha,
      const GenVector<Tx>& x, const GenVector<Ty>& y, 
      const int beta, const MatrixView<T>& A);
       
  template <class T, class Ta> void ElementProd(
      const T alpha, const GenMatrix<Ta>& A, const MatrixView<T>& B);
  template <class T, class Ta, class Tb> void AddElementProd(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C);

  template <class T> class MatrixComposite : public GenMatrix<T>
  {
    public:

      MatrixComposite(StorageType s) : GenMatrix<T>(s,NonConj,1), itsm(0) {}
      MatrixComposite(const MatrixComposite<T>& rhs) : 
	GenMatrix<T>(rhs), itsm(0) {}
      virtual ~MatrixComposite() {}
      virtual void AssignTo(const MatrixView<T>&) const = 0;
      inline const T* cptr() const
      { 
	if (!itsm.get()) {
	  size_t len = colsize()*rowsize();
	  itsm.reset(new T[len]);
	  MatrixView<T>(itsm.get(),colsize(),rowsize(),stepi(),stepj(),stor(),
	      NonConj,len FIRSTLAST1(itsm.get(),itsm.get()+len) ) = *this;
	}
	return itsm.get();
      }
      using GenMatrix<T>::rowsize;
      using GenMatrix<T>::colsize;
      using GenMatrix<T>::isrm;
      using GenMatrix<T>::stor;
      inline int stepi() const { return isrm() ? rowsize() : 1; }
      inline int stepj() const { return isrm() ? 1 : colsize(); }

    private:

      mutable auto_array<T> itsm;
  };

  template <class T, StorageType S, IndexStyle I, class Tx> 
    inline Matrix<T,S>& operator+=(Matrix<T,S,I>& m, const Tx& x) 
  { m.View() += x; return m; }

  template <class T, StorageType S, IndexStyle I, class Tx> 
    inline Matrix<T,S>& operator-=(Matrix<T,S,I>& m, const Tx& x) 
  { m.View() -= x; return m; }

  template <class T, StorageType S, IndexStyle I, class Tx> 
    inline Matrix<T,S>& operator*=(Matrix<T,S,I>& m, const Tx& x) 
  { m.View() *= x; return m; }

  template <class T, StorageType S, IndexStyle I, class Tx> 
    inline Matrix<T,S>& operator/=(Matrix<T,S,I>& m, const Tx& x) 
  { m.View() /= x; return m; }

  template <class T, StorageType S, IndexStyle I, class Tx> 
    inline Matrix<T,S>& operator%=(Matrix<T,S,I>& m, const Tx& x) 
  { m.View() %= x; return m; }

  //
  // Scalar * Matrix
  //

  template <class T, class Tm> class ProdXM : 
    public MatrixComposite<T> 
  {
    // x*m
    public:
      ProdXM(const T _x, const GenMatrix<Tm>& _m) :
	MatrixComposite<T>(BaseStorOf(_m)),
	x(_x), m(&_m) {}
      inline size_t colsize() const { return m->colsize(); }
      inline size_t rowsize() const { return m->rowsize(); }
      inline T GetX() const { return x; }
      inline const GenMatrix<Tm>& GetM() const { return *m; }
      inline void AssignTo(const MatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	m0 = *m;
	if (x != T(1)) MultXM(x,m0);
      }
    private:
      const T x;
      const GenMatrix<Tm>*const m;
  };

  // m*=x
  template <class T> inline const MatrixView<T>& operator*=(
      const MatrixView<T>& m, T x) 
  { MultXM(x,m); return m; }

  template <class T> inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m, T x) 
  { MultXM(CT(x),m); return m; }

  template <class T> inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m, CCT x) 
  { MultXM(CT(x),m); return m; }

  template <class T> inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m, VCT x) 
  { MultXM(CT(x),m); return m; }

  // m/=x
  template <class T> inline const MatrixView<T>& operator/=(
      const MatrixView<T>& m, T x) 
  { MultXM(T(1)/x,m); return m; }

  template <class T> inline const MatrixView<CT>& operator/=(
      const MatrixView<CT>& m, T x) 
  { MultXM(CT(T(1)/x),m); return m; }

  template <class T> inline const MatrixView<CT>& operator/=(
      const MatrixView<CT>& m, CCT x) 
  { MultXM(T(1)/CT(x),m); return m; }

  template <class T> inline const MatrixView<CT>& operator/=(
      const MatrixView<CT>& m, VCT x) 
  { MultXM(T(1)/CT(x),m); return m; }

#define GENMATRIX GenMatrix
#define PRODXM ProdXM
#include "TMV_AuxProdXM.h"
  // Defines things like -m, x*m, m*x, x*(x*m), etc.

  //
  // Matrix + Scalar
  //

  template <class T, class Tm> class SumMX : 
    public MatrixComposite<T>
  {
    // x1*m + x2
    public:
      SumMX(T _x1, const GenMatrix<Tm>& _m, T _x2) :
	MatrixComposite<T>(BaseStorOf(_m)),
	x1(_x1), m(&_m), x2(_x2)
	{ TMVAssert(m->IsSquare()); }
      inline size_t colsize() const { return m->colsize(); }
      inline size_t rowsize() const { return m->rowsize(); }
      inline T GetX1() const { return x1; }
      inline const GenMatrix<Tm>& GetM() const { return *m; }
      inline T GetX2() const { return x2; }
      inline void AssignTo(const MatrixView<T>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	if (x1 == T(1)) m0 = *m;
	else m0 = x1*(*m);
	if (x2 != T(0)) m0.diag().AddToAll(x2);
      }
    private:
      const T x1;
      const GenMatrix<Tm>* m;
      const T x2;
  };

  // m+=x
  template <class T> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, T x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(x); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, T x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(x); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, CCT x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(x)); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, VCT x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(x)); 
    return m;
  }

  // m-=x
  template <class T> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, T x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(-x); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, T x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(-x); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, CCT x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(-x)); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, VCT x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(-x)); 
    return m; 
  }

#define SUMMX SumMX
#include "TMV_AuxSumMX.h"
  // Defines things like m+x, x+m, x-m, m-x, x+x*m, x*(x+m), etc.
#undef SUMMX

  //
  // Vector ^ Vector (OuterProduct)
  //

  template <class T, class T1, class T2> class OProdVV : 
    public MatrixComposite<T>
  {
    public:
      OProdVV(const T _x, const GenVector<T1>& _v1,
	  const GenVector<T2>& _v2) :
	MatrixComposite<T>(RowMajor), // arbitrary
	x(_x), v1(&_v1), v2(&_v2) {}
      inline size_t colsize() const { return v1->size(); }
      inline size_t rowsize() const { return v2->size(); }
      inline T GetX() const { return x; }
      inline const GenVector<T1>& GetV1() const { return *v1; }
      inline const GenVector<T2>& GetV2() const { return *v2; }
      inline void AssignTo(const MatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	Rank1Update(x, *v1, *v2, 0, m0);
      }
    private:
      T x;
      const GenVector<T1>*const v1;
      const GenVector<T2>*const v2;
  };

  // m+=(x*v^v)
  template <class T, class T1, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m0, const OProdVV<T,T1,T2>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update(opvv.GetX(), opvv.GetV1(), opvv.GetV2(), 1, m0);
    return m0;
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m0, const OProdVV<T,T,T>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update(CT(opvv.GetX()), opvv.GetV1(), opvv.GetV2(), 1, m0);
    return m0;
  }

  // m-=(x*v^v)
  template <class T, class T1, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m0, const OProdVV<T,T1,T2>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update(-opvv.GetX(), opvv.GetV1(), opvv.GetV2(), 1, m0);
    return m0;
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m0, const OProdVV<T,T,T>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update(CT(-opvv.GetX()), opvv.GetV1(), opvv.GetV2(), 1, m0);
    return m0;
  }

#define OPRODVV OProdVV
#define GENVECTOR1 GenVector
#define GENVECTOR2 GenVector
#define PRODXV1 ProdXV
#define PRODXV2 ProdXV
#include "TMV_AuxOProdVV.h"
#undef OPRODVV
#undef GENVECTOR1
#undef GENVECTOR2
#undef PRODXV1
#undef PRODXV2


  // Now use some helper files for combinations which are basically the same 
  // for some of the special matrices.
  
  // First things which return a MatrixComposite:

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXM
#define SUMMM SumMM
#define PRODMM ProdMM
#define QUOTXM QuotXM
#define QUOTMM QuotMM
#define RQUOTMM RQuotMM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM

#include "TMV_AuxMatComposite1.h"
  // Defines m+m, m*m, m/m, m%m
#include "TMV_AuxMatComposite3.h"
  // Defines x/m

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM

  // Next things which return a VectorComposite:
  // (Needs to be second so QuotXM is defined - in TMV_AuxMatComposite3.h)
  
#define PRODMV ProdMV
#define QUOTVM QuotVM
#define RQUOTVM RQuotVM
#include "TMV_AuxVecComposite.h"
#undef PRODMV
#undef QUOTVM
#undef RQUOTVM
#undef PRODXM
#undef QUOTXM
#undef GENMATRIX


} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif
