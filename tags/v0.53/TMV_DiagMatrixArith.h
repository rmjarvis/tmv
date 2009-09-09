#ifndef TMV_DiagMatrixArith_H
#define TMV_DiagMatrixArith_H

#include "TMV_DiagMatrix.h"

#define CT complex<T>
#define CCT ConjRef<complex<T> >
#define VCT VarConjRef<complex<T> >

namespace tmv {

  // y = alpha * A * x + beta * y
  template <class T, class Ta,  class Tx> void MultMV(
      const T alpha, const GenDiagMatrix<Ta>& A,
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y);

  // B = alpha * A + beta * B
  template <class T, class Ta> void AddMM(
      const T alpha, const GenDiagMatrix<Ta>& A, 
      const T beta, const MatrixView<T>& B);

  // C = alpha * A * B + beta * C
  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C);
  template <class T, class Ta, class Tb> inline void MultMM(const T alpha,
	const GenMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
	const T beta, const MatrixView<T>& C)
  { MultMM(alpha,B,A.Transpose(),beta,C.Transpose()); }
  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
      const T beta, const DiagMatrixView<T>& C)
  { MultMV(alpha,A,B.diag(),beta,C.diag()); }

  template <class T> class DiagMatrixComposite : 
    public GenDiagMatrix<T>
  {
    public:
      DiagMatrixComposite() : inst(0) {}
      virtual ~DiagMatrixComposite() { if (inst) delete inst; };
      virtual size_t size() const = 0;
      virtual void AssignTo(const DiagMatrixView<T>&) const = 0;
    private:
      mutable const DiagMatrix<T>* inst;
      inline ConstVectorView<T> cdiag() const
      { if (!inst) inst = new DiagMatrix<T>(*this); return inst->diag(); }
      inline int step() const { return 1; }
  };

  template <class T, IndexStyle I, class Tx> inline DiagMatrix<T,I>& operator+=(
      DiagMatrix<T,I>& m, const Tx& x)
  { m.View() += x; return m; }

  template <class T, IndexStyle I, class Tx> inline DiagMatrix<T,I>& operator-=(
      DiagMatrix<T,I>& m, const Tx& x)
  { m.View() -= x; return m; }

  template <class T, IndexStyle I, class Tx> inline DiagMatrix<T,I>& operator*=(
      DiagMatrix<T,I>& m, const Tx& x)
  { m.View() *= x; return m; }

  template <class T, IndexStyle I, class Tx> inline DiagMatrix<T,I>& operator/=(
      DiagMatrix<T,I>& m, const Tx& x)
  { m.View() /= x; return m; }

  template <class T, IndexStyle I, class Tx> inline DiagMatrix<T,I>& operator%=(
      DiagMatrix<T,I>& m, const Tx& x)
  { m.View() %= x; return m; }

  template <class T, StorageType S, IndexStyle I, class T2> 
    inline Matrix<T,S,I>& operator+=(Matrix<T,S,I>& m1, const GenDiagMatrix<T2>& m2)
    { DiagMatrixViewOf(m1) += m2; return m1; }

  template <class T, StorageType S, IndexStyle I, class T2> 
    inline Matrix<T,S,I>& operator-=(Matrix<T,S,I>& m1, const GenDiagMatrix<T2>& m2)
    { DiagMatrixViewOf(m1) -= m2; return m1; }


  //
  // DiagMatrix + Scalar
  //

  template <class T, class T2> class SumDX : 
    public DiagMatrixComposite<T>
  {
    public:
      SumDX(T _x1, const GenDiagMatrix<T2>& _m, T _x2) :
	x1(_x1), m(&_m), x2(_x2) {}
      inline size_t size() const { return m->size(); }
      inline T GetX1() const { return x1; }
      inline const GenDiagMatrix<T2>& GetM() const { return *m; }
      inline T GetX2() const { return x2; }
      inline void AssignTo(const DiagMatrixView<T>& m0) const
      { 
	TMVAssert(m0.size() == size());
	if (x1 == T(1)) m0.diag() = m->diag();
	else m0.diag() = x1 * m->diag();
	m0.diag().AddToAll(x2);
      }
    private:
      const T x1;
      const GenDiagMatrix<T2>* m;
      const T x2;
  };

  template <class T> inline const DiagMatrixView<T>& operator+=(
      const DiagMatrixView<T>& m, T x) 
  { m.diag().AddToAll(x); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
      const DiagMatrixView<CT>& m, T x) 
  { m.diag().AddToAll(x); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
      const DiagMatrixView<CT>& m, CCT x) 
  { m.diag().AddToAll(CT(x)); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
      const DiagMatrixView<CT>& m, VCT x) 
  { m.diag().AddToAll(CT(x)); return m; }

  template <class T> inline const DiagMatrixView<T>& operator-=(
      const DiagMatrixView<T>& m, T x) 
  { m.diag().AddToAll(-x); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
      const DiagMatrixView<CT>& m, T x) 
  { m.diag().AddToAll(-x); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
      const DiagMatrixView<CT>& m, CCT x) 
  { m.diag().AddToAll(-CT(x)); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
      const DiagMatrixView<CT>& m, VCT x) 
  { m.diag().AddToAll(-CT(x)); return m; }

  //
  // Scalar * DiagMatrix
  //

  template <class T, class T2> class ProdXD : 
    public DiagMatrixComposite<T> 
  {
    public:
      ProdXD(T _x, const GenDiagMatrix<T2>& _m) : x(_x), m(&_m) {}
      inline size_t size() const { return m->size(); }
      inline T GetX() const { return x; }
      inline const GenDiagMatrix<T2>& GetM() const { return *m; }
      inline void AssignTo(const DiagMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	m0.diag() = x * m->diag();
      }
    private:
      const T x;
      const GenDiagMatrix<T2>*const m;
  };

  template <class T> inline const DiagMatrixView<T>& operator*=(
      const DiagMatrixView<T>& m, T x) 
  { m.diag() *= x; return m; }

  template <class T> inline const DiagMatrixView<CT>& operator*=(
      const DiagMatrixView<CT>& m, T x) 
  { m.diag() *= x; return m; }

  template <class T> inline const DiagMatrixView<CT>& operator*=(
      const DiagMatrixView<CT>& m, CCT x) 
  { m.diag() *= CT(x); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator*=(
      const DiagMatrixView<CT>& m, VCT x) 
  { m.diag() *= CT(x); return m; }

  template <class T> inline const DiagMatrixView<T>& operator/=(
      const DiagMatrixView<T>& m, T x) 
  { m.diag() /= x; return m; }

  template <class T> inline const DiagMatrixView<CT>& operator/=(
      const DiagMatrixView<CT>& m, T x) 
  { m.diag() /= x; return m; }

  template <class T> inline const DiagMatrixView<CT>& operator/=(
      const DiagMatrixView<CT>& m, CCT x) 
  { m.diag() /= CT(x); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator/=(
      const DiagMatrixView<CT>& m, VCT x) 
  { m.diag() /= CT(x); return m; }

  //
  // DiagMatrix + DiagMatrix
  //
  
  template <class T, class T1, class T2> class SumDD : 
    public DiagMatrixComposite<T> 
  {
    public:
      SumDD(T _x1, const GenDiagMatrix<T1>& _m1, 
	  T _x2, const GenDiagMatrix<T2>& _m2) :
	x1(_x1),m1(&_m1),x2(_x2),m2(&_m2)
	{ TMVAssert(m1->size() == m2->size()); }
      inline size_t size() const { return m1->size(); }
      inline T GetX1() const { return x1; }
      inline const GenDiagMatrix<T1>& GetM1() const { return *m1; }
      inline T GetX2() const { return x2; }
      inline const GenDiagMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const DiagMatrixView<T>& m0) const
      { 
	TMVAssert(m0.size() == size());
	m0.diag() = x1*m1->diag() + x2*m2->diag();
      }
    private:
      T x1;
      const GenDiagMatrix<T1>* m1;
      T x2;
      const GenDiagMatrix<T2>* m2;
  };

  template <class T> inline const DiagMatrixView<T>& operator+=(
      const DiagMatrixView<T>& m1, const GenDiagMatrix<T>& m2) 
  { m1.diag() += m2.diag(); return m1; }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
      const DiagMatrixView<CT>& m1, const GenDiagMatrix<T>& m2) 
  { m1.diag() += m2.diag(); return m1; }

  template <class T> inline const DiagMatrixView<T>& operator-=(
      const DiagMatrixView<T>& m1, const GenDiagMatrix<T>& m2) 
  { m1.diag() -= m2.diag(); return m1; }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
      const DiagMatrixView<CT>& m1, const GenDiagMatrix<T>& m2) 
  { m1.diag() -= m2.diag(); return m1; }

  //
  // DiagMatrix * DiagMatrix
  //

  template <class T, class T1, class T2> class ProdDD : 
    public DiagMatrixComposite<T>
  {
    public:
      ProdDD(T _x, const GenDiagMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) :
	x(_x), m1(&_m1), m2(&_m2)
      { TMVAssert( m1->size() == m2->size()); }
      inline size_t size() const { return m1->size(); }
      inline T GetX() const { return x; }
      inline const GenDiagMatrix<T1>& GetM1() const { return *m1; }
      inline const GenDiagMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const DiagMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	MultMM(x,*m1,*m2,T(0),m0);
      }
    private:
      T x;
      const GenDiagMatrix<T1>*const m1;
      const GenDiagMatrix<T2>*const m2;
  };

  template <class T> inline const DiagMatrixView<T>& operator*=(
      const DiagMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  { MultMM(T(1),m1,m2,T(0),m1); return m1; }

  template <class T> inline const DiagMatrixView<CT>& operator*=(
      const DiagMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  { MultMM(CT(1),m1,m2,CT(0),m1); return m1; }

  //
  // DiagMatrix / % DiagMatrix
  //

  template <class T, class T1, class T2> class QuotDD : 
    public DiagMatrixComposite<T>
  {
    public:
      QuotDD(const T _x, const GenDiagMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) :
	x(_x), m1(&_m1), m2(&_m2)
      { TMVAssert( m1->size() == m2->size() ); }
      inline size_t size() const { return m1->size(); }
      inline T GetX() const { return x; }
      inline const GenDiagMatrix<T1>& GetM1() const { return *m1; }
      inline const GenDiagMatrix<T2>& GetM2() const { return *m2; }
      inline void AssignTo(const DiagMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	if (m0.SameAs(*m1)) m2->DivEq(m0);
	else m2->Div(*m1,m0);
	if (x != T(1)) m0 *= x;
      }
    protected:
      const T x;
      const GenDiagMatrix<T1>*const m1;
      const GenDiagMatrix<T2>*const m2;
  };

  template <class T, class T1, class T2> class RQuotDD : 
    public QuotDD<T,T1,T2>
  {
    public:
      RQuotDD(const T _x, const GenDiagMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) : QuotDD<T,T1,T2>(_x,_m1,_m2) {}
  };

  template <class T> inline const DiagMatrixView<T>& operator/=(
      const DiagMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  { m2.DivEq(m1); return m1; }

  template <class T> inline const DiagMatrixView<CT>& operator/=(
      const DiagMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  { m2.DivEq(m1); return m1; }

  template <class T> inline const DiagMatrixView<T>& operator%=(
      const DiagMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  { m2.DivEq(m1); return m1; }

  template <class T> inline const DiagMatrixView<CT>& operator%=(
      const DiagMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  { m2.DivEq(m1); return m1; }

  //
  // Scalar / DiagMatrix
  //

  template <class T, class Tm> class QuotXD : 
    public DiagMatrixComposite<T>
  {
    public:
      QuotXD(T _x, const GenDiagMatrix<Tm>& _m) :
	x(_x), m(&_m) {}
      inline size_t size() const { return m->size(); }
      inline T GetX() const { return x; }
      inline const GenDiagMatrix<Tm>& GetM() const { return *m; }
      inline void AssignTo(const DiagMatrixView<T>& m0) const
      {
	TMVAssert(m0.size() == size());
	if (x == T(0)) m0.Zero();
	else {
	  // Need temporary, since T, Tm are not the same
	  // T = Tm overridden below.
	  DiagMatrix<Tm> temp(m0.size());
	  m->Inverse(temp.View());
	  if (x != T(1)) m0 = temp*x;
	  else m0 = temp;
	}
      }
    protected:
      T x;
      const GenDiagMatrix<Tm>*const m;
  };

  template <class T> class QuotXD<T,T> : 
    public DiagMatrixComposite<T>
    {
      public:
	QuotXD(T _x, const GenDiagMatrix<T>& _m) :
	  x(_x), m(&_m) {}
	inline size_t size() const { return m->size(); }
	inline T GetX() const { return x; }
	inline const GenDiagMatrix<T>& GetM() const { return *m; }
	inline void AssignTo(const DiagMatrixView<T>& m0) const
	{
	  TMVAssert(m0.size() == size());
	  if (x == T(0)) m0.Zero();
	  else {
	    m->Inverse(m0);
	    if (x != T(1)) m0 *= x;
	  }
	}
      protected:
	T x;
	const GenDiagMatrix<T>*const m;
    };

  // First all the (M) op (D) possibilities

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXD

#define GENVECTOR1 GenVector
#define GENVECTOR2 GenVector
#define PRODXV1 ProdXV
#define PRODXV2 ProdXV

#define SUMMM SumMD
#define SUMMX SumDX
#define PRODXM ProdXD
#define PRODMV ProdDV
#define PRODMM ProdMD
#define QUOTVM QuotVD
#define RQUOTVM RQuotVD
#define QUOTMM QuotMD
#define RQUOTMM RQuotMD
#define TQUOTMM TransientQuotMD
#define TRQUOTMM TransientRQuotMD
#define QUOTXM QuotXD

#include "TMV_AuxVecComposite.h"

  // transm is unimportant for the calculation of MV for DiagMatrices,
  // so specialize the trans = true version to be the same as the
  // trans = false version.
  template <class T, class T1, class T2> class ProdDV<T,T1,true,T2> :
    public ProdDV<T,T1,false,T2>
    {
      public:
	ProdDV(T _x, const GenDiagMatrix<T1>& _m, const GenVector<T2>& _v) : 
	  ProdDV<T,T1,false,T2>(_x,_m,_v) {}
    };

#include "TMV_AuxMatComposite1.h"

#include "TMV_AuxSumMM.h"
#include "TMV_AuxSumXM.h"
#include "TMV_AuxProdXM.h"
#include "TMV_AuxProdVM.h"
#include "TMV_AuxProdMM.h"
#include "TMV_AuxQuotVM.h"
#include "TMV_AuxQuotMM.h"
#include "TMV_AuxQuotXM.h"

  // Next (D) op (D)
  
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenDiagMatrix
#define PRODXM1 ProdXD

#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef RQUOTMM
#define SUMMM SumDD
#define PRODMM ProdDD
#define QUOTMM QuotDD
#define RQUOTMM RQuotDD

#include "TMV_AuxSumMM.h"
#include "TMV_AuxProdMM.h"
#include "TMV_AuxQuotMM.h"

  // Next  (M) op (D)
  
#undef GENMATRIX2
#undef PRODXM2
#define GENMATRIX2 GenMatrix
#define PRODXM2 ProdXM

#undef SUMMM
#undef PRODMM
#undef TQUOTMM
#undef TRQUOTMM
#define SUMMM1 SumMD
#define SUMMM SumDM
#define PRODMM ProdDM
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
