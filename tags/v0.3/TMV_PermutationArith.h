//---------------------------------------------------------------------------
#ifndef TMV_PermutationArith_H
#define TMV_PermutationArith_H

#include "TMV_Permutation.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"

namespace tmv {

  void MultPP(const Permutation& p1, const bool inv1, const Permutation& p2,
      const bool inv2, Permutation& p0);
  template <class T> void MultPV(const Permutation& p1, const bool inv1,
      const GenVector<T>& v2, const VectorView<T>& v0);
  template <class T> void MultPM(const Permutation& p1, const bool inv1,
      const GenMatrix<T>& m2, const MatrixView<T>& m0);

  //
  // P^n
  //

  class PowerPN : public PermutationComposite
  {
    public:
      PowerPN(const Permutation& _p, int _n) : p(&_p), n(_n) {}
      size_t size() const { return p->size(); }
      void AssignTo(Permutation& p0) const
      { 
	TMVAssert(p0.size() == p->size());
	p0 = *p;
	p0 ^= n; 
      }
    private:
      const Permutation*const p;
      int n;
  };

  inline PowerPN operator^(const Permutation& p, int n)
  { return PowerPN(p,n); }

  //
  // Permutation of Permutation
  //

  class ProdPP : public PermutationComposite
  {
    public :

      ProdPP(const Permutation& _p1, const bool _inv1,
	  const Permutation& _p2, const bool _inv2) :
	p1(&_p1), inv1(_inv1), p2(&_p2), inv2(_inv2)
	{ TMVAssert(_p1.size() == _p2.size()); }
      size_t size() const { return p1->size(); }

      void AssignTo(Permutation& p0) const
      {
	TMVAssert(p0.size() == p1->size());
	MultPP(*p1,inv1,*p2,inv2,p0);
      }
    private:
      const Permutation*const p1;
      const bool inv1;
      const Permutation*const p2;
      const bool inv2;
  };

  inline ProdPP operator*(
      const Permutation& p1, const Permutation& p2)
  { return ProdPP(p1,false,p2,false); }

  inline ProdPP operator/(
      const Permutation& p1, const Permutation& p2)
  { return ProdPP(p2,true,p1,false); }

  inline ProdPP operator%(
      const Permutation& p1, const Permutation& p2)
  { return ProdPP(p1,false,p2,true); }


  //
  // Permutations of Vectors
  //

  template <class T> class ProdPV :
    public VectorComposite<T> 
  {
    public :
      ProdPV(const Permutation& _p1, const bool _inv1,
	  const GenVector<T>& _v2) :
	p1(&_p1), inv1(_inv1), v2(&_v2)
	{ TMVAssert(p1->size() == v2->size()); }
      size_t size() const { return v2->size(); }
      void AssignTo(const VectorView<T>& v0) const
      {
	TMVAssert(v0.size() == v2->size()); 
	MultPV(*p1,inv1,*v2,v0);
      }
    private :
      const Permutation*const p1;
      const bool inv1;
      const GenVector<T>*const v2;
  };

  template <class T> inline ProdPV<T> operator*(
      const Permutation& p, const GenVector<T>& v)
    // p * v
  { return ProdPV<T>(p,false,v); }

  template <class T> inline ProdPV<T> operator*(
      const GenVector<T>& v, const Permutation& p)
    // v * p => pt * v
  { return ProdPV<T>(p,true,v); }

  template <class T> inline const VectorView<T>& operator*=(
      const VectorView<T>& v, const Permutation& p)
  { if (p.IsInverse()) p.LMultEq(v); else p.LDivEq(v); return v; }

  template <class T> inline ProdPV<T> operator/(
      const GenVector<T>& v, const Permutation& p)
    // v / p => pt * v
  { return ProdPV<T>(p,true,v); }

  template <class T> inline const VectorView<T>& operator/=(
      const VectorView<T>& v, const Permutation& p)
  { if (p.IsInverse()) p.LMultEq(v); else p.LDivEq(v); return v; }

  template <class T> inline ProdPV<T> operator%(
      const GenVector<T>& v, const Permutation& p)
    // v % p => v * pt => p * v
  { return ProdPV<T>(p,false,v); }

  template <class T> inline const VectorView<T>& operator%=(
      const VectorView<T>& v, const Permutation& p)
  { if (!p.IsInverse()) p.LMultEq(v); else p.LDivEq(v); return v; }


  //
  // Permutations of Matrices
  //

  template <class T> class ProdPM : public MatrixComposite<T>
  {
    public :
      ProdPM(const Permutation& _p1, const bool _inv1, 
	  const GenMatrix<T>& _m2) :
	MatrixComposite<T>(BaseStorOf(_m2)),
	p1(&_p1), inv1(_inv1), m2(&_m2)
	{ TMVAssert(p1->size() == m2->colsize()); }
      size_t colsize() const { return m2->colsize(); }
      size_t rowsize() const { return m2->rowsize(); }
      void AssignTo(const MatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == m2->colsize()); 
	TMVAssert(m0.rowsize() == m2->rowsize()); 
	MultPM(*p1,inv1,*m2,m0);
      }
    private :
      const Permutation*const p1;
      const bool inv1;
      const GenMatrix<T>*const m2;
  };

  template <class T> class ProdMP : public MatrixComposite<T>
  {
    public :
      ProdMP(const GenMatrix<T>& _m1, const Permutation& _p2, 
	  const bool _inv2) :
	MatrixComposite<T>(BaseStorOf(_m1)),
	m1(&_m1), p2(&_p2), inv2(_inv2)
	{ TMVAssert(m1->rowsize() == p2->size()); }
      size_t colsize() const { return m1->colsize(); }
      size_t rowsize() const { return m1->rowsize(); }
      void AssignTo(const MatrixView<T>& m0) const
      {
	TMVAssert(m0.colsize() == m1->colsize()); 
	TMVAssert(m0.rowsize() == m1->rowsize()); 
	// m0 = m1 * p
	// => m0t = pt * m1t
	MultPM(*p2,!inv2,m1->QuickTranspose(),m0.QuickTranspose());
      }
    private :
      const GenMatrix<T>*const m1;
      const Permutation*const p2;
      const bool inv2;
  };

  template <class T> inline ProdPM<T> operator*(
      const Permutation& p, const GenMatrix<T>& m)
    // p * m
  { return ProdPM<T>(p,false,m); }

  template <class T> inline ProdMP<T> operator*(
      const GenMatrix<T>& m, const Permutation& p)
    // m * p
  { return ProdMP<T>(m,p,false); }

  template <class T> inline ProdPM<T> operator/(
      const GenMatrix<T>& m, const Permutation& p)
    // pt * m
  { return ProdPM<T>(p,true,m); }

  template <class T> inline ProdMP<T> operator%(
      const GenMatrix<T>& m, const Permutation& p)
    // m * pt
  { return ProdMP<T>(m,p,true); }

  template <class T> inline const MatrixView<T>& operator*=(
      const MatrixView<T>& m, const Permutation& p)
  { 
    if (p.IsInverse()) p.LMultEq(m.QuickTranspose());
    else p.LDivEq(m.QuickTranspose()); 
    return m; 
  }

  template <class T> inline const MatrixView<T>& operator/=(
      const MatrixView<T>& m, const Permutation& p)
  { if (p.IsInverse()) p.LMultEq(m); else p.LDivEq(m); return m; }

  template <class T> inline const MatrixView<T>& operator%=(
      const MatrixView<T>& m, const Permutation& p)
  { 
    if (!p.IsInverse()) p.LMultEq(m.QuickTranspose());
    else p.LDivEq(m.QuickTranspose()); 
    return m; 
  }


} // namespace tmv

#endif
