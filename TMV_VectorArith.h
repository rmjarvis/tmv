//-----------------------------------------------------------------------------
//
// This file contains the code to do Vector arithmetic (including mixing
// real and complex types) :
//
// Operators:
//
//    -v
//
//    v = v
//
//    v += v
//    v + v
//
//    v -= v
//    v - v
//
//    v *= x
//    x * v
//    v * x
//    v * v   (inner product)
//
//    v /= x
//    v / x
//
//


#ifndef TMV_VectorArith_H
#define TMV_VectorArith_H

#include "TMV_Vector.h"

#define CT complex<T>
#define CCT ConjRef<complex<T> >
#define VCT VarConjRef<complex<T> >

namespace tmv {

  // v = x * v
  template <class T> void MultXV(const T x, const VectorView<T>& v2);
  // v2 = x * v1
  template <class T, class T1> void MultXV(const T x,
      const GenVector<T1>& v1, const VectorView<T>& v2);

  // v2 = x * v1 + v2
  template <class T, class T1> void AddVV(const T x, 
      const GenVector<T1>& v1, const VectorView<T>& v2);
  // v2 = x1*v1 + x2*v2
  template <class T, class T1> void AddVV(const T x1, const GenVector<T1>& v1,
      const T x2, const VectorView<T>& v2);

  // v1 * v2 (dot product)
  // Note: the return type is the type of the first vector
  // This is important for mixing complex and real vectors
  template <class T, class T2> T MultVV(
      const GenVector<T>& v1, const GenVector<T2>& v2);
  template <class T> inline complex<T> MultVV(
      const GenVector<T>& v1, const GenVector<complex<T> >& v2)
  { return MultVV(v2,v1); }

  template <class T, class Tx> void ElementProd(
      const T alpha, const GenVector<Tx>& x, const VectorView<T>& y);
  template <class T, class Tx, class Ty> void AddElementProd(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const VectorView<T>& z);

  template <class T> class VectorComposite : 
    public GenVector<T>
  {
    public:

      VectorComposite() : GenVector<T>(NonConj), itsv(0) {}
      VectorComposite(const VectorComposite<T>& rhs) : 
	GenVector<T>(rhs), itsv(0) {}
      inline virtual ~VectorComposite() {}
      virtual void AssignTo(const VectorView<T>& v0) const = 0;

      inline const T* cptr() const
      {
	if (!itsv.get()) {
	  size_t len = size();
	  itsv.reset(new T[len]);
	  VectorView<T>(itsv.get(),len,step(),NonConj 
	      FIRSTLAST1(itsv.get(),itsv.get()+len) ) = *this;
	}
	return itsv.get();
      }

      using GenVector<T>::size;
      inline int step() const { return 1; }

    private:

      void operator=(const VectorComposite<T>& rhs) const 
      { TMVAssert(FALSE); }

      mutable auto_array<T> itsv;
  };

  // These are what we want to do no matter what type Tx is:
  template <class T, IndexStyle I, class Tx> inline Vector<T>& operator+=(
      Vector<T,I>& v, const Tx& x)
  { v.View() += x; return v; }

  template <class T, IndexStyle I, class Tx> inline Vector<T>& operator-=(
      Vector<T,I>& v, const Tx& x)
  { v.View() -= x; return v; }

  template <class T, IndexStyle I, class Tx> inline Vector<T>& operator*=(
      Vector<T,I>& v, const Tx& x)
  { v.View() *= x; return v; }

  template <class T, IndexStyle I, class Tx> inline Vector<T>& operator/=(
      Vector<T,I>& v, const Tx& x)
  { v.View() /= x; return v; }

  template <class T, IndexStyle I, class Tx> inline Vector<T>& operator%=(
      Vector<T,I>& v, const Tx& x)
  { v.View() %= x; return v; }


  //
  // Vector * / Scalar
  //

  template <class T, class Tv> class ProdXV : public VectorComposite<T> 
  {
    public:
      ProdXV(T _x, const GenVector<Tv>& _v) : x(_x), v(_v) {}
      virtual ~ProdXV() {}
      virtual inline size_t size() const { return v.size(); }
      inline T GetX() const { return x; }
      inline const GenVector<Tv>& GetV() const { return v; }
      inline void AssignTo(const VectorView<T>& v0) const
      {
	TMVAssert(v0.size() == size());
	MultXV(x,v,v0);
      }
    private:
      const T x;
      const GenVector<Tv>& v;
  };

  template <class T> inline const VectorView<T>& operator*=(
      const VectorView<T>& v1, T x2)
  { MultXV(x2,v1); return v1; }

  template <class T> inline const VectorView<CT>& operator*=(
      const VectorView<CT>& v1, T x2)
  { MultXV(CT(x2),v1); return v1; }

  template <class T> inline const VectorView<CT>& operator*=(
      const VectorView<CT>& v1, CCT x2)
  { MultXV(CT(x2),v1); return v1; }

  template <class T> inline const VectorView<CT>& operator*=(
      const VectorView<CT>& v1, VCT x2)
  { MultXV(CT(x2),v1); return v1; }

  template <class T> inline const VectorView<T>& operator/=(
      const VectorView<T>& v1, T x2)
  { MultXV(T(1)/x2,v1); return v1; }

  template <class T> inline const VectorView<CT>& operator/=(
      const VectorView<CT>& v1, T x2)
  { MultXV(CT(T(1)/x2),v1); return v1; }

  template <class T> inline const VectorView<CT>& operator/=(
      const VectorView<CT>& v1, CCT x2)
  { MultXV(T(1)/CT(x2),v1); return v1; }

  template <class T> inline const VectorView<CT>& operator/=(
      const VectorView<CT>& v1, VCT x2)
  { MultXV(T(1)/CT(x2),v1); return v1; }

  //
  // Vector + Vector
  //

  template <class T, class T1, class T2> class SumVV : 
    public VectorComposite<T> 
  {
    public:
      SumVV(T _x1, const GenVector<T1>& _v1, T _x2, const GenVector<T2>& _v2) :
	x1(_x1),v1(_v1),x2(_x2), v2(_v2)
      { TMVAssert(v1.size() == v2.size()); }
      inline ~SumVV() {}
      inline size_t size() const { return v1.size(); }
      inline T GetX1() const { return x1; }
      inline const GenVector<T1>& GetV1() const { return v1; }
      inline T GetX2() const { return x2; }
      inline const GenVector<T2>& GetV2() const { return v2; }
      inline void AssignTo(const VectorView<T>& v0) const
      {
	if (v0.SameStorageAs(v1)) {
	  v0 = v1; 
	  if (x1 == T(1))
	    AddVV(x2,v2,v0);
	  else
	    AddVV(x2,v2,x1,v0);
	} else if (v0.SameStorageAs(v2) || x1 != T(1)) {
	  v0 = v2;
	  if (x2 == T(1))
	    AddVV(x1,v1,v0);
	  else
	    AddVV(x1,v1,x2,v0);
	} else {
	  // x1 == T(1) and no same storage
	  v0 = v1;
	  AddVV(x2,v2,v0);
	}
      }
    private:
      const T x1;
      const GenVector<T1>& v1;
      const T x2;
      const GenVector<T2>& v2;
  };

  // v+=v
  template <class T> inline const VectorView<T>& operator+=(
      const VectorView<T>& v1, const GenVector<T>& v2) 
  { 
    TMVAssert(v1.size() == v2.size());
    AddVV(T(1),v2,v1); 
    return v1; 
  }

  template <class T> inline const VectorView<CT>& operator+=(
      const VectorView<CT>& v1, const GenVector<T>& v2) 
  {
    TMVAssert(v1.size() == v2.size());
    AddVV(CT(1),v2,v1); 
    return v1; 
  }

  // v-=v
  template <class T> inline const VectorView<T>& operator-=(
      const VectorView<T>& v1, const GenVector<T>& v2)
  {
    TMVAssert(v1.size() == v2.size());
    AddVV(T(-1),v2,v1); 
    return v1; 
  }

  template <class T> inline const VectorView<CT>& operator-=(
      const VectorView<CT>& v1, const GenVector<T>& v2) 
  {
    TMVAssert(v1.size() == v2.size());
    AddVV(CT(-1),v2,v1); 
    return v1; 
  }

  // v+=(x*v)
  template <class T, class T2> inline const VectorView<T>& operator+=(
      const VectorView<T>& v, const ProdXV<T,T2>& pxv)
  {
    TMVAssert(v.size() == pxv.size());
    AddVV(pxv.GetX(),pxv.GetV(),v);
    return v; 
  }

  template <class T> inline const VectorView<CT>& operator+=(
      const VectorView<CT>& v, const ProdXV<T,T>& pxv)
  {
    TMVAssert(v.size() == pxv.size());
    AddVV(CT(pxv.GetX()),pxv.GetV(),v);
    return v; 
  }

  // v-=(x*v)
  template <class T, class T2> inline const VectorView<T>& operator-=(
	const VectorView<T>& v, const ProdXV<T,T2>& pxv)
  {
    TMVAssert(v.size() == pxv.size());
    AddVV(-pxv.GetX(),pxv.GetV(),v);
    return v; 
  }

  template <class T> inline const VectorView<CT>& operator-=(
      const VectorView<CT>& v, const ProdXV<T,T>& pxv)
  {
    TMVAssert(v.size() == pxv.size());
    AddVV(CT(-pxv.GetX()),pxv.GetV(),v);
    return v; 
  }

  //
  // Vector * Vector
  //

  template <class T> inline T operator*(
      const GenVector<T>& v1, const GenVector<T>& v2) 
  { 
    TMVAssert(v1.size() == v2.size());
    return MultVV(v1,v2); 
  }

  template <class T> inline CT operator*(
      const GenVector<CT>& v1, const GenVector<T>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    return MultVV(v1,v2); 
  }

  template <class T> inline CT operator*(
      const GenVector<T>& v1, const GenVector<CT>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    return MultVV(v2,v1); 
  }

  // v * (x*v)

  template <class T, class T2> inline T operator*(
      const GenVector<T>& v1, const ProdXV<T,T2>& v2) 
  { 
    TMVAssert(v1.size() == v2.size());
    return v2.GetX()*MultVV(v1,v2.GetV()); 
  }

  template <class T> inline CT operator*(
      const GenVector<CT>& v1, const ProdXV<T,T>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    return v2.GetX()*MultVV(v1,v2.GetV()); 
  }

  template <class T, class T2> inline CT operator*(
      const GenVector<T>& v1, const ProdXV<CT,T2>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    return v2.GetX()*MultVV(v2.GetV(),v1); 
  }

  // (x*v) * v

  template <class T, class T2> inline T operator*(
      const ProdXV<T,T2>& v1, const GenVector<T>& v2)
  {
    TMVAssert(v1.size() == v2.size());
    return v1.GetX()*MultVV(v1.GetV(),v2); 
  }

  template <class T> inline CT operator*(
      const ProdXV<T,T>& v1, const GenVector<CT>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    return v1.GetX()*MultVV(v2,v1.GetV()); 
  }

  template <class T, class T2> inline CT operator*(
      const ProdXV<CT,T2>& v1, const GenVector<T>& v2)
  { 
    TMVAssert(v1.size() == v2.size());
    return v1.GetX()*MultVV(v1.GetV(),v2); 
  }

  // (x*v) * (x*v)

  template <class T, class T1, class T2> inline T operator*(
      const ProdXV<T,T1>& v1, const ProdXV<T,T2>& v2)
  {
    TMVAssert(v1.size() == v2.size());
    return v1.GetX()*v2.GetX()*MultVV(v1.GetV(),v2.GetV()); 
  }

  template <class T, class T1> inline CT operator*(
      const ProdXV<CT,T1>& v1, const ProdXV<T,T>& v2)
  {
    TMVAssert(v1.size() == v2.size());
    return v1.GetX()*v2.GetX()*MultVV(v1.GetV(),v2.GetV()); 
  }

  template <class T, class T2> inline CT operator*(
      const ProdXV<T,T>& v1, const ProdXV<CT,T2>& v2)
  {
    TMVAssert(v1.size() == v2.size());
    return v1.GetX()*v2.GetX()*MultVV(v1.GetV(),v2.GetV()); 
  }


#define GENVECTOR GenVector
#define GENVECTOR1 GenVector
#define GENVECTOR2 GenVector
#define SUMVV SumVV
#define PRODXV ProdXV
#define PRODXV1 ProdXV
#define PRODXV2 ProdXV

#include "TMV_AuxSumVV.h"
#include "TMV_AuxProdXV.h"

#undef GENVECTOR
#undef GENVECTOR1
#undef GENVECTOR2
#undef SUMVV
#undef PRODXV
#undef PRODXV1
#undef PRODXV2

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif 
