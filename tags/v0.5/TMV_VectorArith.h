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
  template <class T, class Tx> void MultXV(const Tx x, 
      const VectorView<T>& v2);

  // v2 = x * v1 + v2
  template <class T, class Tx, class T1> void AddVV(const Tx x, 
      const GenVector<T1>& v1, const VectorView<T>& v2);
  // v0 = x1*v1 + x2*v2
  template <class T, class T1> void AddVV(const T x1, const GenVector<T1>& v1,
      const T x2, const GenVector<T>& v2, const VectorView<T>& v0);
  template <class T> inline void AddVV(const CT x1, const GenVector<CT>& v1,
      const CT x2, const GenVector<T>& v2, const VectorView<CT>& v0)
  { AddVV(x2,v2,x1,v1,v0); }
  template <class T> inline void AddVV(const CT x1, const GenVector<T>& v1,
      const CT x2, const GenVector<T>& v2, const VectorView<CT>& v0)
  { 
    AddVV(REAL(x1),v1,REAL(x2),v2,v0.Real()); 
    AddVV(IMAG(x1),v1,IMAG(x2),v2,v0.Imag()); 
  }

  // v1 * v2 (dot product)
  // Note: the return type is the type of the first vector
  // This is important for mixing complex and real vectors
  template <class T, class T2> T MultVV(
      const GenVector<T>& v1, const GenVector<T2>& v2);
  template <class T> complex<T> MultVV(
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
      VectorComposite() : GenVector<T>(NonConj), inst(0) {}
      virtual ~VectorComposite() { if (inst) delete inst; };
      virtual void AssignTo(const VectorView<T>& v0) const = 0;
    protected:
      mutable const Vector<T>* inst;
      inline const T* cptr() const
      {
	if (!inst) inst = new Vector<T>(*this); 
	return inst->cptr(); 
      }
      inline int step() const { return 1; }
  };

  // These are what we want to do no matter what type T2 is:
  template <class T1, class T2> inline Vector<T1>& operator+=(
      Vector<T1>& v1, const T2& x2)
  { v1.View() += x2; return v1; }

  template <class T1, class T2> inline Vector<T1>& operator-=(
      Vector<T1>& v1, const T2& x2)
  { v1.View() -= x2; return v1; }

  template <class T1, class T2> inline Vector<T1>& operator*=(
      Vector<T1>& v1, const T2& x2)
  { v1.View() *= x2; return v1; }

  template <class T1, class T2> inline Vector<T1>& operator/=(
      Vector<T1>& v1, const T2& x2)
  { v1.View() /= x2; return v1; }

  template <class T1, class T2> inline Vector<T1>& operator%=(
      Vector<T1>& v1, const T2& x2)
  { v1.View() %= x2; return v1; }


  //
  // Vector * / Scalar
  //

  template <class T, class T2> class ProdXV : public VectorComposite<T> 
  {
    public:
      ProdXV(T _x1, const GenVector<T2>& _v2) : x1(_x1), v2(&_v2) {}
      inline size_t size() const { return v2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenVector<T2>& GetV2() const { return *v2; }
      inline void AssignTo(const VectorView<T>& v0) const
      { DoAssignTo(v0); }
      inline void DoAssignTo(const VectorView<T>& v0) const
      {
	TMVAssert(v0.size() == size());
	if (!(v2->SameAs(v0))) v0 = *v2;
	MultXV(x1,v0);
      }
    private:
      const T x1;
      const GenVector<T2>*const v2;
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

  template <class T, class T2, class T4> class SumVV : public VectorComposite<T> 
  {
    public:
      SumVV(T _x1, const GenVector<T2>& _v2, T _x3, const GenVector<T4>& _v4) :
	x1(_x1),v2(&_v2),x3(_x3), v4(&_v4)
      { TMVAssert(v2->size() == v4->size()); }
      inline size_t size() const { return v2->size(); }
      inline T GetX1() const { return x1; }
      inline const GenVector<T2>& GetV2() const { return *v2; }
      inline T GetX3() const { return x3; }
      inline const GenVector<T4>& GetV4() const { return *v4; }
      inline void AssignTo(const VectorView<T>& v0) const
      { DoAssignTo(v0); }
      inline void DoAssignTo(const VectorView<T>& v0) const
      {
	TMVAssert(v0.size() == size());
	AddVV(x1,*v2,x3,*v4,v0);
      }
    private:
      const T x1;
      const GenVector<T2>*const v2;
      const T x3;
      const GenVector<T4>*const v4;
  };

  // v+=v
  template <class T> inline const VectorView<T>& operator+=(
      const VectorView<T>& v1, const GenVector<T>& v2) 
  { AddVV(T(1),v2,v1); return v1; }

  template <class T> inline const VectorView<CT>& operator+=(
      const VectorView<CT>& v1, const GenVector<T>& v2) 
  { AddVV(CT(1),v2,v1); return v1; }

  // v-=v
  template <class T> inline const VectorView<T>& operator-=(
      const VectorView<T>& v1, const GenVector<T>& v2)
  { AddVV(T(-1),v2,v1); return v1; }

  template <class T> inline const VectorView<CT>& operator-=(
      const VectorView<CT>& v1, const GenVector<T>& v2) 
  { AddVV(CT(-1),v2,v1); return v1; }

  // v+=(x*v)
  template <class T, class T2> inline const VectorView<T>& operator+=(
      const VectorView<T>& v, const ProdXV<T,T2>& pxv)
  {
    AddVV(pxv.GetX1(),pxv.GetV2(),v);
    return v; 
  }

  template <class T> inline const VectorView<CT>& operator+=(
      const VectorView<CT>& v1, const ProdXV<T,T>& pxv2)
  {
    AddVV(CT(pxv2.GetX1()),pxv2.GetV2(),v1);
    return v1; 
  }

  // v-=(x*v)
  template <class T, class T2> inline const VectorView<T>& operator-=(
	const VectorView<T>& v1, const ProdXV<T,T2>& pxv2)
  {
    AddVV(-pxv2.GetX1(),pxv2.GetV2(),v1);
    return v1; 
  }

  template <class T> inline const VectorView<CT>& operator-=(
      const VectorView<CT>& v1, const ProdXV<T,T>& pxv2)
  {
    AddVV(CT(-pxv2.GetX1()),pxv2.GetV2(),v1);
    return v1; 
  }

  //
  // Vector * Vector
  //

  template <class T> inline T operator*(
      const GenVector<T>& v1, const GenVector<T>& v2) 
  { return MultVV(v1,v2); }

  template <class T> inline CT operator*(
      const GenVector<CT>& v1, const GenVector<T>& v2)
  { return MultVV(v1,v2); }

  template <class T> inline CT operator*(
      const GenVector<T>& v1, const GenVector<CT>& v2)
  { return MultVV(v2,v1); }

  // v * (x*v)

  template <class T, class T2> inline T operator*(
      const GenVector<T>& v1, const ProdXV<T,T2>& v2) 
  { return v2.GetX1()*MultVV(v1,v2.GetV2()); }

  template <class T> inline CT operator*(
      const GenVector<CT>& v1, const ProdXV<T,T>& v2)
  { return v2.GetX1()*MultVV(v1,v2.GetV2()); }

  template <class T, class T2> inline CT operator*(
      const GenVector<T>& v1, const ProdXV<CT,T2>& v2)
  { return v2.GetX1()*MultVV(v2.GetV2(),v1); }

  // (x*v) * v

  template <class T, class T2> inline T operator*(
      const ProdXV<T,T2>& v1, const GenVector<T>& v2)
  { return v1.GetX1()*MultVV(v1.GetV2(),v2); }

  template <class T> inline CT operator*(
      const ProdXV<T,T>& v1, const GenVector<CT>& v2)
  { return v1.GetX1()*MultVV(v2,v1.GetV2()); }

  template <class T, class T2> inline CT operator*(
      const ProdXV<CT,T2>& v1, const GenVector<T>& v2)
  { return v1.GetX1()*MultVV(v1.GetV2(),v2); }

  // (x*v) * (x*v)

  template <class T, class T2, class T3> inline T operator*(
      const ProdXV<T,T2>& v1, const ProdXV<T,T3>& v2)
  { return v1.GetX1()*v2.GetX1()*MultVV(v1.GetV2(),v2.GetV2()); }

  template <class T, class T2> inline CT operator*(
      const ProdXV<CT,T2>& v1, const ProdXV<T,T>& v2)
  { return v1.GetX1()*v2.GetX1()*MultVV(v1.GetV2(),v2.GetV2()); }

  template <class T, class T3> inline CT operator*(
      const ProdXV<T,T>& v1, const ProdXV<CT,T3>& v2)
  { return v1.GetX1()*v2.GetX1()*MultVV(v1.GetV2(),v2.GetV2()); }


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
