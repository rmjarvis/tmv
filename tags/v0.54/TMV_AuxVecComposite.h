// This file sets up the Composite classes for all operations with a 
// (sparse) matrix that returns a vector
//
// Need to define the following with #define statements.
// (The given definition is for a Band Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX GenBandMatrix
// #define PRODMV ProdBV
// #define QUOTVM QuotVB
// #define RQUOTVM RQuotVB

//
// Matrix * Vector
//

template <class T, class T1, bool transm, class T2> 
class PRODMV : public VectorComposite<T>
{
  public:
    PRODMV(const T _x, const GENMATRIX<T1>& _m,
	const GenVector<T2>& _v) :
      x(_x), m(&_m), v(&_v)
    { TMVAssert(v->size()==(transm ?  m->colsize() : m->rowsize())); }
    inline size_t size() const 
    { return transm ? m->rowsize() : m->colsize(); }
    inline T GetX() const { return x; }
    inline const GENMATRIX<T1>& GetM() const { return *m; }
    inline const GenVector<T2>& GetV() const { return *v; }
    inline void AssignTo(const VectorView<T>& v0) const
    {
      TMVAssert(v0.size() == size());
      if (transm)
	MultMV(x,m->Transpose(),*v,T(0),v0);
      else
	MultMV(x,*m,*v,T(0),v0);
    }
  private:
    const T x;
    const GENMATRIX<T1>*const m;
    const GenVector<T2>*const v;
};

template <class T> inline const VectorView<T>& operator*=(
    const VectorView<T>& v, const GENMATRIX<T>& m)
{ 
  MultMV(T(1),m.Transpose(),v,T(0),v);
  return v;
}

template <class T> inline const VectorView<CT>& operator*=(
    const VectorView<CT>& v, const GENMATRIX<T>& m)
{
  MultMV(CT(1),m.Transpose(),v,CT(0),v);
  return v;
}

template <class T, class T1, bool transm, class T2> 
inline const VectorView<T>& operator+=(
    const VectorView<T>& v, const PRODMV<T,T1,transm,T2>& pmv)
{ 
  if (transm)
    MultMV(pmv.GetX(),pmv.GetM().Transpose(),pmv.GetV(),T(1),v);
  else
    MultMV(pmv.GetX(),pmv.GetM(),pmv.GetV(),T(1),v);
  return v;
}

template <class T, bool transm> inline const VectorView<CT>& operator+=(
    const VectorView<CT>& v, const PRODMV<T,T,transm,T>& pmv)
{ 
  if (transm)
    MultMV(CT(pmv.GetX()),pmv.GetM().Transpose(),pmv.GetV(),CT(1),v);
  else
    MultMV(CT(pmv.GetX()),pmv.GetM(),pmv.GetV(),CT(1),v);
  return v;
}

template <class T, class T1, bool transm, class T2> 
inline const VectorView<T>& operator-=(
    const VectorView<T>& v, const PRODMV<T,T1,transm,T2>& pmv)
{ 
  if (transm)
    MultMV(-pmv.GetX(), pmv.GetM().Transpose(), pmv.GetV(), T(1), v);
  else
    MultMV(-pmv.GetX(), pmv.GetM(), pmv.GetV(), T(1), v);
  return v;
}

template <class T, bool transm> inline const VectorView<CT>& operator-=(
    const VectorView<CT>& v, const PRODMV<T,T,transm,T>& pmv)
{ 
  if (transm)
    MultMV(CT(-pmv.GetX()),pmv.GetM().Transpose(),pmv.GetV(),CT(1),v);
  else
    MultMV(CT(-pmv.GetX()),pmv.GetM(),pmv.GetV(),CT(1),v);
  return v;
}

#define GENVECTOR GenVector
#define PRODXV ProdXV
#include "TMV_AuxProdMV.h"
// Defines things like v*m, m*v, (x*m)*(x*v), etc.
#include "TMV_AuxProdMVa.h"
// Defines things like -(m*v), x*(m*v), (m*v)/x, etc.

//
// Vector / % Matrix 
// v/m is the solution (x) of mx = v
// ie. / is really division from the left: x = m^-1 v
// Use % if you want division from the right (v m^-1)
//

template <class T, class T1, class T2> class QUOTVM : 
  public VectorComposite<T>
{
  public:
    QUOTVM(const T _x, const GenVector<T1>& _v, const GENMATRIX<T2>& _m) :
      x(_x), v(&_v), m(&_m)
    { TMVAssert(v->size()==m->colsize()); }
    inline size_t size() const { return m->rowsize(); }
    inline T GetX() const { return x; }
    inline const GenVector<T1>& GetV() const { return *v; }
    inline const GENMATRIX<T2>& GetM() const { return *m; }
    inline void AssignTo(const VectorView<T>& v0) const
    {
      TMVAssert(v0.size() == size());
      if (v0.SameAs(*v)) m->LDivEq(v0);
      else m->LDiv(*v,v0);
      if (x != T(1)) v0 *= x;
    }
  private:
    const T x;
    const GenVector<T1>*const v;
    const GENMATRIX<T2>*const m;
};

template <class T, class T1, class T2> class RQUOTVM : 
  public VectorComposite<T>
{
  public:
    RQUOTVM(const T _x, const GenVector<T1>& _v, const GENMATRIX<T2>& _m) :
      x(_x), v(&_v), m(&_m)
    { TMVAssert(v->size()==m->rowsize()); }
    inline size_t size() const { return m->colsize(); }
    inline T GetX() const { return x; }
    inline const GenVector<T1>& GetV() const { return *v; }
    inline const GENMATRIX<T2>& GetM() const { return *m; }
    inline void AssignTo(const VectorView<T>& v0) const
    {
      TMVAssert(v0.size() == size());
      if (v0.SameAs(*v)) m->RDivEq(v0);
      else m->RDiv(*v,v0);
      if (x != T(1)) v0 *= x;
    }
  private:
    const T x;
    const GenVector<T1>*const v;
    const GENMATRIX<T2>*const m;
};

template <class T> inline const VectorView<T>& operator/=(
    const VectorView<T>& v, const GENMATRIX<T>& m)
{ 
  TMVAssert(m.IsSquare());
  TMVAssert(m.rowsize() == v.size());
  m.LDivEq(v); 
  return v; 
}

template <class T> inline const VectorView<CT>& operator/=(
    const VectorView<CT>& v, const GENMATRIX<T>& m)
{
  TMVAssert(m.IsSquare());
  TMVAssert(m.rowsize() == v.size());
  m.LDivEq(v); 
  return v; 
}

template <class T> inline const VectorView<T>& operator%=(
    const VectorView<T>& v, const GENMATRIX<T>& m)
{
  TMVAssert(m.IsSquare());
  TMVAssert(m.rowsize() == v.size());
  m.RDivEq(v); 
  return v; 
}

template <class T> inline const VectorView<CT>& operator%=(
    const VectorView<CT>& v, const GENMATRIX<T>& m)
{ 
  TMVAssert(m.IsSquare());
  TMVAssert(m.rowsize() == v.size());
  m.RDivEq(v); 
  return v; 
}

template <class T, class Tm> inline const VectorView<T>& operator*=(
    const VectorView<T>& v, const QUOTXM<T,Tm>& qxm)
{
  TMVAssert(qxm.GetM().IsSquare());
  TMVAssert(qxm.GetM().rowsize() == v.size());
  qxm.GetM().RDivEq(v); 
  v *= qxm.GetX(); 
  return v; 
}

template <class T> inline const VectorView<CT>& operator*=(
    const VectorView<CT>& v, const QUOTXM<T,T>& qxm)
{
  TMVAssert(qxm.GetM().IsSquare());
  TMVAssert(qxm.GetM().rowsize() == v.size());
  qxm.GetM().RDivEq(v); 
  v *= qxm.GetX(); 
  return v; 
}

#include "TMV_AuxQuotVM.h"

#undef GENVECTOR
#undef PRODXV


