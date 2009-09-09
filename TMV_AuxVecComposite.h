// This file sets up the Composite classes for all operations with a 
// (sparse) matrix that returns a vector
//
// Need to define the following with #define statements.
// (The given definition is for a Band Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX2 GenBandMatrix
// #define PRODMV ProdBV
// #define QUOTVM QuotVB
// #define RQUOTVM RQuotVB

//
// Matrix * Vector
//

template <class T, class T2, bool transm2, class T3> 
class PRODMV : public VectorComposite<T>
{
  public:
    PRODMV(const T _x1, const GENMATRIX2<T2>& _m2,
	const GenVector<T3>& _v3) :
      x1(_x1), m2(&_m2), v3(&_v3)
    { TMVAssert(v3->size()==(transm2 ?  m2->colsize() : m2->rowsize())); }
    inline size_t size() const 
    { return transm2 ? m2->rowsize() : m2->colsize(); }
    inline T GetX1() const { return x1; }
    inline const GENMATRIX2<T2>& GetM2() const { return *m2; }
    inline const GenVector<T3>& GetV3() const { return *v3; }
    inline void AssignTo(const VectorView<T>& v0) const
    {
      TMVAssert(v0.size() == size());
      if (transm2)
	MultMV(x1,m2->QuickTranspose(),*v3,T(0),v0);
      else
	MultMV(x1,*m2,*v3,T(0),v0);
    }
  private:
    const T x1;
    const GENMATRIX2<T2>*const m2;
    const GenVector<T3>*const v3;
};

template <class T> inline const VectorView<T>& operator*=(
    const VectorView<T>& v1, const GENMATRIX2<T>& m2)
{ 
  MultMV(T(1),m2.QuickTranspose(),v1,T(0),v1);
  return v1;
}

template <class T> inline const VectorView<CT>& operator*=(
    const VectorView<CT>& v1, const GENMATRIX2<T>& m2)
{
  MultMV(CT(1),m2.QuickTranspose(),v1,CT(0),v1);
  return v1;
}

template <class T, class T2, bool transm2, class T3> 
inline const VectorView<T>& operator+=(
    const VectorView<T>& v, const PRODMV<T,T2,transm2,T3>& pmv)
{ 
  if (transm2)
    MultMV(pmv.GetX1(),pmv.GetM2().QuickTranspose(),pmv.GetV3(),T(1),v);
  else
    MultMV(pmv.GetX1(),pmv.GetM2(),pmv.GetV3(),T(1),v);
  return v;
}

template <class T, bool transm2> inline const VectorView<CT>& operator+=(
    const VectorView<CT>& v, const PRODMV<T,T,transm2,T>& pmv)
{ 
  if (transm2)
    MultMV(CT(pmv.GetX1()),pmv.GetM2().QuickTranspose(),pmv.GetV3(),CT(1),v);
  else
    MultMV(CT(pmv.GetX1()),pmv.GetM2(),pmv.GetV3(),CT(1),v);
  return v;
}

template <class T, class T2, bool transm2, class T3> 
inline const VectorView<T>& operator-=(
    const VectorView<T>& v, const PRODMV<T,T2,transm2,T3>& pmv)
{ 
  if (transm2)
    MultMV(-pmv.GetX1(), pmv.GetM2().QuickTranspose(), pmv.GetV3(), T(1), v);
  else
    MultMV(-pmv.GetX1(), pmv.GetM2(), pmv.GetV3(), T(1), v);
  return v;
}

template <class T, bool transm2> inline const VectorView<CT>& operator-=(
    const VectorView<CT>& v, const PRODMV<T,T,transm2,T>& pmv)
{ 
  if (transm2)
    MultMV(CT(-pmv.GetX1()),pmv.GetM2().QuickTranspose(),pmv.GetV3(),CT(1),v);
  else
    MultMV(CT(-pmv.GetX1()),pmv.GetM2(),pmv.GetV3(),CT(1),v);
  return v;
}

//
// Vector / % Matrix 
// v/m is the solution (x) of mx = v
// ie. / is really division from the left: x = m^-1 v
// Use % if you want division from the right (v m^-1)
//

template <class T0, class T1, class T2> class QUOTVM : 
  public VectorComposite<T0>
{
  public:
    QUOTVM(const GenVector<T1>& _v1, const GENMATRIX2<T2>& _m2) :
      v1(&_v1), m2(&_m2)
    { TMVAssert(v1->size()==m2->colsize()); }
    inline size_t size() const { return m2->rowsize(); }
    inline const GenVector<T1>& GetV1() const { return *v1; }
    inline const GENMATRIX2<T2>& GetM2() const { return *m2; }
    inline void AssignTo(const VectorView<T0>& v0) const
    {
      TMVAssert(v0.size() == size());
      m2->LDiv(*v1,v0);
    }
  private:
    const GenVector<T1>*const v1;
    const GENMATRIX2<T2>*const m2;
};

template <class T0, class T1, class T2> class RQUOTVM : 
  public VectorComposite<T0>
{
  public:
    RQUOTVM(const GenVector<T1>& _v1, const GENMATRIX2<T2>& _m2) :
      v1(&_v1), m2(&_m2)
    { TMVAssert(v1->size()==m2->rowsize()); }
    inline size_t size() const { return m2->colsize(); }
    inline const GenVector<T1>& GetV1() const { return *v1; }
    inline const GENMATRIX2<T2>& GetM2() const { return *m2; }
    inline void AssignTo(const VectorView<T0>& v0) const
    {
      TMVAssert(v0.size() == size());
      m2->RDiv(*v1,v0);
    }
  private:
    const GenVector<T1>*const v1;
    const GENMATRIX2<T2>*const m2;
};

template <class T> inline const VectorView<T>& operator/=(
    const VectorView<T>& v1, const GENMATRIX2<T>& m2)
{ m2.LDivEq(v1); return v1; }

template <class T> inline const VectorView<CT>& operator/=(
    const VectorView<CT>& v1, const GENMATRIX2<T>& m2)
{ m2.LDivEq(v1); return v1; }

template <class T> inline const VectorView<T>& operator%=(
    const VectorView<T>& v1, const GENMATRIX2<T>& m2)
{ m2.RDivEq(v1); return v1; }

template <class T> inline const VectorView<CT>& operator%=(
    const VectorView<CT>& v1, const GENMATRIX2<T>& m2)
{ m2.RDivEq(v1); return v1; }

