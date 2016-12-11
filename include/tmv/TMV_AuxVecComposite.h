///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


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

template <typename T, typename T1, typename T2>
class PRODMV : public VectorComposite<T>
{
public:
    inline PRODMV(
        const T _x, const GENMATRIX<T1>& _m, const GenVector<T2>& _v) :
        x(_x), m(_m), v(_v)
    { TMVAssert(v.size()==m.rowsize()); }
    inline ptrdiff_t size() const { return m.colsize(); }
    inline T getX() const { return x; }
    inline const GENMATRIX<T1>& getM() const { return m; }
    inline const GenVector<T2>& getV() const { return v; }
    inline void assignToV(VectorView<TMV_RealType(T)> v0) const
    {
        TMVAssert(isReal(T()));
        TMVAssert(v0.size() == size());
        MultMV<false>(x,m,v,v0);
    }
    inline void assignToV(VectorView<TMV_ComplexType(T)> v0) const
    {
        TMVAssert(v0.size() == size());
        MultMV<false>(TMV_ComplexType(T)(x),m,v,v0);
    }
private:
    const T x;
    const GENMATRIX<T1>& m;
    const GenVector<T2>& v;
};

template <typename T, typename T1, typename T2>
inline VectorView<T> operator+=(
    VectorView<T> v, const PRODMV<T,T1,T2>& pmv)
{ MultMV<true>(pmv.getX(),pmv.getM(),pmv.getV(),v); return v; }

template <typename T>
inline VectorView<CT> operator+=(
    VectorView<CT> v, const PRODMV<T,T,T>& pmv)
{ MultMV<true>(CT(pmv.getX()),pmv.getM(),pmv.getV(),v); return v; }

template <typename T, typename T1, typename T2>
inline VectorView<T> operator-=(
    VectorView<T> v, const PRODMV<T,T1,T2>& pmv)
{ MultMV<true>(-pmv.getX(), pmv.getM(), pmv.getV(), v); return v; }

template <typename T>
inline VectorView<CT> operator-=(
    VectorView<CT> v, const PRODMV<T,T,T>& pmv)
{ MultMV<true>(CT(-pmv.getX()),pmv.getM(),pmv.getV(),v); return v; }


template <typename T, typename T1, typename T2>
class PRODVM : public VectorComposite<T>
{
public:
    inline PRODVM(
        const T _x, const GenVector<T1>& _v, const GENMATRIX<T2>& _m) :
        x(_x), v(_v), m(_m)
    { TMVAssert(v.size()==m.colsize()); }
    inline ptrdiff_t size() const { return m.rowsize(); }
    inline T getX() const { return x; }
    inline const GenVector<T1>& getV() const { return v; }
    inline const GENMATRIX<T2>& getM() const { return m; }
    inline void assignToV(VectorView<TMV_RealType(T)> v0) const
    {
        TMVAssert(isReal(T()));
        TMVAssert(v0.size() == size());
        MultMV<false>(x,m.transpose(),v,v0);
    }
    inline void assignToV(VectorView<TMV_ComplexType(T)> v0) const
    {
        TMVAssert(v0.size() == size());
        MultMV<false>(TMV_ComplexType(T)(x),m.transpose(),v,v0);
    }
private:
    const T x;
    const GenVector<T1>& v;
    const GENMATRIX<T2>& m;
};

template <typename T>
inline VectorView<T> operator*=(
    VectorView<T> v, const GENMATRIX<T>& m)
{ MultMV<false>(T(1),m.transpose(),v,v); return v; }

template <typename T>
inline VectorView<CT> operator*=(
    VectorView<CT> v, const GENMATRIX<T>& m)
{ MultMV<false>(T(1),m.transpose(),v,v); return v; }

template <typename T, typename T1, typename T2>
inline VectorView<T> operator+=(
    VectorView<T> v, const PRODVM<T,T1,T2>& pmv)
{ MultMV<true>(pmv.getX(),pmv.getM().transpose(),pmv.getV(),v); return v; }

template <typename T>
inline VectorView<CT> operator+=(
    VectorView<CT> v, const PRODVM<T,T,T>& pmv)
{ MultMV<true>(pmv.getX(),pmv.getM().transpose(),pmv.getV(),v); return v; }

template <typename T, typename T1, typename T2>
inline VectorView<T> operator-=(
    VectorView<T> v, const PRODVM<T,T1,T2>& pmv)
{ MultMV<true>(-pmv.getX(), pmv.getM().transpose(), pmv.getV(), v); return v; }

template <typename T>
inline VectorView<CT> operator-=(
    VectorView<CT> v, const PRODVM<T,T,T>& pmv)
{ MultMV<true>(-pmv.getX(),pmv.getM().transpose(),pmv.getV(),v); return v; }


#define GENMATRIX1 GENMATRIX
#define GENMATRIX2 GenVector
#define PRODMM PRODMV
#define PRODXM1 PRODXM
#define PRODXM2 ProdXV
#define GETM1 .getM()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .getM()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2
#define GENMATRIX1 GenVector
#define GENMATRIX2 GENMATRIX
#define PRODMM PRODVM
#define PRODXM1 ProdXV
#define PRODXM2 PRODXM
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

//
// Vector / % Matrix
// v/m is the solution (x) of mx = v
// ie. / is really division from the left: x = m^-1 v
// Use % if you want division from the right (v m^-1)
//

template <typename T, typename T1, typename T2>
class QUOTVM : public VectorComposite<T>
{
public:
    inline QUOTVM(
        const T _x, const GenVector<T1>& _v, const GENMATRIX<T2>& _m) :
        x(_x), v(_v), m(_m)
    { TMVAssert(v.size()==m.colsize()); }
    inline ptrdiff_t size() const { return m.rowsize(); }
    inline T getX() const { return x; }
    inline const GenVector<T1>& getV() const { return v; }
    inline const GENMATRIX<T2>& getM() const { return m; }
    inline void assignToV(VectorView<TMV_RealType(T)> v0) const
    {
        TMVAssert(isReal(T()));
        TMVAssert(v0.size() == size());
        m.LDiv(v,v0);
        if (x != T(1)) v0 *= TMV_REAL(x);
    }
    inline void assignToV(VectorView<TMV_ComplexType(T)> v0) const
    {
        TMVAssert(v0.size() == size());
        m.LDiv(v,v0);
        if (x != T(1)) v0 *= x;
    }
private:
    const T x;
    const GenVector<T1>& v;
    const GENMATRIX<T2>& m;
};

template <typename T, typename T1, typename T2>
class RQUOTVM : public VectorComposite<T>
{
public:
    inline RQUOTVM(
        const T _x, const GenVector<T1>& _v, const GENMATRIX<T2>& _m) :
        x(_x), v(_v), m(_m)
    { TMVAssert(v.size()==m.rowsize()); }
    inline ptrdiff_t size() const { return m.colsize(); }
    inline T getX() const { return x; }
    inline const GenVector<T1>& getV() const { return v; }
    inline const GENMATRIX<T2>& getM() const { return m; }
    inline void assignToV(VectorView<TMV_RealType(T)> v0) const
    {
        TMVAssert(isReal(T()));
        TMVAssert(v0.size() == size());
        m.RDiv(v,v0);
        if (x != T(1)) v0 *= TMV_REAL(x);
    }
    inline void assignToV(VectorView<TMV_ComplexType(T)> v0) const
    {
        TMVAssert(v0.size() == size());
        m.RDiv(v,v0);
        if (x != T(1)) v0 *= x;
    }
private:
    const T x;
    const GenVector<T1>& v;
    const GENMATRIX<T2>& m;
};

template <typename T>
inline VectorView<T> operator/=(
    VectorView<T> v, const GENMATRIX<T>& m)
{
    TMVAssert(m.isSquare());
    TMVAssert(m.rowsize() == v.size());
    m.LDivEq(v);
    return v;
}

template <typename T>
inline VectorView<CT> operator/=(
    VectorView<CT> v, const GENMATRIX<T>& m)
{
    TMVAssert(m.isSquare());
    TMVAssert(m.rowsize() == v.size());
    m.LDivEq(v);
    return v;
}

template <typename T>
inline VectorView<T> operator%=(
    VectorView<T> v, const GENMATRIX<T>& m)
{
    TMVAssert(m.isSquare());
    TMVAssert(m.rowsize() == v.size());
    m.RDivEq(v);
    return v;
}

template <typename T>
inline VectorView<CT> operator%=(
    VectorView<CT> v, const GENMATRIX<T>& m)
{
    TMVAssert(m.isSquare());
    TMVAssert(m.rowsize() == v.size());
    m.RDivEq(v);
    return v;
}

template <typename T, typename Tm>
inline VectorView<T> operator*=(
    VectorView<T> v, const QUOTXM<T,Tm>& qxm)
{
    TMVAssert(qxm.getM().isSquare());
    TMVAssert(qxm.getM().rowsize() == v.size());
    qxm.getM().RDivEq(v);
    v *= qxm.getX();
    return v;
}

template <typename T>
inline VectorView<CT> operator*=(
    VectorView<CT> v, const QUOTXM<T,T>& qxm)
{
    TMVAssert(qxm.getM().isSquare());
    TMVAssert(qxm.getM().rowsize() == v.size());
    qxm.getM().RDivEq(v);
    v *= qxm.getX();
    return v;
}

#define GENMATRIX1 GenVector
#define GENMATRIX2 GENMATRIX
#define PRODXM1 ProdXV
#define PRODXM2 PRODXM
#define QUOTMM QUOTVM
#define RQUOTMM RQUOTVM
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxQuotMM.h"
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxQuotMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTMM
#undef RQUOTMM

