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
// (sparse) matrix that returns a Matrix and is of the form:
// (M) op (S)
// (ie. the first operand is the normal Matrix.)
//
// Need to define the following with #define statements.
// (The given definition is for a Band Matrix.  Modify as
// appropriate for the various other matrices.)
//
// #define GENMATRIX2 GenBandMatrix
//
// #define PRODXM ProdXB
// #define SUMMM SumMB
// #define PRODMM ProdMB
// #define QUOTMM QuotMB
// #define RQUOTMM RQuotMB
// #define TQUOTMM TransientQuotMB
// #define TRQUOTMM TransientRQuotMB


//
// Matrix + Matrix
//

template <typename T, typename T1, typename T2>
class SUMMM :
    public MatrixComposite<T>
{
public:
    inline SUMMM(
        const T _x1, const GenMatrix<T1>& _m1,
        const T _x2, const GENMATRIX2<T2>& _m2) :
        x1(_x1), m1(_m1), x2(_x2), m2(_m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
    }
    inline ptrdiff_t colsize() const { return m1.colsize(); }
    inline ptrdiff_t rowsize() const { return m1.rowsize(); }
    inline T getX1() const { return x1; }
    inline const GenMatrix<T1>& getM1() const { return m1; }
    inline T getX2() const { return x2; }
    inline const GENMATRIX2<T2>& getM2() const { return m2; }
    inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
    {
        TMVAssert(isReal(T()));
        TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
        AddMM(x1,m1,x2,m2,m0);
    }
    inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
    {
        TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
        AddMM(x1,m1,x2,m2,m0);
    }
private:
    const T x1;
    const GenMatrix<T1>& m1;
    const T x2;
    const GENMATRIX2<T2>& m2;
};

template <typename T>
inline MatrixView<T> operator+=(
    MatrixView<T> m1, const GENMATRIX2<T>& m2)
{
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    AddMM(T(1),m2,m1);
    return m1;
}

template <typename T>
inline MatrixView<CT> operator+=(
    MatrixView<CT> m1, const GENMATRIX2<T>& m2)
{
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    AddMM(CT(1),m2,m1);
    return m1;
}

template <typename T>
inline MatrixView<T> operator-=(
    MatrixView<T> m1, const GENMATRIX2<T>& m2)
{
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    AddMM(T(-1),m2,m1);
    return m1;
}

template <typename T>
inline MatrixView<CT> operator-=(
    MatrixView<CT> m1, const GENMATRIX2<T>& m2)
{
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    AddMM(CT(-1),m2,m1);
    return m1;
}

template <typename T, typename T2>
inline MatrixView<T> operator+=(
    MatrixView<T> m, const PRODXM2<T,T2>& pxm)
{
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    AddMM(pxm.getX(),pxm.getM(),m);
    return m;
}

template <typename T>
inline MatrixView<CT> operator+=(
    MatrixView<CT> m, const PRODXM2<T,T>& pxm)
{
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    AddMM(CT(pxm.getX()),pxm.getM(),m);
    return m;
}

template <typename T, typename T2>
inline MatrixView<T> operator-=(
    MatrixView<T> m, const PRODXM2<T,T2>& pxm)
{
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    AddMM(-pxm.getX(),pxm.getM(),m);
    return m;
}

template <typename T>
inline MatrixView<CT> operator-=(
    MatrixView<CT> m, const PRODXM2<T,T>& pxm)
{
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    AddMM(CT(-pxm.getX()),pxm.getM(),m);
    return m;
}

#include "tmv/TMV_AuxSumMM.h"
// Defines things like m+m, m-m, (x*m)-(x*m), etc.
#include "tmv/TMV_AuxSumMMa.h"
// Defines things like -(m+m), x*(m+m), (m+m)/x, etc.

//
// Matrix * Matrix
//

template <typename T, typename T1, typename T2>
class PRODMM :
    public MatrixComposite<T>
{
public:
    inline PRODMM(
        const T _x, const GenMatrix<T1>& _m1,
        const GENMATRIX2<T2>& _m2) : x(_x), m1(_m1), m2(_m2)
    { TMVAssert(m1.rowsize() == m2.colsize()) ; }
    inline ptrdiff_t colsize() const { return m1.colsize(); }
    inline ptrdiff_t rowsize() const { return m2.rowsize(); }
    inline T getX() const { return x; }
    inline const GenMatrix<T1>& getM1() const { return m1; }
    inline const GENMATRIX2<T2>& getM2() const { return m2; }
    inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
    {
        TMVAssert(isReal(T()));
        TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
        MultMM<false>(x, m1, m2, m0);
    }
    inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
    {
        TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
        MultMM<false>(x, m1, m2, m0);
    }
private:
    const T x;
    const GenMatrix<T1>& m1;
    const GENMATRIX2<T2>& m2;
};

template <typename T>
inline MatrixView<T> operator*=(
    MatrixView<T> m1, const GENMATRIX2<T>& m2)
{
    TMVAssert(m2.colsize()==m2.rowsize());
    TMVAssert(m1.rowsize()==m2.rowsize());
    MultMM<false>(T(1),m1,m2,m1);
    return m1;
}

template <typename T>
inline MatrixView<CT> operator*=(
    MatrixView<CT> m1, const GENMATRIX2<T>& m2)
{
    TMVAssert(m2.colsize()==m2.rowsize());
    TMVAssert(m1.rowsize()==m2.rowsize());
    MultMM<false>(CT(1),m1,m2,m1);
    return m1;
}

template <typename T, typename T2>
inline MatrixView<T> operator*=(
    MatrixView<T> m1, const PRODXM2<T,T2>& pxm)
{
    TMVAssert(pxm.colsize()==pxm.rowsize());
    TMVAssert(m1.rowsize()==pxm.rowsize());
    MultMM<false>(pxm.getX(),m1,pxm.getM(),m1);
    return m1;
}

template <typename T>
inline MatrixView<CT> operator*=(
    MatrixView<CT> m1, const PRODXM2<T,T>& pxm)
{
    TMVAssert(pxm.colsize()==pxm.rowsize());
    TMVAssert(m1.rowsize()==pxm.rowsize());
    MultMM<false>(CT(pxm.getX()),m1,pxm.getM(),m1);
    return m1;
}

template <typename T, typename T2, typename T3>
inline MatrixView<T> operator+=(
    MatrixView<T> m, const PRODMM<T,T2,T3>& pmm)
{
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
    return m;
}

template <typename T>
inline MatrixView<CT> operator+=(
    MatrixView<CT> m, const PRODMM<T,T,T>& pmm)
{
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(CT(pmm.getX()),pmm.getM1(),pmm.getM2(),m);
    return m;
}

template <typename T, typename T2, typename T3>
inline MatrixView<T> operator-=(
    MatrixView<T> m, const PRODMM<T,T2,T3>& pmm)
{
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
    return m;
}

template <typename T>
inline MatrixView<CT> operator-=(
    MatrixView<CT> m, const PRODMM<T,T,T>& pmm)
{
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(CT(-pmm.getX()),pmm.getM1(),pmm.getM2(),m);
    return m;
}

#include "tmv/TMV_AuxProdMM.h"
// Defines things like m*m, m*(x*m), etc.
#include "tmv/TMV_AuxProdMMa.h"
// Defines things like -(m*m), x*(m*m), etc.

//
// Matrix / % Matrix
//

template <typename T, typename T1, typename T2>
class QUOTMM : public MatrixComposite<T>
{
public:
    inline QUOTMM(
        const T _x, const GenMatrix<T1>& _m1, const GENMATRIX2<T2>& _m2) :
        x(_x), m1(_m1), m2(_m2)
    { TMVAssert( m1.colsize() == m2.colsize() ); }
    inline ptrdiff_t colsize() const { return m2.rowsize(); }
    inline ptrdiff_t rowsize() const { return m1.rowsize(); }
    inline T getX() const { return x; }
    inline const GenMatrix<T1>& getM1() const { return m1; }
    inline const GENMATRIX2<T2>& getM2() const { return m2; }
    inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
    {
        TMVAssert(isReal(T()));
        TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
        m2.LDiv(m1,m0);
        MultXM(x,m0);
    }
    inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
    {
        TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
        m2.LDiv(m1,m0);
        MultXM(x,m0);
    }
protected:
    const T x;
    const GenMatrix<T1>& m1;
    const GENMATRIX2<T2>& m2;
};

template <typename T, typename T1, typename T2>
class TQUOTMM : public QUOTMM<T,T1,T2>
{
public :
    inline TQUOTMM(
        const T x, auto_ptr<Matrix<T1,ColMajor> > m1,
        const GENMATRIX2<T2>& m2) : QUOTMM<T,T1,T2>(x,*m1,m2), m1p(m1.release()) {}
    inline TQUOTMM(const TQUOTMM<T,T1,T2>& rhs) :
        QUOTMM<T,T1,T2>(rhs), m1p(rhs.m1p.release()) {}
    inline ~TQUOTMM() {}
    inline auto_ptr<Matrix<T1,ColMajor> > getP() const { return m1p; }

private :
    mutable auto_ptr<Matrix<T1,ColMajor> > m1p;
};

template <typename T, typename T1, typename T2>
class RQUOTMM :
    public MatrixComposite<T>
{
public:
    inline RQUOTMM(
        const T _x, const GenMatrix<T1>& _m1, const GENMATRIX2<T2>& _m2) :
        x(_x), m1(_m1), m2(_m2)
    { TMVAssert( m1.rowsize() == m2.rowsize() ); }
    inline ptrdiff_t colsize() const { return m1.colsize(); }
    inline ptrdiff_t rowsize() const { return m2.colsize(); }
    inline T getX() const { return x; }
    inline const GenMatrix<T1>& getM1() const { return m1; }
    inline const GENMATRIX2<T2>& getM2() const { return m2; }
    inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
    {
        TMVAssert(isReal(T()));
        TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
        m2.RDiv(m1,m0);
        if (x != T(1)) m0 *= TMV_REAL(x);
    }
    inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
    {
        TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
        m2.RDiv(m1,m0);
        if (x != T(1)) m0 *= x;
    }
protected:
    const T x;
    const GenMatrix<T1>& m1;
    const GENMATRIX2<T2>& m2;
};

template <typename T, typename T1, typename T2>
class TRQUOTMM : public RQUOTMM<T,T1,T2>
{
public :
    inline TRQUOTMM(
        const T x, auto_ptr<Matrix<T1,RowMajor> > m1,
        const GENMATRIX2<T2>& m2) : RQUOTMM<T,T1,T2>(x,*m1,m2), m1p(m1.release()) {}
    inline TRQUOTMM(const TRQUOTMM<T,T1,T2>& rhs) :
        RQUOTMM<T,T1,T2>(rhs), m1p(rhs.m1p.release()) {}
    inline ~TRQUOTMM() {}
    inline auto_ptr<Matrix<T1,RowMajor> > getP() const { return m1p; }

private :
    mutable auto_ptr<Matrix<T1,RowMajor> > m1p;
};

template <typename T>
inline MatrixView<T> operator/=(
    MatrixView<T> m1, const GENMATRIX2<T>& m2)
{
    TMVAssert(m2.colsize() == m2.rowsize());
    TMVAssert(m1.colsize() == m2.rowsize());
    m2.LDivEq(m1);
    return m1;
}

template <typename T>
inline MatrixView<CT> operator/=(
    MatrixView<CT> m1, const GENMATRIX2<T>& m2)
{
    TMVAssert(m2.colsize() == m2.rowsize());
    TMVAssert(m1.colsize() == m2.rowsize());
    m2.LDivEq(m1);
    return m1;
}

template <typename T>
inline MatrixView<T> operator%=(
    MatrixView<T> m1, const GENMATRIX2<T>& m2)
{
    TMVAssert(m2.colsize() == m2.rowsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    m2.RDivEq(m1);
    return m1;
}

template <typename T>
inline MatrixView<CT> operator%=(
    MatrixView<CT> m1, const GENMATRIX2<T>& m2)
{
    TMVAssert(m2.colsize() == m2.rowsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    m2.RDivEq(m1);
    return m1;
}


template <typename T, typename Tm>
inline MatrixView<T> operator*=(
    MatrixView<T> m1, const QUOTXM<T,Tm>& qxm)
{
    TMVAssert(qxm.getM().isSquare());
    TMVAssert(m1.rowsize() == qxm.getM().rowsize());
    qxm.getM().RDivEq(m1);
    m1 *= qxm.getX();
    return m1;
}

template <typename T>
inline MatrixView<CT> operator*=(
    MatrixView<CT> m1, const QUOTXM<T,T>& qxm)
{
    TMVAssert(qxm.getM().isSquare());
    TMVAssert(m1.rowsize() == qxm.getM().rowsize());
    qxm.getM().RDivEq(m1);
    m1 *= qxm.getX();
    return m1;
}

#include "tmv/TMV_AuxQuotMM.h"
// Defines things like m/m, m%m, (x*m)/m, m/(x*m), etc.
#include "tmv/TMV_AuxQuotMMa.h"
#include "tmv/TMV_AuxTQuotMMa.h"
// Defines things like (m/m)*x, -(m/m), etc.
