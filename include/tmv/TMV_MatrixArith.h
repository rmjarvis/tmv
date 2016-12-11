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


#ifndef TMV_MatrixArith_H
#define TMV_MatrixArith_H

#include "tmv/TMV_BaseMatrix.h"
#include "tmv/TMV_VectorArithFunc.h"
#include "tmv/TMV_MatrixArithFunc.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

    template <typename T, typename Tv>
    class ProdXV;

    template <typename T, int A, typename Tx>
    inline Matrix<T,A>& operator+=(Matrix<T,A>& m, const Tx& x)
    { m.view() += x; return m; }

    template <typename T, int A, typename Tx>
    inline Matrix<T,A>& operator-=(Matrix<T,A>& m, const Tx& x)
    { m.view() -= x; return m; }

    template <typename T, int A, typename Tx>
    inline Matrix<T,A>& operator*=(Matrix<T,A>& m, const Tx& x)
    { m.view() *= x; return m; }

    template <typename T, int A, typename Tx>
    inline Matrix<T,A>& operator/=(Matrix<T,A>& m, const Tx& x)
    { m.view() /= x; return m; }

    template <typename T, int A, typename Tx>
    inline Matrix<T,A>& operator%=(Matrix<T,A>& m, const Tx& x)
    { m.view() %= x; return m; }

    //
    // Scalar * Matrix
    //

    template <typename T, typename Tm>
    class ProdXM : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdXM(const T _x, const GenMatrix<Tm>& _m) : x(_x), m(_m) {}
        inline ptrdiff_t colsize() const { return m.colsize(); }
        inline ptrdiff_t rowsize() const { return m.rowsize(); }
        inline T getX() const { return x; }
        inline const GenMatrix<Tm>& getM() const { return m; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            MultXM(x,m0=m);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            MultXM(x,m0=m);
        }
    private:
        const T x;
        const GenMatrix<Tm>& m;
    };

    // m*=x
    template <typename T>
    inline MatrixView<T> operator*=(MatrixView<T> m, T x)
    { MultXM(x,m); return m; }

    template <typename T>
    inline MatrixView<CT> operator*=(MatrixView<CT> m, T x)
    { MultXM(T(x),m); return m; }

    template <typename T>
    inline MatrixView<CT> operator*=(MatrixView<CT> m, CCT x)
    { MultXM(CT(x),m); return m; }

    template <typename T>
    inline MatrixView<CT> operator*=(MatrixView<CT> m, VCT x)
    { MultXM(CT(x),m); return m; }

    // m/=x
    template <typename T>
    inline MatrixView<T> operator/=(MatrixView<T> m, T x)
    { MultXM(TMV_InverseOf(x),m); return m; }

    template <typename T>
    inline MatrixView<CT> operator/=(MatrixView<CT> m, T x)
    { MultXM(TMV_InverseOf(x),m); return m; }

    template <typename T>
    inline MatrixView<CT> operator/=(MatrixView<CT> m, CCT x)
    { MultXM(TMV_InverseOf(CT(x)),m); return m; }

    template <typename T>
    inline MatrixView<CT> operator/=(MatrixView<CT> m, VCT x)
    { MultXM(TMV_InverseOf(CT(x)),m); return m; }

#define GENMATRIX GenMatrix
#define PRODXM ProdXM
#include "tmv/TMV_AuxProdXM.h"
    // Defines things like -m, x*m, m*x, x*(x*m), etc.
#undef GENMATRIX
#undef PRODXM

    //
    // Matrix + Scalar
    //

    template <typename T, typename Tm>
    class SumMX : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumMX(T _x1, const GenMatrix<Tm>& _m, T _x2) :
            x1(_x1), m(_m), x2(_x2)
        { TMVAssert(m.isSquare()); }
        inline ptrdiff_t colsize() const { return m.colsize(); }
        inline ptrdiff_t rowsize() const { return m.rowsize(); }
        inline T getX1() const { return x1; }
        inline const GenMatrix<Tm>& getM() const { return m; }
        inline T getX2() const { return x2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            MultXM(x1,m0=m);
            m0.diag().addToAll(TMV_REAL(x2));
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            MultXM(x1,m0=m);
            m0.diag().addToAll(x2);
        }
    private:
        const T x1;
        const GenMatrix<Tm>& m;
        const T x2;
    };

    // m+=x
    template <typename T>
    inline MatrixView<T> operator+=(MatrixView<T> m, T x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(x);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(MatrixView<CT> m, T x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(CT(x));
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(MatrixView<CT> m, CCT x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(CT(x));
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(MatrixView<CT> m, VCT x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(CT(x));
        return m;
    }

    // m-=x
    template <typename T>
    inline MatrixView<T> operator-=(MatrixView<T> m, T x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(-x);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(MatrixView<CT> m, T x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(CT(-x));
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(MatrixView<CT> m, CCT x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(-CT(x));
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(MatrixView<CT> m, VCT x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(-CT(x));
        return m;
    }

#define GENMATRIX GenMatrix
#define PRODXM ProdXM
#define SUMMX SumMX
#include "tmv/TMV_AuxSumMX.h"
    // Defines things like m+x, x+m, x-m, m-x, x+x*m, x*(x+m), etc.
#undef SUMMX
#undef GENMATRIX
#undef PRODXM

    //
    // Vector ^ Vector (OuterProduct)
    //

    template <typename T, typename T1, typename T2>
    class OProdVV : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline OProdVV(
            const T _x, const GenVector<T1>& _v1, const GenVector<T2>& _v2) :
            x(_x), v1(_v1), v2(_v2) {}
        inline ptrdiff_t colsize() const { return v1.size(); }
        inline ptrdiff_t rowsize() const { return v2.size(); }
        inline T getX() const { return x; }
        inline const GenVector<T1>& getV1() const { return v1; }
        inline const GenVector<T2>& getV2() const { return v2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            Rank1Update<false>(x, v1, v2, m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            Rank1Update<false>(x, v1, v2, m0);
        }

    private:
        T x;
        const GenVector<T1>& v1;
        const GenVector<T2>& v2;
    };

    // m+=(x*v^v)
    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m0, const OProdVV<T,T1,T2>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(opvv.getX(), opvv.getV1(), opvv.getV2(), m0);
        return m0;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m0, const OProdVV<T,T,T>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(opvv.getX(), opvv.getV1(), opvv.getV2(), m0);
        return m0;
    }

    // m-=(x*v^v)
    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m0, const OProdVV<T,T1,T2>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(-opvv.getX(), opvv.getV1(), opvv.getV2(), m0);
        return m0;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m0, const OProdVV<T,T,T>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(-opvv.getX(), opvv.getV1(), opvv.getV2(), m0);
        return m0;
    }

#define PRODMM OProdVV
#define GENMATRIX1 GenVector
#define GENMATRIX2 GenVector
#define PRODXM1 ProdXV
#define PRODXM2 ProdXV
#define OP operator^
#define GETM1 .getV()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .getV1()
#define GETM2 .getV2()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // Matrix + Matrix
    //

    template <typename T, typename T1, typename T2>
    class SumMM : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumMM(
            const T _x1, const GenMatrix<T1>& _m1,
            const T _x2, const GenMatrix<T2>& _m2) :
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
        inline const GenMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            AddMM(x1,m1,x2,m2,m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            AddMM(x1,m1,x2,m2,m0);
        }
    private:
        const T x1;
        const GenMatrix<T1>& m1;
        const T x2;
        const GenMatrix<T2>& m2;
    };

    template <typename T>
    inline MatrixView<T> operator+=(
        MatrixView<T> m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        AddMM(T(1),m2,m1);
        return m1;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        AddMM(T(1),m2,m1);
        return m1;
    }

    template <typename T>
    inline MatrixView<T> operator-=(
        MatrixView<T> m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        AddMM(T(-1),m2,m1);
        return m1;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        AddMM(T(-1),m2,m1);
        return m1;
    }

    template <typename T, typename T2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ProdXM<T,T2>& pxm)
    {
        TMVAssert(m.colsize() == pxm.colsize());
        TMVAssert(m.rowsize() == pxm.rowsize());
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ProdXM<T,T>& pxm)
    {
        TMVAssert(m.colsize() == pxm.colsize());
        TMVAssert(m.rowsize() == pxm.rowsize());
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T, typename T2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ProdXM<T,T2>& pxm)
    {
        TMVAssert(m.colsize() == pxm.colsize());
        TMVAssert(m.rowsize() == pxm.rowsize());
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ProdXM<T,T>& pxm)
    {
        TMVAssert(m.colsize() == pxm.colsize());
        TMVAssert(m.rowsize() == pxm.rowsize());
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXM
#define SUMMM SumMM
#include "tmv/TMV_AuxSumMM.h"
    // Defines things like m+m, m-m, (x*m)-(x*m), etc.
#include "tmv/TMV_AuxSumMMa.h"
    // Defines things like -(m+m), x*(m+m), (m+m)/x, etc.
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM

    //
    // Matrix * Vector
    //

    template <typename T, typename T1, typename T2>
    class ProdMV : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdMV(
            const T _x, const GenMatrix<T1>& _m, const GenVector<T2>& _v) :
            x(_x), m(_m), v(_v)
        { TMVAssert(v.size()==m.rowsize()); }
        inline ptrdiff_t size() const { return m.colsize(); }
        inline T getX() const { return x; }
        inline const GenMatrix<T1>& getM() const { return m; }
        inline const GenVector<T2>& getV() const { return v; }
        inline void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            MultMV<false>(x,m,v,v0);
        }
        inline void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == size());
            MultMV<false>(x,m,v,v0);
        }
    private:
        const T x;
        const GenMatrix<T1>& m;
        const GenVector<T2>& v;
    };

    template <typename T, typename T1, typename T2>
    inline VectorView<T> operator+=(
        VectorView<T> v, const ProdMV<T,T1,T2>& pmv)
    {
        TMVAssert(v.size() == pmv.size());
        MultMV<true>(pmv.getX(),pmv.getM(),pmv.getV(),v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator+=(
        VectorView<CT> v, const ProdMV<T,T,T>& pmv)
    {
        TMVAssert(v.size() == pmv.size());
        MultMV<true>(pmv.getX(),pmv.getM(),pmv.getV(),v);
        return v;
    }

    template <typename T, typename T1, typename T2>
    inline VectorView<T> operator-=(
        VectorView<T> v, const ProdMV<T,T1,T2>& pmv)
    {
        TMVAssert(v.size() == pmv.size());
        MultMV<true>(-pmv.getX(), pmv.getM(), pmv.getV(), v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator-=(
        VectorView<CT> v, const ProdMV<T,T,T>& pmv)
    {
        TMVAssert(v.size() == pmv.size());
        MultMV<true>(-pmv.getX(),pmv.getM(),pmv.getV(),v);
        return v;
    }

    template <typename T, typename T1, typename T2>
    class ProdVM : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdVM(
            const T _x, const GenVector<T1>& _v, const GenMatrix<T2>& _m) :
            x(_x), v(_v), m(_m)
        { TMVAssert(v.size()==m.colsize()); }
        inline ptrdiff_t size() const { return m.rowsize(); }
        inline T getX() const { return x; }
        inline const GenVector<T1>& getV() const { return v; }
        inline const GenMatrix<T2>& getM() const { return m; }
        inline void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            MultMV<false>(x,m.transpose(),v,v0);
        }
        inline void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == size());
            MultMV<false>(x,m.transpose(),v,v0);
        }
    private:
        const T x;
        const GenVector<T1>& v;
        const GenMatrix<T2>& m;
    };

    template <typename T>
    inline VectorView<T> operator*=(
        VectorView<T> v, const GenMatrix<T>& m)
    {
        TMVAssert(v.size() == m.colsize());
        TMVAssert(v.size() == m.rowsize());
        MultMV<false>(T(1),m.transpose(),v,v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator*=(
        VectorView<CT> v, const GenMatrix<T>& m)
    {
        TMVAssert(v.size() == m.colsize());
        TMVAssert(v.size() == m.rowsize());
        MultMV<false>(T(1),m.transpose(),v,v);
        return v;
    }

    template <typename T, typename T2>
    inline VectorView<T> operator*=(
        VectorView<T> v, const ProdXM<T,T2>& pxm)
    {
        TMVAssert(v.size() == pxm.colsize());
        TMVAssert(v.size() == pxm.rowsize());
        MultMV<false>(pxm.getX(),pxm.getM().transpose(),v,v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator*=(
        VectorView<CT> v, const ProdXM<T,T>& pxm)
    {
        TMVAssert(v.size() == pxm.colsize());
        TMVAssert(v.size() == pxm.rowsize());
        MultMV<false>(pxm.getX(),pxm.getM().transpose(),v,v);
        return v;
    }

    template <typename T, typename T1, typename T2>
    inline VectorView<T> operator+=(
        VectorView<T> v, const ProdVM<T,T1,T2>& pvm)
    {
        TMVAssert(v.size() == pvm.size());
        MultMV<true>(pvm.getX(),pvm.getM().transpose(),pvm.getV(),v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator+=(
        VectorView<CT> v, const ProdVM<T,T,T>& pvm)
    {
        TMVAssert(v.size() == pvm.size());
        MultMV<true>(pvm.getX(),pvm.getM().transpose(),pvm.getV(),v);
        return v;
    }

    template <typename T, typename T1, typename T2>
    inline VectorView<T> operator-=(
        VectorView<T> v, const ProdVM<T,T1,T2>& pvm)
    {
        TMVAssert(v.size() == pvm.size());
        MultMV<true>(-pvm.getX(), pvm.getM().transpose(), pvm.getV(), v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator-=(
        VectorView<CT> v, const ProdVM<T,T,T>& pvm)
    {
        TMVAssert(v.size() == pvm.size());
        MultMV<true>(-pvm.getX(),pvm.getM().transpose(),pvm.getV(),v);
        return v;
    }

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenVector
#define PRODXM1 ProdXM
#define PRODXM2 ProdXV
#define PRODMM ProdMV
#define GETM1 .getM()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMM.h"
    // Defines things like m*v, (x*m)*(x*v), etc.
#define GETM1 .getM()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMMa.h"
    // Defines things like -(m*v), x*(m*v), (m*v)/x, etc.
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef PRODMM
#define GENMATRIX1 GenVector
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXV
#define PRODXM2 ProdXM
#define PRODMM ProdVM
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef PRODMM

    //
    // Matrix * Matrix
    //

    template <typename T, typename T1, typename T2>
    class ProdMM : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdMM(
            const T _x, const GenMatrix<T1>& _m1, const GenMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.rowsize() == m2.colsize()) ; }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m2.rowsize(); }
        inline T getX() const { return x; }
        inline const GenMatrix<T1>& getM1() const { return m1; }
        inline const GenMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            MultMM<false>(x, m1, m2, m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            MultMM<false>(x, m1, m2, m0);
        }

    private:
        const T x;
        const GenMatrix<T1>& m1;
        const GenMatrix<T2>& m2;
    };

    template <typename T>
    inline MatrixView<T> operator*=(
        MatrixView<T> m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.colsize()==m2.rowsize());
        TMVAssert(m1.rowsize()==m2.rowsize());
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T>
    inline MatrixView<CT> operator*=(
        MatrixView<CT> m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.colsize()==m2.rowsize());
        TMVAssert(m1.rowsize()==m2.rowsize());
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T, typename T2>
    inline MatrixView<T> operator*=(
        MatrixView<T> m1, const ProdXM<T,T2>& pxm)
    {
        TMVAssert(pxm.colsize()==pxm.rowsize());
        TMVAssert(m1.rowsize()==pxm.rowsize());
        MultMM<false>(pxm.getX(),m1,pxm.getM(),m1);
        return m1;
    }

    template <typename T>
    inline MatrixView<CT> operator*=(
        MatrixView<CT> m1, const ProdXM<T,T>& pxm)
    {
        TMVAssert(pxm.colsize()==pxm.rowsize());
        TMVAssert(m1.rowsize()==pxm.rowsize());
        MultMM<false>(pxm.getX(),m1,pxm.getM(),m1);
        return m1;
    }

    template <typename T, typename T2, typename T3>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ProdMM<T,T2,T3>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ProdMM<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T2, typename T3>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ProdMM<T,T2,T3>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ProdMM<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXM
#define PRODMM ProdMM
#include "tmv/TMV_AuxProdMM.h"
    // Defines things like m*m, m*(x*m), etc.
#include "tmv/TMV_AuxProdMMa.h"
    // Defines things like -(m*m), x*(m*m), etc.
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef PRODMM


    //
    // Element Product Matrix * Matrix
    //

    template <typename T, typename T1, typename T2>
    class ElemProdMM : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ElemProdMM(
            const T _x, const GenMatrix<T1>& _m1, const GenMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        {
            TMVAssert(m1.colsize() == m2.colsize());
            TMVAssert(m1.rowsize() == m2.rowsize());
        }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline T getX() const { return x; }
        inline const GenMatrix<T1>& getM1() const { return m1; }
        inline const GenMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            ElemMultMM<false>(x, m1, m2, m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            ElemMultMM<false>(x, m1, m2, m0);
        }

    private:
        const T x;
        const GenMatrix<T1>& m1;
        const GenMatrix<T2>& m2;
    };

    template <typename T, typename T2, typename T3>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ElemProdMM<T,T2,T3>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        ElemMultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ElemProdMM<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        ElemMultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T2, typename T3>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ElemProdMM<T,T2,T3>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        ElemMultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ElemProdMM<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        ElemMultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

#define PRODMM ElemProdMM
#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXM
#define OP ElemProd
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // Scalar / Matrix
    //

    template <typename T, typename Tm>
    class QuotXM : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotXM(const T _x, const GenMatrix<Tm>& _m) : x(_x), m(_m) {}
        inline ptrdiff_t colsize() const { return m.rowsize(); }
        inline ptrdiff_t rowsize() const { return m.colsize(); }
        inline T getX() const { return x; }
        inline const GenMatrix<Tm>& getM() const { return m; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            m.makeInverse(m0);
            MultXM(x,m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            m.makeInverse(m0);
            MultXM(x,m0);
        }
    private:
        const T x;
        const GenMatrix<Tm>& m;
    };

#define GENMATRIX GenMatrix
#define PRODXM ProdXM
#define QUOTXM QuotXM
#include "tmv/TMV_AuxQuotXM.h"
    // Defines x/m, x%m
#include "tmv/TMV_AuxQuotXMa.h"
    // Defines x*(x/m), -(x/m), etc.
#undef GENMATRIX
#undef PRODXM
#undef QUOTXM


    //
    // Vector / % Matrix
    // v/m is the solution (x) of mx = v
    // ie. / is really division from the left: x = m^-1 v
    // Use % if you want division from the right (v m^-1)
    //

    template <typename T, typename T1, typename T2>
    class QuotVM : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotVM(
            const T _x, const GenVector<T1>& _v, const GenMatrix<T2>& _m) :
            x(_x), v(_v), m(_m)
        { TMVAssert(v.size()==m.colsize()); }
        inline ptrdiff_t size() const { return m.rowsize(); }
        inline T getX() const { return x; }
        inline const GenVector<T1>& getV() const { return v; }
        inline const GenMatrix<T2>& getM() const { return m; }
        inline void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            m.LDiv(v,v0);
            MultXV(x,v0);
        }
        inline void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == size());
            m.LDiv(v,v0);
            MultXV(x,v0);
        }
    private:
        const T x;
        const GenVector<T1>& v;
        const GenMatrix<T2>& m;
    };

    template <typename T, typename T1, typename T2>
    class RQuotVM : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline RQuotVM(
            const T _x, const GenVector<T1>& _v, const GenMatrix<T2>& _m) :
            x(_x), v(_v), m(_m)
        { TMVAssert(v.size()==m.rowsize()); }
        inline ptrdiff_t size() const { return m.colsize(); }
        inline T getX() const { return x; }
        inline const GenVector<T1>& getV() const { return v; }
        inline const GenMatrix<T2>& getM() const { return m; }
        inline void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            m.RDiv(v,v0);
            MultXV(x,v0);
        }
        inline void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == size());
            m.RDiv(v,v0);
            MultXV(x,v0);
        }
    private:
        const T x;
        const GenVector<T1>& v;
        const GenMatrix<T2>& m;
    };

    template <typename T>
    inline VectorView<T> operator/=(
        VectorView<T> v, const GenMatrix<T>& m)
    {
        TMVAssert(m.isSquare());
        TMVAssert(m.rowsize() == v.size());
        m.LDivEq(v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator/=(
        VectorView<CT> v, const GenMatrix<T>& m)
    {
        TMVAssert(m.isSquare());
        TMVAssert(m.rowsize() == v.size());
        m.LDivEq(v);
        return v;
    }

    template <typename T>
    inline VectorView<T> operator%=(
        VectorView<T> v, const GenMatrix<T>& m)
    {
        TMVAssert(m.isSquare());
        TMVAssert(m.rowsize() == v.size());
        m.RDivEq(v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator%=(
        VectorView<CT> v, const GenMatrix<T>& m)
    {
        TMVAssert(m.isSquare());
        TMVAssert(m.rowsize() == v.size());
        m.RDivEq(v);
        return v;
    }

    template <typename T, typename Tm>
    inline VectorView<T> operator*=(
        VectorView<T> v, const QuotXM<T,Tm>& qxm)
    {
        TMVAssert(qxm.getM().isSquare());
        TMVAssert(qxm.getM().rowsize() == v.size());
        qxm.getM().RDivEq(v);
        MultXV(qxm.getX(),v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator*=(
        VectorView<CT> v, const QuotXM<T,T>& qxm)
    {
        TMVAssert(qxm.getM().isSquare());
        TMVAssert(qxm.getM().rowsize() == v.size());
        qxm.getM().RDivEq(v);
        MultXV(qxm.getX(),v);
        return v;
    }

#define GENMATRIX1 GenVector
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXV
#define PRODXM2 ProdXM
#define QUOTMM QuotVM
#define RQUOTMM RQuotVM
#define QUOTXM QuotXM
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxQuotMM.h"
    // Defines things like v/m, v%m, (x*v)/m, v/(x*m), etc.
#define PRODMM QuotVM
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
    // Defines things like (v/m)*x, -(v/m), etc.
#undef PRODMM
#define PRODMM RQuotVM
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTMM
#undef RQUOTMM
#undef QUOTXM


    //
    // Matrix / % Matrix
    //

    template <typename T, typename T1, typename T2>
    class QuotMM : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotMM(
            const T _x, const GenMatrix<T1>& _m1, const GenMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.colsize() == m2.colsize() ); }
        inline ptrdiff_t colsize() const { return m2.rowsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline T getX() const { return x; }
        inline const GenMatrix<T1>& getM1() const { return m1; }
        inline const GenMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            m2.LDiv(m1,m0);
            MultXM(x,m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            m2.LDiv(m1,m0);
            MultXM(x,m0);
        }
    protected:
        const T x;
        const GenMatrix<T1>& m1;
        const GenMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class TransientQuotMM : public QuotMM<T,T1,T2>
    {
    public :
        inline TransientQuotMM(
            const T x, auto_ptr<Matrix<T1,ColMajor> > m1,
            const GenMatrix<T2>& m2) :
            QuotMM<T,T1,T2>(x,*m1,m2), m1p(m1.release()) {}
        inline TransientQuotMM(const TransientQuotMM<T,T1,T2>& rhs) :
            QuotMM<T,T1,T2>(rhs), m1p(rhs.m1p.release()) {}
        inline ~TransientQuotMM() {}
        inline auto_ptr<Matrix<T1,ColMajor> > getP() const { return m1p; }

    private :
        mutable auto_ptr<Matrix<T1,ColMajor> > m1p;
    };

    template <typename T, typename T1, typename T2>
    class RQuotMM : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline RQuotMM(
            const T _x, const GenMatrix<T1>& _m1, const GenMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.rowsize() == m2.rowsize() ); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m2.colsize(); }
        inline T getX() const { return x; }
        inline const GenMatrix<T1>& getM1() const { return m1; }
        inline const GenMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            m2.RDiv(m1,m0);
            MultXM(x,m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            m2.RDiv(m1,m0);
            MultXM(x,m0);
        }
    protected:
        const T x;
        const GenMatrix<T1>& m1;
        const GenMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class TransientRQuotMM : public RQuotMM<T,T1,T2>
    {
    public :
        inline TransientRQuotMM(
            const T x, auto_ptr<Matrix<T1,RowMajor> > m1,
            const GenMatrix<T2>& m2) :
            RQuotMM<T,T1,T2>(x,*m1,m2), m1p(m1.release()) {}
        inline TransientRQuotMM(const TransientRQuotMM<T,T1,T2>& rhs) :
            RQuotMM<T,T1,T2>(rhs), m1p(rhs.m1p.release()) {}
        inline ~TransientRQuotMM() {}
        inline auto_ptr<Matrix<T1,RowMajor> > getP() const { return m1p; }

    private :
        mutable auto_ptr<Matrix<T1,RowMajor> > m1p;
    };

    template <typename T>
    inline MatrixView<T> operator/=(
        MatrixView<T> m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m2.rowsize());
        TMVAssert(m1.colsize() == m2.rowsize());
        m2.LDivEq(m1);
        return m1;
    }

    template <typename T>
    inline MatrixView<CT> operator/=(
        MatrixView<CT> m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m2.rowsize());
        TMVAssert(m1.colsize() == m2.rowsize());
        m2.LDivEq(m1);
        return m1;
    }

    template <typename T>
    inline MatrixView<T> operator%=(
        MatrixView<T> m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        m2.RDivEq(m1);
        return m1;
    }

    template <typename T>
    inline MatrixView<CT> operator%=(
        MatrixView<CT> m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        m2.RDivEq(m1);
        return m1;
    }

    template <typename T, typename Tm>
    inline MatrixView<T> operator*=(
        MatrixView<T> m1, const QuotXM<T,Tm>& qxm)
    {
        TMVAssert(qxm.getM().isSquare());
        TMVAssert(m1.rowsize() == qxm.getM().rowsize());
        qxm.getM().RDivEq(m1);
        MultXM(qxm.getX(),m1);
        return m1;
    }

    template <typename T>
    inline MatrixView<CT> operator*=(
        MatrixView<CT> m1, const QuotXM<T,T>& qxm)
    {
        TMVAssert(qxm.getM().isSquare());
        TMVAssert(m1.rowsize() == qxm.getM().rowsize());
        qxm.getM().RDivEq(m1);
        MultXM(qxm.getX(),m1);
        return m1;
    }

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXM
#define QUOTXM QuotXM
#define QUOTMM QuotMM
#define RQUOTMM RQuotMM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM
#include "tmv/TMV_AuxQuotMM.h"
    // Defines things like m/m, m%m, (x*m)/m, m/(x*m), etc.
#define PRODMM QuotMM
#include "tmv/TMV_AuxProdMMa.h"
    // Defines things like (m/m)*x, -(m/m), etc.
#undef PRODMM
#define PRODMM RQuotMM
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM TransientQuotMM
#define GETM1 .getP()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM TransientRQuotMM
#define GETM1 .getP()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif
