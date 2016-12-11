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


#ifndef TMV_BandMatrixArith_H
#define TMV_BandMatrixArith_H

#include "tmv/TMV_BaseBandMatrix.h"
#include "tmv/TMV_BandMatrixArithFunc.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

    template <typename T, typename Tv>
    class ProdXV;
    template <typename T, typename Tv>
    class ProdXM;

    template <typename T, int A, typename Tx>
    inline BandMatrix<T,A>& operator+=(BandMatrix<T,A>& m, const Tx& x)
    { m.view() += x; return m; }

    template <typename T, int A, typename Tx>
    inline BandMatrix<T,A>& operator-=(BandMatrix<T,A>& m, const Tx& x)
    { m.view() -= x; return m; }

    template <typename T, int A, typename Tx>
    inline BandMatrix<T,A>& operator*=(BandMatrix<T,A>& m, const Tx& x)
    { m.view() *= x; return m; }

    template <typename T, int A, typename Tx>
    inline BandMatrix<T,A>& operator/=(BandMatrix<T,A>& m, const Tx& x)
    { m.view() /= x; return m; }

    template <typename T, int A, typename Tx>
    inline BandMatrix<T,A>& operator%=(BandMatrix<T,A>& m, const Tx& x)
    { m.view() %= x; return m; }

    template <typename T, typename T2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m1, const GenBandMatrix<T2>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        BandMatrixView<T>(m1,m2.nlo(),m2.nhi()) += m2;
        return m1;
    }

    template <typename T, typename T2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m1, const GenBandMatrix<T2>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        BandMatrixView<T>(m1,m2.nlo(),m2.nhi()) -= m2;
        return m1;
    }

    //
    // Scalar * BandMatrix
    //

    template <typename T, typename Tm>
    class ProdXB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdXB(const T _x, const GenBandMatrix<Tm>& _m) :
            x(_x), m(_m) {}
        inline ptrdiff_t colsize() const { return m.colsize(); }
        inline ptrdiff_t rowsize() const { return m.rowsize(); }
        inline ptrdiff_t nlo() const { return m.nlo(); }
        inline ptrdiff_t nhi() const { return m.nhi(); }
        inline T getX() const { return x; }
        inline const GenBandMatrix<Tm>& getM() const { return m; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            MultXM(x,m0=m);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            MultXM(x,m0=m);
        }
    private:
        const T x;
        const GenBandMatrix<Tm>& m;
    };

    // m*=x
    template <typename T>
    inline BandMatrixView<T> operator*=(
        BandMatrixView<T> m, T x)
    { MultXM(x,m); return m; }

    template <typename T>
    inline BandMatrixView<CT> operator*=(
        BandMatrixView<CT> m, T x)
    { MultXM(x,m); return m; }

    template <typename T>
    inline BandMatrixView<CT> operator*=(
        BandMatrixView<CT> m, CCT x)
    { MultXM(CT(x),m); return m; }

    template <typename T>
    inline BandMatrixView<CT> operator*=(
        BandMatrixView<CT> m, VCT x)
    { MultXM(CT(x),m); return m; }

    // m/=x
    template <typename T>
    inline BandMatrixView<T> operator/=(
        BandMatrixView<T> m, T x)
    { MultXM(TMV_InverseOf(x),m); return m; }

    template <typename T>
    inline BandMatrixView<CT> operator/=(
        BandMatrixView<CT> m, T x)
    { MultXM(TMV_InverseOf(x),m); return m; }

    template <typename T>
    inline BandMatrixView<CT> operator/=(
        BandMatrixView<CT> m, CCT x)
    { MultXM(TMV_InverseOf(CT(x)),m); return m; }

    template <typename T>
    inline BandMatrixView<CT> operator/=(
        BandMatrixView<CT> m, VCT x)
    { MultXM(TMV_InverseOf(CT(x)),m); return m; }

#define GENMATRIX GenBandMatrix
#define PRODXM ProdXB
#include "tmv/TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM


    //
    // BandMatrix + Scalar
    //

    template <typename T, typename Tm>
    class SumBX : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumBX(T _x1, const GenBandMatrix<Tm>& _m, T _x2) :
            x1(_x1), m(_m), x2(_x2)
        { TMVAssert(m.isSquare()); }
        inline ptrdiff_t colsize() const { return m.colsize(); }
        inline ptrdiff_t rowsize() const { return m.rowsize(); }
        inline ptrdiff_t nlo() const { return m.nlo(); }
        inline ptrdiff_t nhi() const { return m.nhi(); }
        inline T getX1() const { return x1; }
        inline const GenBandMatrix<Tm>& getM() const { return m; }
        inline T getX2() const { return x2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            MultXM(x1,m0=m);
            m0.diag().addToAll(TMV_REAL(x2));
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            MultXM(x1,m0=m);
            m0.diag().addToAll(complex_type(x2));
        }
    private:
        const T x1;
        const GenBandMatrix<Tm>& m;
        const T x2;
    };

    // m+=x
    template <typename T>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, T x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(x); return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, T x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(CT(x)); return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, CCT x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(CT(x)); return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, VCT x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(CT(x)); return m;
    }

    // m-=x
    template <typename T>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, T x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(-x); return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, T x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(CT(-x)); return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, CCT x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(-CT(x)); return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, VCT x)
    {
        TMVAssert(m.isSquare());
        m.diag().addToAll(-CT(x)); return m;
    }

#define SUMMX SumBX
#define GENMATRIX GenBandMatrix
#define PRODXM ProdXB
#include "tmv/TMV_AuxSumMX.h"
#undef SUMMX
#undef GENMATRIX
#undef PRODXM


    //
    // BandMatrix + BandMatrix
    //

    template <typename T, typename T1, typename T2>
    class SumBB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumBB(
            T _x1, const GenBandMatrix<T1>& _m1,
            T _x2, const GenBandMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        {
            TMVAssert(m1.rowsize() == m2.rowsize());
            TMVAssert(m1.colsize() == m2.colsize());
        }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline ptrdiff_t nlo() const { return TMV_MAX(m1.nlo(),m2.nlo()); }
        inline ptrdiff_t nhi() const { return TMV_MAX(m1.nhi(),m2.nhi()); }
        inline T getX1() const { return x1; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            AddMM(x1,m1,x2,m2,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            AddMM(x1,m1,x2,m2,m0);
        }
    private:
        T x1;
        const GenBandMatrix<T1>& m1;
        T x2;
        const GenBandMatrix<T2>& m2;
    };

    template <typename T>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() >= m2.nlo());
        TMVAssert(m1.nhi() >= m2.nhi());
        AddMM(T(1),m2,m1); return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() >= m2.nlo());
        TMVAssert(m1.nhi() >= m2.nhi());
        AddMM(T(1),m2,m1); return m1;
    }

    template <typename T>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() >= m2.nlo());
        TMVAssert(m1.nhi() >= m2.nhi());
        AddMM(T(-1),m2,m1); return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() >= m2.nlo());
        TMVAssert(m1.nhi() >= m2.nhi());
        AddMM(T(-1),m2,m1); return m1;
    }

    template <typename T, typename T2>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, const ProdXB<T,T2>& pxm)
    {
        TMVAssert(m.colsize() == pxm.colsize());
        TMVAssert(m.rowsize() == pxm.rowsize());
        TMVAssert(m.nlo() >= pxm.nlo());
        TMVAssert(m.nhi() >= pxm.nhi());
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, const ProdXB<T,T>& pxm)
    {
        TMVAssert(m.colsize() == pxm.colsize());
        TMVAssert(m.rowsize() == pxm.rowsize());
        TMVAssert(m.nlo() >= pxm.nlo());
        TMVAssert(m.nhi() >= pxm.nhi());
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T, typename T2>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, const ProdXB<T,T2>& pxm)
    {
        TMVAssert(m.colsize() == pxm.colsize());
        TMVAssert(m.rowsize() == pxm.rowsize());
        TMVAssert(m.nlo() >= pxm.nlo());
        TMVAssert(m.nhi() >= pxm.nhi());
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, const ProdXB<T,T>& pxm)
    {
        TMVAssert(m.colsize() == pxm.colsize());
        TMVAssert(m.rowsize() == pxm.rowsize());
        TMVAssert(m.nlo() >= pxm.nlo());
        TMVAssert(m.nhi() >= pxm.nhi());
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

#define SUMMM SumBB
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXB
#include "tmv/TMV_AuxSumMM.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // BandMatrix * BandMatrix
    //

    template <typename T, typename T1, typename T2>
    class ProdBB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdBB(
            T _x, const GenBandMatrix<T1>& _m1,
            const GenBandMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.rowsize() == m2.colsize()); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m2.rowsize(); }
        inline ptrdiff_t nlo() const
        { return TMV_MIN(colsize()-1,m1.nlo()+m2.nlo()); }
        inline ptrdiff_t nhi() const
        { return TMV_MIN(rowsize()-1,m1.nhi()+m2.nhi()); }
        inline T getX() const { return x; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline const GenBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            MultMM<false>(x,m1,m2,m0);
        }
    private:
        T x;
        const GenBandMatrix<T1>& m1;
        const GenBandMatrix<T2>& m2;
    };

    template <typename T>
    inline BandMatrixView<T> operator*=(
        BandMatrixView<T> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.nlo() == 0 || m1.nlo() == m1.colsize()-1);
        TMVAssert(m2.nhi() == 0 || m1.nhi() == m1.rowsize()-1);
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator*=(
        BandMatrixView<CT> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.nlo() == 0 || m1.nlo() == m1.colsize()-1);
        TMVAssert(m2.nhi() == 0 || m1.nhi() == m1.rowsize()-1);
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, const ProdBB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, const ProdBB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, const ProdBB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, const ProdBB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2, typename T3>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ProdBB<T1,T2,T3>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixView<T>(m,pmm.nlo(),pmm.nhi()) += pmm;
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ProdBB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixView<CT>(m,pmm.nlo(),pmm.nhi()) += pmm;
        return m;
    }

    template <typename T, typename T1, typename T2, typename T3>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ProdBB<T1,T2,T3>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixView<T>(m,pmm.nlo(),pmm.nhi()) -= pmm;
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ProdBB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixView<CT>(m,pmm.nlo(),pmm.nhi()) -= pmm;
        return m;
    }

#define PRODMM ProdBB
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXB
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // Element Product BandMatrix * BandMatrix
    //

    template <typename T, typename T1, typename T2>
    class ElemProdBB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ElemProdBB(
            T _x, const GenBandMatrix<T1>& _m1, const GenBandMatrix<T2>& _m2) :
            x(_x),m1(_m1),m2(_m2)
        {
            TMVAssert(m1.rowsize() == m2.rowsize());
            TMVAssert(m1.colsize() == m2.colsize());
        }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline ptrdiff_t nlo() const { return TMV_MIN(m1.nlo(),m2.nlo()); }
        inline ptrdiff_t nhi() const { return TMV_MIN(m1.nhi(),m2.nhi()); }
        inline T getX() const { return x; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline const GenBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            ElemMultMM<false>(x,m1,m2,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            ElemMultMM<false>(x,m1,m2,m0);
        }
    private:
        T x;
        const GenBandMatrix<T1>& m1;
        const GenBandMatrix<T2>& m2;
    };

    template <typename T, typename T2, typename T3>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, const ElemProdBB<T,T2,T3>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        ElemMultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, const ElemProdBB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        ElemMultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T2, typename T3>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, const ElemProdBB<T,T2,T3>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        ElemMultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, const ElemProdBB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        ElemMultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

#define PRODMM ElemProdBB
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXB
#define OP ElemProd
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


#define GENMATRIX GenBandMatrix
#define PRODXM ProdXB
#define PRODMV ProdBV
#define PRODVM ProdVB
#define QUOTVM QuotVB
#define QUOTXM QuotXB
#define RQUOTVM RQuotVB
#include "tmv/TMV_AuxVecComposite.h"
#undef GENMATRIX
#undef PRODXM
#undef PRODMV
#undef PRODVM
#undef QUOTVM
#undef QUOTXM
#undef RQUOTVM

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenBandMatrix
#define GENMATRIX GenBandMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXB
#define PRODXM ProdXB
#define SUMMM SumMB
#define PRODMM ProdMB
#define QUOTXM QuotXB
#define QUOTMM QuotMB
#define RQUOTMM RQuotMB
#define TQUOTMM TransientQuotMB
#define TRQUOTMM TransientRQuotMB

#include "tmv/TMV_AuxMatComposite1.h"
#include "tmv/TMV_AuxSumMMb.h"
#include "tmv/TMV_AuxMatComposite3.h"


    // B/B -> Matrix.  Use TransientQuot
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenBandMatrix
#define PRODXM1 ProdXB
#include "tmv/TMV_AuxTQuotMM.h"


#undef GENMATRIX2
#undef PRODXM2
#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM
#define GENMATRIX2 GenMatrix
#define PRODXM2 ProdXM
#define SUMMMa SumMB
#define SUMMM SumBM
#define PRODMM ProdBM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM
#define QUOTXM QuotXM

#include "tmv/TMV_AuxMatComposite2.h"
#include "tmv/TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef GENMATRIX
#undef PRODXM1
#undef PRODXM2
#undef PRODXM
#undef SUMMM
#undef SUMMMa
#undef PRODMM
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#ifdef TMV_DiagMatrixArith_H
#include "tmv/TMV_DiagBandArith.h"
#endif

#ifdef TMV_TriMatrixArith_H
#include "tmv/TMV_TriBandArith.h"
#endif

#ifdef TMV_SymMatrixArith_H
#include "tmv/TMV_BandSymArith.h"
#endif

#ifdef TMV_SymBandMatrixArith_H
#include "tmv/TMV_BandSymBandArith.h"
#endif

#endif
