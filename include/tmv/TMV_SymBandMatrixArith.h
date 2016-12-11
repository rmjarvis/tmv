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


#ifndef TMV_SymBandMatrixArith_H
#define TMV_SymBandMatrixArith_H

#include "tmv/TMV_BaseSymBandMatrix.h"
#include "tmv/TMV_SymBandMatrixArithFunc.h"
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
    inline SymBandMatrix<T,A>& operator+=(SymBandMatrix<T,A>& m, const Tx& x)
    { m.view() += x; return m; }

    template <typename T, int A, typename Tx>
    inline SymBandMatrix<T,A>& operator-=(SymBandMatrix<T,A>& m, const Tx& x)
    { m.view() -= x; return m; }

    template <typename T, int A, typename Tx>
    inline SymBandMatrix<T,A>& operator*=(SymBandMatrix<T,A>& m, const Tx& x)
    { m.view() *= x; return m; }

    template <typename T, int A, typename Tx>
    inline SymBandMatrix<T,A>& operator/=(SymBandMatrix<T,A>& m, const Tx& x)
    { m.view() /= x; return m; }

    template <typename T, int A, typename Tx>
    inline SymBandMatrix<T,A>& operator%=(SymBandMatrix<T,A>& m, const Tx& x)
    { m.view() %= x; return m; }

    template <typename T, int A, typename Tx>
    inline HermBandMatrix<T,A>& operator+=(HermBandMatrix<T,A>& m, const Tx& x)
    { m.view() += x; return m; }

    template <typename T, int A, typename Tx>
    inline HermBandMatrix<T,A>& operator-=(HermBandMatrix<T,A>& m, const Tx& x)
    { m.view() -= x; return m; }

    template <typename T, int A, typename Tx>
    inline HermBandMatrix<T,A>& operator*=(HermBandMatrix<T,A>& m, const Tx& x)
    { m.view() *= x; return m; }

    template <typename T, int A, typename Tx>
    inline HermBandMatrix<T,A>& operator/=(HermBandMatrix<T,A>& m, const Tx& x)
    { m.view() /= x; return m; }

    template <typename T, int A, typename Tx>
    inline HermBandMatrix<T,A>& operator%=(HermBandMatrix<T,A>& m, const Tx& x)
    { m.view() %= x; return m; }


    //
    // Scalar * SymBandMatrix
    //

    template <typename T, typename Tm>
    class ProdXsB : public SymBandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdXsB(const T _x, const GenSymBandMatrix<Tm>& _m) :
            x(_x), m(_m) {}
        inline ptrdiff_t size() const { return m.size(); }
        inline SymType sym() const { return m.issym() ? Sym : Herm; }
        inline ptrdiff_t nlo() const { return m.nlo(); }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<Tm>& getM() const { return m; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nlo());
            m0.diagRange(0,m0.nhi()+1) = m.upperBand();
            if (m.nlo() > 0) m0.diagRange(-m0.nlo(),0) = m.lowerBandOff();
            else if (m0.nlo() > 0) m0.diagRange(-m0.nlo(),0).setZero();
            MultXM(x,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nlo());
            m0.diagRange(0,m0.nhi()+1) = m.upperBand();
            if (m.nlo() > 0) m0.diagRange(-m0.nlo(),0) = m.lowerBandOff();
            else if (m0.nlo() > 0) m0.diagRange(-m0.nlo(),0).setZero();
            MultXM(x,m0);
        }
        inline void assignTosB(SymBandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m.issym() || TMV_IMAG(x)==real_type(0) );
            TMVAssert(m0.size() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(isReal(Tm()) || m0.issym() == m.issym());
            MultXM(x,m0=m);
        }
        inline void assignTosB(SymBandMatrixView<complex_type> m0) const
        {
            TMVAssert(m.issym() || TMV_IMAG(x)==real_type(0) );
            TMVAssert(m0.size() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(isReal(Tm()) || m0.issym() == m.issym());
            MultXM(x,m0=m);
        }
    private:
        const T x;
        const GenSymBandMatrix<Tm>& m;
    };

    // m*=x
    template <typename T>
    inline SymBandMatrixView<T> operator*=(
        SymBandMatrixView<T> m, T x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        MultXM(x,m);
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator*=(
        SymBandMatrixView<CT> m, T x)
    {
        MultXM(x,m);
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator*=(
        SymBandMatrixView<CT> m, CCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        MultXM(CT(x),m);
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator*=(
        SymBandMatrixView<CT> m, VCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        MultXM(CT(x),m);
        return m;
    }

    // m/=x
    template <typename T>
    inline SymBandMatrixView<T> operator/=(
        SymBandMatrixView<T> m, T x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        MultXM(TMV_InverseOf(x),m);
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator/=(
        SymBandMatrixView<CT> m, T x)
    {
        MultXM(TMV_InverseOf(x),m);
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator/=(
        SymBandMatrixView<CT> m, CCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        MultXM(TMV_InverseOf(CT(x)),m);
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator/=(
        SymBandMatrixView<CT> m, VCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        MultXM(TMV_InverseOf(CT(x)),m);
        return m;
    }

#define GENMATRIX GenSymBandMatrix
#define PRODXM ProdXsB
#include "tmv/TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM

    //
    // SymBandMatrix + Scalar
    //

    template <typename T, typename Tm>
    class SumsBX : public SymBandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumsBX(T _x1, const GenSymBandMatrix<Tm>& _m, T _x2) :
            x1(_x1), m(_m), x2(_x2) {}
        inline ptrdiff_t size() const { return m.size(); }
        inline ptrdiff_t nlo() const { return m.nlo(); }
        inline SymType sym() const { return m.issym() ? Sym : Herm; }
        inline T getX1() const { return x1; }
        inline const GenSymBandMatrix<Tm>& getM() const { return m; }
        inline T getX2() const { return x2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nlo());
            MultXM(x1,m0=m);
            m0.diag().addToAll(TMV_REAL(x2));
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nlo());
            MultXM(x1,m0=m);
            m0.diag().addToAll(complex_type(x2));
        }
        inline void assignTosB(SymBandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(isReal(Tm()) || m0.issym() == m.issym());
            TMVAssert(m.issym() || TMV_IMAG(x1)==real_type(0));
            TMVAssert(m.issym() || TMV_IMAG(x2)==real_type(0));
            TMVAssert(m0.size() == size());
            TMVAssert(m0.nlo() >= nlo());
            MultXM(x1,m0=m);
            m0.diag().addToAll(TMV_REAL(x2));
        }
        inline void assignTosB(SymBandMatrixView<complex_type> m0) const
        {
            TMVAssert(isReal(Tm()) || m0.issym() == m.issym());
            TMVAssert(m.issym() || TMV_IMAG(x1)==real_type(0));
            TMVAssert(m.issym() || TMV_IMAG(x2)==real_type(0));
            TMVAssert(m0.size() == size());
            TMVAssert(m0.nlo() >= nlo());
            MultXM(x1,m0=m);
            m0.diag().addToAll(complex_type(x2));
        }
    private:
        const T x1;
        const GenSymBandMatrix<Tm>& m;
        const T x2;
    };

    // m+=x
    template <typename T>
    inline SymBandMatrixView<T> operator+=(
        SymBandMatrixView<T> m, T x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        m.diag().addToAll(x);
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator+=(
        SymBandMatrixView<CT> m, T x)
    {
        m.diag().addToAll(CT(x));
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator+=(
        SymBandMatrixView<CT> m, CCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == T(0));
        m.diag().addToAll(CT(x));
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator+=(
        SymBandMatrixView<CT> m, VCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == T(0));
        m.diag().addToAll(CT(x));
        return m;
    }

    // m-=x
    template <typename T>
    inline SymBandMatrixView<T> operator-=(
        SymBandMatrixView<T> m, T x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        m.diag().addToAll(-x);
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator-=(
        SymBandMatrixView<CT> m, T x)
    {
        m.diag().addToAll(CT(-x));
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator-=(
        SymBandMatrixView<CT> m, CCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == T(0));
        m.diag().addToAll(-CT(x));
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator-=(
        SymBandMatrixView<CT> m, VCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == T(0));
        m.diag().addToAll(-CT(x));
        return m;
    }

#define SUMMX SumsBX
#define GENMATRIX GenSymBandMatrix
#define PRODXM ProdXsB
#include "tmv/TMV_AuxSumMX.h"
#undef SUMMX
#undef GENMATRIX
#undef PRODXM

    //
    // SymBandMatrix + SymBandMatrix
    //

    template <typename T, typename T1, typename T2>
    class SumsBsB : public SymBandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumsBsB(
            T _x1, const GenSymBandMatrix<T1>& _m1,
            T _x2, const GenSymBandMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline ptrdiff_t nlo() const { return TMV_MAX(m1.nlo(),m2.nlo()); }
        inline SymType sym() const
        { return isReal(T1()) ? m2.sym() : m1.sym(); }
        inline T getX1() const { return x1; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenSymBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nlo());
            AddMM(x1,m1,x2,m2,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nlo());
            AddMM(x1,m1,x2,m2,m0);
        }
        inline void assignTosB(SymBandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(isReal(T1()) || isReal(T2()) || m1.sym() == m2.sym());
            TMVAssert(m0.size() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(isReal(T1()) || m0.issym() == m1.issym());
            TMVAssert(isReal(T2()) || m0.issym() == m2.issym());
            TMVAssert(m0.issym() || TMV_IMAG(x1) == real_type(0));
            TMVAssert(m0.issym() || TMV_IMAG(x2) == real_type(0));
            AddMM(x1,m1,x2,m2,m0);
        }
        inline void assignTosB(SymBandMatrixView<complex_type> m0) const
        {
            TMVAssert(isReal(T1()) || isReal(T2()) || m1.sym() == m2.sym());
            TMVAssert(m0.size() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(isReal(T1()) || m0.issym() == m1.issym());
            TMVAssert(isReal(T2()) || m0.issym() == m2.issym());
            TMVAssert(m0.issym() || TMV_IMAG(x1) == real_type(0));
            TMVAssert(m0.issym() || TMV_IMAG(x2) == real_type(0));
            AddMM(x1,m1,x2,m2,m0);
        }
    private:
        T x1;
        const GenSymBandMatrix<T1>& m1;
        T x2;
        const GenSymBandMatrix<T2>& m2;
    };

    template <typename T>
    inline SymBandMatrixView<T> operator+=(
        SymBandMatrixView<T> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.nlo() >= m2.nlo());
        TMVAssert(isReal(T()) || m1.sym() == m2.sym());
        AddMM(T(1),m2,m1); return m1;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator+=(
        SymBandMatrixView<CT> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.nlo() >= m2.nlo());
        AddMM(T(1),m2,m1); return m1;
    }

    template <typename T>
    inline SymBandMatrixView<T> operator-=(
        SymBandMatrixView<T> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.nlo() >= m2.nlo());
        TMVAssert(isReal(T()) || m1.sym() == m2.sym());
        AddMM(T(-1),m2,m1); return m1;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator-=(
        SymBandMatrixView<CT> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.nlo() >= m2.nlo());
        AddMM(T(-1),m2,m1); return m1;
    }

    template <typename T, typename T2>
    inline SymBandMatrixView<T> operator+=(
        SymBandMatrixView<T> m, const ProdXsB<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(m.nlo() >= pxm.nlo());
        TMVAssert(isReal(T2()) || m.issym() == pxm.getM().issym());
        TMVAssert(m.issym() || TMV_IMAG(pxm.getX()) == TMV_RealType(T)(0));
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator+=(
        SymBandMatrixView<CT> m, const ProdXsB<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(m.nlo() >= pxm.nlo());
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T, typename T2>
    inline SymBandMatrixView<T> operator-=(
        SymBandMatrixView<T> m, const ProdXsB<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(m.nlo() >= pxm.nlo());
        TMVAssert(isReal(T2()) || m.issym() == pxm.getM().issym());
        TMVAssert(m.issym() || TMV_IMAG(pxm.getX()) == TMV_RealType(T)(0));
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator-=(
        SymBandMatrixView<CT> m, const ProdXsB<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(m.nlo() >= pxm.nlo());
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

#define SUMMM SumsBsB
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXsB
#include "tmv/TMV_AuxSumMM.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // SymBandMatrix * SymBandMatrix
    //

    template <typename T, typename T1, typename T2>
    class ProdsBsB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdsBsB(
            T _x, const GenSymBandMatrix<T1>& _m1,
            const GenSymBandMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline ptrdiff_t nlo() const
        { return TMV_MIN(colsize()-1,m1.nlo()+m2.nlo()); }
        inline ptrdiff_t nhi() const
        { return TMV_MIN(rowsize()-1,m1.nhi()+m2.nhi()); }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline const GenSymBandMatrix<T2>& getM2() const { return m2; }
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
        const GenSymBandMatrix<T1>& m1;
        const GenSymBandMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, const ProdsBsB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());;
        TMVAssert(m.nhi() >= pmm.nhi());;
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, const ProdsBsB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());;
        TMVAssert(m.nhi() >= pmm.nhi());;
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, const ProdsBsB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());;
        TMVAssert(m.nhi() >= pmm.nhi());;
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, const ProdsBsB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());;
        TMVAssert(m.nhi() >= pmm.nhi());;
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ProdsBsB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixViewOf(m,pmm.nlo(),pmm.nhi()) += pmm;
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ProdsBsB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixViewOf(m,pmm.nlo(),pmm.nhi()) += pmm;
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ProdsBsB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixViewOf(m,pmm.nlo(),pmm.nhi()) -= pmm;
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ProdsBsB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixViewOf(m,pmm.nlo(),pmm.nhi()) -= pmm;
        return m;
    }

#define PRODMM ProdsBsB
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXsB
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2



    //
    // Element Product SymBandMatrix * SymBandMatrix
    //

    template <typename T, typename T1, typename T2>
    class ElemProdsBsB : public SymBandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ElemProdsBsB(
            T _x, const GenSymBandMatrix<T1>& _m1,
            const GenSymBandMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        {
            TMVAssert(m1.size() == m2.size());
            TMVAssert(m1.issym() == m2.issym());
            TMVAssert(m1.isherm() == m2.isherm());
        }
        inline ptrdiff_t size() const { return m1.size(); }
        inline ptrdiff_t nlo() const
        { return TMV_MIN(m1.nlo(),m2.nlo()); }
        inline SymType sym() const
        { return isReal(T1()) ? m2.sym() : m1.sym(); }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline const GenSymBandMatrix<T2>& getM2() const { return m2; }

        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nlo());
            ElemMultMM<false>(x,m1.upperBand(),m2.upperBand(),m0.upperBand());
            if (m0.nlo() > 0) {
                if (m1.nlo() > 0 && m2.nlo() > 0) {
                    ElemMultMM<false>(
                        x,m1.lowerBandOff(),m2.lowerBandOff(),
                        m0.lowerBandOff());
                } else {
                    m0.lowerBandOff().setZero();
                }
            }
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nlo());
            ElemMultMM<false>(x,m1.upperBand(),m2.upperBand(),m0.upperBand());
            if (m0.nlo() > 0) {
                if (m1.nlo() > 0 && m2.nlo() > 0) {
                    ElemMultMM<false>(
                        x,m1.lowerBandOff(),m2.lowerBandOff(),
                        m0.lowerBandOff());
                } else {
                    m0.lowerBandOff().setZero();
                }
            }
        }
        inline void assignTosB(SymBandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(isReal(T1()) || isReal(T2()) || m1.sym() == m2.sym());
            TMVAssert(m0.size() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(isReal(T1()) || m0.issym() == m1.issym());
            TMVAssert(isReal(T2()) || m0.issym() == m2.issym());
            TMVAssert(m0.issym() || TMV_IMAG(x) == real_type(0));
            ElemMultMM<false>(x,m1.upperBand(),m2.upperBand(),m0.upperBand());
        }
        inline void assignTosB(SymBandMatrixView<complex_type> m0) const
        {
            TMVAssert(isReal(T1()) || isReal(T2()) || m1.sym() == m2.sym());
            TMVAssert(m0.size() == size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(isReal(T1()) || m0.issym() == m1.issym());
            TMVAssert(isReal(T2()) || m0.issym() == m2.issym());
            TMVAssert(m0.issym() || TMV_IMAG(x) == real_type(0));
            ElemMultMM<false>(x,m1.upperBand(),m2.upperBand(),m0.upperBand());
        }
    private:
        T x;
        const GenSymBandMatrix<T1>& m1;
        const GenSymBandMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    inline SymBandMatrixView<T> operator+=(
        SymBandMatrixView<T> m, const ElemProdsBsB<T,T1,T2>& pmm)
    {
        TMVAssert(isReal(T1()) || isReal(T2()) ||
                  pmm.getM1().sym() == pmm.getM2().sym());
        TMVAssert(m.size() == pmm.size());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(isReal(T1()) || m.issym() == pmm.getM1().issym());
        TMVAssert(isReal(T2()) || m.issym() == pmm.getM2().issym());
        TMVAssert(m.issym() ||
                  TMV_IMAG(pmm.getX()) == typename Traits<T>::real_type(0));
        ElemMultMM<true>(
            pmm.getX(),pmm.getM1().upperBand(),pmm.getM2().upperBand(),
            m.upperBand());
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator+=(
        SymBandMatrixView<CT> m, const ElemProdsBsB<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.issym() || TMV_IMAG(pmm.getX()) == T(0));
        ElemMultMM<true>(
            pmm.getX(),pmm.getM1().upperBand(),pmm.getM2().upperBand(),
            m.upperBand());
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline SymBandMatrixView<T> operator-=(
        SymBandMatrixView<T> m, const ElemProdsBsB<T,T1,T2>& pmm)
    {
        TMVAssert(isReal(T1()) || isReal(T2()) ||
                  pmm.getM1().sym() == pmm.getM2().sym());
        TMVAssert(m.size() == pmm.size());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(isReal(T1()) || m.issym() == pmm.getM1().issym());
        TMVAssert(isReal(T2()) || m.issym() == pmm.getM2().issym());
        TMVAssert(m.issym() ||
                  TMV_IMAG(pmm.getX()) == typename Traits<T>::real_type(0));
        ElemMultMM<true>(
            -pmm.getX(),pmm.getM1().upperBand(),pmm.getM2().upperBand(),
            m.upperBand());
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator-=(
        SymBandMatrixView<CT> m, const ElemProdsBsB<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.issym() || TMV_IMAG(pmm.getX()) == T(0));
        ElemMultMM<true>(
            -pmm.getX(),pmm.getM1().upperBand(),pmm.getM2().upperBand(),
            m.upperBand());
        return m;
    }

#define PRODMM ElemProdsBsB
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXsB
#define OP ElemProd
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


#define GENMATRIX GenSymBandMatrix
#define PRODXM ProdXsB
#define PRODMV ProdsBV
#define PRODVM ProdVsB
#define QUOTVM QuotVsB
#define RQUOTVM RQuotVsB
#define QUOTXM QuotXsB
#include "tmv/TMV_AuxVecComposite.h"
#undef GENMATRIX
#undef PRODXM
#undef PRODMV
#undef PRODVM
#undef QUOTVM
#undef RQUOTVM
#undef QUOTXM


#define GENMATRIX GenSymBandMatrix
#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXsB
#define SUMMM SumMsB
#define PRODMM ProdMsB
#define QUOTXM QuotXsB
#define QUOTMM QuotMsB
#define RQUOTMM RQuotMsB
#define TQUOTMM TransientQuotMsB
#define TRQUOTMM TransientRQuotMsB

#include "tmv/TMV_AuxMatComposite1.h"
#include "tmv/TMV_AuxSumMMb.h"
#include "tmv/TMV_AuxMatComposite3.h"

    // sB/sB -> Matrix.  Use TransientQuot
#undef GENMATRIX
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenSymBandMatrix
#define PRODXM1 ProdXsB
#include "tmv/TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM


#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXM
#define SUMMMa SumMsB
#define SUMMM SumsBM
#define PRODMM ProdsBM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM
#define QUOTXM QuotXM

#include "tmv/TMV_AuxMatComposite2.h"
#include "tmv/TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
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
#include "tmv/TMV_DiagSymBandArith.h"
#endif

#ifdef TMV_TriMatrixArith_H
#include "tmv/TMV_TriSymBandArith.h"
#endif

#ifdef TMV_BandMatrixArith_H
#include "tmv/TMV_BandSymBandArith.h"
#endif

#ifdef TMV_SymMatrixArith_H
#include "tmv/TMV_SymSymBandArith.h"
#endif

#endif
