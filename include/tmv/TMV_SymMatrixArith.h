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


#ifndef TMV_SymMatrixArith_H
#define TMV_SymMatrixArith_H

#include "tmv/TMV_BaseSymMatrix.h"
#include "tmv/TMV_SymMatrixArithFunc.h"
#include "tmv/TMV_VIt.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

    template <typename T, typename Tv>
    class ProdXV;

    template <typename T, typename Tv>
    class ProdXM;

    template <typename T, typename Tv>
    class QuotXM;

    template <typename T, typename T1, typename T2>
    class OProdVV;

    template <typename T, typename T1, typename T2>
    class ProdMM;

    template <typename T, typename T1, typename T2>
    class ProdLU;

    template <typename T, typename T1, typename T2>
    class ProdUL;

    template <typename T, int A, typename Tx>
    inline SymMatrix<T,A>& operator+=(SymMatrix<T,A>& m, const Tx& x)
    { m.view() += x; return m; }

    template <typename T, int A, typename Tx>
    inline SymMatrix<T,A>& operator-=(SymMatrix<T,A>& m, const Tx& x)
    { m.view() -= x; return m; }

    template <typename T, int A, typename Tx>
    inline SymMatrix<T,A>& operator*=(SymMatrix<T,A>& m, const Tx& x)
    { m.view() *= x; return m; }

    template <typename T, int A, typename Tx>
    inline SymMatrix<T,A>& operator/=(SymMatrix<T,A>& m, const Tx& x)
    { m.view() /= x; return m; }

    template <typename T, int A, typename Tx>
    inline SymMatrix<T,A>& operator%=(SymMatrix<T,A>& m, const Tx& x)
    { m.view() %= x; return m; }

    template <typename T, int A, typename Tx>
    inline HermMatrix<T,A>& operator+=(HermMatrix<T,A>& m, const Tx& x)
    { m.view() += x; return m; }

    template <typename T, int A, typename Tx>
    inline HermMatrix<T,A>& operator-=(HermMatrix<T,A>& m, const Tx& x)
    { m.view() -= x; return m; }

    template <typename T, int A, typename Tx>
    inline HermMatrix<T,A>& operator*=(HermMatrix<T,A>& m, const Tx& x)
    { m.view() *= x; return m; }

    template <typename T, int A, typename Tx>
    inline HermMatrix<T,A>& operator/=(HermMatrix<T,A>& m, const Tx& x)
    { m.view() /= x; return m; }

    template <typename T, int A, typename Tx>
    inline HermMatrix<T,A>& operator%=(HermMatrix<T,A>& m, const Tx& x)
    { m.view() %= x; return m; }


    //
    // Scalar * SymMatrix
    //

    template <typename T, typename Tm>
    class ProdXS : public SymMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdXS(const T _x, const GenSymMatrix<Tm>& _m) :
            x(_x), m(_m) {}
        inline ptrdiff_t size() const { return m.size(); }
        inline SymType sym() const { return m.issym() ? Sym : Herm; }
        inline T getX() const { return x; }
        inline const GenSymMatrix<Tm>& getM() const { return m; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == m.size());
            TMVAssert(m0.rowsize() == m.size());
            m0.upperTri() = m.upperTri();
            if (m.size() > 0)
                m0.lowerTri().offDiag() = m.lowerTri().offDiag();
            MultXM(x,m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == m.size());
            TMVAssert(m0.rowsize() == m.size());
            m0.upperTri() = m.upperTri();
            if (m.size() > 0)
                m0.lowerTri().offDiag() = m.lowerTri().offDiag();
            MultXM(x,m0);
        }
        inline void assignToS(SymMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m.issym() || TMV_IMAG(x)==real_type(0) );
            TMVAssert(m0.size() == m.size());
            TMVAssert(isReal(Tm()) || m0.issym() == m.issym());
            MultXM(x,m0=m);
        }
        inline void assignToS(SymMatrixView<complex_type> m0) const
        {
            TMVAssert(m.issym() || TMV_IMAG(x)==real_type(0) );
            TMVAssert(m0.size() == m.size());
            TMVAssert(isReal(Tm()) || m0.issym() == m.issym());
            MultXM(x,m0=m);
        }
    private:
        const T x;
        const GenSymMatrix<Tm>& m;
    };

    // m*=x
    template <typename T>
    inline SymMatrixView<T> operator*=(SymMatrixView<T> m, T x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        MultXM(x,m);
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator*=(SymMatrixView<CT> m, T x)
    {
        MultXM(x,m);
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator*=(SymMatrixView<CT> m, CCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        MultXM(CT(x),m);
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator*=(SymMatrixView<CT> m, VCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        MultXM(CT(x),m);
        return m;
    }

    // m/=x
    template <typename T>
    inline SymMatrixView<T> operator/=(SymMatrixView<T> m, T x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        MultXM(TMV_InverseOf(x),m);
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator/=(SymMatrixView<CT> m, T x)
    {
        MultXM(TMV_InverseOf(x),m);
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator/=(
        SymMatrixView<CT> m, CCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        MultXM(TMV_InverseOf(CT(x)),m);
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator/=(
        SymMatrixView<CT> m, VCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        MultXM(TMV_InverseOf(CT(x)),m);
        return m;
    }

#define GENMATRIX GenSymMatrix
#define PRODXM ProdXS
#include "tmv/TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM

    //
    // SymMatrix + Scalar
    //

    template <typename T, typename Tm>
    class SumSX : public SymMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumSX(T _x1, const GenSymMatrix<Tm>& _m, T _x2) :
            x1(_x1), m(_m), x2(_x2) {}
        inline ptrdiff_t size() const { return m.size(); }
        inline SymType sym() const { return m.issym() ? Sym : Herm; }
        inline T getX1() const { return x1; }
        inline const GenSymMatrix<Tm>& getM() const { return m; }
        inline T getX2() const { return x2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == m.size());
            TMVAssert(m0.rowsize() == m.size());
            MultXM(x1,m0=m);
            m0.diag().addToAll(TMV_REAL(x2));
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == m.size());
            TMVAssert(m0.rowsize() == m.size());
            MultXM(x1,m0=m);
            m0.diag().addToAll(complex_type(x2));
        }
        inline void assignToS(SymMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(isReal(Tm()) || m0.issym() == m.issym());
            TMVAssert(m.issym() || TMV_IMAG(x1)==real_type(0));
            TMVAssert(m.issym() || TMV_IMAG(x2)==real_type(0));
            TMVAssert(m0.size() == m.size());
            MultXM(x1,m0=m);
            m0.diag().addToAll(TMV_REAL(x2));
        }
        inline void assignToS(SymMatrixView<complex_type> m0) const
        {
            TMVAssert(isReal(Tm()) || m0.issym() == m.issym());
            TMVAssert(m.issym() || TMV_IMAG(x1)==real_type(0));
            TMVAssert(m.issym() || TMV_IMAG(x2)==real_type(0));
            TMVAssert(m0.size() == m.size());
            MultXM(x1,m0=m);
            m0.diag().addToAll(complex_type(x2));
        }
    private:
        const T x1;
        const GenSymMatrix<Tm>& m;
        const T x2;
    };

    // m+=x
    template <typename T>
    inline SymMatrixView<T> operator+=(SymMatrixView<T> m, T x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        m.diag().addToAll(x);
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(SymMatrixView<CT> m, T x)
    {
        m.diag().addToAll(CT(x));
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(
        SymMatrixView<CT> m, CCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == T(0));
        m.diag().addToAll(CT(x));
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(
        SymMatrixView<CT> m, VCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == T(0));
        m.diag().addToAll(CT(x));
        return m;
    }

    // m-=x
    template <typename T>
    inline SymMatrixView<T> operator-=(SymMatrixView<T> m, T x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == TMV_RealType(T)(0));
        m.diag().addToAll(-x);
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(SymMatrixView<CT> m, T x)
    {
        m.diag().addToAll(CT(-x));
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(
        SymMatrixView<CT> m, CCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == T(0));
        m.diag().addToAll(-CT(x));
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(
        SymMatrixView<CT> m, VCT x)
    {
        TMVAssert(m.issym() || TMV_IMAG(x) == T(0));
        m.diag().addToAll(-CT(x));
        return m;
    }

#define SUMMX SumSX
#define GENMATRIX GenSymMatrix
#define PRODXM ProdXS
#include "tmv/TMV_AuxSumMX.h"
#undef SUMMX
#undef GENMATRIX
#undef PRODXM

    //
    // SymMatrix + SymMatrix
    //

    template <typename T, typename T1, typename T2>
    class SumSS : public SymMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumSS(
            T _x1, const GenSymMatrix<T1>& _m1,
            T _x2, const GenSymMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline SymType sym() const
        { return isReal(T1()) ? m2.sym() : m1.sym(); }
        inline T getX1() const { return x1; }
        inline const GenSymMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenSymMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            AddMM(x1,m1,x2,m2,m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            AddMM(x1,m1,x2,m2,m0);
        }
        inline void assignToS(SymMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(isReal(T1()) || isReal(T2()) || m1.sym() == m2.sym());
            TMVAssert(m0.size() == m1.size());
            TMVAssert(isReal(T1()) || m0.issym() == m1.issym());
            TMVAssert(isReal(T2()) || m0.issym() == m2.issym());
            TMVAssert(m0.issym() || TMV_IMAG(x1) == real_type(0));
            TMVAssert(m0.issym() || TMV_IMAG(x2) == real_type(0));
            AddMM(x1,m1,x2,m2,m0);
        }
        inline void assignToS(SymMatrixView<complex_type> m0) const
        {
            TMVAssert(isReal(T1()) || isReal(T2()) || m1.sym() == m2.sym());
            TMVAssert(m0.size() == m1.size());
            TMVAssert(isReal(T1()) || m0.issym() == m1.issym());
            TMVAssert(isReal(T2()) || m0.issym() == m2.issym());
            TMVAssert(m0.issym() || TMV_IMAG(x1) == real_type(0));
            TMVAssert(m0.issym() || TMV_IMAG(x2) == real_type(0));
            AddMM(x1,m1,x2,m2,m0);
        }
    private:
        T x1;
        const GenSymMatrix<T1>& m1;
        T x2;
        const GenSymMatrix<T2>& m2;
    };

    template <typename T>
    inline SymMatrixView<T> operator+=(
        SymMatrixView<T> m1, const GenSymMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(isReal(T()) || m1.sym() == m2.sym());
        AddMM(T(1),m2,m1); return m1;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(
        SymMatrixView<CT> m1, const GenSymMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        AddMM(T(1),m2,m1); return m1;
    }

    template <typename T>
    inline SymMatrixView<T> operator-=(
        SymMatrixView<T> m1, const GenSymMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(isReal(T()) || m1.sym() == m2.sym());
        AddMM(T(-1),m2,m1); return m1;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(
        SymMatrixView<CT> m1, const GenSymMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        AddMM(T(-1),m2,m1); return m1;
    }

    template <typename T, typename T2>
    inline SymMatrixView<T> operator+=(
        SymMatrixView<T> m, const ProdXS<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(isReal(T2()) || m.issym() == pxm.getM().issym());
        TMVAssert(m.issym() || TMV_IMAG(pxm.getX()) == TMV_RealType(T)(0));
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(
        SymMatrixView<CT> m, const ProdXS<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T, typename T2>
    inline SymMatrixView<T> operator-=(
        SymMatrixView<T> m, const ProdXS<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(isReal(T2()) || m.issym() == pxm.getM().issym());
        TMVAssert(m.issym() || TMV_IMAG(pxm.getX()) == TMV_RealType(T)(0));
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(
        SymMatrixView<CT> m, const ProdXS<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

#define SUMMM SumSS
#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXS
#define PRODXM2 ProdXS
#include "tmv/TMV_AuxSumMM.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // SymMatrix * SymMatrix
    //

    template <typename T, typename T1, typename T2>
    class ProdSS : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdSS(
            T _x, const GenSymMatrix<T1>& _m1,
            const GenSymMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline T getX() const { return x; }
        inline const GenSymMatrix<T1>& getM1() const { return m1; }
        inline const GenSymMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            MultMM<false>(x,m1,m2,m0);
        }
    private:
        T x;
        const GenSymMatrix<T1>& m1;
        const GenSymMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ProdSS<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ProdSS<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ProdSS<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ProdSS<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

#define PRODMM ProdSS
#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXS
#define PRODXM2 ProdXS
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // Element Product SymMatrix * SymMatrix
    //

    template <typename T, typename T1, typename T2>
    class ElemProdSS : public SymMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ElemProdSS(
            T _x, const GenSymMatrix<T1>& _m1, const GenSymMatrix<T2>& _m2) :
            x(_x),m1(_m1),m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline SymType sym() const
        { return isReal(T1()) ? m2.sym() : m1.sym(); }
        inline T getX() const { return x; }
        inline const GenSymMatrix<T1>& getM1() const { return m1; }
        inline const GenSymMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            ElemMultMM<false>(x,m1.upperTri(),m2.upperTri(),m0.upperTri());
            if (m1.size() > 1) {
                ElemMultMM<false>(
                    x,m1.lowerTri().offDiag(),m2.lowerTri().offDiag(),
                    m0.lowerTri().offDiag());
            }
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            ElemMultMM<false>(x,m1.upperTri(),m2.upperTri(),m0.upperTri());
            if (m1.size() > 1) {
                ElemMultMM<false>(
                    x,m1.lowerTri().offDiag(),m2.lowerTri().offDiag(),
                    m0.lowerTri().offDiag());
            }
        }
        inline void assignToS(SymMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(isReal(T1()) || isReal(T2()) || m1.sym() == m2.sym());
            TMVAssert(m0.size() == m1.size());
            TMVAssert(isReal(T1()) || m0.issym() == m1.issym());
            TMVAssert(isReal(T2()) || m0.issym() == m2.issym());
            TMVAssert(m0.issym() || TMV_IMAG(x) == real_type(0));
            ElemMultMM<false>(x,m1.upperTri(),m2.upperTri(),m0.upperTri());
        }
        inline void assignToS(SymMatrixView<complex_type> m0) const
        {
            TMVAssert(isReal(T1()) || isReal(T2()) || m1.sym() == m2.sym());
            TMVAssert(m0.size() == m1.size());
            TMVAssert(isReal(T1()) || m0.issym() == m1.issym());
            TMVAssert(isReal(T2()) || m0.issym() == m2.issym());
            TMVAssert(m0.issym() || TMV_IMAG(x) == real_type(0));
            ElemMultMM<false>(x,m1.upperTri(),m2.upperTri(),m0.upperTri());
        }
    private:
        T x;
        const GenSymMatrix<T1>& m1;
        const GenSymMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    inline SymMatrixView<T> operator+=(
        SymMatrixView<T> m, const ElemProdSS<T,T1,T2>& pmm)
    {
        TMVAssert(isReal(T1()) || isReal(T2()) ||
                  pmm.getM1().sym() == pmm.getM2().sym());
        TMVAssert(m.size() == pmm.size());
        TMVAssert(isReal(T1()) || m.issym() == pmm.getM1().issym());
        TMVAssert(isReal(T2()) || m.issym() == pmm.getM2().issym());
        TMVAssert(m.issym() ||
                  TMV_IMAG(pmm.getX()) == typename Traits<T>::real_type(0));
        ElemMultMM<true>(
            pmm.getX(),pmm.getM1().upperTri(),pmm.getM2().upperTri(),
            m.upperTri());
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(
        SymMatrixView<CT> m, const ElemProdSS<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(m.issym() || TMV_IMAG(pmm.getX()) == T(0));
        ElemMultMM<true>(
            pmm.getX(),pmm.getM1().upperTri(),pmm.getM2().upperTri(),
            m.upperTri());
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline SymMatrixView<T> operator-=(
        SymMatrixView<T> m, const ElemProdSS<T,T1,T2>& pmm)
    {
        TMVAssert(isReal(T1()) || isReal(T2()) ||
                  pmm.getM1().sym() == pmm.getM2().sym());
        TMVAssert(m.size() == pmm.size());
        TMVAssert(isReal(T1()) || m.issym() == pmm.getM1().issym());
        TMVAssert(isReal(T2()) || m.issym() == pmm.getM2().issym());
        TMVAssert(m.issym() ||
                  TMV_IMAG(pmm.getX()) == typename Traits<T>::real_type(0));
        ElemMultMM<true>(
            -pmm.getX(),pmm.getM1().upperTri(),pmm.getM2().upperTri(),
            m.upperTri());
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(
        SymMatrixView<CT> m, const ElemProdSS<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(m.issym() || TMV_IMAG(pmm.getX()) == T(0));
        ElemMultMM<true>(
            -pmm.getX(),pmm.getM1().upperTri(),pmm.getM2().upperTri(),
            m.upperTri());
        return m;
    }

#define PRODMM ElemProdSS
#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXS
#define PRODXM2 ProdXS
#define OP ElemProd
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2




    //
    // Vector ^ Vector
    //

    // m += (x*v^v)
    template <typename T, typename Tv>
    inline SymMatrixView<T> operator+=(
        SymMatrixView<T> m0, const OProdVV<T,Tv,Tv>& opvv)
    {
        TMVAssert(m0.size() == opvv.colsize());
        TMVAssert(m0.size() == opvv.rowsize());
        TMVAssert(opvv.getV1().isSameAs(
                m0.isherm() ? opvv.getV2().conjugate() : opvv.getV2().view()));
        Rank1Update<true>(opvv.getX(), opvv.getV1(), m0);
        return m0;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(
        SymMatrixView<CT> m0, const OProdVV<T,T,T>& opvv)
    {
        TMVAssert(m0.size() == opvv.colsize());
        TMVAssert(m0.size() == opvv.rowsize());
        TMVAssert(opvv.getV1().isSameAs(opvv.getV2()));
        Rank1Update<true>(opvv.getX(), opvv.getV1(), m0);
        return m0;
    }

    // m -= (x*v^v)
    template <typename T, typename Tv>
    inline SymMatrixView<T> operator-=(
        SymMatrixView<T> m0, const OProdVV<T,Tv,Tv>& opvv)
    {
        TMVAssert(m0.size() == opvv.colsize());
        TMVAssert(m0.size() == opvv.rowsize());
        TMVAssert(opvv.getV1().isSameAs(
                m0.isherm() ? opvv.getV2().conjugate() : opvv.getV2().view()));
        Rank1Update<true>(-opvv.getX(), opvv.getV1(), m0);
        return m0;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(
        SymMatrixView<CT> m0, const OProdVV<T,T,T>& opvv)
    {
        TMVAssert(m0.size() == opvv.colsize());
        TMVAssert(m0.size() == opvv.rowsize());
        TMVAssert(opvv.getV1().isSameAs(opvv.getV2()));
        Rank1Update<true>(-opvv.getX(), opvv.getV1(), m0);
        return m0;
    }

    //
    // Matrix * Matrix.transpose()
    //

    // m += (x*m*mt)
    template <typename T, typename Tm>
    inline SymMatrixView<T> operator+=(
        SymMatrixView<T> m0, const ProdMM<T,Tm,Tm>& opmm)
    {
        TMVAssert(m0.size() == opmm.colsize());
        TMVAssert(m0.size() == opmm.rowsize());
        TMVAssert(opmm.getM1().isSameAs(
                m0.isherm() ?
                opmm.getM2().adjoint() :
                opmm.getM2().transpose()));
        RankKUpdate<true>(opmm.getX(), opmm.getM1(), m0);
        return m0;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(
        SymMatrixView<CT> m0, const ProdMM<T,T,T>& opmm)
    {
        TMVAssert(m0.size() == opmm.colsize());
        TMVAssert(m0.size() == opmm.rowsize());
        TMVAssert(opmm.getM1().isSameAs(opmm.getM2().transpose()));
        RankKUpdate<true>(opmm.getX(), opmm.getM1(), m0);
        return m0;
    }

    template <typename T, typename Tm>
    inline SymMatrixView<T> operator+=(
        SymMatrixView<T> m0, const ProdLU<T,Tm,Tm>& opmm)
    {
        TMVAssert(m0.size() == opmm.colsize());
        TMVAssert(m0.size() == opmm.rowsize());
        TMVAssert(opmm.getM1().isSameAs(
                m0.isherm() ?
                opmm.getM2().adjoint() :
                opmm.getM2().transpose()));
        RankKUpdate<true>(opmm.getX(), opmm.getM1(), m0);
        return m0;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(
        SymMatrixView<CT> m0, const ProdLU<T,T,T>& opmm)
    {
        TMVAssert(m0.size() == opmm.colsize());
        TMVAssert(m0.size() == opmm.rowsize());
        TMVAssert(opmm.getM1().isSameAs(opmm.getM2().transpose()));
        RankKUpdate<true>(opmm.getX(), opmm.getM1(), m0);
        return m0;
    }

    template <typename T, typename Tm>
    inline SymMatrixView<T> operator+=(
        SymMatrixView<T> m0, const ProdUL<T,Tm,Tm>& opmm)
    {
        TMVAssert(m0.size() == opmm.colsize());
        TMVAssert(m0.size() == opmm.rowsize());
        TMVAssert(opmm.getM1().isSameAs(
                m0.isherm() ?
                opmm.getM2().adjoint() :
                opmm.getM2().transpose()));
        RankKUpdate<true>(opmm.getX(), opmm.getM1(), m0);
        return m0;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(
        SymMatrixView<CT> m0, const ProdUL<T,T,T>& opmm)
    {
        TMVAssert(m0.size() == opmm.colsize());
        TMVAssert(m0.size() == opmm.rowsize());
        TMVAssert(opmm.getM1().isSameAs(opmm.getM2().transpose()));
        RankKUpdate<true>(opmm.getX(), opmm.getM1(), m0);
        return m0;
    }

    // m -= (x*m*mt)
    template <typename T, typename Tm>
    inline SymMatrixView<T> operator-=(
        SymMatrixView<T> m0, const ProdMM<T,Tm,Tm>& opmm)
    {
        TMVAssert(m0.size() == opmm.colsize());
        TMVAssert(m0.size() == opmm.rowsize());
        TMVAssert(opmm.getM1().isSameAs(
                m0.isherm() ?
                opmm.getM2().adjoint() :
                opmm.getM2().transpose()));
        RankKUpdate<true>(-opmm.getX(), opmm.getM1(), m0);
        return m0;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(
        SymMatrixView<CT> m0, const ProdMM<T,T,T>& opmm)
    {
        TMVAssert(m0.size() == opmm.colsize());
        TMVAssert(m0.size() == opmm.rowsize());
        TMVAssert(opmm.getM1().isSameAs(opmm.getM2().transpose()));
        RankKUpdate<true>(-opmm.getX(), opmm.getM1(), m0);
        return m0;
    }

    template <typename T, typename Tm>
    inline SymMatrixView<T> operator-=(
        SymMatrixView<T> m0, const ProdUL<T,Tm,Tm>& opmm)
    {
        TMVAssert(m0.size() == opmm.colsize());
        TMVAssert(m0.size() == opmm.rowsize());
        TMVAssert(opmm.getM1().isSameAs(
                m0.isherm() ?
                opmm.getM2().adjoint() :
                opmm.getM2().transpose()));
        RankKUpdate<true>(-opmm.getX(), opmm.getM1(), m0);
        return m0;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(
        SymMatrixView<CT> m0, const ProdUL<T,T,T>& opmm)
    {
        TMVAssert(m0.size() == opmm.colsize());
        TMVAssert(m0.size() == opmm.rowsize());
        TMVAssert(opmm.getM1().isSameAs(opmm.getM2().transpose()));
        RankKUpdate<true>(-opmm.getX(), opmm.getM1(), m0);
        return m0;
    }

    template <typename T, typename Tm>
    inline SymMatrixView<T> operator-=(
        SymMatrixView<T> m0, const ProdLU<T,Tm,Tm>& opmm)
    {
        TMVAssert(m0.size() == opmm.colsize());
        TMVAssert(m0.size() == opmm.rowsize());
        TMVAssert(opmm.getM1().isSameAs(
                m0.isherm() ?
                opmm.getM1().adjoint() :
                opmm.getM2().transpose()));
        RankKUpdate<true>(-opmm.getX(), opmm.getM1(), m0);
        return m0;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(
        SymMatrixView<CT> m0, const ProdLU<T,T,T>& opmm)
    {
        TMVAssert(m0.size() == opmm.colsize());
        TMVAssert(m0.size() == opmm.rowsize());
        TMVAssert(opmm.getM1().isSameAs(opmm.getM2().transpose()));
        RankKUpdate<true>(-opmm.getX(), opmm.getM1(), m0);
        return m0;
    }



#define GENMATRIX GenSymMatrix
#define PRODXM ProdXS
#define QUOTXM QuotXS
#define PRODMV ProdSV
#define PRODVM ProdVS
#define QUOTVM QuotVS
#define RQUOTVM RQuotVS
#include "tmv/TMV_AuxVecComposite.h"
#undef GENMATRIX
#undef PRODXM
#undef QUOTXM
#undef PRODMV
#undef PRODVM
#undef QUOTVM
#undef RQUOTVM

#define GENMATRIX GenSymMatrix
#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXS
#define SUMMM SumMS
#define PRODMM ProdMS
#define QUOTXM QuotXS
#define QUOTMM QuotMS
#define RQUOTMM RQuotMS
#define TQUOTMM TransientQuotMS
#define TRQUOTMM TransientRQuotMS

#include "tmv/TMV_AuxMatComposite1.h"
#include "tmv/TMV_AuxSumMMb.h"
#include "tmv/TMV_AuxMatComposite3.h"
    // S/S -> Matrix.  Use TransientQuot
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenSymMatrix
#define PRODXM1 ProdXS
#include "tmv/TMV_AuxTQuotMM.h"

#undef GENMATRIX
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
#define SUMMMa SumMS
#define SUMMM SumSM
#define PRODMM ProdSM
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
#include "tmv/TMV_DiagSymArith.h"
#endif

#ifdef TMV_TriMatrixArith_H
#include "tmv/TMV_TriSymArith.h"
#endif

#ifdef TMV_BandMatrixArith_H
#include "tmv/TMV_BandSymArith.h"
#endif

#ifdef TMV_SymBandMatrixArith_H
#include "tmv/TMV_SymSymBandArith.h"
#endif

#endif
