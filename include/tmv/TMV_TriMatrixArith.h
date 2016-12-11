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


#ifndef TMV_TriMatrixArith_H
#define TMV_TriMatrixArith_H

#include "tmv/TMV_TriMatrixArithFunc.h"
#include "tmv/TMV_VectorArithFunc.h"
#include "tmv/TMV_MatrixArithFunc.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

    template <typename T, typename Tv>
    class ProdXV;
    template <typename T, typename Tv>
    class ProdXM;

    //
    // First do everything which returns an UpperTriMatrixComposite
    // or only deals with an UpperTriMatrix, and not a LowerTriMatrix
    //

    template <typename T, int A, typename Tx>
    inline UpperTriMatrix<T,A>& operator+=(
        UpperTriMatrix<T,A>& m, const Tx& x)
    { m.view() += x; return m; }

    template <typename T, int A, typename Tx>
    inline UpperTriMatrix<T,A>& operator-=(
        UpperTriMatrix<T,A>& m, const Tx& x)
    { m.view() -= x; return m; }

    template <typename T, int A, typename Tx>
    inline UpperTriMatrix<T,A>& operator*=(
        UpperTriMatrix<T,A>& m, const Tx& x)
    { m.view() *= x; return m; }

    template <typename T, int A, typename Tx>
    inline UpperTriMatrix<T,A>& operator/=(
        UpperTriMatrix<T,A>& m, const Tx& x)
    { m.view() /= x; return m; }

    template <typename T, int A, typename Tx>
    inline UpperTriMatrix<T,A>& operator%=(
        UpperTriMatrix<T,A>& m, const Tx& x)
    { m.view() %= x; return m; }

    template <typename T, typename T2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m1, const GenUpperTriMatrix<T2>& m2)
    {
        TMVAssert(m1.isSquare());
        TMVAssert(m1.colsize() == m2.size());
        m1.upperTri() += m2;
        return m1;
    }

    template <typename T, typename T2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m1, const GenUpperTriMatrix<T2>& m2)
    {
        TMVAssert(m1.isSquare());
        TMVAssert(m1.colsize() == m2.size());
        m1.upperTri() -= m2;
        return m1;
    }

    //
    // Scalar * TriMatrix
    //

    template <typename T, typename Tm>
    class ProdXU : public UpperTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdXU(const T _x, const GenUpperTriMatrix<Tm>& _m) :
            x(_x), m(_m) {}
        inline ptrdiff_t size() const { return m.size(); }
        inline DiagType dt() const
        { return x==T(1) ? m.dt() : NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenUpperTriMatrix<Tm>& getM() const { return m; }
        inline void assignToU(UpperTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            MultXM(x,m0=m);
        }
        inline void assignToU(UpperTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            MultXM(x,m0=m);
        }
    private:
        const T x;
        const GenUpperTriMatrix<Tm>& m;
    };

    // m*=x
    template <typename T>
    inline UpperTriMatrixView<T> operator*=(UpperTriMatrixView<T> m, T x)
    {
        TMVAssert(!m.isunit());
        MultXM(x,m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator*=(UpperTriMatrixView<CT> m, T x)
    {
        TMVAssert(!m.isunit());
        MultXM(x,m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator*=(UpperTriMatrixView<CT> m, CCT x)
    {
        TMVAssert(!m.isunit());
        MultXM(CT(x),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator*=(UpperTriMatrixView<CT> m, VCT x)
    {
        TMVAssert(!m.isunit());
        MultXM(CT(x),m);
        return m;
    }

    // m/=x
    template <typename T>
    inline UpperTriMatrixView<T> operator/=(UpperTriMatrixView<T> m, T x)
    {
        TMVAssert(!m.isunit());
        MultXM(TMV_InverseOf(x),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator/=(UpperTriMatrixView<CT> m, T x)
    {
        TMVAssert(!m.isunit());
        MultXM(TMV_InverseOf(x),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator/=(UpperTriMatrixView<CT> m, CCT x)
    {
        TMVAssert(!m.isunit());
        MultXM(TMV_InverseOf(CT(x)),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator/=(UpperTriMatrixView<CT> m, VCT x)
    {
        TMVAssert(!m.isunit());
        MultXM(TMV_InverseOf(CT(x)),m);
        return m;
    }

#define GENMATRIX GenUpperTriMatrix
#define PRODXM ProdXU
#include "tmv/TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM

    //
    // TriMatrix + Scalar
    //

    template <typename T, typename Tm>
    class SumUX : public UpperTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumUX(T _x1, const GenUpperTriMatrix<Tm>& _m, T _x2) :
            x1(_x1), m(_m), x2(_x2) {}
        inline ptrdiff_t size() const { return m.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX1() const { return x1; }
        inline const GenUpperTriMatrix<Tm>& getM() const { return m; }
        inline T getX2() const { return x2; }
        inline void assignToU(UpperTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit());
            MultXM(x1,m0=m);
            m0.diag().addToAll(TMV_REAL(x2));
        }
        inline void assignToU(UpperTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit());
            MultXM(x1,m0=m);
            m0.diag().addToAll(complex_type(x2));
        }
    private:
        const T x1;
        const GenUpperTriMatrix<Tm>& m;
        const T x2;
    };

    // m+=x
    template <typename T>
    inline UpperTriMatrixView<T> operator+=(UpperTriMatrixView<T> m, T x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(x);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(UpperTriMatrixView<CT> m, T x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(CT(x));
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(UpperTriMatrixView<CT> m, CCT x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(CT(x));
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(UpperTriMatrixView<CT> m, VCT x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(CT(x));
        return m;
    }

    // m-=x
    template <typename T>
    inline UpperTriMatrixView<T> operator-=(UpperTriMatrixView<T> m, T x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(-x);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(UpperTriMatrixView<CT> m, T x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(CT(-x));
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(UpperTriMatrixView<CT> m, CCT x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(-CT(x));
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(UpperTriMatrixView<CT> m, VCT x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(-CT(x));
        return m;
    }

#define SUMMX SumUX
#define GENMATRIX GenUpperTriMatrix
#define PRODXM ProdXU
#include "tmv/TMV_AuxSumMX.h"
#undef SUMMX
#undef GENMATRIX
#undef PRODXM

    //
    // TriMatrix + TriMatrix
    //

    template <typename T, typename T1, typename T2>
    class SumUU : public UpperTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumUU(T _x1, const GenUpperTriMatrix<T1>& _m1,
                     T _x2, const GenUpperTriMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX1() const { return x1; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToU(UpperTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            AddMM(x1,m1,x2,m2,m0);
        }
        inline void assignToU(UpperTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            AddMM(x1,m1,x2,m2,m0);
        }
    private:
        const T x1;
        const GenUpperTriMatrix<T1>& m1;
        const T x2;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T>
    inline UpperTriMatrixView<T> operator+=(
        UpperTriMatrixView<T> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddMM(T(1),m2,m1);
        return m1;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(
        UpperTriMatrixView<CT> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddMM(T(1),m2,m1);
        return m1;
    }

    template <typename T>
    inline UpperTriMatrixView<T> operator-=(
        UpperTriMatrixView<T> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddMM(T(-1),m2,m1);
        return m1;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(
        UpperTriMatrixView<CT> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddMM(T(-1),m2,m1);
        return m1;
    }

    template <typename T, typename T2>
    inline UpperTriMatrixView<T> operator+=(
        UpperTriMatrixView<T> m, const ProdXU<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(
        UpperTriMatrixView<CT> m, const ProdXU<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T, typename T2>
    inline UpperTriMatrixView<T> operator-=(
        UpperTriMatrixView<T> m, const ProdXU<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(
        UpperTriMatrixView<CT> m, const ProdXU<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

#define SUMMM SumUU
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXU
#include "tmv/TMV_AuxSumMM.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // TriMatrix * TriMatrix
    //

    template <typename T, typename T1, typename T2>
    class ProdUU : public UpperTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdUU(
            T _x, const GenUpperTriMatrix<T1>& _m1,
            const GenUpperTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const
        { return x==T(1) && m1.dt()==m2.dt() ? m1.dt() : NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToU(UpperTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToU(UpperTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            MultMM<false>(x,m1,m2,m0);
        }
    protected:
        T x;
        const GenUpperTriMatrix<T1>& m1;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T>
    inline UpperTriMatrixView<T> operator*=(
        UpperTriMatrixView<T> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit() || m2.isunit());
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator*=(
        UpperTriMatrixView<CT> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit() || m2.isunit());
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T, typename T1, typename T2>
    inline UpperTriMatrixView<T> operator+=(
        UpperTriMatrixView<T> m, const ProdUU<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(
        UpperTriMatrixView<CT> m, const ProdUU<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline UpperTriMatrixView<T> operator-=(
        UpperTriMatrixView<T> m, const ProdUU<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(
        UpperTriMatrixView<CT> m, const ProdUU<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }


#define PRODMM ProdUU
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXU
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // Element Product TriMatrix * TriMatrix
    //

    template <typename T, typename T1, typename T2>
    class ElemProdUU : public UpperTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ElemProdUU(
            T _x, const GenUpperTriMatrix<T1>& _m1,
            const GenUpperTriMatrix<T2>& _m2) :
            x(_x),m1(_m1),m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const
        { return x==T(1) && m1.dt()==m2.dt() ? m1.dt() : NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToU(UpperTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            ElemMultMM<false>(x,m1,m2,m0);
        }
        inline void assignToU(UpperTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            ElemMultMM<false>(x,m1,m2,m0);
        }
    private:
        const T x;
        const GenUpperTriMatrix<T1>& m1;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    inline UpperTriMatrixView<T> operator+=(
        UpperTriMatrixView<T> m, const ElemProdUU<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit());
        ElemMultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(
        UpperTriMatrixView<CT> m, const ElemProdUU<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit());
        ElemMultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline UpperTriMatrixView<T> operator-=(
        UpperTriMatrixView<T> m, const ElemProdUU<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit());
        ElemMultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(
        UpperTriMatrixView<CT> m, const ElemProdUU<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit());
        ElemMultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }


#define PRODMM ElemProdUU
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXU
#define OP ElemProd
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2



    //
    // Scalar / TriMatrix
    //

    template <typename T, typename Tm>
    class QuotXU : public UpperTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotXU(const T _x, const GenUpperTriMatrix<Tm>& _m) :
            x(_x), m(_m)  {}
        inline ptrdiff_t size() const { return m.size(); }
        inline DiagType dt() const
        { return x==T(1) ? m.dt() : NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenUpperTriMatrix<Tm>& getM() const { return m; }
        inline void assignToU(UpperTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            (m0=m).invertSelf();
            MultXM(x,m0);
        }
        inline void assignToU(UpperTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            (m0=m).invertSelf();
            MultXM(x,m0);
        }
    private:
        const T x;
        const GenUpperTriMatrix<Tm>& m;
    };

#define QUOTXM QuotXU
#define GENMATRIX GenUpperTriMatrix
#define PRODXM ProdXU
#include "tmv/TMV_AuxQuotXM.h"
#include "tmv/TMV_AuxQuotXMa.h"
#undef QUOTXM
#undef GENMATRIX
#undef PRODXM

    //
    // TriMatrix / % TriMatrix
    //

    template <typename T, typename T1, typename T2>
    class QuotUU : public UpperTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotUU(
            const T _x, const GenUpperTriMatrix<T1>& _m1,
            const GenUpperTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const
        { return x==T(1) && m1.dt()==m2.dt() ? m1.dt() : NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToU(UpperTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
            m2.LDiv(m1,m0);
            MultXM(x,m0);
        }
        inline void assignToU(UpperTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
            m2.LDiv(m1,m0);
            MultXM(x,m0);
        }
    protected:
        const T x;
        const GenUpperTriMatrix<T1>& m1;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class RQuotUU : public UpperTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline RQuotUU(
            const T _x, const GenUpperTriMatrix<T1>& _m1,
            const GenUpperTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const
        { return x==T(1) && m1.dt()==m2.dt() ? m1.dt() : NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToU(UpperTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
            m2.RDiv(m1,m0);
            MultXM(x,m0);
        }
        inline void assignToU(UpperTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
            m2.RDiv(m1,m0);
            MultXM(x,m0);
        }
    protected:
        const T x;
        const GenUpperTriMatrix<T1>& m1;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T>
    inline UpperTriMatrixView<T> operator/=(
        UpperTriMatrixView<T> m1, const GenUpperTriMatrix<T>& m2)
    { m2.LDivEq(m1); return m1; }

    template <typename T>
    inline UpperTriMatrixView<CT> operator/=(
        UpperTriMatrixView<CT> m1, const GenUpperTriMatrix<T>& m2)
    { m2.LDivEq(m1); return m1; }

    template <typename T>
    inline UpperTriMatrixView<T> operator%=(
        UpperTriMatrixView<T> m1, const GenUpperTriMatrix<T>& m2)
    { m2.RDivEq(m1); return m1; }

    template <typename T>
    inline UpperTriMatrixView<CT> operator%=(
        UpperTriMatrixView<CT> m1, const GenUpperTriMatrix<T>& m2)
    { m2.RDivEq(m1); return m1; }

#define QUOTMM QuotUU
#define QUOTXM QuotXU
#define RQUOTMM RQuotUU
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXU
#include "tmv/TMV_AuxQuotMM.h"
#include "tmv/TMV_AuxQuotMMa.h"
#undef QUOTMM
#undef QUOTXM
#undef RQUOTMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX GenUpperTriMatrix
#define PRODXM ProdXU
#define PRODMV ProdUV
#define PRODVM ProdVU
#define QUOTVM QuotVU
#define QUOTXM QuotXU
#define RQUOTVM RQuotVU
#include "tmv/TMV_AuxVecComposite.h"
#undef GENMATRIX
#undef PRODXM
#undef PRODMV
#undef PRODVM
#undef QUOTVM
#undef QUOTXM
#undef RQUOTVM


#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXU
#define SUMMM SumMU
#define PRODMM ProdMU
#define QUOTMM QuotMU
#define QUOTXM QuotXU
#define RQUOTMM RQuotMU
#define TQUOTMM TransientQuotMU
#define TRQUOTMM TransientRQuotMU

#include "tmv/TMV_AuxMatComposite1.h"
#include "tmv/TMV_AuxSumMMb.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef QUOTXM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM


#define GENMATRIX1 GenUpperTriMatrix
#define PRODXM1 ProdXU
#define GENMATRIX2 GenMatrix
#define PRODXM2 ProdXM
#define SUMMMa SumMU
#define SUMMM SumUM
#define PRODMM ProdUM
#define QUOTXM QuotXM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM

#include "tmv/TMV_AuxMatComposite2.h"
#include "tmv/TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef SUMMMa
#undef PRODMM
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM

    //
    // Now the things which return a LowerTriMatrixComposite
    // or only deal with a LowerTriMatrix
    //

    template <typename T, int A, typename Tx>
    inline LowerTriMatrix<T,A>& operator+=(
        LowerTriMatrix<T,A>& m, const Tx& x)
    { m.view() += x; return m; }

    template <typename T, int A, typename Tx>
    inline LowerTriMatrix<T,A>& operator-=(
        LowerTriMatrix<T,A>& m, const Tx& x)
    { m.view() -= x; return m; }

    template <typename T, int A, typename Tx>
    inline LowerTriMatrix<T,A>& operator*=(
        LowerTriMatrix<T,A>& m, const Tx& x)
    { m.view() *= x; return m; }

    template <typename T, int A, typename Tx>
    inline LowerTriMatrix<T,A>& operator/=(
        LowerTriMatrix<T,A>& m, const Tx& x)
    { m.view() /= x; return m; }

    template <typename T, int A, typename Tx>
    inline LowerTriMatrix<T,A>& operator%=(
        LowerTriMatrix<T,A>& m, const Tx& x)
    { m.view() %= x; return m; }

    template <typename T, typename T2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m1, const GenLowerTriMatrix<T2>& m2)
    {
        TMVAssert(m1.isSquare());
        TMVAssert(m1.colsize() == m2.size());
        m1.lowerTri() += m2;
        return m1;
    }

    template <typename T, typename T2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m1, const GenLowerTriMatrix<T2>& m2)
    {
        TMVAssert(m1.isSquare());
        TMVAssert(m1.colsize() == m2.size());
        m1.lowerTri() -= m2;
        return m1;
    }

    //
    // Scalar * TriMatrix
    //

    template <typename T, typename Tm>
    class ProdXL : public LowerTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdXL(const T _x, const GenLowerTriMatrix<Tm>& _m) :
            x(_x), m(_m) {}
        inline ptrdiff_t size() const { return m.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenLowerTriMatrix<Tm>& getM() const { return m; }
        inline void assignToL(LowerTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            MultXM(x,m0=m);
        }
        inline void assignToL(LowerTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            MultXM(x,m0=m);
        }
    private:
        const T x;
        const GenLowerTriMatrix<Tm>& m;
    };

    // m*=x
    template <typename T>
    inline LowerTriMatrixView<T> operator*=(LowerTriMatrixView<T> m, T x)
    {
        TMVAssert(!m.isunit());
        MultXM(x,m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator*=(LowerTriMatrixView<CT> m, T x)
    {
        TMVAssert(!m.isunit());
        MultXM(x,m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator*=(LowerTriMatrixView<CT> m, CCT x)
    {
        TMVAssert(!m.isunit());
        MultXM(CT(x),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator*=(LowerTriMatrixView<CT> m, VCT x)
    {
        TMVAssert(!m.isunit());
        MultXM(CT(x),m);
        return m;
    }

    // m/=x
    template <typename T>
    inline LowerTriMatrixView<T> operator/=(LowerTriMatrixView<T> m, T x)
    {
        TMVAssert(!m.isunit());
        MultXM(TMV_InverseOf(x),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator/=(LowerTriMatrixView<CT> m, T x)
    {
        TMVAssert(!m.isunit());
        MultXM(TMV_InverseOf(x),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator/=(LowerTriMatrixView<CT> m, CCT x)
    {
        TMVAssert(!m.isunit());
        MultXM(TMV_InverseOf(CT(x)),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator/=(LowerTriMatrixView<CT> m, VCT x)
    {
        TMVAssert(!m.isunit());
        MultXM(TMV_InverseOf(CT(x)),m);
        return m;
    }

#define GENMATRIX GenLowerTriMatrix
#define PRODXM ProdXL
#include "tmv/TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODX


    //
    // TriMatrix + Scalar
    //

    template <typename T, typename Tm>
    class SumLX : public LowerTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumLX(T _x1, const GenLowerTriMatrix<Tm>& _m, T _x2) :
            x1(_x1), m(_m), x2(_x2) {}
        inline ptrdiff_t size() const { return m.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX1() const { return x1; }
        inline const GenLowerTriMatrix<Tm>& getM() const { return m; }
        inline T getX2() const { return x2; }
        inline void assignToL(LowerTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit());
            MultXM(x1,m0=m);
            m0.diag().addToAll(TMV_REAL(x2));
        }
        inline void assignToL(LowerTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit());
            MultXM(x1,m0=m);
            m0.diag().addToAll(complex_type(x2));
        }
    private:
        const T x1;
        const GenLowerTriMatrix<Tm>& m;
        const T x2;
    };

    // m+=x
    template <typename T>
    inline LowerTriMatrixView<T> operator+=(LowerTriMatrixView<T> m, T x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(x);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(LowerTriMatrixView<CT> m, T x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(CT(x));
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(LowerTriMatrixView<CT> m, CCT x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(CT(x));
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(LowerTriMatrixView<CT> m, VCT x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(CT(x));
        return m;
    }

    // m-=x
    template <typename T>
    inline LowerTriMatrixView<T> operator-=(LowerTriMatrixView<T> m, T x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(-x);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(LowerTriMatrixView<CT> m, T x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(CT(-x));
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(LowerTriMatrixView<CT> m, CCT x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(-CT(x));
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(LowerTriMatrixView<CT> m, VCT x)
    {
        TMVAssert(!m.isunit());
        m.diag().addToAll(-CT(x));
        return m;
    }

#define SUMMX SumLX
#define GENMATRIX GenLowerTriMatrix
#define PRODXM ProdXL
#include "tmv/TMV_AuxSumMX.h"
#undef SUMMX
#undef GENMATRIX
#undef PRODX


    //
    // TriMatrix + TriMatrix
    //

    template <typename T, typename T1, typename T2>
    class SumLL : public LowerTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumLL(
            T _x1, const GenLowerTriMatrix<T1>& _m1,
            T _x2, const GenLowerTriMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m2.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX1() const { return x1; }
        inline const GenLowerTriMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToL(LowerTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            AddMM(x1,m1,x2,m2,m0);
        }
        inline void assignToL(LowerTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            AddMM(x1,m1,x2,m2,m0);
        }
    private:
        T x1;
        const GenLowerTriMatrix<T1>& m1;
        T x2;
        const GenLowerTriMatrix<T2>& m2;
    };

    template <typename T>
    inline LowerTriMatrixView<T> operator+=(
        LowerTriMatrixView<T> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddMM(T(1),m2,m1);
        return m1;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(
        LowerTriMatrixView<CT> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddMM(T(1),m2,m1);
        return m1;
    }

    template <typename T>
    inline LowerTriMatrixView<T> operator-=(
        LowerTriMatrixView<T> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddMM(T(-1),m2,m1);
        return m1;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(
        LowerTriMatrixView<CT> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddMM(T(-1),m2,m1);
        return m1;
    }

    template <typename T, typename T2>
    inline LowerTriMatrixView<T> operator+=(
        LowerTriMatrixView<T> m, const ProdXL<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(
        LowerTriMatrixView<CT> m, const ProdXL<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T, typename T2>
    inline LowerTriMatrixView<T> operator-=(
        LowerTriMatrixView<T> m, const ProdXL<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(
        LowerTriMatrixView<CT> m, const ProdXL<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

#define SUMMM SumLL
#define GENMATRIX1 GenLowerTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXL
#define PRODXM2 ProdXL
#include "tmv/TMV_AuxSumMM.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // TriMatrix * TriMatrix
    //

    template <typename T, typename T1, typename T2>
    class ProdLL : public LowerTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdLL(
            T _x, const GenLowerTriMatrix<T1>& _m1,
            const GenLowerTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m2.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const
        { return x==T(1) && m1.dt()==m2.dt() ? m1.dt() : NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenLowerTriMatrix<T1>& getM1() const { return m1; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToL(LowerTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToL(LowerTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            MultMM<false>(x,m1,m2,m0);
        }
    protected:
        T x;
        const GenLowerTriMatrix<T1>& m1;
        const GenLowerTriMatrix<T2>& m2;
    };

    template <typename T>
    inline LowerTriMatrixView<T> operator*=(
        LowerTriMatrixView<T> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit() || m2.isunit());
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator*=(
        LowerTriMatrixView<CT> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit() || m2.isunit());
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T, typename T1, typename T2>
    inline LowerTriMatrixView<T> operator+=(
        LowerTriMatrixView<T> m, const ProdLL<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit() || pmm.isunit());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(
        LowerTriMatrixView<CT> m, const ProdLL<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit() || pmm.isunit());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline LowerTriMatrixView<T> operator-=(
        LowerTriMatrixView<T> m, const ProdLL<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit() || pmm.isunit());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(
        LowerTriMatrixView<CT> m, const ProdLL<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit() || pmm.isunit());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

#define PRODMM ProdLL
#define GENMATRIX1 GenLowerTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXL
#define PRODXM2 ProdXL
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // Element Product TriMatrix * TriMatrix
    //

    template <typename T, typename T1, typename T2>
    class ElemProdLL : public LowerTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ElemProdLL(
            T _x, const GenLowerTriMatrix<T1>& _m1,
            const GenLowerTriMatrix<T2>& _m2) :
            x(_x),m1(_m1),m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const
        { return x==T(1) && m1.dt()==m2.dt() ? m1.dt() : NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenLowerTriMatrix<T1>& getM1() const { return m1; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToL(LowerTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            ElemMultMM<false>(x,m1,m2,m0);
        }
        inline void assignToL(LowerTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            ElemMultMM<false>(x,m1,m2,m0);
        }
    private:
        const T x;
        const GenLowerTriMatrix<T1>& m1;
        const GenLowerTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    inline LowerTriMatrixView<T> operator+=(
        LowerTriMatrixView<T> m, const ElemProdLL<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit());
        ElemMultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(
        LowerTriMatrixView<CT> m, const ElemProdLL<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit());
        ElemMultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline LowerTriMatrixView<T> operator-=(
        LowerTriMatrixView<T> m, const ElemProdLL<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit());
        ElemMultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(
        LowerTriMatrixView<CT> m, const ElemProdLL<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.size());
        TMVAssert(!m.isunit());
        ElemMultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }


#define PRODMM ElemProdLL
#define GENMATRIX1 GenLowerTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXL
#define PRODXM2 ProdXL
#define OP ElemProd
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2



    //
    // Scalar / TriMatrix
    //

    template <typename T, typename Tm>
    class QuotXL : public LowerTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotXL(const T _x, const GenLowerTriMatrix<Tm>& _m) :
            x(_x), m(_m)  {}
        inline ptrdiff_t size() const { return m.size(); }
        inline DiagType dt() const
        { return x==T(1) ? m.dt() : NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenLowerTriMatrix<Tm>& getM() const { return m; }
        inline void assignToL(LowerTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            (m0=m).invertSelf();
            MultXM(x,m0);
        }
        inline void assignToL(LowerTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || m0.dt() == dt());
            (m0=m).invertSelf();
            MultXM(x,m0);
        }
    private:
        const T x;
        const GenLowerTriMatrix<Tm>& m;
    };

#define QUOTXM QuotXL
#define GENMATRIX GenLowerTriMatrix
#define PRODXM ProdXL
#include "tmv/TMV_AuxQuotXM.h"
#include "tmv/TMV_AuxQuotXMa.h"
#undef QUOTXM
#undef GENMATRIX
#undef PRODX

    //
    // TriMatrix / TriMatrix
    //

    template <typename T, typename T1, typename T2>
    class QuotLL : public LowerTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotLL(
            const T _x, const GenLowerTriMatrix<T1>& _m1,
            const GenLowerTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const
        { return x==T(1) && m1.dt()==m2.dt() ? m1.dt() : NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenLowerTriMatrix<T1>& getM1() const { return m1; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToL(LowerTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
            m2.LDiv(m1,m0);
            MultXM(x,m0);
        }
        inline void assignToL(LowerTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
            m2.LDiv(m1,m0);
            MultXM(x,m0);
        }
    protected:
        const T x;
        const GenLowerTriMatrix<T1>& m1;
        const GenLowerTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class RQuotLL : public LowerTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline RQuotLL(
            const T _x, const GenLowerTriMatrix<T1>& _m1,
            const GenLowerTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const
        { return x==T(1) && m1.dt()==m2.dt() ? m1.dt() : NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenLowerTriMatrix<T1>& getM1() const { return m1; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToL(LowerTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
            m2.RDiv(m1,m0);
            MultXM(x,m0);
        }
        inline void assignToL(LowerTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
            m2.RDiv(m1,m0);
            MultXM(x,m0);
        }
    protected:
        const T x;
        const GenLowerTriMatrix<T1>& m1;
        const GenLowerTriMatrix<T2>& m2;
    };

    template <typename T>
    inline LowerTriMatrixView<T> operator/=(
        LowerTriMatrixView<T> m1, const GenLowerTriMatrix<T>& m2)
    { m2.LDivEq(m1); return m1; }

    template <typename T>
    inline LowerTriMatrixView<CT> operator/=(
        LowerTriMatrixView<CT> m1, const GenLowerTriMatrix<T>& m2)
    { m2.LDivEq(m1); return m1; }

    template <typename T>
    inline LowerTriMatrixView<T> operator%=(
        LowerTriMatrixView<T> m1, const GenLowerTriMatrix<T>& m2)
    { m2.RDivEq(m1); return m1; }

    template <typename T>
    inline LowerTriMatrixView<CT> operator%=(
        LowerTriMatrixView<CT> m1, const GenLowerTriMatrix<T>& m2)
    { m2.RDivEq(m1); return m1; }

#define QUOTMM QuotLL
#define RQUOTMM RQuotLL
#define QUOTXM QuotXL
#define GENMATRIX1 GenLowerTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXL
#define PRODXM2 ProdXL
#include "tmv/TMV_AuxQuotMM.h"
#include "tmv/TMV_AuxQuotMMa.h"
#undef QUOTMM
#undef RQUOTMM
#undef QUOTXM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX GenLowerTriMatrix
#define PRODXM ProdXL
#define PRODMV ProdLV
#define PRODVM ProdVL
#define QUOTVM QuotVL
#define QUOTXM QuotXL
#define RQUOTVM RQuotVL
#include "tmv/TMV_AuxVecComposite.h"
#undef PRODMV
#undef PRODVM
#undef QUOTVM
#undef QUOTXM
#undef RQUOTVM
#undef PRODXM
#undef GENMATRIX

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXL
#define SUMMM SumML
#define PRODMM ProdML
#define QUOTMM QuotML
#define QUOTXM QuotXL
#define RQUOTMM RQuotML
#define TQUOTMM TransientQuotML
#define TRQUOTMM TransientRQuotML

#include "tmv/TMV_AuxMatComposite1.h"
#include "tmv/TMV_AuxSumMMb.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef QUOTXM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM

#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL
#define GENMATRIX2 GenMatrix
#define PRODXM2 ProdXM
#define SUMMMa SumML
#define SUMMM SumLM
#define PRODMM ProdLM
#define QUOTXM QuotXM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM

#include "tmv/TMV_AuxMatComposite2.h"
#include "tmv/TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef SUMMMa
#undef PRODMM
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM

    //
    // Finally, the things which mix Upper and Lower
    //

    template <typename T, typename T1, typename T2>
    class SumUL : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumUL(
            T _x1, const GenUpperTriMatrix<T1>& _m1,
            T _x2, const GenLowerTriMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        {  TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline T getX1() const { return x1; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            AddMM(x2,m2,x1,m1,m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            AddMM(x2,m2,x1,m1,m0);
        }
    private:
        T x1;
        const GenUpperTriMatrix<T1>& m1;
        T x2;
        const GenLowerTriMatrix<T2>& m2;
    };

#define SUMMM SumUL
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXL
#include "tmv/TMV_AuxSumMM.h"
#include "tmv/TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    template <typename T, typename T1, typename T2>
    class ProdUL : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdUL(
            T _x, const GenUpperTriMatrix<T1>& _m1,
            const GenLowerTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline T getX() const { return x; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
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

    protected:
        T x;
        const GenUpperTriMatrix<T1>& m1;
        const GenLowerTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ProdUL<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ProdUL<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ProdUL<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ProdUL<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

#define PRODMM ProdUL
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXL
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define TQUOTMM TransientQuotML
#define TRQUOTMM TransientRQuotML
#define QUOTXM QuotXL
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXL
#include "tmv/TMV_AuxTQuotMM.h"
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    template <typename T, typename T1, typename T2>
    class ProdLU : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdLU(
            T _x, const GenLowerTriMatrix<T1>& _m1,
            const GenUpperTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline T getX() const { return x; }
        inline const GenLowerTriMatrix<T1>& getM1() const { return m1; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
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

    protected:
        T x;
        const GenLowerTriMatrix<T1>& m1;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ProdLU<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ProdLU<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ProdLU<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ProdLU<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

#define PRODMM ProdLU
#define GENMATRIX1 GenLowerTriMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXL
#define PRODXM2 ProdXU
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define TQUOTMM TransientQuotMU
#define TRQUOTMM TransientRQuotMU
#define QUOTXM QuotXU
#define GENMATRIX1 GenLowerTriMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXL
#define PRODXM2 ProdXU
#include "tmv/TMV_AuxTQuotMM.h"
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#ifdef TMV_DiagMatrixArith_H
#include "tmv/TMV_DiagTriArith.h"
#endif

#ifdef TMV_BandMatrixArith_H
#include "tmv/TMV_TriBandArith.h"
#endif

#ifdef TMV_SymMatrixArith_H
#include "tmv/TMV_TriSymArith.h"
#endif

#ifdef TMV_SymBandMatrixArith_H
#include "tmv/TMV_TriSymBandArith.h"
#endif

#endif
