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


#ifndef TMV_DiagBandArith_H
#define TMV_DiagBandArith_H

#define CT std::complex<T>

namespace tmv {

    //
    // DiagMatrix + BandMatrix
    //

    template <typename T, typename T1, typename T2>
    class SumDB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumDB(
            T _x1, const GenDiagMatrix<T1>& _m1,
            T _x2, const GenBandMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        {
            TMVAssert(m1.size() == m2.colsize());
            TMVAssert(m1.size() == m2.rowsize());
        }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline ptrdiff_t nlo() const { return m2.nlo(); }
        inline ptrdiff_t nhi() const { return m2.nhi(); }
        inline T getX1() const { return x1; }
        inline const GenDiagMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            AddMM(x1,BandMatrixViewOf(m1),x2,m2,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            AddMM(x1,BandMatrixViewOf(m1),x2,m2,m0);
        }
    private:
        const T x1;
        const GenDiagMatrix<T1>& m1;
        const T x2;
        const GenBandMatrix<T2>& m2;
    };

    template <typename T>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.rowsize() == m2.size());
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.rowsize() == m2.size());
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.rowsize() == m2.size());
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.rowsize() == m2.size());
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T, typename T2>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, const ProdXD<T,T2>& pxm)
    {
        TMVAssert(m.colsize() == pxm.size());
        TMVAssert(m.rowsize() == pxm.size());
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, const ProdXD<T,T>& pxm)
    {
        TMVAssert(m.colsize() == pxm.size());
        TMVAssert(m.rowsize() == pxm.size());
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T, typename T2>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, const ProdXD<T,T2>& pxm)
    {
        TMVAssert(m.colsize() == pxm.size());
        TMVAssert(m.rowsize() == pxm.size());
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, const ProdXD<T,T>& pxm)
    {
        TMVAssert(m.colsize() == pxm.size());
        TMVAssert(m.rowsize() == pxm.size());
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline DiagMatrixView<T> operator+=(
        DiagMatrixView<T> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nlo() == 0 && m2.nhi() == 0);
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator+=(
        DiagMatrixView<CT> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nlo() == 0 && m2.nhi() == 0);
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline DiagMatrixView<T> operator-=(
        DiagMatrixView<T> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nlo() == 0 && m2.nhi() == 0);
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator-=(
        DiagMatrixView<CT> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nlo() == 0 && m2.nhi() == 0);
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T, typename T2>
    inline DiagMatrixView<T> operator+=(
        DiagMatrixView<T> m, const ProdXB<T,T2>& pxm)
    {
        TMVAssert(pxm.colsize() == m.size());
        TMVAssert(pxm.rowsize() == m.size());
        TMVAssert(pxm.nlo() == 0 && pxm.nhi() == 0);
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator+=(
        DiagMatrixView<CT> m, const ProdXB<T,T>& pxm)
    {
        TMVAssert(pxm.colsize() == m.size());
        TMVAssert(pxm.rowsize() == m.size());
        TMVAssert(pxm.nlo() == 0 && pxm.nhi() == 0);
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T, typename T2>
    inline DiagMatrixView<T> operator-=(
        DiagMatrixView<T> m, const ProdXB<T,T2>& pxm)
    {
        TMVAssert(pxm.colsize() == m.size());
        TMVAssert(pxm.rowsize() == m.size());
        TMVAssert(pxm.nlo() == 0 && pxm.nhi() == 0);
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator-=(
        DiagMatrixView<CT> m, const ProdXB<T,T>& pxm)
    {
        TMVAssert(pxm.colsize() == m.size());
        TMVAssert(pxm.rowsize() == m.size());
        TMVAssert(pxm.nlo() == 0 && pxm.nhi() == 0);
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

#define SUMMM SumDB
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXB
#include "tmv/TMV_AuxSumMM.h"
#include "tmv/TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // DiagMatrix * BandMatrix
    //

    template <typename T, typename T1, typename T2>
    class ProdDB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdDB(
            T _x, const GenDiagMatrix<T1>& _m1,
            const GenBandMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.colsize()); }
        inline ptrdiff_t colsize() const { return m2.colsize(); }
        inline ptrdiff_t rowsize() const { return m2.rowsize(); }
        inline ptrdiff_t nlo() const { return m2.nlo(); }
        inline ptrdiff_t nhi() const { return m2.nhi(); }
        inline T getX() const { return x; }
        inline const GenDiagMatrix<T1>& getM1() const { return m1; }
        inline const GenBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,BandMatrixViewOf(m1),m2,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,BandMatrixViewOf(m1),m2,m0);
        }
    protected:
        T x;
        const GenDiagMatrix<T1>& m1;
        const GenBandMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class ProdBD : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdBD(
            T _x, const GenBandMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.rowsize() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline ptrdiff_t nlo() const { return m1.nlo(); }
        inline ptrdiff_t nhi() const { return m1.nhi(); }
        inline T getX() const { return x; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,m1,BandMatrixViewOf(m2),m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,m1,BandMatrixViewOf(m2),m0);
        }
    protected:
        T x;
        const GenBandMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

    template <typename T>
    inline BandMatrixView<T> operator*=(
        BandMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.rowsize() == m2.size());
        MultMM<false>(T(1),m1,BandMatrixViewOf(m2),m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator*=(
        BandMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.rowsize() == m2.size());
        MultMM<false>(T(1),m1,BandMatrixViewOf(m2),m1);
        return m1;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, const ProdDB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(
            pmm.getX(),BandMatrixViewOf(pmm.getM1()),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, const ProdDB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(
            pmm.getX(),BandMatrixViewOf(pmm.getM1()),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, const ProdDB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(
            -pmm.getX(),BandMatrixViewOf(pmm.getM1()),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, const ProdDB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(
            -pmm.getX(),BandMatrixViewOf(pmm.getM1()),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, const ProdBD<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(
            pmm.getX(),pmm.getM1(),BandMatrixViewOf(pmm.getM2()),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, const ProdBD<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(
            pmm.getX(),pmm.getM1(),BandMatrixViewOf(pmm.getM2()),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, const ProdBD<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(
            -pmm.getX(),pmm.getM1(),BandMatrixViewOf(pmm.getM2()),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, const ProdBD<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(
            -pmm.getX(),pmm.getM1(),BandMatrixViewOf(pmm.getM2()),m);
        return m;
    }

    template <typename T>
    inline DiagMatrixView<T> operator*=(
        DiagMatrixView<T> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.nlo() == 0);
        TMVAssert(m2.nhi() == 0);
        MultMM<false>(
            T(1),BandMatrixViewOf(m1),m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator*=(
        DiagMatrixView<CT> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.nlo() == 0);
        TMVAssert(m2.nhi() == 0);
        MultMM<false>(
            T(1),BandMatrixViewOf(m1),m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T, typename T1, typename T2>
    inline DiagMatrixView<T> operator+=(
        DiagMatrixView<T> m, const ProdDB<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(pmm.nlo() == 0);
        MultMM<true>(
            pmm.getX(),BandMatrixViewOf(pmm.getM1()),pmm.getM2(),
            BandMatrixViewOf(m));
        return m;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator+=(
        DiagMatrixView<CT> m, const ProdDB<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(pmm.nlo() == 0);
        MultMM<true>(
            T(pmm.getX()),BandMatrixViewOf(pmm.getM1()),pmm.getM2(),
            BandMatrixViewOf(m));
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline DiagMatrixView<T> operator-=(
        DiagMatrixView<T> m, const ProdDB<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(pmm.nlo() == 0);
        MultMM<true>(
            -pmm.getX(),BandMatrixViewOf(pmm.getM1()),pmm.getM2(),
            BandMatrixViewOf(m));
        return m;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator-=(
        DiagMatrixView<CT> m, const ProdDB<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(pmm.nlo() == 0);
        MultMM<true>(
            -pmm.getX(),BandMatrixViewOf(pmm.getM1()),pmm.getM2(),
            BandMatrixViewOf(m));
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline DiagMatrixView<T> operator+=(
        DiagMatrixView<T> m, const ProdBD<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(pmm.nlo() == 0);
        MultMM<true>(
            pmm.getX(),pmm.getM1(),BandMatrixViewOf(pmm.getM2()),
            BandMatrixViewOf(m));
        return m;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator+=(
        DiagMatrixView<CT> m, const ProdBD<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(pmm.nlo() == 0);
        MultMM<true>(
            pmm.getX(),pmm.getM1(),BandMatrixViewOf(pmm.getM2()),
            BandMatrixViewOf(m));
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline DiagMatrixView<T> operator-=(
        DiagMatrixView<T> m, const ProdBD<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(pmm.nlo() == 0);
        MultMM<true>(
            -pmm.getX(),pmm.getM1(),BandMatrixViewOf(pmm.getM2()),
            BandMatrixViewOf(m));
        return m;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator-=(
        DiagMatrixView<CT> m, const ProdBD<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(pmm.nlo() == 0);
        MultMM<true>(
            -pmm.getX(),pmm.getM1(),BandMatrixViewOf(pmm.getM2()),
            BandMatrixViewOf(m));
        return m;
    }

#define PRODMM ProdBD
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXD
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM ProdDB
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXB
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // BandMatrix / % DiagMatrix
    //

    template <typename T, typename T1, typename T2>
    class QuotBD : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotBD(
            const T _x, const GenBandMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.colsize() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline ptrdiff_t nlo() const { return m1.nlo(); }
        inline ptrdiff_t nhi() const { return m1.nhi(); }
        inline T getX() const { return x; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            MultMM<false>(
                x,BandMatrixViewOf(DiagMatrix<T2>(m2.inverse())),m1,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            MultMM<false>(
                x,BandMatrixViewOf(DiagMatrix<T2>(m2.inverse())),m1,m0);
        }
    protected:
        const T x;
        const GenBandMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class RQuotBD : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline RQuotBD(
            const T _x, const GenBandMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.rowsize() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline ptrdiff_t nlo() const { return m1.nlo(); }
        inline ptrdiff_t nhi() const { return m1.nhi(); }
        inline T getX() const { return x; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            MultMM<false>(
                x,m1,BandMatrixViewOf(DiagMatrix<T2>(m2.inverse())),m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            MultMM<false>(
                x,m1,BandMatrixViewOf(DiagMatrix<T2>(m2.inverse())),m0);
        }
    protected:
        const T x;
        const GenBandMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

    template <typename T>
    inline BandMatrixView<T> operator/=(
        BandMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        MultMM<false>(
            T(1),BandMatrixViewOf(DiagMatrix<T>(m2.inverse())),m1,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator/=(
        BandMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        MultMM<false>(
            T(1),BandMatrixViewOf(DiagMatrix<T>(m2.inverse())),m1,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<T> operator%=(
        BandMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.rowsize() == m2.size());
        MultMM<false>(
            T(1),m1,BandMatrixViewOf(DiagMatrix<T>(m2.inverse())),m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator%=(
        BandMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.rowsize() == m2.size());
        MultMM<false>(
            T(1),m1,BandMatrixViewOf(DiagMatrix<T>(m2.inverse())),m1);
        return m1;
    }

#define QUOTMM QuotBD
#define QUOTXM QuotXD
#define RQUOTMM RQuotBD
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXD
#include "tmv/TMV_AuxQuotMM.h"
#include "tmv/TMV_AuxQuotMMa.h"
#undef QUOTMM
#undef QUOTXM
#undef RQUOTMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXB
#define QUOTXM QuotXB
#define TQUOTMM TransientQuotMB
#define TRQUOTMM TransientRQuotMB

#include "tmv/TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM

} // namespace tmv

#undef CT

#endif
