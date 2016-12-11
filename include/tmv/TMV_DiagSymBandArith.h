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


#ifndef TMV_DiagSymBandArith_H
#define TMV_DiagSymBandArith_H

#define CT std::complex<T>

namespace tmv {

    //
    // DiagMatrix + SymBandMatrix
    //

    template <typename T, typename T1, typename T2>
    class SumDsB : public SymBandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumDsB(
            T _x1, const GenDiagMatrix<T1>& _m1,
            T _x2, const GenSymBandMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        {
            TMVAssert(m1.size() == m2.colsize());
            TMVAssert(m1.size() == m2.rowsize());
        }
        inline ptrdiff_t size() const { return m1.size(); }
        inline ptrdiff_t nlo() const { return m2.nlo(); }
        inline SymType sym() const { return m2.sym(); }
        inline T getX1() const { return x1; }
        inline const GenDiagMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenSymBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(m0.nhi() >= nlo());
            TMVAssert(m0.nlo() >= nlo());
            AddMM(
                x1,BandMatrixViewOf(m1),x2,m2.upperBand(),
                m0.diagRange(0,m0.nhi()+1));
            if (m2.nlo() > 0)
                MultXM(x2,m0.diagRange(-m0.nlo(),0)=m2.lowerBandOff());
            else if (m0.nlo() > 0)
                m0.diagRange(-m0.nlo(),0).setZero();
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(m0.nhi() >= nlo());
            TMVAssert(m0.nlo() >= nlo());
            AddMM(
                x1,BandMatrixViewOf(m1),x2,m2.upperBand(),
                m0.diagRange(0,m0.nhi()+1));
            if (m2.nlo() > 0)
                MultXM(x2,m0.diagRange(-m0.nlo(),0)=m2.lowerBandOff());
            else if (m0.nlo() > 0)
                m0.diagRange(-m0.nlo(),0).setZero();
        }
        inline void assignTosB(SymBandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(m0.nlo() >= nlo());
            AddMM(
                x1,BandMatrixViewOf(m1),x2,m2.upperBand(),
                  m0.diagRange(0,m0.nhi()+1));
        }
        inline void assignTosB(SymBandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(m0.nlo() >= nlo());
            if (!m2.issym()) {
                TMVAssert(isReal(T1()) &&
                          TMV_IMAG(x1) == T1(0) &&
                          TMV_IMAG(x2) == T1(0));
            }
            AddMM(
                x1,BandMatrixViewOf(m1),x2,m2.upperBand(),
                  m0.diagRange(0,m0.nhi()+1));
        }
    private:
        const T x1;
        const GenDiagMatrix<T1>& m1;
        const T x2;
        const GenSymBandMatrix<T2>& m2;
    };

    template <typename T>
    inline SymBandMatrixView<T> operator+=(
        SymBandMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.issym());
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator+=(
        SymBandMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline SymBandMatrixView<T> operator-=(
        SymBandMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.issym());
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator-=(
        SymBandMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T, typename T2>
    inline SymBandMatrixView<T> operator+=(
        SymBandMatrixView<T> m, const ProdXD<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(m.issym());
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator+=(
        SymBandMatrixView<CT> m, const ProdXD<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T, typename T2>
    inline SymBandMatrixView<T> operator-=(
        SymBandMatrixView<T> m, const ProdXD<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(m.issym());
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline SymBandMatrixView<CT> operator-=(
        SymBandMatrixView<CT> m, const ProdXD<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline DiagMatrixView<T> operator+=(
        DiagMatrixView<T> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m2.nlo() == 0);
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator+=(
        DiagMatrixView<CT> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m2.nlo() == 0);
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline DiagMatrixView<T> operator-=(
        DiagMatrixView<T> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m2.nlo() == 0);
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator-=(
        DiagMatrixView<CT> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m2.nlo() == 0);
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T, typename T2>
    inline DiagMatrixView<T> operator+=(
        DiagMatrixView<T> m, const ProdXsB<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(pxm.nlo() == 0);
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator+=(
        DiagMatrixView<CT> m, const ProdXsB<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(pxm.nlo() == 0);
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T, typename T2>
    inline DiagMatrixView<T> operator-=(
        DiagMatrixView<T> m, const ProdXsB<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(pxm.nlo() == 0);
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator-=(
        DiagMatrixView<CT> m, const ProdXsB<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(pxm.nlo() == 0);
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

#define SUMMM SumDsB
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXsB
#include "tmv/TMV_AuxSumMM.h"
#include "tmv/TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // DiagMatrix * SymBandMatrix
    //

    template <typename T, typename T1, typename T2>
    class ProdDsB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdDsB(
            T _x, const GenDiagMatrix<T1>& _m1,
            const GenSymBandMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline ptrdiff_t nlo() const { return m2.nlo(); }
        inline ptrdiff_t nhi() const { return m2.nlo(); }
        inline T getX() const { return x; }
        inline const GenDiagMatrix<T1>& getM1() const { return m1; }
        inline const GenSymBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            if (SameStorage(m0,m1.diag())) {
                DiagMatrix<T1> m1x = m1;
                MultXM(x,m0=m2);
                MultMM<false>(T(1),BandMatrixViewOf(m1x),m0,m0);
            } else {
                MultXM(x,m0=m2);
                MultMM<false>(T(1),BandMatrixViewOf(m1),m0,m0);
            }
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            if (SameStorage(m0,m1.diag())) {
                DiagMatrix<T1> m1x = m1;
                MultXM(x,m0=m2);
                MultMM<false>(T(1),BandMatrixViewOf(m1x),m0,m0);
            } else {
                MultXM(x,m0=m2);
                MultMM<false>(T(1),BandMatrixViewOf(m1),m0,m0);
            }
        }
    protected:
        T x;
        const GenDiagMatrix<T1>& m1;
        const GenSymBandMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class ProdsBD : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdsBD(
            T _x, const GenSymBandMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline ptrdiff_t nlo() const { return m1.nlo(); }
        inline ptrdiff_t nhi() const { return m1.nlo(); }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            if (SameStorage(m0,m2.diag())) {
                DiagMatrix<T2> m2x = m2;
                MultXM(x,m0=m1);
                MultMM<false>(T(1),m0,BandMatrixViewOf(m2x),m0);
            } else {
                MultXM(x,m0=m1);
                MultMM<false>(T(1),m0,BandMatrixViewOf(m2),m0);
            }
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            if (SameStorage(m0,m2.diag())) {
                DiagMatrix<T2> m2x = m2;
                MultXM(x,m0=m1);
                MultMM<false>(T(1),m0,BandMatrixViewOf(m2x),m0);
            } else {
                MultXM(x,m0=m1);
                MultMM<false>(T(1),m0,BandMatrixViewOf(m2),m0);
            }
        }
    protected:
        T x;
        const GenSymBandMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

    template <typename T>
    inline DiagMatrixView<T> operator*=(
        DiagMatrixView<T> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m2.nlo() == 0);
        MultMM<false>(T(1),m1,DiagMatrixViewOf(m2.diag()),m1);
        return m1;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator*=(
        DiagMatrixView<CT> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m2.nlo() == 0);
        MultMM<false>(T(1),m1,DiagMatrixViewOf(m2.diag()),m1);
        return m1;
    }

    template <typename T, typename T2>
    inline DiagMatrixView<T> operator*=(
        DiagMatrixView<T> m, const ProdXsB<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(pxm.nlo() == 0);
        MultMM<false>(pxm.getX(),m,DiagMatrixViewOf(pxm.getM().diag()),m);
        return m;
    }

    template <typename T>
    inline DiagMatrixView<CT> operator*=(
        DiagMatrixView<CT> m, const ProdXsB<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(pxm.nlo() == 0);
        MultMM<false>(pxm.getX(),m,DiagMatrixViewOf(pxm.getM().diag()),m);
        return m;
    }

#define PRODMM ProdsBD
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXD
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM ProdDsB
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXsB
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // SymBandMatrix / % DiagMatrix
    //

    template <typename T, typename T1, typename T2>
    class QuotsBD : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotsBD(
            const T _x, const GenSymBandMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline ptrdiff_t nlo() const { return m1.nlo(); }
        inline ptrdiff_t nhi() const { return m1.nlo(); }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == m1.size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            MultXM(x,m0=m1);
            MultMM<false>(
                T(1),BandMatrixViewOf(DiagMatrix<T2>(m2.inverse())),m0,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == m1.size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            MultXM(x,m0=m1);
            MultMM<false>(
                T(1),BandMatrixViewOf(DiagMatrix<T2>(m2.inverse())),m0,m0);
        }
    protected:
        const T x;
        const GenSymBandMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class RQuotsBD : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline RQuotsBD(
            const T _x, const GenSymBandMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline ptrdiff_t nlo() const { return m1.nlo(); }
        inline ptrdiff_t nhi() const { return m1.nlo(); }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == m1.size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nlo());
            MultXM(x,m0=m1);
            MultMM<false>(
                T(1),m0,BandMatrixViewOf(DiagMatrix<T2>(m2.inverse())),m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == m1.size());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nlo());
            MultXM(x,m0=m1);
            MultMM<false>(
                T(1),m0,BandMatrixViewOf(DiagMatrix<T2>(m2.inverse())),m0);
        }
    protected:
        const T x;
        const GenSymBandMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

#define QUOTMM QuotsBD
#define QUOTXM QuotXD
#define RQUOTMM RQuotsBD
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXsB
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
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXsB
#define QUOTXM QuotXsB
#define TQUOTMM TransientQuotMsB
#define TRQUOTMM TransientRQuotMsB
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
