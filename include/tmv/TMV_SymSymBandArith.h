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


#ifndef TMV_SymSymBandArith_H
#define TMV_SymSymBandArith_H

#define CT std::complex<T>

namespace tmv {

    //
    // SymBandMatrix + SymMatrix
    //

    template <typename T, typename T1, typename T2>
    class SumsBS : public SymMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumsBS(
            T _x1, const GenSymBandMatrix<T1>& _m1,
            T _x2, const GenSymMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m2.size(); }
        inline SymType sym() const
        { return isReal(T1()) ? m2.sym() : m1.sym(); }
        inline T getX1() const { return x1; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenSymMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            if (SameStorage(m0,m1)) {
                SymBandMatrix<T1> m1x = m1;
                MultXM(x2,m0=m2);
                AddMM(x1,m1x,m0);
            } else {
                MultXM(x2,m0=m2);
                AddMM(x1,m1,m0);
            }
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            if (SameStorage(m0,m1)) {
                if (m1.issym()) {
                    SymBandMatrix<T1> m1x = m1;
                    MultXM(x2,m0=m2);
                    AddMM(x1,m1x,m0);
                } else {
                    HermBandMatrix<T1> m1x = m1;
                    MultXM(x2,m0=m2);
                    AddMM(x1,m1x,m0);
                }
            } else {
                MultXM(x2,m0=m2);
                AddMM(x1,m1,m0);
            }
        }
        inline void assignToS(SymMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(isReal(T1()) || m0.issym() == m1.issym());
            TMVAssert(isReal(T2()) || m0.issym() == m2.issym());
            TMVAssert(m0.issym() || TMV_IMAG(x1) == real_type(0));
            TMVAssert(m0.issym() || TMV_IMAG(x2) == real_type(0));
            if (SameStorage(m0,m1)) {
                SymBandMatrix<T1> m1x = m1;
                MultXM(x2,m0=m2);
                AddMM(x1,m1x,SymBandMatrixViewOf(m0,m1.nlo()));
            } else {
                MultXM(x2,m0=m2);
                AddMM(x1,m1,SymBandMatrixViewOf(m0,m1.nlo()));
            }
        }
        inline void assignToS(SymMatrixView<complex_type> m0) const
        {
            TMVAssert(isReal(T1()) || isReal(T2()) || m1.sym() == m2.sym());
            TMVAssert(m0.size() == size());
            TMVAssert(isReal(T1()) || m0.issym() == m1.issym());
            TMVAssert(isReal(T2()) || m0.issym() == m2.issym());
            TMVAssert(m0.issym() || TMV_IMAG(x1) == real_type(0));
            TMVAssert(m0.issym() || TMV_IMAG(x2) == real_type(0));
            if (SameStorage(m0,m1)) {
                if (m1.issym()) {
                    SymBandMatrix<T1> m1x = m1;
                    MultXM(x2,m0=m2);
                    AddMM(x1,m1x,SymBandMatrixViewOf(m0,m1.nlo()));
                } else {
                    HermBandMatrix<T1> m1x = m1;
                    MultXM(x2,m0=m2);
                    AddMM(x1,m1x,SymBandMatrixViewOf(m0,m1.nlo()));
                }
            } else {
                MultXM(x2,m0=m2);
                AddMM(x1,m1,SymBandMatrixViewOf(m0,m1.nlo()));
            }
        }

    private:
        const T x1;
        const GenSymBandMatrix<T1>& m1;
        const T x2;
        const GenSymMatrix<T2>& m2;
    };

    template <typename T>
    inline SymMatrixView<T> operator+=(
        SymMatrixView<T> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(isReal(T()) || m1.sym() == m2.sym());
        AddMM(T(1),m2,SymBandMatrixViewOf(m1,m2.nlo()));
        return m1;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(
        SymMatrixView<CT> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        AddMM(T(1),m2,SymBandMatrixViewOf(m1,m2.nlo()));
        return m1;
    }

    template <typename T>
    inline SymMatrixView<T> operator-=(
        SymMatrixView<T> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(isReal(T()) || m1.sym() == m2.sym());
        AddMM(T(-1),m2,SymBandMatrixViewOf(m1,m2.nlo()));
        return m1;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(
        SymMatrixView<CT> m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        AddMM(T(-1),m2,SymBandMatrixViewOf(m1,m2.nlo()));
        return m1;
    }

    template <typename T, typename T2>
    inline SymMatrixView<T> operator+=(
        SymMatrixView<T> m, const ProdXsB<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(isReal(T2()) || m.issym() == pxm.getM().issym());
        TMVAssert(m.issym() || TMV_IMAG(pxm.getX()) == TMV_RealType(T)(0));
        AddMM(pxm.getX(),pxm.getM(),SymBandMatrixViewOf(m,pxm.nlo()));
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(
        SymMatrixView<CT> m, const ProdXsB<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        AddMM(pxm.getX(),pxm.getM(),SymBandMatrixViewOf(m,pxm.nlo()));
        return m;
    }

    template <typename T, typename T2>
    inline SymMatrixView<T> operator-=(
        SymMatrixView<T> m, const ProdXsB<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(isReal(T2()) || m.issym() == pxm.getM().issym());
        TMVAssert(m.issym() || TMV_IMAG(pxm.getX()) == TMV_RealType(T)(0));
        AddMM(-pxm.getX(),pxm.getM(),SymBandMatrixViewOf(m,pxm.nlo()));
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(
        SymMatrixView<CT> m, const ProdXsB<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        AddMM(-pxm.getX(),pxm.getM(),SymBandMatrixViewOf(m,pxm.nlo()));
        return m;
    }


#define SUMMM SumsBS
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXS
#include "tmv/TMV_AuxSumMM.h"
#include "tmv/TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // SymBandMatrix * SymMatrix
    //

    template <typename T, typename T1, typename T2>
    class ProdsBS : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdsBS(
            T _x, const GenSymBandMatrix<T1>& _m1,
            const GenSymMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m2.size(); }
        inline ptrdiff_t rowsize() const { return m2.size(); }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline const GenSymMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            Matrix<T> m2x = m2;
            MultMM<false>(x,m1,m2x,m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            Matrix<T> m2x = m2;
            MultMM<false>(x,m1,m2x,m0);
        }
    protected:
        T x;
        const GenSymBandMatrix<T1>& m1;
        const GenSymMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class ProdSsB : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdSsB(
            T _x, const GenSymMatrix<T1>& _m1,
            const GenSymBandMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        {
            TMVAssert(m2.colsize() == m1.size());
            TMVAssert(m2.rowsize() == m1.size());
        }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline T getX() const { return x; }
        inline const GenSymMatrix<T1>& getM1() const { return m1; }
        inline const GenSymBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            Matrix<T> m1x = m1;
            MultMM<false>(x,m1x,m2,m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            Matrix<T> m1x = m1;
            MultMM<false>(x,m1x,m2,m0);
        }
    protected:
        T x;
        const GenSymMatrix<T1>& m1;
        const GenSymBandMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ProdSsB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        AddMM(T(1),Matrix<T>(pmm),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ProdSsB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        AddMM(T(1),Matrix<T>(pmm),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ProdSsB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        AddMM(T(-1),Matrix<T>(pmm),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ProdSsB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        AddMM(T(-1),Matrix<T>(pmm),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ProdsBS<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        AddMM(T(1),Matrix<T>(pmm),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ProdsBS<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        AddMM(T(1),Matrix<T>(pmm),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ProdsBS<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        AddMM(T(-1),Matrix<T>(pmm),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ProdsBS<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        AddMM(T(-1),Matrix<T>(pmm),m);
        return m;
    }


#define PRODMM ProdSsB
#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXS
#define PRODXM2 ProdXsB
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM ProdsBS
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXS
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // SymMatrix / % SymBandMatrix
    //

#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXS
#define QUOTXM QuotXS
#define TQUOTMM TransientQuotMS
#define TRQUOTMM TransientRQuotMS
#include "tmv/TMV_AuxTQuotMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM

#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXS
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
