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


#ifndef TMV_DiagSymArith_H
#define TMV_DiagSymArith_H

#define CT std::complex<T>

namespace tmv {

    //
    // DiagMatrix + SymMatrix
    //

    template <typename T, typename T1, typename T2>
    class SumDS : public SymMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumDS(
            T _x1, const GenDiagMatrix<T1>& _m1,
            T _x2, const GenSymMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline SymType sym() const { return m2.sym(); }
        inline T getX1() const { return x1; }
        inline const GenDiagMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenSymMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            if (SameStorage(m0,m1.diag())) {
                DiagMatrix<T1> m1x = m1;
                MultXM(x2,m0=m2);
                AddVV(x1,m1x.diag(),m0.diag());
            } else {
                MultXM(x2,m0=m2);
                AddVV(x1,m1.diag(),m0.diag());
            }
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            if (SameStorage(m0,m1.diag())) {
                DiagMatrix<T1> m1x = m1;
                MultXM(x2,m0=m2);
                AddVV(x1,m1.diag(),m0.diag());
            } else {
                MultXM(x2,m0=m2);
                AddVV(x1,m1.diag(),m0.diag());
            }
        }
        inline void assignToS(SymMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            if (SameStorage(m0,m1.diag())) {
                DiagMatrix<T1> m1x = m1;
                MultXM(x2,m0=m2);
                AddVV(x1,m1x.diag(),m0.diag());
            } else {
                MultXM(x2,m0=m2);
                AddVV(x1,m1.diag(),m0.diag());
            }
        }
        inline void assignToS(SymMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == size());
            TMVAssert(m0.rowsize() == size());
            if (!m2.issym()) {
                TMVAssert(isReal(T1()) &&
                          TMV_IMAG(x1) == T1(0) &&
                          TMV_IMAG(x2) == T1(0));
            }
            if (SameStorage(m0,m1.diag())) {
                DiagMatrix<T1> m1x = m1;
                MultXM(x2,m0=m2);
                AddVV(x1,m1x.diag(),m0.diag());
            } else {
                MultXM(x2,m0=m2);
                AddVV(x1,m1.diag(),m0.diag());
            }
        }
    private:
        const T x1;
        const GenDiagMatrix<T1>& m1;
        const T x2;
        const GenSymMatrix<T2>& m2;
    };

    template <typename T>
    inline SymMatrixView<T> operator+=(
        SymMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.issym());
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(
        SymMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline SymMatrixView<T> operator-=(
        SymMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.issym());
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(
        SymMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T, typename T2>
    inline SymMatrixView<T> operator+=(
        SymMatrixView<T> m, const ProdXD<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(m.issym());
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator+=(
        SymMatrixView<CT> m, const ProdXD<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T, typename T2>
    inline SymMatrixView<T> operator-=(
        SymMatrixView<T> m, const ProdXD<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(m.issym());
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline SymMatrixView<CT> operator-=(
        SymMatrixView<CT> m, const ProdXD<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

#define SUMMM SumDS
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXS
#include "tmv/TMV_AuxSumMM.h"
#include "tmv/TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // DiagMatrix * SymMatrix
    //

    template <typename T, typename T1, typename T2>
    class ProdDS : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdDS(
            T _x, const GenDiagMatrix<T1>& _m1,
            const GenSymMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m2.size(); }
        inline ptrdiff_t rowsize() const { return m2.size(); }
        inline T getX() const { return x; }
        inline const GenDiagMatrix<T1>& getM1() const { return m1; }
        inline const GenSymMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m1.diag())) {
                DiagMatrix<T1> m1x = m1;
                MultXM(x,m0=m2);
                MultMM<false>(T(1),m1x,m0,m0);
            } else {
                MultXM(x,m0=m2);
                MultMM<false>(T(1),m1,m0,m0);
            }
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m1.diag())) {
                DiagMatrix<T1> m1x = m1;
                MultXM(x,m0=m2);
                MultMM<false>(T(1),m1x,m0,m0);
            } else {
                MultXM(x,m0=m2);
                MultMM<false>(T(1),m1,m0,m0);
            }
        }
    protected:
        T x;
        const GenDiagMatrix<T1>& m1;
        const GenSymMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class ProdSD : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdSD(
            T _x, const GenSymMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline T getX() const { return x; }
        inline const GenSymMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m2.diag())) {
                DiagMatrix<T2> m2x = m2;
                MultXM(x,m0=m1);
                MultMM<false>(T(1),m0,m2x,m0);
            } else {
                MultXM(x,m0=m1);
                MultMM<false>(T(1),m0,m2,m0);
            }
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m2.diag())) {
                DiagMatrix<T2> m2x = m2;
                MultXM(x,m0=m1);
                MultMM<false>(T(1),m0,m2x,m0);
            } else {
                MultXM(x,m0=m1);
                MultMM<false>(T(1),m0,m2,m0);
            }
        }
    protected:
        T x;
        const GenSymMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

#define PRODMM ProdSD
#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXS
#define PRODXM2 ProdXD
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM ProdDS
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXS
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // SymMatrix / % DiagMatrix
    //

    template <typename T, typename T1, typename T2>
    class QuotSD : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotSD(
            const T _x, const GenSymMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline T getX() const { return x; }
        inline const GenSymMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            MultXM(x,m0=m1);
            m2.LDivEq(m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            MultXM(x,m0=m1);
            m2.LDivEq(m0);
        }
    protected:
        const T x;
        const GenSymMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class RQuotSD : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline RQuotSD(
            const T _x, const GenSymMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline T getX() const { return x; }
        inline const GenSymMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            MultXM(x,m0=m1);
            m2.RDivEq(m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            MultXM(x,m0=m1);
            m2.RDivEq(m0);
        }
    protected:
        const T x;
        const GenSymMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

#define QUOTMM QuotSD
#define QUOTXM QuotXD
#define RQUOTMM RQuotSD
#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXS
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
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXD
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

} // namespace tmv

#undef CT

#endif
