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


#ifndef TMV_TriBandArith_H
#define TMV_TriBandArith_H

#include "tmv/TMV_TriBandArithFunc.h"

#define CT std::complex<T>

namespace tmv {

    //
    // TriMatrix + BandMatrix
    //

    template <typename T, typename T1, typename T2>
    class SumUB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumUB(
            T _x1, const GenUpperTriMatrix<T1>& _m1,
            T _x2, const GenBandMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        {
            TMVAssert(m1.size() == m2.colsize());
            TMVAssert(m1.size() == m2.rowsize());
        }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline ptrdiff_t nlo() const { return m2.nlo(); }
        inline ptrdiff_t nhi() const { return m1.size()-1; }
        inline T getX1() const { return x1; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            if (SameStorage(m0,m2)) {
                if (SameStorage(m0,m1)) {
                    BandMatrix<T2> m2x = m2;
                    MultXM(x1,m0=m1);
                    AddMM(x2,m2x,m0);
                } else {
                    MultXM(x2,m0=m2);
                    AddMM(x1,m1,m0);
                }
            } else {
                MultXM(x1,m0=m1);
                AddMM(x2,m2,m0);
            }
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            if (SameStorage(m0,m2)) {
                if (SameStorage(m0,m1)) {
                    BandMatrix<T2> m2x = m2;
                    MultXM(x1,m0=m1);
                    AddMM(x2,m2x,m0);
                } else {
                    MultXM(x2,m0=m2);
                    AddMM(x1,m1,m0);
                }
            } else {
                MultXM(x1,m0=m1);
                AddMM(x2,m2,m0);
            }
        }
    private:
        const T x1;
        const GenUpperTriMatrix<T1>& m1;
        const T x2;
        const GenBandMatrix<T2>& m2;
    };

    template <typename T>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nhi() == m2.size()-1);
        AddMM(T(1),m2,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nhi() == m2.size()-1);
        AddMM(T(1),m2,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nhi() == m2.size()-1);
        AddMM(T(-1),m2,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nhi() == m2.size()-1);
        AddMM(T(-1),m2,m1);
        return m1;
    }

    template <typename T, typename T2>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, const ProdXU<T,T2>& pxm)
    {
        TMVAssert(m.colsize() == pxm.size());
        TMVAssert(m.rowsize() == pxm.size());
        TMVAssert(m.nhi() == pxm.size()-1);
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, const ProdXU<T,T>& pxm)
    {
        TMVAssert(m.colsize() == pxm.size());
        TMVAssert(m.rowsize() == pxm.size());
        TMVAssert(m.nhi() == pxm.size()-1);
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T, typename T2>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, const ProdXU<T,T2>& pxm)
    {
        TMVAssert(m.colsize() == pxm.size());
        TMVAssert(m.rowsize() == pxm.size());
        TMVAssert(m.nhi() == pxm.size()-1);
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, const ProdXU<T,T>& pxm)
    {
        TMVAssert(m.colsize() == pxm.size());
        TMVAssert(m.rowsize() == pxm.size());
        TMVAssert(m.nhi() == pxm.size()-1);
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<T> operator+=(
        UpperTriMatrixView<T> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nlo() == 0);
        TMVAssert(!m1.isunit());
        AddMM(T(1),m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(
        UpperTriMatrixView<CT> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nlo() == 0);
        TMVAssert(!m1.isunit());
        AddMM(T(1),m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T>
    inline UpperTriMatrixView<T> operator-=(
        UpperTriMatrixView<T> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nlo() == 0);
        TMVAssert(!m1.isunit());
        AddMM(T(-1),m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(
        UpperTriMatrixView<CT> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nlo() == 0);
        TMVAssert(!m1.isunit());
        AddMM(T(-1),m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T, typename T2>
    inline UpperTriMatrixView<T> operator+=(
        UpperTriMatrixView<T> m, const ProdXB<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.colsize());
        TMVAssert(m.size() == pxm.rowsize());
        TMVAssert(pxm.nlo() == 0);
        TMVAssert(!m.isunit());
        BandMatrixViewOf(m) += pxm;
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(
        UpperTriMatrixView<CT> m, const ProdXB<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.colsize());
        TMVAssert(m.size() == pxm.rowsize());
        TMVAssert(pxm.nlo() == 0);
        TMVAssert(!m.isunit());
        BandMatrixViewOf(m) += pxm;
        return m;
    }

    template <typename T, typename T2>
    inline UpperTriMatrixView<T> operator-=(
        UpperTriMatrixView<T> m, const ProdXB<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.colsize());
        TMVAssert(m.size() == pxm.rowsize());
        TMVAssert(pxm.nlo() == 0);
        TMVAssert(!m.isunit());
        BandMatrixViewOf(m) -= pxm;
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(
        UpperTriMatrixView<CT> m, const ProdXB<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.colsize());
        TMVAssert(m.size() == pxm.rowsize());
        TMVAssert(pxm.nlo() == 0);
        TMVAssert(!m.isunit());
        BandMatrixViewOf(m) -= pxm;
        return m;
    }

    template <typename T, typename T1, typename T2>
    class SumLB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumLB(
            T _x1, const GenLowerTriMatrix<T1>& _m1,
            T _x2, const GenBandMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        {
            TMVAssert(m1.size() == m2.colsize());
            TMVAssert(m1.size() == m2.rowsize());
        }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline ptrdiff_t nlo() const { return m1.size()-1; }
        inline ptrdiff_t nhi() const { return m2.nhi(); }
        inline T getX1() const { return x1; }
        inline const GenLowerTriMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            if (SameStorage(m0,m2)) {
                if (SameStorage(m0,m1)) {
                    BandMatrix<T2> m2x = m2;
                    MultXM(x1,m0=m1);
                    AddMM(x2,m2x,m0);
                } else {
                    MultXM(x2,m0=m2);
                    AddMM(x1,m1,m0);
                }
            } else {
                MultXM(x1,m0=m1);
                AddMM(x2,m2,m0);
            }
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            if (SameStorage(m0,m2)) {
                if (SameStorage(m0,m1)) {
                    BandMatrix<T2> m2x = m2;
                    MultXM(x1,m0=m1);
                    AddMM(x2,m2x,m0);
                } else {
                    MultXM(x2,m0=m2);
                    AddMM(x1,m1,m0);
                }
            } else {
                MultXM(x1,m0=m1);
                AddMM(x2,m2,m0);
            }
        }
    private:
        const T x1;
        const GenLowerTriMatrix<T1>& m1;
        const T x2;
        const GenBandMatrix<T2>& m2;
    };

    template <typename T>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nlo() == m2.size()-1);
        AddMM(T(1),m2,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nlo() == m2.size()-1);
        AddMM(T(1),m2,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nlo() == m2.size()-1);
        AddMM(T(-1),m2,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nlo() == m2.size()-1);
        AddMM(T(-1),m2,m1);
        return m1;
    }

    template <typename T, typename T2>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, const ProdXL<T,T2>& pxm)
    {
        TMVAssert(m.colsize() == pxm.size());
        TMVAssert(m.rowsize() == pxm.size());
        TMVAssert(m.nlo() == pxm.size()-1);
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, const ProdXL<T,T>& pxm)
    {
        TMVAssert(m.colsize() == pxm.size());
        TMVAssert(m.rowsize() == pxm.size());
        TMVAssert(m.nlo() == pxm.size()-1);
        AddMM(pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T, typename T2>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, const ProdXL<T,T2>& pxm)
    {
        TMVAssert(m.colsize() == pxm.size());
        TMVAssert(m.rowsize() == pxm.size());
        TMVAssert(m.nlo() == pxm.size()-1);
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, const ProdXL<T,T>& pxm)
    {
        TMVAssert(m.colsize() == pxm.size());
        TMVAssert(m.rowsize() == pxm.size());
        TMVAssert(m.nlo() == pxm.size()-1);
        AddMM(-pxm.getX(),pxm.getM(),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<T> operator+=(
        LowerTriMatrixView<T> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nhi() == 0);
        TMVAssert(!m1.isunit());
        AddMM(T(1),m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(
        LowerTriMatrixView<CT> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nhi() == 0);
        TMVAssert(!m1.isunit());
        AddMM(T(1),m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T>
    inline LowerTriMatrixView<T> operator-=(
        LowerTriMatrixView<T> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nhi() == 0);
        TMVAssert(!m1.isunit());
        AddMM(T(-1),m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(
        LowerTriMatrixView<CT> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nhi() == 0);
        TMVAssert(!m1.isunit());
        AddMM(T(-1),m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T, typename T2>
    inline LowerTriMatrixView<T> operator+=(
        LowerTriMatrixView<T> m, const ProdXB<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.colsize());
        TMVAssert(m.size() == pxm.rowsize());
        TMVAssert(pxm.nhi() == 0);
        TMVAssert(!m.isunit());
        BandMatrixViewOf(m) += pxm;
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(
        LowerTriMatrixView<CT> m, const ProdXB<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.colsize());
        TMVAssert(m.size() == pxm.rowsize());
        TMVAssert(pxm.nhi() == 0);
        TMVAssert(!m.isunit());
        BandMatrixViewOf(m) += pxm;
        return m;
    }

    template <typename T, typename T2>
    inline LowerTriMatrixView<T> operator-=(
        LowerTriMatrixView<T> m, const ProdXB<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.colsize());
        TMVAssert(m.size() == pxm.rowsize());
        TMVAssert(pxm.nhi() == 0);
        TMVAssert(!m.isunit());
        BandMatrixViewOf(m) -= pxm;
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(
        LowerTriMatrixView<CT> m, const ProdXB<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.colsize());
        TMVAssert(m.size() == pxm.rowsize());
        TMVAssert(pxm.nhi() == 0);
        TMVAssert(!m.isunit());
        BandMatrixViewOf(m) -= pxm;
        return m;
    }

#define SUMMM SumUB
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXB
#include "tmv/TMV_AuxSumMM.h"
#include "tmv/TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef PRODXM1
#define SUMMM SumLB
#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL
#include "tmv/TMV_AuxSumMM.h"
#include "tmv/TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // TriMatrix * BandMatrix
    //

    template <typename T, typename T1, typename T2>
    class ProdUB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdUB(
            T _x, const GenUpperTriMatrix<T1>& _m1,
            const GenBandMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.colsize()); }
        inline ptrdiff_t colsize() const { return m2.colsize(); }
        inline ptrdiff_t rowsize() const { return m2.rowsize(); }
        inline ptrdiff_t nlo() const { return m2.nlo(); }
        inline ptrdiff_t nhi() const
        { return TMV_MIN(rowsize()-1,m1.size()-1+m2.nhi()); }
        inline T getX() const { return x; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline const GenBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,m1,m2,m0);
        }
    protected:
        T x;
        const GenUpperTriMatrix<T1>& m1;
        const GenBandMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class ProdBU : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdBU(
            T _x, const GenBandMatrix<T1>& _m1,
            const GenUpperTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.rowsize() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline ptrdiff_t nlo() const { return m1.nlo(); }
        inline ptrdiff_t nhi() const
        { return TMV_MIN(rowsize()-1,m2.size()-1+m1.nhi()); }
        inline T getX() const { return x; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,m1,m2,m0);
        }
    protected:
        T x;
        const GenBandMatrix<T1>& m1;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T>
    inline BandMatrixView<T> operator*=(
        BandMatrixView<T> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nhi() == m2.size()-1);
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator*=(
        BandMatrixView<CT> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nhi() == m2.size()-1);
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, const ProdUB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, const ProdUB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, const ProdUB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, const ProdUB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, const ProdBU<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, const ProdBU<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, const ProdBU<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, const ProdBU<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<T> operator*=(
        UpperTriMatrixView<T> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nlo() == 0);
        TMVAssert(!m1.isunit());
        MultMM<false>(T(1),m1,m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator*=(
        UpperTriMatrixView<CT> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nlo() == 0);
        TMVAssert(!m1.isunit());
        MultMM<false>(T(1),m1,m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T, typename T1, typename T2>
    inline UpperTriMatrixView<T> operator+=(
        UpperTriMatrixView<T> m, const ProdUB<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nlo() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(
        UpperTriMatrixView<CT> m, const ProdUB<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nlo() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline UpperTriMatrixView<T> operator-=(
        UpperTriMatrixView<T> m, const ProdUB<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nlo() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            -pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(
        UpperTriMatrixView<CT> m, const ProdUB<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nlo() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            -pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline UpperTriMatrixView<T> operator+=(
        UpperTriMatrixView<T> m, const ProdBU<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nlo() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(
        UpperTriMatrixView<CT> m, const ProdBU<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nlo() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline UpperTriMatrixView<T> operator-=(
        UpperTriMatrixView<T> m, const ProdBU<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nlo() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            -pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(
        UpperTriMatrixView<CT> m, const ProdBU<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nlo() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            -pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T, typename T1, typename T2>
    class ProdLB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdLB(
            T _x, const GenLowerTriMatrix<T1>& _m1,
            const GenBandMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.colsize()); }
        inline ptrdiff_t colsize() const { return m2.colsize(); }
        inline ptrdiff_t rowsize() const { return m2.rowsize(); }
        inline ptrdiff_t nlo() const
        { return TMV_MIN(colsize()-1,m1.size()-1+m2.nlo()); }
        inline ptrdiff_t nhi() const { return m2.nhi(); }
        inline T getX() const { return x; }
        inline const GenLowerTriMatrix<T1>& getM1() const { return m1; }
        inline const GenBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,m1,m2,m0);
        }
    protected:
        T x;
        const GenLowerTriMatrix<T1>& m1;
        const GenBandMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class ProdBL : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdBL(
            T _x, const GenBandMatrix<T1>& _m1,
            const GenLowerTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.rowsize() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline ptrdiff_t nlo() const
        { return TMV_MIN(colsize()-1,m2.size()-1+m1.nlo()); }
        inline ptrdiff_t nhi() const { return m1.nhi(); }
        inline T getX() const { return x; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,m1,m2,m0);
        }
    protected:
        T x;
        const GenBandMatrix<T1>& m1;
        const GenLowerTriMatrix<T2>& m2;
    };

    template <typename T>
    inline BandMatrixView<T> operator*=(
        BandMatrixView<T> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nlo() == m2.size()-1);
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator*=(
        BandMatrixView<CT> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nlo() == m2.size()-1);
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, const ProdLB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, const ProdLB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, const ProdLB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, const ProdLB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator+=(
        BandMatrixView<T> m, const ProdBL<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator+=(
        BandMatrixView<CT> m, const ProdBL<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline BandMatrixView<T> operator-=(
        BandMatrixView<T> m, const ProdBL<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline BandMatrixView<CT> operator-=(
        BandMatrixView<CT> m, const ProdBL<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nhi() >= pmm.nhi());
        TMVAssert(m.nlo() >= pmm.nlo());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<T> operator*=(
        LowerTriMatrixView<T> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nhi() == 0);
        TMVAssert(!m1.isunit());
        MultMM<false>(T(1),m1,m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator*=(
        LowerTriMatrixView<CT> m1, const GenBandMatrix<T>& m2)
    {
        TMVAssert(m2.colsize() == m1.size());
        TMVAssert(m2.rowsize() == m1.size());
        TMVAssert(m2.nhi() == 0);
        TMVAssert(!m1.isunit());
        MultMM<false>(T(1),m1,m2,BandMatrixViewOf(m1));
        return m1;
    }

    template <typename T, typename T1, typename T2>
    inline LowerTriMatrixView<T> operator+=(
        LowerTriMatrixView<T> m, const ProdLB<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(
        LowerTriMatrixView<CT> m, const ProdLB<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline LowerTriMatrixView<T> operator-=(
        LowerTriMatrixView<T> m, const ProdLB<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            -pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(
        LowerTriMatrixView<CT> m, const ProdLB<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            -pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline LowerTriMatrixView<T> operator+=(
        LowerTriMatrixView<T> m, const ProdBL<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(
        LowerTriMatrixView<CT> m, const ProdBL<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline LowerTriMatrixView<T> operator-=(
        LowerTriMatrixView<T> m, const ProdBL<T,T1,T2>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            -pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(
        LowerTriMatrixView<CT> m, const ProdBL<T,T,T>& pmm)
    {
        TMVAssert(m.size() == pmm.colsize());
        TMVAssert(m.size() == pmm.rowsize());
        TMVAssert(pmm.nhi() == 0);
        TMVAssert(!m.isunit());
        MultMM<true>(
            -pmm.getX(),pmm.getM1(),pmm.getM2(),BandMatrixViewOf(m));
        return m;
    }


#define PRODMM ProdBU
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXU
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX2
#undef PRODXM2
#define PRODMM ProdBL
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM2 ProdXL
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM ProdUB
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXB
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef PRODXM1
#define PRODMM ProdLB
#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // BandMatrix / % TriMatrix
    //

    template <typename T, typename T1, typename T2>
    class QuotBU : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotBU(
            const T _x, const GenBandMatrix<T1>& _m1,
            const GenUpperTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.colsize() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline ptrdiff_t nlo() const { return m1.nlo(); }
        inline ptrdiff_t nhi() const
        { return TMV_MIN(rowsize()-1,m2.size()-1+m1.nhi()); }
        inline T getX() const { return x; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            UpperTriMatrix<T2> m2inv(m2.inverse());
            MultMM<false>(x,BandMatrixViewOf(m2inv),m1,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            UpperTriMatrix<T2> m2inv(m2.inverse());
            MultMM<false>(x,BandMatrixViewOf(m2inv),m1,m0);
        }
    protected:
        const T x;
        const GenBandMatrix<T1>& m1;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class RQuotBU : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline RQuotBU(
            const T _x, const GenBandMatrix<T1>& _m1,
            const GenUpperTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.rowsize() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline ptrdiff_t nlo() const { return m1.nlo(); }
        inline ptrdiff_t nhi() const
        { return TMV_MIN(rowsize()-1,m2.size()-1+m1.nhi()); }
        inline T getX() const { return x; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            UpperTriMatrix<T2> m2inv(m2.inverse());
            MultMM<false>(x,m1,BandMatrixViewOf(m2inv),m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            UpperTriMatrix<T2> m2inv(m2.inverse());
            MultMM<false>(x,m1,BandMatrixViewOf(m2inv),m0);
        }
    protected:
        const T x;
        const GenBandMatrix<T1>& m1;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T>
    inline BandMatrixView<T> operator/=(
        BandMatrixView<T> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.nhi() == m2.size()-1);
        MultMM<false>(
            T(1),BandMatrixViewOf(UpperTriMatrix<T>(m2.inverse())),m1,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator/=(
        BandMatrixView<CT> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.nhi() == m2.size()-1);
        MultMM<false>(
            T(1),BandMatrixViewOf(UpperTriMatrix<T>(m2.inverse())),m1,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<T> operator%=(
        BandMatrixView<T> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nhi() == m2.size()-1);
        MultMM<false>(
            T(1),m1,BandMatrixViewOf(UpperTriMatrix<T>(m2.inverse())),m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator%=(
        BandMatrixView<CT> m1, const GenUpperTriMatrix<T>& m2)
    {
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nhi() == m2.size()-1);
        MultMM<false>(
            T(1),m1,BandMatrixViewOf(UpperTriMatrix<T>(m2.inverse())),m1);
        return m1;
    }

    template <typename T, typename T1, typename T2>
    class QuotBL : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotBL(
            const T _x, const GenBandMatrix<T1>& _m1,
            const GenLowerTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.colsize() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline ptrdiff_t nlo() const
        { return TMV_MIN(colsize()-1,m2.size()-1+m1.nlo()); }
        inline ptrdiff_t nhi() const { return m1.nhi(); }
        inline T getX() const { return x; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            LowerTriMatrix<T2> m2inv(m2.inverse());
            MultMM<false>(x,BandMatrixViewOf(m2inv),m1,m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            LowerTriMatrix<T2> m2inv(m2.inverse());
            MultMM<false>(x,BandMatrixViewOf(m2inv),m1,m0);
        }
    protected:
        const T x;
        const GenBandMatrix<T1>& m1;
        const GenLowerTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class RQuotBL : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline RQuotBL(
            const T _x, const GenBandMatrix<T1>& _m1,
            const GenLowerTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.rowsize() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline ptrdiff_t nlo() const
        { return TMV_MIN(colsize()-1,m2.size()-1+m1.nlo()); }
        inline ptrdiff_t nhi() const { return m1.nhi(); }
        inline T getX() const { return x; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(BandMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            LowerTriMatrix<T2> m2inv(m2.inverse());
            MultMM<false>(x,m1,BandMatrixViewOf(m2inv),m0);
        }
        inline void assignToB(BandMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nlo() >= nlo());
            TMVAssert(m0.nhi() >= nhi());
            LowerTriMatrix<T2> m2inv(m2.inverse());
            MultMM<false>(x,m1,BandMatrixViewOf(m2inv),m0);
        }
    protected:
        const T x;
        const GenBandMatrix<T1>& m1;
        const GenLowerTriMatrix<T2>& m2;
    };

    template <typename T>
    inline BandMatrixView<T> operator/=(
        BandMatrixView<T> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.nlo() == m2.size()-1);
        MultMM<false>(
            T(1),BandMatrixViewOf(LowerTriMatrix<T>(m2.inverse())),m1,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator/=(
        BandMatrixView<CT> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.nlo() == m2.size()-1);
        MultMM<false>(
            T(1),BandMatrixViewOf(LowerTriMatrix<T>(m2.inverse())),m1,m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<T> operator%=(
        BandMatrixView<T> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nlo() == m2.size()-1);
        MultMM<false>(
            T(1),m1,BandMatrixViewOf(LowerTriMatrix<T>(m2.inverse())),m1);
        return m1;
    }

    template <typename T>
    inline BandMatrixView<CT> operator%=(
        BandMatrixView<CT> m1, const GenLowerTriMatrix<T>& m2)
    {
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m1.nlo() == m2.size()-1);
        MultMM<false>(
            T(1),m1,BandMatrixViewOf(LowerTriMatrix<T>(m2.inverse())),m1);
        return m1;
    }

#define QUOTMM QuotBL
#define QUOTXM QuotXL
#define RQUOTMM RQuotBL
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXL
#include "tmv/TMV_AuxQuotMM.h"
#include "tmv/TMV_AuxQuotMMa.h"
#undef QUOTMM
#undef QUOTXM
#undef RQUOTMM
#undef GENMATRIX2
#undef PRODXM2
#define QUOTMM QuotBU
#define QUOTXM QuotXU
#define RQUOTMM RQuotBU
#define GENMATRIX2 GenUpperTriMatrix
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

#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXB
#define QUOTXM QuotXB
#define TQUOTMM TransientQuotMB
#define TRQUOTMM TransientRQuotMB
#include "tmv/TMV_AuxTQuotMM.h"
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL
#include "tmv/TMV_AuxTQuotMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM

}

#undef CT

#endif
