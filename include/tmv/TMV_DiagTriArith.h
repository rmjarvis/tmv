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


#ifndef TMV_DiagTriArith_H
#define TMV_DiagTriArith_H

#include "tmv/TMV_DiagTriArithFunc.h"

#define CT std::complex<T>

namespace tmv {

    //
    // DiagMatrix + TriMatrix
    //

    template <typename T, typename T1, typename T2>
    class SumDU : public UpperTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumDU(
            T _x1, const GenDiagMatrix<T1>& _m1,
            T _x2, const GenUpperTriMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX1() const { return x1; }
        inline const GenDiagMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToU(UpperTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(m0.dt() == dt());
            if (SameStorage(m0,m1)) {
                DiagMatrix<T1> temp = m1;
                MultXM(x2,m0=m2);
                AddVV(x1,temp.diag(),m0.diag());
            } else {
                MultXM(x2,m0=m2);
                AddVV(x1,m1.diag(),m0.diag());
            }
        }
        inline void assignToU(UpperTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(m0.dt() == dt());
            if (SameStorage(m0,m1)) {
                DiagMatrix<T1> temp = m1;
                MultXM(x2,m0=m2);
                AddVV(x1,temp.diag(),m0.diag());
            } else {
                MultXM(x2,m0=m2);
                AddVV(x1,m1.diag(),m0.diag());
            }
        }
    private:
        const T x1;
        const GenDiagMatrix<T1>& m1;
        const T x2;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T>
    inline UpperTriMatrixView<T> operator+=(
        UpperTriMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(
        UpperTriMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline UpperTriMatrixView<T> operator-=(
        UpperTriMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(
        UpperTriMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T, typename T2>
    inline UpperTriMatrixView<T> operator+=(
        UpperTriMatrixView<T> m, const ProdXD<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(
        UpperTriMatrixView<CT> m, const ProdXD<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T, typename T2>
    inline UpperTriMatrixView<T> operator-=(
        UpperTriMatrixView<T> m, const ProdXD<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(
        UpperTriMatrixView<CT> m, const ProdXD<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T, typename T1, typename T2>
    class SumDL : public LowerTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumDL(
            T _x1, const GenDiagMatrix<T1>& _m1,
            T _x2, const GenLowerTriMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX1() const { return x1; }
        inline const GenDiagMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToL(LowerTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(m0.dt() == dt());
            if (SameStorage(m0,m1)) {
                DiagMatrix<T1> temp = m1;
                MultXM(x2,m0=m2);
                AddVV(x1,temp.diag(),m0.diag());
            } else {
                MultXM(x2,m0=m2);
                AddVV(x1,m1.diag(),m0.diag());
            }
        }
        inline void assignToL(LowerTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(m0.dt() == dt());
            if (SameStorage(m0,m1)) {
                DiagMatrix<T1> temp = m1;
                MultXM(x2,m0=m2);
                AddVV(x1,temp.diag(),m0.diag());
            } else {
                MultXM(x2,m0=m2);
                AddVV(x1,m1.diag(),m0.diag());
            }
        }
    private:
        const T x1;
        const GenDiagMatrix<T1>& m1;
        const T x2;
        const GenLowerTriMatrix<T2>& m2;
    };

    template <typename T>
    inline LowerTriMatrixView<T> operator+=(
        LowerTriMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(
        LowerTriMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddVV(T(1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline LowerTriMatrixView<T> operator-=(
        LowerTriMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(
        LowerTriMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m1.isunit());
        AddVV(T(-1),m2.diag(),m1.diag());
        return m1;
    }

    template <typename T, typename T2>
    inline LowerTriMatrixView<T> operator+=(
        LowerTriMatrixView<T> m, const ProdXD<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(
        LowerTriMatrixView<CT> m, const ProdXD<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddVV(pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T, typename T2>
    inline LowerTriMatrixView<T> operator-=(
        LowerTriMatrixView<T> m, const ProdXD<T,T2>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(
        LowerTriMatrixView<CT> m, const ProdXD<T,T>& pxm)
    {
        TMVAssert(m.size() == pxm.size());
        TMVAssert(!m.isunit());
        AddVV(-pxm.getX(),pxm.getM().diag(),m.diag());
        return m;
    }

#define SUMMM SumDU
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXU
#include "tmv/TMV_AuxSumMM.h"
#include "tmv/TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX2
#undef PRODXM2
#define SUMMM SumDL
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM2 ProdXL
#include "tmv/TMV_AuxSumMM.h"
#include "tmv/TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // DiagMatrix * TriMatrix
    //

    template <typename T, typename T1, typename T2>
    class ProdDU : public UpperTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdDU(
            T _x, const GenDiagMatrix<T1>& _m1,
            const GenUpperTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenDiagMatrix<T1>& getM1() const { return m1; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToU(UpperTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(m0.dt() == dt());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToU(UpperTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(m0.dt() == dt());
            MultMM<false>(x,m1,m2,m0);
        }
    protected:
        T x;
        const GenDiagMatrix<T1>& m1;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class ProdUD : public UpperTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdUD(
            T _x, const GenUpperTriMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToU(UpperTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(m0.dt() == dt());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToU(UpperTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(m0.dt() == dt());
            MultMM<false>(x,m1,m2,m0);
        }
    protected:
        T x;
        const GenUpperTriMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

    template <typename T>
    inline UpperTriMatrixView<T> operator*=(
        UpperTriMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(!m1.isunit());
        TMVAssert(m1.size() == m2.size());
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator*=(
        UpperTriMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(!m1.isunit());
        TMVAssert(m1.size() == m2.size());
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T, typename T1, typename T2>
    inline UpperTriMatrixView<T> operator+=(
        UpperTriMatrixView<T> m, const ProdDU<T,T1,T2>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(
        UpperTriMatrixView<CT> m, const ProdDU<T,T,T>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline UpperTriMatrixView<T> operator-=(
        UpperTriMatrixView<T> m, const ProdDU<T,T1,T2>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(
        UpperTriMatrixView<CT> m, const ProdDU<T,T,T>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline UpperTriMatrixView<T> operator+=(
        UpperTriMatrixView<T> m, const ProdUD<T,T1,T2>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator+=(
        UpperTriMatrixView<CT> m, const ProdUD<T,T,T>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline UpperTriMatrixView<T> operator-=(
        UpperTriMatrixView<T> m, const ProdUD<T,T1,T2>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline UpperTriMatrixView<CT> operator-=(
        UpperTriMatrixView<CT> m, const ProdUD<T,T,T>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    class ProdDL : public LowerTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdDL(
            T _x, const GenDiagMatrix<T1>& _m1,
            const GenLowerTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenDiagMatrix<T1>& getM1() const { return m1; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToL(LowerTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(m0.dt() == dt());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToL(LowerTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(m0.dt() == dt());
            MultMM<false>(x,m1,m2,m0);
        }
    protected:
        T x;
        const GenDiagMatrix<T1>& m1;
        const GenLowerTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class ProdLD : public LowerTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdLD(
            T _x, const GenLowerTriMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenLowerTriMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToL(LowerTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(m0.dt() == dt());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToL(LowerTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(m0.dt() == dt());
            MultMM<false>(x,m1,m2,m0);
        }
    protected:
        T x;
        const GenLowerTriMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

    template <typename T>
    inline LowerTriMatrixView<T> operator*=(
        LowerTriMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(!m1.isunit());
        TMVAssert(m1.size() == m2.size());
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator*=(
        LowerTriMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    {
        TMVAssert(!m1.isunit());
        TMVAssert(m1.size() == m2.size());
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <typename T, typename T1, typename T2>
    inline LowerTriMatrixView<T> operator+=(
        LowerTriMatrixView<T> m, const ProdDL<T,T1,T2>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(
        LowerTriMatrixView<CT> m, const ProdDL<T,T,T>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline LowerTriMatrixView<T> operator-=(
        LowerTriMatrixView<T> m, const ProdDL<T,T1,T2>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(
        LowerTriMatrixView<CT> m, const ProdDL<T,T,T>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline LowerTriMatrixView<T> operator+=(
        LowerTriMatrixView<T> m, const ProdLD<T,T1,T2>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator+=(
        LowerTriMatrixView<CT> m, const ProdLD<T,T,T>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T, typename T1, typename T2>
    inline LowerTriMatrixView<T> operator-=(
        LowerTriMatrixView<T> m, const ProdLD<T,T1,T2>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <typename T>
    inline LowerTriMatrixView<CT> operator-=(
        LowerTriMatrixView<CT> m, const ProdLD<T,T,T>& pmm)
    {
        TMVAssert(!m.isunit());
        TMVAssert(m.size() == pmm.size());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

#define PRODMM ProdUD
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXD
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef PRODXM1
#define PRODMM ProdLD
#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM ProdDU
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXU
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX2
#undef PRODXM2
#define PRODMM ProdDL
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM2 ProdXL
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // TriMatrix / % DiagMatrix
    //

    template <typename T, typename T1, typename T2>
    class QuotUD : public UpperTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotUD(
            const T _x, const GenUpperTriMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToU(UpperTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit());
            MultMM<false>(x,DiagMatrix<T2>(m2.inverse()),m1,m0);
        }
        inline void assignToU(UpperTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit());
            MultMM<false>(x,DiagMatrix<T2>(m2.inverse()),m1,m0);
        }
    protected:
        const T x;
        const GenUpperTriMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class RQuotUD : public UpperTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline RQuotUD(
            const T _x, const GenUpperTriMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToU(UpperTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit());
            MultMM<false>(x,m1,DiagMatrix<T2>(m2.inverse()),m0);
        }
        inline void assignToU(UpperTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit());
            MultMM<false>(x,m1,DiagMatrix<T2>(m2.inverse()),m0);
        }
    protected:
        const T x;
        const GenUpperTriMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

    template <typename T>
    inline UpperTriMatrixView<T> operator/=(
        UpperTriMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    { MultMM<false>(T(1),DiagMatrix<T>(m2.inverse()),m1,m1); return m1; }

    template <typename T>
    inline UpperTriMatrixView<CT> operator/=(
        UpperTriMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    { MultMM<false>(T(1),DiagMatrix<T>(m2.inverse()),m1,m1); return m1; }

    template <typename T>
    inline UpperTriMatrixView<T> operator%=(
        UpperTriMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    { MultMM<false>(T(1),m1,DiagMatrix<T>(m2.inverse()),m1); return m1; }

    template <typename T>
    inline UpperTriMatrixView<CT> operator%=(
        UpperTriMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    { MultMM<false>(T(1),m1,DiagMatrix<T>(m2.inverse()),m1); return m1; }


    template <typename T, typename T1, typename T2>
    class QuotLD : public LowerTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotLD(
            const T _x, const GenLowerTriMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenLowerTriMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToL(LowerTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit());
            MultMM<false>(x,DiagMatrix<T2>(m2.inverse()),m1,m0);
        }
        inline void assignToL(LowerTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit());
            MultMM<false>(x,DiagMatrix<T2>(m2.inverse()),m1,m0);
        }
    protected:
        const T x;
        const GenLowerTriMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class RQuotLD : public LowerTriMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline RQuotLD(
            const T _x, const GenLowerTriMatrix<T1>& _m1,
            const GenDiagMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t size() const { return m1.size(); }
        inline DiagType dt() const { return NonUnitDiag; }
        inline T getX() const { return x; }
        inline const GenLowerTriMatrix<T1>& getM1() const { return m1; }
        inline const GenDiagMatrix<T2>& getM2() const { return m2; }
        inline void assignToL(LowerTriMatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit());
            MultMM<false>(x,m1,DiagMatrix<T2>(m2.inverse()),m0);
        }
        inline void assignToL(LowerTriMatrixView<complex_type> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(!m0.isunit());
            MultMM<false>(x,m1,DiagMatrix<T2>(m2.inverse()),m0);
        }
    protected:
        const T x;
        const GenLowerTriMatrix<T1>& m1;
        const GenDiagMatrix<T2>& m2;
    };

    template <typename T>
    inline LowerTriMatrixView<T> operator/=(
        LowerTriMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    { MultMM<false>(T(1),DiagMatrix<T>(m2.inverse()),m1,m1); return m1; }

    template <typename T>
    inline LowerTriMatrixView<CT> operator/=(
        LowerTriMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    { MultMM<false>(T(1),DiagMatrix<T>(m2.inverse()),m1,m1); return m1; }

    template <typename T>
    inline LowerTriMatrixView<T> operator%=(
        LowerTriMatrixView<T> m1, const GenDiagMatrix<T>& m2)
    { MultMM<false>(T(1),m1,DiagMatrix<T>(m2.inverse()),m1); return m1; }

    template <typename T>
    inline LowerTriMatrixView<CT> operator%=(
        LowerTriMatrixView<CT> m1, const GenDiagMatrix<T>& m2)
    { MultMM<false>(T(1),m1,DiagMatrix<T>(m2.inverse()),m1); return m1; }


#define QUOTMM QuotUD
#define QUOTXM QuotXD
#define RQUOTMM RQuotUD
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXD
#include "tmv/TMV_AuxQuotMM.h"
#include "tmv/TMV_AuxQuotMMa.h"
#undef QUOTMM
#undef RQUOTMM
#undef GENMATRIX1
#undef PRODXM1
#define QUOTMM QuotLD
#define RQUOTMM RQuotLD
#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL
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
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXU
#define QUOTXM QuotXU
#define TQUOTMM TransientQuotMU
#define TRQUOTMM TransientRQuotMU
#include "tmv/TMV_AuxTQuotMM.h"
#undef GENMATRIX2
#undef PRODXM2
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM2 ProdXL
#define QUOTXM QuotXL
#define TQUOTMM TransientQuotML
#define TRQUOTMM TransientRQuotML
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
