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


#ifndef TMV_TriSymBandArith_H
#define TMV_TriSymBandArith_H

#include "tmv/TMV_TriBandArithFunc.h"

#define CT std::complex<T>

namespace tmv {

    //
    // TriMatrix + SymBandMatrix
    //

    template <typename T, typename T1, typename T2>
    class SumUsB : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumUsB(
            T _x1, const GenUpperTriMatrix<T1>& _m1,
            T _x2, const GenSymBandMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline T getX1() const { return x1; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenSymBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m2)) {
                BandMatrix<T2> m2x = m2;
                MultXM(x1,m0=m1);
                AddMM(x2,m2x,m0);
            } else {
                MultXM(x1,m0=m1);
                AddMM(x2,m2,m0);
            }
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m2)) {
                BandMatrix<T2> m2x = m2;
                MultXM(x1,m0=m1);
                AddMM(x2,m2x,m0);
            } else {
                MultXM(x1,m0=m1);
                AddMM(x2,m2,m0);
            }
        }
    private:
        const T x1;
        const GenUpperTriMatrix<T1>& m1;
        const T x2;
        const GenSymBandMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class SumLsB : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumLsB(
            T _x1, const GenLowerTriMatrix<T1>& _m1,
            T _x2, const GenSymBandMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline T getX1() const { return x1; }
        inline const GenLowerTriMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenSymBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m2)) {
                BandMatrix<T2> m2x = m2;
                MultXM(x1,m0=m1);
                AddMM(x2,m2x,m0);
            } else {
                MultXM(x1,m0=m1);
                AddMM(x2,m2,m0);
            }
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m2)) {
                BandMatrix<T2> m2x = m2;
                MultXM(x1,m0=m1);
                AddMM(x2,m2x,m0);
            } else {
                MultXM(x1,m0=m1);
                AddMM(x2,m2,m0);
            }
        }
    private:
        const T x1;
        const GenLowerTriMatrix<T1>& m1;
        const T x2;
        const GenSymBandMatrix<T2>& m2;
    };


#define SUMMM SumUsB
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXsB
#include "tmv/TMV_AuxSumMM.h"
#include "tmv/TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef PRODXM1
#define SUMMM SumLsB
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
    // TriMatrix * SymBandMatrix
    //

    template <typename T, typename T1, typename T2>
    class ProdUsB : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdUsB(
            T _x, const GenUpperTriMatrix<T1>& _m1,
            const GenSymBandMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m2.size(); }
        inline ptrdiff_t rowsize() const { return m2.size(); }
        inline T getX() const { return x; }
        inline const GenUpperTriMatrix<T1>& getM1() const { return m1; }
        inline const GenSymBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m1)) {
                BandMatrix<T1> m1x = m1;
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
            if (SameStorage(m0,m1)) {
                BandMatrix<T1> m1x = m1;
                MultXM(x,m0=m2);
                MultMM<false>(T(1),m1x,m0,m0);
            } else {
                MultXM(x,m0=m2);
                MultMM<false>(T(1),m1,m0,m0);
            }
        }
    protected:
        T x;
        const GenUpperTriMatrix<T1>& m1;
        const GenSymBandMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class ProdsBU : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdsBU(
            T _x, const GenSymBandMatrix<T1>& _m1,
            const GenUpperTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m2)) {
                BandMatrix<T2> m2x = m2;
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
            if (SameStorage(m0,m2)) {
                BandMatrix<T2> m2x = m2;
                MultXM(x,m0=m1);
                MultMM<false>(T(1),m0,m2x,m0);
            } else {
                MultXM(x,m0=m1);
                MultMM<false>(T(1),m0,m2,m0);
            }
        }
    protected:
        T x;
        const GenSymBandMatrix<T1>& m1;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class ProdLsB : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdLsB(
            T _x, const GenLowerTriMatrix<T1>& _m1,
            const GenSymBandMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m2.size(); }
        inline ptrdiff_t rowsize() const { return m2.size(); }
        inline T getX() const { return x; }
        inline const GenLowerTriMatrix<T1>& getM1() const { return m1; }
        inline const GenSymBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m1)) {
                BandMatrix<T1> m1x = m1;
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
            if (SameStorage(m0,m1)) {
                BandMatrix<T1> m1x = m1;
                MultXM(x,m0=m2);
                MultMM<false>(T(1),m1x,m0,m0);
            } else {
                MultXM(x,m0=m2);
                MultMM<false>(T(1),m1,m0,m0);
            }
        }
    protected:
        T x;
        const GenLowerTriMatrix<T1>& m1;
        const GenSymBandMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class ProdsBL : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdsBL(
            T _x, const GenSymBandMatrix<T1>& _m1,
            const GenLowerTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.size() == m2.size()); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m2)) {
                BandMatrix<T2> m2x = m2;
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
            if (SameStorage(m0,m2)) {
                BandMatrix<T2> m2x = m2;
                MultXM(x,m0=m1);
                MultMM<false>(T(1),m0,m2x,m0);
            } else {
                MultXM(x,m0=m1);
                MultMM<false>(T(1),m0,m2,m0);
            }
        }
    protected:
        T x;
        const GenSymBandMatrix<T1>& m1;
        const GenLowerTriMatrix<T2>& m2;
    };


#define PRODMM ProdsBU
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXU
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX2
#undef PRODXM2
#define PRODMM ProdsBL
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM2 ProdXL
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM ProdUsB
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXsB
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef PRODXM1
#define PRODMM ProdLsB
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
    // SymBandMatrix / % TriMatrix
    //

    template <typename T, typename T1, typename T2>
    class QuotsBU : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotsBU(
            const T _x, const GenSymBandMatrix<T1>& _m1,
            const GenUpperTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
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
        const GenSymBandMatrix<T1>& m1;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class RQuotsBU : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline RQuotsBU(
            const T _x, const GenSymBandMatrix<T1>& _m1,
            const GenUpperTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.rowsize() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline const GenUpperTriMatrix<T2>& getM2() const { return m2; }
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
        const GenSymBandMatrix<T1>& m1;
        const GenUpperTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class QuotsBL : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline QuotsBL(
            const T _x, const GenSymBandMatrix<T1>& _m1,
            const GenLowerTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.size() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.size(); }
        inline ptrdiff_t rowsize() const { return m1.size(); }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
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
        const GenSymBandMatrix<T1>& m1;
        const GenLowerTriMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class RQuotsBL : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline RQuotsBL(
            const T _x, const GenSymBandMatrix<T1>& _m1,
            const GenLowerTriMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert( m1.rowsize() == m2.size() ); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m1.rowsize(); }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline const GenLowerTriMatrix<T2>& getM2() const { return m2; }
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
        const GenSymBandMatrix<T1>& m1;
        const GenLowerTriMatrix<T2>& m2;
    };

#define QUOTMM QuotsBL
#define QUOTXM QuotXL
#define RQUOTMM RQuotsBL
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXL
#include "tmv/TMV_AuxQuotMM.h"
#include "tmv/TMV_AuxQuotMMa.h"
#undef QUOTMM
#undef QUOTXM
#undef RQUOTMM
#undef GENMATRIX2
#undef PRODXM2
#define QUOTMM QuotsBU
#define QUOTXM QuotXU
#define RQUOTMM RQuotsBU
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
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXsB
#define QUOTXM QuotXsB
#define TQUOTMM TransientQuotMsB
#define TRQUOTMM TransientRQuotMsB
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
