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


#ifndef TMV_BandSymArith_H
#define TMV_BandSymArith_H

#define CT std::complex<T>

namespace tmv {

    //
    // BandMatrix + SymMatrix
    //

    template <typename T, typename T1, typename T2>
    class SumBS : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumBS(
            T _x1, const GenBandMatrix<T1>& _m1,
            T _x2, const GenSymMatrix<T2>& _m2) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        {
            TMVAssert(m1.colsize() == m2.size());
            TMVAssert(m1.rowsize() == m2.size());
        }
        inline ptrdiff_t colsize() const { return m2.size(); }
        inline ptrdiff_t rowsize() const { return m2.size(); }
        inline T getX1() const { return x1; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenSymMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m1)) {
                BandMatrix<T1> m1x = m1;
                MultXM(x2,m0=m2);
                AddMM(x1,m1x,m0);
            } else {
                MultXM(x2,m0=m2);
                AddMM(x1,m1,m0);
            }
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            if (SameStorage(m0,m1)) {
                BandMatrix<T1> m1x = m1;
                MultXM(x2,m0=m2);
                AddMM(x1,m1x,m0);
            } else {
                MultXM(x2,m0=m2);
                AddMM(x1,m1,m0);
            }
        }
    private:
        const T x1;
        const GenBandMatrix<T1>& m1;
        const T x2;
        const GenSymMatrix<T2>& m2;
    };

#define SUMMM SumBS
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXS
#include "tmv/TMV_AuxSumMM.h"
#include "tmv/TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // BandMatrix * SymMatrix
    //

    template <typename T, typename T1, typename T2>
    class ProdBS : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdBS(
            T _x, const GenBandMatrix<T1>& _m1,
            const GenSymMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.rowsize() == m2.colsize()); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m2.rowsize(); }
        inline T getX() const { return x; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline const GenSymMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            Matrix<T> xm2 = m2;
            MultMM<false>(x,m1,xm2,m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            Matrix<T> xm2 = m2;
            MultMM<false>(x,m1,xm2,m0);
        }
    protected:
        T x;
        const GenBandMatrix<T1>& m1;
        const GenSymMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2>
    class ProdSB : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdSB(
            T _x, const GenSymMatrix<T1>& _m1,
            const GenBandMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.rowsize() == m2.colsize()); }
        inline ptrdiff_t colsize() const { return m1.colsize(); }
        inline ptrdiff_t rowsize() const { return m2.rowsize(); }
        inline T getX() const { return x; }
        inline const GenSymMatrix<T1>& getM1() const { return m1; }
        inline const GenBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            Matrix<T> xm1 = m1;
            MultMM<false>(x,xm1,m2,m0);
        }
        inline void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            Matrix<T> xm1 = m1;
            MultMM<false>(x,xm1,m2,m0);
        }
    protected:
        T x;
        const GenSymMatrix<T1>& m1;
        const GenBandMatrix<T2>& m2;
    };

#define PRODMM ProdSB
#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXS
#define PRODXM2 ProdXB
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM ProdBS
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXS
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // SymMatrix / % BandMatrix
    //

#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXB
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
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXS
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
