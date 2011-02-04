///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef TMV_BandSymBandArith_H
#define TMV_BandSymBandArith_H

#define CT std::complex<T>

namespace tmv {

    //
    // BandMatrix + SymBandMatrix
    //

    template <class T, class T1, class T2> 
    class SumBsB : public BandMatrixComposite<T> 
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumBsB(
            T _x1, const GenBandMatrix<T1>& _m1, 
            T _x2, const GenSymBandMatrix<T2>& _m2
        ) :
            x1(_x1),m1(_m1),x2(_x2),m2(_m2)
        { 
            TMVAssert(m1.colsize() == m2.size()); 
            TMVAssert(m1.rowsize() == m2.size()); 
        }
        inline size_t colsize() const { return m2.size(); }
        inline size_t rowsize() const { return m2.size(); }
        inline int nlo() const { return TMV_MAX(m1.nlo(),m2.nlo()); }
        inline int nhi() const { return TMV_MAX(m1.nhi(),m2.nhi()); }
        inline StorageType stor() const 
        { return m1.stor() == m2.stor() ? BaseStorOf(m2) : DiagMajor; }
        inline T getX1() const { return x1; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline T getX2() const { return x2; }
        inline const GenSymBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(const BandMatrixView<real_type>& m0) const
        { 
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            if (SameStorage(m0,m1)) {
                BandMatrix<T1> m1x = m1;
                MultXM(x2,m0=m2);
                AddMM(x1,m1x,m0);
            } else {
                MultXM(x2,m0=m2);
                AddMM(x1,m1,m0);
            }
        }
        inline void assignToB(const BandMatrixView<complex_type>& m0) const
        { 
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
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
        const GenSymBandMatrix<T2>& m2;
    };

    template <class T> 
    inline const BandMatrixView<T>& operator+=(
        const BandMatrixView<T>& m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() >= m2.nlo());
        TMVAssert(m1.nhi() >= m2.nhi());
        AddMM(T(1),m2.upperBand(),m1.diagRange(0,m1.nhi()+1));
        if (m2.nlo() > 0) 
            AddMM(T(1),m2.lowerBandOff(),m1.diagRange(-m1.nlo(),0));
        return m1;
    }

    template <class T> 
    inline const BandMatrixView<CT>& operator+=(
        const BandMatrixView<CT>& m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() >= m2.nlo());
        TMVAssert(m1.nhi() >= m2.nhi());
        AddMM(T(1),m2.upperBand(),m1.diagRange(0,m1.nhi()+1));
        if (m2.nlo() > 0) 
            AddMM(T(1),m2.lowerBandOff(),m1.diagRange(-m1.nlo(),0));
        return m1;
    }

    template <class T> 
    inline const BandMatrixView<T>& operator-=(
        const BandMatrixView<T>& m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() >= m2.nlo());
        TMVAssert(m1.nhi() >= m2.nhi());
        AddMM(T(-1),m2.upperBand(),m1.diagRange(0,m1.nhi()+1));
        if (m2.nlo() > 0) 
            AddMM(T(-1),m2.lowerBandOff(),m1.diagRange(-m1.nlo(),0));
        return m1;
    }

    template <class T> 
    inline const BandMatrixView<CT>& operator-=(
        const BandMatrixView<CT>& m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() >= m2.nlo());
        TMVAssert(m1.nhi() >= m2.nhi());
        AddMM(T(-1),m2.upperBand(),m1.diagRange(0,m1.nhi()+1));
        if (m2.nlo() > 0) 
            AddMM(T(-1),m2.lowerBandOff(),m1.diagRange(-m1.nlo(),0));
        return m1;
    }

    template <class T, class T2> 
    inline const BandMatrixView<T>& operator+=(
        const BandMatrixView<T>& m, const ProdXsB<T,T2>& pxm)
    {
        TMVAssert(m.colsize() == pxm.colsize());
        TMVAssert(m.rowsize() == pxm.rowsize());
        TMVAssert(m.nlo() >= pxm.nlo());
        TMVAssert(m.nhi() >= pxm.nhi());
        AddMM(pxm.getX(),pxm.getM().upperBand(),m.diagRange(0,m.nhi()+1));
        if (pxm.nlo() > 0)
            AddMM(pxm.getX(),pxm.getM().lowerBandOff(),
                  m.diagRange(-m.nlo(),0));
        return m;
    }

    template <class T> 
    inline const BandMatrixView<CT>& operator+=(
        const BandMatrixView<CT>& m, const ProdXsB<T,T>& pxm)
    {
        TMVAssert(m.colsize() == pxm.colsize());
        TMVAssert(m.rowsize() == pxm.rowsize());
        TMVAssert(m.nlo() >= pxm.nlo());
        TMVAssert(m.nhi() >= pxm.nhi());
        AddMM(pxm.getX(),pxm.getM().upperBand(),m.diagRange(0,m.nhi()+1));
        if (pxm.nlo() > 0)
            AddMM(pxm.getX(),pxm.getM().lowerBandOff(),
                  m.diagRange(-m.nlo(),0));
        return m;
    }

    template <class T, class T2> 
    inline const BandMatrixView<T>& operator-=(
        const BandMatrixView<T>& m, const ProdXsB<T,T2>& pxm)
    {
        TMVAssert(m.colsize() == pxm.colsize());
        TMVAssert(m.rowsize() == pxm.rowsize());
        TMVAssert(m.nlo() >= pxm.nlo());
        TMVAssert(m.nhi() >= pxm.nhi());
        AddMM(-pxm.getX(),pxm.getM().upperBand(),m.diagRange(0,m.nhi()+1));
        if (pxm.nlo() > 0)
            AddMM(-pxm.getX(),pxm.getM().lowerBandOff(),
                  m.diagRange(-m.nlo(),0));
        return m;
    }

    template <class T> 
    inline const BandMatrixView<CT>& operator-=(
        const BandMatrixView<CT>& m, const ProdXsB<T,T>& pxm)
    {
        TMVAssert(m.colsize() == pxm.colsize());
        TMVAssert(m.rowsize() == pxm.rowsize());
        TMVAssert(m.nlo() >= pxm.nlo());
        TMVAssert(m.nhi() >= pxm.nhi());
        AddMM(-pxm.getX(),pxm.getM().upperBand(),m.diagRange(0,m.nhi()+1));
        if (pxm.nlo() > 0)
            AddMM(-pxm.getX(),pxm.getM().lowerBandOff(),
                  m.diagRange(-m.nlo(),0));
        return m;
    }


#define SUMMM SumBsB
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXsB
#include "tmv/TMV_AuxSumMM.h"
#include "tmv/TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // BandMatrix * SymBandMatrix
    //

    template <class T, class T1, class T2> 
    class ProdBsB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdBsB(
            T _x, const GenBandMatrix<T1>& _m1,
            const GenSymBandMatrix<T2>& _m2
        ) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.rowsize() == m2.colsize()); }
        inline size_t colsize() const { return m1.colsize(); }
        inline size_t rowsize() const { return m2.rowsize(); }
        inline int nlo() const 
        { return TMV_MIN(int(colsize()-1),m1.nlo()+m2.nlo()); }
        inline int nhi() const 
        { return TMV_MIN(int(rowsize()-1),m1.nhi()+m2.nhi()); }
        inline StorageType stor() const
        {
            return m1.isrm() ? RowMajor : m2.iscm() ? ColMajor :
                (m1.iscm() && m2.isrm()) ? ColMajor : DiagMajor;
        }
        inline T getX() const { return x; }
        inline const GenBandMatrix<T1>& getM1() const { return m1; }
        inline const GenSymBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(const BandMatrixView<real_type>& m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToB(const BandMatrixView<complex_type>& m0) const
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
        const GenSymBandMatrix<T2>& m2;
    };

    template <class T, class T1, class T2> 
    class ProdsBB : public BandMatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdsBB(
            T _x, const GenSymBandMatrix<T1>& _m1,
            const GenBandMatrix<T2>& _m2
        ) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.rowsize() == m2.colsize()); }
        inline size_t colsize() const { return m1.colsize(); }
        inline size_t rowsize() const { return m2.rowsize(); }
        inline int nlo() const 
        { return TMV_MIN(int(colsize()-1),m1.nlo()+m2.nlo()); }
        inline int nhi() const
        { return TMV_MIN(int(rowsize()-1),m1.nhi()+m2.nhi()); }
        inline StorageType stor() const 
        {
            return m1.isrm() ? RowMajor : m2.iscm() ? ColMajor :
                (m1.iscm() && m2.isrm()) ? ColMajor : DiagMajor;
        }
        inline T getX() const { return x; }
        inline const GenSymBandMatrix<T1>& getM1() const { return m1; }
        inline const GenBandMatrix<T2>& getM2() const { return m2; }
        inline void assignToB(const BandMatrixView<real_type>& m0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,m1,m2,m0);
        }
        inline void assignToB(const BandMatrixView<complex_type>& m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(m0.nhi() >= nhi());
            TMVAssert(m0.nlo() >= nlo());
            MultMM<false>(x,m1,m2,m0);
        }
    protected:
        T x;
        const GenSymBandMatrix<T1>& m1;
        const GenBandMatrix<T2>& m2;
    };

    template <class T> 
    inline const BandMatrixView<T>& operator*=(
        const BandMatrixView<T>& m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m2.nlo() == 0 || m1.nlo() == int(m1.colsize())-1);
        TMVAssert(m2.nhi() == 0 || m1.nhi() == int(m1.rowsize())-1);
        MultMM<false>(T(1),m1,m2,m1); 
        return m1;
    }

    template <class T> 
    inline const BandMatrixView<CT>& operator*=(
        const BandMatrixView<CT>& m1, const GenSymBandMatrix<T>& m2)
    {
        TMVAssert(m2.nlo() == 0 || m1.nlo() == int(m1.colsize())-1);
        TMVAssert(m2.nhi() == 0 || m1.nhi() == int(m1.rowsize())-1);
        MultMM<false>(T(1),m1,m2,m1);
        return m1;
    }

    template <class T, class T1, class T2>
    inline const BandMatrixView<T>& operator+=(
        const BandMatrixView<T>& m, const ProdBsB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <class T> 
    inline const BandMatrixView<CT>& operator+=(
        const BandMatrixView<CT>& m, const ProdBsB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <class T, class T1, class T2>
    inline const BandMatrixView<T>& operator-=(
        const BandMatrixView<T>& m, const ProdBsB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <class T> 
    inline const BandMatrixView<CT>& operator-=(
        const BandMatrixView<CT>& m, const ProdBsB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <class T, class T1, class T2>
    inline const BandMatrixView<T>& operator+=(
        const BandMatrixView<T>& m, const ProdsBB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <class T> 
    inline const BandMatrixView<CT>& operator+=(
        const BandMatrixView<CT>& m, const ProdsBB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <class T, class T1, class T2>
    inline const BandMatrixView<T>& operator-=(
        const BandMatrixView<T>& m, const ProdsBB<T,T1,T2>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <class T> 
    inline const BandMatrixView<CT>& operator-=(
        const BandMatrixView<CT>& m, const ProdsBB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        TMVAssert(m.nlo() >= pmm.nlo());
        TMVAssert(m.nhi() >= pmm.nhi());
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m);
        return m;
    }

    template <class T, class T1, class T2, class T3>
    inline const MatrixView<T>& operator+=(
        const MatrixView<T>& m, const ProdBsB<T1,T2,T3>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixView<T>(m,pmm.nlo(),pmm.nhi()) += pmm;
        return m;
    }

    template <class T> 
    inline const MatrixView<CT>& operator+=(
        const MatrixView<CT>& m, const ProdBsB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixView<CT>(m,pmm.nlo(),pmm.nhi()) += pmm;
        return m;
    }

    template <class T, class T1, class T2, class T3>
    inline const MatrixView<T>& operator-=(
        const MatrixView<T>& m, const ProdBsB<T1,T2,T3>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixView<T>(m,pmm.nlo(),pmm.nhi()) -= pmm;
        return m;
    }

    template <class T> 
    inline const MatrixView<CT>& operator-=(
        const MatrixView<CT>& m, const ProdBsB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixView<CT>(m,pmm.nlo(),pmm.nhi()) -= pmm;
        return m;
    }
    template <class T, class T1, class T2, class T3>
    inline const MatrixView<T>& operator+=(
        const MatrixView<T>& m, const ProdsBB<T1,T2,T3>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixView<T>(m,pmm.nlo(),pmm.nhi()) += pmm;
        return m;
    }

    template <class T> 
    inline const MatrixView<CT>& operator+=(
        const MatrixView<CT>& m, const ProdsBB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixView<CT>(m,pmm.nlo(),pmm.nhi()) += pmm;
        return m;
    }

    template <class T, class T1, class T2, class T3>
    inline const MatrixView<T>& operator-=(
        const MatrixView<T>& m, const ProdsBB<T1,T2,T3>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixView<T>(m,pmm.nlo(),pmm.nhi()) -= pmm;
        return m;
    }

    template <class T> 
    inline const MatrixView<CT>& operator-=(
        const MatrixView<CT>& m, const ProdsBB<T,T,T>& pmm)
    {
        TMVAssert(m.colsize() == pmm.colsize());
        TMVAssert(m.rowsize() == pmm.rowsize());
        BandMatrixView<CT>(m,pmm.nlo(),pmm.nhi()) -= pmm;
        return m;
    }


#define PRODMM ProdsBB
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXB
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM ProdBsB
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXsB
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // SymBandMatrix / % BandMatrix
    //

#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXsB
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

#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXB
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

}

#undef CT

#endif
