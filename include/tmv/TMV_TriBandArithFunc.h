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


#ifndef TMV_TriBandArithFunc_H
#define TMV_TriBandArithFunc_H

#include "tmv/TMV_BaseTriMatrix.h"
#include "tmv/TMV_BaseBandMatrix.h"
#include "tmv/TMV_BandMatrixArithFunc.h"

#define CT std::complex<T>

namespace tmv {

    template <typename T, typename T2>
    inline void AddMM(
        const T x2, const GenUpperTriMatrix<T2>& m2,
        BandMatrixView<T> m1)
    {
        if (m2.isunit()) {
            if (m2.size() > 1)
                AddMM(x2,BandMatrixViewOf(m2.offDiag()),
                      m1.diagRange(1,m2.size()));
            m1.diag().addToAll(x2);
        } else {
            AddMM(x2,BandMatrixViewOf(m2),m1);
        }
    }

    template <typename T, typename T2>
    inline void AddMM(
        const T x2, const GenLowerTriMatrix<T2>& m2,
        BandMatrixView<T> m1)
    {
        if (m2.isunit()) {
            if (m2.size() > 1)
                AddMM(x2,BandMatrixViewOf(m2.offDiag()),
                      m1.diagRange(-m2.size()+1,0));
            m1.diag().addToAll(x2);
        } else {
            AddMM(x2,BandMatrixViewOf(m2),m1);
        }
    }

    template <bool add, typename T, typename T1, typename T2>
    inline void MultMM(
        const T x, const GenUpperTriMatrix<T1>& m1,
        const GenBandMatrix<T2>& m2, BandMatrixView<T> m0)
    {
        if (m1.isunit()) {
            UpperTriMatrix<T1,NonUnitDiag|RowMajor> m1x = m1;
            MultMM<add>(x,BandMatrixViewOf(m1x),m2,m0);
        } else {
            MultMM<add>(x,BandMatrixViewOf(m1),m2,m0);
        }
    }
    template <bool add, typename T, typename T1, typename T2>
    inline void MultMM(
        const T x, const GenBandMatrix<T1>& m1,
        const GenUpperTriMatrix<T2>& m2, BandMatrixView<T> m0)
    {
        if (m2.isunit()) {
            UpperTriMatrix<T2,NonUnitDiag|RowMajor> m2x = m2;
            MultMM<add>(x,m1,BandMatrixViewOf(m2x),m0);
        } else {
            MultMM<add>(x,m1,BandMatrixViewOf(m2),m0);
        }
    }
    template <bool add, typename T, typename T1, typename T2>
    inline void MultMM(
        const T x, const GenLowerTriMatrix<T1>& m1,
        const GenBandMatrix<T2>& m2, BandMatrixView<T> m0)
    {
        if (m1.isunit()) {
            LowerTriMatrix<T1,NonUnitDiag|RowMajor> m1x = m1;
            MultMM<add>(x,BandMatrixViewOf(m1x),m2,m0);
        } else {
            MultMM<add>(x,BandMatrixViewOf(m1),m2,m0);
        }
    }
    template <bool add, typename T, typename T1, typename T2>
    inline void MultMM(
        const T x, const GenBandMatrix<T1>& m1,
        const GenLowerTriMatrix<T2>& m2, BandMatrixView<T> m0)
    {
        if (m2.isunit()) {
            LowerTriMatrix<T2,NonUnitDiag|RowMajor> m2x = m2;
            MultMM<add>(x,m1,BandMatrixViewOf(m2x),m0);
        } else {
            MultMM<add>(x,m1,BandMatrixViewOf(m2),m0);
        }
    }

    template <typename T, typename T2>
    inline void AddMM(
        const T x2, const GenUpperTriMatrix<T2>& m2,
        BandMatrixView<CT> m1)
    { AddMM(CT(x2),m2,m1); }
    template <typename T, typename T2>
    inline void AddMM(
        const T x2, const GenLowerTriMatrix<T2>& m2,
        BandMatrixView<CT> m1)
    { AddMM(CT(x2),m2,m1); }

    template <typename T, typename T2>
    inline void AddMM(
        const CT , const GenUpperTriMatrix<T2>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename T2>
    inline void AddMM(
        const CT , const GenLowerTriMatrix<T2>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename T1, typename T2>
    inline void MultMM(
        const T x, const GenUpperTriMatrix<T1>& m1,
        const GenBandMatrix<T2>& m2, BandMatrixView<CT> m0)
    { MultMM<add>(CT(x),m1,m2,m0); }
    template <bool add, typename T, typename T1, typename T2>
    inline void MultMM(
        const T x, const GenBandMatrix<T1>& m1,
        const GenUpperTriMatrix<T2>& m2, BandMatrixView<CT> m0)
    { MultMM<add>(CT(x),m1,m2,m0); }
    template <bool add, typename T, typename T1, typename T2>
    inline void MultMM(
        const T x, const GenLowerTriMatrix<T1>& m1,
        const GenBandMatrix<T2>& m2, BandMatrixView<CT> m0)
    { MultMM<add>(CT(x),m1,m2,m0); }
    template <bool add, typename T, typename T1, typename T2>
    inline void MultMM(
        const T x, const GenBandMatrix<T1>& m1,
        const GenLowerTriMatrix<T2>& m2, BandMatrixView<CT> m0)
    { MultMM<add>(CT(x),m1,m2,m0); }

    template <bool add, typename T, typename T1, typename T2>
    inline void MultMM(
        const CT , const GenUpperTriMatrix<T1>& ,
        const GenBandMatrix<T2>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename T1, typename T2>
    inline void MultMM(
        const CT , const GenBandMatrix<T1>& ,
        const GenUpperTriMatrix<T2>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename T1, typename T2>
    inline void MultMM(
        const CT , const GenLowerTriMatrix<T1>& ,
        const GenBandMatrix<T2>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename T1, typename T2>
    inline void MultMM(
        const CT , const GenBandMatrix<T1>& ,
        const GenLowerTriMatrix<T2>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

}

#undef CT

#endif
