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


#ifndef TMV_BaseBandMatrix_H
#define TMV_BaseBandMatrix_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    template <typename T>
    class GenBandMatrix;

    template <typename T, int A=0>
    class ConstBandMatrixView;

    template <typename T, int A=0>
    class BandMatrixView;

    template <typename T, int A=0>
    class BandMatrix;

    template <typename T>
    class BandLUDiv;

    template <typename T>
    class BandSVDiv;

    template <typename T>
    class BandQRDiv;

    template <typename T, typename Tm>
    class QuotXB;

    ptrdiff_t BandStorageLength(StorageType s, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, ptrdiff_t hi);
    ptrdiff_t BandNumElements(ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, ptrdiff_t hi);

    template <typename T1, typename T2>
    void Copy(const GenBandMatrix<T1>& m1, BandMatrixView<T2> m2);

    template <typename T>
    struct AssignableToBandMatrix : virtual public AssignableToMatrix<T>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        virtual ptrdiff_t nlo() const = 0;
        virtual ptrdiff_t nhi() const = 0;
        virtual void assignToB(BandMatrixView<RT> m) const = 0;
        virtual void assignToB(BandMatrixView<CT> m) const = 0;
        virtual inline ~AssignableToBandMatrix() {}
    };

    template <typename T, int A>
    inline std::string TMV_Text(const BandMatrix<T,A>& )
    {
        return std::string("BandMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <typename T>
    inline std::string TMV_Text(const GenBandMatrix<T>& )
    {
        return std::string("GenBandMatrix<") + TMV_Text(T()) + ">";
    }
    template <typename T, int A>
    inline std::string TMV_Text(const ConstBandMatrixView<T,A>& )
    {
        return std::string("ConstBandMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }
    template <typename T, int A>
    inline std::string TMV_Text(const BandMatrixView<T,A>& )
    {
        return std::string("BandMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

} // namespace tmv

#endif
