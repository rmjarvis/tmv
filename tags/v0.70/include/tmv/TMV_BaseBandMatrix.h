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


#ifndef TMV_BaseBandMatrix_H
#define TMV_BaseBandMatrix_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    template <class T> 
    class GenBandMatrix;

    template <class T, int A=0>
    class ConstBandMatrixView;

    template <class T, int A=0>
    class BandMatrixView;

    template <class T, int A=0>
    class BandMatrix;

    template <class T> 
    class BandLUDiv;

    template <class T> 
    class BandSVDiv;

    template <class T> 
    class BandQRDiv;

    template <class T, class Tm> 
    class QuotXB;

    ptrdiff_t BandStorageLength(StorageType s, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, ptrdiff_t hi);
    ptrdiff_t BandNumElements(ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, ptrdiff_t hi);

    template <class T1, class T2> 
    void Copy(const GenBandMatrix<T1>& m1, BandMatrixView<T2> m2);

    template <class T> 
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

    template <class T, int A>
    inline std::string TMV_Text(const BandMatrix<T,A>& )
    { 
        return std::string("BandMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <class T> 
    inline std::string TMV_Text(const GenBandMatrix<T>& )
    {
        return std::string("GenBandMatrix<") + TMV_Text(T()) + ">";
    }
    template <class T, int A>
    inline std::string TMV_Text(const ConstBandMatrixView<T,A>& )
    {
        return std::string("ConstBandMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }
    template <class T, int A>
    inline std::string TMV_Text(const BandMatrixView<T,A>& )
    {
        return std::string("BandMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

} // namespace tmv

#endif
