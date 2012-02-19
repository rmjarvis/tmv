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


#ifndef TMV_BaseSymBandMatrix_H
#define TMV_BaseSymBandMatrix_H

#include "tmv/TMV_BaseSymMatrix.h"
#include "tmv/TMV_BaseBandMatrix.h"

namespace tmv {

    template <class T> 
    class GenSymBandMatrix;

    template <class T, int A=0>
    class ConstSymBandMatrixView;

    template <class T, int A=0>
    class SymBandMatrixView;

    template <class T, int A=0>
    class SymBandMatrix;

    template <class T, int A=0>
    class HermBandMatrix;

    template <class T> 
    class HermBandCHDiv;

    template <class T> 
    class HermBandSVDiv;

    template <class T> 
    class SymBandSVDiv;

    template <class T, class Tm> 
    class QuotXsB;

    template <class T, class T1, class T2> 
    class ProdBB;

    template <class T> 
    struct AssignableToSymBandMatrix :
        virtual public AssignableToSymMatrix<T>,
        virtual public AssignableToBandMatrix<T>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        virtual void assignTosB(SymBandMatrixView<RT> m) const = 0;
        virtual void assignTosB(SymBandMatrixView<CT> m) const = 0;
        virtual inline ~AssignableToSymBandMatrix() {}
    };


    template <class T, int A>
    inline std::string TMV_Text(const SymBandMatrix<T,A>& )
    {
        return std::string("SymBandMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <class T, int A>
    inline std::string TMV_Text(const HermBandMatrix<T,A>& )
    {
        return std::string("HermBandMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <class T> 
    inline std::string TMV_Text(const GenSymBandMatrix<T>& )
    {
        return std::string("GenSymBandMatrix<") + TMV_Text(T()) + ">";
    }

    template <class T, int A>
    inline std::string TMV_Text(const ConstSymBandMatrixView<T,A>& )
    {
        return std::string("ConstSymBandMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <class T, int A>
    inline std::string TMV_Text(const SymBandMatrixView<T,A>& )
    {
        return std::string("SymBandMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

}

#endif
