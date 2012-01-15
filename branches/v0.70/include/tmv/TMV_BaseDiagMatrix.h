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


#ifndef TMV_BaseDiagMatrix_H
#define TMV_BaseDiagMatrix_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    template <class T> 
    class GenDiagMatrix;

    template <class T, int A=0>
    class ConstDiagMatrixView;

    template <class T, int A=0>
    class DiagMatrixView;

    template <class T, int A=0>
    class DiagMatrix;

    template <class T, class Tm> 
    class QuotXD;

    template <class T> 
    struct AssignableToDiagMatrix :
        virtual public AssignableToMatrix<T>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        virtual int size() const = 0;
        virtual void assignToD(const DiagMatrixView<RT>& m) const = 0;
        virtual void assignToD(const DiagMatrixView<CT>& m) const = 0;
        virtual inline ~AssignableToDiagMatrix() {}
    };

    template <class T, int A>
    inline std::string TMV_Text(const DiagMatrix<T,A>& )
    {
        return std::string("DiagMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::vtext() + ">";
    }

    template <class T> 
    inline std::string TMV_Text(const GenDiagMatrix<T>& )
    {
        return std::string("GenDiagMatrix<") + TMV_Text(T()) + ">";
    }

    template <class T, int A>
    inline std::string TMV_Text(const ConstDiagMatrixView<T,A>& )
    {
        return std::string("ConstDiagMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::vtext() + ">";
    }

    template <class T, int A>
    inline std::string TMV_Text(const DiagMatrixView<T,A>& )
    {
        return std::string("DiagMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::vtext() + ">";
    }

} // namespace tmv

#endif
