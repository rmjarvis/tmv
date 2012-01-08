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


#ifndef TMV_BaseTriMatrix_H
#define TMV_BaseTriMatrix_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    template <class T> 
    class GenUpperTriMatrix;

    template <class T> 
    class GenLowerTriMatrix;

    template <class T, int A=0>
    class ConstUpperTriMatrixView;

    template <class T, int A=0>
    class ConstLowerTriMatrixView;

    template <class T, int A=0>
    class UpperTriMatrixView;

    template <class T, int A=0>
    class LowerTriMatrixView;

    template <class T, int A=0>
    class UpperTriMatrix;

    template <class T, int A=0>
    class LowerTriMatrix;

    template <class T> 
    class UpperTriDiv;

    template <class T> 
    class LowerTriDiv;

    template <class T, class Tm> 
    class QuotXU;

    template <class T, class Tm> 
    class QuotXL;

    template <class T1, class T2> 
    void Copy(
        const GenUpperTriMatrix<T1>& m1, const UpperTriMatrixView<T2>& m2);

    template <class T> 
    struct AssignableToUpperTriMatrix :
        virtual public AssignableToMatrix<T>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        virtual int size() const = 0;
        virtual DiagType dt() const = 0;
        virtual void assignToU(const UpperTriMatrixView<RT>& m) const = 0;
        virtual void assignToU(const UpperTriMatrixView<CT>& m) const = 0;
        virtual inline ~AssignableToUpperTriMatrix() {}
    };

    template <class T> 
    struct AssignableToLowerTriMatrix :
        virtual public AssignableToMatrix<T>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        virtual int size() const = 0;
        virtual DiagType dt() const = 0;
        virtual void assignToL(const LowerTriMatrixView<RT>& m) const = 0;
        virtual void assignToL(const LowerTriMatrixView<CT>& m) const = 0;
        virtual inline ~AssignableToLowerTriMatrix() {}
    };

    template <class T, int A>
    inline std::string TMV_Text(const UpperTriMatrix<T,A>& )
    {
        return std::string("UpperTriMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }
    template <class T, int A>
    inline std::string TMV_Text(const LowerTriMatrix<T,A>& )
    {
        return std::string("LowerTriMatrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <class T> 
    inline std::string TMV_Text(const GenUpperTriMatrix<T>& m)
    {
        return std::string("GenUpperTriMatrix<") + TMV_Text(T()) + ">";
    }
    template <class T> 
    inline std::string TMV_Text(const GenLowerTriMatrix<T>& m)
    {
        return std::string("GenLowerTriMatrix<") + TMV_Text(T()) + ">";
    }

    template <class T, int A>
    inline std::string TMV_Text(
        const ConstUpperTriMatrixView<T,A>& m)
    {
        return std::string("ConstUpperTriMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }
    template <class T, int A>
    inline std::string TMV_Text(
        const ConstLowerTriMatrixView<T,A>& m)
    {
        return std::string("ConstLowerTriMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

    template <class T, int A>
    inline std::string TMV_Text(
        const UpperTriMatrixView<T,A>& m)
    {
        return std::string("UpperTriMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }
    template <class T, int A>
    inline std::string TMV_Text(
        const LowerTriMatrixView<T,A>& m)
    {
        return std::string("LowerTriMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

} // namespace tmv

#endif
