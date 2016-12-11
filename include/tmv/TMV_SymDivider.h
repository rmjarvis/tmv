//////////////////////////////////////////////////////////////////////////////
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


//---------------------------------------------------------------------------
//
// This file defines the TMV SymDivider class.
//
// It adds Inverse(SymMatrix) capability to the regular Divier class.
//

#ifndef TMV_SymDivider_H
#define TMV_SymDivider_H

#include "tmv/TMV_Divider.h"
#include "tmv/TMV_BaseSymMatrix.h"

namespace tmv {

#define RT TMV_RealType(T)
#define CT TMV_ComplexType(T)

    template <typename T>
    class SymDivider : public Divider<T>
    {

    public :

        SymDivider() {}
        virtual ~SymDivider() {}

        using Divider<T>::makeInverse;
        virtual void makeInverse(SymMatrixView<RT> sinv) const =0;
        virtual void makeInverse(SymMatrixView<CT> sinv) const =0;
    };

    template <typename T>
    inline std::string TMV_Text(const SymDivider<T>& d)
    { return std::string("SymDivider<")+tmv::TMV_Text(T())+">"; }

#undef RT
#undef CT

} // namespace tmv

#endif
