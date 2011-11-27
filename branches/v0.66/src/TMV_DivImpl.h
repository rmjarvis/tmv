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


#ifndef TMV_DIVIMPL_H
#define TMV_DIVIMPL_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    template <class T> 
    struct DivHelper<T>::DivImpl 
    {
    public :

        DivImpl(const BaseMatrix<T>& _m) : 
            m(_m), div(0), dt(m.isSquare() ? LU : QR),
            inplace(false), cache(false) {}
        ~DivImpl() {}

        const BaseMatrix<T>& m;
        mutable std::auto_ptr<Divider<T> > div;
        DivType dt;
        bool inplace;
        bool cache;

    private :

        DivImpl(const DivImpl&);
        DivImpl& operator=(const DivImpl&);
    };

}

#endif

