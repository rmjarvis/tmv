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


#ifndef TMV_ProdVV_H
#define TMV_ProdVV_H

#include "TMV_ProdXV.h"

namespace tmv {

    //
    // Vector * Vector
    //

#define PT typename ProdType<V1,V2>::type

    // These are defined in TMV_MultVV.h
    template <class V1, class V2>
    inline PT MultVV(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2);
    template <class V1, class V2>
    inline PT NoAliasMultVV(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2);
    template <class V1, class V2>
    inline PT InlineMultVV(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2);
    template <class V1, class V2>
    inline PT AliasMultVV(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2);

    // v * v
    template <class V1, class V2>
    inline PT operator*(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2) 
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        return MultVV(v1.calc(),v2.calc()); 
    }

#define PT2 typename Traits2<Tx,PT>::type
    // v * (x*v)
    template <class V1, int ix2, class Tx, class V2>
    inline PT2 operator*(const BaseVector<V1>& v1, const ProdXV<ix2,Tx,V2>& v2) 
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        return v2.getX() * (v1 * v2.getV());
    }

    // (x*v) * v
    template <int ix1, class Tx, class V1, class V2>
    inline PT2 operator*(const ProdXV<ix1,Tx,V1>& v1, const BaseVector<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        return v1.getX() * (v1.getV() * v2);
    }
#undef PT2

#define PT2 typename Traits2<Tx1,typename Traits2<Tx2,PT>::type>::type
    // (x*v) * (x*v)
    template <int ix1, class Tx1, class V1, int ix2, class Tx2, class V2>
    inline PT2 operator*(
        const ProdXV<ix1,Tx1,V1>& v1, const ProdXV<ix2,Tx2,V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        return v1.getX() * v2.getX() * (v1.getV() * v2.getV());
    }
#undef PT2

#undef PT

} // namespace tmv

#endif 
