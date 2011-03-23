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


#ifndef TMV_MultXD_H
#define TMV_MultXD_H

#include "TMV_BaseMatrix_Diag.h"
#include "TMV_MultXV_Funcs.h"

namespace tmv {

    //
    // D *= x
    //

    template <int ix, class T, class M>
    static inline void Scale(
        const Scaling<ix,T>& x, BaseMatrix_Diag_Mutable<M>& m)
    {
        typename M::diag_type md = m.diag();
        Scale(x,md);
    }


    //
    // D (+)= x * D
    //

    template <bool add, int ix1, class T1, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Diag_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        MultXV<add>(x1,m1d,m2d);
    }

    template <bool add, int ix1, class T1, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Diag_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        NoAliasMultXV<add>(x1,m1d,m2d);
    }

    template <bool add, int ix1, class T1, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Diag_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        AliasMultXV<add>(x1,m1d,m2d);
    }

    //
    // M (+)= x * D
    //

    template <bool add, int ix1, class T1, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        if (SameStorage(m1,m2)) {
            AliasMultXV<add>(x1,m1d,m2d);
            Maybe<!add>::zero_offdiag2(m2.upperTri());
            Maybe<!add>::zero_offdiag2(m2.lowerTri());
        } else {
            Maybe<!add>::zero(m2);
            NoAliasMultXV<add>(x1,m1d,m2d);
        }
    }

    template <bool add, int ix1, class T1, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        Maybe<!add>::zero(m2);
        NoAliasMultXV<add>(x1,m1d,m2d);
    }

    template <bool add, int ix1, class T1, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2)
    { MultXM<add>(x1,m1,m2); }

} // namespace tmv

#endif
