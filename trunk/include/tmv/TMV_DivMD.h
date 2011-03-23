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
// Boston, MA  02320-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef TMV_DivMD_H
#define TMV_DivMD_H

#include "TMV_MultMM_Funcs.h"

namespace tmv {

    template <int ix, class T, class M1, class M2, class M3>
    static inline void LDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { MultMM<false>(x,m2.inverse().calc(),m1.mat(),m3.mat()); }
    template <int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { NoAliasMultMM<false>(x,m2.inverse().calc(),m1.mat(),m3.mat()); }
    template <int ix, class T, class M1, class M2, class M3>
    static inline void AliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { AliasMultMM<false>(x,m2.inverse().calc(),m1.mat(),m3.mat()); }

    template <class M1, class M2>
    static inline void LDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Diag<M2>& m2)
    { LDiv(Scaling<1,typename M2::real_type>(),m1.mat(),m2,m1); }
    template <class M1, class M2>
    static inline void NoAliasLDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Diag<M2>& m2)
    { NoAliasLDiv(Scaling<1,typename M2::real_type>(),m1.mat(),m2,m1); }
    template <class M1, class M2>
    static inline void AliasLDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Diag<M2>& m2)
    { AliasLDiv(Scaling<1,typename M2::real_type>(),m1.mat(),m2,m1); }

    template <int ix, class T, class M1, class M2, class M3>
    static inline void RDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { MultMM<false>(x,m1.mat(),m2.inverse().calc(),m3.mat()); }
    template <int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { NoAliasMultMM<false>(x,m1.mat(),m2.inverse().calc(),m3.mat()); }
    template <int ix, class T, class M1, class M2, class M3>
    static inline void AliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { AliasMultMM<false>(x,m1.mat(),m2.inverse().calc(),m3.mat()); }

    template <class M1, class M2>
    static inline void RDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Diag<M2>& m2)
    { RDiv(Scaling<1,typename M2::real_type>(),m1.mat(),m2,m1); }
    template <class M1, class M2>
    static inline void NoAliasRDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Diag<M2>& m2)
    { NoAliasRDiv(Scaling<1,typename M2::real_type>(),m1.mat(),m2,m1); }
    template <class M1, class M2>
    static inline void AliasRDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Diag<M2>& m2)
    { AliasRDiv(Scaling<1,typename M2::real_type>(),m1.mat(),m2,m1); }

} // namespace tmv

#undef TMV_ZeroIX

#endif 
