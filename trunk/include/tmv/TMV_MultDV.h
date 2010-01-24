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


#ifndef TMV_MultDV_H
#define TMV_MultDV_H

#include "TMV_BaseMatrix_Diag.h"
#include "TMV_ElemMultVV.h"

namespace tmv {

    // 
    // MultDV
    //

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void MultMV(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseVector_Calc<V2>& v2,
        BaseVector_Mutable<V3>& v3)
    { ElemMultVV<add>(x,m1.diag(),v2.vec(),v3.vec()); }

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void NoAliasMultMV(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseVector_Calc<V2>& v2,
        BaseVector_Mutable<V3>& v3)
    { NoAliasElemMultVV<add>(x,m1.diag(),v2.vec(),v3.vec()); }

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void AliasMultMV(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseVector_Calc<V2>& v2,
        BaseVector_Mutable<V3>& v3)
    { AliasElemMultVV<add>(x,m1.diag(),v2.vec(),v3.vec()); }

    template <bool add, int ix, class T, class V1, class M2, class V3>
    inline void MultVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1, const BaseMatrix_Diag<M2>& m2,
        BaseVector_Mutable<V3>& v3)
    { ElemMultVV<add>(x,m2.diag(),v1.vec(),v3.vec()); }

    template <bool add, int ix, class T, class V1, class M2, class V3>
    inline void NoAliasMultVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1, const BaseMatrix_Diag<M2>& m2,
        BaseVector_Mutable<V3>& v3)
    { NoAliasElemMultVV<add>(x,m2.diag(),v1.vec(),v3.vec()); }

    template <bool add, int ix, class T, class V1, class M2, class V3>
    inline void AliasMultVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1, const BaseMatrix_Diag<M2>& m2,
        BaseVector_Mutable<V3>& v3)
    { AliasElemMultVV<add>(x,m2.diag(),v1.vec(),v3.vec()); }

    template <class V1, int ix, class T, class M2>
    inline void MultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { ElemMultVV<false>(x,m2.diag(),v1.vec(),v1.vec()); }

    template <class V1, int ix, class T, class M2>
    inline void NoAliasMultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { NoAliasElemMultVV<false>(x,m2.diag(),v1.vec(),v1.vec()); }

    template <class V1, int ix, class T, class M2>
    inline void AliasMultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { AliasElemMultVV<false>(x,m2.diag(),v1.vec(),v1.vec()); }


    //
    // MultDD
    //

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Diag_Mutable<M3>& m3)
    { 
        typename M3::diag_type m3d = m3.diag();
        ElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Diag_Mutable<M3>& m3)
    { 
        typename M3::diag_type m3d = m3.diag();
        NoAliasElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Diag_Mutable<M3>& m3)
    { 
        typename M3::diag_type m3d = m3.diag();
        AliasElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
    }

    template <class M1, int ix, class T, class M2>
    inline void MultEqMM(
        BaseMatrix_Diag_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { 
        typename M1::diag_type m1d = m1.diag();
        ElemMultVV<false>(x,m1.diag(),m2.diag(),m1d); 
    }

    template <class M1, int ix, class T, class M2>
    inline void NoAliasMultEqMM(
        BaseMatrix_Diag_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { 
        typename M1::diag_type m1d = m1.diag();
        NoAliasElemMultVV<false>(x,m1.diag(),m2.diag(),m1d); 
    }

    template <class M1, int ix, class T, class M2>
    inline void AliasMultEqMM(
        BaseMatrix_Diag_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { 
        typename M1::diag_type m1d = m1.diag();
        AliasElemMultVV<false>(x,m1.diag(),m2.diag(),m1d); 
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M3::diag_type m3d = m3.diag();
        if (!add && (SameStorage(m1,m3) || SameStorage(m2,m3))) {
            ElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
            m3.upperTri().offDiag().setZero();
            m3.lowerTri().offDiag().setZero();
        } else {
            Maybe<!add>::zero(m3);
            ElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
        }
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M3::diag_type m3d = m3.diag();
        Maybe<!add>::zero(m3);
        ElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M3::diag_type m3d = m3.diag();
        if (!add && (SameStorage(m1,m3) || SameStorage(m2,m3))) {
            AliasElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
            m3.upperTri().offDiag().setZero();
            m3.lowerTri().offDiag().setZero();
        } else {
            Maybe<!add>::zero(m3);
            ElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
        }
    }

} // namespace tmv

#endif
