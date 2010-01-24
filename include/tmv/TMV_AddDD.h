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


#ifndef TMV_AddDD_H
#define TMV_AddDD_H

#include "TMV_BaseMatrix_Diag.h"
#include "TMV_AddVV.h"
#include "TMV_MultXM.h"
#include "TMV_MultXV.h"

namespace tmv {

    //
    // D = x * D + x * D
    //

    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Diag_Mutable<M3>& m3)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::const_diag_type m2d = m2.diag();
        typename M3::diag_type m3d = m3.diag();
        AddVV(x1,m1d,x2,m2d,m3d); 
    }
    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Diag_Mutable<M3>& m3)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::const_diag_type m2d = m2.diag();
        typename M3::diag_type m3d = m3.diag();
        NoAliasAddVV(x1,m1d,x2,m2d,m3d); 
    }
    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Diag_Mutable<M3>& m3)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::const_diag_type m2d = m2.diag();
        typename M3::diag_type m3d = m3.diag();
        AliasAddVV(x1,m1d,x2,m2d,m3d); 
    }

    //
    // M = x * D + x * D
    //

    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::const_diag_type m2d = m2.diag();
        typename M3::diag_type m3d = m3.diag();
        AddVV(x1,m1d,x2,m2d,m3d); 
        m3.upperTri().offDiag().setZero();
        m3.lowerTri().offDiag().setZero();
    }
    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::const_diag_type m2d = m2.diag();
        typename M3::diag_type m3d = m3.diag();
        m3.setZero();
        NoAliasAddVV(x1,m1d,x2,m2d,m3d); 
    }
    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::const_diag_type m2d = m2.diag();
        typename M3::diag_type m3d = m3.diag();
        AliasAddVV(x1,m1d,x2,m2d,m3d); 
        m3.upperTri().offDiag().setZero();
        m3.lowerTri().offDiag().setZero();
    }

    //
    // M = x * D + x * M
    //

    // TODO: This function is supposed to only check for aliasing if 
    // the sizes are UNKNOWN.
    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Calc<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (SameStorage(m1,m3)) {
            typename M1::copy_type m1c = m1;
            typename M1::copy_type::const_diag_type m1d = m1c.diag();
            typename M3::diag_type m3d = m3.diag();
            MultXM<false>(x2,m2.mat(),m3.mat());
            NoAliasMultXV<true>(x1,m1d,m3d);
        } else {
            typename M1::const_diag_type m1d = m1.diag();
            typename M3::diag_type m3d = m3.diag();
            MultXM<false>(x2,m2.mat(),m3.mat());
            NoAliasMultXV<true>(x1,m1d,m3d);
        } 
    }
    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Calc<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M3::diag_type m3d = m3.diag();
        NoAliasMultXM<false>(x2,m2.mat(),m3.mat());
        NoAliasMultXV<true>(x1,m1d,m3d);
    }
    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Calc<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (SameStorage(m1,m3)) {
            typename M1::copy_type m1c = m1;
            typename M1::copy_type::const_diag_type m1d = m1c.diag();
            typename M3::diag_type m3d = m3.diag();
            MultXM<false>(x2,m2.mat(),m3.mat());
            NoAliasMultXV<true>(x1,m1d,m3d);
        } else {
            typename M1::const_diag_type m1d = m1.diag();
            typename M3::diag_type m3d = m3.diag();
            MultXM<false>(x2,m2.mat(),m3.mat());
            NoAliasMultXV<true>(x1,m1d,m3d);
        } 
    }

    //
    // M = x * M + x * D
    //

    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Calc<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    { AddMM(x2,m2,x1,m1,m3); }
    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Calc<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    { NoAliasAddMM(x2,m2,x1,m1,m3); }
    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Calc<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    { AliasAddMM(x2,m2,x1,m1,m3); }

} // namespace tmv

#endif
