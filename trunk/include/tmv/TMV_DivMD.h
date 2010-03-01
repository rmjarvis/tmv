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

#include "TMV_InvertD.h"
#include "TMV_MultMD.h"

namespace tmv {

    // m3 = x * m2^-1 * m1  (or m3 = x * m1 / m2)
    template <int algo, int ix, class T, class M1, class M2, class M3>
    inline void DoLDivMD(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
        TMVStaticAssert((Sizes<M3::mcolsize,M1::mcolsize>::same));
        TMVStaticAssert((Sizes<M3::mrowsize,M1::mrowsize>::same));
        TMVStaticAssert((Sizes<M3::mcolsize,M2::msize>::same));
        TMVAssert(m3.colsize() == m1.colsize());
        TMVAssert(m3.rowsize() == m1.rowsize());
        TMVAssert(m3.colsize() == m2.size());
        const int cs1 = Sizes<M3::mcolsize,M1::mcolsize>::size;
        const int cs = Sizes<cs1,M2::msize>::size;
        const int N = cs == UNKNOWN ? m3.colsize() : cs;
        const int rs = Sizes<M3::mrowsize,M1::mrowsize>::size;

        // First calculate x * m2^-1
        typedef typename M2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        typedef typename MCopyHelper<PT2,Diag,cs,cs,false,false>::type M2c;
        M2c m2c(N);
        NoAliasInvert(x,m2,m2c);

        typedef typename M1::const_transpose_type::const_cview_type M1t;
        typedef typename M2c::const_cview_type M2cv;
        typedef typename M3::const_transpose_type::cview_type M3t;
        M1t m1t = m1.transpose().cView();
        M2cv m2cv = m2c.cView();
        M3t m3t = m3.transpose().cView();
        typedef typename Traits<T>::real_type RT;
        Scaling<1,RT> one;
        MultMD_Helper<algo,cs,rs,false,1,RT,M1t,M2cv,M3t>::call(one,m1t,m2cv,m3t);
    }

    template <int ix, class T, class M1, class M2, class M3>
    inline void LDivMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoLDivMD<-1>(x,m1,m2,m3); }

    template <int ix, class T, class M1, class M2, class M3>
    inline void NoAliasLDivMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoLDivMD<-2>(x,m1,m2,m3); }

    template <int ix, class T, class M1, class M2, class M3>
    inline void AliasLDivMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoLDivMD<99>(x,m1,m2,m3); }

    template <class M1, int ix, class T, class M2>
    inline void LDivEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { DoLDivMD<-1>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    inline void NoAliasLDivEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { DoLDivMD<-2>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    inline void AliasLDivEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { DoLDivMD<99>(x,m1.mat(),m2.mat(),m1.mat()); }

    // m3 = x * m1 * m2^-1 (or m3 = x * m1 % m2)
    template <int algo, int ix, class T, class M1, class M2, class M3>
    inline void DoRDivMD(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
        TMVStaticAssert((Sizes<M3::mcolsize,M1::mcolsize>::same));
        TMVStaticAssert((Sizes<M3::mrowsize,M1::mrowsize>::same));
        TMVStaticAssert((Sizes<M3::mrolsize,M2::msize>::same));
        TMVAssert(m3.colsize() == m1.colsize());
        TMVAssert(m3.rowsize() == m1.rowsize());
        TMVAssert(m3.rolsize() == m2.size());
        const int rs1 = Sizes<M3::mrowsize,M1::mrowsize>::size;
        const int rs = Sizes<rs1,M2::msize>::size;
        const int cs = Sizes<M3::mcolsize,M1::mcolsize>::size;
        const int N = rs == UNKNOWN ? m3.rowsize() : rs;

        // First calculate x * m2^-1
        typedef typename M2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        typedef typename MCopyHelper<PT2,Diag,rs,rs,false,false>::type M2c;
        M2c m2c(N);
        NoAliasInvert(x,m2,m2c);

        typedef typename M1::const_cview_type M1v;
        typedef typename M2c::const_cview_type M2cv;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2cv m2cv = m2c.cView();
        M3v m3v = m3.cView();
        typedef typename Traits<T>::real_type RT;
        Scaling<1,RT> one;
        MultMD_Helper<algo,cs,rs,false,1,RT,M1v,M2cv,M3v>::call(
            one,m1v,m2cv,m3v);
    }

    template <int ix, class T, class M1, class M2, class M3>
    inline void RDivMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoRDivMD<-1>(x,m1,m2,m3); }

    template <int ix, class T, class M1, class M2, class M3>
    inline void NoAliasRDivMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoRDivMD<-2>(x,m1,m2,m3); }

    template <int ix, class T, class M1, class M2, class M3>
    inline void AliasRDivMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoRDivMD<99>(x,m1,m2,m3); }

    template <class M1, int ix, class T, class M2>
    inline void RDivEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { DoRDivMD<-1>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    inline void NoAliasRDivEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { DoRDivMD<-2>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    inline void AliasRDivEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { DoRDivMD<99>(x,m1.mat(),m2.mat(),m1.mat()); }

} // namespace tmv

#undef TMV_ZeroIX

#endif 
