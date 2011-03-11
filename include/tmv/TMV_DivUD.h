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


#ifndef TMV_DivUD_H
#define TMV_DivUD_H

#include "TMV_InvertD.h"
#include "TMV_MultUD.h"

namespace tmv {

    // m3 = x * m2^-1 * m1  (or m3 = x * m1 / m2)
    template <int algo, int ix, class T, class M1, class M2, class M3>
    static void DoLDivUD(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
        TMVStaticAssert((Sizes<M3::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_colsize,M2::_size>::same));
        TMVAssert(m3.colsize() == m1.colsize());
        TMVAssert(m3.rowsize() == m1.rowsize());
        TMVAssert(m3.colsize() == m2.size());
        const int cs1 = Sizes<M3::_colsize,M1::_colsize>::size;
        const int cs = Sizes<cs1,M2::_size>::size;
        const int N = cs == UNKNOWN ? m3.colsize() : cs;

        // First calculate x * m2^-1
        typedef typename M2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        typedef typename MCopyHelper<PT2,Diag,cs,cs,false,false>::type M2c;
        M2c m2c(N);
        typedef typename M2::const_diag_type::const_cview_type M2d;
        typedef typename M2c::diag_type::cview_type M2cd;
        M2d m2d = m2.diag().cView();
        M2cd m2cd = m2c.diag().cView();
        ElemInvert_Helper<-2,cs,ix,T,M2d,M2cd>::call(x,m2d,m2cd);

        typedef typename M1::const_transpose_type::const_cview_type M1v;
        typedef typename M2c::const_cview_type M2cv;
        typedef typename M3::const_transpose_type::cview_type M3v;
        M1t m1t = m1.transpose().cView();
        M2cv m2cv = m2c.cView();
        M3t m3t = m3.transpose().cView();
        typedef typename Traits<T>::real_type RT;
        Scaling<1,RT> one;
        MultUD_Helper<algo,cs,rs,false,1,RT,M1t,M2cv,M3t>::call(one,m1t,m2cv,m3t);
    }

    template <int ix, class T, class M1, class M2, class M3>
    static void LDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    { DoLDivUD<-1>(x,m1,m2,m3); }

    template <int ix, class T, class M1, class M2, class M3>
    static void NoAliasLDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    { DoLDivUD<-2>(x,m1,m2,m3); }

    template <int ix, class T, class M1, class M2, class M3>
    static void AliasLDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    { DoLDivUD<99>(x,m1,m2,m3); }

    template <int ix, class T, class M1, class M2, class M3>
    static void LDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M3::uppertri_type,
                typename M3::lowertri_type>::type M3u;
        M3u m3u = Maybe<upper>::uppertri(m3);
        LDiv(x,m1,m2,m3u)
        Maybe<!upper>::uppertri(m3).offDiag().setZero();
    }

    template <int ix, class T, class M1, class M2, class M3>
    static void NoAliasLDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { 
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M3::uppertri_type,
                typename M3::lowertri_type>::type M3u;
        M3u m3u = Maybe<upper>::uppertri(m3);
        m3.setZero();
        NoAliasLDiv(x,m1,m2,m3u)
    }

    template <int ix, class T, class M1, class M2, class M3>
    static void AliasLDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M3::uppertri_type,
                typename M3::lowertri_type>::type M3u;
        M3u m3u = Maybe<upper>::uppertri(m3);
        AliasLDiv(x,m1,m2,m3u)
        Maybe<!upper>::uppertri(m3).offDiag().setZero();
    }


    template <class M1, int ix, class T, class M2>
    static void LDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { DoLDivUD<-1>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    static void NoAliasLDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { DoLDivUD<-2>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    static void AliasLDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { DoLDivUD<99>(x,m1.mat(),m2.mat(),m1.mat()); }



    // m3 = x * m1 * m2^-1 (or m3 = x * m1 % m2)
    template <int algo, int ix, class T, class M1, class M2, class M3>
    static void DoRDivUD(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
        TMVStaticAssert((Sizes<M3::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_rolsize,M2::_size>::same));
        TMVAssert(m3.colsize() == m1.colsize());
        TMVAssert(m3.rowsize() == m1.rowsize());
        TMVAssert(m3.rolsize() == m2.size());
        const int rs1 = Sizes<M3::_rowsize,M1::_rowsize>::size;
        const int rs = Sizes<rs1,M2::_size>::size;
        const int N = rs == UNKNOWN ? m3.rowsize() : rs;

        // First calculate x * m2^-1
        typedef typename M2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        typedef typename MCopyHelper<PT2,Diag,rs,rs,false,false>::type M2c;
        M2c m2c(N);
        typedef typename M2::const_diag_type::const_cview_type M2d;
        typedef typename M2c::diag_type::cview_type M2cd;
        M2d m2d = m2.diag().cView();
        M2cd m2cd = m2c.diag().cView();
        ElemInvert_Helper<-2,rs,ix,T,M2d,M2cd>::call(x,m2d,m2cd);

        typedef typename M1::const_cview_type M1v;
        typedef typename M2c::const_cview_type M2cv;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2cv m2cv = m2c.cView();
        M3v m3v = m3.cView();
        typedef typename Traits<T>::real_type RT;
        Scaling<1,RT> one;
        MultUD_Helper<algo,cs,rs,false,1,RT,M1v,M2cv,M3v>::call(one,m1v,m2cv,m3v);
    }

    template <int ix, class T, class M1, class M2, class M3>
    static void RDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    { DoRDivUD<-1>(x,m1,m2,m3); }

    template <int ix, class T, class M1, class M2, class M3>
    static void NoAliasRDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    { DoRDivUD<-2>(x,m1,m2,m3); }

    template <int ix, class T, class M1, class M2, class M3>
    static void AliasRDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    { DoRDivUD<99>(x,m1,m2,m3); }

    template <int ix, class T, class M1, class M2, class M3>
    static void RDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M3::uppertri_type,
                typename M3::lowertri_type>::type M3u;
        M3u m3u = Maybe<upper>::uppertri(m3);
        RDiv(x,m1,m2,m3u)
        Maybe<!upper>::uppertri(m3).offDiag().setZero();
    }

    template <int ix, class T, class M1, class M2, class M3>
    static void NoAliasRDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { 
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M3::uppertri_type,
                typename M3::lowertri_type>::type M3u;
        M3u m3u = Maybe<upper>::uppertri(m3);
        m3.setZero();
        NoAliasRDiv(x,m1,m2,m3u)
    }

    template <int ix, class T, class M1, class M2, class M3>
    static void AliasRDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M3::uppertri_type,
                typename M3::lowertri_type>::type M3u;
        M3u m3u = Maybe<upper>::uppertri(m3);
        AliasRDiv(x,m1,m2,m3u)
        Maybe<!upper>::uppertri(m3).offDiag().setZero();
    }


    template <class M1, int ix, class T, class M2>
    static void RDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { DoRDivUD<-1>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    static void NoAliasRDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { DoRDivUD<-2>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    static void AliasRDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { DoRDivUD<99>(x,m1.mat(),m2.mat(),m1.mat()); }

} // namespace tmv

#undef TMV_ZeroIX

#endif 
