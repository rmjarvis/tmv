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


#ifndef TMV_MultMM_Funcs_H
#define TMV_MultMM_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //


    // From TMV_MultMM.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <class M1, int ix, class T, class M2>
    static void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static void NoAliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static void AliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M2>& m2);

    // From TMV_MultMM_Block.h
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM_RecursiveBlock(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM_Block(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    // From TMV_MultMM_Winograd.h
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM_Winograd(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    // From TMV_MultMM_OpenMP
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM_OpenMP(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    // From TMV_MultMD.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <class M1, int ix, class T, class M2>
    static void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static void NoAliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static void AliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);

    // From TMV_MultDV.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Diag_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Diag_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Diag_Mutable<M3>& m3);

    template <class M1, int ix, class T, class M2>
    static void MultEqMM(
        BaseMatrix_Diag_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static void NoAliasMultEqMM(
        BaseMatrix_Diag_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static void AliasMultEqMM(
        BaseMatrix_Diag_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    // From TMV_MultUM.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <class M1, int ix, class T, class M2>
    static void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static void NoAliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static void AliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);

    // From TMV_MultUU.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3);

    template <class M1, int ix, class T, class M2>
    static void MultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static void NoAliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static void AliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);

    // From TMV_MultUL.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    // From TMV_MultUD.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3);

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Mutable<M3>& m3);

    template <class M1, int ix, class T, class M2>
    static void MultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static void NoAliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static void AliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);

} // namespace tmv

#endif
