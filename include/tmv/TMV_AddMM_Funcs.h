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


#ifndef TMV_AddMM_Funcs_H
#define TMV_AddMM_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //

    // From TMV_AddMM.h:
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);

    // From TMV_AddDD.h:
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Diag_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Diag_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Diag_Mutable<M3>& m3);

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Mutable<M3>& m3);

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Calc<M2>& m2,
        BaseMatrix_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Calc<M2>& m2,
        BaseMatrix_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Calc<M2>& m2,
        BaseMatrix_Mutable<M3>& m3);

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Calc<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Calc<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Calc<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Mutable<M3>& m3);

    // From TMV_AddUU.h:
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3);

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);

} // namespace tmv

#endif 
