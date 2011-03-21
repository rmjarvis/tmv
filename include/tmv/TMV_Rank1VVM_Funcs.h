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


#ifndef TMV_Rank1VVM_Funcs_H
#define TMV_Rank1VVM_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //

    // From TMV_Rank1_VVM.h:
    template <bool add, int ix, class T, class V1, class V2, class M3>
    static void Rank1Update(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class V1, class V2, class M3>
    static void NoAliasRank1Update(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class V1, class V2, class M3>
    static void AliasRank1Update(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseMatrix_Rec_Mutable<M3>& m3);

} // namespace tmv

#endif 
