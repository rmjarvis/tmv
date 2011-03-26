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


#ifndef TMV_DivVM_Funcs_H
#define TMV_DivVM_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //

    // From TMV_DivM.h:
    template <class V1, class M2>
    static void LDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2);
    template <class V1, class M2>
    static void NoAliasLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2);
    template <class V1, class M2>
    static void AliasLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2);

    template <class V1, class M2>
    static void RDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2);
    template <class V1, class M2>
    static void NoAliasRDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2);
    template <class V1, class M2>
    static void AliasRDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2);

    template <int ix, class T, class V1, class M2, class V3>
    static void LDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    static void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    static void AliasLDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);

    template <int ix, class T, class V1, class M2, class V3>
    static void RDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    static void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    static void AliasRDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);

    // From TMV_DivVU.h:
    template <class V1, class M2>
    static void TriLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Tri<M2>& m2);
    template <class V1, class M2>
    static void NoAliasTriLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Tri<M2>& m2);
    template <class V1, class M2>
    static void AliasTriLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Tri<M2>& m2);

    // From TMV_MultPV.h:
    template <int ix, class T, class V1, class V3>
    static void LDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V3>
    static void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V3>
    static void AliasLDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3);

    template <class V1>
    static void LDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2);
    template <class V1>
    static void NoAliasLDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2);
    template <class V1>
    static void AliasLDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2);

    template <int ix, class T, class V1, class V3>
    static void RDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V3>
    static void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V3>
    static void AliasRDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3);

    template <class V1>
    static void RDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2);
    template <class V1>
    static void NoAliasRDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2);
    template <class V1>
    static void AliasRDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2);

} // namespace tmv

#endif 
