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


#ifndef TMV_MultMV_Funcs_H
#define TMV_MultMV_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //


    // From TMV_MultMV.h:
    template <bool add, int ix, class T, class M1, class V2, class V3>
    static inline void MultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class M1, class V2, class V3>
    static inline void NoAliasMultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class M1, class V2, class V3>
    static inline void AliasMultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3);

    template <bool add, int ix, class T, class V1, class M2, class V3>
    static inline void MultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Rec<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V1, class M2, class V3>
    static inline void NoAliasMultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Rec<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V1, class M2, class V3>
    static inline void AliasMultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Rec<M2>& m2, BaseVector_Mutable<V3>& v3);

    template <class V1, int ix, class T, class M2>
    static inline void MultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M2>& m2);
    template <class V1, int ix, class T, class M2>
    static inline void NoAliasMultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M2>& m2);
    template <class V1, int ix, class T, class M2>
    static inline void AliasMultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M2>& m2);

    // From TMV_MultDV.h:
    template <bool add, int ix, class T, class M1, class V2, class V3>
    static inline void MultMV(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class M1, class V2, class V3>
    static inline void NoAliasMultMV(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class M1, class V2, class V3>
    static inline void AliasMultMV(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);

    template <bool add, int ix, class T, class V1, class M2, class V3>
    static inline void MultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V1, class M2, class V3>
    static inline void NoAliasMultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V1, class M2, class V3>
    static inline void AliasMultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3);

    template <class V1, int ix, class T, class M2>
    static inline void MultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class V1, int ix, class T, class M2>
    static inline void NoAliasMultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class V1, int ix, class T, class M2>
    static inline void AliasMultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);

    // From TMV_MultUV.h:
    template <bool add, int ix, class T, class M1, class V2, class V3>
    static inline void MultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class M1, class V2, class V3>
    static inline void NoAliasMultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class M1, class V2, class V3>
    static inline void AliasMultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3);

    template <bool add, int ix, class T, class V1, class M2, class V3>
    static inline void MultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V1, class M2, class V3>
    static inline void NoAliasMultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V1, class M2, class V3>
    static inline void AliasMultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3);
    
    template <class V1, int ix, class T, class M2>
    static inline void MultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class V1, int ix, class T, class M2>
    static inline void NoAliasMultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class V1, int ix, class T, class M2>
    static inline void AliasMultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);

    // From TMV_MultPV.h:
    template <bool add, int ix, class T, class V2, class V3>
    static inline void MultMV(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V2, class V3>
    static inline void NoAliasMultMV(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V2, class V3>
    static inline void AliasMultMV(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3);

    template <bool add, int ix, class T, class V1, class V3>
    static inline void MultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const Permutation& m2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V1, class V3>
    static inline void NoAliasMultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const Permutation& m2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V1, class V3>
    static inline void AliasMultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const Permutation& m2, BaseVector_Mutable<V3>& v3);
    
    template <int ix, class T, class V1>
    static inline void MultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const Permutation& m2);
    template <int ix, class T, class V1>
    static inline void NoAliasMultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const Permutation& m2);
    template <int ix, class T, class V1>
    static inline void AliasMultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const Permutation& m2);

} // namespace tmv

#endif

