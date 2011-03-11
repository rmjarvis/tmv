
#ifndef TMV_MultPV_H
#define TMV_MultPV_H

#include "TMV_Permutation.h"
#include "TMV_MultXV.h"

namespace tmv {

    //
    // P * v
    // TODO: break out add option at compile time rather than runtime if
    //

    template <bool add, int ix, class T, class V2, class V3>
    static void MultMV(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        if (add) {
            v3 += (x*m1*v2).calc();
        } else {
            MultXV<false>(x,v2.vec(),v3.vec());
            m1.applyOnLeft(v3);
        }
    }
    template <bool add, int ix, class T, class V2, class V3>
    static void NoAliasMultMV(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        if (add) {
            v3 += (x*m1*v2).calc();
        } else {
            NoAliasMultXV<false>(x,v2.vec(),v3.vec());
            m1.applyOnLeft(v3);
        }
    }
    template <bool add, int ix, class T, class V2, class V3>
    static void AliasMultMV(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        if (add) {
            v3 += (x*m1*v2).calc();
        } else {
            AliasMultXV<false>(x,v2.vec(),v3.vec());
            m1.applyOnLeft(v3);
        }
    }

    
    //
    // v * P
    //

    template <bool add, int ix, class T, class V1, class V3>
    static void MultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    { MultMV<add>(x,m2.transpose(),v1,v3); }
    template <bool add, int ix, class T, class V1, class V3>
    static void NoAliasMultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    { NoAliasMultMV<add>(x,m2.transpose(),v1,v3); }
    template <bool add, int ix, class T, class V1, class V3>
    static void AliasMultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    { AliasMultMV<add>(x,m2.transpose(),v1,v3); }
 

    //
    // v *= P
    //

    template <int ix, class T, class V1>
    static void MultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const Permutation& m2)
    { m2.applyOnRight(v1); Scale(x,v1); }
    template <int ix, class T, class V1>
    static void NoAliasMultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const Permutation& m2)
    { m2.applyOnRight(v1); Scale(x,v1); }
    template <int ix, class T, class V1>
    static void AliasMultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const Permutation& m2)
    { m2.applyOnRight(v1); Scale(x,v1); }


    //
    // v / P
    //

    template <int ix, class T, class V1, class V3>
    static void LDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    { 
        MultXV<false>(x,v1.vec(),v3.vec());
        m2.applyOnRight(v3);
    }
    template <int ix, class T, class V1, class V3>
    static void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    { 
        NoAliasMultXV<false>(x,v1.vec(),v3.vec());
        m2.applyOnRight(v3);
    }
    template <int ix, class T, class V1, class V3>
    static void AliasLDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    { 
        AliasMultXV<false>(x,v1.vec(),v3.vec());
        m2.applyOnRight(v3);
    }


    //
    // v /= P
    //

    template <class V1>
    static void LDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnRight(v1); }
    template <class V1>
    static void NoAliasLDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnRight(v1); }
    template <class V1>
    static void AliasLDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnRight(v1); }


    // 
    // v % P
    //

    template <int ix, class T, class V1, class V3>
    static void RDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    { 
        MultXV<false>(x,v1.vec(),v3.vec());
        m2.applyOnLeft(v3);
    }
    template <int ix, class T, class V1, class V3>
    static void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    { 
        NoAliasMultXV<false>(x,v1.vec(),v3.vec());
        m2.applyOnLeft(v3);
    }
    template <int ix, class T, class V1, class V3>
    static void AliasRDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    { 
        AliasMultXV<false>(x,v1.vec(),v3.vec());
        m2.applyOnLeft(v3);
    }


    // 
    // v %= P
    //

    template <class V1>
    static void RDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnLeft(v1); }
    template <class V1>
    static void NoAliasRDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnLeft(v1); }
    template <class V1>
    static void AliasRDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnLeft(v1); }

} // namespace tmv

#endif
