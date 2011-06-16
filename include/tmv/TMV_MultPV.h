
#ifndef TMV_MultPV_H
#define TMV_MultPV_H

#include "TMV_Permutation.h"
#include "TMV_BaseVector.h"
#include "TMV_MultXV_Funcs.h"

namespace tmv {

    //
    // P * v
    // TODO: break out add option at compile time rather than runtime if
    //

    template <bool add, int ix, class T, class V2, class V3>
    static inline void MultMV(
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
    static inline void NoAliasMultMV(
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
    static inline void AliasMultMV(
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
    static TMV_INLINE void MultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    { MultMV<add>(x,m2.transpose(),v1,v3); }
    template <bool add, int ix, class T, class V1, class V3>
    static TMV_INLINE void NoAliasMultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    { NoAliasMultMV<add>(x,m2.transpose(),v1,v3); }
    template <bool add, int ix, class T, class V1, class V3>
    static TMV_INLINE void AliasMultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    { AliasMultMV<add>(x,m2.transpose(),v1,v3); }
 

    //
    // v *= P
    //

    template <int ix, class T, class V1>
    static TMV_INLINE void MultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const Permutation& m2)
    { m2.applyOnRight(v1); Scale(x,v1); }
    template <int ix, class T, class V1>
    static TMV_INLINE void NoAliasMultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const Permutation& m2)
    { m2.applyOnRight(v1); Scale(x,v1); }
    template <int ix, class T, class V1>
    static TMV_INLINE void AliasMultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const Permutation& m2)
    { m2.applyOnRight(v1); Scale(x,v1); }


    //
    // v / P
    //

    template <int ix, class T, class V1, class V3>
    static inline void LDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    {
        MultXV<false>(x,v1.vec(),v3.vec());
        m2.applyOnRight(v3);
    }
    template <int ix, class T, class V1, class V3>
    static inline void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    {
        NoAliasMultXV<false>(x,v1.vec(),v3.vec());
        m2.applyOnRight(v3);
    }
    template <int ix, class T, class V1, class V3>
    static inline void AliasLDiv(
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
    static TMV_INLINE void LDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnRight(v1); }
    template <class V1>
    static TMV_INLINE void NoAliasLDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnRight(v1); }
    template <class V1>
    static TMV_INLINE void AliasLDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnRight(v1); }


    // 
    // v % P
    //

    template <int ix, class T, class V1, class V3>
    static inline void RDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    {
        MultXV<false>(x,v1.vec(),v3.vec());
        m2.applyOnLeft(v3);
    }
    template <int ix, class T, class V1, class V3>
    static inline void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    {
        NoAliasMultXV<false>(x,v1.vec(),v3.vec());
        m2.applyOnLeft(v3);
    }
    template <int ix, class T, class V1, class V3>
    static inline void AliasRDiv(
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
    static TMV_INLINE void RDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnLeft(v1); }
    template <class V1>
    static TMV_INLINE void NoAliasRDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnLeft(v1); }
    template <class V1>
    static TMV_INLINE void AliasRDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnLeft(v1); }

} // namespace tmv

#endif
