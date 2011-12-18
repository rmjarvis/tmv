
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
    inline void MultMV(
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

    
    //
    // v * P
    //

    template <bool add, int ix, class T, class V1, class V3>
    TMV_INLINE void MultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    { MultMV<add>(x,m2.transpose(),v1,v3); }
 

    //
    // v *= P
    //

    template <int ix, class T, class V1>
    TMV_INLINE void MultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const Permutation& m2)
    { m2.applyOnRight(v1); Scale(x,v1); }


    //
    // v / P
    //

    template <int ix, class T, class V1, class V3>
    inline void LDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    {
        MultXV<false>(x,v1.vec(),v3.vec());
        m2.applyOnRight(v3);
    }


    //
    // v /= P
    //

    template <class V1>
    TMV_INLINE void LDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnRight(v1); }


    // 
    // v % P
    //

    template <int ix, class T, class V1, class V3>
    inline void RDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3)
    {
        MultXV<false>(x,v1.vec(),v3.vec());
        m2.applyOnLeft(v3);
    }

    // 
    // v %= P
    //

    template <class V1>
    TMV_INLINE void RDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2)
    { m2.applyOnLeft(v1); }

} // namespace tmv

#endif
