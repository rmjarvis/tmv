

#ifndef TMV_ProdVV_H
#define TMV_ProdVV_H

#include "TMV_BaseVector.h"
#include "TMV_ProdXV.h"
#include "TMV_MultVV_Funcs.h"

namespace tmv {

    //
    // Vector * Vector
    //

#define PT typename ProdType<V1,V2>::type

    // These first few are for when an argument is a composite vector
    // and needs to be calculated before running MultVV.
    template <class V1, class V2>
    static TMV_INLINE PT MultVV(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    { return MultVV(v1.calc(),v2.calc()); }
    template <class V1, class V2>
    static TMV_INLINE PT NoAliasMultVV(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    { return NoAliasMultVV(v1.calc(),v2.calc()); }
    template <class V1, class V2>
    static TMV_INLINE PT AliasMultVV(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    { return AliasMultVV(v1.calc(),v2.calc()); }

    // v * v
    template <class V1, class V2>
    static TMV_INLINE_ND PT operator*(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2) 
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVAssert(v1.size() == v2.size());
        return MultVV(v1.vec(),v2.vec()); 
    }

#define PT2 typename Traits2<Tx,PT>::type
    // v * (x*v)
    template <class V1, int ix2, class Tx, class V2>
    static TMV_INLINE_ND PT2 operator*(
        const BaseVector<V1>& v1, const ProdXV<ix2,Tx,V2>& v2) 
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVAssert(v1.size() == v2.size());
        return v2.getX() * (v1.vec() * v2.getV().vec());
    }

    // (x*v) * v
    template <int ix1, class Tx, class V1, class V2>
    static TMV_INLINE_ND PT2 operator*(
        const ProdXV<ix1,Tx,V1>& v1, const BaseVector<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVAssert(v1.size() == v2.size());
        return v1.getX() * (v1.getV().vec() * v2.vec());
    }
#undef PT2

#define PT2 typename Traits2<Tx1,typename Traits2<Tx2,PT>::type>::type
    // (x*v) * (x*v)
    template <int ix1, class Tx1, class V1, int ix2, class Tx2, class V2>
    static TMV_INLINE_ND PT2 operator*(
        const ProdXV<ix1,Tx1,V1>& v1, const ProdXV<ix2,Tx2,V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVAssert(v1.size() == v2.size());
        return v1.getX() * v2.getX() * (v1.getV().vec() * v2.getV().vec());
    }
#undef PT2

#undef PT

} // namespace tmv

#endif 
