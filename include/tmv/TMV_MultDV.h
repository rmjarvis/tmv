

#ifndef TMV_MultDV_H
#define TMV_MultDV_H

#include "TMV_BaseMatrix_Diag.h"
#include "TMV_BaseVector.h"
#include "TMV_MultVV_Funcs.h"

namespace tmv {

    // 
    // MultDV
    //

    template <bool add, int ix, class T, class M1, class V2, class V3>
    TMV_INLINE void MultMV(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    { ElemMultVV<add>(x,m1.diag(),v2.vec(),v3.vec()); }

    template <bool add, int ix, class T, class V1, class M2, class V3>
    TMV_INLINE void MultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    { ElemMultVV<add>(x,m2.diag(),v1.vec(),v3.vec()); }

    template <class V1, int ix, class T, class M2>
    TMV_INLINE void MultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { ElemMultVV<false>(x,m2.diag(),v1.vec(),v1.vec()); }


    //
    // MultDD
    //

    template <bool add, int ix, class T, class M1, class M2, class M3>
    TMV_INLINE void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Diag_Mutable<M3>& m3)
    {
        typename M3::diag_type m3d = m3.diag();
        ElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
    }

    template <class M1, int ix, class T, class M2>
    TMV_INLINE void MultEqMM(
        BaseMatrix_Diag_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    {
        typename M1::diag_type m1d = m1.diag();
        ElemMultVV<false>(x,m1.diag(),m2.diag(),m1d); 
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M3::diag_type m3d = m3.diag();
        if (!add && (SameStorage(m1,m3) || SameStorage(m2,m3))) {
            ElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
            if (m1.size() > 1) {
                m3.upperTri().offDiag().setZero();
                m3.lowerTri().offDiag().setZero();
            }
        } else {
            Maybe<!add>::zero(m3);
            ElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
        }
    }

} // namespace tmv

#endif
