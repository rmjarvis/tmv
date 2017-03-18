

#ifndef TMV_MultXD_H
#define TMV_MultXD_H

#include "TMV_BaseMatrix_Diag.h"
#include "TMV_MultXV_Funcs.h"

namespace tmv {

    //
    // D *= x
    //

    template <int ix, class T, class M>
    inline void Scale(
        const Scaling<ix,T>& x, BaseMatrix_Diag_Mutable<M>& m)
    {
        typename M::diag_type md = m.diag();
        Scale(x,md);
    }


    //
    // D (+)= x * D
    //

    template <bool add, int ix1, class T1, class M1, class M2>
    inline void MultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Diag_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        MultXV<add>(x1,m1d,m2d);
    }


    //
    // M (+)= x * D
    //

    template <bool add, int ix1, class T1, class M1, class M2>
    inline void MultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        if (SameStorage(m1,m2)) {
            MultXV<add>(x1,m1d,m2d);
            Maybe<!add>::zero_offdiag2(m2.upperTri());
            Maybe<!add>::zero_offdiag2(m2.lowerTri());
        } else {
            Maybe<!add>::zero(m2);
            typename M2::diag_type::noalias_type m2dna = m2d.noAlias();
            MultXV<add>(x1,m1d,m2dna);
        }
    }

} // namespace tmv

#endif
