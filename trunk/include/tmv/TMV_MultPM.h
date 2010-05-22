
#ifndef TMV_MultPM_H
#define TMV_MultPM_H

#include "TMV_Permutation.h"
#include "TMV_MultXM.h"

namespace tmv {

    //
    // P * m
    // TODO: break out add option at compile time rather than using a 
    //       runtime if statement.
    //

    template <bool add, int ix, class T, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const BaseMatrix<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += (x*m1*m2).calc();
        } else {
            MultXM<false>(x,m2.mat(),m3.mat());
            m1.applyOnLeft(m3);
        }
    }
    template <bool add, int ix, class T, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const BaseMatrix<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += (x*m1*m2).calc();
        } else {
            NoAliasMultXM<false>(x,m2.mat(),m3.mat());
            m1.applyOnLeft(m3);
        }
    }
    template <bool add, int ix, class T, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const BaseMatrix<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += (x*m1*m2).calc();
        } else {
            AliasMultXM<false>(x,m2.mat(),m3.mat());
            m1.applyOnLeft(m3);
        }
    }

    template <bool add, int ix, class T, class M1, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += (x*m1*m2).calc();
        } else {
            MultXM<false>(x,m1.mat(),m3.mat());
            m2.applyOnRight(m3);
        }
    }
    template <bool add, int ix, class T, class M1, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += (x*m1*m2).calc();
        } else {
            NoAliasMultXM<false>(x,m1.mat(),m3.mat());
            m2.applyOnRight(m3);
        }
    }
    template <bool add, int ix, class T, class M1, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += (x*m1*m2).calc();
        } else {
            AliasMultXM<false>(x,m1.mat(),m3.mat());
            m2.applyOnRight(m3);
        }
    }

    template <bool add, int ix, class T, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += (x*m1*m2).calc();
        } else {
            m3 = x*m1*m2.calc();
        }
    }
    template <bool add, int ix, class T, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    { MultMM<add>(x,m1,m2,m3); }
    template <bool add, int ix, class T, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    { MultMM<add>(x,m1,m2,m3); }

} // namespace tmv

#endif
