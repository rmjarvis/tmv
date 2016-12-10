
#ifndef TMV_MultPM_H
#define TMV_MultPM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_Permutation.h"
#include "TMV_MultXM_Funcs.h"

namespace tmv {

    //
    // P * m
    // TODO: break out add option at compile time rather than using a 
    //       runtime if statement.
    //

    template <int ix, class T, class M1, class M2> class ProdMM;

    template <bool add, int ix, class T, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const BaseMatrix<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += ProdMM<ix,T,Permutation,M2>(x,m1,m2).calc();
        } else {
            MultXM<false>(x,m2.calc(),m3.mat());
            m1.applyOnLeft(m3);
        }
    }

    //
    // m * P
    //
    
    template <bool add, int ix, class T, class M1, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += ProdMM<ix,T,M1,Permutation>(x,m1,m2).calc();
        } else {
            MultXM<false>(x,m1.calc(),m3.mat());
            m2.applyOnRight(m3);
        }
    }

    //
    // P * P
    //

    template <bool add, int ix, class T, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += ProdMM<ix,T,Permutation,Permutation>(x,m1,m2).calc();
        } else {
            m3 = ProdMM<ix,T,Permutation,Permutation::calc_type>(x,m1,m2.calc());
        }
    }


    //
    // m / P
    // 

    template <int ix, class T, class M1, class M3>
    TMV_INLINE void LDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    { MultMM<false>(x,m2.inverse(),m1.calc(),m3.mat()); }


    //
    // m % P
    // 

    template <int ix, class T, class M1, class M3>
    TMV_INLINE void RDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    { MultMM<false>(x,m1.calc(),m2.inverse(),m3.mat()); }


    // 
    // m *= P
    //

    template <class M1, int ix, class T>
    inline void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const Permutation& m2)
    {
        m2.applyOnRight(m1); 
        Scale(x,m1.mat());
    }

    //
    // m /= P
    //

    template <class M1>
    TMV_INLINE void LDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& m2)
    { m2.inverse().applyOnLeft(m1); }

    //
    // m %= P
    //

    template <class M1>
    TMV_INLINE void RDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& m2)
    { m2.inverse().applyOnRight(m1); }



} // namespace tmv

#endif
