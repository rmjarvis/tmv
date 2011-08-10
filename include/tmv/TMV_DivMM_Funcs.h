

#ifndef TMV_DivMM_Funcs_H
#define TMV_DivMM_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //


    // From TMV_DivM.h:
    template <class M1, class M2>
    static inline void LDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2);
    template <class M1, class M2>
    static inline void NoAliasLDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2);
    template <class M1, class M2>
    static inline void AliasLDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2);

    template <class M1, class M2>
    static inline void RDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2);
    template <class M1, class M2>
    static inline void NoAliasRDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2);
    template <class M1, class M2>
    static inline void AliasRDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2);

    template <int ix, class T, class M1, class M2, class M3>
    static inline void LDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3);
    template <int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3);
    template <int ix, class T, class M1, class M2, class M3>
    static inline void AliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3);

    template <int ix, class T, class M1, class M2, class M3>
    static inline void RDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3);
    template <int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3);
    template <int ix, class T, class M1, class M2, class M3>
    static inline void AliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3);

    // From TMV_DivMU.h:
    template <class M1, class M2>
    static inline void TriLDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2);
    template <class M1, class M2>
    static inline void NoAliasTriLDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2);
    template <class M1, class M2>
    static inline void AliasTriLDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2);

    // From TMV_DivUU.h:
    template <class M1, class M2>
    static inline void TriLDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2);
    template <class M1, class M2>
    static inline void NoAliasTriLDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2);
    template <class M1, class M2>
    static inline void AliasTriLDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2);

    // From TMV_MultPM.h
    template <int ix, class T, class M1, class M3>
    static inline void LDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <int ix, class T, class M1, class M2, class M3>
    static inline void AliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <int ix, class T, class M1, class M3>
    static inline void RDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <int ix, class T, class M1, class M2, class M3>
    static inline void AliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <class M1>
    static inline void LDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& m2);
    template <class M1>
    static inline void NoAliasLDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& m2);
    template <class M1>
    static inline void AliasLDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& m2);

    template <class M1>
    static inline void RDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& m2);
    template <class M1>
    static inline void NoAliasRDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& m2);
    template <class M1>
    static inline void AliasRDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& m2);

} // namespace tmv

#endif 
