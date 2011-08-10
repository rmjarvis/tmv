

#ifndef TMV_DivVM_Funcs_H
#define TMV_DivVM_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //

    // From TMV_DivM.h:
    template <class V1, class M2>
    static inline void LDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2);
    template <class V1, class M2>
    static inline void NoAliasLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2);
    template <class V1, class M2>
    static inline void AliasLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2);

    template <class V1, class M2>
    static inline void RDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2);
    template <class V1, class M2>
    static inline void NoAliasRDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2);
    template <class V1, class M2>
    static inline void AliasRDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2);

    template <int ix, class T, class V1, class M2, class V3>
    static inline void LDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    static inline void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    static inline void AliasLDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);

    template <int ix, class T, class V1, class M2, class V3>
    static inline void RDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    static inline void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    static inline void AliasRDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);

    // From TMV_DivVU.h:
    template <class V1, class M2>
    static inline void TriLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Tri<M2>& m2);
    template <class V1, class M2>
    static inline void NoAliasTriLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Tri<M2>& m2);
    template <class V1, class M2>
    static inline void AliasTriLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Tri<M2>& m2);

    // From TMV_MultPV.h:
    template <int ix, class T, class V1, class V3>
    static inline void LDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V3>
    static inline void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V3>
    static inline void AliasLDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3);

    template <class V1>
    static inline void LDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2);
    template <class V1>
    static inline void NoAliasLDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2);
    template <class V1>
    static inline void AliasLDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2);

    template <int ix, class T, class V1, class V3>
    static inline void RDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V3>
    static inline void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V3>
    static inline void AliasRDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Permutation& m2, BaseVector_Mutable<V3>& v3);

    template <class V1>
    static inline void RDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2);
    template <class V1>
    static inline void NoAliasRDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2);
    template <class V1>
    static inline void AliasRDivEq(
        BaseVector_Mutable<V1>& v1, const Permutation& m2);

} // namespace tmv

#endif 
