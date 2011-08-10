

#ifndef TMV_Rank1VVM_Funcs_H
#define TMV_Rank1VVM_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //

    // From TMV_Rank1_VVM.h:
    template <bool add, int ix, class T, class V1, class V2, class M3>
    static inline void Rank1Update(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class V1, class V2, class M3>
    static inline void NoAliasRank1Update(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class V1, class V2, class M3>
    static inline void AliasRank1Update(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseMatrix_Rec_Mutable<M3>& m3);

} // namespace tmv

#endif 
