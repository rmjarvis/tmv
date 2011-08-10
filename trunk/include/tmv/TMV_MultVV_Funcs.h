

#ifndef TMV_MultVV_Funcs_H
#define TMV_MultVV_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //


    // From TMV_MultVV.h:
    template <class V1, class V2>
    static inline typename ProdType<V1,V2>::type MultVV(
        const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2);

    // From TMV_ElemMultVV.h:
    template <bool add, int ix, class T, class V1, class V2, class V3>
    static inline void ElemMultVV(
        const Scaling<ix,T>& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V1, class V2, class V3>
    static inline void NoAliasElemMultVV(
        const Scaling<ix,T>& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V1, class V2, class V3>
    static inline void AliasElemMultVV(
        const Scaling<ix,T>& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);

} // namespace tmv

#endif 
