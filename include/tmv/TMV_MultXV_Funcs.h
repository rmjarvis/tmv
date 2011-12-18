

#ifndef TMV_MultXV_Funcs_H
#define TMV_MultXV_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //

    // From TMV_ScaleV.h:
    template <int ix, class T, class V>
    inline void Scale(const Scaling<ix,T>& x, BaseVector_Mutable<V>& v);

    // From TMV_MultXV.h:
    template <bool add, int ix, class T, class V1, class V2>
    inline void MultXV(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        BaseVector_Mutable<V2>& v2);

} // namespace tmv

#endif

