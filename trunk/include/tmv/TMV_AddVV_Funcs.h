

#ifndef TMV_AddVV_Funcs_H
#define TMV_AddVV_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //

    // From TMV_AddVV.h:
    template <int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    inline void AddVV(
        const Scaling<ix1,T1>& x1, const BaseVector_Calc<V1>& v1,
        const Scaling<ix2,T2>& x2, const BaseVector_Calc<V2>& v2,
        BaseVector_Mutable<V3>& v3);

} // namespace tmv

#endif 
