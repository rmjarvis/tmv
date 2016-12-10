

#ifndef TMV_MultVV_Funcs_H
#define TMV_MultVV_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //

    // From TMV_InvertM.h:
    template <int ix, class T, class M1, class M2>
    inline void MakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& m2);

    template <class M1, class M2>
    inline void MakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2);

    // From TMV_InvertD.h:
    template <class M>
    inline void InvertSelf(BaseMatrix_Diag_Mutable<M>& m);

    // From TMV_InvertU.h:
    template <class M>
    inline void InvertSelf(BaseMatrix_Tri_Mutable<M>& m);

} // namespace tmv

#endif 
