

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
    inline void LDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2);

    template <class M1, class M2>
    inline void RDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2);

    template <int ix, class T, class M1, class M2, class M3>
    inline void LDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3);

    template <int ix, class T, class M1, class M2, class M3>
    inline void RDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3);

    // From TMV_DivMU.h:
    template <class M1, class M2>
    inline void TriLDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2);

    // From TMV_DivUU.h:
    template <class M1, class M2>
    inline void TriLDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2);

    // From TMV_MultPM.h
    template <int ix, class T, class M1, class M3>
    inline void LDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <int ix, class T, class M1, class M3>
    inline void RDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <class M1>
    inline void LDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& m2);

    template <class M1>
    inline void RDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& m2);

} // namespace tmv

#endif 
