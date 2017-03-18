

#ifndef TMV_AddMM_Funcs_H
#define TMV_AddMM_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //

    // From TMV_AddMM.h:
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);

    // From TMV_AddDD.h:
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Diag_Mutable<M3>& m3);

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Mutable<M3>& m3);

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Calc<M2>& m2,
        BaseMatrix_Mutable<M3>& m3);

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Calc<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
        BaseMatrix_Mutable<M3>& m3);

    // From TMV_AddUU.h:
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3);

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);


    // From TMV_AddBB.h:
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Band<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Band<M2>& m2, 
        BaseMatrix_Band_Mutable<M3>& m3);

} // namespace tmv

#endif 
