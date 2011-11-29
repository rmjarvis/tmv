

#ifndef TMV_MultXM_Funcs_H
#define TMV_MultXM_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //

    // From TMV_MultXM.h:
    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        BaseMatrix_Rec_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        BaseMatrix_Rec_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        BaseMatrix_Rec_Mutable<M2>& m2);

    // From TMV_MultXD.h:
    template <bool add, int ix1, class T1, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Diag_Mutable<M2>& m2);
    template <bool add, int ix1, class T1, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Diag_Mutable<M2>& m2);
    template <bool add, int ix1, class T1, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Diag_Mutable<M2>& m2);

    template <bool add, int ix1, class T1, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2);
    template <bool add, int ix1, class T1, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2);
    template <bool add, int ix1, class T1, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2);

    // From TMV_MultXU.h:
    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1, 
        BaseMatrix_Tri_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1, 
        BaseMatrix_Tri_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1, 
        BaseMatrix_Tri_Mutable<M2>& m2);

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2);

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Tri_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Tri_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Tri_Mutable<M2>& m2);

    // From TMV_MultXB.h:
    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1, 
        BaseMatrix_Band_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1, 
        BaseMatrix_Band_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1, 
        BaseMatrix_Band_Mutable<M2>& m2);

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2);

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        BaseMatrix_Tri_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        BaseMatrix_Tri_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        BaseMatrix_Tri_Mutable<M2>& m2);

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Band_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Band_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Band_Mutable<M2>& m2);

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        BaseMatrix_Band_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        BaseMatrix_Band_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        BaseMatrix_Band_Mutable<M2>& m2);


} // namespace tmv

#endif
