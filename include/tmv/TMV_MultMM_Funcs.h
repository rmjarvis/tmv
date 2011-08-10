

#ifndef TMV_MultMM_Funcs_H
#define TMV_MultMM_Funcs_H

namespace tmv {

    //
    // This is just a list of all declarations for all the functions
    // that are defined elsewhere to help get the overloading right
    // so you don't have to worry so much about which order header 
    // files are included.
    //


    // From TMV_MultMM.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <class M1, int ix, class T, class M2>
    static inline void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static inline void NoAliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static inline void AliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M2>& m2);

    // From TMV_MultMM_Block.h
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM_RecursiveBlock(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM_Block(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    // From TMV_MultMM_Winograd.h
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM_Winograd(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    // From TMV_MultMM_OpenMP
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM_OpenMP(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    // From TMV_MultMD.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <class M1, int ix, class T, class M2>
    static inline void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static inline void NoAliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static inline void AliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);

    // From TMV_MultDV.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Diag_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Diag_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Diag_Mutable<M3>& m3);

    template <class M1, int ix, class T, class M2>
    static inline void MultEqMM(
        BaseMatrix_Diag_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static inline void NoAliasMultEqMM(
        BaseMatrix_Diag_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static inline void AliasMultEqMM(
        BaseMatrix_Diag_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    // From TMV_MultUM.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <class M1, int ix, class T, class M2>
    static inline void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static inline void NoAliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static inline void AliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);

    // From TMV_MultUU.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3);

    template <class M1, int ix, class T, class M2>
    static inline void MultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static inline void NoAliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static inline void AliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);

    // From TMV_MultUL.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    // From TMV_MultUD.h:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3);

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Mutable<M3>& m3);

    template <class M1, int ix, class T, class M2>
    static inline void MultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static inline void NoAliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    static inline void AliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);


    // From TMV_MultPM.h
    template <bool add, int ix, class T, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const BaseMatrix<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const BaseMatrix<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const BaseMatrix<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <bool add, int ix, class T, class M1, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <bool add, int ix, class T, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const Permutation& m1, 
        const Permutation& m2, BaseMatrix_Rec_Mutable<M3>& m3);

    template <class M1, int ix, class T>
    static inline void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const Permutation& m2);
    template <class M1, int ix, class T>
    static inline void NoAliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const Permutation& m2);
    template <class M1, int ix, class T>
    static inline void AliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const Permutation& m2);

} // namespace tmv

#endif
