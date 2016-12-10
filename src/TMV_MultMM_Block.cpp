
#include "TMV_Blas.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_MultMM_Block.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"

namespace tmv {

    template <bool add, class M1, class M2, class T>
    static inline void DoMultMM_Block(
        const T& x, const M1& m1, const M2& m2, MatrixView<T> m3)
    { InlineMultMM_Block<add>(Scaling<0,T>(x),m1,m2,m3); }

    template <bool add, class M1, class M2, class T>
    static inline void DoMultMM_Block(
        const std::complex<T>& x, const M1& m1, const M2& m2,
        MatrixView<std::complex<T> > m3)
    {
        typedef typename Traits<T>::complex_type CT;
        InlineMultMM_Block<add>(Scaling<0,CT>(x),m1,m2,m3);
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM_Block(
        const T3 x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { 
#if !defined(BLAS) && TMV_OPT >= 3
        DoMultMM_Block<false>(x,m1,m2,m3); 
#else
        InstMultMM(x,m1,m2,m3); 
#endif
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM_Block(
        const T3 x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
#if !defined(BLAS) && TMV_OPT >= 1
        DoMultMM_Block<true>(x,m1,m2,m3); 
#else
        InstAddMultMM(x,m1,m2,m3); 
#endif
    }

    template <bool add, class M1, class M2, class T>
    static inline void DoMultMM_RecursiveBlock(
        const T& x, const M1& m1, const M2& m2, MatrixView<T> m3)
    { InlineMultMM_RecursiveBlock<add>(Scaling<0,T>(x),m1,m2,m3); }

    template <bool add, class M1, class M2, class T>
    static inline void DoMultMM_RecursiveBlock(
        const std::complex<T>& x, const M1& m1, const M2& m2,
        MatrixView<std::complex<T> > m3)
    {
        typedef typename Traits<T>::complex_type CT;
        InlineMultMM_RecursiveBlock<add>(Scaling<0,CT>(x),m1,m2,m3);
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM_RecursiveBlock(
        const T3 x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
#if !defined(BLAS) && (TMV_OPT >= 3) && defined(TMV_MM_USE_RECURSIVE_BLOCK)
        DoMultMM_RecursiveBlock<false>(x,m1,m2,m3); 
#else
        InstMultMM(x,m1,m2,m3); 
#endif
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM_RecursiveBlock(
        const T3 x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
#if !defined(BLAS) && (TMV_OPT >= 1) && defined(TMV_MM_USE_RECURSIVE_BLOCK)
        DoMultMM_RecursiveBlock<true>(x,m1,m2,m3); 
#else
        InstAddMultMM(x,m1,m2,m3); 
#endif
    }

#define InstFile "TMV_MultMM_Block.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmg


