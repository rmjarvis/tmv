
#include "TMV_Blas.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_MultMM_Winograd.h"
#include "tmv/TMV_MultMM_Block.h"
#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"

namespace tmv {

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM_Winograd(
        const T3 x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
#if !defined(BLAS) && TMV_OPT >= 3 && defined(TMV_MM_USE_WINOGRAD)
        InlineMultMM_Winograd<false>(Scaling<0,T3>(x),m1,m2,m3); 
#else
        InstMultMM(x,m1,m2,m3); 
#endif
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM_Winograd(
        const T3 x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
#if !defined(BLAS) && defined(TMV_MM_USE_WINOGRAD)
        InlineMultMM_Winograd<true>(Scaling<0,T3>(x),m1,m2,m3); 
#else
        InstAddMultMM(x,m1,m2,m3); 
#endif
    }

#define InstFile "TMV_MultMM_Winograd.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


