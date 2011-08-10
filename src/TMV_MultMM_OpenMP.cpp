
#include "TMV_Blas.h"
#include "tmv/TMV_MultMM_OpenMP.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_Matrix.h"

namespace tmv {

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM_OpenMP(
        const T3 x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
#if !defined(BLAS) && TMV_OPT > 2 && defined(_OPENMP)
        InlineMultMM_OpenMP<false>(Scaling<0,T3>(x),m1,m2,m3); 
#else
        InstMultMM(x,m1,m2,m3);
#endif
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM_OpenMP(
        const T3 x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
#if !defined(BLAS) && TMV_OPT > 0 && defined(_OPENMP)
        InlineMultMM_OpenMP<true>(Scaling<0,T3>(x),m1,m2,m3); 
#else
        InstAddMultMM(x,m1,m2,m3);
#endif
    }

#define InstFile "TMV_MultMM_OpenMP.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


