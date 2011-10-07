
//#define XDEBUG_MM

#include "TMV_Blas.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_MultMM_OpenMP.h"
#include "tmv/TMV_MultMM_Block.h"
#include "tmv/TMV_MultMM_Winograd.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"

namespace tmv {

    template <class T1, int C1, class T2, int C2, class T3>
    void DoInstMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3,ColMajor> m3)
    {
#if TMV_OPT > 2 && !defined(BLAS)
        typedef typename Traits<T3>::real_type RT;
        if (TMV_IMAG(x) == RT(0)) {
            if (x == RT(1))
                InlineMultMM<false>(Scaling<1,RT>(),m1,m2,m3);
            else if (x == RT(-1))
                InlineMultMM<false>(Scaling<-1,RT>(),m1,m2,m3);
            else if (x == RT(0))
                m3.setZero();
            else
                InlineMultMM<false>(Scaling<0,RT>(TMV_REAL(x)),m1,m2,m3);
        } else {
            InlineMultMM<false>(Scaling<0,T3>(x),m1,m2,m3);
        }
#else
        InstMultMM(x,m1,m2,m3.xView());
#endif
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void DoInstAddMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3,ColMajor> m3)
    {
#if !defined(BLAS)
        typedef typename Traits<T3>::real_type RT;
        if (TMV_IMAG(x) == RT(0)) {
            if (x == RT(1))
                InlineMultMM<true>(Scaling<1,RT>(),m1,m2,m3);
            else if (x == RT(-1))
                InlineMultMM<true>(Scaling<-1,RT>(),m1,m2,m3);
            else if (x != RT(0))
                InlineMultMM<true>(Scaling<0,RT>(TMV_REAL(x)),m1,m2,m3);
        } else {
            InlineMultMM<true>(Scaling<0,T3>(x),m1,m2,m3);
        }
#else
        InstAddMultMM(x,m1,m2,m3.xView());
#endif
    }

#define InstFile "TMV_MultMM_CRC.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


