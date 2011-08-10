
#include "TMV_Blas.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_MultMM_OpenMP.h"
#include "tmv/TMV_MultMM_Block.h"
#include "tmv/TMV_MultMM_Winograd.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_SumMM.h"
#include "tmv/TMV_Vector.h"

#ifdef BLAS
#include "TMV_MultMM_Blas.h"
#endif

namespace tmv {

    template <bool add, class M1, class M2, class T3>
    static void DoMultMM(
        const T3 x, const M1& m1, const M2& m2,
        MatrixView<T3,ColMajor> m3)
    {
        typedef typename Traits<T3>::real_type RT;
        if (TMV_IMAG(x) == RT(0)) {
            if (x == RT(1))
                InlineMultMM<add>(Scaling<1,RT>(),m1,m2,m3);
            else if (x == RT(-1))
                InlineMultMM<add>(Scaling<-1,RT>(),m1,m2,m3);
            else if (x == RT(0))
                Maybe<!add>::zero(m3); 
            else
                InlineMultMM<add>(Scaling<0,RT>(TMV_REAL(x)),m1,m2,m3);
        } else {
            InlineMultMM<add>(Scaling<0,T3>(x),m1,m2,m3);
        }
    }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
    template <bool add, class M1, class M2>
    static void DoMultMM(
        double x, const M1& m1, const M2& m2, MatrixView<double,ColMajor> m3)
    { BlasMultMM(x,m1,m2,add?1:0,m3); }
    template <bool add, class M1, class M2>
    static void DoMultMM(
        std::complex<double> x, const M1& m1, const M2& m2,
        MatrixView<std::complex<double>,ColMajor> m3)
    { BlasMultMM(x,m1,m2,add?1:0,m3); }
#endif // TMV_INST_DOUBLE
#ifdef TMV_INST_FLOAT
    template <bool add, class M1, class M2>
    static void DoMultMM(
        float x, const M1& m1, const M2& m2, MatrixView<float,ColMajor> m3)
    { BlasMultMM(x,m1,m2,add?1:0,m3); }
    template <bool add, class M1, class M2>
    static void DoMultMM(
        std::complex<float> x, const M1& m1, const M2& m2,
        MatrixView<std::complex<float>,ColMajor> m3)
    { BlasMultMM(x,m1,m2,add?1:0,m3); }
#endif // TMV_INST_FLOAT
#endif // BLAS

    template <bool add, class T1, int C1, class T2, int C2, class T3>
    void DoInstMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3,ColMajor> m3)
    {
        if (m3.colsize() > 0 && m3.rowsize() > 0) {
            if (m1.rowsize() == 0) Maybe<!add>::zero(m3);
            else DoMultMM<add>(x,m1,m2,m3); 
        }
    }
 

} // namespace tmv


