
#include "TMV_Blas.h"
#include "tmv/TMV_AddBB.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_ScaleB.h"

namespace tmv {

    // Just do the most common cases: m+m, m-m
    // For anything else, do it in two function calls to MultXM.
    template <class M1, class M2, class T3, int A3>
    static void DoAddMM(
        const T3 x1, const M1& m1, const T3 x2, const M2& m2,
        BandMatrixView<T3,A3> m3)
    {
        typedef typename Traits<T3>::real_type RT;
        if (x1 == RT(1)) {
            if (x2 == RT(1)) {
                InlineAddMM(Scaling<1,RT>(),m1,Scaling<1,RT>(),m2,m3);
            } else if (x2 == RT(-1)) {
                InlineAddMM(Scaling<1,RT>(),m1,Scaling<-1,RT>(),m2,m3);
            } else if (x2 == RT(0)) {
                InstMultXM(x1,m1.xView(),m3.xView());
            } else {
                InstMultXM(x1,m1.xView(),m3.xView());
                InstAddMultXM(x2,m2.xView(),m3.xView());
            }
        } else if (x1 == RT(0)) {
            InstMultXM(x2,m2.xView(),m3.xView());
        } else {
            InstMultXM(x1,m1.xView(),m3.xView());
            InstAddMultXM(x2,m2.xView(),m3.xView());
        }
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMM(
        const T3 x1, const ConstBandMatrixView<T1,C1>& m1,
        const T3 x2, const ConstBandMatrixView<T2,C2>& m2,
        BandMatrixView<T3> m3)
    {
#ifdef BLAS
        InstMultXM(x1,m1,m3);
        InstAddMultXM(x2,m2,m3);
#else
        if (m1.iscm() && m2.iscm() && m3.iscm()) 
            DoAddMM(x1,m1.cmView(),x2,m2.cmView(),m3.cmView());
        else if (m1.isrm() && m2.isrm() && m3.isrm()) 
            DoAddMM(x1,m1.rmView(),x2,m2.rmView(),m3.rmView());
        else if (m1.isdm() && m2.isdm() && m3.isdm()) 
            DoAddMM(x1,m1.dmView(),x2,m2.dmView(),m3.dmView());
        else if (x1 == T3(0)) 
            InstMultXM(x2,m2,m3);
        else if (x2 == T3(0)) 
            InstMultXM(x1,m1,m3);
        else {
            InstMultXM(x1,m1,m3);
            InstAddMultXM(x2,m2,m3);
        }
#endif
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMM(
        const T3 x1, const ConstBandMatrixView<T1,C1>& m1,
        const T3 x2, const ConstBandMatrixView<T2,C2>& m2,
        BandMatrixView<T3> m3)
    { InlineAliasAddMM(Scaling<0,T3>(x1),m1,Scaling<0,T3>(x2),m2,m3); }

#define InstFile "TMV_AddBB.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


