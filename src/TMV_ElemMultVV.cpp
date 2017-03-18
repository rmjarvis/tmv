

#include "tmv/TMV_ElemMultVV.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_MultXV.h"

namespace tmv {

    //
    // ElementProd
    //

    template <bool add, class T, class V1, class V2, class V3> 
    static void DoElemMultVV(
        const T x, const V1& v1, const V2& v2, V3& v3)
    {
        if (x == T(1))
            InlineElemMultVV<add>(Scaling<1,T>(x),v1,v2,v3);
        else if (x == T(-1))
            InlineElemMultVV<add>(Scaling<-1,T>(x),v1,v2,v3);
        else if (x == T(0)) 
            Maybe<!add>::zero(v3);
        else
            InlineElemMultVV<add>(Scaling<0,T>(x),v1,v2,v3);
    }

    template <bool add, class T, class V1, class V2, class V3> 
    static void DoElemMultVV(
        const std::complex<T> x, const V1& v1, const V2& v2, V3& v3)
    {
        if (imag(x) == T(0)) {
            if (real(x) == T(1))
                InlineElemMultVV<add>(Scaling<1,T>(real(x)),v1,v2,v3);
            else if (real(x) == T(-1))
                InlineElemMultVV<add>(Scaling<-1,T>(real(x)),v1,v2,v3);
            else if (real(x) == T(0)) 
                Maybe<!add>::zero(v3);
            else
                InlineElemMultVV<add>(Scaling<0,T>(real(x)),v1,v2,v3);
        } else
            InlineElemMultVV<add>(Scaling<0,std::complex<T> >(x),v1,v2,v3);
    }

    template <class T1, int C1, class T2, int C2, class T3> 
    void InstElemMultVV(
        const T3 x, const ConstVectorView<T1,C1>& v1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    {
        if (v1.step() == 1 && v2.step() == 1 && v3.step() == 1) {
            ConstVectorView<T1,C1|Unit> v1unit = v1.unitView();
            ConstVectorView<T2,C2|Unit> v2unit = v2.unitView();
            VectorView<T3,Unit> v3unit = v3.unitView();
            DoElemMultVV<false>(x,v1unit,v2unit,v3unit);
        } else 
            DoElemMultVV<false>(x,v1,v2,v3);
    }

    template <class T1, int C1, class T2, int C2, class T3> 
    void InstAddElemMultVV(
        const T3 x, const ConstVectorView<T1,C1>& v1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    {
        if (v1.step() == 1 && v2.step() == 1 && v3.step() == 1) {
            ConstVectorView<T1,C1|Unit> v1unit = v1.unitView();
            ConstVectorView<T2,C2|Unit> v2unit = v2.unitView();
            VectorView<T3,Unit> v3unit = v3.unitView();
            DoElemMultVV<true>(x,v1unit,v2unit,v3unit);
        } else 
            DoElemMultVV<true>(x,v1,v2,v3);
    }

    template <class T1, int C1, class T2, int C2, class T3> 
    void InstAliasElemMultVV(
        const T3 x, const ConstVectorView<T1,C1>& v1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    { InlineAliasElemMultVV<false>(Scaling<0,T3>(x),v1,v2,v3); }

    template <class T1, int C1, class T2, int C2, class T3> 
    void InstAliasAddElemMultVV(
        const T3 x, const ConstVectorView<T1,C1>& v1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    { InlineAliasElemMultVV<true>(Scaling<0,T3>(x),v1,v2,v3); }

#define InstFile "TMV_ElemMultVV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


