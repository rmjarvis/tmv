
//#define PRINTALGO_DIVVD

#include "tmv/TMV_DivVD.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_ScaleV.h"
#include "tmv/TMV_CopyV.h"

namespace tmv {

    //
    // Element-wise division
    //

    template <class T1, int C1, class T2, int C2, class T3>
    void InstElemDivVV(
        const T3 x, const ConstVectorView<T1,C1>& v1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    {
        typedef typename Traits<T3>::real_type RT;
        const Scaling<1,RT> one;
        if (v1.step() == 1 && v2.step() == 1 && v3.step() == 1) {
            ConstVectorView<T1,C1|Unit> v1unit = v1.unitView();
            ConstVectorView<T2,C2|Unit> v2unit = v2.unitView();
            VectorView<T3,Unit> v3unit = v3.unitView();
            InlineElemDivVV(one,v1unit,v2unit,v3unit);
        } else 
            InlineElemDivVV(one,v1,v2,v3);
        InstScale(x,v3);
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasElemDivVV(
        const T3 x, const ConstVectorView<T1,C1>& v1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    { InlineAliasElemDivVV(Scaling<0,T3>(x),v1,v2,v3); }

#define InstFile "TMV_DivVD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


