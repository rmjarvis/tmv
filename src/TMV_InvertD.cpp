
#include "tmv/TMV_InvertD.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_ScaleV.h"
#include "tmv/TMV_MultXD.h"
#include "tmv/TMV_Det.h"

namespace tmv {

    //
    // Element-wise division
    //

    template <class T>
    void InstElemInvert(VectorView<T> v)
    {
        if (v.step() == 1) {
            VectorView<T,Unit> vunit = v.unitView();
            InlineElemInvert(vunit);
        } else {
            InlineElemInvert(v);
        }
    }

#define InstFile "TMV_InvertD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


