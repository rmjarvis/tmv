
//#define PRINTALGO_Det

#include "tmv/TMV_Det.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_MinMax.h"

namespace tmv {

    //
    // ProdElements
    //

    template <class T>
    T InstProdElements(const ConstVectorView<T>& v)
    {
        if (v.step() == 1) {
            ConstVectorView<T,Unit> vunit = v.unitView();
            return InlineProdElements(vunit);
        } else {
            return InlineProdElements(v);
        }
    }

    template <class T>
    typename ConstVectorView<T>::float_type InstLogProdElements(
        const ConstVectorView<T>& v,
        typename ConstVectorView<T>::zfloat_type* sign)
    {
        if (v.step() == 1) {
            ConstVectorView<T,Unit> vunit = v.unitView();
            return InlineLogProdElements(vunit,sign);
        } else {
            return InlineLogProdElements(v,sign);
        }
    }

    template <class T>
    bool InstHasZeroElement(const ConstVectorView<T>& v)
    {
        if (v.step() == 1) {
            ConstVectorView<T,Unit> vunit = v.unitView();
            return InlineHasZeroElement(vunit);
        } else {
            return InlineHasZeroElement(v);
        }
    }

#define InstFile "TMV_Det.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


