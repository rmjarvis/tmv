
//#define PRINTALGO_Det

#include "tmv/TMV_IntegerDet.h"

namespace tmv {

    template <class T>
    T InstIntegerDet(const ConstMatrixView<T>& m)
    { return InlineIntegerDet(m); }

#define InstFile "TMV_IntegerDet.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


