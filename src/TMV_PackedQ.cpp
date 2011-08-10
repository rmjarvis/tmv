
//#define PRINTALGO_QR

#include "TMV_Blas.h"
#include "tmv/TMV_PackedQ.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_ConjugateV.h"

namespace tmv {

    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_MultEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        MatrixView<T2> m2)
    { InlinePackedQ_MultEq(Q,beta,m2); }
    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_LDivEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        MatrixView<T2> m2)
    { InlinePackedQ_LDivEq(Q,beta,m2); }
    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_MultEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        VectorView<T2> v2)
    { InlinePackedQ_MultEq(Q,beta,v2); }
    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_LDivEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        VectorView<T2> v2)
    { InlinePackedQ_LDivEq(Q,beta,v2); }


#define InstFile "TMV_PackedQ.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


