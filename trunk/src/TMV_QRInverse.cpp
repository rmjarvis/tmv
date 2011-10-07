
//#define PRINTALGO_QR

#include "TMV_Blas.h"
#include "tmv/TMV_QRInverse.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_DivVU.h"
#include "tmv/TMV_DivMU.h"
#include "tmv/TMV_PermuteM.h"
#include "tmv/TMV_MultUL.h"
#include "tmv/TMV_MultPM.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_ScaleM.h"

#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_QuotMM.h"
#include "tmv/TMV_DivM.h"

namespace tmv {

    //
    // Inverse
    //
    
    template <class T1, int C1, class RT1, class T2>
    void InstQR_Inverse(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, MatrixView<T2> m2)
    {
        InlineQR_Inverse(QR,beta,P,N1,m2);
    }

    //
    // InverseATA
    //

    template <class T1, int C1, class RT1, class T2>
    void InstQR_InverseATA(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, MatrixView<T2> m2)
    {
        InlineQR_InverseATA(QR,beta,P,N1,m2);
    }

#define InstFile "TMV_QRInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


