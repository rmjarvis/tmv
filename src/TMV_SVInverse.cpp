
//#define PRINTALGO_SVD

#include "TMV_Blas.h"
#include "tmv/TMV_SVInverse.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Permutation.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_DivM.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_MultMD.h"
#include "tmv/TMV_InvertD.h"
#include "tmv/TMV_InvertM.h"
#include "tmv/TMV_QuotXM.h"
#include "tmv/TMV_QuotMM.h"

namespace tmv {

    //
    // Inverse
    //
    
    template <class T1, int C1, class RT1, class T2>
    void InstSV_Inverse(
        const ConstMatrixView<T1,C1>& U, const ConstDiagMatrixView<RT1>& S,
        const ConstMatrixView<T1,C1>& V, MatrixView<T2> minv)
    { InlineSV_Inverse(U,S,V,minv); }

    //
    // InverseATA
    //

    template <class T1, int C1, class RT1, class T2>
    void InstSV_InverseATA(
        const ConstMatrixView<T1,C1>& U, const ConstDiagMatrixView<RT1>& S,
        const ConstMatrixView<T1,C1>& V, MatrixView<T2> ata)
    { InlineSV_InverseATA(U,S,V,ata); }



#define InstFile "TMV_SVInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


