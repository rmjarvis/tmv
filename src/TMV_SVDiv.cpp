
//#define PRINTALGO_SVD

#include "TMV_Blas.h"
#include "tmv/TMV_SVDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_SmallDiagMatrix.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_TransposeM.h"
#include "tmv/TMV_ConjugateV.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_QuotXM.h"
#include "tmv/TMV_DivM.h"
#include "tmv/TMV_DivVD.h"
#include "tmv/TMV_InvertM.h"
#include "tmv/TMV_InvertD.h"
#include "tmv/TMV_QuotXM.h"
#include "tmv/TMV_MultMD.h"
#include "tmv/TMV_Det.h"


namespace tmv {

    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstSV_Solve(
        const ConstMatrixView<T1,C1>& U, const ConstDiagMatrixView<RT1>& S,
        const ConstMatrixView<T1,C1>& V, 
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineSV_Solve(U,S,V,m2,m3); }
    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstSV_Solve(
        const ConstMatrixView<T1,C1>& U, const ConstDiagMatrixView<RT1>& S,
        const ConstMatrixView<T1,C1>& V, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    { InlineSV_Solve(U,S,V,v2,v3); }

#define InstFile "TMV_SVDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


