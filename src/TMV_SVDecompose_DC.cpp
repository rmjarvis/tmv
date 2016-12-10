
//#undef NDEBUG
//#define PRINTALGO_SVD
//#define XDEBUG_SVD
//#include "TMV.h"

#include "tmv/TMV_SVDecompose_DC.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_ProdMV.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_SVDecompose_QR.h"
#include "tmv/TMV_MinMax.h"
#include "tmv/TMV_SwapV.h"
#include "tmv/TMV_SortV.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_NormV.h"
#include "tmv/TMV_MultPV.h"
#include "tmv/TMV_MultPM.h"
#include "tmv/TMV_SumVV.h"
#include "tmv/TMV_SumMM.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_PermuteM.h"

namespace tmv {

    template <class Tu, class RT, class Tv>
    void InstSV_DecomposeFromBidiagonal_DC(
        MatrixView<Tu> U, VectorView<RT> D, VectorView<RT> E,
        MatrixView<Tv> V, bool UisI, bool VisI)
    {
        TMVAssert(U.iscm());
        TMVAssert(V.iscm() || V.isrm());
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        if (D.size() > 0) {
            MatrixView<Tu,ColMajor> Ucm = U;
            VectorView<RT,Unit> Du = D;
            VectorView<RT,Unit> Eu = E;
            if (V.iscm()) {
                MatrixView<Tv,ColMajor> Vcm = V;
                InlineSV_DecomposeFromBidiagonal_DC(Ucm,Du,Eu,Vcm,UisI,VisI);
            } else {
                MatrixView<Tv,RowMajor> Vrm = V;
                InlineSV_DecomposeFromBidiagonal_DC(Ucm,Du,Eu,Vrm,UisI,VisI);
            }
        }
    }


#define InstFile "TMV_SVDecompose_DC.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


