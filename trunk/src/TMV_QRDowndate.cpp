
//#define PRINTALGO_QR

#include "tmv/TMV_QRDowndate.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SmallTriMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_SwapV.h"
#include "tmv/TMV_MinMax.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_Rank1VVM.h"
#include "tmv/TMV_MultMV.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_MultMM_Block.h"
#include "tmv/TMV_MultMM_Winograd.h"
#include "tmv/TMV_MultMM_OpenMP.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultUM.h"
#include "tmv/TMV_DivM.h"
#include "tmv/TMV_DivMU.h"
#include "tmv/TMV_DivVU.h"
#include "tmv/TMV_SumMX.h"
#include "tmv/TMV_ScaleM.h"

#include "tmv/TMV_ProdVV.h"
#include "tmv/TMV_SumVV.h"
#include "tmv/TMV_SumMM.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_ProdMV.h"
#include "tmv/TMV_QuotMM.h"
#include "tmv/TMV_OProdVV.h"

namespace tmv {

    template <class T, int C>
    void InstQR_Downdate(UpperTriMatrixView<T> R, MatrixView<T,C> A)
    {
        TMVAssert(!R.isunit());
        UpperTriMatrixView<T,NonUnitDiag> Rn = R;
        InlineQR_Downdate(Rn,A); 
    }

#define InstFile "TMV_QRDowndate.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


