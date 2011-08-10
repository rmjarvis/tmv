
//#define PRINTALGO_QR

#include "TMV_Blas.h"
#include "tmv/TMV_QRDecompose.h"
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
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultUM.h"

namespace tmv {

#ifdef LAP
    template <class T, int S, class RT> 
    static inline void LapQR_Decompose(
        MatrixView<T,S> A, VectorView<RT,Unit> beta)
    { InlineQR_Decompose(A,beta); }
#endif // LAP

    template <class T, class RT> 
    void InstQR_Decompose(MatrixView<T> A, VectorView<RT> beta)
    {
        if (A.rowsize() > 0) {
            if (beta.step() == 1) {
                VectorView<RT,Unit> beta1 = beta;
                if (A.iscm()) {
                    MatrixView<T,ColMajor> Acm = A;
#ifdef LAP
                    LapQR_Decompose(Acm,beta1);
#else
                    InlineQR_Decompose(Acm,beta1);
#endif
                } else if (A.isrm()) {
                    MatrixView<T,RowMajor> Arm = A;
#ifdef LAP
                    LapQR_Decompose(Arm,beta1);
#else
                    InlineQR_Decompose(Arm,beta1);
#endif
                } else {
                    Matrix<T,ColMajor|NoDivider> Ac = A;
                    InstQR_Decompose(Ac.xView(),beta);
                    InstCopy(Ac.constView().xView(),A);
                }
            } else {
                Vector<RT> betac = beta;
                InstQR_Decompose(A,betac.xView());
                InstCopy(betac.constView().xView(),beta);
            }
        }
    }

#define InstFile "TMV_QRDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


