
//#define PRINTALGO_QR

#include "TMV_Blas.h"
#include "tmv/TMV_UnpackQ.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_SwapV.h"
#include "tmv/TMV_MinMax.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultUM.h"
#include "tmv/TMV_MultUU.h"
#include "tmv/TMV_MultUL.h"

namespace tmv {

#ifdef LAP
    template <class T, int S, class RT> 
    static inline void LapUnpackQ(
        MatrixView<T,S> A, const ConstVectorView<RT,Unit>& beta)
    { InlineUnpackQ(A,beta); }
#endif // LAP

    template <class T, class RT> 
    void InstUnpackQ(MatrixView<T> A, const ConstVectorView<RT>& beta)
    {
        if (A.rowsize() > 0) {
            if (beta.step() == 1) {
                ConstVectorView<RT,Unit> beta1 = beta;
                if (A.iscm()) {
                    MatrixView<T,ColMajor> Acm = A;
#ifdef LAP
                    LapUnpackQ(Acm,beta1);
#else
                    InlineUnpackQ(Acm,beta1);
#endif
                } else if (A.isrm()) {
                    MatrixView<T,RowMajor> Arm = A;
#ifdef LAP
                    LapUnpackQ(Arm,beta1);
#else
                    InlineUnpackQ(Arm,beta1);
#endif
                } else {
                    Matrix<T,ColMajor|NoDivider> Ac = A;
                    InstUnpackQ(Ac.xView(),beta);
                    InstCopy(Ac.constView().xView(),A);
                }
            } else {
                Vector<RT> betac = beta;
                InstUnpackQ(A,betac.constView().xView());
            }
        }
    }

#define InstFile "TMV_UnpackQ.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


