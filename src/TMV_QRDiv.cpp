
//#define PRINTALGO_QR
//#define XDEBUG_QR

#include "TMV_Blas.h"
#include "tmv/TMV_QRDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SmallTriMatrix.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_DivVU.h"
#include "tmv/TMV_DivMU.h"
#include "tmv/TMV_PermuteM.h"
#include "tmv/TMV_TransposeM.h"
#include "tmv/TMV_ConjugateV.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_MultMM_Block.h"
#include "tmv/TMV_MultMM_Winograd.h"
#include "tmv/TMV_MultMM_OpenMP.h"
#include "tmv/TMV_Det.h"
#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_MultXU.h"


namespace tmv {

    template <class T1, int C1, class RT1, class T2>
    void InstQR_SolveInPlace(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, MatrixView<T2> m2)
    { InlineQR_SolveInPlace(QR,beta,P,N1,m2); }
    template <class T1, int C1, class RT1, class T2>
    void InstQR_SolveTransposeInPlace(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, MatrixView<T2> m2)
    { InlineQR_SolveTransposeInPlace(QR,beta,P,N1,m2); }
    template <class T1, int C1, class RT1, class T2>
    void InstQR_SolveInPlace(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, VectorView<T2> v2)
    { InlineQR_SolveInPlace(QR,beta,P,N1,v2); }
    template <class T1, int C1, class RT1, class T2>
    void InstQR_SolveTransposeInPlace(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, VectorView<T2> v2)
    { InlineQR_SolveTransposeInPlace(QR,beta,P,N1,v2); }



    //
    // The Solve functions:
    //

    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstQR_Solve(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineQR_Solve(QR,beta,P,N1,m2,m3); }
    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstQR_SolveTranspose(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineQR_SolveTranspose(QR,beta,P,N1,m2,m3); }
    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstQR_Solve(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    { InlineQR_Solve(QR,beta,P,N1,v2,v3); }
    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstQR_SolveTranspose(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    { InlineQR_SolveTranspose(QR,beta,P,N1,v2,v3); }

#define InstFile "TMV_QRDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


