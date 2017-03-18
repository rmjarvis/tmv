
//#define PRINTALGO_QR

#include "tmv/TMV_QRPD.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SmallTriMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_MultUL.h"
#include "tmv/TMV_MultPM.h"
#include "tmv/TMV_MultPV.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_SumMM.h"
#include "tmv/TMV_AddMM.h"
#include "tmv/TMV_QRInverse.h"
#include "tmv/TMV_QRPDecompose.h"
#include "tmv/TMV_QRDiv.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_ConjugateV.h"
#include "tmv/TMV_Det.h"
#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_MultMM_Block.h"
#include "tmv/TMV_MultMM_Winograd.h"
#include "tmv/TMV_MultMM_OpenMP.h"
#include "tmv/TMV_PermuteM.h"
#include "tmv/TMV_TransposeM.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_DivM.h"
#include "tmv/TMV_DivMU.h"
#include "tmv/TMV_MultXU.h"

namespace tmv {


    template <class T> template <int C>
    InstQRPD<T>::InstQRPD(const ConstMatrixView<T,C>& A, bool _inplace) : 
        base(A,_inplace) {}
    template <class T> 
    InstQRPD<T>::InstQRPD(const InstQRPD<T>& rhs) : base(rhs) {}
    template <class T> 
    InstQRPD<T>::~InstQRPD() {}

    template <class T> 
    void InstQRPD<T>::doSolveInPlace(MatrixView<RT> m2) const 
    { base::solveInPlace(m2);  }
    template <class T> 
    void InstQRPD<T>::doSolveInPlace(MatrixView<CT> m2) const 
    { base::solveInPlace(m2);  }
    template <class T> 
    void InstQRPD<T>::doSolveInPlace(MatrixView<CT,Conj> m2) const 
    { base::solveInPlace(m2);  }
    template <class T> 
    void InstQRPD<T>::doSolveInPlace(VectorView<RT> v2) const 
    { base::solveInPlace(v2);  }
    template <class T> 
    void InstQRPD<T>::doSolveInPlace(VectorView<CT> v2) const 
    { base::solveInPlace(v2);  }
    template <class T> 
    void InstQRPD<T>::doSolveInPlace(VectorView<CT,Conj> v2) const 
    { base::solveInPlace(v2);  }

    template <class T> 
    void InstQRPD<T>::doSolveTransposeInPlace(MatrixView<RT> m2) const 
    { base::solveTransposeInPlace(m2); } 
    template <class T> 
    void InstQRPD<T>::doSolveTransposeInPlace(MatrixView<CT> m2) const 
    { base::solveTransposeInPlace(m2); } 
    template <class T> 
    void InstQRPD<T>::doSolveTransposeInPlace(MatrixView<CT,Conj> m2) const 
    { base::solveTransposeInPlace(m2); } 
    template <class T> 
    void InstQRPD<T>::doSolveTransposeInPlace(VectorView<RT> v2) const 
    { base::solveTransposeInPlace(v2); } 
    template <class T> 
    void InstQRPD<T>::doSolveTransposeInPlace(VectorView<CT> v2) const 
    { base::solveTransposeInPlace(v2); } 
    template <class T> 
    void InstQRPD<T>::doSolveTransposeInPlace(VectorView<CT,Conj> v2) const 
    { base::solveTransposeInPlace(v2); } 

    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstMatrixView<RT>& m1, MatrixView<CT,Conj> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstMatrixView<CT>& m1, MatrixView<CT,Conj> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstMatrixView<CT,Conj>& m1, MatrixView<CT> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstMatrixView<CT,Conj>& m1, MatrixView<CT,Conj> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstVectorView<RT>& v1, VectorView<RT> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstVectorView<RT>& v1, VectorView<CT> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstVectorView<RT>& v1, VectorView<CT,Conj> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstVectorView<CT>& v1, VectorView<CT> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstVectorView<CT>& v1, VectorView<CT,Conj> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstVectorView<CT,Conj>& v1, VectorView<CT> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstQRPD<T>::doSolve(
        const ConstVectorView<CT,Conj>& v1, VectorView<CT,Conj> v2) const 
    { base::solve(v1,v2); } 

    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstMatrixView<RT>& m1, MatrixView<CT,Conj> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstMatrixView<CT>& m1, MatrixView<CT,Conj> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstMatrixView<CT,Conj>& m1, MatrixView<CT> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstMatrixView<CT,Conj>& m1, MatrixView<CT,Conj> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstVectorView<RT>& v1, VectorView<RT> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstVectorView<RT>& v1, VectorView<CT> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstVectorView<RT>& v1, VectorView<CT,Conj> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstVectorView<CT>& v1, VectorView<CT> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstVectorView<CT>& v1, VectorView<CT,Conj> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstVectorView<CT,Conj>& v1, VectorView<CT> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstQRPD<T>::doSolveTranspose(
        const ConstVectorView<CT,Conj>& v1, VectorView<CT,Conj> v2) const 
    { base::solveTranspose(v1,v2); }

    template <class T> 
    T InstQRPD<T>::det() const
    { return base::det(); }
    template <class T> 
    typename InstQRPD<T>::FT InstQRPD<T>::logDet(
        typename InstQRPD<T>::ZFT* sign) const
    { return base::logDet(sign); }
    template <class T> 
    bool InstQRPD<T>::isSingular() const
    { return base::isSingular(); }

    template <class T> 
    void InstQRPD<T>::doMakeInverse(MatrixView<RT> minv) const 
    { base::makeInverse(minv); } 
    template <class T> 
    void InstQRPD<T>::doMakeInverse(MatrixView<CT> minv) const 
    { base::makeInverse(minv); } 
    template <class T> 
    void InstQRPD<T>::doMakeInverse(MatrixView<CT,Conj> minv) const 
    { base::makeInverse(minv); } 

    template <class T> 
    void InstQRPD<T>::doMakeInverseATA(MatrixView<RT> ata) const
    { base::makeInverseATA(ata); }
    template <class T> 
    void InstQRPD<T>::doMakeInverseATA(MatrixView<CT> ata) const
    { base::makeInverseATA(ata); }
    template <class T> 
    void InstQRPD<T>::doMakeInverseATA(MatrixView<CT,Conj> ata) const
    { base::makeInverseATA(ata); }

    template <class T> 
    typename InstQRPD<T>::RT InstQRPD<T>::condition(
        typename InstQRPD<T>::RT normInf) const
    { return base::condition(normInf); }

    template <class T> 
    bool InstQRPD<T>::preferInPlace() const
    { return base::preferInPlace(); }


#define InstFile "TMV_QRPD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


