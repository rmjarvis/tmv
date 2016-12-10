
//#define PRINTALGO_SVD

#include "tmv/TMV_SVD.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_SmallDiagMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_SumMM.h"
#include "tmv/TMV_AddMM.h"
#include "tmv/TMV_SVInverse.h"
#include "tmv/TMV_SVDecompose.h"
#include "tmv/TMV_SVDiv.h"
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
#include "tmv/TMV_TransposeM.h"
#include "tmv/TMV_DivM.h"
#include "tmv/TMV_VectorIO.h"

namespace tmv {


    template <class T> template <int C>
    InstSVD<T>::InstSVD(const ConstMatrixView<T,C>& A, bool _inplace) : 
        base(A,_inplace) {}
    template <class T> 
    InstSVD<T>::InstSVD(const InstSVD<T>& rhs) : base(rhs) {}
    template <class T> 
    InstSVD<T>::~InstSVD() {}

    template <class T> 
    void InstSVD<T>::doSolveInPlace(MatrixView<RT> m2) const 
    { base::solveInPlace(m2);  }
    template <class T> 
    void InstSVD<T>::doSolveInPlace(MatrixView<CT> m2) const 
    { base::solveInPlace(m2);  }
    template <class T> 
    void InstSVD<T>::doSolveInPlace(MatrixView<CT,Conj> m2) const 
    { base::solveInPlace(m2);  }
    template <class T> 
    void InstSVD<T>::doSolveInPlace(VectorView<RT> v2) const 
    { base::solveInPlace(v2);  }
    template <class T> 
    void InstSVD<T>::doSolveInPlace(VectorView<CT> v2) const 
    { base::solveInPlace(v2);  }
    template <class T> 
    void InstSVD<T>::doSolveInPlace(VectorView<CT,Conj> v2) const 
    { base::solveInPlace(v2);  }

    template <class T> 
    void InstSVD<T>::doSolveTransposeInPlace(MatrixView<RT> m2) const 
    { base::solveTransposeInPlace(m2); } 
    template <class T> 
    void InstSVD<T>::doSolveTransposeInPlace(MatrixView<CT> m2) const 
    { base::solveTransposeInPlace(m2); } 
    template <class T> 
    void InstSVD<T>::doSolveTransposeInPlace(MatrixView<CT,Conj> m2) const 
    { base::solveTransposeInPlace(m2); } 
    template <class T> 
    void InstSVD<T>::doSolveTransposeInPlace(VectorView<RT> v2) const 
    { base::solveTransposeInPlace(v2); } 
    template <class T> 
    void InstSVD<T>::doSolveTransposeInPlace(VectorView<CT> v2) const 
    { base::solveTransposeInPlace(v2); } 
    template <class T> 
    void InstSVD<T>::doSolveTransposeInPlace(VectorView<CT,Conj> v2) const 
    { base::solveTransposeInPlace(v2); } 

    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstMatrixView<RT>& m1, MatrixView<CT,Conj> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstMatrixView<CT>& m1, MatrixView<CT,Conj> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstMatrixView<CT,Conj>& m1, MatrixView<CT> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstMatrixView<CT,Conj>& m1, MatrixView<CT,Conj> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstVectorView<RT>& v1, VectorView<RT> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstVectorView<RT>& v1, VectorView<CT> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstVectorView<RT>& v1, VectorView<CT,Conj> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstVectorView<CT>& v1, VectorView<CT> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstVectorView<CT>& v1, VectorView<CT,Conj> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstVectorView<CT,Conj>& v1, VectorView<CT> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstSVD<T>::doSolve(
        const ConstVectorView<CT,Conj>& v1, VectorView<CT,Conj> v2) const 
    { base::solve(v1,v2); } 

    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstMatrixView<RT>& m1, MatrixView<CT,Conj> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstMatrixView<CT>& m1, MatrixView<CT,Conj> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstMatrixView<CT,Conj>& m1, MatrixView<CT> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstMatrixView<CT,Conj>& m1, MatrixView<CT,Conj> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstVectorView<RT>& v1, VectorView<RT> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstVectorView<RT>& v1, VectorView<CT> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstVectorView<RT>& v1, VectorView<CT,Conj> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstVectorView<CT>& v1, VectorView<CT> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstVectorView<CT>& v1, VectorView<CT,Conj> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstVectorView<CT,Conj>& v1, VectorView<CT> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstSVD<T>::doSolveTranspose(
        const ConstVectorView<CT,Conj>& v1, VectorView<CT,Conj> v2) const 
    { base::solveTranspose(v1,v2); }

    template <class T> 
    T InstSVD<T>::det() const
    { return base::det(); }
    template <class T> 
    typename InstSVD<T>::FT InstSVD<T>::logDet(
        typename InstSVD<T>::ZFT* sign) const
    { return base::logDet(sign); }
    template <class T> 
    bool InstSVD<T>::isSingular() const
    { return base::isSingular(); }

    template <class T> 
    void InstSVD<T>::doMakeInverse(MatrixView<RT> minv) const 
    { base::makeInverse(minv); } 
    template <class T> 
    void InstSVD<T>::doMakeInverse(MatrixView<CT> minv) const 
    { base::makeInverse(minv); } 
    template <class T> 
    void InstSVD<T>::doMakeInverse(MatrixView<CT,Conj> minv) const 
    { base::makeInverse(minv); } 

    template <class T> 
    void InstSVD<T>::doMakeInverseATA(MatrixView<RT> ata) const
    { base::makeInverseATA(ata); }
    template <class T> 
    void InstSVD<T>::doMakeInverseATA(MatrixView<CT> ata) const
    { base::makeInverseATA(ata); }
    template <class T> 
    void InstSVD<T>::doMakeInverseATA(MatrixView<CT,Conj> ata) const
    { base::makeInverseATA(ata); }

    template <class T> 
    typename InstSVD<T>::RT InstSVD<T>::condition() const
    { return base::condition(); }

    template <class T> 
    bool InstSVD<T>::preferInPlace() const
    { return base::preferInPlace(); }


#define InstFile "TMV_SVD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


