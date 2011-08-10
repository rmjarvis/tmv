
//#define PRINTALGO_SVD
//#define XDEBUG_SVD

#include "tmv/TMV_SVDecompose_Bidiag.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_MultMM.h"

namespace tmv {

    template <class T, class RT>
    void InstBidiagonalize(
        MatrixView<T> A, VectorView<RT> Ubeta, VectorView<RT> Vbeta,
        VectorView<T> D, VectorView<T> E)
    {
        TMVAssert(A.iscm());
        TMVAssert(Ubeta.step() == 1);
        TMVAssert(Vbeta.step() == 1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        if (D.size() > 0) {
            MatrixView<T,ColMajor> Acm = A;
            VectorView<RT,Unit> Ubu = Ubeta;
            VectorView<RT,Unit> Vbu = Vbeta;
            VectorView<T,Unit> Du = D;
            VectorView<T,Unit> Eu = E;
            InlineBidiagonalize(Acm,Ubu,Vbu,Du,Eu);
        }
    }

#define InstFile "TMV_SVDecompose_Bidiag.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


