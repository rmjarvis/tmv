
//#define PRINTALGO_SVD
//#define XDEBUG_SVD

#include "TMV_Blas.h"
#include "tmv/TMV_SVDecompose.h"
#include "tmv/TMV_SVDecompose_Bidiag.h"
#include "tmv/TMV_SVDecompose_QR.h"
#include "tmv/TMV_SVDecompose_DC.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_ProdXV.h"
#include "tmv/TMV_ScaleV.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MinMax.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_SortV.h"
#include "tmv/TMV_MultPM.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_PermuteM.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_UnpackQ.h"
#include "tmv/TMV_SumMX.h"
#include "tmv/TMV_QRDecompose.h"

namespace tmv {

#ifdef LAP
    template <class Tu, class RT, class Tv> 
    static inline void LapSV_DecomposeFromBidiagonal(
        MatrixView<T> U, VectorView<RT>& D, VectorView<RT>& E,
        MatrixView<T> V, bool setUV)
    { InlineSV_DecomposeFromBidiagonal(U,D,E,V,setUV); }
#endif // LAP

    template <class Tu, class RT, class Tv>
    void InstSV_DecomposeFromBidiagonal(
        MatrixView<Tu> U, VectorView<RT> D, VectorView<RT> E,
        MatrixView<Tv> V, bool setUV)
    {
        TMVAssert(U.iscm());
        TMVAssert(V.iscm() || V.isrm());
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        if (D.size() > 0) {
#ifdef LAP
            LapSV_DecomposeFromBidiagonal(U,D,E,V,setUV);
#else
            MatrixView<Tu,ColMajor> Ucm = U;
            VectorView<RT,Unit> Du = D;
            VectorView<RT,Unit> Eu = E;
            if (V.iscm()) {
                MatrixView<Tv,ColMajor> Vcm = V;
                InlineSV_DecomposeFromBidiagonal(Ucm,Du,Eu,Vcm,setUV);
            } else {
                MatrixView<Tv,RowMajor> Vrm = V;
                InlineSV_DecomposeFromBidiagonal(Ucm,Du,Eu,Vrm,setUV);
            }
#endif
        }
    }

    template <class T, class RT>
    void InstSV_Decompose(
        MatrixView<T> U, DiagMatrixView<RT> S,
        MatrixView<T> V, T& signdet, RT& logdet, bool StoreU)
    {
        if (S.size() > 0) {
            if (U.iscm()) {
                if (V.isrm() || V.iscm() || !V.cptr()) {
                    MatrixView<T,ColMajor> Ucm = U;
                    if (V.iscm() || !V.cptr()) {
                        MatrixView<T,ColMajor> Vcm = V;
                        InlineSV_Decompose(Ucm,S,Vcm,signdet,logdet,StoreU); 
                    } else {
                        MatrixView<T,RowMajor> Vrm = V;
                        InlineSV_Decompose(Ucm,S,Vrm,signdet,logdet,StoreU); 
                    }
                } else {
                    InstSV_Decompose(U,S,V.copy().xView(),signdet,logdet,StoreU);
                } 
            } else {
                InstSV_Decompose(U.copy().xView(),S,V,signdet,logdet,StoreU);
            }
        }
    }


#define InstFile "TMV_SVDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


