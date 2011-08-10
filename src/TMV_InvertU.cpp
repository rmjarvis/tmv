
#include "tmv/TMV_InvertU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Det.h" // For isSingular().
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_Vector.h"

namespace tmv {

    template <class T>
    void InstInvertSelf(UpperTriMatrixView<T> m)
    {
        if (m.iscm()) {
            UpperTriMatrixView<T,ColMajor> mcm = m.cmView();
            InlineInvertSelf(mcm);
        } else if (m.isrm()) {
            UpperTriMatrixView<T,RowMajor> mrm = m.rmView();
            InlineInvertSelf(mrm);
        } else if (m.isunit()) {
            UpperTriMatrix<T,UnitDiag|ColMajor> mc = m;
            UpperTriMatrixView<T,ColMajor> mcm = mc.cmView();
            InlineInvertSelf(mcm);
            InstCopy(mc.constView().xView(),m);
        } else {
            UpperTriMatrix<T,NonUnitDiag|ColMajor> mc = m;
            UpperTriMatrixView<T,ColMajor> mcm = mc.cmView();
            InlineInvertSelf(mcm);
            InstCopy(mc.constView().xView(),m);
        }
    }

    template <class T>
    void InstInvertSelf(LowerTriMatrixView<T> m)
    { InstInvertSelf(m.transpose()); }

#define InstFile "TMV_InvertU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


