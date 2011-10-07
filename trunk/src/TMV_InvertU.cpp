
#include "tmv/TMV_InvertU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Det.h" // For isSingular().
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_Vector.h"

namespace tmv {

    template <class T>
    void DoInvertSelf(UpperTriMatrixView<T> m)
    {
        TMVAssert(m.iscm() || m.isrm());
        if (m.iscm()) {
            UpperTriMatrixView<T,ColMajor> mcm = m.cmView();
            InlineInvertSelf(mcm);
        } else {
            UpperTriMatrixView<T,RowMajor> mrm = m.rmView();
            InlineInvertSelf(mrm);
        }
    }

#ifdef LAP
#ifdef TMV_INST_DOUBLE
    void DoInvertSelf(const UpperTriMatrixView<double>& m)
    {
        int n = m.size();
        int lda = m.iscm() ? m.stepj() : m.stepi();
        LAPNAME(dtrtri) (
            LAPCM m.iscm()?LAPCH_UP:LAPCH_LO,
            m.isunit()?LAPCH_U:LAPCH_NU,LAPV(n),LAPP(m.ptr()),LAPV(lda) LAPINFO
            LAP1 LAP1);
        LAP_Results("dtrtri");
    }
    void DoInvertSelf(
        const UpperTriMatrixView<std::complex<double> >& m)
    {
        int n = m.size();
        int lda = m.iscm() ? m.stepj() : m.stepi();
        LAPNAME(ztrtri) (
            LAPCM m.iscm()?LAPCH_UP:LAPCH_LO,
            m.isunit()?LAPCH_U:LAPCH_NU,LAPV(n),LAPP(m.ptr()),LAPV(lda) LAPINFO
            LAP1 LAP1);
        LAP_Results("ztrtri");
    }
#endif
#ifdef TMV_INST_FLOAT
    void DoInvertSelf(const UpperTriMatrixView<float>& m)
    {
        int n = m.size();
        int lda = m.iscm() ? m.stepj() : m.stepi();
        LAPNAME(strtri) (
            LAPCM m.iscm()?LAPCH_UP:LAPCH_LO,
            m.isunit()?LAPCH_U:LAPCH_NU,LAPV(n),LAPP(m.ptr()),LAPV(lda) LAPINFO
            LAP1 LAP1);
        LAP_Results("strtri");
    }
    void DoInvertSelf(
        const UpperTriMatrixView<std::complex<float> >& m)
    {
        int n = m.size();
        int lda = m.iscm() ? m.stepj() : m.stepi();
        LAPNAME(ctrtri) (
            LAPCM m.iscm()?LAPCH_UP:LAPCH_LO,
            m.isunit()?LAPCH_U:LAPCH_NU,LAPV(n),LAPP(m.ptr()),LAPV(lda) LAPINFO
            LAP1 LAP1);
        LAP_Results("ctrtri");
    }
#endif // FLOAT
#endif // ALAP

    template <class T>
    void InstInvertSelf(UpperTriMatrixView<T> m)
    {
        if (m.size() > 0) {
            if (m.iscm() || m.isrm()) {
                DoInvertSelf(m);
            } else if (m.isunit()) {
                UpperTriMatrix<T,UnitDiag|ColMajor> mc = m;
                DoInvertSelf(mc.xView());
                InstCopy(mc.constView().xView(),m);
            } else {
                UpperTriMatrix<T,NonUnitDiag|ColMajor> mc = m;
                DoInvertSelf(mc.xView());
                InstCopy(mc.constView().xView(),m);
            }
        }
    }

    template <class T>
    void InstInvertSelf(LowerTriMatrixView<T> m)
    { InstInvertSelf(m.transpose()); }

#define InstFile "TMV_InvertU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


