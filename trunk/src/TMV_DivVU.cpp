
//#define PRINTALGO_DIVVU

#include "tmv/TMV_DivVU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_Det.h"

namespace tmv {

    template <class T1, class M2>
    static void DoTriLDivEq(VectorView<T1,Unit> v1, const M2& m2)
    {
        const ptrdiff_t xx = Unknown;
        if (m2.iscm()) 
            InlineTriLDivEq(v1,m2.cmView());
        else if (m2.isrm())
            InlineTriLDivEq(v1,m2.rmView());
        else {
            typedef typename M2::value_type T;
            const ptrdiff_t N = m2.size();
            if (m2.isunit()) {
                const ptrdiff_t s = ShapeTraits<M2::_shape>::unit_shape;
                typename MCopyHelper<T,s,xx,xx>::type mc(N);
                InstCopy(m2,mc.xView());
                InlineTriLDivEq(v1,mc.xView().constView().cmView());
            } else  {
                const ptrdiff_t s = ShapeTraits<M2::_shape>::nonunit_shape;
                typename MCopyHelper<T,s,xx,xx>::type mc(N);
                InstCopy(m2,mc.xView());
                InlineTriLDivEq(v1,mc.xView().constView().cmView());
            }
        }
    }

#ifdef BLAS
    template <class V1, class M2, class T2>
    static void BlasTriLDivEq(const V1& b, M2& A, T2)
    { DoTriLDivEq(b.unitView(),A); }
#ifdef TMV_INST_DOUBLE
    template <class M2> 
    static void BlasTriLDivEq(VectorView<double> b, const M2& A, double)
    {
        int n=A.size();
        int lda = A.iscm()?A.stepj():A.stepi();
        int bs = b.step();
        double* bp = b.ptr();
        if (bs < 0) bp += (n-1)*bs;
        BLASNAME(dtrsv) (
            BLASCM A.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
    }
    template <class M2> 
    static void BlasTriLDivEq(
        VectorView<std::complex<double> > b, const M2& A, std::complex<double>)
    {
        int n=A.size();
        int lda = A.iscm()?A.stepj():A.stepi();
        int bs = b.step();
        std::complex<double>* bp = b.ptr();
        if (bs < 0) bp += (n-1)*bs;
        if (A.iscm() && A.isconj()) {
            b.conjugateSelf();
            BLASNAME(ztrsv) (
                BLASCM M2::mupper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
            b.conjugateSelf();
        } else {
            BLASNAME(ztrsv) (
                BLASCM A.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO,
                A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
        }
    }
    template <class M2> 
    static void BlasTriLDivEq(
        VectorView<std::complex<double> > b, const M2& A, double)
    {
        int n=A.size();
        int lda = A.iscm()?A.stepj():A.stepi();
        int bs = 2*b.step();
        double* bp = (double*)(b.ptr());
        if (bs < 0) bp += (n-1)*bs;
        BLASNAME(dtrsv) (
            BLASCM A.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
        BLASNAME(dtrsv) (
            BLASCM A.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp+1),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
    }
#endif
#ifdef TMV_INST_FLOAT
    template <class M2> 
    static void BlasTriLDivEq(VectorView<float> b, const M2& A, float)
    {
        int n=A.size();
        int lda = A.iscm()?A.stepj():A.stepi();
        int bs = b.step();
        float* bp = b.ptr();
        if (bs < 0) bp += (n-1)*bs;
        BLASNAME(strsv) (
            BLASCM A.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
    }
    template <class M2> 
    static void BlasTriLDivEq(
        VectorView<std::complex<float> > b, const M2& A, std::complex<float>)
    {
        int n=A.size();
        int lda = A.iscm()?A.stepj():A.stepi();
        int bs = b.step();
        std::complex<float>* bp = b.ptr();
        if (bs < 0) bp += (n-1)*bs;
        if (A.iscm() && A.isconj()) {
            b.conjugateSelf();
            BLASNAME(ctrsv) (
                BLASCM M2::mupper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
            b.conjugateSelf();
        } else {
            BLASNAME(ctrsv) (
                BLASCM A.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO,
                A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
        }
    }
    template <class M2> 
    static void BlasTriLDivEq(
        VectorView<std::complex<float> > b, const M2& A, float)
    {
        int n=A.size();
        int lda = A.iscm()?A.stepj():A.stepi();
        int bs = 2*b.step();
        float* bp = (float*)(b.ptr());
        if (bs < 0) bp += (n-1)*bs;
        BLASNAME(strsv) (
            BLASCM A.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
        BLASNAME(strsv) (
            BLASCM A.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp+1),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
    }
#endif
#endif // BLAS

    template <class T1, class M2>
    static inline void CallTriLDivEq(VectorView<T1> v1, const M2& m2)
    {
#ifdef BLAS
        const typename M2::value_type t2(0);
        if ( (m2.iscm() && m2.stepj()>0) || (m2.isrm() && m2.stepi()>0) ) {
            BlasTriLDivEq(v1,m2,t2);
        } else {
            if (m2.isunit()) {
                BlasTriLDivEq(
                    v1,m2.copy().viewAsUnitDiag().constView().xView(),t2);
            } else {
                BlasTriLDivEq(v1,m2.copy().constView().xView(),t2);
            }
        }
#else
        if (v1.step() == 1) {
            DoTriLDivEq(v1.unitView(),m2);
        } else {
            Vector<T1> v1c(v1.size());
            InstCopy(v1.constView(),v1c.xView());
            DoTriLDivEq(v1c.view(),m2);
            InstCopy(v1c.xView().constView(),v1);
        }
#endif
    }

    template <class T1, class T2, int C2>
    void InstTriLDivEq(
        VectorView<T1> v1, const ConstUpperTriMatrixView<T2,C2>& m2)
    { CallTriLDivEq(v1,m2); }
    template <class T1, class T2, int C2>
    void InstTriLDivEq(
        VectorView<T1> v1, const ConstLowerTriMatrixView<T2,C2>& m2)
    { CallTriLDivEq(v1,m2); }

    template <class T1, class T2, int C2>
    void InstAliasTriLDivEq(
        VectorView<T1> v1, const ConstUpperTriMatrixView<T2,C2>& m2)
    { InlineAliasTriLDivEq(v1,m2); }
    template <class T1, class T2, int C2>
    void InstAliasTriLDivEq(
        VectorView<T1> v1, const ConstLowerTriMatrixView<T2,C2>& m2)
    { InlineAliasTriLDivEq(v1,m2); }


#define InstFile "TMV_DivVU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


