
//#define PRINTALGO_QR
//#define XDEBUG_QR
//#define XDEBUG_MM
//#define XDEBUG_MV

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

    template <class M, class V> 
    static inline void DoUnpackQ(M& Q, const V& beta)
    {
        TMVAssert(Q.iscm());
        TMVAssert(beta.step() == 1);
        typename M::cmview_type Qcm = Q.cmView();
        typename V::const_unitview_type betau = beta.unitView();
        InlineUnpackQ(Qcm,betau);
    }

#ifdef LAP
#ifdef TMV_INST_DOUBLE
    void DoUnpackQ(
        MatrixView<double> Q, const ConstVectorView<double>& beta)
    {
        int m = Q.colsize();
        int n = Q.rowsize();
        if (Q.iscm()) {
            int ldq = Q.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            LAPNAME(dorgqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(dorgqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"dorgqr");
#else
            LAP_Results(Lap_info,int(work[0]),m,n,lwork,"dorgqr");
#endif
        } else {
            int ldq = Q.stepi();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            LAPNAME(dorglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(dorglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"dorglq");
#else
            LAP_Results(Lap_info,int(work[0]),m,n,lwork,"dorglq");
#endif
        }
    }
    void DoUnpackQ(
        MatrixView<std::complex<double> > Q,
        const ConstVectorView<double>& beta)
    {
        int m = Q.colsize();
        int n = Q.rowsize();
        Vector<std::complex<double> > cbeta = beta;
        if (Q.iscm()) {
            int ldq = Q.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            work.get()[0] = 0.;
            LAPNAME(zungqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(cbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(zungqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(cbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"zungqr");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"zungqr");
#endif
        } else {
            int ldq = Q.stepi();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            work.get()[0] = 0.;
            LAPNAME(zunglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(cbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(zunglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(cbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"zunglq");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"zunglq");
#endif
        }
    }
#endif
#ifdef TMV_INST_FLOAT
    void DoUnpackQ(
        MatrixView<float> Q, const ConstVectorView<float>& beta)
    {
        int m = Q.colsize();
        int n = Q.rowsize();
        if (Q.iscm()) {
            int ldq = Q.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.F;
            LAPNAME(sorgqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(sorgqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"sorgqr");
#else
            LAP_Results(Lap_info,int(work[0]),m,n,lwork,"sorgqr");
#endif
        } else {
            int ldq = Q.stepi();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.F;
            LAPNAME(sorglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(sorglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"sorglq");
#else
            LAP_Results(Lap_info,int(work[0]),m,n,lwork,"sorglq");
#endif
        }
    }
    void DoUnpackQ(
        MatrixView<std::complex<float> > Q,
        const ConstVectorView<float>& beta)
    {
        int m = Q.colsize();
        int n = Q.rowsize();
        Vector<std::complex<float> > cbeta = beta;
        if (Q.iscm()) {
            int ldq = Q.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            work.get()[0] = 0.F;
            LAPNAME(cungqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(cbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(cungqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(cbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"cungqr");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"cungqr");
#endif
        } else {
            int ldq = Q.stepi();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            work.get()[0] = 0.F;
            LAPNAME(cunglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(cbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(cunglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(cbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"cunglq");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"cunglq");
#endif
        }
    }
#endif
#endif // LAP

    template <class T, class RT> 
    void InstUnpackQ(MatrixView<T> Q, const ConstVectorView<RT>& beta)
    {
        if (beta.step() == 1) {
            if (Q.iscm()) {
                DoUnpackQ(Q,beta);
            } else {
                Matrix<T,ColMajor|NoDivider> Qc = Q;
                InstUnpackQ(Qc.xView(),beta);
                InstCopy(Qc.constView().xView(),Q);
            }
        } else {
            Vector<RT> betac = beta;
            InstUnpackQ(Q,betac.constView().xView());
        }
    }

#define InstFile "TMV_UnpackQ.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


