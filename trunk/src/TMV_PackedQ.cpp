
//#define PRINTALGO_QR
//#define XDEBUG_QR

#include "TMV_Blas.h"
#include "tmv/TMV_PackedQ.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_ConjugateV.h"

#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_ProdMV.h"
#include "tmv/TMV_SumVV.h"
#include "tmv/TMV_OProdVV.h"

#ifdef XDEBUG_QR
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#include "TMV_Vector.h"
#include "TMV_Matrix.h"
#endif

namespace tmv {

    template <class M1, class V1, class M2>
    static void DoPackedQ_MultEq(
        const M1& Q, const V1& beta, M2& m2)
    { InlinePackedQ_MultEq(Q,beta.unitView(),m2); }

    template <class M1, class V1, class M2>
    static void DoPackedQ_LDivEq(
        const M1& Q, const V1& beta, M2& m2)
    { InlinePackedQ_LDivEq(Q,beta.unitView(),m2); }

#ifdef LAP
#ifdef TMV_INST_DOUBLE
    static void DoPackedQ_LDivEq(
        const ConstMatrixView<double>& Q,
        const ConstVectorView<double>& beta,
        MatrixView<double> x)
    {
        int ldx = x.isrm() ? x.stepi() : x.stepj();
        int m = x.isrm() ? x.rowsize() : x.colsize();
        int n = x.isrm() ? x.colsize() : x.rowsize();
        int k = Q.rowsize();
        if (Q.isrm()) {
            int ldq = Q.stepi();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            LAPNAME(dormlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_T:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(dormlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_T:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("dormlq");
#else
            LAP_Results(int(work[0]),m,n,lwork,"dormlq");
#endif
        } else {
            int ldq = Q.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            LAPNAME(dormqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
                x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(dormqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
                x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("dormqr");
#else
            LAP_Results(int(work[0]),m,n,lwork,"dormqr");
#endif
        }
    }
    template <int C>
    void DoPackedQ_LDivEq(
        const ConstMatrixView<std::complex<double>,C>& Q,
        const ConstVectorView<double>& beta,
        MatrixView<std::complex<double> > x)
    {
        int k = Q.rowsize();
        Vector<std::complex<double> > cbeta = beta;
        if (Q.isrm()) {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepi();
            if (x.isrm() == Q.isconj()) x.conjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            work.get()[0] = 0.;
            LAPNAME(zunmlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(zunmlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() == Q.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("zunmlq");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"zunmlq");
#endif
        } else {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepj();
            if (x.isrm() != Q.isconj()) x.conjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            work.get()[0] = 0.;
            LAPNAME(zunmqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(zunmqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() != Q.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("zunmqr");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"zunmqr");
#endif
        }
    }
    void DoPackedQ_MultEq(
        const ConstMatrixView<double>& Q,
        const ConstVectorView<double>& beta,
        MatrixView<double> x)
    {
        int ldx = x.isrm() ? x.stepi() : x.stepj();
        int m = x.isrm() ? x.rowsize() : x.colsize();
        int n = x.isrm() ? x.colsize() : x.rowsize();
        int k = Q.rowsize();
        if (Q.isrm()) {
            int ldq = Q.stepi();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            LAPNAME(dormlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_T,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(dormlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_T,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("dormlq");
#else
            LAP_Results(int(work[0]),m,n,lwork,"dormlq");
#endif
        } else {
            int ldq = Q.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            LAPNAME(dormqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
                x.isrm()?LAPCH_T:LAPCH_NT, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(dormqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
                x.isrm()?LAPCH_T:LAPCH_NT, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("dormqr");
#else
            LAP_Results(int(work[0]),m,n,lwork,"dormqr");
#endif
        }
    }
    template <int C>
    void DoPackedQ_MultEq(
        const ConstMatrixView<std::complex<double>,C>& Q,
        const ConstVectorView<double>& beta,
        MatrixView<std::complex<double> > x)
    {
        int k = Q.rowsize();
        Vector<std::complex<double> > cbeta = beta;
        if (Q.isrm()) {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepi();
            if (x.isrm() == Q.isconj()) x.conjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            work.get()[0] = 0.;
            LAPNAME(zunmlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(zunmlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() == Q.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("zunmlq");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"zunmlq");
#endif
        } else {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepj();
            if (x.isrm() != Q.isconj()) x.conjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            work.get()[0] = 0.;
            LAPNAME(zunmqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(zunmqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() != Q.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("zunmqr");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"zunmqr");
#endif
        }
    }
    static void DoPackedQ_LDivEq(
        const ConstMatrixView<double>& Q,
        const ConstVectorView<double>& beta,
        VectorView<double> x)
    { DoPackedQ_LDivEq(Q,beta,ColVectorViewOf(x).xView()); }
    template <int C>
    void DoPackedQ_LDivEq(
        const ConstMatrixView<std::complex<double>,C>& Q,
        const ConstVectorView<double>& beta,
        VectorView<std::complex<double> > x)
    { DoPackedQ_LDivEq(Q,beta,ColVectorViewOf(x).xView()); }
    static void DoPackedQ_MultEq(
        const ConstMatrixView<double>& Q,
        const ConstVectorView<double>& beta,
        VectorView<double> x)
    { DoPackedQ_MultEq(Q,beta,ColVectorViewOf(x).xView()); }
    template <int C>
    void DoPackedQ_MultEq(
        const ConstMatrixView<std::complex<double>,C>& Q,
        const ConstVectorView<double>& beta,
        VectorView<std::complex<double> > x)
    { DoPackedQ_MultEq(Q,beta,ColVectorViewOf(x).xView()); }
#endif
#ifdef TMV_INST_FLOAT
    static void DoPackedQ_LDivEq(
        const ConstMatrixView<float>& Q,
        const ConstVectorView<float>& beta,
        MatrixView<float> x)
    {
        int ldx = x.isrm() ? x.stepi() : x.stepj();
        int m = x.isrm() ? x.rowsize() : x.colsize();
        int n = x.isrm() ? x.colsize() : x.rowsize();
        int k = Q.rowsize();
        if (Q.isrm()) {
            int ldq = Q.stepi();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.;
            LAPNAME(sormlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_T:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(sormlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_T:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("sormlq");
#else
            LAP_Results(int(work[0]),m,n,lwork,"sormlq");
#endif
        } else {
            int ldq = Q.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.;
            LAPNAME(sormqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
                x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(sormqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
                x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("sormqr");
#else
            LAP_Results(int(work[0]),m,n,lwork,"sormqr");
#endif
        }
    }
    template <int C>
    void DoPackedQ_LDivEq(
        const ConstMatrixView<std::complex<float>,C>& Q,
        const ConstVectorView<float>& beta,
        MatrixView<std::complex<float> > x)
    {
        int k = Q.rowsize();
        Vector<std::complex<float> > cbeta = beta;
        if (Q.isrm()) {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepi();
            if (x.isrm() == Q.isconj()) x.conjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            work.get()[0] = 0.;
            LAPNAME(cunmlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(cunmlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() == Q.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("cunmlq");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"cunmlq");
#endif
        } else {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepj();
            if (x.isrm() != Q.isconj()) x.conjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            work.get()[0] = 0.;
            LAPNAME(cunmqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(cunmqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() != Q.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("cunmqr");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"cunmqr");
#endif
        }
    }
    void DoPackedQ_MultEq(
        const ConstMatrixView<float>& Q,
        const ConstVectorView<float>& beta,
        MatrixView<float> x)
    {
        int ldx = x.isrm() ? x.stepi() : x.stepj();
        int m = x.isrm() ? x.rowsize() : x.colsize();
        int n = x.isrm() ? x.colsize() : x.rowsize();
        int k = Q.rowsize();
        if (Q.isrm()) {
            int ldq = Q.stepi();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.;
            LAPNAME(sormlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_T,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(sormlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_T,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("sormlq");
#else
            LAP_Results(int(work[0]),m,n,lwork,"sormlq");
#endif
        } else {
            int ldq = Q.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.;
            LAPNAME(sormqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
                x.isrm()?LAPCH_T:LAPCH_NT, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(sormqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
                x.isrm()?LAPCH_T:LAPCH_NT, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("sormqr");
#else
            LAP_Results(int(work[0]),m,n,lwork,"sormqr");
#endif
        }
    }
    template <int C>
    void DoPackedQ_MultEq(
        const ConstMatrixView<std::complex<float>,C>& Q,
        const ConstVectorView<float>& beta,
        MatrixView<std::complex<float> > x)
    {
        int k = Q.rowsize();
        Vector<std::complex<float> > cbeta = beta;
        if (Q.isrm()) {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepi();
            if (x.isrm() == Q.isconj()) x.conjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            work.get()[0] = 0.;
            LAPNAME(cunmlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(cunmlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() == Q.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("cunmlq");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"cunmlq");
#endif
        } else {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepj();
            if (x.isrm() != Q.isconj()) x.conjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            work.get()[0] = 0.;
            LAPNAME(cunmqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(cunmqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(cbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() != Q.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("cunmqr");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"cunmqr");
#endif
        }
    }
    static void DoPackedQ_LDivEq(
        const ConstMatrixView<float>& Q,
        const ConstVectorView<float>& beta,
        VectorView<float> x)
    { DoPackedQ_LDivEq(Q,beta,ColVectorViewOf(x).xView()); }
    template <int C>
    void DoPackedQ_LDivEq(
        const ConstMatrixView<std::complex<float>,C>& Q,
        const ConstVectorView<float>& beta,
        VectorView<std::complex<float> > x)
    { DoPackedQ_LDivEq(Q,beta,ColVectorViewOf(x).xView()); }
    static void DoPackedQ_MultEq(
        const ConstMatrixView<float>& Q,
        const ConstVectorView<float>& beta,
        VectorView<float> x)
    { DoPackedQ_MultEq(Q,beta,ColVectorViewOf(x).xView()); }
    template <int C>
    void DoPackedQ_MultEq(
        const ConstMatrixView<std::complex<float>,C>& Q,
        const ConstVectorView<float>& beta,
        VectorView<std::complex<float> > x)
    { DoPackedQ_MultEq(Q,beta,ColVectorViewOf(x).xView()); }
#endif
#endif // LAP

    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_MultEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        MatrixView<T2> m2)
    {
#ifdef XDEBUG_QR
        std::cout<<"InstPackedQ_MultEq matrix\n";
        std::cout<<"Q = "<<TMV_Text(Q)<<"  "<<Q<<std::endl;
        std::cout<<"beta = "<<TMV_Text(beta)<<"  "<<beta<<std::endl;
        std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
        Matrix<T2> m2x = m2;
        InlinePackedQ_MultEq(Q,beta,m2x);
#endif
        if (Q.isrm() || Q.iscm()) {
            if (beta.step() == 1) {
                if (m2.iscm() || m2.isrm()) {
                    DoPackedQ_MultEq(Q,beta,m2); 
                } else {
                    Matrix<T2,ColMajor|NoDivider> m2c = m2;
                    InstPackedQ_MultEq(Q,beta,m2c.xView());
                    InstCopy(m2c.constView().xView(),m2);
                }
            } else {
                Vector<RT1> betac = beta;
                InstPackedQ_MultEq(Q,betac.constView().xView(),m2);
            }
        } else {
            Matrix<T1,ColMajor|NoDivider> Qc = Q;
            InstPackedQ_MultEq(Qc.constView().xView(),beta,m2);
        }
#ifdef XDEBUG_QR
        std::cout<<"m2 => "<<m2<<std::endl;
        std::cout<<"Norm(m2-m2x) = "<<Norm(m2-m2x)<<std::endl;
        if (!(Norm(m2-m2x) < 1.e-3*Norm(Q)*Norm(m2x))) {
            abort();
        }
#endif
    }

    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_LDivEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        MatrixView<T2> m2)
    {
#ifdef XDEBUG_QR
        std::cout<<"InstPackedQ_LDivEq matrix\n";
        std::cout<<"Q = "<<TMV_Text(Q)<<"  "<<Q<<std::endl;
        std::cout<<"beta = "<<TMV_Text(beta)<<"  "<<beta<<std::endl;
        std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
        Matrix<T2> m2x = m2;
        InlinePackedQ_LDivEq(Q,beta,m2x);
#endif
        if (Q.isrm() || Q.iscm()) {
            if (beta.step() == 1) {
                if (m2.iscm() || m2.isrm()) {
                    DoPackedQ_LDivEq(Q,beta,m2); 
                } else {
                    Matrix<T2,ColMajor|NoDivider> m2c = m2;
                    InstPackedQ_LDivEq(Q,beta,m2c.xView());
                    InstCopy(m2c.constView().xView(),m2);
                }
            } else {
                Vector<RT1> betac = beta;
                InstPackedQ_LDivEq(Q,betac.constView().xView(),m2);
            }
        } else {
            Matrix<T1,ColMajor|NoDivider> Qc = Q;
            InstPackedQ_LDivEq(Qc.constView().xView(),beta,m2);
        }
#ifdef XDEBUG_QR
        std::cout<<"m2 => "<<m2<<std::endl;
        std::cout<<"Norm(m2-m2x) = "<<Norm(m2-m2x)<<std::endl;
        if (!(Norm(m2-m2x) < 1.e-3*Norm(Q)*Norm(m2x))) {
            abort();
        }
#endif
    }
 
    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_MultEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        VectorView<T2> v2)
    {
#ifdef XDEBUG_QR
        std::cout<<"InstPackedQ_MultEq vector\n";
        std::cout<<"Q = "<<TMV_Text(Q)<<"  "<<Q<<std::endl;
        std::cout<<"beta = "<<TMV_Text(beta)<<"  "<<beta<<std::endl;
        std::cout<<"v2 = "<<TMV_Text(v2)<<"  "<<v2<<std::endl;
        Vector<T2> v2x = v2;
        InlinePackedQ_MultEq(Q,beta,v2x);
#endif
        if (Q.isrm() || Q.iscm()) {
            if (beta.step() == 1) {
                if (v2.step() == 1) {
                    DoPackedQ_MultEq(Q,beta,v2); 
                } else {
                    Vector<T2> v2c = v2;
                    InstPackedQ_MultEq(Q,beta,v2c.xView());
                    InstCopy(v2c.constView().xView(),v2);
                }
            } else {
                Vector<RT1> betac = beta;
                InstPackedQ_MultEq(Q,betac.constView().xView(),v2);
            }
        } else {
            Matrix<T1,ColMajor|NoDivider> Qc = Q;
            InstPackedQ_MultEq(Qc.constView().xView(),beta,v2);
        }
#ifdef XDEBUG_QR
        std::cout<<"v2 => "<<v2<<std::endl;
        std::cout<<"Norm(v2-v2x) = "<<Norm(v2-v2x)<<std::endl;
        if (!(Norm(v2-v2x) < 1.e-3*Norm(Q)*Norm(v2x))) {
            abort();
        }
#endif
    }

    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_LDivEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        VectorView<T2> v2)
    {
#ifdef XDEBUG_QR
        std::cout<<"InstPackedQ_LDivEq vector\n";
        std::cout<<"Q = "<<TMV_Text(Q)<<"  "<<Q<<std::endl;
        std::cout<<"beta = "<<TMV_Text(beta)<<"  "<<beta<<std::endl;
        std::cout<<"v2 = "<<TMV_Text(v2)<<"  "<<v2<<std::endl;
        Vector<T2> v2x = v2;
        InlinePackedQ_LDivEq(Q,beta,v2x);
#endif
        if (Q.isrm() || Q.iscm()) {
            if (beta.step() == 1) {
                if (v2.step() == 1) {
                    DoPackedQ_LDivEq(Q,beta,v2); 
                } else {
                    Vector<T2> v2c = v2;
                    InstPackedQ_LDivEq(Q,beta,v2c.xView());
                    InstCopy(v2c.constView().xView(),v2);
                }
            } else {
                Vector<RT1> betac = beta;
                InstPackedQ_LDivEq(Q,betac.constView().xView(),v2);
            }
        } else {
            Matrix<T1,ColMajor|NoDivider> Qc = Q;
            InstPackedQ_LDivEq(Qc.constView().xView(),beta,v2);
        }
#ifdef XDEBUG_QR
        std::cout<<"v2 => "<<v2<<std::endl;
        std::cout<<"Norm(v2-v2x) = "<<Norm(v2-v2x)<<std::endl;
        if (!(Norm(v2-v2x) < 1.e-3*Norm(Q)*Norm(v2x))) {
            abort();
        }
#endif
    }


#define InstFile "TMV_PackedQ.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


