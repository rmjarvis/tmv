
//#define PRINTALGO_MV
//#define XDEBUG_MV

#include "TMV_Blas.h"
#include "tmv/TMV_MultMV.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_MultXV.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_ProdXV.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_ProdMV.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_ConjugateV.h"

#ifdef BLAS
#include "tmv/TMV_AddVV.h"
#include "tmv/TMV_SumVV.h"
#include "tmv/TMV_ProdXV.h"
#include "tmv/TMV_ProdXM.h"
#endif

// The CBLAS trick of using RowMajor with ConjTrans when we have a 
// case of A.conjugate() * x doesn't seem to be working with MKL 10.2.2.
// I haven't been able to figure out why.  (e.g. Is it a bug in the MKL
// code, or am I doing something wrong?)  So for now, just disable it.
#ifdef CBLAS
#undef CBLAS
#endif

namespace tmv {

    template <bool add, class T, class M1, class V2>
    static void CallInlineMultMV(
        const M1& m1, const V2& v2, VectorView<T,Unit> v3)
    {
        TMVAssert(v2.step() == 1);
        typedef typename Traits<T>::real_type RT;
        const Scaling<1,RT> one;
        if (m1.iscm())
            InlineMultMV<add>(one,m1.cmView(),v2.unitView(),v3);
        else if (m1.isrm())
            InlineMultMV<add>(one,m1.rmView(),v2.unitView(),v3);
        else
            InlineMultMV<add>(one,m1,v2.unitView(),v3);
    }
    
    template <class T, class M1, class V2>
    static void DoMultMV(
        const T x, const M1& m1, const V2& v2, VectorView<T>& v3)
    {
        // Check for non-unit step and x != 1, and do the necessary copies here,
        // rather than in the InlineMultMV function.  
        // This is faster to compile, since it keeps the InlineMultMV
        // algo path to the ones that have vstep == 1.

        typedef typename Traits<T>::real_type RT;
        const ptrdiff_t M = m1.colsize();
        const ptrdiff_t N = m1.rowsize();

        if (v3.step() == 1) {
            if (v2.step() == 1) {
                if (x == RT(1))
                    CallInlineMultMV<false>(m1,v2,v3.unitView());
                else if (M > 4*N) {
                    Vector<T> xv2(N);
                    InstMultXV(x,v2,xv2.xView());
                    CallInlineMultMV<false>(m1,xv2,v3.unitView());
                } else {
                    CallInlineMultMV<false>(m1,v2,v3.unitView());
                    InstScale(x,v3);
                }
            } else {
                Vector<T> xv2(N);
                InstMultXV(x,v2,xv2.xView());
                CallInlineMultMV<false>(m1,xv2,v3.unitView());
            }
        } else {
            Vector<T> v3c(M);
            if (v2.step() == 1) {
                CallInlineMultMV<false>(m1,v2,v3c.unitView());
                InstMultXV(x,v3c.constView().xView(),v3);
            } else if (M > N) {
                Vector<T> xv2(N);
                InstMultXV(x,v2,xv2.xView());
                CallInlineMultMV<false>(m1,xv2,v3c.unitView());
                InstCopy(v3c.constView().xView(),v3);
            } else {
                CallInlineMultMV<false>(m1,v2.copy(),v3c.unitView());
                InstMultXV(x,v3c.constView().xView(),v3);
            }
        }
    }

    template <class T, class M1, class V2>
    static void DoAddMultMV(
        const T x, const M1& m1, const V2& v2, VectorView<T>& v3)
    {
        // Check for non-unit step and x != 1, and do the necessary copies here,
        // rather than in the InlineMultMV function.  
        // This is faster to compile, since it keeps the InlineMultMV
        // algo path to the ones that have vstep == 1.

        typedef typename Traits<T>::real_type RT;
        const ptrdiff_t M = m1.colsize();
        const ptrdiff_t N = m1.rowsize();

        if (v3.step() == 1) {
            if (v2.step() == 1) {
                if (x == RT(1))
                    CallInlineMultMV<true>(m1,v2,v3.unitView());
                else if (N > 4*M) {
                    Vector<T> v3c(M);
                    CallInlineMultMV<false>(m1,v2,v3c.unitView());
                    InstAddMultXV(x,v3c.constView().xView(),v3);
                } else {
                    Vector<T> xv2(N);
                    InstMultXV(x,v2,xv2.xView());
                    CallInlineMultMV<true>(m1,xv2,v3.unitView());
                }
            } else {
                Vector<T> xv2(N);
                InstMultXV(x,v2,xv2.xView());
                CallInlineMultMV<true>(m1,xv2,v3.unitView());
            }
        } else {
#if TMV_OPT > 2
            Vector<T> v3c(v3.size());
            const bool add2 = false;
#else
            Vector<T> v3c(v3.size(),T(0));
            const bool add2 = true;
#endif
            if (v2.step() == 1) {
                CallInlineMultMV<add2>(m1,v2,v3c.unitView());
                InstAddMultXV(x,v3c.constView().xView(),v3);
            } else if (M > N) {
                Vector<T> xv2(N);
                InstMultXV(x,v2,xv2.xView());
                CallInlineMultMV<add2>(m1,xv2,v3c.unitView());
                InstAddMultXV(T(1),v3c.constView().xView(),v3);
            } else {
                CallInlineMultMV<add2>(m1,v2.copy(),v3c.unitView());
                InstAddMultXV(x,v3c.constView().xView(),v3);
            }
        }
    }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
    static void DoAddMultMV(
        double alpha,
        const ConstMatrixView<double>& A,
        const ConstVectorView<double>& x,
        VectorView<double> y)
    {
        //std::cout<<"dddd\n";
        int m = A.iscm() ? A.colsize() : A.rowsize();
        int n = A.iscm() ? A.rowsize() : A.colsize();
        int lda = A.iscm() ? A.stepj() : A.stepi();
        if (lda < m) { TMVAssert(n==1); lda = m; }
        int xs = x.step();
        int ys = y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const double* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        double* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        double beta = 1.;
        //std::cout<<"m,n = "<<m<<','<<n<<std::endl;
        //std::cout<<"alpha,Ap,lda = "<<alpha<<','<<A.cptr()<<','<<lda<<std::endl;
        //std::cout<<"xp,xs = "<<xp<<','<<xs<<std::endl;
        //std::cout<<"beta,yp,ys = "<<beta<<','<<yp<<','<<ys<<std::endl;
        BLASNAME(dgemv) (
            BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASV(beta),BLASP(yp),BLASV(ys)
            BLAS1);
    }
    template <int C1, int C2>
    static void DoAddMultMV(
        std::complex<double> alpha,
        const ConstMatrixView<std::complex<double>,C1>& A,
        const ConstVectorView<std::complex<double>,C2>& x,
        VectorView<std::complex<double> > y)
    {
        //std::cout<<"zzzz\n";
        if (x.isconj()
#ifndef CBLAS
            && !(A.isconj() && A.iscm()) 
#endif
        ) {
            const Vector<std::complex<double> > xx = alpha*x;
            return DoAddMultMV(std::complex<double>(1.),A,xx.xView(),y);
        } 

        int m = A.iscm() ? A.colsize() : A.rowsize();
        int n = A.iscm() ? A.rowsize() : A.colsize();
        int lda = A.iscm() ? A.stepj() : A.stepi();
        if (lda < m) { TMVAssert(n==1); lda = m; }
        int xs = x.step();
        int ys = y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const std::complex<double>* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        std::complex<double>* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        std::complex<double> xbeta = 1.;
        if (A.isconj() && A.iscm()) {
#ifdef CBLAS
            TMV_SWAP(m,n);
            BLASNAME(zgemv) (
                BLASRM BLASCH_CT,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                BLAS1);
#else
            std::complex<double> ca = TMV_CONJ(alpha);
            if (x.isconj()) {
                y.conjugateSelf();
                BLASNAME(zgemv) (
                    BLASCM BLASCH_NT,
                    BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                    BLAS1);
                y.conjugateSelf();
            } else {
                Vector<std::complex<double> > xx = ca*x.conjugate();
                ca = std::complex<double>(1.);
                xs = 1;
                xp = xx.cptr();
                y.conjugateSelf();
                BLASNAME(zgemv) (
                    BLASCM BLASCH_NT,
                    BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                    BLAS1);
                y.conjugateSelf();
            }
#endif
        } else {
            BLASNAME(zgemv) (
                BLASCM A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                BLAS1);
        }
    }
#ifdef TMV_INST_MIX
    template <int C1>
    static void DoAddMultMV(
        std::complex<double> alpha,
        const ConstMatrixView<std::complex<double>,C1>& A,
        const ConstVectorView<double>& x,
        VectorView<std::complex<double> > y)
    {
        //std::cout<<"zzdz\n";
        double beta = 1.;
        if (A.iscm()) {
            if (y.step() != 1) {
                Vector<std::complex<double> > yy(y.size(),0.);
                DoAddMultMV(1.,A,x,yy.xView());
                y += alpha*yy;
            } else {
                if (A.isconj()) {
                    Vector<std::complex<double> > yy(y.size(),0.);
                    DoAddMultMV(1.,A.conjugate(),x,yy.xView());
                    y += alpha*yy.conjugate();
                } else if (TMV_IMAG(alpha) == 0.) {
                    int m = 2*A.colsize();
                    int n = A.rowsize();
                    int lda = 2*A.stepj();
                    if (lda < m) { TMVAssert(n==1); lda = m; }
                    int xs = x.step();
                    int ys = 1;
                    const double* xp = x.cptr();
                    if (xs < 0) xp += (x.size()-1)*xs;
                    double* yp = (double*) y.ptr();
                    if (ys < 0) yp += (y.size()-1)*ys;
                    double xalpha(TMV_REAL(alpha));
                    BLASNAME(dgemv) (
                        BLASCM BLASCH_NT,
                        BLASV(m),BLASV(n),BLASV(xalpha),
                        BLASP((double*)A.cptr()),BLASV(lda),
                        BLASP(xp),BLASV(xs),BLASV(beta),
                        BLASP(yp),BLASV(ys) BLAS1);
                } else {
                    Vector<std::complex<double> > yy(y.size(),0.);
                    DoAddMultMV(std::complex<double>(1.),A,x,yy.xView());
                    y += alpha*yy;
                }
            }
        } else { // A.isrm
            const Vector<std::complex<double> > xx = x;
            DoAddMultMV(alpha,A,xx.xView(),y);
        }
    }
    template <int C2>
    static void DoAddMultMV(  
        std::complex<double> alpha,
        const ConstMatrixView<double>& A,
        const ConstVectorView<std::complex<double>,C2>& x,
        VectorView<std::complex<double> > y)
    {
        //std::cout<<"zdzz\n";
        double beta = 1.;
        if (TMV_IMAG(alpha) == 0. && !x.isconj()) {
            int m = A.iscm() ? A.colsize() : A.rowsize();
            int n = A.iscm() ? A.rowsize() : A.colsize();
            int lda = A.iscm() ? A.stepj() : A.stepi();
            if (lda < m) { TMVAssert(n==1); lda = m; }
            int xs = 2*x.step();
            int ys = 2*y.step();
            if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
            if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
            const double* xp = (const double*) x.cptr();
            if (xs < 0) xp += (x.size()-1)*xs;
            double* yp = (double*) y.ptr();
            if (ys < 0) yp += (y.size()-1)*ys;
            double xalpha(TMV_REAL(alpha));
            BLASNAME(dgemv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(dgemv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp+1),BLASV(xs),BLASV(beta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        } else {
            const Vector<std::complex<double> > xx = alpha*x;
            DoAddMultMV(std::complex<double>(1.),A,xx.xView(),y);
        }
    }
    static void DoAddMultMV(
        std::complex<double> alpha,
        const ConstMatrixView<double>& A,
        const ConstVectorView<double>& x,
        VectorView<std::complex<double> > y)
    {
        //std::cout<<"zddz\n";
        int m = A.iscm() ? A.colsize() : A.rowsize();
        int n = A.iscm() ? A.rowsize() : A.colsize();
        int lda = A.iscm() ? A.stepj() : A.stepi();
        if (lda < m) { TMVAssert(n==1); lda = m; }
        int xs = x.step();
        int ys = 2*y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const double* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        double* yp = (double*) y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        double ar(TMV_REAL(alpha));
        double ai(TMV_IMAG(alpha));
        double beta = 1.;
        if (ar != 0.) {
            BLASNAME(dgemv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(ar),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp),BLASV(ys) BLAS1);
        }
        if (ai != 0.) {
            BLASNAME(dgemv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(ai),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        }
    }
#endif
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
    static void DoAddMultMV(
        float alpha,
        const ConstMatrixView<float>& A,
        const ConstVectorView<float>& x,
        VectorView<float> y)
    {
        int m = A.iscm() ? A.colsize() : A.rowsize();
        int n = A.iscm() ? A.rowsize() : A.colsize();
        int lda = A.iscm() ? A.stepj() : A.stepi();
        if (lda < m) { TMVAssert(n==1); lda = m; }
        int xs = x.step();
        int ys = y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const float* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        float* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        float beta = 1.F;
        BLASNAME(sgemv) (
            BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASV(beta),BLASP(yp),BLASV(ys)
            BLAS1);
    }
    template <int C1, int C2>
    static void DoAddMultMV(
        std::complex<float> alpha,
        const ConstMatrixView<std::complex<float>,C1>& A,
        const ConstVectorView<std::complex<float>,C2>& x,
        VectorView<std::complex<float> > y)
    {
        if (x.isconj()
#ifndef CBLAS
            && !(A.isconj() && A.iscm()) 
#endif
        ) {
            const Vector<std::complex<float> > xx = alpha*x;
            return DoAddMultMV(std::complex<float>(1.F),A,xx.xView(),y);
        } 

        int m = A.iscm() ? A.colsize() : A.rowsize();
        int n = A.iscm() ? A.rowsize() : A.colsize();
        int lda = A.iscm() ? A.stepj() : A.stepi();
        if (lda < m) { TMVAssert(n==1); lda = m; }
        int xs = x.step();
        int ys = y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const std::complex<float>* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        std::complex<float>* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        std::complex<float> xbeta = 1.F;
        if (A.isconj() && A.iscm()) {
#ifdef CBLAS
            TMV_SWAP(m,n);
            BLASNAME(cgemv) (
                BLASRM BLASCH_CT,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                BLAS1);
#else
            std::complex<float> ca = TMV_CONJ(alpha);
            if (x.isconj()) {
                y.conjugateSelf();
                BLASNAME(cgemv) (
                    BLASCM BLASCH_NT,
                    BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                    BLAS1);
                y.conjugateSelf();
            } else {
                Vector<std::complex<float> > xx = ca*x.conjugate();
                ca = std::complex<float>(1.F);
                xs = 1;
                xp = xx.cptr();
                y.conjugateSelf();
                BLASNAME(cgemv) (
                    BLASCM BLASCH_NT,
                    BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                    BLAS1);
                y.conjugateSelf();
            }
#endif
        } else {
            BLASNAME(cgemv) (
                BLASCM A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                BLAS1);
        }
    }
#ifdef TMV_INST_MIX
    template <int C1>
    static void DoAddMultMV(
        std::complex<float> alpha,
        const ConstMatrixView<std::complex<float>,C1>& A,
        const ConstVectorView<float>& x,
        VectorView<std::complex<float> > y)
    {
        if (A.iscm()) {
            if (y.step() != 1) {
                Vector<std::complex<float> > yy(y.size(),0.F);
                DoAddMultMV(1.F,A,x,yy.xView());
                y += alpha*yy;
            } else {
                if (A.isconj()) {
                    Vector<std::complex<float> > yy(y.size(),0.F);
                    DoAddMultMV(1.F,A.conjugate(),x,yy.xView());
                    y += alpha*yy.conjugate();
                } else if (TMV_IMAG(alpha) == 0.F) {
                    int m = 2*A.colsize();
                    int n = A.rowsize();
                    int lda = 2*A.stepj();
                    if (lda < m) { TMVAssert(n==1); lda = m; }
                    int xs = x.step();
                    int ys = 1;
                    const float* xp = x.cptr();
                    if (xs < 0) xp += (x.size()-1)*xs;
                    float* yp = (float*) y.ptr();
                    if (ys < 0) yp += (y.size()-1)*ys;
                    float xalpha(TMV_REAL(alpha));
                    float beta = 1.;
                    BLASNAME(sgemv) (
                        BLASCM BLASCH_NT,
                        BLASV(m),BLASV(n),BLASV(xalpha),
                        BLASP((float*)A.cptr()),BLASV(lda),
                        BLASP(xp),BLASV(xs),BLASV(beta),
                        BLASP(yp),BLASV(ys) BLAS1);
                } else {
                    Vector<std::complex<float> > yy(y.size(),0.F);
                    DoAddMultMV(std::complex<float>(1.F),A,x,yy.xView());
                    y += alpha*yy;
                }
            }
        } else { // A.isrm
            const Vector<std::complex<float> > xx = x;
            DoAddMultMV(alpha,A,xx.xView(),y);
        }
    }
    template <int C2>
    static void DoAddMultMV(  
        std::complex<float> alpha,
        const ConstMatrixView<float>& A,
        const ConstVectorView<std::complex<float>,C2>& x,
        VectorView<std::complex<float> > y)
    {
        if (TMV_IMAG(alpha) == 0.F && !x.isconj()) {
            int m = A.iscm() ? A.colsize() : A.rowsize();
            int n = A.iscm() ? A.rowsize() : A.colsize();
            int lda = A.iscm() ? A.stepj() : A.stepi();
            if (lda < m) { TMVAssert(n==1); lda = m; }
            int xs = 2*x.step();
            int ys = 2*y.step();
            if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
            if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
            const float* xp = (const float*) x.cptr();
            if (xs < 0) xp += (x.size()-1)*xs;
            float* yp = (float*) y.ptr();
            if (ys < 0) yp += (y.size()-1)*ys;
            float xalpha(TMV_REAL(alpha));
            float beta = 1.;
            BLASNAME(sgemv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(sgemv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp+1),BLASV(xs),BLASV(beta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        } else {
            const Vector<std::complex<float> > xx = alpha*x;
            DoAddMultMV(std::complex<float>(1.F),A,xx.xView(),y);
        }
    }
    static void DoAddMultMV(
        std::complex<float> alpha,
        const ConstMatrixView<float>& A,
        const ConstVectorView<float>& x,
        VectorView<std::complex<float> > y)
    {
        int m = A.iscm() ? A.colsize() : A.rowsize();
        int n = A.iscm() ? A.rowsize() : A.colsize();
        int lda = A.iscm() ? A.stepj() : A.stepi();
        if (lda < m) { TMVAssert(n==1); lda = m; }
        int xs = x.step();
        int ys = 2*y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const float* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        float* yp = (float*) y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        float ar(TMV_REAL(alpha));
        float ai(TMV_IMAG(alpha));
        float beta = 1.F;
        if (ar != 0.F) {
            BLASNAME(sgemv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(ar),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp),BLASV(ys) BLAS1);
        }
        if (ai != 0.F) {
            BLASNAME(sgemv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(ai),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        }
    }
#endif
#endif // FLOAT
#endif // BLAS

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMV(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    {
        //std::cout<<"InstMultMV:\n";
        //std::cout<<"x = "<<x<<std::endl;
        //std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
        //std::cout<<"v2 = "<<TMV_Text(v2)<<"  "<<v2<<std::endl;
        //std::cout<<"v3 = "<<TMV_Text(v3)<<"  "<<v3<<std::endl;
#if TMV_OPT <= 2 || defined(BLAS) 
        // Some BLAS implementations seem to have trouble if the 
        // input y has a nan in it and beta = 0.
        // They propagate the nan into the  output.
        // I guess they strictly interpret y = beta*y + alpha*m*x,
        // so if beta = 0, then beta*nan = nan.
        // This is explicitly wrong behavior according to the standard,
        // but we still need to guard against it, since some implementations
        // do it this way.
        // Anyway, to fix this problem, we always use beta=1, and just
        // zero out y before calling the blas function if beta is 0.
        // In other words, only have Blas for the AddMultMV path.
        v3.setZero();
        InstAddMultMV(x,m1,v2,v3);
#else
        if (v3.size() > 0) {
            if (v2.size() > 0) {
                DoMultMV(x,m1,v2,v3);
            } else {
                v3.setZero();
            }
        }
#endif
        //std::cout<<"InstMultMV: v3 -> "<<v3<<std::endl;
    }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMV(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    {
        //std::cout<<"InstAddMultMV:\n";
        //std::cout<<"x = "<<x<<std::endl;
        //std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
        //std::cout<<"v2 = "<<TMV_Text(v2)<<"  "<<v2<<std::endl;
        //std::cout<<"v3 = "<<TMV_Text(v3)<<"  "<<v3<<std::endl;
        if (v3.size() > 0 && v2.size() > 0) {
#ifdef BLAS
            if ( BlasIsCM(m1) || BlasIsRM(m1) ) 
                DoAddMultMV(x,m1,v2,v3);
            else {
                Matrix<T1,ColMajor|NoDivider> m1c = m1;
                DoAddMultMV(x,m1c.constView(),v2,v3);
            }
#else
            DoAddMultMV(x,m1,v2,v3);
#endif
        }
        //std::cout<<"InstAddMultMV: v3 -> "<<v3<<std::endl;
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMV(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    { InlineAliasMultMV<false>(Scaling<0,T3>(x),m1,v2,v3); }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMV(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    { InlineAliasMultMV<true>(Scaling<0,T3>(x),m1,v2,v3); }

#define InstFile "TMV_MultMV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


