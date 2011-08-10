
//#define PRINTALGO_MV

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

    template <class T, class M1, class V2>
    static void DoMultMV(
        const T x, const M1& m1, const V2& v2, VectorView<T>& v3)
    {
        // Check for non-unit step and x != 1, and do the necessary copies here,
        // rather than in the InlineMultMV function.  
        // This is faster to compile, since it keeps the InlineMultMV
        // algo path to the ones that have vstep == 1.

        typedef typename Traits<T>::real_type RT;
        const Scaling<1,RT> one;

        const int M = m1.colsize();
        const int N = m1.rowsize();

        if (v3.step() == 1) {
            VectorView<T,Unit> v3u = v3.unitView();
            if (v2.step() == 1) {
                if (x == RT(1))
                    InlineMultMV<false>(one,m1,v2.unitView(),v3u);
                else if (M > 4*N) {
                    Vector<T> xv2(N);
                    InstMultXV(x,v2,xv2.xView());
                    InlineMultMV<false>(one,m1,xv2,v3u);
                } else {
                    InlineMultMV<false>(one,m1,v2.unitView(),v3u);
                    InstScale(x,v3);
                }
            } else {
                Vector<T> xv2(N);
                InstMultXV(x,v2,xv2.xView());
                InlineMultMV<false>(one,m1,xv2,v3u);
            }
        } else {
            Vector<T> v3c(M);
            VectorView<T,Unit> v3u = v3c.unitView();
            if (v2.step() == 1) {
                InlineMultMV<false>(one,m1,v2.unitView(),v3u);
                InstMultXV(x,v3c.constView().xView(),v3);
            } else if (M > N) {
                Vector<T> xv2(N);
                InstMultXV(x,v2,xv2.xView());
                InlineMultMV<false>(one,m1,xv2,v3u);
                InstCopy(v3c.constView().xView(),v3);
            } else {
                InlineMultMV<false>(one,m1,v2.copy(),v3u);
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
        const Scaling<1,RT> one;

        const int M = m1.colsize();
        const int N = m1.rowsize();

        if (v3.step() == 1) {
            VectorView<T,Unit> v3u = v3.unitView();
            if (v2.step() == 1) {
                if (x == RT(1))
                    InlineMultMV<true>(one,m1,v2.unitView(),v3u);
                else if (N > 4*M) {
                    Vector<T> v3c(M);
                    InlineMultMV<false>(one,m1,v2.unitView(),v3c);
                    InstAddMultXV(x,v3c.constView().xView(),v3);
                } else {
                    Vector<T> xv2(N);
                    InstMultXV(x,v2,xv2.xView());
                    InlineMultMV<true>(one,m1,xv2,v3u);
                }
            } else {
                Vector<T> xv2(N);
                InstMultXV(x,v2,xv2.xView());
                InlineMultMV<true>(one,m1,xv2,v3u);
            }
        } else {
            Vector<T> v3c(v3.size(),T(0));
            VectorView<T,Unit> v3u = v3c.unitView();
            if (v2.step() == 1) {
                InlineMultMV<true>(one,m1,v2.unitView(),v3u);
                InstAddMultXV(x,v3c.constView().xView(),v3);
            } else if (M > N) {
                Vector<T> xv2(N);
                InstMultXV(x,v2,xv2.xView());
                InlineMultMV<true>(one,m1,xv2,v3u);
                InstAddMultXV(T(1),v3c.constView().xView(),v3);
            } else {
                InlineMultMV<true>(one,m1,v2.copy(),v3u);
                InstAddMultXV(x,v3c.constView().xView(),v3);
            }
        }
    }

#ifdef BLAS
    template <bool add, int dummy>
    struct BlasMultMV_Helper;
    template <int dummy>
    struct BlasMultMV_Helper<false,dummy>
    {
        template <class T, class M1, class V2>
        static void call(
            const T x, const M1& m1, const V2& v2, VectorView<T>& v3)
        { DoMultMV(x,m1,v2,v3); }
    };
    template <int dummy>
    struct BlasMultMV_Helper<true,dummy>
    {
        template <class T, class M1, class V2>
        static void call(
            const T x, const M1& m1, const V2& v2, VectorView<T>& v3)
        { DoAddMultMV(x,m1,v2,v3); }
    };
    template <bool add, class T, class M1, class V2>
    static void BlasMultMV(
        const T x, const M1& m1, const V2& v2, VectorView<T>& v3)
    { BlasMultMV_Helper<add,1>::call(x,m1,v2,v3); }
#ifdef TMV_INST_DOUBLE
    static void BlasMultMV(
        double alpha,
        const ConstMatrixView<double>& A,
        const ConstVectorView<double>& x,
        double beta, VectorView<double> y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(!SameStorage(A,y));

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
        BLASNAME(dgemv) (
            BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASV(beta),BLASP(yp),BLASV(ys)
            BLAS1);
    }
    template <int C1, int C2>
    static void BlasMultMV(
        std::complex<double> alpha,
        const ConstMatrixView<std::complex<double>,C1>& A,
        const ConstVectorView<std::complex<double>,C2>& x,
        double beta, VectorView<std::complex<double> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(!SameStorage(A,y));

        if (x.isconj()
#ifndef CBLAS
            && !(A.isconj() && A.iscm()) 
#endif
        ) {
            const Vector<std::complex<double> > xx = alpha*x;
            return BlasMultMV(std::complex<double>(1.),A,xx.xView(),beta,y);
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
        std::complex<double> xbeta = beta;
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
                BLASCM A.isrm()?A.isconj()?BLASCH_CT:BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                BLAS1);
        }
    }
#ifdef TMV_INST_MIX
    template <int C1>
    static void BlasMultMV(
        std::complex<double> alpha,
        const ConstMatrixView<std::complex<double>,C1>& A,
        const ConstVectorView<double>& x,
        double beta, VectorView<std::complex<double> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));

        if (A.iscm()) {
            if (y.step() != 1) {
                Vector<std::complex<double> > yy(y.size());
                BlasMultMV(1.,A,x,0.,yy.xView());
                if (beta == 0.) y = alpha*yy;
                else y += alpha*yy;
            } else {
                if (beta == 0.) {
                    int m = 2*A.colsize();
                    int n = A.rowsize();
                    int lda = 2*A.stepj();
                    if (lda < m) { TMVAssert(n==1); lda = m; }
                    int xs = x.step();
                    int ys = 1;
                    const double* xp = x.cptr();
                    if (xs < 0) xp += (x.size()-1)*xs;
                    double* yp = (double*) y.ptr();
                    double xalpha(1);
                    BLASNAME(dgemv) (
                        BLASCM BLASCH_NT,
                        BLASV(m),BLASV(n),BLASV(xalpha),
                        BLASP((double*)A.cptr()),BLASV(lda),
                        BLASP(xp),BLASV(xs),BLASV(beta),
                        BLASP(yp),BLASV(ys) BLAS1);
                    if (A.isconj()) y.conjugateSelf();
                    y *= alpha;
                } else if (A.isconj()) {
                    Vector<std::complex<double> > yy(y.size());
                    BlasMultMV(1.,A.conjugate(),x,0.,yy.xView());
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
                    Vector<std::complex<double> > yy(y.size());
                    BlasMultMV(std::complex<double>(1.),A,x,0.,yy.xView());
                    y += alpha*yy;
                }
            }
        } else { // A.isrm
            const Vector<std::complex<double> > xx = x;
            BlasMultMV(alpha,A,xx.xView(),beta,y);
        }
    }
    template <int C2>
    static void BlasMultMV(  
        std::complex<double> alpha,
        const ConstMatrixView<double>& A,
        const ConstVectorView<std::complex<double>,C2>& x,
        double beta, VectorView<std::complex<double> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));

        if (beta == 0.) {
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
            double xalpha(1);
            BLASNAME(dgemv) (
                BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(dgemv) (
                BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp+1),BLASV(xs),BLASV(beta),
                BLASP(yp+1),BLASV(ys) BLAS1);
            if (x.isconj()) y.conjugateSelf();
            y *= alpha;
        } else if (TMV_IMAG(alpha) == 0. && !x.isconj()) {
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
                BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(dgemv) (
                BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp+1),BLASV(xs),BLASV(beta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        } else {
            const Vector<std::complex<double> > xx = alpha*x;
            BlasMultMV(std::complex<double>(1.),A,xx.xView(),1.,y);
        }
    }
    static void BlasMultMV(
        std::complex<double> alpha,
        const ConstMatrixView<double>& A,
        const ConstVectorView<double>& x,
        double beta, VectorView<std::complex<double> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));

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
        if (ar == 0.) {
            if (beta == 0.) y.realPart().setZero();
        } else {
            BLASNAME(dgemv) (
                BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASV(ar),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp),BLASV(ys) BLAS1);
        }
        if (ai == 0.) {
            if (beta == 0.) y.imagPart().setZero();
        } else {
            BLASNAME(dgemv) (
                BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASV(ai),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        }
    }
#endif
    template <bool add, class M1, class V2>
    static void BlasMultMV(
        double alpha, const M1& A, const V2& x, VectorView<double> y)
    { BlasMultMV(alpha,A,x,add?1:0,y); }
    template <bool add, class M1, class V2>
    static void BlasMultMV(
        std::complex<double> alpha, const M1& A, const V2& x,
        VectorView<std::complex<double> > y)
    { BlasMultMV(alpha,A,x,add?1:0,y); }
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
    static void BlasMultMV(
        float alpha,
        const ConstMatrixView<float>& A,
        const ConstVectorView<float>& x,
        float beta, VectorView<float> y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(!SameStorage(A,y));

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
        BLASNAME(sgemv) (
            BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASV(beta),BLASP(yp),BLASV(ys)
            BLAS1);
    }
    template <int C1, int C2>
    static void BlasMultMV(
        std::complex<float> alpha,
        const ConstMatrixView<std::complex<float>,C1>& A,
        const ConstVectorView<std::complex<float>,C2>& x,
        float beta, VectorView<std::complex<float> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(!SameStorage(A,y));

        if (x.isconj()
#ifndef CBLAS
            && !(A.isconj() && A.iscm()) 
#endif
        ) {
            const Vector<std::complex<float> > xx = alpha*x;
            return BlasMultMV(std::complex<float>(1.F),A,xx.xView(),beta,y);
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
        std::complex<float> xbeta = beta;
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
                BLASCM A.isrm()?A.isconj()?BLASCH_CT:BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
                BLAS1);
        }
    }
#ifdef TMV_INST_MIX
    template <int C1>
    static void BlasMultMV(
        std::complex<float> alpha,
        const ConstMatrixView<std::complex<float>,C1>& A,
        const ConstVectorView<float>& x,
        float beta, VectorView<std::complex<float> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));

        if (A.iscm()) {
            if (y.step() != 1) {
                Vector<std::complex<float> > yy(y.size());
                BlasMultMV(1.F,A,x,0.F,yy.xView());
                if (beta == 0.F) y = alpha*yy;
                else y += alpha*yy;
            } else {
                if (beta == 0.F) {
                    int m = 2*A.colsize();
                    int n = A.rowsize();
                    int lda = 2*A.stepj();
                    if (lda < m) { TMVAssert(n==1); lda = m; }
                    int xs = x.step();
                    int ys = 1;
                    const float* xp = x.cptr();
                    if (xs < 0) xp += (x.size()-1)*xs;
                    float* yp = (float*) y.ptr();
                    float xalpha(1);
                    BLASNAME(sgemv) (
                        BLASCM BLASCH_NT,
                        BLASV(m),BLASV(n),BLASV(xalpha),
                        BLASP((float*)A.cptr()),BLASV(lda),
                        BLASP(xp),BLASV(xs),BLASV(beta),
                        BLASP(yp),BLASV(ys) BLAS1);
                    if (A.isconj()) y.conjugateSelf();
                    y *= alpha;
                } else if (A.isconj()) {
                    Vector<std::complex<float> > yy(y.size());
                    BlasMultMV(1.F,A.conjugate(),x,0.F,yy.xView());
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
                    BLASNAME(sgemv) (
                        BLASCM BLASCH_NT,
                        BLASV(m),BLASV(n),BLASV(xalpha),
                        BLASP((float*)A.cptr()),BLASV(lda),
                        BLASP(xp),BLASV(xs),BLASV(beta),
                        BLASP(yp),BLASV(ys) BLAS1);
                } else {
                    Vector<std::complex<float> > yy(y.size());
                    BlasMultMV(std::complex<float>(1.F),A,x,0.F,yy.xView());
                    y += alpha*yy;
                }
            }
        } else { // A.isrm
            const Vector<std::complex<float> > xx = x;
            BlasMultMV(alpha,A,xx.xView(),beta,y);
        }
    }
    template <int C2>
    static void BlasMultMV(  
        std::complex<float> alpha,
        const ConstMatrixView<float>& A,
        const ConstVectorView<std::complex<float>,C2>& x,
        float beta, VectorView<std::complex<float> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));

        if (beta == 0.F) {
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
            float xalpha(1);
            BLASNAME(sgemv) (
                BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(sgemv) (
                BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp+1),BLASV(xs),BLASV(beta),
                BLASP(yp+1),BLASV(ys) BLAS1);
            if (x.isconj()) y.conjugateSelf();
            y *= alpha;
        } else if (TMV_IMAG(alpha) == 0.F && !x.isconj()) {
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
            BLASNAME(sgemv) (
                BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(sgemv) (
                BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp+1),BLASV(xs),BLASV(beta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        } else {
            const Vector<std::complex<float> > xx = alpha*x;
            BlasMultMV(std::complex<float>(1.F),A,xx.xView(),1.F,y);
        }
    }
    static void BlasMultMV(
        std::complex<float> alpha,
        const ConstMatrixView<float>& A,
        const ConstVectorView<float>& x,
        float beta, VectorView<std::complex<float> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));

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
        if (ar == 0.F) {
            if (beta == 0.F) y.realPart().setZero();
        } else {
            BLASNAME(sgemv) (
                BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASV(ar),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp),BLASV(ys) BLAS1);
        }
        if (ai == 0.F) {
            if (beta == 0.F) y.imagPart().setZero();
        } else {
            BLASNAME(sgemv) (
                BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
                BLASV(m),BLASV(n),BLASV(ai),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs),BLASV(beta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        }
    }
#endif
#endif // FLOAT
    template <bool add, class M1, class V2>
    static void BlasMultMV(
        float alpha, const M1& A, const V2& x, VectorView<float> y)
    { BlasMultMV(alpha,A,x,add?1:0,y); }
    template <bool add, class M1, class V2>
    static void BlasMultMV(
        std::complex<float> alpha, const M1& A, const V2& x,
        VectorView<std::complex<float> > y)
    { BlasMultMV(alpha,A,x,add?1:0,y); }
#endif // BLAS

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMV(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    {
#ifdef BLAS
        BlasMultMV<false>(x,m1,v2,v3);
#elif TMV_OPT <= 2
        v3.setZero();
        InstAddMultMV(x,m1,v2,v3);
#else
        if (m1.isrm())
            DoMultMV(x,m1.rmView(),v2,v3);
        else if (m1.iscm())
            DoMultMV(x,m1.cmView(),v2,v3);
        else {
            Matrix<T1,ColMajor|NoDivider> m1c = m1;
            DoMultMV(x,m1c.constView(),v2,v3);
        }
#endif
    }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMV(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    {
#ifdef BLAS
        BlasMultMV<true>(x,m1,v2,v3);
#else
        if (m1.isrm())
            DoAddMultMV(x,m1.rmView(),v2,v3);
        else if (m1.iscm())
            DoAddMultMV(x,m1.cmView(),v2,v3);
        else {
            Matrix<T1,ColMajor|NoDivider> m1c = m1;
            DoAddMultMV(x,m1c.constView(),v2,v3);
        }
#endif
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


