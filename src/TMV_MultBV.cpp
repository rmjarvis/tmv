
//#define PRINTALGO_BV
//#define XDEBUG_BV

#include "TMV_Blas.h"
#include "tmv/TMV_MultBV.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_MultXB.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_ProdMV.h"

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
        else if (m1.isdm())
            InlineMultMV<add>(one,m1.dmView(),v2.unitView(),v3);
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
        const int M = m1.colsize();
        const int N = m1.rowsize();

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
        const int M = m1.colsize();
        const int N = m1.rowsize();

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
        const ConstBandMatrixView<double>& A,
        const ConstVectorView<double>& x,
        VectorView<double> y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(!SameStorage(A,y));
        TMVAssert(A.isrm() || A.iscm());
        if (A.isrm()) TMVAssert(A.stepi() >= A.nlo()+A.nhi());
        else TMVAssert(A.stepj() >= A.nlo()+A.nhi());

        int m = A.iscm() ? A.colsize() : A.rowsize();
        int n = A.iscm() ? A.rowsize() : A.colsize();
        int lo = A.iscm() ? A.nlo() : A.nhi();
        int hi = A.iscm() ? A.nhi() : A.nlo();
        int ds = A.diagstep();
        int xs = x.step();
        int ys = y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const double* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        double* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        double beta = 1.;
        BLASNAME(dgbmv) (
            BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
            BLASV(alpha),BLASP(A.cptr()-hi),BLASV(ds),
            BLASP(xp),BLASV(xs),BLASV(beta),
            BLASP(yp),BLASV(ys) BLAS1);
    }
    template <int C1, int C2>
    static void DoAddMultMV(
        std::complex<double> alpha,
        const ConstBandMatrixView<std::complex<double>,C1>& A,
        const ConstVectorView<std::complex<double>,C2>& x,
        VectorView<std::complex<double> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(!SameStorage(A,y));
        if (A.isrm()) TMVAssert(A.stepi() >= A.nlo()+A.nhi());
        else TMVAssert(A.stepj() >= A.nlo()+A.nhi());

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
        int lo = A.iscm() ? A.nlo() : A.nhi();
        int hi = A.iscm() ? A.nhi() : A.nlo();
        int ds = A.diagstep();
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
            TMV_SWAP(lo,hi);
            BLASNAME(zgbmv) (
                BLASRM BLASCH_CT,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASP(&alpha),BLASP(A.cptr()-lo),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
#else
            std::complex<double> ca = TMV_CONJ(alpha);
            if (x.isconj()) {
                y.conjugateSelf();
                BLASNAME(zgbmv) (
                    BLASCM BLASCH_NT, BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                    BLASP(&ca),BLASP(A.cptr()-hi),BLASV(ds),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),
                    BLASP(yp),BLASV(ys) BLAS1);
                y.conjugateSelf();
            } else {
                Vector<std::complex<double> > xx = ca*x.conjugate();
                ca = std::complex<double>(1.);
                xs = 1;
                xp = xx.cptr();
                y.conjugateSelf();
                BLASNAME(zgbmv) (
                    BLASCM BLASCH_NT, BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                    BLASP(&ca),BLASP(A.cptr()-hi),BLASV(ds),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),
                    BLASP(yp),BLASV(ys) BLAS1);
                y.conjugateSelf();
            }
#endif
        } else {
            BLASNAME(zgbmv) (
                BLASCM A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASP(&alpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
        }
    }
#ifdef TMV_INST_MIX
    template <int C1>
    static void DoAddMultMV(
        std::complex<double> alpha,
        const ConstBandMatrixView<std::complex<double>,C1>& A,
        const ConstVectorView<double>& x,
        VectorView<std::complex<double> > y)
    { DoAddMultMV(alpha,A,Vector<std::complex<double> >(x),beta,y); }
    template <int C2>
    static void DoAddMultMV(  
        std::complex<double> alpha,
        const ConstBandMatrixView<double>& A,
        const ConstVectorView<std::complex<double>,C2>& x,
        VectorView<std::complex<double> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));
        if (A.isrm()) TMVAssert(A.stepi() >= A.nlo()+A.nhi());
        else TMVAssert(A.stepj() >= A.nlo()+A.nhi());

        double beta = 1.;
        if (TMV_IMAG(alpha) == 0. && !x.isconj()) {
            int m = A.iscm() ? A.colsize() : A.rowsize();
            int n = A.iscm() ? A.rowsize() : A.colsize();
            int lo = A.iscm() ? A.nlo() : A.nhi();
            int hi = A.iscm() ? A.nhi() : A.nlo();
            int ds = A.diagstep();
            int xs = 2*x.step();
            int ys = 2*y.step();
            if (xs == 0) { TMVAssert(x.size() == 1); xs = 2; }
            if (ys == 0) { TMVAssert(y.size() == 1); ys = 2; }
            const double* xp = (const double*) x.cptr();
            if (xs < 0) xp += (x.size()-1)*xs;
            double* yp = (double*) y.ptr();
            if (ys < 0) yp += (y.size()-1)*ys;
            double xalpha(TMV_REAL(alpha));
            BLASNAME(dgbmv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(xalpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(dgbmv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(xalpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp+1),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        } else {
            const Vector<std::complex<double> > xx = alpha*x;
            DoAddMultMV(std::complex<double>(1.),A,xx.xView(),y);
        }
    }
    static void DoAddMultMV(
        std::complex<double> alpha,
        const ConstBandMatrixView<double>& A,
        const ConstVectorView<double>& x,
        VectorView<std::complex<double> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));
        if (A.isrm()) TMVAssert(A.stepi() >= A.nlo()+A.nhi());
        else TMVAssert(A.stepj() >= A.nlo()+A.nhi());

        int m = A.iscm() ? A.colsize() : A.rowsize();
        int n = A.iscm() ? A.rowsize() : A.colsize();
        int lo = A.iscm() ? A.nlo() : A.nhi();
        int hi = A.iscm() ? A.nhi() : A.nlo();
        int ds = A.diagstep();
        int xs = x.step();
        int ys = 2*y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 2; }
        const double* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        double* yp = (double*) y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        double ar(TMV_REAL(alpha));
        double ai(TMV_IMAG(alpha));
        double beta = 1.;
        if (ar != 0.) {
            BLASNAME(dgbmv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(ar),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
        }
        if (ai != 0.) {
            BLASNAME(dgbmv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(ai),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        }
    }
#endif
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
    static void DoAddMultMV(
        float alpha,
        const ConstBandMatrixView<float>& A,
        const ConstVectorView<float>& x,
        VectorView<float> y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(!SameStorage(A,y));
        TMVAssert(A.isrm() || A.iscm());
        if (A.isrm()) TMVAssert(A.stepi() >= A.nlo()+A.nhi());
        else TMVAssert(A.stepj() >= A.nlo()+A.nhi());

        int m = A.iscm() ? A.colsize() : A.rowsize();
        int n = A.iscm() ? A.rowsize() : A.colsize();
        int lo = A.iscm() ? A.nlo() : A.nhi();
        int hi = A.iscm() ? A.nhi() : A.nlo();
        int ds = A.diagstep();
        int xs = x.step();
        int ys = y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
        const float* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        float* yp = y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        float beta = 1.F;
        BLASNAME(sgbmv) (
            BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
            BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
            BLASV(alpha),BLASP(A.cptr()-hi),BLASV(ds),
            BLASP(xp),BLASV(xs),BLASV(beta),
            BLASP(yp),BLASV(ys) BLAS1);
    }
    template <int C1, int C2>
    static void DoAddMultMV(
        std::complex<float> alpha,
        const ConstBandMatrixView<std::complex<float>,C1>& A,
        const ConstVectorView<std::complex<float>,C2>& x,
        VectorView<std::complex<float> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));
        TMVAssert(!SameStorage(A,y));
        if (A.isrm()) TMVAssert(A.stepi() >= A.nlo()+A.nhi());
        else TMVAssert(A.stepj() >= A.nlo()+A.nhi());

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
        int lo = A.iscm() ? A.nlo() : A.nhi();
        int hi = A.iscm() ? A.nhi() : A.nlo();
        int ds = A.diagstep();
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
            TMV_SWAP(lo,hi);
            BLASNAME(cgbmv) (
                BLASRM BLASCH_CT,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASP(&alpha),BLASP(A.cptr()-lo),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
#else
            std::complex<float> ca = TMV_CONJ(alpha);
            if (x.isconj()) {
                y.conjugateSelf();
                BLASNAME(cgbmv) (
                    BLASCM BLASCH_NT, BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                    BLASP(&ca),BLASP(A.cptr()-hi),BLASV(ds),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),
                    BLASP(yp),BLASV(ys) BLAS1);
                y.conjugateSelf();
            } else {
                Vector<std::complex<float> > xx = ca*x.conjugate();
                ca = std::complex<float>(1.F);
                xs = 1;
                xp = xx.cptr();
                y.conjugateSelf();
                BLASNAME(cgbmv) (
                    BLASCM BLASCH_NT, BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                    BLASP(&ca),BLASP(A.cptr()-hi),BLASV(ds),
                    BLASP(xp),BLASV(xs),BLASP(&xbeta),
                    BLASP(yp),BLASV(ys) BLAS1);
                y.conjugateSelf();
            }
#endif
        } else {
            BLASNAME(cgbmv) (
                BLASCM A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASP(&alpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASP(&xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
        }
    }
#ifdef TMV_INST_MIX
    template <int C1>
    static void DoAddMultMV(
        std::complex<float> alpha,
        const ConstBandMatrixView<std::complex<float>,C1>& A,
        const ConstVectorView<float>& x,
        VectorView<std::complex<float> > y)
    { DoAddMultMV(alpha,A,Vector<std::complex<float> >(x),beta,y); }
    template <int C2>
    static void DoAddMultMV(  
        std::complex<float> alpha,
        const ConstBandMatrixView<float>& A,
        const ConstVectorView<std::complex<float>,C2>& x,
        VectorView<std::complex<float> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));
        if (A.isrm()) TMVAssert(A.stepi() >= A.nlo()+A.nhi());
        else TMVAssert(A.stepj() >= A.nlo()+A.nhi());

        float beta = 1.F;
        if (TMV_IMAG(alpha) == 0. && !x.isconj()) {
            int m = A.iscm() ? A.colsize() : A.rowsize();
            int n = A.iscm() ? A.rowsize() : A.colsize();
            int lo = A.iscm() ? A.nlo() : A.nhi();
            int hi = A.iscm() ? A.nhi() : A.nlo();
            int ds = A.diagstep();
            int xs = 2*x.step();
            int ys = 2*y.step();
            if (xs == 0) { TMVAssert(x.size() == 1); xs = 2; }
            if (ys == 0) { TMVAssert(y.size() == 1); ys = 2; }
            const float* xp = (const float*) x.cptr();
            if (xs < 0) xp += (x.size()-1)*xs;
            float* yp = (float*) y.ptr();
            if (ys < 0) yp += (y.size()-1)*ys;
            float xalpha(TMV_REAL(alpha));
            BLASNAME(sgbmv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(xalpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
            BLASNAME(sgbmv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(xalpha),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp+1),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        } else {
            const Vector<std::complex<float> > xx = alpha*x;
            DoAddMultMV(std::complex<float>(1.F),A,xx.xView(),y);
        }
    }
    static void DoAddMultMV(
        std::complex<float> alpha,
        const ConstBandMatrixView<float>& A,
        const ConstVectorView<float>& x,
        VectorView<std::complex<float> > y)
    {
        TMVAssert(A.rowsize() == x.size());
        TMVAssert(A.colsize() == y.size());
        TMVAssert(x.size() > 0);
        TMVAssert(y.size() > 0);
        TMVAssert(!SameStorage(x,y));
        if (A.isrm()) TMVAssert(A.stepi() >= A.nlo()+A.nhi());
        else TMVAssert(A.stepj() >= A.nlo()+A.nhi());

        int m = A.iscm() ? A.colsize() : A.rowsize();
        int n = A.iscm() ? A.rowsize() : A.colsize();
        int lo = A.iscm() ? A.nlo() : A.nhi();
        int hi = A.iscm() ? A.nhi() : A.nlo();
        int ds = A.diagstep();
        int xs = x.step();
        int ys = 2*y.step();
        if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
        if (ys == 0) { TMVAssert(y.size() == 1); ys = 2; }
        const float* xp = x.cptr();
        if (xs < 0) xp += (x.size()-1)*xs;
        float* yp = (float*) y.ptr();
        if (ys < 0) yp += (y.size()-1)*ys;
        float ar(TMV_REAL(alpha));
        float ai(TMV_IMAG(alpha));
        float beta = 1.F;
        if (ar != 0.) {
            BLASNAME(sgbmv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(ar),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp),BLASV(ys) BLAS1);
        }
        if (ai != 0.) {
            BLASNAME(sgbmv) (
                BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
                BLASV(m),BLASV(n),BLASV(lo),BLASV(hi),
                BLASV(ai),BLASP(A.cptr()-hi),BLASV(ds),
                BLASP(xp),BLASV(xs),BLASV(xbeta),
                BLASP(yp+1),BLASV(ys) BLAS1);
        }
    }
#endif
#endif // DOUBLE
#endif // BLAS

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMV(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
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
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    {
        //std::cout<<"InstAddMultMV:\n";
        //std::cout<<"x = "<<x<<std::endl;
        //std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
        //std::cout<<"v2 = "<<TMV_Text(v2)<<"  "<<v2<<std::endl;
        //std::cout<<"v3 = "<<TMV_Text(v3)<<"  "<<v3<<std::endl;
        if (v3.size() > 0 && v2.size() > 0) {
#ifdef BLAS
            if ( BlasIsCMBand(m1) || BlasIsRMBand(m1) ) {
                // Check for BandMatrixes that have more diagonals
                // than either colsize or rowsize.  Then stepi or stepj
                // is not large enough for the BLAS function.
                // Need to break up the problem into parts.
                if ((m1.isrm() && m1.stepi() < m1.nlo()+m1.nhi()) ||
                    (m1.iscm() && m1.stepj() < m1.nlo()+m1.nhi())) {
                    if (m1.nlo()+1 == int(m1.colsize())) {
                        if (m1.nhi()+1 == int(m1.rowsize())) {
                            ConstMatrixView<T1,C1> m1a =
                                m1.subMatrix(0,m1.colsize(),0,m1.rowsize());
                            v3 += x * m1a * v2;
                        } else {
                            ConstMatrixView<T1,C1> m1a =
                                m1.subMatrix(0,m1.colsize(),0,m1.nhi());
                            v3 += x * m1a * v2.subVector(0,m1.nhi());
                            ConstBandMatrixView<T1,C1> m1b =
                                m1.colRange(m1.nhi(),m1.rowsize());
                            DoAddMultMV(
                                x,m1b,v2.subVector(m1.nhi(),m1.rowsize()),1,v3);
                        }
                    } else {
                        TMVAssert(m1.nlo()>0);
                        if (m1.nhi()+1 == int(m1.rowsize())) {
                            ConstMatrixView<T1,C1> m1a =
                                m1.subMatrix(0,m1.nlo(),0,m1.rowsize());
                            v3.subVector(0,m1.nlo()) += x * m1a * v2;
                        } else {
                            ConstBandMatrixView<T1,C1> m1a =
                                m1.rowRange(0,m1.nlo());
                            DoAddMultMV(
                                x,m1a,v2.subVector(0,m1a.rowsize()),
                                v3.subVector(0,m1.nlo()));
                        }
                        ConstBandMatrixView<T1,C1> m1b =
                            m1.rowRange(m1.nlo(),m1.colsize());
                        DoAddMultMV(
                            x,m1b,v2,v3.subVector(m1.nlo(),m1.colsize()));
                    }
                } else {
                    DoAddMultMV(x,m1,v2,v3);
                }
            } else {
                BandMatrix<T1,ColMajor|NoDivider> m1c = m1;
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
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    { InlineAliasMultMV<false>(Scaling<0,T3>(x),m1,v2,v3); }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMV(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    { InlineAliasMultMV<true>(Scaling<0,T3>(x),m1,v2,v3); }

#define InstFile "TMV_MultBV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


