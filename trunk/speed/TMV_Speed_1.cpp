
//#undef NDEBUG
//#define PRINTALGO
//#define TMV_INLINE

#include <iostream>
#include "TMV.h"

// How big do you want the vectors to be?
const int N = 1533;
//const int N = 2;

// Define the type to use:
#define TISFLOAT
//#define TISCOMPLEX

// Define whether you want to include the error checks.
#define ERRORCHECK

// Define which versions you want to test:
// Note: often it seems that accurate timings only happen when doing one
// version at a time.  I don't really undrestand why.
#define DOREG
//#define DOSMALL
//#define DOBLAS
//#define DOEIGEN
//#define DOEIGENSMALL

// Define which batches of functions you want to test:
//#define DOMULTXV
//#define DOADDVV
//#define DOMULTVV
#define DOMINMAX
//#define DONORM
//#define DOELEMMULT
//#define DOELEMDIV
//#define DOSWAP

// Set this if you only want to do a single loop.
// Not so useful for timing, but useful for debugging.
//#define ONELOOP

// Set up the target number of operations and memory to use for testing
const int targetnflops = 100000000; // in real ops
const int targetmem = 1000000; // in bytes

// Include the BLAS library you want to test TMV against:
//#include "mkl.h"
extern "C" {
#include "util/fblas.h"
}

// Set these as appropriate given the BLAS library being included above.
//#define CBLAS
#define FBLAS



// The rest of this is based on the above values.
// You shouldn't need to change anything else here.

#include <sys/time.h>
#include <algorithm>
#include <numeric>
#include <iostream>

#if (defined(DOEIGEN) || defined(DOEIGENSMALL))
#include "Eigen/Core"
#include "Eigen/Array"
#endif

#ifdef TISFLOAT
typedef float RT;
#else
typedef double RT;
#endif

#ifdef TISCOMPLEX
typedef std::complex<RT> T;
#define XTWO 2
#else
typedef RT T;
#define XTWO 1
#endif

#ifdef ONELOOP
const int nloops1 = 1;
const int nloops2 = 1;
#else
const int nloops2 = (
    (N * sizeof(T) > targetmem) ? 1 : targetmem / N / (sizeof(T)));
const int nloops1 = targetnflops / (N * nloops2) / XTWO;
#endif

#ifdef DOBLAS

#ifdef TISCOMPLEX
#ifdef TISFLOAT
#define BLAS_LETTER c
#else
#define BLAS_LETTER z
#endif
#else
#ifdef TISFLOAT
#define BLAS_LETTER s
#else
#define BLAS_LETTER d
#endif
#endif

#if defined(FBLAS) && defined(TISCOMPLEX) && defined(TISFLOAT)
#define BT cfloat
#define BP(x) (cfloat*) x
#elif defined(FBLAS) && defined(TISCOMPLEX) && !defined(TISFLOAT)
#define BT cdouble
#define BP(x) (cdouble*) x
#else
#define BP(x) x
#endif

#ifdef CBLAS
#define BLASNAME1(c,x) cblas_ ## c ## x
#define BLASINAME1(c,x) cblas_i ## c ## x
#define BLASM1
#elif defined FBLAS
#define BLASNAME1(c,x) c ## x ## _
#define BLASINAME1(c,x) i ## c ## x ## _
#define BLASM1 -1
#else
UNKNOWN BLAS definition....  // Give compile error
#define BLASNAME1(c,x) c ## x
#define BLASINAME1(c,x) i ## c ## x
#endif
#define BLASNAME2(c,x) BLASNAME1(c,x)
#define BLASINAME2(c,x) BLASINAME1(c,x)

#define BLASNAME(x) BLASNAME2(BLAS_LETTER,x)
#define BLASINAME(x) BLASINAME2(BLAS_LETTER,x)

#ifdef TISCOMPLEX
#ifdef TISFLOAT
#define BLASDZNAME(x) BLASNAME2(sc,x)
#define BLASDNAME(x) BLASNAME2(s,x)
#else
#define BLASDZNAME(x) BLASNAME2(dz,x)
#define BLASDNAME(x) BLASNAME2(d,x)
#endif
#else
#define BLASDZNAME(x) BLASNAME(x)
#define BLASDNAME(x) BLASNAME(x)
#endif

#ifdef TISCOMPLEX
const T xzero(0);
const T xone(1);
const T xseven(7);
const T xeight(8);
const T xoneeighth(1./8.);
const T xmone(-1);
const T xmtwelve(-12);
const T xz89(8,9);
const T xz71(7,1);

const BT* zero = (const BT*) &xzero;
const BT* one = (const BT*) &xone;
const BT* seven = (const BT*) &xseven;
const BT* eight = (const BT*) &xeight;
const BT* oneeighth = (const BT*) &xoneeighth;
const BT* mone = (const BT*) &xmone;
const BT* mtwelve = (const BT*) &xmtwelve;
const BT* z89 = (const BT*) &xz89;
const BT* z71 = (const BT*) &xz71;

const RT dmone(-1);
#else
const T zero(0);
const T one(1);
const T seven(7);
const T eight(8);
const T oneeighth(1./8.);
const T mone(-1);
const T mtwelve(-12);
const T dmone(-1);
#endif

#endif

#ifdef DOEIGEN

#ifdef TISCOMPLEX
#ifdef TISFLOAT
#define EIGENV VectorXcf
#else
#define EIGENV VectorXcd
#endif
#else
#ifdef TISFLOAT
#define EIGENV VectorXf
#else
#define EIGENV VectorXd
#endif
#endif

#endif

#ifdef DOEIGENSMALL
const int nloops2x = (
    (N*sizeof(T) < 1536 * 1024 && N!=Eigen::Dynamic) ? nloops2 : 0 );
#endif

static void ClearCache()
{
    tmv::Vector<double> X(1000000,8.);
    tmv::Vector<double> Y(1000000,8.);
    tmv::Vector<double> Z(1000000);
    Z = X + Y;
    if (Norm(Z) < 5.) exit(1);
}

#ifdef DOMULTXV
static void MultXV(
    const std::vector<tmv::Vector<T> >& A1,
    std::vector<tmv::Vector<T> >& B1)
{
#ifdef ERRORCHECK
    std::vector<tmv::Vector<T> > A0 = A1;
    std::vector<tmv::Vector<T> > B0 = B1;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallVector<T,N> > A2(nloops2);
    std::vector<tmv::SmallVector<T,N> > B2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; B2[k] = B1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Vector<T> > A3 = A1;
    std::vector<tmv::Vector<T> > B3 = B1;
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENV> A4(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENV> B4(nloops2,Eigen::EIGENV(N));
    for(int k=0;k<nloops2;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        const tmv::Vector<T>& B1k = B1[k];
        Eigen::EIGENV& A4k = A4[k];
        Eigen::EIGENV& B4k = B4[k];
        for(int i=0;i<N;++i) A4k(i) = A1k(i);
        for(int i=0;i<N;++i) B4k(i) = B1k(i);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,N,1> > A5;
    std::vector<Eigen::Matrix<T,N,1> > B5;
    if (nloops2x) { A5.resize(nloops2x); B5.resize(nloops2x); }
    for(int k=0;k<nloops2x;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        const tmv::Vector<T>& B1k = B1[k];
        Eigen::Matrix<T,N,1>& A5k = A5[k];
        Eigen::Matrix<T,N,1>& B5k = B5[k];
        for(int i=0;i<N;++i) A5k(i) = A1k(i);
        for(int i=0;i<N;++i) B5k(i) = B1k(i);
    }
#endif

    timeval tp;

    double t1_reg=0., t1_small=0., t1_blas=0., t1_eigen=0., t1_smalleigen=0.;
    double t2_reg=0., t2_small=0., t2_blas=0., t2_eigen=0., t2_smalleigen=0.;
    double t3_reg=0., t3_small=0., t3_blas=0., t3_eigen=0., t3_smalleigen=0.;
    double t4_reg=0., t4_small=0., t4_blas=0., t4_eigen=0., t4_smalleigen=0.;
    double t5_reg=0., t5_small=0., t5_blas=0., t5_eigen=0., t5_smalleigen=0.;
    double t6_reg=0., t6_small=0., t6_blas=0., t6_eigen=0., t6_smalleigen=0.;
    double t7_reg=0., t7_small=0., t7_blas=0., t7_eigen=0., t7_smalleigen=0.;
    double t8_reg=0., t8_small=0., t8_blas=0., t8_eigen=0., t8_smalleigen=0.;
    double ta,tb;

#ifdef ERRORCHECK
    double e1_reg=0., e1_small=0., e1_blas=0., e1_eigen=0., e1_smalleigen=0.;
    double e2_reg=0., e2_small=0., e2_blas=0., e2_eigen=0., e2_smalleigen=0.;
    double e3_reg=0., e3_small=0., e3_blas=0., e3_eigen=0., e3_smalleigen=0.;
    double e4_reg=0., e4_small=0., e4_blas=0., e4_eigen=0., e4_smalleigen=0.;
    double e5_reg=0., e5_small=0., e5_blas=0., e5_eigen=0., e5_smalleigen=0.;
    double e6_reg=0., e6_small=0., e6_blas=0., e6_eigen=0., e6_smalleigen=0.;
    double e7_reg=0., e7_small=0., e7_blas=0., e7_eigen=0., e7_smalleigen=0.;
    double e8_reg=0., e8_small=0., e8_blas=0., e8_eigen=0., e8_smalleigen=0.;
#endif

    for (int i=0; i<nloops1; ++i) {

#if 1 // B = A
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& B0k = B0[k];
                for(int i=0;i<N;++i) B0k(i) = A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B1[k] = A1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(B1[k]-B0[k]);
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B2[k] = A2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(B2[k]-B0[k]);
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(copy) (N,BP(A3[k].cptr()),1,BP(B3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(B3[k]-B0[k]);
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B4[k] = A4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e1_eigen += tmv::TMV_NORM(B4[k](i)-B0[k](i));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            B5[k] = A5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e1_smalleigen += tmv::TMV_NORM(B5[k](i)-B0[k](i));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // B = 8 * A
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& B0k = B0[k];
                for(int i=0;i<N;++i) B0k(i) = RT(8) * A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B1[k] = 8 * A1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) e2_reg += NormSq(B1[k]-B0[k]);
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= 8.*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B2[k] = 8 * A2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) e2_small += NormSq(B2[k]-B0[k]);
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= 8.*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(copy) (N,BP(A3[k].cptr()),1,BP(B3[k].ptr()),1);
            BLASNAME(scal) (N,eight,BP(B3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) e2_blas += NormSq(B3[k]-B0[k]);
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= 8.*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B4[k] = 8 * A4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e2_eigen += tmv::TMV_NORM(B4[k](i)-B0[k](i));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= 8.*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            B5[k] = 8 * A5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e2_smalleigen += tmv::TMV_NORM(B5[k](i)-B0[k](i));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= 8.*Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // B *= 7
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                tmv::Vector<T>& B0k = B0[k];
                for(int i=0;i<N;++i) B0k(i) *= RT(7);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B1[k] *= 7;

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) e3_reg += NormSq(B1[k]-B0[k]);
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B2[k] *= 7;

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) e3_small += NormSq(B2[k]-B0[k]);
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(scal) (N,seven,BP(B3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) e3_blas += NormSq(B3[k]-B0[k]);
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B4[k] *= 7;

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e3_eigen += tmv::TMV_NORM(B4[k](i)-B0[k](i));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            B5[k] *= 7;

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e3_smalleigen += tmv::TMV_NORM(B5[k](i)-B0[k](i));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // B = -A
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& B0k = B0[k];
                for(int i=0;i<N;++i) B0k(i) = -A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B1[k] = -A1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) e4_reg += NormSq(B1[k]-B0[k]);
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B2[k] = -A2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) e4_small += NormSq(B2[k]-B0[k]);
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(copy) (N,BP(A3[k].cptr()),1,BP(B3[k].ptr()),1);
            BLASNAME(scal) (N,mone,BP(B3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) e4_blas += NormSq(B3[k]-B0[k]);
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B4[k] = -A4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e4_eigen += tmv::TMV_NORM(B4[k](i)-B0[k](i));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            B5[k] = -A5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e4_smalleigen += tmv::TMV_NORM(B5[k](i)-B0[k](i));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // B = -B
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& B0k = B0[k];
                for(int i=0;i<N;++i) B0k(i) = -B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B1[k] = -B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) e5_reg += NormSq(B1[k]-B0[k]);
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B2[k] = -B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) e5_small += NormSq(B2[k]-B0[k]);
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(scal) (N,mone,BP(B3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) e5_blas += NormSq(B3[k]-B0[k]);
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B4[k] = -B4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e5_eigen += tmv::TMV_NORM(B4[k](i)-B0[k](i));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            B5[k] = -B5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e5_smalleigen += tmv::TMV_NORM(B5[k](i)-B0[k](i));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // B *= complex(8,9)
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                tmv::Vector<T>& B0k = B0[k];
                for(int i=0;i<N;++i) B0k(i) *= T(8,9);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B1[k] *= T(8,9);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) e6_reg += NormSq(B1[k]-B0[k]);
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B2[k] *= T(8,9);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) e6_small += NormSq(B2[k]-B0[k]);
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(scal) (N,z89,BP(B3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) e6_blas += NormSq(B3[k]-B0[k]);
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B4[k] *= T(8,9);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e6_eigen += tmv::TMV_NORM(B4[k](i)-B0[k](i));
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            B5[k] *= T(8,9);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e6_smalleigen += tmv::TMV_NORM(B5[k](i)-B0[k](i));
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(B0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // B = complex(8,9) * A
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& B0k = B0[k];
                for(int i=0;i<N;++i) B0k(i) = T(8,9) * A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B1[k] = T(8,9) * A1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) e7_reg += NormSq(B1[k]-B0[k]);
            e7_reg = sqrt(e7_reg/nloops2);
            e7_reg /= 12.*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B2[k] = T(8,9) * A2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) e7_small += NormSq(B2[k]-B0[k]);
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= 12.*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(copy) (N,BP(A3[k].cptr()),1,BP(B3[k].ptr()),1);
            BLASNAME(scal) (N,z89,BP(B3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) e7_blas += NormSq(B3[k]-B0[k]);
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= 12.*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B4[k] = T(8,9) * A4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e7_eigen += tmv::TMV_NORM(B4[k](i)-B0[k](i));
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= 12.*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            B5[k] = T(8,9) * A5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e7_smalleigen += tmv::TMV_NORM(B5[k](i)-B0[k](i));
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= 12.*Norm(A0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // B = complex(8,9) * A*
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& B0k = B0[k];
                for(int i=0;i<N;++i) B0k(i) = T(8,9) * std::conj(A0k(i));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B1[k] = T(8,9) * A1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) e8_reg += NormSq(B1[k]-B0[k]);
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= 12.*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B2[k] = T(8,9) * A2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) e8_small += NormSq(B2[k]-B0[k]);
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= 12.*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(copy) (N,BP(A3[k].cptr()),1,BP(B3[k].ptr()),1);
            BLASDNAME(scal) (N,dmone,B3[k].realPart().ptr()+1,2);
            BLASNAME(scal) (N,z89,BP(B3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) e8_blas += NormSq(B3[k]-B0[k]);
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= 12.*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            B4[k] = T(8,9) * A4[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e8_eigen += tmv::TMV_NORM(B4[k](i)-B0[k](i));
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= 12.*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            B5[k] = T(8,9) * A5[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e8_smalleigen += tmv::TMV_NORM(B5[k](i)-B0[k](i));
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= 12.*Norm(A0[0]);
        }
#endif
#endif
#endif
#endif
    }

    std::cout<<"B = A                   "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"B = 8*A                 "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"B *= 7                  "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"B = -A                  "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"B = -B                  "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"B *= (8,9)              "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"B = (8,9) * A           "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"B = (8,9) * A*          "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"B = A                   "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"B = 8 * A               "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"B *= 7                  "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"B = -A                  "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"B = -B                  "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"B *= (8,9)              "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
    std::cout<<"B = (8,9) * A           "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
    std::cout<<"B = (8,9) * A*          "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
#endif
    std::cout<<"\n\n";
#endif

}
#endif

#ifdef DOADDVV
static void AddVV(
    const std::vector<tmv::Vector<T> >& A1,
    const std::vector<tmv::Vector<T> >& B1,
    std::vector<tmv::Vector<T> >& C1)
{
#ifdef ERRORCHECK
    std::vector<tmv::Vector<T> > A0 = A1;
    std::vector<tmv::Vector<T> > B0 = B1;
    std::vector<tmv::Vector<T> > C0 = C1;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallVector<T,N> > A2(nloops2);
    std::vector<tmv::SmallVector<T,N> > B2(nloops2);
    std::vector<tmv::SmallVector<T,N> > C2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; B2[k] = B1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Vector<T> > A3 = A1;
    std::vector<tmv::Vector<T> > B3 = B1;
    std::vector<tmv::Vector<T> > C3 = C1;
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENV> A4(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENV> B4(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENV> C4(nloops2,Eigen::EIGENV(N));
    for(int k=0;k<nloops2;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        const tmv::Vector<T>& B1k = B1[k];
        Eigen::EIGENV& A4k = A4[k];
        Eigen::EIGENV& B4k = B4[k];
        for(int i=0;i<N;++i) A4k(i) = A1k(i);
        for(int i=0;i<N;++i) B4k(i) = B1k(i);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,N,1> > A5;
    std::vector<Eigen::Matrix<T,N,1> > B5;
    std::vector<Eigen::Matrix<T,N,1> > C5;
    if (nloops2x) 
    { A5.resize(nloops2x); B5.resize(nloops2x); C5.resize(nloops2x); }
    for(int k=0;k<nloops2x;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        const tmv::Vector<T>& B1k = B1[k];
        Eigen::Matrix<T,N,1>& A5k = A5[k];
        Eigen::Matrix<T,N,1>& B5k = B5[k];
        for(int i=0;i<N;++i) A5k(i) = A1k(i);
        for(int i=0;i<N;++i) B5k(i) = B1k(i);
    }
#endif

    timeval tp;

    double t1_reg=0., t1_small=0., t1_blas=0., t1_eigen=0., t1_smalleigen=0.;
    double t2_reg=0., t2_small=0., t2_blas=0., t2_eigen=0., t2_smalleigen=0.;
    double t3_reg=0., t3_small=0., t3_blas=0., t3_eigen=0., t3_smalleigen=0.;
    double t4_reg=0., t4_small=0., t4_blas=0., t4_eigen=0., t4_smalleigen=0.;
    double t5_reg=0., t5_small=0., t5_blas=0., t5_eigen=0., t5_smalleigen=0.;
    double t6_reg=0., t6_small=0., t6_blas=0., t6_eigen=0., t6_smalleigen=0.;
    double t7_reg=0., t7_small=0., t7_blas=0., t7_eigen=0., t7_smalleigen=0.;
    double t8_reg=0., t8_small=0., t8_blas=0., t8_eigen=0., t8_smalleigen=0.;
    double t9_reg=0., t9_small=0., t9_blas=0., t9_eigen=0., t9_smalleigen=0.;
    double t10_reg=0., t10_small=0., t10_blas=0., t10_eigen=0., t10_smalleigen=0.;
    double t11_reg=0., t11_small=0., t11_blas=0., t11_eigen=0., t11_smalleigen=0.;
    double ta,tb;

#ifdef ERRORCHECK
    double e1_reg=0., e1_small=0., e1_blas=0., e1_eigen=0., e1_smalleigen=0.;
    double e2_reg=0., e2_small=0., e2_blas=0., e2_eigen=0., e2_smalleigen=0.;
    double e3_reg=0., e3_small=0., e3_blas=0., e3_eigen=0., e3_smalleigen=0.;
    double e4_reg=0., e4_small=0., e4_blas=0., e4_eigen=0., e4_smalleigen=0.;
    double e5_reg=0., e5_small=0., e5_blas=0., e5_eigen=0., e5_smalleigen=0.;
    double e6_reg=0., e6_small=0., e6_blas=0., e6_eigen=0., e6_smalleigen=0.;
    double e7_reg=0., e7_small=0., e7_blas=0., e7_eigen=0., e7_smalleigen=0.;
    double e8_reg=0., e8_small=0., e8_blas=0., e8_eigen=0., e8_smalleigen=0.;
    double e9_reg=0., e9_small=0., e9_blas=0., e9_eigen=0., e9_smalleigen=0.;
    double e10_reg=0., e10_small=0., e10_blas=0., e10_eigen=0., e10_smalleigen=0.;
    double e11_reg=0., e11_small=0., e11_blas=0., e11_eigen=0., e11_smalleigen=0.;
#endif

    for (int i=0; i<nloops1; ++i) {

#if 1 // C = A + B
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = A0k(i) + B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = A1[k] + B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(C1[k]-C0[k]);
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = A2[k] + B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(C2[k]-C0[k]);
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(copy) (N,BP(B3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(axpy) (N,one,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(C3[k]-C0[k]);
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] = A4[k] + B4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e1_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] = A5[k] + B5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e1_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C += A
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) += A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] += A1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) e2_reg += NormSq(C1[k]-C0[k]);
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] += A2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) e2_small += NormSq(C2[k]-C0[k]);
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(axpy) (N,one,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) e2_blas += NormSq(C3[k]-C0[k]);
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] += A4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e2_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] += A5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e2_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // C -= A
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) -= A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] -= A1[k];


        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) e3_reg += NormSq(C1[k]-C0[k]);
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] -= A2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) e3_small += NormSq(C2[k]-C0[k]);
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(axpy) (N,mone,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) e3_blas += NormSq(C3[k]-C0[k]);
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] -= A4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e3_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] -= A5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e3_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // C += 8 * A
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) += RT(8) * A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] += 8 * A1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) e4_reg += NormSq(C1[k]-C0[k]);
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] += 8 * A2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) e4_small += NormSq(C2[k]-C0[k]);
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(axpy) (N,eight,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) e4_blas += NormSq(C3[k]-C0[k]);
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] += 8 * A4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e4_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] += 8 * A5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e4_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // C = A - B
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = A0k(i) - B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = A1[k] - B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) e5_reg += NormSq(C1[k]-C0[k]);
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = A2[k] - B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) e5_small += NormSq(C2[k]-C0[k]);
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(copy) (N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(axpy) (N,mone,BP(B3[k].cptr()),1,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) e5_blas += NormSq(C3[k]-C0[k]);
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] = A4[k] - B4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e5_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] = A5[k] - B5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e5_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = 7*A - 12*B
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = RT(7) * A0k(i) - RT(12) * B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = 7 * A1[k] - 12 * B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) e6_reg += NormSq(C1[k]-C0[k]);
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = 7 * A2[k] - 12 * B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) e6_small += NormSq(C2[k]-C0[k]);
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(copy) (N,BP(B3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(scal) (N,mtwelve,BP(C3[k].ptr()),1);
            BLASNAME(axpy) (N,seven,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) e6_blas += NormSq(C3[k]-C0[k]);
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] = 7 * A4[k] - 12 * B4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e6_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] = 7 * A5[k] - 12 * B5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e6_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = 7*C - 12*B
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = RT(7) * C0k(i) - RT(12) * B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = 7 * C1[k] - 12 * B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) e7_reg += NormSq(C1[k]-C0[k]);
            e7_reg = sqrt(e7_reg/nloops2);
            e7_reg /= 7.*Norm(C0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = 7 * C2[k] - 12 * B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) e7_small += NormSq(C2[k]-C0[k]);
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= 7.*Norm(C0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(scal) (N,seven,BP(C3[k].ptr()),1);
            BLASNAME(axpy) (N,mtwelve,BP(B3[k].cptr()),1,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) e7_blas += NormSq(C3[k]-C0[k]);
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= 7.*Norm(C0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] = 7 * C4[k] - 12 * B4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e7_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= 7.*Norm(C0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] = 7 * C5[k] - 12 * B5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e7_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= 7.*Norm(C0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C += complex(8,9) * A
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) += T(8,9) * A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] += T(8,9) * A1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) e8_reg += NormSq(C1[k]-C0[k]);
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= Norm(C0[0]);
            //std::cout<<"e8_reg = "<<e8_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] += T(8,9) * A2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) e8_small += NormSq(C2[k]-C0[k]);
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(axpy) (N,z89,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) e8_blas += NormSq(C3[k]-C0[k]);
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] += T(8,9) * A4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e8_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] += T(8,9) * A5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e8_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // C = 7*A + (8,9)*B
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = RT(7) * A0k(i) + T(8,9) * B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = 7 * A1[k] + T(8,9) * B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_reg = 0.;
            for (int k=0; k<nloops2; ++k) e9_reg += NormSq(C1[k]-C0[k]);
            e9_reg = sqrt(e9_reg/nloops2);
            e9_reg /= 7.*Norm(C0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = 7 * A2[k] + T(8,9) * B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_small = 0.;
            for (int k=0; k<nloops2; ++k) e9_small += NormSq(C2[k]-C0[k]);
            e9_small = sqrt(e9_small/nloops2);
            e9_small /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(copy) (N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(scal) (N,seven,BP(C3[k].ptr()),1);
            BLASNAME(axpy) (N,z89,BP(B3[k].cptr()),1,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_blas = 0.;
            for (int k=0; k<nloops2; ++k) e9_blas += NormSq(C3[k]-C0[k]);
            e9_blas = sqrt(e9_blas/nloops2);
            e9_blas /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] = 7 * A4[k] + T(8,9) * B4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e9_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e9_eigen = sqrt(e9_eigen/nloops2);
            e9_eigen /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] = 7 * A5[k] + T(8,9) * B5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e9_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e9_smalleigen = sqrt(e9_smalleigen/nloops2);
            e9_smalleigen /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // C = (7,1)*A + (8,9)*B
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = T(7,1) * A0k(i) + T(8,9) * B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = T(7,1) * A1[k] + T(8,9) * B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e10_reg = 0.;
            for (int k=0; k<nloops2; ++k) e10_reg += NormSq(C1[k]-C0[k]);
            e10_reg = sqrt(e10_reg/nloops2);
            e10_reg /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = T(7,1) * A2[k] + T(8,9) * B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e10_small = 0.;
            for (int k=0; k<nloops2; ++k) e10_small += NormSq(C2[k]-C0[k]);
            e10_small = sqrt(e10_small/nloops2);
            e10_small /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(copy) (N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(scal) (N,z71,BP(C3[k].ptr()),1);
            BLASNAME(axpy) (N,z89,BP(B3[k].cptr()),1,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e10_blas = 0.;
            for (int k=0; k<nloops2; ++k) e10_blas += NormSq(C3[k]-C0[k]);
            e10_blas = sqrt(e10_blas/nloops2);
            e10_blas /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] = T(7,1) * A4[k] + T(8,9) * B4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e10_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e10_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e10_eigen = sqrt(e10_eigen/nloops2);
            e10_eigen /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] = T(7,1) * A5[k] + T(8,9) * B5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e10_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e10_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e10_smalleigen = sqrt(e10_smalleigen/nloops2);
            e10_smalleigen /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // C = (7,1)*A - B
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = T(7,1) * A0k(i) - B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = T(7,1) * A1[k] - B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e11_reg = 0.;
            for (int k=0; k<nloops2; ++k) e11_reg += NormSq(C1[k]-C0[k]);
            e11_reg = sqrt(e11_reg/nloops2);
            e11_reg /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = T(7,1) * A2[k] - B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e11_small = 0.;
            for (int k=0; k<nloops2; ++k) e11_small += NormSq(C2[k]-C0[k]);
            e11_small = sqrt(e11_small/nloops2);
            e11_small /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(copy) (N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(scal) (N,z71,BP(C3[k].ptr()),1);
            BLASNAME(axpy) (N,mone,BP(B3[k].cptr()),1,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e11_blas = 0.;
            for (int k=0; k<nloops2; ++k) e11_blas += NormSq(C3[k]-C0[k]);
            e11_blas = sqrt(e11_blas/nloops2);
            e11_blas /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] = T(7,1) * A4[k] - B4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e11_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e11_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e11_eigen = sqrt(e11_eigen/nloops2);
            e11_eigen /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] = T(7,1) * A5[k] - B5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e11_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e11_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e11_smalleigen = sqrt(e11_smalleigen/nloops2);
            e11_smalleigen /= 7.*Norm(A0[0]) + 12.*Norm(B0[0]);
        }
#endif
#endif
#endif
#endif
    }

    std::cout<<"C = A + B               "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"C += A                  "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"C -= A                  "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"C += 8*A                "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"C = A - B               "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
    std::cout<<"C = 7*A - 12*B          "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"C = 7*C - 12*B          "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C += (8,9)*A            "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
    std::cout<<"C = 7*A + (8,9)*B       "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
    std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
    std::cout<<"C = (7,1)*A + (8,9)*B   "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
    std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
    std::cout<<"C = (7,1)*A - B         "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
    std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"C = A + B               "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"C += A                  "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"C -= A                  "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"C += 8*A                "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"C = A - B               "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
    std::cout<<"C = 7*A - 12*B          "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
    std::cout<<"C = 7*C - 12*B          "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C += (8,9)*A            "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
    std::cout<<"C = 7*A + (8,9)*B       "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
    std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
    std::cout<<"C = (7,1)*A + (8,9)*B   "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
    std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
    std::cout<<"C = (7,1)*A - B         "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
    std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
#endif
    std::cout<<"\n\n";
#endif
}
#endif

#ifdef DOMULTVV
static void MultVV(
    const std::vector<tmv::Vector<T> >& A1,
    const std::vector<tmv::Vector<T> >& B1)
{
#ifdef ERRORCHECK
    std::vector<tmv::Vector<T> > A0 = A1;
    std::vector<tmv::Vector<T> > B0 = B1;
    std::vector<T> D0(nloops2);
#endif

    std::vector<T> D1(nloops2);

#ifdef DOSMALL
    std::vector<tmv::SmallVector<T,N> > A2(nloops2);
    std::vector<tmv::SmallVector<T,N> > B2(nloops2);
    std::vector<T> D2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; B2[k] = B1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Vector<T> > A3 = A1;
    std::vector<tmv::Vector<T> > B3 = B1;
    std::vector<T> D3(nloops2);
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENV> A4(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENV> B4(nloops2,Eigen::EIGENV(N));
    std::vector<T> D4(nloops2);
    for(int k=0;k<nloops2;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        const tmv::Vector<T>& B1k = B1[k];
        Eigen::EIGENV& A4k = A4[k];
        Eigen::EIGENV& B4k = B4[k];
        for(int i=0;i<N;++i) A4k(i) = A1k(i);
        for(int i=0;i<N;++i) B4k(i) = B1k(i);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,N,1> > A5;
    std::vector<Eigen::Matrix<T,N,1> > B5;
    if (nloops2x) { A5.resize(nloops2x); B5.resize(nloops2x); }
    std::vector<T> D5(nloops2x);
    for(int k=0;k<nloops2x;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        const tmv::Vector<T>& B1k = B1[k];
        Eigen::Matrix<T,N,1>& A5k = A5[k];
        Eigen::Matrix<T,N,1>& B5k = B5[k];
        for(int i=0;i<N;++i) A5k(i) = A1k(i);
        for(int i=0;i<N;++i) B5k(i) = B1k(i);
    }
#endif

    timeval tp;

    double t1_reg=0., t1_small=0., t1_blas=0., t1_eigen=0., t1_smalleigen=0.;
    double t2_reg=0., t2_small=0., t2_blas=0., t2_eigen=0., t2_smalleigen=0.;
    double t3_reg=0., t3_small=0., t3_blas=0., t3_eigen=0., t3_smalleigen=0.;
    double t4_reg=0., t4_small=0., t4_blas=0., t4_eigen=0., t4_smalleigen=0.;
    double t5_reg=0., t5_small=0., t5_blas=0., t5_eigen=0., t5_smalleigen=0.;
    double t6_reg=0., t6_small=0., t6_blas=0., t6_eigen=0., t6_smalleigen=0.;
    double ta,tb;

#ifdef ERRORCHECK
    double e1_reg=0., e1_small=0., e1_blas=0., e1_eigen=0., e1_smalleigen=0.;
    double e2_reg=0., e2_small=0., e2_blas=0., e2_eigen=0., e2_smalleigen=0.;
    double e3_reg=0., e3_small=0., e3_blas=0., e3_eigen=0., e3_smalleigen=0.;
    double e4_reg=0., e4_small=0., e4_blas=0., e4_eigen=0., e4_smalleigen=0.;
    double e5_reg=0., e5_small=0., e5_blas=0., e5_eigen=0., e5_smalleigen=0.;
    double e6_reg=0., e6_small=0., e6_blas=0., e6_eigen=0., e6_smalleigen=0.;
#endif

    for (int i=0; i<nloops1; ++i) {

#if 1 // D = A * B
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                T& D0k = D0[k];
                D0k = T(0);
                for(int i=0;i<N;++i) D0k += A0k(i) * B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = A1[k] * B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_reg += std::abs((D1[k]-D0[k])*(D1[k]-D0[k]));
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = A2[k] * B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_small += std::abs((D2[k]-D0[k])*(D2[k]-D0[k]));
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
#ifdef CBLAS
            BLASNAME(dotu)_sub (N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1,BP(&D3[k]));
#else
            BLASNAME(dotu) (BP(&D3[k]),N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1);
#endif
#else
            D3[k] = BLASNAME(dot) (N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_blas += std::abs((D3[k]-D0[k])*(D3[k]-D0[k]));
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
#if 1
#ifdef TISCOMPLEX
            // Apparently Eigen automatically makes the dot product
            // use the conjugate of the second vector.  
            // So we need to undo that.
            D4[k] = A4[k].dot(B4[k].conjugate());
#else
        D4[k] = A4[k].dot(B4[k]);
#endif
#else
        // Or it can be done this way.
        // But this is much, much slower.
        D4[k] = (A4[k].transpose() * B4[k])(0);
#endif

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_eigen += std::abs((D4[k]-D0[k])*(D4[k]-D0[k]));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
#ifdef TISCOMPLEX
            D5[k] = A5[k].dot(B5[k].conjugate());
#else
        D5[k] = A5[k].dot(B5[k]);
#endif

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e1_smalleigen += std::abs((D5[k]-D0[k])*(D5[k]-D0[k]));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // D = 8 * A * B
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                T& D0k = D0[k];
                D0k = T(0);
                for(int i=0;i<N;++i) D0k += A0k(i) * B0k(i);
                D0k *= RT(8);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = 8 * A1[k] * B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_reg += std::abs((D1[k]-D0[k])*(D1[k]-D0[k]));
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= 8.*(Norm(A0[0]) * Norm(B0[0]));
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = 8 * A2[k] * B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_small += std::abs((D2[k]-D0[k])*(D2[k]-D0[k]));
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= 8.*(Norm(A0[0]) * Norm(B0[0]));
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
#ifdef CBLAS
            BLASNAME(dotu)_sub (N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1,BP(&D3[k]));
#else
            BLASNAME(dotu) (BP(&D3[k]),N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1);
#endif
            D3[k] *= 8;
#else
            D3[k] = 8 * BLASNAME(dot) (N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_blas += std::abs((D3[k]-D0[k])*(D3[k]-D0[k]));
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= 8.*(Norm(A0[0]) * Norm(B0[0]));
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
#ifdef TISCOMPLEX
            D4[k] = RT(8) * A4[k].dot(B4[k].conjugate());
#else
        D4[k] = RT(8) * A4[k].dot(B4[k]);
#endif

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_eigen += std::abs((D4[k]-D0[k])*(D4[k]-D0[k]));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= 8.*(Norm(A0[0]) * Norm(B0[0]));
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
#ifdef TISCOMPLEX
            D5[k] = RT(8) * A5[k].dot(B5[k].conjugate());
#else
        D5[k] = RT(8) * A5[k].dot(B5[k]);
#endif

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e2_smalleigen += std::abs((D5[k]-D0[k])*(D5[k]-D0[k]));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= 8.*(Norm(A0[0]) * Norm(B0[0]));
        }
#endif
#endif
#endif

#if 1 // D = 7*A * 12*B
        // This is really just showing off a feature of TMV
        // that it automatically extracts the two scaling factors and 
        // multiplies them outside of the vector product.
        // It's probably a useless feature in that nobody probably ever 
        // uses it, but there it is.
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                T& D0k = D0[k];
                D0k = T(0);
                for(int i=0;i<N;++i) D0k += A0k(i) * B0k(i);
                D0k *= RT(7) * RT(12);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = (7 * A1[k]) * (12 * B1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_reg += std::abs((D1[k]-D0[k])*(D1[k]-D0[k]));
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= 7.*12.*(Norm(A0[0]) * Norm(B0[0]));
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = (7 * A2[k]) * (12 * B2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_small += std::abs((D2[k]-D0[k])*(D2[k]-D0[k]));
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= 7.*12.*(Norm(A0[0]) * Norm(B0[0]));
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
#ifdef CBLAS
            BLASNAME(dotu)_sub (N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1,BP(&D3[k]));
#else
            BLASNAME(dotu) (BP(&D3[k]),N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1);
#endif
            D3[k] *= 7 * 12;
#else
            D3[k] = 7 * 12 * BLASNAME(dot) (N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_blas += std::abs((D3[k]-D0[k])*(D3[k]-D0[k]));
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= 7.*12.*(Norm(A0[0]) * Norm(B0[0]));
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
#ifdef TISCOMPLEX
            D4[k] = RT(7) * RT(12) * A4[k].dot(B4[k].conjugate());
#else
        D4[k] = RT(7) * RT(12) * A4[k].dot(B4[k]);
#endif

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_eigen += std::abs((D4[k]-D0[k])*(D4[k]-D0[k]));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= 7.*12.*(Norm(A0[0]) * Norm(B0[0]));
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
#ifdef TISCOMPLEX
            D5[k] = RT(7) * RT(12) * A5[k].dot(B5[k].conjugate());
#else
        D5[k] = RT(7) * RT(12) * A5[k].dot(B5[k]);
#endif

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e3_smalleigen += std::abs((D5[k]-D0[k])*(D5[k]-D0[k]));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= 7.*12.*(Norm(A0[0]) * Norm(B0[0]));
        }
#endif
#endif
#endif

#if 1 // D = A * B.conjugate()
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                T& D0k = D0[k];
                D0k = T(0);
                for(int i=0;i<N;++i) D0k += A0k(i) * std::conj(B0k(i));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = A1[k] * B1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_reg += std::abs((D1[k]-D0[k])*(D1[k]-D0[k]));
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = A2[k] * B2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_small += std::abs((D2[k]-D0[k])*(D2[k]-D0[k]));
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef CBLAS
            BLASNAME(dotc)_sub (N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1,BP(&D3[k]));
#else
            BLASNAME(dotc) (BP(&D3[k]),N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1);
#endif
            D3[k] = std::conj(D3[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_blas += std::abs((D3[k]-D0[k])*(D3[k]-D0[k]));
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            // Here the automatic conjugate of B is correct.
            D4[k] = A4[k].dot(B4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_eigen += std::abs((D4[k]-D0[k])*(D4[k]-D0[k]));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = A5[k].dot(B5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e4_smalleigen += std::abs((D5[k]-D0[k])*(D5[k]-D0[k]));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // D = A.conjugate() * B
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                T& D0k = D0[k];
                D0k = T(0);
                for(int i=0;i<N;++i) D0k += std::conj(A0k(i)) * B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = A1[k].conjugate() * B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e5_reg += std::abs((D1[k]-D0[k])*(D1[k]-D0[k]));
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = A2[k].conjugate() * B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e5_small += std::abs((D2[k]-D0[k])*(D2[k]-D0[k]));
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef CBLAS
            BLASNAME(dotc)_sub (N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1,BP(&D3[k]));
#else
            BLASNAME(dotc) (BP(&D3[k]),N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e5_blas += std::abs((D3[k]-D0[k])*(D3[k]-D0[k]));
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = std::conj(A4[k].dot(B4[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e5_eigen += std::abs((D4[k]-D0[k])*(D4[k]-D0[k]));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = A5[k].conjugate().dot(B5[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e5_smalleigen += std::abs((D5[k]-D0[k])*(D5[k]-D0[k]));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // D = A.conjugate() * B.conjugate()
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                T& D0k = D0[k];
                D0k = T(0);
                for(int i=0;i<N;++i) D0k += std::conj(A0k(i)) * std::conj(B0k(i));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = A1[k].conjugate() * B1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e6_reg += std::abs((D1[k]-D0[k])*(D1[k]-D0[k]));
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = A2[k].conjugate() * B2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e6_small += std::abs((D2[k]-D0[k])*(D2[k]-D0[k]));
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef CBLAS
            BLASNAME(dotu)_sub (N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1,BP(&D3[k]));
#else
            BLASNAME(dotu) (BP(&D3[k]),N,BP(A3[k].cptr()),1,BP(B3[k].cptr()),1);
#endif
            D3[k] = std::conj(D3[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e6_blas += std::abs((D3[k]-D0[k])*(D3[k]-D0[k]));
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = A4[k].conjugate().dot(B4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e6_eigen += std::abs((D4[k]-D0[k])*(D4[k]-D0[k]));
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = A5[k].conjugate().dot(B5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e6_smalleigen += std::abs((D5[k]-D0[k])*(D5[k]-D0[k]));
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(A0[0]) * Norm(B0[0]);
        }
#endif
#endif
#endif
#endif
    }

    std::cout<<"D = A * B               "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"D = 8 * A * B           "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"D = 7*A * 12*B          "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"D = A * conj(B)         "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"D = conj(A) * B         "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
    std::cout<<"D = conj(A) * conj(B)   "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"D = A * B               "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"D = 8 * A * B           "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"D = 7*A * 12*B          "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"D = A * conj(B)         "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"D = conj(A) * B         "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
    std::cout<<"D = conj(A) * conj(B)   "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
#endif
    std::cout<<"\n\n";
#endif
}
#endif

// Used for Non-Blas SumAbsElements function
#ifdef TISCOMPLEX
struct Summer
{
    T operator()(const T& sum, const T& el)
    { return sum + std::abs(el); }
};
#endif

#ifdef DONORM
static void NormV(
    const std::vector<tmv::Vector<T> >& A1)
{
#ifdef ERRORCHECK
    std::vector<tmv::Vector<T> > A0 = A1;
    std::vector<T> D0(nloops2);
    std::vector<RT> F0(nloops2);
#endif

    std::vector<T> D1(nloops2);
    std::vector<RT> F1(nloops2);

#ifdef DOSMALL
    std::vector<tmv::SmallVector<T,N> > A2(nloops2);
    std::vector<T> D2(nloops2);
    std::vector<RT> F2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Vector<T> > A3 = A1;
    std::vector<T> D3(nloops2);
    std::vector<RT> F3(nloops2);
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENV> A4(nloops2,Eigen::EIGENV(N));
    std::vector<T> D4(nloops2);
    std::vector<RT> F4(nloops2);
    for(int k=0;k<nloops2;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        Eigen::EIGENV& A4k = A4[k];
        for(int i=0;i<N;++i) A4k(i) = A1k(i);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,N,1> > A5;
    if (nloops2x) { A5.resize(nloops2x); }
    std::vector<T> D5(nloops2x);
    std::vector<RT> F5(nloops2x);
    for(int k=0;k<nloops2x;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        Eigen::Matrix<T,N,1>& A5k = A5[k];
        for(int i=0;i<N;++i) A5k(i) = A1k(i);
    }
#endif

    timeval tp;

    double t1_reg=0., t1_small=0., t1_blas=0., t1_eigen=0., t1_smalleigen=0.;
    double t2_reg=0., t2_small=0., t2_blas=0., t2_eigen=0., t2_smalleigen=0.;
    double t3_reg=0., t3_small=0., t3_blas=0., t3_eigen=0., t3_smalleigen=0.;
    double t4_reg=0., t4_small=0., t4_blas=0., t4_eigen=0., t4_smalleigen=0.;
    double t5_reg=0., t5_small=0., t5_blas=0., t5_eigen=0., t5_smalleigen=0.;
    double t6_reg=0., t6_small=0., t6_blas=0., t6_eigen=0., t6_smalleigen=0.;
    double t7_reg=0., t7_small=0., t7_blas=0., t7_eigen=0., t7_smalleigen=0.;
    double t8_reg=0., t8_small=0., t8_blas=0., t8_eigen=0., t8_smalleigen=0.;
    double ta,tb;

#ifdef ERRORCHECK
    double e1_reg=0., e1_small=0., e1_blas=0., e1_eigen=0., e1_smalleigen=0.;
    double e2_reg=0., e2_small=0., e2_blas=0., e2_eigen=0., e2_smalleigen=0.;
    double e3_reg=0., e3_small=0., e3_blas=0., e3_eigen=0., e3_smalleigen=0.;
    double e4_reg=0., e4_small=0., e4_blas=0., e4_eigen=0., e4_smalleigen=0.;
    double e5_reg=0., e5_small=0., e5_blas=0., e5_eigen=0., e5_smalleigen=0.;
    double e6_reg=0., e6_small=0., e6_blas=0., e6_eigen=0., e6_smalleigen=0.;
    double e7_reg=0., e7_small=0., e7_blas=0., e7_eigen=0., e7_smalleigen=0.;
    double e8_reg=0., e8_small=0., e8_blas=0., e8_eigen=0., e8_smalleigen=0.;
#endif

    for (int i=0; i<nloops1; ++i) {

#if 1 // F = Norm1(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                RT& F0k = F0[k];
                F0k = RT(0);
                for(int i=0;i<N;++i) F0k += std::abs(A0k(i));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = Norm1(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm1(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = Norm1(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm1(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;


        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
            // BLAS doesn't have this capability
            F3[k] = real(std::accumulate(A3[k].begin(),A3[k].end(),T(0.),Summer()));
#else
            F3[k] = std::abs(BLASNAME(asum) (N,BP(A3[k].cptr()),1));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm1(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = A4[k].lpNorm<1>();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm1(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = A5[k].lpNorm<1>();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e1_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm1(A0[0]);
        }
#endif
#endif
#endif

#if 1 // F = Norm2(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                RT& F0k = F0[k];
                F0k = RT(0);
                for(int i=0;i<N;++i) F0k += std::abs(A0k(i)*A0k(i));
                F0k = std::sqrt(F0k);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = Norm2(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= Norm2(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = Norm2(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= Norm2(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F3[k] = BLASDZNAME(nrm2) (N,BP(A3[k].cptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= Norm2(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = A4[k].lpNorm<2>();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm2(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = A5[k].lpNorm<2>();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e2_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm2(A0[0]);
        }
#endif
#endif
#endif

#if 1 // F = NormInf(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                RT& F0k = F0[k];
                F0k = RT(0);
                for(int i=0;i<N;++i) 
                    if (std::abs(A0k(i)) > F0k) F0k = std::abs(A0k(i));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = NormInf(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= NormInf(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = NormInf(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= NormInf(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
            // BLAS doesn't have this capability
            tmv::Vector<T>::iterator maxi = 
                std::max_element(A3[k].begin(),A3[k].end(),
                                 tmv::Compare<tmv::Ascend,tmv::AbsComp,T>());
            F3[k] = std::abs(*maxi);
#else
            F3[k] = std::abs(A3[k](BLASINAME(amax) (N,BP(A3[k].cptr()),1) BLASM1));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= NormInf(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = A4[k].lpNorm<Eigen::Infinity>();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= NormInf(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = A5[k].lpNorm<Eigen::Infinity>();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e3_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= NormInf(A0[0]);
        }
#endif
#endif
#endif

#if 1 // F = NormSq(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                RT& F0k = F0[k];
                F0k = RT(0);
                for(int i=0;i<N;++i) F0k += std::abs(A0k(i) * A0k(i));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = NormSq(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= NormSq(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = NormSq(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= NormSq(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
#ifdef CBLAS
            BLASNAME(dotc)_sub (N,BP(A3[k].cptr()),1,BP(A3[k].cptr()),1,BP(&D3[k]));
#else
            BLASNAME(dotc) (BP(&D3[k]),N,BP(A3[k].cptr()),1,BP(A3[k].cptr()),1);
#endif
            F3[k] = real(D3[k]);
#else
            D3[k] = BLASNAME(dot) (N,BP(A3[k].cptr()),1,BP(A3[k].cptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= NormSq(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = A4[k].squaredNorm();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= NormSq(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = A5[k].squaredNorm();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e4_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= NormSq(A0[0]);
        }
#endif
#endif
#endif

#if 1 // F = Norm(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                RT& F0k = F0[k];
                F0k = RT(0);
                for(int i=0;i<N;++i) F0k += std::abs(A0k(i)*A0k(i));
                F0k = std::sqrt(F0k);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = Norm(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e5_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = Norm(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e5_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F3[k] = BLASDZNAME(nrm2) (N,BP(A3[k].cptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e5_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = A4[k].norm();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e5_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = A5[k].norm();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e5_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // D = SumElements(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                T& D0k = D0[k];
                D0k = T(0);
                for(int i=0;i<N;++i) D0k += A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = SumElements(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e6_reg += std::abs((D1[k]-D0[k])*(D1[k]-D0[k]));
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= std::abs(SumElements(A0[0]));
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = SumElements(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e6_small += std::abs((D2[k]-D0[k])*(D2[k]-D0[k]));
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= std::abs(SumElements(A0[0]));
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // No BLAS function for this.  Use STL.
            D3[k] = std::accumulate(A3[k].begin(),A3[k].end(),T(0.));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e6_blas += std::abs((D3[k]-D0[k])*(D3[k]-D0[k]));
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= std::abs(SumElements(A0[0]));
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = A4[k].sum();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e6_eigen += std::abs((D4[k]-D0[k])*(D4[k]-D0[k]));
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= std::abs(SumElements(A0[0]));
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = A5[k].sum();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e6_smalleigen += std::abs((D5[k]-D0[k])*(D5[k]-D0[k]));
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= std::abs(SumElements(A0[0]));
        }
#endif
#endif
#endif

#if 1 // F = SumAbsElements(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                RT& F0k = F0[k];
                F0k = RT(0);
                for(int i=0;i<N;++i) F0k += std::abs(A0k(i));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = SumAbsElements(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e7_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e7_reg = sqrt(e7_reg/nloops2);
            e7_reg /= SumAbsElements(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = SumAbsElements(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e7_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= SumAbsElements(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
            // BLAS doesn't have this capability
            F3[k] = real(std::accumulate(A3[k].begin(),A3[k].end(),T(0.),Summer()));
#else
            F3[k] = std::abs(BLASNAME(asum) (N,BP(A3[k].cptr()),1));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e7_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= SumAbsElements(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = A4[k].cwise().abs().sum();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e7_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= SumAbsElements(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = A5[k].cwise().abs().sum();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e7_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= SumAbsElements(A0[0]);
        }
#endif
#endif
#endif

#if 1 // F = SumAbs2Elements(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                RT& F0k = F0[k];
                F0k = RT(0);
                for(int i=0;i<N;++i) {
#ifdef TISCOMPLEX
                    RT temp = std::abs(real(A0k(i))) + std::abs(imag(A0k(i)));
#else
                    RT temp = std::abs(A0k(i));
#endif
                    F0k += temp;
                }
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = SumAbs2Elements(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e8_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= SumAbs2Elements(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = SumAbs2Elements(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e8_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= SumAbs2Elements(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F3[k] = std::abs(BLASDZNAME(asum) (N,BP(A3[k].cptr()),1));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e8_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= SumAbs2Elements(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = (
                A4[k].real().cwise().abs().sum() + 
                A4[k].imag().cwise().abs().sum());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e8_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= SumAbs2Elements(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = (
                A5[k].real().cwise().abs().sum() + 
                A5[k].imag().cwise().abs().sum());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e8_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= SumAbs2Elements(A0[0]);
        }
#endif
#endif
#endif
    }

    std::cout<<"Norm1(A)                "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"Norm2(A)                "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"NormInf(A)              "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"NormSq(A)               "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"Norm(A)                 "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
    std::cout<<"SumElements(A)          "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"SumAbsElements(A)       "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"SumAbs2Elements(A)      "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"Norm1(A)                "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"Norm2(A)                "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"NormInf(A)              "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"NormSq(A)               "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"Norm(A)                 "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
    std::cout<<"SumElements(A)          "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
    std::cout<<"SumAbsElements(A)       "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
    std::cout<<"SumAbs2Elements(A)      "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
    std::cout<<"\n\n";
#endif
}
#endif

#ifdef DOMINMAX
static void MinMaxV(
    const std::vector<tmv::Vector<T> >& A1)
{
#ifdef ERRORCHECK
    std::vector<tmv::Vector<T> > A0 = A1;
    std::vector<T> D0(nloops2);
    std::vector<int> E0(nloops2);
    std::vector<RT> F0(nloops2);
#endif

    std::vector<T> D1(nloops2);
    std::vector<int> E1(nloops2);
    std::vector<RT> F1(nloops2);

#ifdef DOSMALL
    std::vector<tmv::SmallVector<T,N> > A2(nloops2);
    std::vector<T> D2(nloops2);
    std::vector<int> E2(nloops2);
    std::vector<RT> F2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Vector<T> > A3 = A1;
    std::vector<T> D3(nloops2);
    std::vector<int> E3(nloops2);
    std::vector<RT> F3(nloops2);
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENV> A4(nloops2,Eigen::EIGENV(N));
    std::vector<T> D4(nloops2);
    std::vector<int> E4(nloops2);
    std::vector<RT> F4(nloops2);
    for(int k=0;k<nloops2;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        Eigen::EIGENV& A4k = A4[k];
        for(int i=0;i<N;++i) A4k(i) = A1k(i);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,N,1> > A5;
    if (nloops2x) { A5.resize(nloops2x); }
    std::vector<T> D5(nloops2x);
    std::vector<int> E5(nloops2);
    std::vector<RT> F5(nloops2x);
    for(int k=0;k<nloops2x;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        Eigen::Matrix<T,N,1>& A5k = A5[k];
        for(int i=0;i<N;++i) A5k(i) = A1k(i);
    }
#endif

    timeval tp;

    double t1_reg=0., t1_small=0., t1_blas=0., t1_eigen=0., t1_smalleigen=0.;
    double t2_reg=0., t2_small=0., t2_blas=0., t2_eigen=0., t2_smalleigen=0.;
    double t3_reg=0., t3_small=0., t3_blas=0., t3_eigen=0., t3_smalleigen=0.;
    double t4_reg=0., t4_small=0., t4_blas=0., t4_eigen=0., t4_smalleigen=0.;
    double t5_reg=0., t5_small=0., t5_blas=0., t5_eigen=0., t5_smalleigen=0.;
    double t6_reg=0., t6_small=0., t6_blas=0., t6_eigen=0., t6_smalleigen=0.;
    double t7_reg=0., t7_small=0., t7_blas=0., t7_eigen=0., t7_smalleigen=0.;
    double t8_reg=0., t8_small=0., t8_blas=0., t8_eigen=0., t8_smalleigen=0.;
    double t9_reg=0., t9_small=0., t9_blas=0., t9_eigen=0., t9_smalleigen=0.;
    double t10_reg=0., t10_small=0., t10_blas=0., t10_eigen=0., t10_smalleigen=0.;
    double ta,tb;

#ifdef ERRORCHECK
    double e1_reg=0., e1_small=0., e1_blas=0., e1_eigen=0., e1_smalleigen=0.;
    double e2_reg=0., e2_small=0., e2_blas=0., e2_eigen=0., e2_smalleigen=0.;
    double e3_reg=0., e3_small=0., e3_blas=0., e3_eigen=0., e3_smalleigen=0.;
    double e4_reg=0., e4_small=0., e4_blas=0., e4_eigen=0., e4_smalleigen=0.;
    double e5_reg=0., e5_small=0., e5_blas=0., e5_eigen=0., e5_smalleigen=0.;
    double e6_reg=0., e6_small=0., e6_blas=0., e6_eigen=0., e6_smalleigen=0.;
    double e7_reg=0., e7_small=0., e7_blas=0., e7_eigen=0., e7_smalleigen=0.;
    double e8_reg=0., e8_small=0., e8_blas=0., e8_eigen=0., e8_smalleigen=0.;
    double e9_reg=0., e9_small=0., e9_blas=0., e9_eigen=0., e9_smalleigen=0.;
    double e10_reg=0., e10_small=0., e10_blas=0., e10_eigen=0., e10_smalleigen=0.;
#endif

    for (int i=0; i<nloops1; ++i) {

#if 1 // D = MaxElement(A)
#ifndef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                T& D0k = D0[k];
                D0k = A0k(0);
                for(int i=0;i<N;++i) 
                    if (A0k(i) > D0k) D0k = A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) 
            D1[k] = MaxElement(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_reg += std::abs((D1[k]-D0[k])*(D1[k]-D0[k]));
            e1_reg = sqrt(e1_reg/nloops2);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = MaxElement(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_small += std::abs((D2[k]-D0[k])*(D2[k]-D0[k]));
            e1_small = sqrt(e1_small/nloops2);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D3[k] = *std::max_element(A3[k].begin(),A3[k].end());
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_blas += std::abs((D3[k]-D0[k])*(D3[k]-D0[k]));
            e1_blas = sqrt(e1_blas/nloops2);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = A4[k].maxCoeff();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_eigen += std::abs((D4[k]-D0[k])*(D4[k]-D0[k]));
            e1_eigen = sqrt(e1_eigen/nloops2);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = A5[k].maxCoeff();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e1_smalleigen += std::abs((D5[k]-D0[k])*(D5[k]-D0[k]));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
        }
#endif
#endif
#endif
#endif

#if 1 // F = MaxAbsElement(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                RT& F0k = F0[k];
                F0k = std::abs(A0k(0));
                for(int i=0;i<N;++i) 
                    if (std::abs(A0k(i)) > F0k) F0k = std::abs(A0k(i));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = MaxAbsElement(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e2_reg = sqrt(e2_reg/nloops2);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = MaxAbsElement(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e2_small = sqrt(e2_small/nloops2);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
            // BLAS doesn't have this capability
            tmv::Vector<T>::iterator maxi = 
                std::max_element(A3[k].begin(),A3[k].end(),
                                 tmv::Compare<tmv::Ascend,tmv::AbsComp,T>());
            F3[k] = std::abs(*maxi);
#else
            F3[k] = std::abs(A3[k](BLASINAME(amax) (N,BP(A3[k].cptr()),1) BLASM1));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e2_blas = sqrt(e2_blas/nloops2);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = A4[k].cwise().abs().maxCoeff();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e2_eigen = sqrt(e2_eigen/nloops2);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = A5[k].cwise().abs().maxCoeff();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e2_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
        }
#endif
#endif
#endif

#if 1 // D = MinElement(A)
#ifndef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                T& D0k = D0[k];
                D0k = A0k(0);
                for(int i=0;i<N;++i) 
                    if (A0k(i) < D0k) D0k = A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = MinElement(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_reg += std::abs((D1[k]-D0[k])*(D1[k]-D0[k]));
            e3_reg = sqrt(e3_reg/nloops2);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = MinElement(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_small += std::abs((D2[k]-D0[k])*(D2[k]-D0[k]));
            e3_small = sqrt(e3_small/nloops2);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D3[k] = *std::min_element(A3[k].begin(),A3[k].end());
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_blas += std::abs((D3[k]-D0[k])*(D3[k]-D0[k]));
            e3_blas = sqrt(e3_blas/nloops2);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = A4[k].minCoeff();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_eigen += std::abs((D4[k]-D0[k])*(D4[k]-D0[k]));
            e3_eigen = sqrt(e3_eigen/nloops2);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = A5[k].minCoeff();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e3_smalleigen += std::abs((D5[k]-D0[k])*(D5[k]-D0[k]));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
        }
#endif
#endif
#endif
#endif

#if 1 // F = MinAbsElement(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                RT& F0k = F0[k];
                F0k = std::abs(A0k(0));
                for(int i=0;i<N;++i) 
                    if (std::abs(A0k(i)) < F0k) F0k = std::abs(A0k(i));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = MinAbsElement(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e4_reg = sqrt(e4_reg/nloops2);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = MinAbsElement(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e4_small = sqrt(e4_small/nloops2);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F3[k] = std::abs(*std::min_element(
                    A3[k].begin(),A3[k].end(),
                    tmv::Compare<tmv::Ascend,tmv::AbsComp,T>()));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e4_blas = sqrt(e4_blas/nloops2);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = A4[k].cwise().abs().minCoeff();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e4_eigen = sqrt(e4_eigen/nloops2);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = A5[k].cwise().abs().minCoeff();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e4_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
        }
#endif
#endif
#endif

#if 1 // D = MaxElement(A,&i)
#ifndef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                T& D0k = D0[k];
                int& E0k = E0[k];
                D0k = A0k(0);
                E0k = 0;
                for(int i=0;i<N;++i) 
                    if (A0k(i) > D0k) { D0k = A0k(i); E0k = i; }
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = MaxElement(A1[k],&E1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                e5_reg += std::abs((D1[k]-D0[k])*(D1[k]-D0[k]));
                e5_reg += std::abs((E1[k]-E0[k])*(E1[k]-E0[k]));
            }
            e5_reg = sqrt(e5_reg/nloops2);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = MaxElement(A2[k],&E2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                e5_small += std::abs((D2[k]-D0[k])*(D2[k]-D0[k]));
                e5_small += std::abs((E2[k]-E0[k])*(E2[k]-E0[k]));
            }
            e5_small = sqrt(e5_small/nloops2);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            tmv::Vector<T>::iterator maxi = 
                std::max_element(A3[k].begin(),A3[k].end());
            D3[k] = *maxi;
            E3[k] = maxi - A3[k].begin();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                e5_blas += std::abs((D3[k]-D0[k])*(D3[k]-D0[k]));
                e5_blas += std::abs((E3[k]-E0[k])*(E3[k]-E0[k]));
            }
            e5_blas = sqrt(e5_blas/nloops2);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = A4[k].maxCoeff(&E4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e5_eigen += std::abs((D4[k]-D0[k])*(D4[k]-D0[k]));
                e5_eigen += std::abs((E4[k]-E0[k])*(E4[k]-E0[k]));
            }
            e5_eigen = sqrt(e5_eigen/nloops2);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = A5[k].maxCoeff(&E5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                e5_smalleigen += std::abs((D5[k]-D0[k])*(D5[k]-D0[k]));
                e5_smalleigen += std::abs((E5[k]-E0[k])*(E5[k]-E0[k]));
            }
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
        }
#endif
#endif
#endif
#endif

#if 1 // F = MaxAbsElement(A,&i)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                RT& F0k = F0[k];
                int& E0k = E0[k];
                F0k = std::abs(A0k(0));
                E0k = 0;
                for(int i=0;i<N;++i) 
                    if (std::abs(A0k(i)) > F0k) { F0k = std::abs(A0k(i)); E0k = i; }
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = MaxAbsElement(A1[k],&E1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                e6_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
                e6_reg += std::abs((E1[k]-E0[k])*(E1[k]-E0[k]));
            }
            e6_reg = sqrt(e6_reg/nloops2);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = MaxAbsElement(A2[k],&E2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                e6_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
                e6_small += std::abs((E2[k]-E0[k])*(E2[k]-E0[k]));
            }
            e6_small = sqrt(e6_small/nloops2);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
            // BLAS doesn't have this capability
            tmv::Vector<T>::iterator maxi = 
                std::max_element(
                    A3[k].begin(),A3[k].end(),
                    tmv::Compare<tmv::Ascend,tmv::AbsComp,T>());
            F3[k] = std::abs(*maxi);
            E3[k] = maxi - A3[k].begin();
#else
            E3[k] = BLASINAME(amax) (N,BP(A3[k].cptr()),1) BLASM1;
            F3[k] = std::abs(A3[k][E3[k]]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                e6_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
                e6_blas += std::abs((E3[k]-E0[k])*(E3[k]-E0[k]));
            }
            e6_blas = sqrt(e6_blas/nloops2);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = A4[k].cwise().abs().maxCoeff(&E4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e6_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
                e6_eigen += std::abs((E4[k]-E0[k])*(E4[k]-E0[k]));
            }
            e6_eigen = sqrt(e6_eigen/nloops2);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = A5[k].cwise().abs().maxCoeff(&E5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                e6_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
                e6_smalleigen += std::abs((E5[k]-E0[k])*(E5[k]-E0[k]));
            }
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
        }
#endif
#endif
#endif

#if 1 // D = MinElement(A,&i)
#ifndef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                T& D0k = D0[k];
                int& E0k = E0[k];
                D0k = A0k(0);
                E0k = 0;
                for(int i=0;i<N;++i) 
                    if (A0k(i) < D0k) { D0k = A0k(i); E0k = i; }
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = MinElement(A1[k],&E1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                e7_reg += std::abs((D1[k]-D0[k])*(D1[k]-D0[k]));
                e7_reg += std::abs((E1[k]-E0[k])*(E1[k]-E0[k]));
            }
            e7_reg = sqrt(e7_reg/nloops2);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = MinElement(A2[k],&E2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                e7_small += std::abs((D2[k]-D0[k])*(D2[k]-D0[k]));
                e7_small += std::abs((E2[k]-E0[k])*(E2[k]-E0[k]));
            }
            e7_small = sqrt(e7_small/nloops2);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            tmv::Vector<T>::iterator mini = 
                std::min_element(A3[k].begin(),A3[k].end());
            D3[k] = *mini;
            E3[k] = mini - A3[k].begin();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                e7_blas += std::abs((D3[k]-D0[k])*(D3[k]-D0[k]));
                e7_blas += std::abs((E3[k]-E0[k])*(E3[k]-E0[k]));
            }
            e7_blas = sqrt(e7_blas/nloops2);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = A4[k].minCoeff(&E4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e7_eigen += std::abs((D4[k]-D0[k])*(D4[k]-D0[k]));
                e7_eigen += std::abs((E4[k]-E0[k])*(E4[k]-E0[k]));
            }
            e7_eigen = sqrt(e7_eigen/nloops2);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = A5[k].minCoeff(&E5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                e7_smalleigen += std::abs((D5[k]-D0[k])*(D5[k]-D0[k]));
                e7_smalleigen += std::abs((E5[k]-E0[k])*(E5[k]-E0[k]));
            }
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
        }
#endif
#endif
#endif
#endif

#if 1 // F = MinAbsElement(A,&i)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                RT& F0k = F0[k];
                int& E0k = E0[k];
                F0k = std::abs(A0k(0));
                E0k = 0;
                for(int i=0;i<N;++i) {
                    if (std::abs(A0k(i)) < F0k) {
                        F0k = std::abs(A0k(i));
                        E0k = i; 
                    }
                }
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) 
            F1[k] = MinAbsElement(A1[k],&E1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                e8_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
                e8_reg += std::abs((E1[k]-E0[k])*(E1[k]-E0[k]));
            }
            e8_reg = sqrt(e8_reg/nloops2);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = MinAbsElement(A2[k],&E2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                e8_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
                e8_small += std::abs((E2[k]-E0[k])*(E2[k]-E0[k]));
            }
            e8_small = sqrt(e8_small/nloops2);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            tmv::Vector<T>::iterator mini = 
                std::min_element(
                    A3[k].begin(),A3[k].end(),
                    tmv::Compare<tmv::Ascend,tmv::AbsComp,T>());
            F3[k] = std::abs(*mini);
            E3[k] = mini - A3[k].begin();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                e8_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
                e8_blas += std::abs((E3[k]-E0[k])*(E3[k]-E0[k]));
            }
            e8_blas = sqrt(e8_blas/nloops2);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = A4[k].cwise().abs().minCoeff(&E4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e8_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
                e8_eigen += std::abs((E4[k]-E0[k])*(E4[k]-E0[k]));
            }
            e8_eigen = sqrt(e8_eigen/nloops2);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = A5[k].cwise().abs().minCoeff(&E5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                e8_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
                e8_smalleigen += std::abs((E5[k]-E0[k])*(E5[k]-E0[k]));
            }
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
        }
#endif
#endif
#endif

#if 1 // F = MaxAbs2Element(A)
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                RT& F0k = F0[k];
                F0k = std::abs(A0k(0));
                for(int i=0;i<N;++i) {
#ifdef TISCOMPLEX
                    RT temp = std::abs(real(A0k(i))) + std::abs(imag(A0k(i)));
#else
                    RT temp = std::abs(A0k(i));
#endif
                    if (temp > F0k) F0k = temp;
                }
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) 
            F1[k] = MaxAbs2Element(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e9_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e9_reg = sqrt(e9_reg/nloops2);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = MaxAbs2Element(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e9_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e9_small = sqrt(e9_small/nloops2);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D3[k] = A3[k](BLASINAME(amax) (N,BP(A3[k].cptr()),1) BLASM1);
#ifdef TISCOMPLEX
            F3[k] = std::abs(real(D3[k])) + std::abs(imag(D3[k]));
#else
            F3[k] = std::abs(D3[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e9_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e9_blas = sqrt(e9_blas/nloops2);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
        {
            // Eigen doesn't have this capability
            RT& best = F4[k];
            best = 0.;
            for(int i=0;i<N;i++)
            {
                RT temp = std::abs(real(A4[k][i])) + std::abs(imag(A4[k][i]));
                if (temp > best) { best = temp; }
            }
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e9_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e9_eigen = sqrt(e9_eigen/nloops2);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) 
        {
            RT& best = F5[k];
            best = 0.;
            for(int i=0;i<N;i++)
            {
                RT temp = std::abs(real(A5[k][i])) + std::abs(imag(A5[k][i]));
                if (temp > best) { best = temp; }
            }
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e9_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e9_smalleigen = sqrt(e9_smalleigen/nloops2);
        }
#endif
#endif
#endif
#endif

#if 1 // F = MaxAbs2Element(A,&i)
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                RT& F0k = F0[k];
                int& E0k = E0[k];
                F0k = std::abs(A0k(0));
                E0k = 0;
                for(int i=0;i<N;++i) {
#ifdef TISCOMPLEX
                    RT temp = std::abs(real(A0k(i))) + std::abs(imag(A0k(i)));
#else
                    RT temp = std::abs(A0k(i));
#endif
                    if (temp > F0k) { F0k = temp; E0k = i; }
                }
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = MaxAbs2Element(A1[k],&E1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e10_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                e10_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
                e10_reg += std::abs((E1[k]-E0[k])*(E1[k]-E0[k]));
            }
            e10_reg = sqrt(e10_reg/nloops2);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = MaxAbs2Element(A2[k],&E2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e10_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                e10_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
                e10_small += std::abs((E2[k]-E0[k])*(E2[k]-E0[k]));
            }
            e10_small = sqrt(e10_small/nloops2);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            E3[k] = BLASINAME(amax) (N,BP(A3[k].cptr()),1) BLASM1;
            D3[k] = A3[k](E3[k]);
#ifdef TISCOMPLEX
            F3[k] = std::abs(real(D3[k])) + std::abs(imag(D3[k]));
#else
            F3[k] = std::abs(D3[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e10_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                e10_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
                e10_blas += std::abs((E3[k]-E0[k])*(E3[k]-E0[k]));
            }
            e10_blas = sqrt(e10_blas/nloops2);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
        {
            RT& best = F4[k];
            int& ibest = E4[k];
            best = 0.;
            ibest = 0;
            for(int i=0;i<N;i++)
            {
                RT temp = std::abs(real(A4[k][i])) + std::abs(imag(A4[k][i]));
                if (temp > best) { best = temp; ibest = i; }
            }
        }


        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e10_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e10_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
                e10_eigen += std::abs((E4[k]-E0[k])*(E4[k]-E0[k]));
            }
            e10_eigen = sqrt(e10_eigen/nloops2);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
        {
            RT& best = F5[k];
            int& ibest = E5[k];
            best = 0.;
            ibest = 0;
            for(int i=0;i<N;i++)
            {
                RT temp = std::abs(real(A5[k][i])) + std::abs(imag(A5[k][i]));
                if (temp > best) { best = temp; ibest = i; }
            }
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e10_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e10_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
                e10_smalleigen += std::abs((E5[k]-E0[k])*(E5[k]-E0[k]));
            }
            e10_smalleigen = sqrt(e10_smalleigen/nloops2);
        }
#endif
#endif
#endif
#endif
    }

#ifndef TISCOMPLEX
    std::cout<<"MaxElement(A)           "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
#endif
    std::cout<<"MaxAbsElement(A)        "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
#ifndef TISCOMPLEX
    std::cout<<"MinElement(A)           "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
#endif
    std::cout<<"MinAbsElement(A)        "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
#ifndef TISCOMPLEX
    std::cout<<"MaxElement(A,&E)        "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
#endif
    std::cout<<"MaxAbsElement(A,&E)     "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
#ifndef TISCOMPLEX
    std::cout<<"MinElement(A,&E)        "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
#endif
    std::cout<<"MinAbsElement(A,&E)     "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"MaxAbs2Element(A)       "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
    std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
    std::cout<<"MaxAbs2Element(A,&E)    "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
    std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
#ifndef TISCOMPLEX
    std::cout<<"MaxElements(A)          "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
#endif
    std::cout<<"MaxAbsElement(A)        "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
#ifndef TISCOMPLEX
    std::cout<<"MinElement(A)           "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
#endif
    std::cout<<"MinAbsElement(A)        "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
#ifndef TISCOMPLEX
    std::cout<<"MaxElement(A,&E)        "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
#endif
    std::cout<<"MaxAbsElement(A,&E)     "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
#ifndef TISCOMPLEX
    std::cout<<"MinElement(A,&E))       "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
#endif
    std::cout<<"MinAbsElement(A,&E)     "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"MaxAbs2Element(A)       "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
    std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
    std::cout<<"MaxAbs2Element(A,&E)    "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
    std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
#endif
    std::cout<<"\n\n";
#endif
}
#endif

#ifdef DOSWAP
static void SwapV(
    std::vector<tmv::Vector<T> >& A1,
    std::vector<tmv::Vector<T> >& B1)
{
    // Make a plausible random permutation:
    tmv::Vector<T> A_temp = A1[0];
    int P[N];
    A_temp.sort(P);

#ifdef ERRORCHECK
    std::vector<tmv::Vector<T> > A0 = A1;
    std::vector<tmv::Vector<T> > B0 = B1;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallVector<T,N> > A2(nloops2);
    std::vector<tmv::SmallVector<T,N> > B2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; B2[k] = B1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Vector<T> > A3 = A1;
    std::vector<tmv::Vector<T> > B3 = B1;
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENV> A4(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENV> B4(nloops2,Eigen::EIGENV(N));
    for(int k=0;k<nloops2;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        const tmv::Vector<T>& B1k = B1[k];
        Eigen::EIGENV& A4k = A4[k];
        Eigen::EIGENV& B4k = B4[k];
        for(int i=0;i<N;++i) A4k(i) = A1k(i);
        for(int i=0;i<N;++i) B4k(i) = B1k(i);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,N,1> > A5;
    std::vector<Eigen::Matrix<T,N,1> > B5;
    if (nloops2x) { A5.resize(nloops2x); B5.resize(nloops2x); }
    for(int k=0;k<nloops2x;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        const tmv::Vector<T>& B1k = B1[k];
        Eigen::Matrix<T,N,1>& A5k = A5[k];
        Eigen::Matrix<T,N,1>& B5k = B5[k];
        for(int i=0;i<N;++i) A5k(i) = A1k(i);
        for(int i=0;i<N;++i) B5k(i) = B1k(i);
    }
#endif

    timeval tp;

    double t1_reg=0., t1_small=0., t1_blas=0., t1_eigen=0., t1_smalleigen=0.;
    double t2_reg=0., t2_small=0., t2_blas=0., t2_eigen=0., t2_smalleigen=0.;
    double t3_reg=0., t3_small=0., t3_blas=0., t3_eigen=0., t3_smalleigen=0.;
    double t4_reg=0., t4_small=0., t4_blas=0., t4_eigen=0., t4_smalleigen=0.;
    double t5_reg=0., t5_small=0., t5_blas=0., t5_eigen=0., t5_smalleigen=0.;
    double t6_reg=0., t6_small=0., t6_blas=0., t6_eigen=0., t6_smalleigen=0.;
    double t7_reg=0., t7_small=0., t7_blas=0., t7_eigen=0., t7_smalleigen=0.;
    double t8_reg=0., t8_small=0., t8_blas=0., t8_eigen=0., t8_smalleigen=0.;
    double ta,tb;

#ifdef ERRORCHECK
    double e1_reg=0., e1_small=0., e1_blas=0., e1_eigen=0., e1_smalleigen=0.;
    double e2_reg=0., e2_small=0., e2_blas=0., e2_eigen=0., e2_smalleigen=0.;
    double e3_reg=0., e3_small=0., e3_blas=0., e3_eigen=0., e3_smalleigen=0.;
    double e4_reg=0., e4_small=0., e4_blas=0., e4_eigen=0., e4_smalleigen=0.;
    double e5_reg=0., e5_small=0., e5_blas=0., e5_eigen=0., e5_smalleigen=0.;
    double e6_reg=0., e6_small=0., e6_blas=0., e6_eigen=0., e6_smalleigen=0.;
    double e7_reg=0., e7_small=0., e7_blas=0., e7_eigen=0., e7_smalleigen=0.;
    double e8_reg=0., e8_small=0., e8_blas=0., e8_eigen=0., e8_smalleigen=0.;
#endif

    for (int i=0; i<nloops1; ++i) {

#if 1 // Swap(A,B);
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& B0k = B0[k];
                for(int i=0;i<N;++i) std::swap(A0k(i),B0k(i));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) 
            Swap(A1[k],B1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(B1[k]-B0[k]);
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(A1[k]-A0[k]);
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            Swap(A2[k],B2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(B2[k]-B0[k]);
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(A2[k]-A0[k]);
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(swap) (N,BP(A3[k].ptr()),1,BP(B3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(B3[k]-B0[k]);
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(A3[k]-A0[k]);
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            A4[k].swap(B4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e1_eigen += tmv::TMV_NORM(A4[k](i)-A0[k](i));
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e1_eigen += tmv::TMV_NORM(B4[k](i)-B0[k](i));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            A5[k].swap(B5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e1_smalleigen += tmv::TMV_NORM(A5[k](i)-A0[k](i));
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e1_smalleigen += tmv::TMV_NORM(B5[k](i)-B0[k](i));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // Swap(A[0::N/2],B[0::N/2])
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& B0k = B0[k];
                for(int i=0;i<N/2;++i) std::swap(A0k(i),B0k(i));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) 
            Swap(A1[k].subVector(0,N/2),B1[k].subVector(0,N/2));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) e2_reg += NormSq(A1[k]-A0[k]);
            for (int k=0; k<nloops2; ++k) e2_reg += NormSq(B1[k]-B0[k]);
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            Swap(A2[k].subVector(0,N/2),B2[k].subVector(0,N/2));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) e2_small += NormSq(A2[k]-A0[k]);
            for (int k=0; k<nloops2; ++k) e2_small += NormSq(B2[k]-B0[k]);
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(swap) (N/2,BP(A3[k].ptr()),1,BP(B3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) e2_blas += NormSq(A3[k]-A0[k]);
            for (int k=0; k<nloops2; ++k) e2_blas += NormSq(B3[k]-B0[k]);
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            A4[k].segment(0,N/2).swap(B4[k].segment(0,N/2));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e2_eigen += tmv::TMV_NORM(A4[k](i)-A0[k](i));
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e2_eigen += tmv::TMV_NORM(B4[k](i)-B0[k](i));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        // Eigen has a bug where it gets this wrong for N == 1
        if (N > 1) 
            for (int k=0; k<nloops2x; ++k)
                A5[k].segment(0,N/2).swap(B5[k].segment(0,N/2));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e2_smalleigen += tmv::TMV_NORM(A5[k](i)-A0[k](i));
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e2_smalleigen += tmv::TMV_NORM(B5[k](i)-B0[k](i));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // A.reverseSelf()
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                tmv::Vector<T>& A0k = A0[k];
                for(int i=0;i<N/2;++i) std::swap(A0k(i),A0k(N-i-1));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            A1[k].reverseSelf();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) e3_reg += NormSq(A1[k]-A0[k]);
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            A2[k].reverseSelf();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) e3_small += NormSq(A2[k]-A0[k]);
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(swap) (N/2,BP(A3[k].ptr()),1,BP(A3[k].ptr()+N-N/2),-1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) e3_blas += NormSq(A3[k]-A0[k]);
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // I couldn't find this capability in Eigen...
            for(int i=0;i<N/2;i++) std::swap(A4[k][i],A4[k][N-i-1]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e3_eigen += tmv::TMV_NORM(A4[k](i)-A0[k](i));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            for(int i=0;i<N/2;i++) std::swap(A5[k][i],A5[k][N-i-1]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e3_smalleigen += tmv::TMV_NORM(A5[k](i)-A0[k](i));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // A.sort()
#ifndef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                tmv::Vector<T>& A0k = A0[k];
                std::sort(A0k.ptr(),A0k.ptr()+N);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            A1[k].sort();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) e4_reg += NormSq(A1[k]-A0[k]);
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            A2[k].sort();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) e4_small += NormSq(A2[k]-A0[k]);
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::sort(A3[k].begin(),A3[k].end());
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) e4_blas += NormSq(A3[k]-A0[k]);
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            std::sort(A4[k].data(),A4[k].data()+N);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e4_eigen += tmv::TMV_NORM(A4[k](i)-A0[k](i));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            std::sort(A5[k].data(),A5[k].data()+N);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e4_smalleigen += tmv::TMV_NORM(A5[k](i)-A0[k](i));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // A.conjugateSelf()
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                tmv::Vector<T>& A0k = A0[k];
                for(int i=0;i<N;++i) A0k(i) = std::conj(A0k(i));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            A1[k].conjugateSelf();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) e5_reg += NormSq(A1[k]-A0[k]);
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            A2[k].conjugateSelf();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) e5_small += NormSq(A2[k]-A0[k]);
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASDNAME(scal) (N,dmone,A3[k].realPart().ptr()+1,2);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) e5_blas += NormSq(A3[k]-A0[k]);
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            A4[k] = A4[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e5_eigen += tmv::TMV_NORM(A4[k](i)-A0[k](i));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            A5[k] = A5[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e5_smalleigen += tmv::TMV_NORM(A5[k](i)-A0[k](i));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // Swap(A,B.conjugate);
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& B0k = B0[k];
                for(int i=0;i<N;++i) {
                    T temp = conj(B0k(i));
                    B0k(i) = conj(A0k(i));
                    A0k(i) = temp;
                }
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) 
            Swap(A1[k],B1[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) e6_reg += NormSq(B1[k]-B0[k]);
            for (int k=0; k<nloops2; ++k) e6_reg += NormSq(A1[k]-A0[k]);
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            Swap(A2[k],B2[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) e6_small += NormSq(B2[k]-B0[k]);
            for (int k=0; k<nloops2; ++k) e6_small += NormSq(A2[k]-A0[k]);
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(swap) (N,BP(A3[k].ptr()),1,BP(B3[k].ptr()),1);
            BLASDNAME(scal) (N,dmone,A3[k].realPart().ptr()+1,2);
            BLASDNAME(scal) (N,dmone,B3[k].realPart().ptr()+1,2);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) e6_blas += NormSq(B3[k]-B0[k]);
            for (int k=0; k<nloops2; ++k) e6_blas += NormSq(A3[k]-A0[k]);
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // Neither of these work.
            //A4[k].swap(B4[k].conjugate());
            //A4[k].conjugate().swap(B4[k]);
            A4[k].swap(B4[k]);
            A4[k] = A4[k].conjugate();
            B4[k] = B4[k].conjugate();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e6_eigen += tmv::TMV_NORM(A4[k](i)-A0[k](i));
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e6_eigen += tmv::TMV_NORM(B4[k](i)-B0[k](i));
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            A5[k].swap(B5[k]);
            A5[k] = A5[k].conjugate();
            B5[k] = B5[k].conjugate();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e6_smalleigen += tmv::TMV_NORM(A5[k](i)-A0[k](i));
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e6_smalleigen += tmv::TMV_NORM(B5[k](i)-B0[k](i));
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // A.permute(P);
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                tmv::Vector<T>& A0k = A0[k];
                for(int j=0;j<N;++j) std::swap(A0k(j),A0k(P[j]));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            A1[k].permute(P);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) e7_reg += NormSq(A1[k]-A0[k]);
            e7_reg = sqrt(e7_reg/nloops2);
            e7_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            A2[k].permute(P);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) e7_small += NormSq(A2[k]-A0[k]);
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // Blas doesn't have this. 
            for(int i=0;i<N;++i) std::swap(A3[k](i),A3[k](P[i]));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) e7_blas += NormSq(A3[k]-A0[k]);
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // I couldn't find this capability in Eigen.
            // When they do this as part of LU, they just write out the loop.
            for(int i=0;i<N;++i) std::swap(A4[k](i),A4[k](P[i]));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e7_eigen += tmv::TMV_NORM(A4[k](i)-A0[k](i));
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            for(int i=0;i<N;++i) std::swap(A5[k](i),A5[k](P[i]));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e7_smalleigen += tmv::TMV_NORM(A5[k](i)-A0[k](i));
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // A.reversePermute(P);
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                tmv::Vector<T>& A0k = A0[k];
                for(int j=N-1;j>=0;--j) std::swap(A0k(j),A0k(P[j]));
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            A1[k].reversePermute(P);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) e8_reg += NormSq(A1[k]-A0[k]);
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            A2[k].reversePermute(P);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) e8_small += NormSq(A2[k]-A0[k]);
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // Blas doesn't have this. 
            for(int i=N-1;i>=0;--i) std::swap(A3[k](i),A3[k](P[i]));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) e8_blas += NormSq(A3[k]-A0[k]);
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // I couldn't find this capability in Eigen.
            // When they do this as part of LU, they just write out the loop.
            for(int i=N-1;i>=0;--i) std::swap(A4[k](i),A4[k](P[i]));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e8_eigen += tmv::TMV_NORM(A4[k](i)-A0[k](i));
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            for(int i=N-1;i>=0;--i) std::swap(A5[k](i),A5[k](P[i]));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e8_smalleigen += tmv::TMV_NORM(A5[k](i)-A0[k](i));
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif
    }

    std::cout<<"Swap(A,B)               "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"Swap(A(0:N/2),B(0:N/2)) "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"A.reverseSelf()         "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
#ifndef TISCOMPLEX
    std::cout<<"A.sort()                "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
#endif
#ifdef TISCOMPLEX
    std::cout<<"A.conjugateSelf()       "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
    std::cout<<"Swap(A,B.conjugate())   "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
#endif
    std::cout<<"A.permute(P)            "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"A.reversePermute(P)     "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"Swap(A,B)               "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"Swap(A(0:N/2),B(0:N/2)) "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"A.reverseSelf()         "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
#ifndef TISCOMPLEX
    std::cout<<"A.sort()                "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
#endif
#ifdef TISCOMPLEX
    std::cout<<"A.conjugateSelf()       "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
    std::cout<<"Swap(A,B.conjugate())   "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
#endif
    std::cout<<"A.permute(P)            "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
    std::cout<<"A.reversePermute(P)     "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
    std::cout<<"\n\n";
#endif
}
#endif

#ifdef DOELEMMULT
static void ElemMultVV(
    const std::vector<tmv::Vector<T> >& A1,
    const std::vector<tmv::Vector<T> >& B1,
    std::vector<tmv::Vector<T> >& C1)
{
#ifdef ERRORCHECK
    std::vector<tmv::Vector<T> > A0 = A1;
    std::vector<tmv::Vector<T> > B0 = B1;
    std::vector<tmv::Vector<T> > C0 = C1;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallVector<T,N> > A2(nloops2);
    std::vector<tmv::SmallVector<T,N> > B2(nloops2);
    std::vector<tmv::SmallVector<T,N> > C2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; B2[k] = B1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Vector<T> > A3 = A1;
    std::vector<tmv::Vector<T> > B3 = B1;
    std::vector<tmv::Vector<T> > C3 = C1;
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENV> A4(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENV> B4(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENV> C4(nloops2,Eigen::EIGENV(N));
    for(int k=0;k<nloops2;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        const tmv::Vector<T>& B1k = B1[k];
        Eigen::EIGENV& A4k = A4[k];
        Eigen::EIGENV& B4k = B4[k];
        for(int i=0;i<N;++i) A4k(i) = A1k(i);
        for(int i=0;i<N;++i) B4k(i) = B1k(i);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,N,1> > A5;
    std::vector<Eigen::Matrix<T,N,1> > B5;
    std::vector<Eigen::Matrix<T,N,1> > C5;
    if (nloops2x) 
    { A5.resize(nloops2x); B5.resize(nloops2x); C5.resize(nloops2x); }
    for(int k=0;k<nloops2x;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        const tmv::Vector<T>& B1k = B1[k];
        Eigen::Matrix<T,N,1>& A5k = A5[k];
        Eigen::Matrix<T,N,1>& B5k = B5[k];
        for(int i=0;i<N;++i) A5k(i) = A1k(i);
        for(int i=0;i<N;++i) B5k(i) = B1k(i);
    }
#endif

    timeval tp;

    double t1_reg=0., t1_small=0., t1_blas=0., t1_eigen=0., t1_smalleigen=0.;
    double t2_reg=0., t2_small=0., t2_blas=0., t2_eigen=0., t2_smalleigen=0.;
    double t3_reg=0., t3_small=0., t3_blas=0., t3_eigen=0., t3_smalleigen=0.;
    double t4_reg=0., t4_small=0., t4_blas=0., t4_eigen=0., t4_smalleigen=0.;
    double t5_reg=0., t5_small=0., t5_blas=0., t5_eigen=0., t5_smalleigen=0.;
    double t6_reg=0., t6_small=0., t6_blas=0., t6_eigen=0., t6_smalleigen=0.;
    double t7_reg=0., t7_small=0., t7_blas=0., t7_eigen=0., t7_smalleigen=0.;
    double t8_reg=0., t8_small=0., t8_blas=0., t8_eigen=0., t8_smalleigen=0.;
    double t9_reg=0., t9_small=0., t9_blas=0., t9_eigen=0., t9_smalleigen=0.;
    double t10_reg=0., t10_small=0., t10_blas=0., t10_eigen=0., t10_smalleigen=0.;
    double ta,tb;

#ifdef ERRORCHECK
    double e1_reg=0., e1_small=0., e1_blas=0., e1_eigen=0., e1_smalleigen=0.;
    double e2_reg=0., e2_small=0., e2_blas=0., e2_eigen=0., e2_smalleigen=0.;
    double e3_reg=0., e3_small=0., e3_blas=0., e3_eigen=0., e3_smalleigen=0.;
    double e4_reg=0., e4_small=0., e4_blas=0., e4_eigen=0., e4_smalleigen=0.;
    double e5_reg=0., e5_small=0., e5_blas=0., e5_eigen=0., e5_smalleigen=0.;
    double e6_reg=0., e6_small=0., e6_blas=0., e6_eigen=0., e6_smalleigen=0.;
    double e7_reg=0., e7_small=0., e7_blas=0., e7_eigen=0., e7_smalleigen=0.;
    double e8_reg=0., e8_small=0., e8_blas=0., e8_eigen=0., e8_smalleigen=0.;
    double e9_reg=0., e9_small=0., e9_blas=0., e9_eigen=0., e9_smalleigen=0.;
    double e10_reg=0., e10_small=0., e10_blas=0., e10_eigen=0., e10_smalleigen=0.;
#endif

    for (int i=0; i<nloops1; ++i) {

#if 1 // C = DiagMatrixViewOf(A) * B
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = A0k(i) * B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = DiagMatrixViewOf(A1[k]) * B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(C1[k]-C0[k]);
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = DiagMatrixViewOf(A2[k]) * B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(C2[k]-C0[k]);
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(A3[k].begin(),A3[k].end(),B3[k].begin(),C3[k].begin(),
                           std::multiplies<T>());
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(C3[k]-C0[k]);
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            // This version does ok speed-wise
            C4[k] = A4[k].asDiagonal() * B4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e1_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] = A5[k].asDiagonal() * B5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e1_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C *= 8 * DiagMatrixViewOf(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) *= RT(8) * A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] *= 8 * DiagMatrixViewOf(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) e2_reg += NormSq(C1[k]-C0[k]);
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] *= 8 * DiagMatrixViewOf(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) e2_small += NormSq(C2[k]-C0[k]);
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(A3[k].begin(),A3[k].end(),C3[k].begin(),C3[k].begin(),
                           std::multiplies<T>());
            BLASNAME(scal) (N,eight,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) e2_blas += NormSq(C3[k]-C0[k]);
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // This version is very slow, so unfair comparison when optimizing.
            //C4[k] = 8 * A4[k].asDiagonal() * C4[k];
            // This version doesn't compile.
            //C4[k] = 8 * A4[k].cwise() * C4[k];
            // The following double statement is the only version I've found
            // that is efficient and compiles.
            C4[k] = A4[k].cwise() * C4[k];
            C4[k] *= 8;
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e2_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            C5[k] = A5[k].cwise() * C5[k];
            C5[k] *= 8;
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e2_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // C = 7 * DiagMatrixViewOf(A) * B
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = RT(7) * A0k(i) * B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = 7 * DiagMatrixViewOf(A1[k]) * B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) e3_reg += NormSq(C1[k]-C0[k]);
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = 7 * DiagMatrixViewOf(A2[k]) * B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) e3_small += NormSq(C2[k]-C0[k]);
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(A3[k].begin(),A3[k].end(),B3[k].begin(),C3[k].begin(),
                           std::multiplies<T>());
            BLASNAME(scal) (N,seven,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) e3_blas += NormSq(C3[k]-C0[k]);
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            //C4[k] = 7 * A4[k].asDiagonal() * B4[k];
            C4[k] = A4[k].cwise() * B4[k];
            C4[k] *= 7;
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e3_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            C5[k] = A5[k].cwise() * B5[k];
            C5[k] *= 7;
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e3_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif
#endif

#if 1 // C *= (8,9) * DiagMatrixViewOf(A)
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) *= T(8,9) * A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] *= T(8,9) * DiagMatrixViewOf(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) e4_reg += NormSq(C1[k]-C0[k]);
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] *= T(8,9) * DiagMatrixViewOf(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) e4_small += NormSq(C2[k]-C0[k]);
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(A3[k].begin(),A3[k].end(),C3[k].begin(),C3[k].begin(),
                           std::multiplies<T>());
            BLASNAME(scal) (N,z89,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) e4_blas += NormSq(C3[k]-C0[k]);
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            C4[k] = A4[k].cwise() * C4[k];
            C4[k] *= T(8,9);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e4_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            C5[k] = A5[k].cwise() * C5[k];
            C5[k] *= T(8,9);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e4_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // C = (7,1) * DiagMatrixViewOf(A) * B
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = T(7,1) * A0k(i) * B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = T(7,1) * DiagMatrixViewOf(A1[k]) * B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) e5_reg += NormSq(C1[k]-C0[k]);
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = T(7,1) * DiagMatrixViewOf(A2[k]) * B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) e5_small += NormSq(C2[k]-C0[k]);
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(A3[k].begin(),A3[k].end(),B3[k].begin(),C3[k].begin(),
                           std::multiplies<T>());
            BLASNAME(scal) (N,z71,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) e5_blas += NormSq(C3[k]-C0[k]);
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            //C4[k] = T(7,1) * A4[k].asDiagonal() * B4[k];
            C4[k] = A4[k].cwise() * B4[k];
            C4[k] *= T(7,1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e5_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            C5[k] = A5[k].cwise() * B5[k];
            C5[k] *= T(7,1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e5_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif
#endif
#endif
    }

    std::cout<<"C = diag(A) * B         "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"C *= 8*diag(A)          "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"C = 7*diag(A) * B       "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C *= (8,9)*diag(A)      "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"C = (7,1)*diag(A) * B   "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"C = diag(A) * B         "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"C *= 8*diag(A)          "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"C = 7*diag(A) * B       "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C *= (8,9)*diag(A)      "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"C = (7,1)*diag(A) * B   "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
#endif
    std::cout<<"\n\n";
#endif
}
#endif

#ifdef DOELEMDIV
static void ElemDivVV(
    const std::vector<tmv::Vector<T> >& A1,
    const std::vector<tmv::Vector<T> >& B1,
    std::vector<tmv::Vector<T> >& C1)
{
#ifdef ERRORCHECK
    std::vector<tmv::Vector<T> > A0 = A1;
    std::vector<tmv::Vector<T> > B0 = B1;
    std::vector<tmv::Vector<T> > C0 = C1;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallVector<T,N> > A2(nloops2);
    std::vector<tmv::SmallVector<T,N> > B2(nloops2);
    std::vector<tmv::SmallVector<T,N> > C2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; B2[k] = B1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Vector<T> > A3 = A1;
    std::vector<tmv::Vector<T> > B3 = B1;
    std::vector<tmv::Vector<T> > C3 = C1;
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENV> A4(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENV> B4(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENV> C4(nloops2,Eigen::EIGENV(N));
    for(int k=0;k<nloops2;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        const tmv::Vector<T>& B1k = B1[k];
        Eigen::EIGENV& A4k = A4[k];
        Eigen::EIGENV& B4k = B4[k];
        for(int i=0;i<N;++i) A4k(i) = A1k(i);
        for(int i=0;i<N;++i) B4k(i) = B1k(i);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,N,1> > A5;
    std::vector<Eigen::Matrix<T,N,1> > B5;
    std::vector<Eigen::Matrix<T,N,1> > C5;
    if (nloops2x) 
    { A5.resize(nloops2x); B5.resize(nloops2x); C5.resize(nloops2x); }
    for(int k=0;k<nloops2x;++k) {
        const tmv::Vector<T>& A1k = A1[k];
        const tmv::Vector<T>& B1k = B1[k];
        Eigen::Matrix<T,N,1>& A5k = A5[k];
        Eigen::Matrix<T,N,1>& B5k = B5[k];
        for(int i=0;i<N;++i) A5k(i) = A1k(i);
        for(int i=0;i<N;++i) B5k(i) = B1k(i);
    }
#endif

    timeval tp;

    double t1_reg=0., t1_small=0., t1_blas=0., t1_eigen=0., t1_smalleigen=0.;
    double t2_reg=0., t2_small=0., t2_blas=0., t2_eigen=0., t2_smalleigen=0.;
    double t3_reg=0., t3_small=0., t3_blas=0., t3_eigen=0., t3_smalleigen=0.;
    double t4_reg=0., t4_small=0., t4_blas=0., t4_eigen=0., t4_smalleigen=0.;
    double t5_reg=0., t5_small=0., t5_blas=0., t5_eigen=0., t5_smalleigen=0.;
    double t6_reg=0., t6_small=0., t6_blas=0., t6_eigen=0., t6_smalleigen=0.;
    double t7_reg=0., t7_small=0., t7_blas=0., t7_eigen=0., t7_smalleigen=0.;
    double t8_reg=0., t8_small=0., t8_blas=0., t8_eigen=0., t8_smalleigen=0.;
    double t9_reg=0., t9_small=0., t9_blas=0., t9_eigen=0., t9_smalleigen=0.;
    double t10_reg=0., t10_small=0., t10_blas=0., t10_eigen=0., t10_smalleigen=0.;
    double ta,tb;

#ifdef ERRORCHECK
    double e1_reg=0., e1_small=0., e1_blas=0., e1_eigen=0., e1_smalleigen=0.;
    double e2_reg=0., e2_small=0., e2_blas=0., e2_eigen=0., e2_smalleigen=0.;
    double e3_reg=0., e3_small=0., e3_blas=0., e3_eigen=0., e3_smalleigen=0.;
    double e4_reg=0., e4_small=0., e4_blas=0., e4_eigen=0., e4_smalleigen=0.;
    double e5_reg=0., e5_small=0., e5_blas=0., e5_eigen=0., e5_smalleigen=0.;
    double e6_reg=0., e6_small=0., e6_blas=0., e6_eigen=0., e6_smalleigen=0.;
    double e7_reg=0., e7_small=0., e7_blas=0., e7_eigen=0., e7_smalleigen=0.;
    double e8_reg=0., e8_small=0., e8_blas=0., e8_eigen=0., e8_smalleigen=0.;
    double e9_reg=0., e9_small=0., e9_blas=0., e9_eigen=0., e9_smalleigen=0.;
    double e10_reg=0., e10_small=0., e10_blas=0., e10_eigen=0., e10_smalleigen=0.;
#endif

    for (int i=0; i<nloops1; ++i) {

#if 1 // C = B / DiagMatrixViewOf(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = B0k(i) / A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = B1[k] / DiagMatrixViewOf(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(C1[k]-C0[k]);
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = B2[k] / DiagMatrixViewOf(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(C2[k]-C0[k]);
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(
                B3[k].begin(),B3[k].end(),A3[k].begin(),C3[k].begin(),
                std::divides<T>());
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(C3[k]-C0[k]);
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] = B4[k].cwise() / A4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) 
                    e1_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] = B5[k].cwise() / A5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) 
                    e1_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C /= 8 * DiagMatrixViewOf(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) /= RT(8) * A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] /= 8 * DiagMatrixViewOf(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) e2_reg += NormSq(C1[k]-C0[k]);
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] /= 8 * DiagMatrixViewOf(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) e2_small += NormSq(C2[k]-C0[k]);
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(
                C3[k].begin(),C3[k].end(),A3[k].begin(),C3[k].begin(),
                std::divides<T>());
            BLASNAME(scal) (N,oneeighth,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) e2_blas += NormSq(C3[k]-C0[k]);
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            C4[k].cwise() /= 8 * A4[k];
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) e2_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            C5[k].cwise() /= 8 * A5[k];
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) e2_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // C = 7 * DiagMatrixViewOf(A).inverse() * B
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = RT(7) / A0k(i) * B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = 7 * DiagMatrixViewOf(A1[k]).inverse() * B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) e3_reg += NormSq(C1[k]-C0[k]);
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = 7 * DiagMatrixViewOf(A2[k]).inverse() * B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) e3_small += NormSq(C2[k]-C0[k]);
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(
                B3[k].begin(),B3[k].end(),A3[k].begin(),C3[k].begin(),
                std::divides<T>());
            BLASNAME(scal) (N,seven,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) e3_blas += NormSq(C3[k]-C0[k]);
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            C4[k] = B4[k].cwise() / A4[k];
            C4[k] *= 7;
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) 
                    e3_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            C5[k] = B5[k].cwise() / A5[k];
            C5[k] *= 7;
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) 
                    e3_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif
#endif

#if 1 // C = (8,9) * C % DiagMatrixViewOf(A)
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = T(8,9) * C0k(i) / A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = T(8,9) * C1[k] % DiagMatrixViewOf(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) e4_reg += NormSq(C1[k]-C0[k]);
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = T(8,9) * C2[k] % DiagMatrixViewOf(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) e4_small += NormSq(C2[k]-C0[k]);
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(
                C3[k].begin(),C3[k].end(),A3[k].begin(),C3[k].begin(),
                std::divides<T>());
            BLASNAME(scal) (N,z89,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) e4_blas += NormSq(C3[k]-C0[k]);
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            C4[k].cwise() /= A4[k];
            C4[k] *= T(8,9);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) 
                    e4_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            C5[k].cwise() /= A5[k];
            C5[k] *= T(8,9);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) 
                    e4_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // C = (7,1) * DiagMatrixViewOf(A).inverse() * B
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = T(7,1) / A0k(i) * B0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C1[k] = T(7,1) / DiagMatrixViewOf(A1[k]) * B1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) e5_reg += NormSq(C1[k]-C0[k]);
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C2[k] = T(7,1) / DiagMatrixViewOf(A2[k]) * B2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) e5_small += NormSq(C2[k]-C0[k]);
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(
                B3[k].begin(),B3[k].end(),A3[k].begin(),C3[k].begin(),
                std::divides<T>());
            BLASNAME(scal) (N,z71,BP(C3[k].ptr()),1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) e5_blas += NormSq(C3[k]-C0[k]);
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            C4[k] = B4[k].cwise() / A4[k];
            C4[k] *= T(7,1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) 
                    e5_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            C5[k] = B5[k].cwise() / A5[k];
            C5[k] *= T(7,1);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) 
                    e5_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= 7.*(Norm(A0[0]) + Norm(B0[0]));
        }
#endif
#endif
#endif
#endif

#if 1 // C = DiagMatrixViewOf(A).inverse()
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = T(1) / A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            DiagMatrixViewOf(C1[k]) = DiagMatrixViewOf(A1[k]).inverse();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) e6_reg += NormSq(C1[k]-C0[k]);
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            DiagMatrixViewOf(C2[k]) = DiagMatrixViewOf(A2[k]).inverse();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) e6_small += NormSq(C2[k]-C0[k]);
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(
                A3[k].begin(),A3[k].end(),C3[k].begin(),
                bind1st(std::divides<T>(),T(1)));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) e6_blas += NormSq(C3[k]-C0[k]);
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] = A4[k].cwise().inverse();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) 
                    e6_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] = A5[k].cwise().inverse();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) 
                    e6_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = 7 / DiagMatrixViewOf(A)
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = T(7) / A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            DiagMatrixViewOf(C1[k]) = 7 / DiagMatrixViewOf(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) e7_reg += NormSq(C1[k]-C0[k]);
            e7_reg = sqrt(e7_reg/nloops2);
            e7_reg /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            DiagMatrixViewOf(C2[k]) = 7 / DiagMatrixViewOf(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) e7_small += NormSq(C2[k]-C0[k]);
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(
                A3[k].begin(),A3[k].end(),C3[k].begin(),
                bind1st(std::divides<T>(),T(7)));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) e7_blas += NormSq(C3[k]-C0[k]);
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] = A4[k].cwise().inverse() * 7;

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) 
                    e7_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] = A5[k].cwise().inverse() * 7;

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) 
                    e7_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = (7,1) / DiagMatrixViewOf(A)
#ifdef TISCOMPLEX
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = T(7,1) / A0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            DiagMatrixViewOf(C1[k]) = T(7,1) / DiagMatrixViewOf(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) e8_reg += NormSq(C1[k]-C0[k]);
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            DiagMatrixViewOf(C2[k]) = T(7,1) / DiagMatrixViewOf(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) e8_small += NormSq(C2[k]-C0[k]);
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(
                A3[k].begin(),A3[k].end(),C3[k].begin(),
                bind1st(std::divides<T>(),T(7,1)));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) e8_blas += NormSq(C3[k]-C0[k]);
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] = A4[k].cwise().inverse() * T(7,1);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) 
                    e8_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] = A5[k].cwise().inverse() * T(7,1);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) 
                    e8_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // C = DiagMatrixViewOf(C).inverse()
#ifdef ERRORCHECK
        if (i == 0) {
            for (int k=0; k<nloops2; ++k) {
                const tmv::Vector<T>& A0k = A0[k];
                const tmv::Vector<T>& B0k = B0[k];
                tmv::Vector<T>& C0k = C0[k];
                for(int i=0;i<N;++i) C0k(i) = T(1) / C0k(i);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            DiagMatrixViewOf(C1[k]) = DiagMatrixViewOf(C1[k]).inverse();
        

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_reg = 0.;
            for (int k=0; k<nloops2; ++k) e9_reg += NormSq(C1[k]-C0[k]);
            e9_reg = sqrt(e9_reg/nloops2);
            e9_reg /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            DiagMatrixViewOf(C2[k]) = DiagMatrixViewOf(C2[k]).inverse();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_small = 0.;
            for (int k=0; k<nloops2; ++k) e9_small += NormSq(C2[k]-C0[k]);
            e9_small = sqrt(e9_small/nloops2);
            e9_small /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            std::transform(
                C3[k].begin(),C3[k].end(),C3[k].begin(),
                bind1st(std::divides<T>(),T(1)));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_blas += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_blas = 0.;
            for (int k=0; k<nloops2; ++k) e9_blas += NormSq(C3[k]-C0[k]);
            e9_blas = sqrt(e9_blas/nloops2);
            e9_blas /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            C4[k] = C4[k].cwise().inverse();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) 
                    e9_eigen += tmv::TMV_NORM(C4[k](i)-C0[k](i));
            e9_eigen = sqrt(e9_eigen/nloops2);
            e9_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            C5[k] = C5[k].cwise().inverse();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (i == 0) {
            e9_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<N;++i) 
                    e9_smalleigen += tmv::TMV_NORM(C5[k](i)-C0[k](i));
            e9_smalleigen = sqrt(e9_smalleigen/nloops2);
            e9_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif
    }

    std::cout<<"C = B / diag(A)         "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"C /= 8 * diag(A)        "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"C = 7*diag(A)^-1 * B    "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C = (8,9)*C \% diag(A)  "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"C = (7,1)*diag(A)^-1*B  "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
#endif
    std::cout<<"C = diag(A)^-1          "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"C = 7 / diag(A)         "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C = (7,1) / diag(A)     "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
#endif
    std::cout<<"C = diag(C)^-1          "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
    std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"C = B / diag(A)         "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"C /= 8 * diag(A)        "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"C = 7*diag(A)^-1 * B    "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C = (8,9)*C \% diag(A)  "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"C = (7,1)*diag(A)^-1*B  "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
#endif
    std::cout<<"C = diag(A)^-1          "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
    std::cout<<"C = 7 / diag(A)         "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C = (7,1) / diag(A)     "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
#endif
    std::cout<<"C = diag(C)^-1          "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
    std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
    std::cout<<"\n\n";
#endif
}
#endif

#ifdef TISCOMPLEX
#define RAND ( T(rand(),rand()) / RT(RAND_MAX) )
#else
#define RAND ( T(rand()) / RT(RAND_MAX) )
#endif

int main() try 
{
    srand(518423972);

#ifdef MEMDEBUG
    atexit(&DumpUnfreed);
#endif

    std::vector<tmv::Vector<T> > A(nloops2,tmv::Vector<T>(N));
    std::vector<tmv::Vector<T> > B(nloops2,tmv::Vector<T>(N));
    std::vector<tmv::Vector<T> > C(nloops2,tmv::Vector<T>(N));
    for(int k=0;k<nloops2;++k) {
        for(int i=0;i<N;++i) {
            A[k](i) = RAND * RT(1000) + RAND - RT(512);
            B[k](i) = RAND * RT(1000) + RAND - RT(512);
        }
    }
    std::cout<<"N = "<<N<<std::endl;
    //std::cout<<"A[0] = "<<A[0]<<std::endl;
    std::cout<<"nloops = "<<nloops1<<" x "<<nloops2;
    std::cout<<" = "<<nloops1*nloops2<<std::endl;
#ifdef DOEIGENSMALL
    if (nloops2x == 0)
        std::cout<<"Vector is too big for the stack, so no \"Eigen Known\" tests.\n";
#endif
    std::cout<<"\nTime for:               ";
    std::cout<<"  TMV    "<<" TMV Small"<<"   BLAS   "<<"  Eigen  "<<"Eigen Known"<<std::endl;

#ifdef DOMULTXV
    MultXV(A,B);
#endif

    for(int k=0;k<nloops2;++k) {
        for(int i=0;i<N;++i) {
            B[k](i) = RAND * RT(1000) + RAND - RT(512);
        }
    }


#ifdef DOADDVV
    AddVV(A,B,C);
#endif

#ifdef DOMULTVV
    MultVV(A,B);
#endif

#ifdef DOMINMAX
    MinMaxV(A);
#endif

#ifdef DONORM
    NormV(A);
#endif

#ifdef DOELEMMULT
    ElemMultVV(A,B,C);
#endif

#ifdef DOELEMDIV
    ElemDivVV(A,B,C);
#endif

    // Do this last, since it screws around with A,B.
#ifdef DOSWAP
    SwapV(A,B);
#endif

    return 0;
}
#if 1
catch (int) {}
#else
catch (tmv::Error& e) {
        std::cout<<"Caught error "<<e<<std::endl;
            return 1;
}
#endif

