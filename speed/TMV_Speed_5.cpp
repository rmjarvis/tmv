//#define PRINTALGO_QR
//#define PRINTALGO_DIVM
//#define PRINTALGO_DIVU
//#define PRINTALGO_DIVU_OMP
//#define PRINTALGO_INVU
//#define PRINTALGO_INVM
//#define PRINTALGO_UL
//#define PRINTALGO_XV
//#define PRINTALGO_MM
//#define PRINTALGO_UM
//#define PRINTALGO_MV
//#define PRINTALGO_DET

//#define XDEBUG_PRODMM
//#define XDEBUG_QUOTMM
//#define XDEBUG_HOUSE

//#undef NDEBUG
#include <iostream>

#define TMV_NO_LIB
#include "TMV.h"

// How big do you want the matrices to be?
// (The matrix being inverted is MxN.  
//  K is the other dimension of the matrix being solved.)
const int M = 8;
const int N = 3;
const int K = 5;
//#define AISSQUARE

// Define the type to use:
//#define TISFLOAT
//#define TISCOMPLEX

// Define the parts of the matrices to use
// 1 = all (rectangle matrix)
#define PART1 1

// Which algorithm to use?
// 1 = QR
// 2 = QRP
// 3 = SV
#define ALGO 2
//#define STRICTQRP

// Define whether you want to include the error checks.
//#define ERRORCHECK

// Define which versions you want to test:
#define DOREG
#define DOSMALL
#define DOLAP
#define DOEIGEN
#define DOEIGENSMALL

// Skip the separate decomposition and solution steps?
//#define NO_SEPARATE_DECOMP

// For large N, this can be useful to keep the condition of the matrix
// not too large.
//#define SMALL_OFFDIAG

// Define which batches of functions you want to test:
#define DO_C
//#define DO_R  

// Set this if you only want to do a single loop.
// Not so useful for timing, but useful for debugging.
// Also useful if M,N,K are very large to avoid overflow in product
//#define ONELOOP

// Normally I fill the matrices and vectors with random values, but 
// that can make it difficult to debug the arithmetic.
// This define uses simpler values that I can multiply in my head.
// (If doing det, should also set SMALL_OFFDIAG so det != 0)
//#define SIMPLE_VALUES

// Show progress with dots....
//#define SHOW_DOTS

// Set up the target number of operations and memory to use for testing
const long long targetnmegaflops(1000); // in 10^6 real ops
const long long targetmem(10000); // in Kbytes

// Include the LAPACK library you want to test TMV against:
#if 0
#include "mkl.h"
#define MKL
#define CBLAS
#define CLAPACK
#define XLAP
#elif 0
extern "C" {
#include "cblas.h"
#include "util/flapack.h"
}
#define CBLAS
#define FLAP
#define XLAP
#else
extern "C" {
#include "util/fblas.h"
#include "util/flapack.h"
}
#define FBLAS
#define FLAP
#define XLAP
#endif


// The rest of this is based on the above values.
// You shouldn't need to change anything else here.

#include <complex>

#ifdef TISFLOAT
typedef float RT;
#else
typedef double RT;
#endif

#if ALGO == 1
#define GETDIV(A) A.qrd()
#define TMV_QR tmv::QR
#define TMV_QRD tmv::QRD
#elif ALGO == 2
#define GETDIV(A) A.qrpd()
#define TMV_QR tmv::QRP
#define TMV_QRD tmv::QRPD
#elif ALGO == 3
#define GETDIV(A) A.svd()
#define TMV_QR tmv::SV
#define TMV_QRD tmv::SVD
#endif

#ifdef TMV_NOLIB
#define DIVIDEUSING(A)
#else
#define DIVIDEUSING(A) A.divideUsing(TMV_QR)
#endif

#ifdef TISCOMPLEX
typedef std::complex<RT> T;
#define XFOUR 4
#else
typedef RT T;
#define XFOUR 1
#endif

#ifdef ONELOOP
const long long nloops1 = 1;
const long long nloops2 = 1;
#else
const long long MN = M*N;
const long long MK = M*K;
const long long NK = N*K;
const long long nloops2 = targetmem*1000 / ((2*MN+2*NK+2*MK)*sizeof(T)) + 1;
const long long nloops1 = targetnmegaflops*1000000 / (MN*(N+K)*nloops2*XFOUR) + 1;
#endif

#include <sys/time.h>
#include <algorithm>
#include <numeric>
#include <iostream>

#if (defined(DOEIGEN) || defined(DOEIGENSMALL))

#define ALLOC(X) Eigen::aligned_allocator< X >
#define EIGENP Eigen::VectorXi
#define EIGENQ Eigen::RowVectorXi

// Supposedly, these are required when using vector's of Eigen types,
// but they don't seem to work.  
// And leaving them out does seem to work.
//#define EIGEN_USE_NEW_STDVECTOR
//#include "Eigen/StdVector"
#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/QR"
#include "Eigen/SVD"

#if ALGO == 1
#define EIGENDIV .householderQr()
#define EIGENQR Eigen::HouseholderQR
#elif ALGO == 2
#define EIGENDIV .colPivHouseholderQr()
#define EIGENQR Eigen::ColPivHouseholderQR
#elif ALGO == 3
#define EIGENDIV .jacobiSvd()
#define EIGENQR Eigen::JacobiSVD
#endif

#endif

#if (PART1 == 1) // Full rectangle
#define MPART1(m) m
#define EPART1(m) m
#define EPART1t(m) m
#ifndef NO_SEPARATE_DECOMP
#define DO_SEPARATE_DECOMP
#endif

#endif


#ifdef DOLAP

#define LAP_BLOCKSIZE 64

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
#define BT T
#define BP(x) x
#endif

#ifdef TISCOMPLEX
#define BX(x) BP(&x)
#else
#define BX(x) x
#endif

#ifdef CBLAS
#define BLASNAME1(c,x) cblas_ ## c ## x
#define BLASINAME1(c,x) cblas_i ## c ## x
#define BLASM1
#define BLASCM CblasColMajor,
#define BLAST CblasTrans
#define BLASNT CblasNoTrans
#define BLASCT CblasConjTrans
#define BLASLeft CblasLeft
#define BLASRight CblasRight
#define BLASUpper CblasUpper
#define BLASLower CblasLower
#define BLASUnit CblasUnit
#define BLASNonUnit CblasNonUnit
#define BLAS1
#elif defined FBLAS
const char BlasCh_N = 'N';
const char BlasCh_T = 'T';
const char BlasCh_C = 'C';
const char BlasCh_L = 'L';
const char BlasCh_R = 'R';
const char BlasCh_U = 'U';
#define BLASNAME1(c,x) c ## x ## _
#define BLASINAME1(c,x) i ## c ## x ## _
#define BLASM1 -1
#define BLAS1 ,1
#define BLASCM 
#define BLASNT BlasCh_N
#define BLAST BlasCh_T
#define BLASCT BlasCh_C
#define BLASLeft BlasCh_L
#define BLASRight BlasCh_R
#define BLASUpper BlasCh_U
#define BLASLower BlasCh_L
#define BLASUnit BlasCh_U
#define BLASNonUnit BlasCh_N
#else
#error UNKNOWN BLAS definition....  
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
const T xmone(-1);
const T xmtwelve(-12);
const T xz89(8,9);
const T xz89c(8,-9);
const T xz71(7,1);
const T xz71c(7,-1);

const BT* zero = (const BT*) &xzero;
const BT* one = (const BT*) &xone;
const BT* seven = (const BT*) &xseven;
const BT* eight = (const BT*) &xeight;
const BT* mone = (const BT*) &xmone;
const BT* mtwelve = (const BT*) &xmtwelve;
const BT* z89 = (const BT*) &xz89;
const BT* z89c = (const BT*) &xz89c;
const BT* z71 = (const BT*) &xz71;
const BT* z71c = (const BT*) &xz71c;

const RT dmone(-1);
#else
const T zero(0);
const T one(1);
const T seven(7);
const T eight(8);
const T mone(-1);
const T mtwelve(-12);
const T dmone(-1);
#endif

#if defined(FLAP) && defined(TISCOMPLEX) && defined(TISFLOAT)
#define LP(x) (cfloat*) x
#elif defined(FBLAS) && defined(TISCOMPLEX) && !defined(TISFLOAT)
#define LP(x) (cdouble*) x
#else
#define LP(x) x
#endif

const int lap_worksize = N * 128;
#ifdef CLAP
#define LAPNAME1(c,x) clapack_ ## c ## x
#define LAPCM CblasColMajor,
RT* lap_rwork = new RT[lap_worksize];
T* lap_work = new T[lap_worksize];
#define LAPWORK lap_work
#define LAPRWORK lap_rwork
#elif defined FLAP
const char LapCh_N = 'N';
const char LapCh_T = 'T';
const char LapCh_C = 'C';
const char LapCh_L = 'L';
const char LapCh_R = 'R';
#define LAPNAME1(c,x) c ## x ## _
#define LAPCM 
RT* lap_rwork = new RT[lap_worksize];
T* lap_work = new T[lap_worksize];
#define LAPWORK lap_work
#define LAPRWORK lap_rwork
#define LAPNT LapCh_N
#define LAPT LapCh_T
#define LAPCT LapCh_C
#define LAPLeft LapCh_L
#define LAPRight LapCh_R
#else
#error UNKNOWN LAP definition....  
#define LAPNAME1(c,x) c ## x
#endif

int lapinfo;

#define LAPNAME2(c,x) LAPNAME1(c,x)
#define LAPNAME(x) LAPNAME2(BLAS_LETTER,x)

static void ConvertIndexToPermute(int n, const int* newIndex, int* P)
{
    // newIndex[i]=j means value at original j location needs to go to i.
    // Although newIndex is in FortranStyle, not CStyle, so use newIndex[i]-1.
    std::vector<int> currIndex(n);
    std::vector<int> origIndex(n);
    for(int i=0;i<n;++i) {
        currIndex[i] = i;
        origIndex[i] = i;
    }
    // currIndex[i]=j means value at original i location is currently at j.
    // origIndex[j]=i means value at original i location is currently at j.
    for(int i=0;i<n;++i) {
        int ip = currIndex[newIndex[i]-1];
        P[i] = ip;
        if (i != ip) {
            int origi = origIndex[i];
            int origip = origIndex[ip];
            currIndex[origi] = ip;
            currIndex[origip] = i;
            origIndex[i] = origip;
            origIndex[ip] = origi;
        }
    }
}
                    
#endif // DOLAP

#if defined(DOEIGEN) || defined(DOEIGENSMALL)
#ifdef TISCOMPLEX
#ifdef TISFLOAT
#define EIGENV Eigen::VectorXcf
#define EIGENM Eigen::MatrixXcf
#else
#define EIGENV Eigen::VectorXcd
#define EIGENM Eigen::MatrixXcd
#endif
#else
#ifdef TISFLOAT
#define EIGENV Eigen::VectorXf
#define EIGENM Eigen::MatrixXf
#else
#define EIGENV Eigen::VectorXd
#define EIGENM Eigen::MatrixXd
#endif
#endif

#endif

#ifdef DOEIGENSMALL
const int nloops2x = (
    (M*N*sizeof(T) < 256 * 1024 && 
     N*K*sizeof(T) < 256 * 1024 && 
     M!=Eigen::Dynamic && N!=Eigen::Dynamic && K!=Eigen::Dynamic) ? 
    nloops2 : 0 );

#define EIGENSMA Eigen::Matrix<T,M,N>
#define EIGENSMB Eigen::Matrix<T,N,M>
#define EIGENSMC Eigen::Matrix<T,M,K>
#define EIGENSMD Eigen::Matrix<T,N,K>
#define EIGENSME Eigen::Matrix<T,K,M>
#define EIGENSMF Eigen::Matrix<T,K,N>

#endif

static void ClearCache()
{
    static tmv::Vector<double> X(1000000,8.);
    static tmv::Vector<double> Y(1000000,8.);
    static tmv::Vector<double> Z(1000000);
    Z = X + Y;
    if (Norm(Z) < 5.) exit(1);
}

#ifdef DO_C
static void QR_C(
    const std::vector<tmv::Matrix<T> >& A1,
    std::vector<tmv::Matrix<T> >& B1,
    const std::vector<tmv::Matrix<T> >& C1,
    std::vector<tmv::Matrix<T> >& D1,
    const std::vector<tmv::Matrix<T> >& E1,
    std::vector<tmv::Matrix<T> >& F1)
{
    //std::cout<<"Start QR_C.\n";
    //std::cout<<"A1.size = "<<A1.size()<<std::endl;
#ifdef TMV_NO_LIB
    std::vector<TMV_QRD<tmv::Matrix<T> >*> QR1(nloops2,0);
#endif
    std::vector<T> d1(nloops2);

#ifdef ERRORCHECK
    std::vector<tmv::Matrix<T> > A0 = A1;
    std::vector<tmv::Matrix<T> > C0 = C1;
    std::vector<tmv::Matrix<T> > D0 = D1;
    std::vector<tmv::Matrix<T> > E0 = E1;
    std::vector<tmv::Matrix<T> > F0 = F1;
    std::vector<T> d0(nloops2);
    //std::cout<<"Made A0, etc."<<std::endl;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallMatrix<T,M,N> > A2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,M> > B2(nloops2);
    std::vector<tmv::SmallMatrix<T,M,K> > C2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,K> > D2(nloops2);
    std::vector<tmv::SmallMatrix<T,K,M> > E2(nloops2);
    std::vector<tmv::SmallMatrix<T,K,N> > F2(nloops2);
    std::vector<T> d2(nloops2);
    std::vector<TMV_QRD<tmv::SmallMatrix<T,M,N> >*> QR2(nloops2,0);
    //std::cout<<"Made A2, etc."<<std::endl;

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; C2[k] = C1[k]; E2[k] = E1[k];
    }
    //std::cout<<"Filled A2,C2,E2"<<std::endl;
#endif

#ifdef DOLAP
    std::vector<tmv::Matrix<T> > A3 = A1;
    std::vector<tmv::Matrix<T> > B3 = B1;
    std::vector<tmv::Matrix<T> > C3 = C1;
    std::vector<tmv::Matrix<T> > D3 = D1;
    std::vector<tmv::Matrix<T> > E3 = E1;
    std::vector<tmv::Matrix<T> > F3 = F1;
    std::vector<T> d3(nloops2);
    std::vector<tmv::Matrix<T> > QR3 = A1;
    std::vector<tmv::Vector<T> > Beta3(nloops2);
    std::vector<std::vector<int> > P3(nloops2);
    std::vector<int> P3i(N);
    tmv::Matrix<T> G3(M,K);
    tmv::Matrix<T> Q3(M,N);
    int lwork = M*N;
    tmv::Vector<T> Work3(lwork);
#ifdef TISCOMPLEX
    tmv::Vector<RT> RWork3(lwork);
#endif
    for(int k=0;k<nloops2;++k) {
        P3[k].resize(N);
        Beta3[k].resize(N);
    }
    //std::cout<<"Made A3, etc."<<std::endl;
#endif

#ifdef DOEIGEN
    std::vector<EIGENM,ALLOC(EIGENM) > A4(nloops2,EIGENM(M,N));
    std::vector<EIGENM,ALLOC(EIGENM) > B4(nloops2,EIGENM(N,M));
    std::vector<EIGENM,ALLOC(EIGENM) > C4(nloops2,EIGENM(M,K));
    std::vector<EIGENM,ALLOC(EIGENM) > D4(nloops2,EIGENM(N,K));
    std::vector<EIGENM,ALLOC(EIGENM) > E4(nloops2,EIGENM(K,M));
    std::vector<EIGENM,ALLOC(EIGENM) > F4(nloops2,EIGENM(K,N));
    EIGENM Q4(M,N);
    std::vector<T> d4(nloops2);
    std::vector<EIGENQR<EIGENM>*,ALLOC(EIGENQR<EIGENM>*) > QR4(nloops2);
    for(int k=0;k<nloops2;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) A4[k](i,j) = A1[k](i,j);
        for(int j=0;j<K;++j) for(int i=0;i<M;++i) C4[k](i,j) = C1[k](i,j);
        for(int j=0;j<M;++j) for(int i=0;i<K;++i) E4[k](i,j) = E1[k](i,j);
    }
    //std::cout<<"Made A4, etc."<<std::endl;
#endif

#ifdef DOEIGENSMALL
    std::vector<EIGENSMA,ALLOC(EIGENSMA) > A5;
    std::vector<EIGENSMB,ALLOC(EIGENSMB) > B5;
    std::vector<EIGENSMC,ALLOC(EIGENSMC) > C5;
    std::vector<EIGENSMD,ALLOC(EIGENSMD) > D5;
    std::vector<EIGENSME,ALLOC(EIGENSME) > E5;
    std::vector<EIGENSMF,ALLOC(EIGENSMF) > F5;
    std::vector<EIGENSMA,ALLOC(EIGENSMA) > Q5;
    std::vector<T> d5(nloops2x);
    std::vector<EIGENQR<EIGENSMA>*,ALLOC(EIGENQR<EIGENSMA>*) > QR5;
    //std::cout<<"Made vectors of A5, etc."<<std::endl;
    if (nloops2x) { 
        A5.resize(nloops2x); 
        B5.resize(nloops2x); 
        C5.resize(nloops2x);
        D5.resize(nloops2x);
        E5.resize(nloops2x); 
        F5.resize(nloops2x); 
        if (nloops2x > 0) Q5.resize(1);
        QR5.resize(nloops2x,0); 
    }
    //std::cout<<"Resized vectors of A5, etc."<<std::endl;
    for(int k=0;k<nloops2x;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) A5[k](i,j) = A1[k](i,j);
        for(int j=0;j<K;++j) for(int i=0;i<M;++i) C5[k](i,j) = C1[k](i,j);
        for(int j=0;j<M;++j) for(int i=0;i<K;++i) E5[k](i,j) = E1[k](i,j);
    }
    //std::cout<<"Made A5, etc."<<std::endl;
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
    double t12_reg=0., t12_small=0., t12_blas=0., t12_eigen=0., t12_smalleigen=0.;
    double t13_reg=0., t13_small=0., t13_blas=0., t13_eigen=0., t13_smalleigen=0.;
    double t14_reg=0., t14_small=0., t14_blas=0., t14_eigen=0., t14_smalleigen=0.;
    double t15_reg=0., t15_small=0., t15_blas=0., t15_eigen=0., t15_smalleigen=0.;
    double t16_reg=0., t16_small=0., t16_blas=0., t16_eigen=0., t16_smalleigen=0.;
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
    double e12_reg=0., e12_small=0., e12_blas=0., e12_eigen=0., e12_smalleigen=0.;
    double e13_reg=0., e13_small=0., e13_blas=0., e13_eigen=0., e13_smalleigen=0.;
    double e14_reg=0., e14_small=0., e14_blas=0., e14_eigen=0., e14_smalleigen=0.;
    double e15_reg=0., e15_small=0., e15_blas=0., e15_eigen=0., e15_smalleigen=0.;
    double e16_reg=0., e16_small=0., e16_blas=0., e16_eigen=0., e16_smalleigen=0.;
#endif

    for (int n=0; n<nloops1; ++n) {
        //std::cout<<"start loop n = "<<n<<std::endl;

#ifdef DO_SEPARATE_DECOMP
#if 1 // A -> QR
        ClearCache();
        //std::cout<<"After clear cache"<<std::endl;

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            //std::cout<<"A1[k] = "<<A1[k]<<std::endl;
            QR1[k] = new TMV_QRD<tmv::Matrix<T> >(A1[k]);
#else
            //std::cout<<"A1[k] = "<<A1[k]<<std::endl;
            A1[k].saveDiv();
            //std::cout<<"After saveDiv"<<std::endl;
            DIVIDEUSING(A1[k]);
            //std::cout<<"After divideUsing"<<std::endl;
            A1[k].setDiv();
            //std::cout<<"After setDiv"<<std::endl;
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
#if ALGO == 1
                tmv::Matrix<T> temp = QR1[k]->getQ() * QR1[k]->getR();
#elif ALGO == 2
                tmv::Matrix<T> temp = QR1[k]->getQ() * QR1[k]->getR() * QR1[k]->getP();
#else
#endif
                if (QR1[k]->isTrans()) A0[k] = temp.transpose();
                else A0[k] = temp;
#else // !NO_LIB
#if ALGO == 1
                tmv::Matrix<T> temp = A1[k].qrd().getQ() * A1[k].qrd().getR();
#elif ALGO == 2
                tmv::Matrix<T> temp = A1[k].qrpd().getQ() * A1[k].qrpd().getR() * A1[k].qrpd().getP();
#else 
#endif
                if (GETDIV(A1[k]).isTrans()) A0[k] = temp.transpose();
                else A0[k] = temp;
#endif
                //std::cout<<"A0[k] => "<<A0[k]<<std::endl;
                e1_reg += NormSq(A0[k]-A1[k]);
                //std::cout<<"diff = "<<A0[k]-A1[k]<<std::endl;
                //std::cout<<"NormSq(diff) = "<<e1_reg<<std::endl;
            }
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A1[0]);
            //std::cout<<"e1_reg = "<<e1_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            //std::cout<<"k = "<<k<<std::endl;
            QR2[k] = new TMV_QRD<tmv::SmallMatrix<T,M,N> >(A2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) {
#if ALGO == 1
                tmv::Matrix<T> temp = QR2[k]->getQ() * QR2[k]->getR();
#elif ALGO == 2
                tmv::Matrix<T> temp = QR2[k]->getQ() * QR2[k]->getR() * QR2[k]->getP();
#else
#endif
                if (QR2[k]->isTrans()) A0[k] = temp.transpose();
                else A0[k] = temp;
                e1_small += NormSq(A2[k]-A0[k]);
                //std::cout<<"e1_small => "<<e1_small<<std::endl;
            }
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A1[0]);
            //std::cout<<"e1_small = "<<e1_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            QR3[k] = A3[k];
#if ALGO == 1
            LAPNAME(geqrf) (
                LAPCM M,N,LP(QR3[k].ptr()),M,LP(Beta3[k].ptr()),
                LP(Work3.ptr()),lwork,&lapinfo);
#elif ALGO == 2
            for(int i=0;i<N;++i) P3i[i] = 0;
            LAPNAME(geqp3) (
                LAPCM M,N,LP(QR3[k].ptr()),M,&(P3i[0]),LP(Beta3[k].ptr()),
                LP(Work3.ptr()),lwork,
#ifdef TISCOMPLEX
                RWork3.ptr(),
#endif
                &lapinfo);
            ConvertIndexToPermute(N,&(P3i[0]),&(P3[k][0]));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = QR3[k];
#ifdef TISCOMPLEX
                LAPNAME(ungqr) (
                    LAPCM M,N,N,LP(A0[k].ptr()),M,LP(Beta3[k].cptr()),
                    LP(Work3.ptr()),lwork,&lapinfo);
#else
                LAPNAME(orgqr) (
                    LAPCM M,N,N,LP(A0[k].ptr()),M,LP(Beta3[k].cptr()),
                    LP(Work3.ptr()),lwork,&lapinfo);
#endif
                A0[k] *= QR3[k].upperTri();
#if ALGO == 2
                A0[k].reversePermuteCols(&P3[k][0]);
#endif
                e1_blas += NormSq(A3[k]-A0[k]);
            }
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A1[0]);
            //std::cout<<"e1_blas = "<<e1_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            QR4[k] = new EIGENQR<EIGENM>(A4[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
#if (ALGO == 1) || (ALGO == 2)
                EIGENM temp1 = 
                    EIGENM(QR4[k]->householderQ()) *
                    EIGENM(QR4[k]->matrixQR().triangularView<Eigen::Upper>());
#endif
#if ALGO == 2
                temp1.applyOnTheRight(QR4[k]->colsPermutation().inverse());
#endif
                e1_eigen += (A4[k]-temp1).squaredNorm();
            }
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A1[0]);
            //std::cout<<"e1_eigen = "<<e1_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        //std::cout<<"start EIGENSMALL QRDecompose"<<std::endl;
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            QR5[k] = new EIGENQR<EIGENSMA>(A5[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
#if (ALGO == 1) || (ALGO == 2)
                EIGENM temp1 = 
                    EIGENM(QR5[k]->householderQ()) *
                    EIGENM(QR5[k]->matrixQR().triangularView<Eigen::Upper>());
#endif
#if ALGO == 2
                temp1.applyOnTheRight(QR5[k]->colsPermutation().inverse());
#endif
                e1_smalleigen += (A5[k]-temp1).squaredNorm();
            }
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A1[0]);
            //std::cout<<"e1_smalleigen = "<<e1_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // D /= A (no decomp)
#ifdef AISSQUARE
#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                A0[k] = A1[k];
                D0[k] = (A0[k].adjoint()*C0[k]) / (A0[k].adjoint()*A0[k]);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) D1[k] = C1[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            QR1[k]->solveInPlace(D1[k]);
#else
            D1[k] /= A1[k];
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                e2_reg += NormSq(D0[k]-D1[k]);
            }
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= Norm(D0[0]);
            //std::cout<<"e2_reg = "<<e2_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            QR2[k]->solveInPlace(D2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                e2_small += NormSq(D0[k]-D2[k]);
            }
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= Norm(D0[0]);
            //std::cout<<"e2_small = "<<e2_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        for (int k=0; k<nloops2; ++k) D3[k] = C3[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
            LAPNAME(unmqr) (
                LAPCM LAPLeft, LAPCT, M,K,N,LP(QR3[k].cptr()),M,
                LP(Beta3[k].cptr()),LP(D3[k].ptr()),N,
                LP(Work3.ptr()),lwork,&lapinfo);
#else
            LAPNAME(ormqr) (
                LAPCM LAPLeft, LAPT, M,K,N,LP(QR3[k].cptr()),M,
                LP(Beta3[k].cptr()),LP(D3[k].ptr()),N,
                LP(Work3.ptr()),lwork,&lapinfo);
#endif
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUpper, BLASNT, BLASNonUnit,
                N,K,one,BP(QR3[k].cptr()),N,BP(D3[k].ptr()),N
                BLAS1 BLAS1 BLAS1 BLAS1);
#if ALGO == 2
            D3[k].reversePermuteRows(&P3[k][0]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                e2_blas += NormSq(D0[k]-D3[k]);
            }
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= Norm(D0[0]);
            //std::cout<<"e2_blas = "<<e2_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D4[k] = QR4[k]->solve(D4[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<N;++i) for(int j=0;j<K;++j)
                    e2_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            }
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm(D0[0]);
            //std::cout<<"e2_eigen = "<<e2_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        //std::cout<<"start EIGENSMALL /= "<<std::endl;
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            D5[k] = QR5[k]->solve(D5[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                for(int i=0;i<N;++i) for(int j=0;j<K;++j)
                    e2_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            }
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm(D0[0]);
            //std::cout<<"e2_smalleigen = "<<e2_smalleigen<<std::endl;
        }
#endif
#endif
#endif
#endif

#if 1 // D = C / A (no decomp)
#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                // Solve AD = C
                // AtAD = AtC
                // D = (AtC) / (AtA)
                A0[k] = A1[k];
                D0[k] = (A0[k].adjoint()*C0[k]) / (A0[k].adjoint()*A0[k]);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            QR1[k]->solve(C1[k],D1[k]);
#else
            D1[k] = C1[k] / A1[k];
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"C1 = "<<C1[k]<<std::endl;
                //std::cout<<"D1 = "<<D1[k]<<std::endl;
                //std::cout<<"D0 = "<<D0[k]<<std::endl;
                e3_reg += NormSq(D0[k]-D1[k]);
            }
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= Norm(D0[0]);
            //std::cout<<"e3_reg = "<<e3_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            QR2[k]->solve(C2[k],D2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                e3_small += NormSq(D0[k]-D2[k]);
            }
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= Norm(D0[0]);
            //std::cout<<"e3_small = "<<e3_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // QRD=C
            // RD = QtC
            // D = R^-1QtC
            G3 = C3[k];
#ifdef TISCOMPLEX
            LAPNAME(unmqr) (
                LAPCM LAPLeft, LAPCT, M,K,N,LP(QR3[k].cptr()),M,
                LP(Beta3[k].cptr()),LP(G3.ptr()),M,
                LP(Work3.ptr()),lwork,&lapinfo);
#else
            LAPNAME(ormqr) (
                LAPCM LAPLeft, LAPT, M,K,N,LP(QR3[k].cptr()),M,
                LP(Beta3[k].cptr()),LP(G3.ptr()),M,
                LP(Work3.ptr()),lwork,&lapinfo);
#endif
            D3[k] = G3.rowRange(0,N);
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUpper, BLASNT, BLASNonUnit,
                N,K,one,BP(QR3[k].cptr()),M,BP(D3[k].ptr()),N
                BLAS1 BLAS1 BLAS1 BLAS1);
#if ALGO == 2
            D3[k].reversePermuteRows(&P3[k][0]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                e3_blas += NormSq(D0[k]-D3[k]);
            }
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= Norm(D0[0]);
            //std::cout<<"e3_blas = "<<e3_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D4[k] = QR4[k]->solve(C4[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<N;++i) for(int j=0;j<K;++j)
                    e3_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            }
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(D0[0]);
            //std::cout<<"e3_eigen = "<<e3_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        //std::cout<<"start EIGENSMALL / "<<std::endl;
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            D5[k] = QR5[k]->solve(C5[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                for(int i=0;i<N;++i) for(int j=0;j<K;++j)
                    e3_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            }
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(D0[0]);
            //std::cout<<"e3_smalleigen = "<<e3_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // B = A.inverse() (no decomp)
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            QR1[k]->makeInverse(B1[k]);
#else
            B1[k] = A1[k].inverse();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"A1 = "<<A1[k]<<std::endl;
                //std::cout<<"B1 = "<<B1[k]<<std::endl;
                tmv::Matrix<T> temp = B1[k] * MPART1(A1[k]);
                //std::cout<<"temp = "<<temp<<std::endl;
                e4_reg += NormSq(temp-1);
            }
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= Norm(A1[0]);
            //std::cout<<"e4_reg = "<<e4_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            QR2[k]->makeInverse(B2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                tmv::Matrix<T> temp = B2[k] * MPART1(A2[k]);
                e4_small += NormSq(temp-1);
            }
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(A1[0]);
            //std::cout<<"e4_small = "<<e4_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            Q3 = QR3[k];
#ifdef TISCOMPLEX
            LAPNAME(ungqr) (
                LAPCM M,N,N,LP(Q3.ptr()),M,LP(Beta3[k].cptr()),
                LP(Work3.ptr()),lwork,&lapinfo);
            Q3.conjugateSelf();
#else
            LAPNAME(orgqr) (
                LAPCM M,N,N,LP(Q3.ptr()),M,LP(Beta3[k].cptr()),
                LP(Work3.ptr()),lwork,&lapinfo);
#endif
            BLASNAME(trsm) (
                BLASCM BLASRight, BLASUpper, BLAST, BLASNonUnit,
                M,N,one,BP(QR3[k].cptr()),M,BP(Q3.ptr()),M
                BLAS1 BLAS1 BLAS1 BLAS1);
            B3[k] = Q3.transpose();
#if ALGO == 2
            B3[k].reversePermuteRows(&P3[k][0]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"A3[k] = "<<A3[k]<<std::endl;
                //std::cout<<"B3[k] = "<<B3[k]<<std::endl;
                tmv::Matrix<T> temp = B3[k] * MPART1(A3[k]);
                //std::cout<<"B3[k]*A3[k] = "<<temp<<std::endl;
                e4_blas += NormSq(temp-1);
            }
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(A1[0]);
            //std::cout<<"e4_blas = "<<e4_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (ALGO == 1) || (ALGO == 2)
            // A = QR
            // A^-1 = R^-1 Q^-1
            // B = R^-1 Qt
            // Bt = Q Rt^-1
            Q4.block(0,0,N,N).setIdentity();
            QR4[k]->matrixQR().block(0,0,N,N).triangularView<Eigen::Upper>().adjoint().solveInPlace(Q4.block(0,0,N,N));
            Q4.block(N,0,M-N,N).setZero();
            Q4.applyOnTheLeft(QR4[k]->householderQ());
#if ALGO == 2
            Q4.applyOnTheRight(QR4[k]->colsPermutation().inverse());
#endif
            B4[k] = Q4.adjoint();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"A4 = "<<A4[k]<<std::endl;
                //std::cout<<"B4 = "<<B4[k]<<std::endl;
                EIGENM temp1 = B4[k] * EPART1(A4[k]);
                e4_eigen += (temp1-temp1.Identity(N,N)).squaredNorm();
            }
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(A1[0]);
            //std::cout<<"e4_eigen = "<<e4_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        //std::cout<<"start EIGENSMALL inverse "<<std::endl;
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if (ALGO == 1) || (ALGO == 2)
            Q5[0].block(0,0,N,N).setIdentity();
            QR5[k]->matrixQR().block(0,0,N,N).triangularView<Eigen::Upper>().adjoint().solveInPlace(Q5[0].block(0,0,N,N));
            Q5[0].block(N,0,M-N,N).setZero();
            Q5[0].applyOnTheLeft(QR5[k]->householderQ());
#if ALGO == 2
            Q5[0].applyOnTheRight(QR5[k]->colsPermutation().inverse());
#endif
            B5[k] = Q5[0].adjoint();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                EIGENM temp1 = B5[k] * EPART1(A5[k]);
                e4_smalleigen += (temp1-temp1.Identity(N,N)).squaredNorm();
            }
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(A1[0]);
            //std::cout<<"e4_smalleigen = "<<e4_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // d = A.det() (no decomp)
#ifdef AISSQUARE
#ifdef ERRORCHECK
        for (int k=0; k<nloops2; ++k) {
            A0[k] = A1[k];
            d0[k] = MPART1(A0[k]).det();
        }
        //if (n==0) std::cout<<"det = "<<d0[0]<<std::endl;
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            d1[k] = QR1[k]->det();
#else
            d1[k] = A1[k].det();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"d1 = "<<d1[k]<<std::endl;
                //std::cout<<"d0 = "<<d0[k]<<std::endl;
                e11_reg += tmv::TMV_NORM(d0[k]-d1[k]);
            }
            e11_reg = sqrt(e11_reg/nloops2);
            e11_reg /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e11_reg = "<<e11_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            d2[k] = QR2[k]->det();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                e11_small += tmv::TMV_NORM(d0[k]-d2[k]);
            }
            e11_small = sqrt(e11_small/nloops2);
            e11_small /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e11_small = "<<e11_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            d3[k] = T(1);
            for(int i=0;i<N;++i) d3[k] *= QR3[k](i,i);
            for(int i=0;i<N;++i) if (Beta3[k][i] != T(0)) {
#ifdef TISCOMPLEX
                T b = Beta3[k][i];
                d3[k] *= -(b*b)/std::norm(b);
#else
                d3[k] = -d3[k];
#endif
            }
#if ALGO == 2
            for(int i=0;i<N;++i) if (P3[k][i] != i) d3[k] = -d3[k];
#endif
        }
        //if (n==0) std::cout<<"LAP det = "<<d3[0]<<std::endl;

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                e11_blas += tmv::TMV_NORM(d0[k]-d3[k]);
            }
            e11_blas = sqrt(e11_blas/nloops2);
            e11_blas /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e11_blas = "<<e11_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            d4[k] = QR4[k]->absDeterminant();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e11_eigen += tmv::TMV_NORM(std::abs(d0[k])-d4[k]);
            }
            e11_eigen = sqrt(e11_eigen/nloops2);
            e11_eigen /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e11_eigen = "<<e11_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        //std::cout<<"start EIGENSMALL det "<<std::endl;
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            d5[k] = QR5[k]->absDeterminant();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e11_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                e11_smalleigen += tmv::TMV_NORM(std::abs(d0[k])-d5[k]);
            }
            e11_smalleigen = sqrt(e11_smalleigen/nloops2);
            e11_smalleigen /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e11_smalleigen = "<<e11_smalleigen<<std::endl;
        }
#endif
#endif
#endif
#endif

#ifdef DOREG
#ifndef TMV_NO_LIB
        for (int k=0; k<nloops2; ++k) 
        { A1[k].unsetDiv(); A1[k].unsaveDiv(); }
#endif
#endif
#endif // separate decomposition and use

#if 1 // D /= A
#ifdef AISSQUARE
#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                A0[k] = A1[k];
                D0[k] = (A0[k].adjoint()*C0[k]) / (A0[k].adjoint()*A0[k]);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) D1[k] = C1[k]; 
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            GETDIV(A1[k]).solveInPlace(D1[k]);
#else
            DIVIDEUSING(A1[k]);
            D1[k] /= MPART1(A1[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"C1 = "<<C1[k]<<std::endl;
                //std::cout<<"D1 = "<<D1[k]<<std::endl;
                //std::cout<<"D0 = "<<D0[k]<<std::endl;
                e5_reg += NormSq(D0[k]-D1[k]);
            }
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= Norm(D0[0]);
            //std::cout<<"e5_reg = "<<e5_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ALGO == 1 && !defined(AISSQUARE)
            D2[k] /= MPART1(A2[k]);
#else
            GETDIV(A2[k]).solveInPlace(D2[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                e5_small += NormSq(D0[k]-D2[k]);
            }
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= Norm(D0[0]);
            //std::cout<<"e5_small = "<<e5_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        for (int k=0; k<nloops2; ++k) D3[k] = C3[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            QR3[k] = A3[k];
#if ALGO == 1
            LAPNAME(gels) (
                LAPCM LAPNT,M,N,K,LP(QR3[k].ptr()),M,LP(D3[k].ptr()),N,
                LP(Work3.ptr()),lwork,&lapinfo);
#elif ALGO == 2
            RT rcond = 1.e-6;
            int rank;
            for(int i=0;i<N;++i) P3i[i] = 0;
            LAPNAME(gelsy) (
                LAPCM M,N,K,LP(QR3[k].ptr()),M,LP(D3[k].ptr()),N,&(P3i[0]),
                rcond,&rank,LP(Work3.ptr()),lwork,
#ifdef TISCOMPLEX
                RWork3.ptr(),
#endif
                &lapinfo);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                e5_blas += NormSq(D0[k]-D3[k]);
            }
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= Norm(D0[0]);
            //std::cout<<"e5_blas = "<<e5_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D4[k] = A4[k] EIGENDIV .solve(D4[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<N;++i) for(int j=0;j<K;++j)
                    e5_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            }
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(D0[0]);
            //std::cout<<"e5_eigen = "<<e5_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        //std::cout<<"start EIGENSMALL /= "<<std::endl;
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            D5[k] = A5[k] EIGENDIV .solve(D5[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                for(int i=0;i<N;++i) for(int j=0;j<K;++j)
                    e5_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            }
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(D0[0]);
            //std::cout<<"e5_smalleigen = "<<e5_smalleigen<<std::endl;
        }
#endif
#endif
#endif
#endif

#if 1 // D = C / A
#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                // Solve AD = C
                // AtAD = AtC
                // D = (AtC) / (AtA)
                A0[k] = A1[k];
                D0[k] = (A0[k].adjoint()*C0[k]) / (A0[k].adjoint()*A0[k]);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            GETDIV(A1[k]).solve(C1[k],D1[k]);
#else
            DIVIDEUSING(A1[k]);
            D1[k] = C1[k] / MPART1(A1[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                double x = NormSq(D0[k]-D1[k]);
                e6_reg += NormSq(D0[k]-D1[k]);
            }
            //std::cout<<"Done loop: e6_reg = "<<e6_reg<<std::endl;
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= Norm(D0[0]);
            //std::cout<<"e6_reg = "<<e6_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ALGO == 1 && !defined(AISSQUARE)
            D2[k] = C2[k] / MPART1(A2[k]);
#else
            GETDIV(A2[k]).solve(C2[k],D2[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                e6_small += NormSq(D0[k]-D2[k]);
            }
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= Norm(D0[0]);
            //std::cout<<"e6_small = "<<e6_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            QR3[k] = A3[k];
            G3 = C3[k];
#if ALGO == 1
            LAPNAME(gels) (
                LAPCM LAPNT,M,N,K,LP(QR3[k].ptr()),M,LP(G3.ptr()),M,
                LP(Work3.ptr()),lwork,&lapinfo);
#elif ALGO == 2
            RT rcond = 1.e-6;
            int rank;
            for(int i=0;i<N;++i) P3i[i] = 0;
            LAPNAME(gelsy) (
                LAPCM M,N,K,LP(QR3[k].ptr()),M,LP(G3.ptr()),M,&(P3i[0]),
                rcond,&rank,LP(Work3.ptr()),lwork,
#ifdef TISCOMPLEX
                RWork3.ptr(),
#endif
                &lapinfo);
#endif
            D3[k] = G3.rowRange(0,N);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                e6_blas += NormSq(D0[k]-D3[k]);
            }
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= Norm(D0[0]);
            //std::cout<<"e6_blas = "<<e6_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D4[k] = A4[k] EIGENDIV .solve(C4[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<N;++i) for(int j=0;j<K;++j)
                    e6_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            }
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(D0[0]);
            //std::cout<<"e6_eigen = "<<e6_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        //std::cout<<"start EIGENSMALL / "<<std::endl;
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            D5[k] = A5[k] EIGENDIV .solve(C5[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                for(int i=0;i<N;++i) for(int j=0;j<K;++j)
                    e6_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            }
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(D0[0]);
            //std::cout<<"e6_smalleigen = "<<e6_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // B = A.inverse()
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            GETDIV(A1[k]).makeInverse(B1[k]);
#else
            DIVIDEUSING(A1[k]);
            B1[k] = MPART1(A1[k]).inverse();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"A1 = "<<A1[k]<<std::endl;
                //std::cout<<"B1 = "<<B1[k]<<std::endl;
                tmv::Matrix<T> temp = B1[k] * MPART1(A1[k]);
                //std::cout<<"temp = "<<temp<<std::endl;
                e7_reg += NormSq(temp-1);
            }
            e7_reg = sqrt(e7_reg/nloops2);
            e7_reg /= Norm(A1[0]);
            //std::cout<<"e7_reg = "<<e7_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ALGO == 1 && !defined(AISSQUARE)
            B2[k] = MPART1(A2[k]).inverse();
#else
            GETDIV(A2[k]).makeInverse(B2[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                tmv::Matrix<T> temp = B2[k] * MPART1(A2[k]);
                e7_small += NormSq(temp-1);
            }
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= Norm(A1[0]);
            //std::cout<<"e7_small = "<<e7_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            QR3[k] = A3[k];
#if ALGO == 1
            LAPNAME(geqrf) (
                LAPCM M,N,LP(QR3[k].ptr()),M,LP(Beta3[k].ptr()),
                LP(Work3.ptr()),lwork,&lapinfo);
#elif ALGO == 2
            for(int i=0;i<N;++i) P3i[i] = 0;
            LAPNAME(geqp3) (
                LAPCM M,N,LP(QR3[k].ptr()),M,&(P3i[0]),LP(Beta3[k].ptr()),
                LP(Work3.ptr()),lwork,
#ifdef TISCOMPLEX
                RWork3.ptr(),
#endif
                &lapinfo);
            ConvertIndexToPermute(N,&(P3i[0]),&(P3[k][0]));
#endif
            Q3 = QR3[k];
#ifdef TISCOMPLEX
            LAPNAME(ungqr) (
                LAPCM M,N,N,LP(Q3.ptr()),M,LP(Beta3[k].cptr()),
                LP(Work3.ptr()),lwork,&lapinfo);
            Q3.conjugateSelf();
#else
            LAPNAME(orgqr) (
                LAPCM M,N,N,LP(Q3.ptr()),M,LP(Beta3[k].cptr()),
                LP(Work3.ptr()),lwork,&lapinfo);
#endif
            BLASNAME(trsm) (
                BLASCM BLASRight, BLASUpper, BLAST, BLASNonUnit,
                M,N,one,BP(QR3[k].cptr()),M,BP(Q3.ptr()),M
                BLAS1 BLAS1 BLAS1 BLAS1);
            B3[k] = Q3.transpose();
#if ALGO == 2
            B3[k].reversePermuteRows(&P3[k][0]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                tmv::Matrix<T> temp = B3[k] * MPART1(A3[k]);
                e7_blas += NormSq(temp-1);
            }
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= Norm(A1[0]);
            //std::cout<<"e7_blas = "<<e7_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (ALGO == 1) || (ALGO == 2)
            EIGENQR<EIGENM> qr(A4[k]);
            Q4.block(0,0,N,N).setIdentity();
            qr.matrixQR().block(0,0,N,N).triangularView<Eigen::Upper>().adjoint().solveInPlace(Q4.block(0,0,N,N));
            Q4.block(N,0,M-N,N).setZero();
            Q4.applyOnTheLeft(qr.householderQ());
#if ALGO == 2
            Q4.applyOnTheRight(qr.colsPermutation().inverse());
#endif
            B4[k] = Q4.adjoint();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = B4[k] * EPART1(A4[k]);
                e7_eigen += (temp1-temp1.Identity(N,N)).squaredNorm();
            }
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= Norm(A1[0]);
            //std::cout<<"e7_eigen = "<<e7_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        //std::cout<<"start EIGENSMALL inverse"<<std::endl;
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if (ALGO == 1) || (ALGO == 2)
            EIGENQR<EIGENSMA> qr(A5[k]);
            Q5[0].block(0,0,N,N).setIdentity();
            qr.matrixQR().block(0,0,N,N).triangularView<Eigen::Upper>().adjoint().solveInPlace(Q5[0].block(0,0,N,N));
            Q5[0].block(N,0,M-N,N).setZero();
            Q5[0].applyOnTheLeft(qr.householderQ());
#if ALGO == 2
            Q5[0].applyOnTheRight(qr.colsPermutation().inverse());
#endif
            B5[k] = Q5[0].adjoint();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                EIGENM temp1 = B5[k] * EPART1(A5[k]);
                e7_smalleigen += (temp1-temp1.Identity(N,N)).squaredNorm();
            }
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= Norm(A1[0]);
            //std::cout<<"e7_smalleigen = "<<e7_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // B = B.inverse()
#ifdef AISSQUARE
        ClearCache();

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) B1[k] = A1[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            GETDIV(B1[k]).makeInverse(B1[k]);
#else
            DIVIDEUSING(B1[k]);
            B1[k] = MPART1(B1[k]).inverse();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"A1 = "<<A1[k]<<std::endl;
                //std::cout<<"B1 = "<<B1[k]<<std::endl;
                tmv::Matrix<T> temp = B1[k] * MPART1(A1[k]);
                //std::cout<<"temp = "<<temp<<std::endl;
                e8_reg += NormSq(temp-1);
            }
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= Norm(A1[0]);
            //std::cout<<"e8_reg = "<<e8_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) B2[k] = A2[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ALGO == 1 && !defined(AISSQUARE)
            B2[k] = MPART1(B2[k]).inverse();
#else
            GETDIV(B2[k]).makeInverse(B2[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                tmv::Matrix<T> temp = B2[k] * MPART1(A2[k]);
                e8_small += NormSq(temp-1);
            }
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= Norm(A1[0]);
            //std::cout<<"e8_small = "<<e8_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        for (int k=0; k<nloops2; ++k) B3[k] = A3[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            QR3[k] = B3[k];
#if ALGO == 1
            LAPNAME(geqrf) (
                LAPCM M,N,LP(QR3[k].ptr()),M,LP(Beta3[k].ptr()),
                LP(Work3.ptr()),lwork,&lapinfo);
#elif ALGO == 2
            for(int i=0;i<N;++i) P3i[i] = 0;
            LAPNAME(geqp3) (
                LAPCM M,N,LP(QR3[k].ptr()),M,&(P3i[0]),LP(Beta3[k].ptr()),
                LP(Work3.ptr()),lwork,
#ifdef TISCOMPLEX
                RWork3.ptr(),
#endif
                &lapinfo);
            ConvertIndexToPermute(N,&(P3i[0]),&(P3[k][0]));
#endif
            Q3 = QR3[k];
#ifdef TISCOMPLEX
            LAPNAME(ungqr) (
                LAPCM M,N,N,LP(Q3.ptr()),M,LP(Beta3[k].cptr()),
                LP(Work3.ptr()),lwork,&lapinfo);
            Q3.conjugateSelf();
#else
            LAPNAME(orgqr) (
                LAPCM M,N,N,LP(Q3.ptr()),M,LP(Beta3[k].cptr()),
                LP(Work3.ptr()),lwork,&lapinfo);
#endif
            BLASNAME(trsm) (
                BLASCM BLASRight, BLASUpper, BLAST, BLASNonUnit,
                M,N,one,BP(QR3[k].cptr()),M,BP(Q3.ptr()),M
                BLAS1 BLAS1 BLAS1 BLAS1);
            B3[k] = Q3.transpose();
#if ALGO == 2
            B3[k].reversePermuteRows(&P3[k][0]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                tmv::Matrix<T> temp = B3[k] * MPART1(A3[k]);
                e8_blas += NormSq(temp-1);
            }
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= Norm(A1[0]);
            //std::cout<<"e8_blas = "<<e8_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) B4[k] = A4[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (ALGO == 1) || (ALGO == 2)
            EIGENQR<EIGENM> qr(B4[k]);
            Q4.block(0,0,N,N).setIdentity();
            qr.matrixQR().block(0,0,N,N).triangularView<Eigen::Upper>().adjoint().solveInPlace(Q4.block(0,0,N,N));
            Q4.block(N,0,M-N,N).setZero();
            Q4.applyOnTheLeft(qr.householderQ());
#if ALGO == 2
            Q4.applyOnTheRight(qr.colsPermutation().inverse());
#endif
            B4[k] = Q4.adjoint();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = B4[k] * EPART1(A4[k]);
                e8_eigen += (temp1-temp1.Identity(N,N)).squaredNorm();
            }
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= Norm(A1[0]);
            //std::cout<<"e8_eigen = "<<e8_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        //std::cout<<"start EIGENSMALL inverse 2"<<std::endl;
        for (int k=0; k<nloops2x; ++k) B5[k] = A5[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if (ALGO == 1) || (ALGO == 2)
            EIGENQR<EIGENSMA> qr(B5[k]);
            Q5[0].block(0,0,N,N).setIdentity();
            qr.matrixQR().block(0,0,N,N).triangularView<Eigen::Upper>().adjoint().solveInPlace(Q5[0].block(0,0,N,N));
            Q5[0].block(N,0,M-N,N).setZero();
            Q5[0].applyOnTheLeft(qr.householderQ());
#if ALGO == 2
            Q5[0].applyOnTheRight(qr.colsPermutation().inverse());
#endif
            B5[k] = Q5[0].adjoint();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                EIGENM temp1 = B5[k] * EPART1(A5[k]);
                e8_smalleigen += (temp1-temp1.Identity(N,N)).squaredNorm();
            }
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= Norm(A1[0]);
            //std::cout<<"e8_smalleigen = "<<e8_smalleigen<<std::endl;
        }
#endif
#endif
#endif
#endif

#if 1 // F %= A
#ifdef AISSQUARE
#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                A0[k] = A1[k];
                F0[k] = E0[k] % A0[k];
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) {
            F1[k] = E1[k];
        }
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            GETDIV(A1[k]).solveTransposeInPlace(F1[k].transpose());
#else
            DIVIDEUSING(A1[k]);
            F1[k] %= MPART1(A1[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"E1 = "<<E1[k]<<std::endl;
                //std::cout<<"F1 = "<<F1[k]<<std::endl;
                //std::cout<<"F0 = "<<F0[k]<<std::endl;
                e9_reg += NormSq(F0[k]-F1[k]);
            }
            e9_reg = sqrt(e9_reg/nloops2);
            e9_reg /= Norm(F0[0]);
            //std::cout<<"e9_reg = "<<e9_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) F2[k] = E2[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ALGO == 1 && !defined(AISSQUARE)
            F2[k] %= MPART1(A2[k]);
#else
            GETDIV(A2[k]).solveTransposeInPlace(F2[k].transpose());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                e9_small += NormSq(F0[k]-F2[k]);
            }
            e9_small = sqrt(e9_small/nloops2);
            e9_small /= Norm(F0[0]);
            //std::cout<<"e9_small = "<<e9_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        for (int k=0; k<nloops2; ++k) F3[k] = E3[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            G3 = F3[k].adjoint();
#if ALGO == 1
            QR3[k] = A3[k];
            LAPNAME(gels) (
                LAPCM 
#ifdef TISCOMPLEX
                LAPCT,
#else
                LAPT,
#endif
                M,N,K,LP(QR3[k].ptr()),M,LP(G3.ptr()),M,
                LP(Work3.ptr()),lwork,&lapinfo);
#elif ALGO == 2
            QR3[k] = A3[k].adjoint();
            RT rcond = 1.e-6;
            int rank;
            for(int i=0;i<N;++i) P3i[i] = 0;
            LAPNAME(gelsy) (
                LAPCM M,N,K,LP(QR3[k].ptr()),M,LP(G3.ptr()),M,&(P3i[0]),
                rcond,&rank,LP(Work3.ptr()),lwork,
#ifdef TISCOMPLEX
                RWork3.ptr(),
#endif
                &lapinfo);
#endif
            F3[k].adjoint() = G3.rowRange(0,N);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                e9_blas += NormSq(F0[k]-F3[k]);
            }
            e9_blas = sqrt(e9_blas/nloops2);
            e9_blas /= Norm(F0[0]);
            //std::cout<<"e9_blas = "<<e9_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) F4[k] = E4[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F4[k].transpose() = A4[k].transpose() EIGENDIV .solve(F4[k].transpose());
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<K;++i) for(int j=0;j<N;++j)
                    e9_eigen += tmv::TMV_NORM(F4[k](i,j)-F0[k](i,j));
            }
            e9_eigen = sqrt(e9_eigen/nloops2);
            e9_eigen /= Norm(F0[0]);
            //std::cout<<"e9_eigen = "<<e9_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        //std::cout<<"start EIGENSMALL %="<<std::endl;
        for (int k=0; k<nloops2x; ++k) F5[k] = E5[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            F5[k].transpose() = A5[k].transpose() EIGENDIV .solve(F5[k].transpose());
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e9_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                for(int i=0;i<K;++i) for(int j=0;j<N;++j)
                    e9_smalleigen += tmv::TMV_NORM(F5[k](i,j)-F0[k](i,j));
            }
            e9_smalleigen = sqrt(e9_smalleigen/nloops2);
            e9_smalleigen /= Norm(F0[0]);
            //std::cout<<"e9_smalleigen = "<<e9_smalleigen<<std::endl;
        }
#endif
#endif
#endif
#endif

#if 1 // F = E % AT
#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                // Solve FAt = E
                // FAtA = EA
                // F = (EA) % (AtA)
                A0[k] = A1[k];
                F0[k] = (E0[k]*A0[k]) % (A0[k].adjoint()*A0[k]);
            }
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            GETDIV(A1[k].adjoint()).solveTranspose(
                E1[k].transpose(),F1[k].transpose());
#elif ALGO == 1 && !defined(AISSQUARE)
            F1[k] = E1[k] % MPART1(A1[k].adjoint());
#else
            B1[k] = A1[k].adjoint();
            DIVIDEUSING(B1[k]);
            F1[k] = E1[k] % MPART1(B1[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"E1 = "<<E1[k]<<std::endl;
                //std::cout<<"F1 = "<<F1[k]<<std::endl;
                //std::cout<<"F0 = "<<F0[k]<<std::endl;
                e10_reg += NormSq(F0[k]-F1[k]);
            }
            e10_reg = sqrt(e10_reg/nloops2);
            e10_reg /= Norm(F0[0]);
            //std::cout<<"e10_reg = "<<e10_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ALGO == 1 && !defined(AISSQUARE)
            F2[k] = E2[k] % MPART1(A2[k].adjoint());
#else
            GETDIV(A2[k].adjoint()).solveTranspose(
                E2[k].transpose(),F2[k].transpose());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                e10_small += NormSq(F0[k]-F2[k]);
            }
            e10_small = sqrt(e10_small/nloops2);
            e10_small /= Norm(F0[0]);
            //std::cout<<"e10_small = "<<e10_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            QR3[k] = A3[k];
            G3 = E3[k].adjoint();
#if ALGO == 1
            LAPNAME(gels) (
                LAPCM LAPNT,M,N,K,LP(QR3[k].ptr()),M,LP(G3.ptr()),M,
                LP(Work3.ptr()),lwork,&lapinfo);
#elif ALGO == 2
            RT rcond = 1.e-6;
            int rank;
            for(int i=0;i<N;++i) P3i[i] = 0;
            LAPNAME(gelsy) (
                LAPCM M,N,K,LP(QR3[k].ptr()),M,LP(G3.ptr()),M,&(P3i[0]),
                rcond,&rank,LP(Work3.ptr()),lwork,
#ifdef TISCOMPLEX
                RWork3.ptr(),
#endif
                &lapinfo);
#endif
            F3[k].adjoint() = G3.rowRange(0,N);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                e10_blas += NormSq(F0[k]-F3[k]);
            }
            e10_blas = sqrt(e10_blas/nloops2);
            e10_blas /= Norm(F0[0]);
            //std::cout<<"e10_blas = "<<e10_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // F = E % At
            // F = E * At^-1
            // F At = E
            // A Ft = Et
            F4[k].transpose() = A4[k] EIGENDIV .solve(E4[k].adjoint());
            F4[k] = F4[k].conjugate();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"A4 = "<<A4[k]<<std::endl;
                //std::cout<<"E4 = "<<E4[k]<<std::endl;
                //std::cout<<"F4 = "<<F4[k]<<std::endl;
                //std::cout<<"F0 = "<<F0[k]<<std::endl;
                for(int i=0;i<K;++i) for(int j=0;j<N;++j)
                    e10_eigen += tmv::TMV_NORM(F4[k](i,j)-F0[k](i,j));
            }
            e10_eigen = sqrt(e10_eigen/nloops2);
            e10_eigen /= Norm(F0[0]);
            //std::cout<<"e10_eigen = "<<e10_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        //std::cout<<"start EIGENSMALL %"<<std::endl;
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            F5[k].transpose() = A5[k] EIGENDIV .solve(E5[k].adjoint());
            F5[k] = F5[k].conjugate();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e10_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                for(int i=0;i<K;++i) for(int j=0;j<N;++j)
                    e10_smalleigen += tmv::TMV_NORM(F5[k](i,j)-F0[k](i,j));
            }
            e10_smalleigen = sqrt(e10_smalleigen/nloops2);
            e10_smalleigen /= Norm(F0[0]);
            //std::cout<<"e10_smalleigen = "<<e10_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // d = A.det()
#ifdef AISSQUARE
#ifdef ERRORCHECK
        for (int k=0; k<nloops2; ++k) {
            A0[k] = A1[k];
            d0[k] = MPART1(A0[k]).det();
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            d1[k] = GETDIV(A1[k]).det();
#else
            DIVIDEUSING(A1[k]);
            d1[k] = MPART1(A1[k]).det();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"d1 = "<<d1[k]<<std::endl;
                //std::cout<<"d0 = "<<d0[k]<<std::endl;
                e12_reg += tmv::TMV_NORM(d0[k]-d1[k]);
            }
            e12_reg = sqrt(e12_reg/nloops2);
            e12_reg /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e12_reg = "<<e12_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ALGO == 1 && !defined(AISSQUARE)
            d2[k] = MPART1(A2[k]).det();
#else
            d2[k] = GETDIV(A2[k]).det();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                e12_small += tmv::TMV_NORM(d0[k]-d2[k]);
            }
            e12_small = sqrt(e12_small/nloops2);
            e12_small /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e12_small = "<<e12_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            QR3[k] = A3[k];
#if ALGO == 1
            LAPNAME(geqrf) (
                LAPCM M,N,LP(QR3[k].ptr()),M,LP(Beta3[k].ptr()),
                LP(Work3.ptr()),lwork,&lapinfo);
#elif ALGO == 2
            for(int i=0;i<N;++i) P3i[i] = 0;
            LAPNAME(geqp3) (
                LAPCM M,N,LP(QR3[k].ptr()),M,&(P3i[0]),LP(Beta3[k].ptr()),
                LP(Work3.ptr()),lwork,
#ifdef TISCOMPLEX
                RWork3.ptr(),
#endif
                &lapinfo);
            ConvertIndexToPermute(N,&(P3i[0]),&(P3[k][0]));
#endif
            d3[k] = T(1);
            for(int i=0;i<N;++i) d3[k] *= QR3[k](i,i);
            for(int i=0;i<N;++i) if (Beta3[k][i] != T(0)) {
#ifdef TISCOMPLEX
                T b = Beta3[k][i];
                d3[k] *= -(b*b)/std::norm(b);
#else
                d3[k] = -d3[k];
#endif
            }
#if ALGO == 2
            for(int i=0;i<N;++i) if (P3[k][i] != i) d3[k] = -d3[k];
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                e12_blas += tmv::TMV_NORM(d0[k]-d3[k]);
            }
            e12_blas = sqrt(e12_blas/nloops2);
            e12_blas /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e12_blas = "<<e12_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            d4[k] = A4[k] EIGENDIV .absDeterminant();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e12_eigen += tmv::TMV_NORM(std::abs(d0[k])-d4[k]);
            }
            e12_eigen = sqrt(e12_eigen/nloops2);
            e12_eigen /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e12_eigen = "<<e12_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        //std::cout<<"start EIGENSMALL det"<<std::endl;
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            d5[k] = A5[k] EIGENDIV .absDeterminant();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e12_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                e12_smalleigen += tmv::TMV_NORM(std::abs(d0[k])-d5[k]);
            }
            e12_smalleigen = sqrt(e12_smalleigen/nloops2);
            e12_smalleigen /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e12_smalleigen = "<<e12_smalleigen<<std::endl;
        }
#endif
#endif
#endif
#endif

#ifdef SHOW_DOTS
        std::cout<<"."; std::cout.flush();
#endif

#ifdef DOREG
#ifdef TMV_NO_LIB
        for(int i=0;i<nloops2;++i) if (QR1[i]) delete QR1[i];
#endif
#endif
#ifdef DOSMALL
        for(int i=0;i<nloops2;++i) if (QR2[i]) delete QR2[i];
#endif
#ifdef DOEIGEN
        for(int i=0;i<nloops2;++i) if (QR4[i]) delete QR4[i];
#endif
#ifdef DOEIGENSMALL
        for(int i=0;i<nloops2x;++i) if (QR5[i]) delete QR5[i];
#endif
    }
    std::cout<<"\n";

#if ALGO==1
    const std::string algo_text = "QR ";
#elif ALGO==2
    const std::string algo_text = "QRP";
#elif ALGO==3
    const std::string algo_text = "USV";
#endif

#ifdef DO_SEPARATE_DECOMP
    std::cout<<"A -> "<<algo_text<<"               "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"No decomposition:\n";
    std::cout<<"D = C / A              "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"B = A.inverse()        "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
#ifdef AISSQUARE
    std::cout<<"D /= A                 "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"d = A.det()            "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
    std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
#endif
    std::cout<<"Including decomposition:\n";
#endif
    std::cout<<"D = C / A              "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"B = A.inverse()        "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"F = E % AT             "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
    std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
#ifdef AISSQUARE
    std::cout<<"D /= A                 "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
    std::cout<<"B = B.inverse()        "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
    std::cout<<"F %= A                 "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
    std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
    std::cout<<"d = A.det()            "<<t12_reg<<"  "<<t12_small<<"  "<<t12_blas;
    std::cout<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
#ifdef DO_SEPARATE_DECOMP
    std::cout<<"A -> "<<algo_text<<"               "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"D = C / A ND           "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"B = A.inverse() ND     "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
#ifdef AISSQUARE
    std::cout<<"D /= A ND              "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"d = A.det() ND         "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
    std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
#endif
#endif
    std::cout<<"D = C / A              "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
    std::cout<<"B = A.inverse()        "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
    std::cout<<"F = E % AT             "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
    std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
#ifdef AISSQUARE
    std::cout<<"D /= A                 "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
    std::cout<<"B = B.inverse()        "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
    std::cout<<"F %= A                 "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
    std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
    std::cout<<"d = A.det()            "<<e12_reg<<"  "<<e12_small<<"  "<<e12_blas;
    std::cout<<"  "<<e12_eigen<<"  "<<e12_smalleigen<<std::endl;
#endif
    std::cout<<"\n\n";
#endif
}
#endif

#ifdef DO_R
static void QR_R(
    const std::vector<tmv::Matrix<T,tmv::RowMajor> >& A1,
    std::vector<tmv::Matrix<T> >& B1,
    const std::vector<tmv::Matrix<T> >& C1,
    std::vector<tmv::Matrix<T> >& D1,
    const std::vector<tmv::Matrix<T> >& E1,
    std::vector<tmv::Matrix<T> >& F1)
#endif


#ifdef SIMPLE_VALUES
#ifdef TISCOMPLEX
#define RAND ( T(i+1,j+1) )
#else
#define RAND ( T(i+j+1) )
#endif
#else
#ifdef TISCOMPLEX
#define RAND1 ( T(rand(),rand()) / RT(RAND_MAX) )
#else
#define RAND1 ( T(rand()) / RT(RAND_MAX) )
#endif
#define RAND  ( RAND1 * RT(1000) + RAND1 - RT(500) ) / RT(1000)
#endif


int main() try
{
    srand(572924738);

#ifdef MEMDEBUG
    atexit(&DumpUnfreed);
#endif

    std::vector<tmv::Matrix<T> > AC(nloops2,tmv::Matrix<T>(M,N));
    std::vector<tmv::Matrix<T,tmv::RowMajor> > AR(nloops2,tmv::Matrix<T,tmv::RowMajor>(M,N));
    std::vector<tmv::Matrix<T> > B(nloops2,tmv::Matrix<T>(N,M));
    std::vector<tmv::Matrix<T> > C(nloops2,tmv::Matrix<T>(M,K));
    std::vector<tmv::Matrix<T> > D(nloops2,tmv::Matrix<T>(N,K));
    std::vector<tmv::Matrix<T> > E(nloops2,tmv::Matrix<T>(K,M));
    std::vector<tmv::Matrix<T> > F(nloops2,tmv::Matrix<T>(K,N));
    for(int k=0;k<nloops2;++k) {
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
            AC[k](i,j) = RAND;
#if 1
            if (i+2*j == N || i+2*j == 2*N) {
                AC[k](i,j) *= T(50);
            }
#endif
        }
        AC[k].diag().addToAll(100.);
#ifdef SMALL_OFFDIAG
        AC[k].upperTri().offDiag() /= 10000.;
        AC[k].rowRange(0,N).lowerTri().offDiag() /= 10000.;
        AC[k].rowRange(N,M) /= 10000.;
#endif
        for(int i=0;i<M;++i) for(int j=0;j<K;++j) {
            C[k](i,j) = RAND;
            E[k](j,i) = RAND;
        }
        B[k].setZero();
        D[k].setZero();
        F[k].setZero();
        AR[k] = AC[k];
        //AC[k].col(0).setZero();
        //AC[k](0,0) = 1.e-20;
        //AC[k](0,1) = 3.e-20;
    }

    std::cout<<"M,N,K = "<<M<<" , "<<N<<" , "<<K<<std::endl;
    std::cout<<"nloops = "<<nloops1<<" x "<<nloops2;
    std::cout<<" = "<<nloops1*nloops2<<std::endl;
#ifdef DOEIGENSMALL
    if (nloops2x == 0)
        std::cout<<"Matrix is too big for the stack, so no \"Eigen Known\" tests.\n";
#endif
    std::cout<<"\nTime for:               ";
    std::cout<<" TMV    TMV Small  LAPACK   Eigen Eigen Known"<<std::endl;

#ifdef STRICTQRP
    tmv::UseStrictQRP();
#endif

#ifdef DO_C
    QR_C(AC,B,C,D,E,F);
#endif

    for(int k=0;k<nloops2;++k) { D[k].setZero(); F[k].setZero(); }

#ifdef DO_R
    QR_R(AR,B,C,D,E,F);
#endif

    return 0;
}
#if 0
catch (int) {}
#else
catch (tmv::Error& e) {
    std::cout<<"Caught error "<<e<<std::endl;
    return 1;
}
#endif
