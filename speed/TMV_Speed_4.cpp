//#define PRINTALGO_LU
//#define PRINTALGO_DIVU
//#define PRINTALGO_DIVU_OMP

//#undef NDEBUG

//#define TMV_INLINE
#include "TMV.h"

// How big do you want the matrices to be?
// (The matrix being inverted is NxN.  
//  K is the other dimension of the matrix being solved.)
const int N = 7;
const int K = 77;

// Define the type to use:
//#define TISFLOAT
//#define TISCOMPLEX

// Define the parts of the matrices to use
// 1 = all (rectangle matrix)
// 2 = uppertri
// 3 = unituppertri
// 4 = lowertri
// 5 = unitlowertri
#define PART1 1

// Define whether you want to include the error checks.
#define ERRORCHECK

// Define which versions you want to test:
#define DOREG
#define DOSMALL
#define DOLAP
#define DOEIGEN
//#define DOEIGENSMALL

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
//#define SIMPLE_VALUES

// Set up the target number of operations and memory to use for testing
const long long targetnmegaflops(1000); // in 10^6 real ops
const long long targetmem(10000); // in Kbytes

// Include the LAPACK library you want to test TMV against:
#if 0
#include "mkl.h"
#define CBLAS
#define CLAPACK
#define XLAP
#else
extern "C" {
#ifdef TMV_V063
#include "fblas.h"
#include "flapack.h"
#else
#include "util/fblas.h"
#include "util/flapack.h"
#endif
#define FBLAS
#define FLAP
#define XLAP
}
#endif


// The rest of this is based on the above values.
// You shouldn't need to change anything else here.

#include <complex>

#ifdef TISFLOAT
typedef float RT;
#else
typedef double RT;
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
const long long nloops2 = (
    (((long long)(sizeof(T)))*(N*N+N*K)/1000 > targetmem ? 1 :
     targetmem*1000 / (((long long)(sizeof(T)))*(N*N+N*K))));
const long long nloops1 = (
    (nloops2*N*N*(N+K)*XFOUR/1000000 > targetnmegaflops ? 1 :
     targetnmegaflops*1000000 / (nloops2*N*N*(N+K)) / XFOUR ));
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
#if !EIGEN_VERSION_AT_LEAST(2,91,0)
#include "Eigen/Array"
#endif
#include "Eigen/LU"

#endif

#if (PART1 == 1) // Full rectangle
#define INPART1(i,j) true
#define UNITPART1(i,j) false
#define MPART1(m) m
#define EPART1(m) m
#define EPART1t(m) m
#define COL_LEN1(j) M
#define COL_START1(j) 0
#define COL_END1(j) M
#define ROW_LEN1(i) N
#define ROW_START1(i) 0
#define ROW_END1(i) N

#elif (PART1 == 2) // UpperTri
#define INPART1(i,j) (i<=j)
#define UNITPART1(i,j) false
#define MPART1(m) m.upperTri()
#ifdef EIGEN_VERSION_AT_LEAST
#if EIGEN_VERSION_AT_LEAST(2,91,0)
#define EPART1(m) m.triangularView<Eigen::Upper>()
#define EPART1t(m) m.triangularView<Eigen::Lower>()
#else
#define EPART1(m) m.part<Eigen::UpperTriangular>()
#define EPART1t(m) m.part<Eigen::LowerTriangular>()
#endif
#endif
#define COL_LEN1(j) j+1
#define COL_START1(j) 0
#define COL_END1(j) j+1
#define ROW_LEN1(i) N-i
#define ROW_START1(i) i
#define ROW_END1(i) N
#define SMALL_OFFDIAG

#elif (PART1 == 3) // UnitUpperTri
#define INPART1(i,j) (i<j)
#define UNITPART1(i,j) (i==j)
#define MPART1(m) m.unitUpperTri()
#ifdef EIGEN_VERSION_AT_LEAST
#if EIGEN_VERSION_AT_LEAST(2,91,0)
#define EPART1(m) m.triangularView<Eigen::UnitUpper>()
#define EPART1t(m) m.triangularView<Eigen::UnitLower>()
#else
#define EPART1(m) m.part<Eigen::UnitUpperTriangular>()
#define EPART1t(m) m.part<Eigen::UnitLowerTriangular>()
#endif
#endif
#define COL_LEN1(j) j
#define COL_START1(j) 0
#define COL_END1(j) j
#define ROW_LEN1(i) N-i-1
#define ROW_START1(i) i+1
#define ROW_END1(i) N
#define SMALL_OFFDIAG

#elif (PART1 == 4) // LowerTri
#define INPART1(i,j) (i>=j)
#define UNITPART1(i,j) false
#define MPART1(m) m.lowerTri()
#ifdef EIGEN_VERSION_AT_LEAST
#if EIGEN_VERSION_AT_LEAST(2,91,0)
#define EPART1(m) m.triangularView<Eigen::Lower>()
#define EPART1t(m) m.triangularView<Eigen::Upper>()
#else
#define EPART1(m) m.part<Eigen::LowerTriangular>()
#define EPART1t(m) m.part<Eigen::UpperTriangular>()
#endif
#endif
#define COL_LEN1(j) N-j
#define COL_START1(j) j
#define COL_END1(j) N
#define ROW_LEN1(i) i+1
#define ROW_START1(i) 0
#define ROW_END1(i) i+1
#define SMALL_OFFDIAG

#elif (PART1 == 5) // UnitLowerTri
#define INPART1(i,j) (i>j)
#define UNITPART1(i,j) (i==j)
#define MPART1(m) m.unitLowerTri()
#ifdef EIGEN_VERSION_AT_LEAST
#if EIGEN_VERSION_AT_LEAST(2,91,0)
#define EPART1(m) m.triangularView<Eigen::UnitLower>()
#define EPART1t(m) m.triangularView<Eigen::UnitUpper>()
#else
#define EPART1(m) m.part<Eigen::UnitLowerTriangular>()
#define EPART1t(m) m.part<Eigen::UnitUpperTriangular>()
#endif
#endif
#define COL_LEN1(j) N-j-1
#define COL_START1(j) j+1
#define COL_END1(j) N
#define ROW_LEN1(i) i
#define ROW_START1(i) 0
#define ROW_END1(i) i
#define SMALL_OFFDIAG

#endif


#ifdef DOSMALL
#ifdef TMV_V063
#include "TMV_Small.h"
#endif
#endif

#ifdef DOLAP

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
UNKNOWN BLAS definition....  // Give compile error
#define BLASNAME1(c,x) c ## x
#define BLASINAME1(c,x) i ## c ## x
#endif
#define BLASNAME2(c,x) BLASNAME1(c,x)
#define BLASINAME2(c,x) BLASINAME1(c,x)

#define BLASNAME(x) BLASNAME2(BLAS_LETTER,x)
#define BLASINAME(x) BLASINAME2(BLAS_LETTER,x)

#if (PART1 == 2 || PART1 == 3)
#define BLASUPLO1 BLASUpper
#define BLASUPLO1x BLASLower
#endif
#if (PART2 == 2 || PART2 == 3)
#define BLASUPLO2 BLASUpper
#define BLASUPLO2x BLASLower
#endif
#if (PART1 == 4 || PART1 == 5)
#define BLASUPLO1 BLASLower
#define BLASUPLO1x BLASUpper
#endif
#if (PART2 == 4 || PART2 == 5)
#define BLASUPLO2 BLASLower
#define BLASUPLO2x BLASUpper
#endif
#if (PART1 == 3 || PART1 == 5)
#define BLASDIAG1 BLASUnit
#endif
#if (PART1 == 2 || PART1 == 4)
#define BLASDIAG1 BLASNonUnit
#endif
#if (PART2 == 3 || PART2 == 5)
#define BLASDIAG2 BLASUnit
#endif
#if (PART2 == 2 || PART2 == 4)
#define BLASDIAG2 BLASNonUnit
#endif

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
#define LAPNAME1(c,x) c ## x ## _
#define LAPCM 
RT* lap_rwork = new RT[lap_worksize];
T* lap_work = new T[lap_worksize];
#define LAPWORK lap_work
#define LAPRWORK lap_rwork
#define LAPNT LapCh_N
#define LAPT LapCh_T
#define LAPCT LapCh_C
#else
UNKNOWN LAP definition....  // Give compile error
#define LAPNAME1(c,x) c ## x
#endif

int lapinfo;

#define LAPNAME2(c,x) LAPNAME1(c,x)
#define LAPNAME(x) LAPNAME2(BLAS_LETTER,x)

#endif // DOLAP

#ifdef DOEIGEN
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
    (N*N*sizeof(T) < 256 * 1024 && 
     N*K*sizeof(T) < 256 * 1024 && 
     N!=10000 && K!=10000) ? 
    nloops2 : 0 );

#define EIGENSMA Eigen::Matrix<T,N,N>
#define EIGENSMB Eigen::Matrix<T,N,N>
#define EIGENSMC Eigen::Matrix<T,N,K>
#define EIGENSMD Eigen::Matrix<T,N,K>
#define EIGENSME Eigen::Matrix<T,K,N>
#define EIGENSMF Eigen::Matrix<T,K,N>

#endif

static void ClearCache()
{
    tmv::Vector<double> X(1000000,8.);
    tmv::Vector<double> Y(1000000,8.);
    tmv::Vector<double> Z(1000000);
    Z = X + Y;
    if (Norm(Z) < 5.) exit(1);
}

#ifdef DO_C
static void LU_C(
    const std::vector<tmv::Matrix<T> >& A1,
    const std::vector<tmv::Matrix<T> >& C1,
    std::vector<tmv::Matrix<T> >& D1,
    const std::vector<tmv::Matrix<T> >& E1,
    std::vector<tmv::Matrix<T> >& F1)
{
    std::vector<tmv::Matrix<T> > LU1 = A1;
    std::vector<std::vector<int> > P1(nloops2);
    for(int k=0;k<nloops2;++k) {
        P1[k].resize(N);
    }

#ifdef ERRORCHECK
    std::vector<tmv::Matrix<T> > A0 = A1;
    std::vector<tmv::Matrix<T> > C0 = C1;
    std::vector<tmv::Matrix<T> > E0 = E1;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallMatrix<T,N,N> > A2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,K> > C2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,K> > D2(nloops2);
    std::vector<tmv::SmallMatrix<T,K,N> > E2(nloops2);
    std::vector<tmv::SmallMatrix<T,K,N> > F2(nloops2);
#ifdef TMV_V063
    std::vector<tmv::SmallMatrix<T,N,N> > LU2(nloops2);
    std::vector<std::vector<int> > P2(nloops2);
#else
    std::vector<tmv::LUD<tmv::SmallMatrix<T,N,N> >*> LU2(nloops2);
#endif

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; C2[k] = C1[k]; E2[k] = E1[k];
#ifdef TMV_V063
        P2[k].resize(N);
#endif
    }
#endif

#ifdef DOLAP
    std::vector<tmv::Matrix<T> > A3 = A1;
    std::vector<tmv::Matrix<T> > C3 = C1;
    std::vector<tmv::Matrix<T> > D3 = D1;
    std::vector<tmv::Matrix<T> > E3 = E1;
    std::vector<tmv::Matrix<T> > F3 = F1;
    std::vector<tmv::Matrix<T> > LU3 = A1;
    std::vector<std::vector<int> > P3(nloops2);
    for(int k=0;k<nloops2;++k) {
        P3[k].resize(N);
    }
#endif

#ifdef DOEIGEN
    std::vector<EIGENM,ALLOC(EIGENM) > A4(nloops2,EIGENM(N,N));
    std::vector<EIGENM,ALLOC(EIGENM) > C4(nloops2,EIGENM(N,K));
    std::vector<EIGENM,ALLOC(EIGENM) > D4(nloops2,EIGENM(N,K));
    std::vector<EIGENM,ALLOC(EIGENM) > E4(nloops2,EIGENM(K,N));
    std::vector<EIGENM,ALLOC(EIGENM) > F4(nloops2,EIGENM(K,N));
#if EIGEN_VERSION_AT_LEAST(2,91,0)
    std::vector<Eigen::PartialPivLU<EIGENM>*,ALLOC(Eigen::PartialPivLU<EIGENM>*) > LU4(nloops2);
#else
    std::vector<Eigen::LU<EIGENM>*,ALLOC(Eigen::LU<EIGENM>*) > LU4(nloops2);
#endif
    for(int k=0;k<nloops2;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<N;++i) A4[k](i,j) = A1[k](i,j);
        for(int j=0;j<K;++j) for(int i=0;i<N;++i) C4[k](i,j) = C1[k](i,j);
        for(int j=0;j<K;++j) for(int i=0;i<N;++i) D4[k](i,j) = D1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<K;++i) E4[k](i,j) = E1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<K;++i) F4[k](i,j) = F1[k](i,j);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<EIGENSMA,ALLOC(EIGENSMA) > A5;
    std::vector<EIGENSMC,ALLOC(EIGENSMC) > C5;
    std::vector<EIGENSMD,ALLOC(EIGENSMD) > D5;
    std::vector<EIGENSME,ALLOC(EIGENSME) > E5;
    std::vector<EIGENSMF,ALLOC(EIGENSMF) > F5;
#if EIGEN_VERSION_AT_LEAST(2,91,0)
    std::vector<Eigen::PartialPivLU<EIGENSMA>*,ALLOC(Eigen::PartialPivLU<EIGENSMA>*) > LU5;
#else
    std::vector<Eigen::LU<EIGENSMA>*,ALLOC(Eigen::LU<EIGENSMA>*) > LU5;
#endif
    if (nloops2x) { 
        A5.resize(nloops2x); 
        C5.resize(nloops2x);
        D5.resize(nloops2x);
        E5.resize(nloops2x); 
        F5.resize(nloops2x); 
        LU5.resize(nloops2x); 
    }
    for(int k=0;k<nloops2x;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<N;++i) A5[k](i,j) = A1[k](i,j);
        for(int j=0;j<K;++j) for(int i=0;i<N;++i) C5[k](i,j) = C1[k](i,j);
        for(int j=0;j<K;++j) for(int i=0;i<N;++i) D5[k](i,j) = D1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<K;++i) E5[k](i,j) = E1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<K;++i) F5[k](i,j) = F1[k](i,j);
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

#if 1 // A -> PLU
#if PART1 == 1
        ClearCache();
        for (int k=0; k<nloops2; ++k) A1[k].unsetDiv();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_V063
            LU1[k] = A1[k];
            LU_Decompose(LU1[k].view(),&P1[k][0]);
#else
            //LU1[k] = A1[k];
            //LU_Decompose(LU1[k],&P1[k][0]);
            A1[k].saveDiv();
            A1[k].divideUsing(tmv::LU);
            A1[k].setDiv();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
#ifdef TMV_V063
                A0[k] = LU1[k].lowerTri(tmv::UnitDiag) * LU1[k].upperTri();
                A0[k].reversePermuteRows(&P1[k][0]);
#else
                //A0[k] = LU1[k].unitLowerTri() * LU1[k].upperTri();
                //A0[k].reversePermuteRows(&P1[k][0]);
                A0[k] = A1[k].lud().getL() * A1[k].lud().getU();
                A0[k].reversePermuteRows(A1[k].lud().getP());
#endif
                e1_reg += NormSq(A0[k]-A1[k]);
            }
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_V063
            LU2[k] = A2[k];
            T sign=0;
            DoLUD(LU2[k],&P2[k][0],sign);
#else
            LU2[k] = new tmv::LUD<tmv::SmallMatrix<T,N,N> >(A2[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) {
#ifdef TMV_V063
                A0[k] = LU2[k].lowerTri(tmv::UnitDiag) * LU2[k].upperTri();
                A0[k].reversePermuteRows(&P2[k][0]);
#else
                A0[k] = LU2[k]->getL() * LU2[k]->getU();
                A0[k].reversePermuteRows(LU2[k]->getP());
#endif
                e1_small += NormSq(A2[k]-A0[k]);
            }
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            LU3[k] = A3[k];
            LAPNAME(getrf) (
                LAPCM N,N,LP(LU3[k].ptr()),N,&P3[k][0],&lapinfo);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<N;++i) --P3[k][i];
#ifdef TMV_V063
                A0[k] = LU3[k].lowerTri(tmv::UnitDiag) * LU3[k].upperTri();
#else
                A0[k] = LU3[k].unitLowerTri() * LU3[k].upperTri();
#endif
                A0[k].reversePermuteRows(&P3[k][0]);
                for(int i=0;i<N;++i) ++P3[k][i];
                e1_blas += NormSq(A3[k]-A0[k]);
            }
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if EIGEN_VERSION_AT_LEAST(2,91,0)
            LU4[k] = new Eigen::PartialPivLU<EIGENM>(A4[k]);
#else
            LU4[k] = new Eigen::LU<EIGENM>(A4[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
#if EIGEN_VERSION_AT_LEAST(2,91,0)
                EIGENM temp1 = 
                    EIGENM(LU4[k]->matrixLU().triangularView<Eigen::UnitLower>()) *
                    EIGENM(LU4[k]->matrixLU().triangularView<Eigen::Upper>());
                temp1 = LU4[k]->permutationP().inverse() * temp1;
#else
                EIGENM temp1 = 
                    LU4[k]->matrixLU().part<Eigen::UnitLowerTriangular>() *
                    LU4[k]->matrixLU().part<Eigen::UpperTriangular>();
                EIGENM temp2(temp1.rows(),temp1.cols());
                for(int i=0;i<temp1.rows();++i) 
                    temp2.row(i) = temp1.row(LU4[k]->permutationP()(i));
                for(int j=0;j<temp1.cols();++j) 
                    temp1.col(LU4[k]->permutationQ()(j)) = temp2.col(j);
#endif
                e1_eigen += (A4[k]-temp1).squaredNorm();
            }
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if EIGEN_VERSION_AT_LEAST(2,91,0)
            LU5[k] = new Eigen::PartialPivLU<EIGENSMA>(A5[k]);
#else
            LU5[k] = new Eigen::LU<EIGENSMA>(A5[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
#if EIGEN_VERSION_AT_LEAST(2,91,0)
                EIGENM temp1 = 
                    EIGENM(LU5[k]->matrixLU().triangularView<Eigen::UnitLower>()) *
                    EIGENM(LU5[k]->matrixLU().triangularView<Eigen::Upper>());
                temp1 = LU4[k]->permutationP().inverse() * temp1;
#else
                EIGENM temp1 = 
                    LU5[k]->matrixLU().part<Eigen::UnitLowerTriangular>() *
                    LU5[k]->matrixLU().part<Eigen::UpperTriangular>();
                EIGENM temp2(temp1.rows(),temp1.cols());
                for(int i=0;i<temp1.rows();++i) 
                    temp2.row(i) = temp1.row(LU5[k]->permutationP()(i));
                for(int j=0;j<temp1.cols();++j) 
                    temp1.col(LU5[k]->permutationQ()(j)) = temp2.col(j);
#endif
                e1_smalleigen += (A5[k]-temp1).squaredNorm();
            }
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // D /= A
        ClearCache();

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) D1[k] = C1[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D1[k] /= MPART1(A1[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A1[k]) * D1[k];
                e2_reg += NormSq(C0[k]-C1[k]);
            }
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_V063
            D2[k] /= MPART1(A2[k]);
#else
#if PART1 == 1
            LU2[k]->solveInPlace(D2[k]);
#else
            D2[k] /= MPART1(A2[k]);
#endif
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A2[k]) * D2[k];
                e2_small += NormSq(C0[k]-C2[k]);
            }
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOLAP
        for (int k=0; k<nloops2; ++k) D3[k] = C3[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LAPNAME(getrs) (
                LAPCM LAPNT, N,K,LP(LU3[k].ptr()),N,&P3[k][0],
                LP(D3[k].ptr()),N,&lapinfo);
#elif PART1 >= 2 && PART1 <= 5
            BLASNAME(trsm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
                            N,K,one,BP(A3[k].cptr()),N,BP(D3[k].ptr()),N
                            BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A3[k]) * D3[k];
                e2_blas += NormSq(C0[k]-C3[k]);
            }
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
#if PART1 >= 2 && PART1 <= 5
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];
#endif
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
#if EIGEN_VERSION_AT_LEAST(2,91,0)
            D4[k] = LU4[k]->solve(C4[k]);
#else
            LU4[k]->solve(C4[k],&D4[k]);
#endif
#elif PART1 >= 2 && PART1 <= 5
            EPART1(A4[k]).solveTriangularInPlace(D4[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k]) * D4[k];
                e2_eigen += (C4[k]-temp1).squaredNorm();
            }
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
#if PART1 >= 2 && PART1 <= 5
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];
#endif
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
#if EIGEN_VERSION_AT_LEAST(2,91,0)
            D5[k] = LU5[k]->solve(C5[k]);
            //D5[k] = A5[k].lu().solve(C5[k]);
#else
            //LU5[k]->solve(C5[k],&D5[k]);
            A5[k].lu().solve(C5[k],&D5[k]);
#endif
#elif PART1 >= 2 && PART1 <= 5
            EPART1(A5[k]).solveTriangularInPlace(D5[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k]) * D5[k];
                e2_smalleigen += (C5[k]-temp1).squaredNorm();
            }
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // F %= A
        ClearCache();

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) F1[k] = E1[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F1[k] %= MPART1(A1[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F1[k] * MPART1(A1[k]);
                e3_reg += NormSq(E0[k]-E1[k]);
            }
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) F2[k] = E2[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_V063
            F2[k] %= MPART1(A2[k]);
#else
#if PART1 == 1
            tmv::SmallMatrix<T,K,N>::transpose_type F2t = 
                F2[k].transpose();
            LU2[k]->solveTransposeInPlace(F2t);
#else
            F2[k] %= MPART1(A2[k]);
#endif
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F2[k] * MPART1(A2[k]);
                e3_small += NormSq(E0[k]-E2[k]);
            }
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOLAP
        for (int k=0; k<nloops2; ++k) F3[k] = E3[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            // This only works for RowMajor F3, so need a temporary.
            // Really, I should use the BLAS commands to effect this
            // copy, but it's so tedious...
            tmv::Matrix<T,tmv::RowMajor> temp = F3[k];
            LAPNAME(getrs) (
                LAPCM LAPT,N,K,LP(LU3[k].ptr()),N,&P3[k][0],
                LP(temp.ptr()),N,&lapinfo);
            F3[k] = temp;
#elif PART1 >= 2 && PART1 <= 5
            BLASNAME(trsm) (BLASCM BLASRight, BLASUPLO1, BLASNT, BLASDIAG1,
                            K,N,one,BP(A3[k].cptr()),N,BP(F3[k].ptr()),K
                            BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F3[k] * MPART1(A3[k]);
                e3_blas += NormSq(E0[k]-E3[k]);
            }
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
#if PART1 >= 2 && PART1 <= 5
        for (int k=0; k<nloops2; ++k) F4[k] = E4[k];
#endif
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            EIGENM temp(N,K);
            // I don't think Eigen has any way to do this once the LU
            // decomposition of A is already done.  So the comparison is 
            // a bit unfair, but apparently we are forced to do the 
            // decomposition over again on A.transpose().
#if EIGEN_VERSION_AT_LEAST(2,91,0)
            temp = A4[k].transpose().lu().solve(E4[k].transpose());
#else
            A4[k].transpose().lu().solve(E4[k].transpose(),&temp);
#endif
            F4[k] = temp.transpose();
#elif PART1 >= 2 && PART1 <= 5
            EPART1(A4[k]).transpose().solveTriangularInPlace(F4[k].transpose());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = F4[k] * EPART1(A4[k]);
                e3_eigen += (E4[k]-temp1).squaredNorm();
            }
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
#if PART1 >= 2 && PART1 <= 5
        for (int k=0; k<nloops2x; ++k) F5[k] = E5[k];
#endif
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            EIGENSMC temp;
#if EIGEN_VERSION_AT_LEAST(2,91,0)
            temp = A5[k].transpose().lu().solve(E5[k].transpose());
#else
            A5[k].transpose().lu().solve(E5[k].transpose(),&temp);
#endif
            F5[k] = temp.transpose();
#elif PART1 >= 2 && PART1 <= 5
            EPART1(A5[k]).transpose().solveTriangularInPlace(F5[k].transpose());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = F5[k] * EPART1(A5[k]);
                e3_smalleigen += (E5[k]-temp1).squaredNorm();
            }
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif
    }

#if PART1 == 1
    std::cout<<"A -> PLU               "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
#endif
    std::cout<<"D /= A                 "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"F %= A                 "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
#if PART1 == 1
    std::cout<<"A -> PLU               "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
#endif
    std::cout<<"D /= A                 "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"F %= A                 "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"\n\n";
#endif
}
#endif

#ifdef DO_R
static void LU_R(
    const std::vector<tmv::Matrix<T> >& B1,
    const std::vector<tmv::Matrix<T> >& C1,
    std::vector<tmv::Matrix<T> >& E1)
{
}
#endif


#ifdef SIMPLE_VALUES
#ifdef TISCOMPLEX
//#define RAND ( T(1,2) )
#define RAND ( T(i+1,j+1) )
#else
//#define RAND ( T(1) )
#define RAND ( T(i+j+1) )
#endif
#else
#ifdef TISCOMPLEX
#define RAND1 ( T(rand(),rand()) / RT(RAND_MAX) )
#else
#define RAND1 ( T(rand()) / RT(RAND_MAX) )
#endif
#define RAND  ( RAND1 * RT(1000) + RAND1 - RT(500) )
#endif


int main() try
{
    srand(572924738);

#ifdef MEMDEBUG
    atexit(&DumpUnfreed);
#endif

    std::vector<tmv::Matrix<T> > A(nloops2,tmv::Matrix<T>(N,N));
    std::vector<tmv::Matrix<T> > B(nloops2,tmv::Matrix<T>(N,N));
    std::vector<tmv::Matrix<T> > C(nloops2,tmv::Matrix<T>(N,K));
    std::vector<tmv::Matrix<T> > D(nloops2,tmv::Matrix<T>(N,K));
    std::vector<tmv::Matrix<T> > E(nloops2,tmv::Matrix<T>(K,N));
    std::vector<tmv::Matrix<T> > F(nloops2,tmv::Matrix<T>(K,N));
    for(int k=0;k<nloops2;++k) {
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) {
            A[k](i,j) = RAND;
            B[k](j,i) = RAND;
#if 0
            if (i+2*j == N || 3*i+2*j == 2*N) {
                A[k](i,j) *= T(50);
                B[k](j,i) *= T(60);
            }
#endif
        }
#ifdef SMALL_OFFDIAG
        A[k].upperTri().offDiag() /= 10000.;
        A[k].lowerTri().offDiag() /= 10000.;
        B[k].upperTri().offDiag() /= 10000.;
        B[k].lowerTri().offDiag() /= 10000.;
#endif
        for(int i=0;i<N;++i) for(int j=0;j<K;++j) {
            C[k](i,j) = RAND;
            E[k](j,i) = RAND;
        }
        D[k].setZero();
        F[k].setZero();
    }

    std::cout<<"N,K = "<<N<<" , "<<K<<std::endl;
    std::cout<<"nloops = "<<nloops1<<" x "<<nloops2;
    std::cout<<" = "<<nloops1*nloops2<<std::endl;
#ifdef DOIGENSMALL
    if (nloops2x == 0)
        std::cout<<"Matrix is too big for the stack, so no \"Eigen Known\" tests.\n";
#endif
    std::cout<<"\nTime for:               ";
    std::cout<<" TMV    TMV Small  LAPACK   Eigen Eigen Known"<<std::endl;

#ifdef DO_C
    LU_C(A,C,D,E,F);
#endif

    for(int k=0;k<nloops2;++k) E[k].setZero();

#ifdef DO_R
    LU_R(B,C,D,E,F);
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
