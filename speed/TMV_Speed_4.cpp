//#define PRINTALGO_LU
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

//#undef NDEBUG
#include <iostream>

#define TMV_NO_LIB
#include "TMV.h"

// How big do you want the matrices to be?
// (The matrix being inverted is NxN.  
//  K is the other dimension of the matrix being solved.)
const int N = 6;
const int K = 11;

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
#define DOEIGENSMALL

// Skip the separate decomposition and solution steps?
//#define NO_SEPARATE_DECOMP

// For large N, this can be useful to keep the condition of the matrix
// not too large.
//#define SMALL_OFFDIAG

// Define which batches of functions you want to test:
//#define DO_C
#define DO_R  

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
const long long targetnmegaflops(100); // in 10^6 real ops
const long long targetmem(100000); // in Kbytes

// Include the LAPACK library you want to test TMV against:
#if 0
#include "mkl.h"
#define MKL
#define CBLAS
#define CLAPACK
#define XLAP
#else
extern "C" {
#include "util/fblas.h"
#include "util/flapack.h"
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
const long long nloops2 = targetmem*1000 / ((2*N*N+4*N*K) * sizeof(T)) + 1;
const long long nloops1 = targetnmegaflops*1000000 / (N*N*(N+K)*nloops2*XFOUR) + 1;
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

#endif

#if (PART1 == 1) // Full rectangle
#define MPART1(m) m
#define EPART1(m) m
#define EPART1t(m) m
#define DO_SEPARATE_DECOMP

#elif (PART1 == 2) // UpperTri
#define MPART1(m) m.upperTri()
#define EPART1(m) m.triangularView<Eigen::Upper>()
#define EPART1t(m) m.triangularView<Eigen::Lower>()
#define SMALL_OFFDIAG

#elif (PART1 == 3) // UnitUpperTri
#define MPART1(m) m.unitUpperTri()
#define EPART1(m) m.triangularView<Eigen::UnitUpper>()
#define EPART1t(m) m.triangularView<Eigen::UnitLower>()
#define SMALL_OFFDIAG

#elif (PART1 == 4) // LowerTri
#define MPART1(m) m.lowerTri()
#define EPART1(m) m.triangularView<Eigen::Lower>()
#define EPART1t(m) m.triangularView<Eigen::Upper>()
#define SMALL_OFFDIAG

#elif (PART1 == 5) // UnitLowerTri
#define MPART1(m) m.unitLowerTri()
#define EPART1(m) m.triangularView<Eigen::UnitLower>()
#define EPART1t(m) m.triangularView<Eigen::UnitUpper>()
#define SMALL_OFFDIAG

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
    (N*N*sizeof(T) < 256 * 1024 && 
     N*K*sizeof(T) < 256 * 1024 && 
     N!=Eigen::Dynamic && K!=Eigen::Dynamic) ? 
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
    std::vector<tmv::Matrix<T> >& B1,
    const std::vector<tmv::Matrix<T> >& C1,
    std::vector<tmv::Matrix<T> >& D1,
    const std::vector<tmv::Matrix<T> >& E1,
    std::vector<tmv::Matrix<T> >& F1)
{
    //std::cout<<"Start LU_C.\n";
    //std::cout<<"A1.size = "<<A1.size()<<std::endl;
#if defined(TMV_NO_LIB)
    std::vector<tmv::LUD<tmv::Matrix<T> >*> LU1(nloops2,0);
#endif
    std::vector<T> d1(nloops2);

#ifdef ERRORCHECK
    std::vector<tmv::Matrix<T> > A0 = A1;
    std::vector<tmv::Matrix<T> > C0 = C1;
    std::vector<tmv::Matrix<T> > E0 = E1;
    std::vector<T> d0(nloops2);
    //std::cout<<"Made A0, etc."<<std::endl;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallMatrix<T,N,N> > A2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,N> > B2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,K> > C2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,K> > D2(nloops2);
    std::vector<tmv::SmallMatrix<T,K,N> > E2(nloops2);
    std::vector<tmv::SmallMatrix<T,K,N> > F2(nloops2);
    std::vector<T> d2(nloops2);
    std::vector<tmv::LUD<tmv::SmallMatrix<T,N,N> >*> LU2(nloops2,0);
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
    std::vector<tmv::Matrix<T> > LU3 = A1;
    std::vector<std::vector<int> > P3(nloops2);
    for(int k=0;k<nloops2;++k) {
        P3[k].resize(N);
    }
    //std::cout<<"Made A3, etc."<<std::endl;
#endif

#ifdef DOEIGEN
    std::vector<EIGENM,ALLOC(EIGENM) > A4(nloops2,EIGENM(N,N));
    std::vector<EIGENM,ALLOC(EIGENM) > B4(nloops2,EIGENM(N,N));
    std::vector<EIGENM,ALLOC(EIGENM) > C4(nloops2,EIGENM(N,K));
    std::vector<EIGENM,ALLOC(EIGENM) > D4(nloops2,EIGENM(N,K));
    std::vector<EIGENM,ALLOC(EIGENM) > E4(nloops2,EIGENM(K,N));
    std::vector<EIGENM,ALLOC(EIGENM) > F4(nloops2,EIGENM(K,N));
    std::vector<T> d4(nloops2);
    std::vector<Eigen::PartialPivLU<EIGENM>*,ALLOC(Eigen::PartialPivLU<EIGENM>*) > LU4(nloops2);
    for(int k=0;k<nloops2;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<N;++i) A4[k](i,j) = A1[k](i,j);
        for(int j=0;j<K;++j) for(int i=0;i<N;++i) C4[k](i,j) = C1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<K;++i) E4[k](i,j) = E1[k](i,j);
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
    std::vector<T> d5(nloops2x);
    std::vector<Eigen::PartialPivLU<EIGENSMA>*,ALLOC(Eigen::PartialPivLU<EIGENSMA>*) > LU5;
    if (nloops2x) { 
        A5.resize(nloops2x); 
        B5.resize(nloops2x); 
        C5.resize(nloops2x);
        D5.resize(nloops2x);
        E5.resize(nloops2x); 
        F5.resize(nloops2x); 
        LU5.resize(nloops2x,0); 
    }
    for(int k=0;k<nloops2x;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<N;++i) A5[k](i,j) = A1[k](i,j);
        for(int j=0;j<K;++j) for(int i=0;i<N;++i) C5[k](i,j) = C1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<K;++i) E5[k](i,j) = E1[k](i,j);
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
#if 1 // A -> PLU
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            LU1[k] = new tmv::LUD<tmv::Matrix<T> >(A1[k]);
#else
            //LU1[k] = A1[k];
            //LU_Decompose(LU1[k],&P1[k][0]);
            //std::cout<<"A1[k] = "<<A1[k]<<std::endl;
            A1[k].saveDiv();
            //std::cout<<"After saveDiv"<<std::endl;
            A1[k].divideUsing(tmv::LU);
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
                A0[k] = LU1[k]->getP() * LU1[k]->getL() * LU1[k]->getU();
                if (LU1[k]->isTrans()) A0[k] = A0[k].transpose().copy();
#else
                //A0[k] = LU1[k].unitLowerTri() * LU1[k].upperTri();
                //A0[k].reversePermuteRows(&P1[k][0]);
                //std::cout<<"P = "<<A1[k].lud().getP()<<std::endl;
                //std::cout<<"LU = "<<A1[k].lud().getLU()<<std::endl;
                //std::cout<<"L = "<<A1[k].lud().getL()<<std::endl;
                //std::cout<<"U = "<<A1[k].lud().getU()<<std::endl;
                A0[k] = A1[k].lud().getP() * 
                    A1[k].lud().getL() * A1[k].lud().getU();
                if (A1[k].lud().isTrans()) A0[k] = A0[k].transpose().copy();
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
            LU2[k] = new tmv::LUD<tmv::SmallMatrix<T,N,N> >(A2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"Error check: k = "<<k<<std::endl;
                //std::cout<<"A0[k] = "<<A0[0]<<std::endl;
                //std::cout<<"P = "<<LU2[k]->getP()<<std::endl;
                //std::cout<<"L = "<<LU2[k]->getL()<<std::endl;
                //std::cout<<"U = "<<LU2[k]->getU()<<std::endl;
                A0[k] = LU2[k]->getP() * LU2[k]->getL() * LU2[k]->getU();
                if (LU2[k]->isTrans()) A0[k] = A0[k].transpose().copy();
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
                A0[k] = LU3[k].unitLowerTri() * LU3[k].upperTri();
                A0[k].reversePermuteRows(&P3[k][0]);
                for(int i=0;i<N;++i) ++P3[k][i];
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
            LU4[k] = new Eigen::PartialPivLU<EIGENM>(A4[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = 
                    EIGENM(LU4[k]->matrixLU().triangularView<Eigen::UnitLower>()) *
                    EIGENM(LU4[k]->matrixLU().triangularView<Eigen::Upper>());
                temp1 = LU4[k]->permutationP().inverse() * temp1;
                e1_eigen += (A4[k]-temp1).squaredNorm();
            }
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A1[0]);
            //std::cout<<"e1_eigen = "<<e1_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            LU5[k] = new Eigen::PartialPivLU<EIGENSMA>(A5[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = 
                    EIGENM(LU5[k]->matrixLU().triangularView<Eigen::UnitLower>()) *
                    EIGENM(LU5[k]->matrixLU().triangularView<Eigen::Upper>());
                temp1 = LU5[k]->permutationP().inverse() * temp1;
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
        ClearCache();

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) D1[k] = C1[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            LU1[k]->solveInPlace(D1[k]);
#else
            D1[k] /= MPART1(A1[k]);
#endif
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
            e2_reg /= Norm(A1[0]);
            //std::cout<<"e2_reg = "<<e2_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LU2[k]->solveInPlace(D2[k]);
#else
            D2[k] /= MPART1(A2[k]);
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
            e2_small /= Norm(A1[0]);
            //std::cout<<"e2_small = "<<e2_small<<std::endl;
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
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
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
            e2_blas /= Norm(A1[0]);
            //std::cout<<"e2_blas = "<<e2_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            D4[k] = LU4[k]->solve(D4[k]);
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
            e2_eigen /= Norm(A1[0]);
            //std::cout<<"e2_eigen = "<<e2_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            D5[k] = LU5[k]->solve(D5[k]);
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
            e2_smalleigen /= Norm(A1[0]);
            //std::cout<<"e2_smalleigen = "<<e2_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // D = C / A (no decomp)
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            LU1[k]->solve(C1[k],D1[k]);
#else
            D1[k] = C1[k] / MPART1(A1[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A1[k]) * D1[k];
                e3_reg += NormSq(C0[k]-C1[k]);
            }
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= Norm(A1[0]);
            //std::cout<<"e3_reg = "<<e3_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LU2[k]->solve(C2[k],D2[k]);
#else
            D2[k] = C2[k] / MPART1(A2[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A2[k]) * D2[k];
                e3_small += NormSq(C0[k]-C2[k]);
            }
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= Norm(A1[0]);
            //std::cout<<"e3_small = "<<e3_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            D3[k] = C3[k];
            LAPNAME(getrs) (
                LAPCM LAPNT, N,K,LP(LU3[k].ptr()),N,&P3[k][0],
                LP(D3[k].ptr()),N,&lapinfo);
#elif PART1 >= 2 && PART1 <= 5
            D3[k] = C3[k];
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
                N,K,one,BP(A3[k].cptr()),N,BP(D3[k].ptr()),N
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
                C0[k] = MPART1(A3[k]) * D3[k];
                e3_blas += NormSq(C0[k]-C3[k]);
            }
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= Norm(A1[0]);
            //std::cout<<"e3_blas = "<<e3_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            D4[k] = LU4[k]->solve(C4[k]);
#elif PART1 >= 2 && PART1 <= 5
            D4[k] = EPART1(A4[k]).solveTriangular(C4[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k]) * D4[k];
                e3_eigen += (C4[k]-temp1).squaredNorm();
            }
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(A1[0]);
            //std::cout<<"e3_eigen = "<<e3_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            D5[k] = LU5[k]->solve(C5[k]);
#elif PART1 >= 2 && PART1 <= 5
            D5[k] = EPART1(A5[k]).solveTriangular(C5[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k]) * D5[k];
                e3_smalleigen += (C5[k]-temp1).squaredNorm();
            }
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(A1[0]);
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
            LU1[k]->makeInverse(B1[k]);
#else
            B1[k] = MPART1(A1[k]).inverse();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A1[k]) * B1[k];
                e4_reg += NormSq(A0[k]-1);
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
#if (PART1 == 1)
            LU2[k]->makeInverse(B2[k]);
#else
            B2[k] = MPART1(A2[k]).inverse();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A2[k]) * B2[k];
                e4_small += NormSq(A0[k]-1);
            }
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(A1[0]);
            //std::cout<<"e4_small = "<<e4_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
#if (PART1 != 1) || !defined(TISCOMPLEX) || defined(MKL)
        // Atlas LAPACK gives me a segmentation fault for zgetri
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            B3[k] = LU3[k];
            int lwork = 10*N*LAP_BLOCKSIZE;
            std::vector<RT> work(lwork);
            LAPNAME(getri) (
                LAPCM N,LP(B3[k].ptr()),N,&P3[k][0],
                LP(&work[0]),lwork,&lapinfo);
#elif PART1 >= 2 && PART1 <= 5
            B3[k].setToIdentity();
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
                N,N,one,BP(A3[k].cptr()),N,BP(B3[k].ptr()),N
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A3[k]) * B3[k];
                e4_blas += NormSq(A0[k]-1);
            }
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(A1[0]);
            //std::cout<<"e4_blas = "<<e4_blas<<std::endl;
        }
#endif
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            B4[k] = LU4[k]->solve(B4[k].Identity(N,N));
#elif PART1 >= 2 && PART1 <= 5
            B4[k] = EPART1(A4[k]).solveTriangular(B4[k].Identity(N,N));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k]) * B4[k];
                e4_eigen += (temp1-temp1.Identity(N,N)).squaredNorm();
            }
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(A1[0]);
            //std::cout<<"e4_eigen = "<<e4_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            B5[k] = LU5[k]->solve(B5[k].Identity());
#elif PART1 >= 2 && PART1 <= 5
            B5[k] = EPART1(A5[k]).solveTriangular(B5[k].Identity());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k]) * B5[k];
                e4_smalleigen += (temp1-B5[k].Identity()).squaredNorm();
            }
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(A1[0]);
            //std::cout<<"e4_smalleigen = "<<e4_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // d = A.det() (no decomp)
#ifdef ERRORCHECK
        for (int k=0; k<nloops2; ++k) {
            // I just use the TMV calculation for the error check.
            // So if there is an error there, it only turns up in comparison
            // with the other (Eigen, LAP) answers.
            d0[k] = MPART1(A1[k]).det();
        }
        //if (n==0) std::cout<<"det = "<<d0[0]<<std::endl;
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            d1[k] = LU1[k]->det();
#else
            d1[k] = MPART1(A1[k]).det();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
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
#if (PART1 == 1)
            d2[k] = LU2[k]->det();
#else
            d2[k] = MPART1(A2[k]).det();
#endif
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
#if PART1 == 1
            for(int i=0;i<N;++i) d3[k] *= LU3[k](i,i);
            for(int i=0;i<N;++i) if (P3[k][i] != i+1) d3[k] = -d3[k];
#elif PART1 >= 2 && PART1 <= 5
            for(int i=0;i<N;++i) d3[k] *= A3[k](i,i);
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
#if PART1 == 1
            d4[k] = LU4[k]->determinant();
#elif PART1 >= 2 && PART1 <= 5
            d4[k] = A4[k].determinant();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e11_eigen += tmv::TMV_NORM(d0[k]-d4[k]);
            }
            e11_eigen = sqrt(e11_eigen/nloops2);
            e11_eigen /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e11_eigen = "<<e11_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            d5[k] = LU5[k]->determinant();
#elif PART1 >= 2 && PART1 <= 5
            d5[k] = A5[k].determinant();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e11_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e11_smalleigen += tmv::TMV_NORM(d0[k]-d5[k]);
            }
            e11_smalleigen = sqrt(e11_smalleigen/nloops2);
            e11_smalleigen /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e11_smalleigen = "<<e11_smalleigen<<std::endl;
        }
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
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A1[k]) * D1[k];
                e5_reg += NormSq(C0[k]-C1[k]);
            }
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= Norm(A1[0]);
            //std::cout<<"e5_reg = "<<e5_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D2[k] /= MPART1(A2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A2[k]) * D2[k];
                e5_small += NormSq(C0[k]-C2[k]);
            }
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= Norm(A1[0]);
            //std::cout<<"e5_small = "<<e5_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        for (int k=0; k<nloops2; ++k) D3[k] = C3[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LU3[k] = A3[k];
            LAPNAME(getrf) (
                LAPCM N,N,LP(LU3[k].ptr()),N,&P3[k][0],&lapinfo);
            LAPNAME(getrs) (
                LAPCM LAPNT, N,K,LP(LU3[k].ptr()),N,&P3[k][0],
                LP(D3[k].ptr()),N,&lapinfo);
#elif PART1 >= 2 && PART1 <= 5
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
                N,K,one,BP(A3[k].cptr()),N,BP(D3[k].ptr()),N
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A3[k]) * D3[k];
                e5_blas += NormSq(C0[k]-C3[k]);
            }
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= Norm(A1[0]);
            //std::cout<<"e5_blas = "<<e5_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            D4[k] = A4[k].lu().solve(D4[k]);
#elif PART1 >= 2 && PART1 <= 5
            EPART1(A4[k]).solveTriangularInPlace(D4[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k]) * D4[k];
                e5_eigen += (C4[k]-temp1).squaredNorm();
            }
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A1[0]);
            //std::cout<<"e5_eigen = "<<e5_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            D5[k] = A5[k].inverse() * D5[k];
#elif PART1 >= 2 && PART1 <= 5
            EPART1(A5[k]).solveTriangularInPlace(D5[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k]) * D5[k];
                e5_smalleigen += (C5[k]-temp1).squaredNorm();
            }
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A1[0]);
            //std::cout<<"e5_smalleigen = "<<e5_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // D = C / A
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D1[k] = C1[k] / MPART1(A1[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A1[k]) * D1[k];
                e6_reg += NormSq(C0[k]-C1[k]);
            }
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= Norm(A1[0]);
            //std::cout<<"e6_reg = "<<e6_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D2[k] = C2[k] / MPART1(A2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A2[k]) * D2[k];
                e6_small += NormSq(C0[k]-C2[k]);
            }
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= Norm(A1[0]);
            //std::cout<<"e6_small = "<<e6_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LU3[k] = A3[k];
            LAPNAME(getrf) (
                LAPCM N,N,LP(LU3[k].ptr()),N,&P3[k][0],&lapinfo);
            D3[k] = C3[k];
            LAPNAME(getrs) (
                LAPCM LAPNT, N,K,LP(LU3[k].ptr()),N,&P3[k][0],
                LP(D3[k].ptr()),N,&lapinfo);
#elif PART1 >= 2 && PART1 <= 5
            D3[k] = C3[k];
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
                N,K,one,BP(A3[k].cptr()),N,BP(D3[k].ptr()),N
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A3[k]) * D3[k];
                e6_blas += NormSq(C0[k]-C3[k]);
            }
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= Norm(A1[0]);
            //std::cout<<"e6_blas = "<<e6_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            D4[k] = A4[k].lu().solve(C4[k]);
#elif PART1 >= 2 && PART1 <= 5
            D4[k] = EPART1(A4[k]).solveTriangular(C4[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k]) * D4[k];
                e6_eigen += (C4[k]-temp1).squaredNorm();
            }
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(A1[0]);
            //std::cout<<"e6_eigen = "<<e6_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            D5[k] = A5[k].inverse() * C5[k];
#elif PART1 >= 2 && PART1 <= 5
            D5[k] = EPART1(A5[k]).solveTriangular(C5[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k]) * D5[k];
                e6_smalleigen += (C5[k]-temp1).squaredNorm();
            }
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(A1[0]);
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
            B1[k] = MPART1(A1[k]).inverse();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A1[k]) * B1[k];
                e7_reg += NormSq(A0[k]-1);
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
            B2[k] = MPART1(A2[k]).inverse();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A2[k]) * B2[k];
                e7_small += NormSq(A0[k]-1);
            }
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= Norm(A1[0]);
            //std::cout<<"e7_small = "<<e7_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
#if (PART1 != 1) || !defined(TISCOMPLEX) || defined(MKL)
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            B3[k] = A3[k];
            LAPNAME(getrf) (
                LAPCM N,N,LP(B3[k].ptr()),N,&P3[k][0],&lapinfo);
            int lwork = N*LAP_BLOCKSIZE;
            std::vector<RT> work(lwork);
            LAPNAME(getri) (
                LAPCM N,LP(B3[k].ptr()),N,&P3[k][0],
                LP(&work[0]),lwork,&lapinfo);
#elif PART1 >= 2 && PART1 <= 5
            B3[k].setToIdentity();
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
                N,N,one,BP(A3[k].cptr()),N,BP(B3[k].ptr()),N
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A3[k]) * B3[k];
                e7_blas += NormSq(A0[k]-1);
            }
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= Norm(A1[0]);
            //std::cout<<"e7_blas = "<<e7_blas<<std::endl;
        }
#endif
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            B4[k] = A4[k].inverse();
#elif PART1 >= 2 && PART1 <= 5
            B4[k] = EPART1(A4[k]).solveTriangular(B4[k].Identity(N,N));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k]) * B4[k];
                e7_eigen += (temp1-temp1.Identity(N,N)).squaredNorm();
            }
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= Norm(A1[0]);
            //std::cout<<"e7_eigen = "<<e7_eigen<<std::endl;
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
            B5[k] = A5[k].inverse();
#elif PART1 >= 2 && PART1 <= 5
            B5[k] = EPART1(A5[k]).solveTriangular(B5[k].Identity());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k]) * B5[k];
                e7_smalleigen += (temp1-B5[k].Identity()).squaredNorm();
            }
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= Norm(A1[0]);
            //std::cout<<"e7_smalleigen = "<<e7_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // B = B.inverse()
        ClearCache();

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) { B1[k] = A1[k]; }
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 >= 2 && PART1 <= 5
            MPART1(B1[k]).invertSelf();
#else
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
                A0[k] = MPART1(A1[k]) * MPART1(B1[k]);
                e8_reg += NormSq(A0[k]-1);
            }
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= Norm(A1[0]);
            //std::cout<<"e8_reg = "<<e8_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) { B2[k] = A2[k]; }
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 >= 2 && PART1 <= 5
            MPART1(B2[k]).invertSelf();
#else
            B2[k] = MPART1(B2[k]).inverse();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A2[k]) * MPART1(B2[k]);
                e8_small += NormSq(A0[k]-1);
            }
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= Norm(A1[0]);
            //std::cout<<"e8_small = "<<e8_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
#if (PART1 != 1) || !defined(TISCOMPLEX) || defined(MKL)
        for (int k=0; k<nloops2; ++k) { B3[k] = A3[k]; }
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LAPNAME(getrf) (
                LAPCM N,N,LP(B3[k].ptr()),N,&P3[k][0],&lapinfo);
            int lwork = N*LAP_BLOCKSIZE;
            std::vector<RT> work(lwork);
            LAPNAME(getri) (
                LAPCM N,LP(B3[k].ptr()),N,&P3[k][0],
                LP(&work[0]),lwork,&lapinfo);
#elif PART1 >= 2 && PART1 <= 5
            B3[k].setToIdentity();
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
                N,N,one,BP(A3[k].cptr()),N,BP(B3[k].ptr()),N
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A3[k]) * MPART1(B3[k]);
                e8_blas += NormSq(A0[k]-1);
            }
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= Norm(A1[0]);
            //std::cout<<"e8_blas = "<<e8_blas<<std::endl;
        }
#endif
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) { B4[k] = A4[k]; }
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            B4[k] = B4[k].inverse();
#elif PART1 >= 2 && PART1 <= 5
            B4[k] = EPART1(A4[k]).solveTriangular(B4[k].Identity(N,N));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k]) * EPART1(B4[k]);
                e8_eigen += (temp1-temp1.Identity(N,N)).squaredNorm();
            }
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= Norm(A1[0]);
            //std::cout<<"e8_eigen = "<<e8_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2; ++k) { B5[k] = A5[k]; }
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            B5[k] = B5[k].inverse().eval();
#elif PART1 >= 2 && PART1 <= 5
            B5[k] = EPART1(A5[k]).solveTriangular(B5[k].Identity());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k]) * EPART1(B5[k]);
                e8_smalleigen += (temp1-B5[k].Identity()).squaredNorm();
            }
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= Norm(A1[0]);
            //std::cout<<"e8_smalleigen = "<<e8_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // F %= A
        ClearCache();

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) { F1[k] = E1[k]; }
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F1[k] %= MPART1(A1[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F1[k] * MPART1(A1[k]);
                e9_reg += NormSq(E0[k]-E1[k]);
            }
            e9_reg = sqrt(e9_reg/nloops2);
            e9_reg /= Norm(A1[0]);
            //std::cout<<"e9_reg = "<<e9_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) F2[k] = E2[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F2[k] %= MPART1(A2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F2[k] * MPART1(A2[k]);
                e9_small += NormSq(E0[k]-E2[k]);
            }
            e9_small = sqrt(e9_small/nloops2);
            e9_small /= Norm(A1[0]);
            //std::cout<<"e9_small = "<<e9_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        for (int k=0; k<nloops2; ++k) F3[k] = E3[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LU3[k] = A3[k];
            LAPNAME(getrf) (
                LAPCM N,N,LP(LU3[k].ptr()),N,&P3[k][0],&lapinfo);
            // This only works for RowMajor F3, so need a temporary.
            tmv::Matrix<T,tmv::RowMajor> temp = F3[k];
            LAPNAME(getrs) (
                LAPCM LAPT,N,K,LP(LU3[k].ptr()),N,&P3[k][0],
                LP(temp.ptr()),N,&lapinfo);
            F3[k] = temp;
#elif PART1 >= 2 && PART1 <= 5
            BLASNAME(trsm) (
                BLASCM BLASRight, BLASUPLO1, BLASNT, BLASDIAG1,
                K,N,one,BP(A3[k].cptr()),N,BP(F3[k].ptr()),K
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F3[k] * MPART1(A3[k]);
                e9_blas += NormSq(E0[k]-E3[k]);
            }
            e9_blas = sqrt(e9_blas/nloops2);
            e9_blas /= Norm(A1[0]);
            //std::cout<<"e9_blas = "<<e9_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) F4[k] = E4[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            F4[k].transpose() = A4[k].transpose().lu().solve(F4[k].transpose());
#elif PART1 >= 2 && PART1 <= 5
            F4[k].transpose() = 
                EPART1(A4[k]).transpose().solveTriangular(F4[k].transpose());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = F4[k] * EPART1(A4[k]);
                e9_eigen += (E4[k]-temp1).squaredNorm();
            }
            e9_eigen = sqrt(e9_eigen/nloops2);
            e9_eigen /= Norm(A1[0]);
            //std::cout<<"e9_eigen = "<<e9_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) F5[k] = E5[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            F5[k] *= A5[k].inverse();
#elif PART1 >= 2 && PART1 <= 5
            EPART1(A5[k]).transpose().solveTriangularInPlace(F5[k].transpose());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e9_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = F5[k] * EPART1(A5[k]);
                e9_smalleigen += (E5[k]-temp1).squaredNorm();
            }
            e9_smalleigen = sqrt(e9_smalleigen/nloops2);
            e9_smalleigen /= Norm(A1[0]);
            //std::cout<<"e9_smalleigen = "<<e9_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // F = E % A
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F1[k] = E1[k] % MPART1(A1[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F1[k] * MPART1(A1[k]);
                e10_reg += NormSq(E0[k]-E1[k]);
            }
            e10_reg = sqrt(e10_reg/nloops2);
            e10_reg /= Norm(A1[0]);
            //std::cout<<"e10_reg = "<<e10_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F2[k] = E2[k] % MPART1(A2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F2[k] * MPART1(A2[k]);
                e10_small += NormSq(E0[k]-E2[k]);
            }
            e10_small = sqrt(e10_small/nloops2);
            e10_small /= Norm(A1[0]);
            //std::cout<<"e10_small = "<<e10_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LU3[k] = A3[k];
            LAPNAME(getrf) (
                LAPCM N,N,LP(LU3[k].ptr()),N,&P3[k][0],&lapinfo);
            tmv::Matrix<T,tmv::RowMajor> temp = E3[k];
            LAPNAME(getrs) (
                LAPCM LAPT,N,K,LP(LU3[k].ptr()),N,&P3[k][0],
                LP(temp.ptr()),N,&lapinfo);
            F3[k] = temp;
#elif PART1 >= 2 && PART1 <= 5
            F3[k] = E3[k];
            BLASNAME(trsm) (
                BLASCM BLASRight, BLASUPLO1, BLASNT, BLASDIAG1,
                K,N,one,BP(A3[k].cptr()),N,BP(F3[k].ptr()),K
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F3[k] * MPART1(A3[k]);
                e10_blas += NormSq(E0[k]-E3[k]);
            }
            e10_blas = sqrt(e10_blas/nloops2);
            e10_blas /= Norm(A1[0]);
            //std::cout<<"e10_blas = "<<e10_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            F4[k].transpose() = A4[k].transpose().lu().solve(E4[k].transpose());
#elif PART1 >= 2 && PART1 <= 5
            F4[k].transpose() =
                EPART1(A4[k]).transpose().solveTriangular(E4[k].transpose());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = F4[k] * EPART1(A4[k]);
                e10_eigen += (E4[k]-temp1).squaredNorm();
            }
            e10_eigen = sqrt(e10_eigen/nloops2);
            e10_eigen /= Norm(A1[0]);
            //std::cout<<"e10_eigen = "<<e10_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            F5[k] = E5[k] * A5[k].inverse();
#elif PART1 >= 2 && PART1 <= 5
            F5[k].transpose() = 
                EPART1(A5[k]).transpose().solveTriangular(E5[k].transpose());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e10_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = F5[k] * EPART1(A5[k]);
                e10_smalleigen += (E5[k]-temp1).squaredNorm();
            }
            e10_smalleigen = sqrt(e10_smalleigen/nloops2);
            e10_smalleigen /= Norm(A1[0]);
            //std::cout<<"e10_smalleigen = "<<e10_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // d = A.det()
#ifdef ERRORCHECK
        for (int k=0; k<nloops2; ++k) {
            // I just use the TMV calculation for the error check.
            // So if there is an error there, it only turns up in comparison
            // with the other 4 answers.
            d0[k] = MPART1(A1[k]).det();
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            d1[k] = MPART1(A1[k]).det();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
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
            d2[k] = MPART1(A2[k]).det();
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
            d3[k] = T(1);
#if PART1 == 1
            LU3[k] = A3[k];
            LAPNAME(getrf) (
                LAPCM N,N,LP(LU3[k].ptr()),N,&P3[k][0],&lapinfo);
            for(int i=0;i<N;++i) d3[k] *= LU3[k](i,i);
            for(int i=0;i<N;++i) if (P3[k][i] != i+1) d3[k] = -d3[k];
#elif PART1 >= 2 && PART1 <= 5
            for(int i=0;i<N;++i) d3[k] *= A3[k](i,i);
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
            d4[k] = A4[k].determinant();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e12_eigen += tmv::TMV_NORM(d0[k]-d4[k]);
            }
            e12_eigen = sqrt(e12_eigen/nloops2);
            e12_eigen /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e12_eigen = "<<e12_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            d5[k] = A5[k].determinant();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e12_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e12_smalleigen += tmv::TMV_NORM(d0[k]-d5[k]);
            }
            e12_smalleigen = sqrt(e12_smalleigen/nloops2);
            e12_smalleigen /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e12_smalleigen = "<<e12_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#ifdef SHOW_DOTS
        std::cout<<"."; std::cout.flush();
#endif

#ifdef DOREG
#ifdef TMV_NO_LIB
        for(int i=0;i<nloops2;++i) if (LU1[i]) delete LU1[i];
#endif
#endif
#ifdef DOSMALL
        for(int i=0;i<nloops2;++i) if (LU2[i]) delete LU2[i];
#endif
#ifdef DOEIGEN
        for(int i=0;i<nloops2;++i) if (LU4[i]) delete LU4[i];
#endif
#ifdef DOEIGENSMALL
        for(int i=0;i<nloops2x;++i) if (LU5[i]) delete LU5[i];
#endif
    }
    std::cout<<"\n";

#ifdef DO_SEPARATE_DECOMP
    std::cout<<"A -> PLU               "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"No decomposition:\n";
    std::cout<<"D /= A                 "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"D = C / A              "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"B = A.inverse()        "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"d = A.det()            "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
    std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
    std::cout<<"Including decomposition:\n";
#endif
    std::cout<<"D /= A                 "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
    std::cout<<"D = C / A              "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"B = A.inverse()        "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"B = B.inverse()        "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
    std::cout<<"F %= A                 "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
    std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
    std::cout<<"F = E % A              "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
    std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
    std::cout<<"d = A.det()            "<<t12_reg<<"  "<<t12_small<<"  "<<t12_blas;
    std::cout<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
#ifdef DO_SEPARATE_DECOMP
    std::cout<<"A -> PLU               "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"D /= A ND              "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"D = C / A ND           "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"B = A.inverse() ND     "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"d = A.det() ND         "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
    std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
#endif
    std::cout<<"D /= A                 "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
    std::cout<<"D = C / A              "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
    std::cout<<"B = A.inverse()        "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
    std::cout<<"B = B.inverse()        "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
    std::cout<<"F %= A                 "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
    std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
    std::cout<<"F = E % A              "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
    std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
    std::cout<<"d = A.det()            "<<e12_reg<<"  "<<e12_small<<"  "<<e12_blas;
    std::cout<<"  "<<e12_eigen<<"  "<<e12_smalleigen<<std::endl;
    std::cout<<"\n\n";
#endif
}
#endif

#ifdef DO_R
static void LU_R(
    const std::vector<tmv::Matrix<T,tmv::RowMajor> >& A1,
    std::vector<tmv::Matrix<T> >& B1,
    const std::vector<tmv::Matrix<T> >& C1,
    std::vector<tmv::Matrix<T> >& D1,
    const std::vector<tmv::Matrix<T> >& E1,
    std::vector<tmv::Matrix<T> >& F1)
{
    //std::cout<<"Start LU_R.\n";
    //std::cout<<"A1.size = "<<A1.size()<<std::endl;
#ifdef TMV_NO_LIB
    std::vector<tmv::LUD<tmv::Matrix<T,tmv::RowMajor> >*> LU1(nloops2,0);
#endif
    std::vector<T> d1(nloops2);

#ifdef ERRORCHECK
    std::vector<tmv::Matrix<T,tmv::RowMajor> > A0 = A1;
    std::vector<tmv::Matrix<T> > C0 = C1;
    std::vector<tmv::Matrix<T> > E0 = E1;
    std::vector<T> d0(nloops2);
    //std::cout<<"Made A0, etc."<<std::endl;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallMatrix<T,N,N,tmv::RowMajor> > A2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,N> > B2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,K> > C2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,K> > D2(nloops2);
    std::vector<tmv::SmallMatrix<T,K,N> > E2(nloops2);
    std::vector<tmv::SmallMatrix<T,K,N> > F2(nloops2);
    std::vector<T> d2(nloops2);
    std::vector<tmv::LUD<tmv::SmallMatrix<T,N,N,tmv::RowMajor> >*> LU2(nloops2,0);
    //std::cout<<"Made A2, etc."<<std::endl;

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; C2[k] = C1[k]; E2[k] = E1[k];
    }
    //std::cout<<"Filled A2,C2,E2"<<std::endl;
#endif

#ifdef DOLAP
    std::vector<tmv::Matrix<T,tmv::RowMajor> > A3 = A1;
    std::vector<tmv::Matrix<T> > B3 = B1;
    std::vector<tmv::Matrix<T> > C3 = C1;
    std::vector<tmv::Matrix<T> > D3 = D1;
    std::vector<tmv::Matrix<T> > E3 = E1;
    std::vector<tmv::Matrix<T> > F3 = F1;
    std::vector<T> d3(nloops2);
    std::vector<tmv::Matrix<T> > LU3(nloops2);
    std::vector<std::vector<int> > P3(nloops2);
    for(int k=0;k<nloops2;++k) {
        LU3[k].resize(A1[k].colsize(),A1[k].rowsize());
        LU3[k] = A1[k];
        P3[k].resize(N);
    }
    //std::cout<<"Made A3, etc."<<std::endl;
#endif

#ifdef DOEIGEN
    std::vector<EIGENM,ALLOC(EIGENM) > A4(nloops2,EIGENM(N,N));
    std::vector<EIGENM,ALLOC(EIGENM) > B4(nloops2,EIGENM(N,N));
    std::vector<EIGENM,ALLOC(EIGENM) > C4(nloops2,EIGENM(N,K));
    std::vector<EIGENM,ALLOC(EIGENM) > D4(nloops2,EIGENM(N,K));
    std::vector<EIGENM,ALLOC(EIGENM) > E4(nloops2,EIGENM(K,N));
    std::vector<EIGENM,ALLOC(EIGENM) > F4(nloops2,EIGENM(K,N));
    std::vector<T> d4(nloops2);
    std::vector<Eigen::PartialPivLU<EIGENM>*,ALLOC(Eigen::PartialPivLU<EIGENM>*) > LU4(nloops2);
    for(int k=0;k<nloops2;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<N;++i) A4[k](i,j) = A1[k](i,j);
        for(int j=0;j<K;++j) for(int i=0;i<N;++i) C4[k](i,j) = C1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<K;++i) E4[k](i,j) = E1[k](i,j);
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
    std::vector<T> d5(nloops2x);
    std::vector<Eigen::PartialPivLU<EIGENSMA>*,ALLOC(Eigen::PartialPivLU<EIGENSMA>*) > LU5;
    if (nloops2x) { 
        A5.resize(nloops2x); 
        B5.resize(nloops2x); 
        C5.resize(nloops2x);
        D5.resize(nloops2x);
        E5.resize(nloops2x); 
        F5.resize(nloops2x); 
        LU5.resize(nloops2x,0); 
    }
    for(int k=0;k<nloops2x;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<N;++i) A5[k](i,j) = A1[k](i,j);
        for(int j=0;j<K;++j) for(int i=0;i<N;++i) C5[k](i,j) = C1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<K;++i) E5[k](i,j) = E1[k](i,j);
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
#if 1 // A -> PLU
        ClearCache();
        //std::cout<<"After clear cache"<<std::endl;

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            LU1[k] = new tmv::LUD<tmv::Matrix<T,tmv::RowMajor> >(A1[k]);
#else
            //LU1[k] = A1[k];
            //LU_Decompose(LU1[k],&P1[k][0]);
            //std::cout<<"A1[k] = "<<A1[k]<<std::endl;
            A1[k].saveDiv();
            //std::cout<<"After saveDiv"<<std::endl;
            A1[k].divideUsing(tmv::LU);
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
                A0[k] = LU1[k]->getP() * LU1[k]->getL() * LU1[k]->getU();
                if (LU1[k]->isTrans()) A0[k] = A0[k].transpose().copy();
#else
                //A0[k] = LU1[k].unitLowerTri() * LU1[k].upperTri();
                //A0[k].reversePermuteRows(&P1[k][0]);
                //std::cout<<"P = "<<A1[k].lud().getP()<<std::endl;
                //std::cout<<"LU = "<<A1[k].lud().getLU()<<std::endl;
                //std::cout<<"L = "<<A1[k].lud().getL()<<std::endl;
                //std::cout<<"U = "<<A1[k].lud().getU()<<std::endl;
                A0[k] = A1[k].lud().getP() * 
                    A1[k].lud().getL() * A1[k].lud().getU();
                if (A1[k].lud().isTrans()) A0[k] = A0[k].transpose().copy();
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
            LU2[k] = new tmv::LUD<tmv::SmallMatrix<T,N,N,tmv::RowMajor> >(A2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                //std::cout<<"Error check: k = "<<k<<std::endl;
                //std::cout<<"A2[k] = "<<A2[0]<<std::endl;
                //std::cout<<"P = "<<LU2[k]->getP()<<std::endl;
                //std::cout<<"L = "<<LU2[k]->getL()<<std::endl;
                //std::cout<<"U = "<<LU2[k]->getU()<<std::endl;
                //std::cout<<"P*L = "<<LU2[k]->getP() * LU2[k]->getL()<<std::endl;
                //std::cout<<"(P*L)*U = "<<(LU2[k]->getP() * LU2[k]->getL())*LU2[k]->getU()<<std::endl;
                //std::cout<<"L*U = "<<LU2[k]->getL() * LU2[k]->getU()<<std::endl;
                //std::cout<<"P*(L*U) = "<<LU2[k]->getP() * (LU2[k]->getL()*LU2[k]->getU()).calc()<<std::endl;
                A0[k] = LU2[k]->getP() * LU2[k]->getL() * LU2[k]->getU();
                //std::cout<<"A0[k] = "<<A0[0]<<std::endl;
                if (LU2[k]->isTrans()) A0[k] = A0[k].transpose().copy();
                //std::cout<<"A0[k] = "<<A0[0]<<std::endl;
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
                A0[k] = LU3[k].unitLowerTri() * LU3[k].upperTri();
                A0[k].reversePermuteRows(&P3[k][0]);
                for(int i=0;i<N;++i) ++P3[k][i];
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
            LU4[k] = new Eigen::PartialPivLU<EIGENM>(A4[k].transpose());
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = 
                    EIGENM(LU4[k]->matrixLU().triangularView<Eigen::UnitLower>()) *
                    EIGENM(LU4[k]->matrixLU().triangularView<Eigen::Upper>());
                temp1 = LU4[k]->permutationP().inverse() * temp1;
                EIGENM temp2 = temp1.transpose();
                e1_eigen += (A4[k]-temp2).squaredNorm();
            }
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A1[0]);
            //std::cout<<"e1_eigen = "<<e1_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            LU5[k] = new Eigen::PartialPivLU<EIGENSMA>(A5[k].transpose());
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = 
                    EIGENM(LU5[k]->matrixLU().triangularView<Eigen::UnitLower>()) *
                    EIGENM(LU5[k]->matrixLU().triangularView<Eigen::Upper>());
                temp1 = LU5[k]->permutationP().inverse() * temp1;
                EIGENM temp2 = temp1.transpose();
                e1_smalleigen += (A5[k]-temp2).squaredNorm();
            }
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A1[0]);
            //std::cout<<"e1_smalleigen = "<<e1_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // D /= A (no decomp)
        ClearCache();

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) D1[k] = C1[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            LU1[k]->solveInPlace(D1[k]);
#else
            D1[k] /= MPART1(A1[k]);
#endif
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
            e2_reg /= Norm(A1[0]);
            //std::cout<<"e2_reg = "<<e2_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LU2[k]->solveInPlace(D2[k]);
#else
            D2[k] /= MPART1(A2[k]);
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
            e2_small /= Norm(A1[0]);
            //std::cout<<"e2_small = "<<e2_small<<std::endl;
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
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
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
            e2_blas /= Norm(A1[0]);
            //std::cout<<"e2_blas = "<<e2_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            D4[k] = LU4[k]->solve(D4[k]);
#elif PART1 >= 2 && PART1 <= 5
            EPART1(A4[k].transpose()).solveTriangularInPlace(D4[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k].transpose()) * D4[k];
                e2_eigen += (C4[k]-temp1).squaredNorm();
            }
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm(A1[0]);
            //std::cout<<"e2_eigen = "<<e2_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            D5[k] = LU5[k]->solve(D5[k]);
#elif PART1 >= 2 && PART1 <= 5
            EPART1(A5[k].transpose()).solveTriangularInPlace(D5[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k].transpose()) * D5[k];
                e2_smalleigen += (C5[k]-temp1).squaredNorm();
            }
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm(A1[0]);
            //std::cout<<"e2_smalleigen = "<<e2_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // D = C / A (no decomp)
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            LU1[k]->solve(C1[k],D1[k]);
#else
            D1[k] = C1[k] / MPART1(A1[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A1[k]) * D1[k];
                e3_reg += NormSq(C0[k]-C1[k]);
            }
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= Norm(A1[0]);
            //std::cout<<"e3_reg = "<<e3_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LU2[k]->solve(C2[k],D2[k]);
#else
            D2[k] = C2[k] / MPART1(A2[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A2[k]) * D2[k];
                e3_small += NormSq(C0[k]-C2[k]);
            }
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= Norm(A1[0]);
            //std::cout<<"e3_small = "<<e3_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            D3[k] = C3[k];
            LAPNAME(getrs) (
                LAPCM LAPNT, N,K,LP(LU3[k].ptr()),N,&P3[k][0],
                LP(D3[k].ptr()),N,&lapinfo);
#elif PART1 >= 2 && PART1 <= 5
            D3[k] = C3[k];
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
                N,K,one,BP(A3[k].cptr()),N,BP(D3[k].ptr()),N
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
                C0[k] = MPART1(A3[k]) * D3[k];
                e3_blas += NormSq(C0[k]-C3[k]);
            }
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= Norm(A1[0]);
            //std::cout<<"e3_blas = "<<e3_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            D4[k] = LU4[k]->solve(C4[k]);
#elif PART1 >= 2 && PART1 <= 5
            D4[k] = EPART1(A4[k].transpose()).solveTriangular(C4[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k].transpose()) * D4[k];
                e3_eigen += (C4[k]-temp1).squaredNorm();
            }
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(A1[0]);
            //std::cout<<"e3_eigen = "<<e3_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            D5[k] = LU5[k]->solve(C5[k]);
#elif PART1 >= 2 && PART1 <= 5
            D5[k] = EPART1(A5[k].transpose()).solveTriangular(C5[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k].transpose()) * D5[k];
                e3_smalleigen += (C5[k]-temp1).squaredNorm();
            }
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(A1[0]);
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
            LU1[k]->makeInverse(B1[k]);
#else
            B1[k] = MPART1(A1[k]).inverse();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A1[k]) * B1[k];
                e4_reg += NormSq(A0[k]-1);
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
#if (PART1 == 1)
            LU2[k]->makeInverse(B2[k]);
#else
            B2[k] = MPART1(A2[k]).inverse();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A2[k]) * B2[k];
                e4_small += NormSq(A0[k]-1);
            }
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(A1[0]);
            //std::cout<<"e4_small = "<<e4_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
#if (PART1 != 1) || !defined(TISCOMPLEX) || defined(MKL)
        // Atlas LAPACK gives me a segmentation fault for zgetri
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            B3[k] = LU3[k];
            int lwork = 10*N*LAP_BLOCKSIZE;
            std::vector<RT> work(lwork);
            LAPNAME(getri) (
                LAPCM N,LP(B3[k].ptr()),N,&P3[k][0],
                LP(&work[0]),lwork,&lapinfo);
#elif PART1 >= 2 && PART1 <= 5
            B3[k].setToIdentity();
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
                N,N,one,BP(A3[k].cptr()),N,BP(B3[k].ptr()),N
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A3[k]) * B3[k];
                e4_blas += NormSq(A0[k]-1);
            }
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(A1[0]);
            //std::cout<<"e4_blas = "<<e4_blas<<std::endl;
        }
#endif
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            B4[k] = LU4[k]->solve(B4[k].Identity(N,N));
#elif PART1 >= 2 && PART1 <= 5
            B4[k] = EPART1(A4[k].transpose()).solveTriangular(B4[k].Identity(N,N));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k].transpose()) * B4[k];
                e4_eigen += (temp1-temp1.Identity(N,N)).squaredNorm();
            }
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(A1[0]);
            //std::cout<<"e4_eigen = "<<e4_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            B5[k] = LU5[k]->solve(B5[k].Identity());
#elif PART1 >= 2 && PART1 <= 5
            B5[k] = EPART1(A5[k].transpose()).solveTriangular(B5[k].Identity());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k].transpose()) * B5[k];
                e4_smalleigen += (temp1-B5[k].Identity()).squaredNorm();
            }
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(A1[0]);
            //std::cout<<"e4_smalleigen = "<<e4_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // d = A.det() (no decomp)
#ifdef ERRORCHECK
        for (int k=0; k<nloops2; ++k) {
            // I just use the TMV calculation for the error check.
            // So if there is an error there, it only turns up in comparison
            // with the other (Eigen, LAP) answers.
            d0[k] = MPART1(A1[k]).det();
        }
        //if (n==0) std::cout<<"det = "<<d0[0]<<std::endl;
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TMV_NO_LIB
            d1[k] = LU1[k]->det();
#else
            d1[k] = MPART1(A1[k]).det();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
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
#if (PART1 == 1)
            d2[k] = LU2[k]->det();
#else
            d2[k] = MPART1(A2[k]).det();
#endif
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
#if PART1 == 1
            for(int i=0;i<N;++i) d3[k] *= LU3[k](i,i);
            for(int i=0;i<N;++i) if (P3[k][i] != i+1) d3[k] = -d3[k];
#elif PART1 >= 2 && PART1 <= 5
            for(int i=0;i<N;++i) d3[k] *= A3[k](i,i);
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
#if PART1 == 1
            d4[k] = LU4[k]->determinant();
#elif PART1 >= 2 && PART1 <= 5
            d4[k] = A4[k].transpose().determinant();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e11_eigen += tmv::TMV_NORM(d0[k]-d4[k]);
            }
            e11_eigen = sqrt(e11_eigen/nloops2);
            e11_eigen /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e11_eigen = "<<e11_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            d5[k] = LU5[k]->determinant();
#elif PART1 >= 2 && PART1 <= 5
            d5[k] = A5[k].transpose().determinant();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e11_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e11_smalleigen += tmv::TMV_NORM(d0[k]-d5[k]);
            }
            e11_smalleigen = sqrt(e11_smalleigen/nloops2);
            e11_smalleigen /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e11_smalleigen = "<<e11_smalleigen<<std::endl;
        }
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
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A1[k]) * D1[k];
                e5_reg += NormSq(C0[k]-C1[k]);
            }
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= Norm(A1[0]);
            //std::cout<<"e5_reg = "<<e5_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D2[k] /= MPART1(A2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A2[k]) * D2[k];
                e5_small += NormSq(C0[k]-C2[k]);
            }
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= Norm(A1[0]);
            //std::cout<<"e5_small = "<<e5_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        for (int k=0; k<nloops2; ++k) D3[k] = C3[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LU3[k] = A3[k];
            LAPNAME(getrf) (
                LAPCM N,N,LP(LU3[k].ptr()),N,&P3[k][0],&lapinfo);
            LAPNAME(getrs) (
                LAPCM LAPNT, N,K,LP(LU3[k].ptr()),N,&P3[k][0],
                LP(D3[k].ptr()),N,&lapinfo);
#elif PART1 >= 2 && PART1 <= 5
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
                N,K,one,BP(A3[k].cptr()),N,BP(D3[k].ptr()),N
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A3[k]) * D3[k];
                e5_blas += NormSq(C0[k]-C3[k]);
            }
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= Norm(A1[0]);
            //std::cout<<"e5_blas = "<<e5_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            D4[k] = A4[k].transpose().lu().solve(D4[k]);
#elif PART1 >= 2 && PART1 <= 5
            EPART1(A4[k].transpose()).solveTriangularInPlace(D4[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k].transpose()) * D4[k];
                e5_eigen += (C4[k]-temp1).squaredNorm();
            }
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A1[0]);
            //std::cout<<"e5_eigen = "<<e5_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            D5[k] = A5[k].transpose().inverse() * D5[k];
#elif PART1 >= 2 && PART1 <= 5
            EPART1(A5[k].transpose()).solveTriangularInPlace(D5[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k].transpose()) * D5[k];
                e5_smalleigen += (C5[k]-temp1).squaredNorm();
            }
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A1[0]);
            //std::cout<<"e5_smalleigen = "<<e5_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // D = C / A
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D1[k] = C1[k] / MPART1(A1[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A1[k]) * D1[k];
                e6_reg += NormSq(C0[k]-C1[k]);
            }
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= Norm(A1[0]);
            //std::cout<<"e6_reg = "<<e6_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            D2[k] = C2[k] / MPART1(A2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A2[k]) * D2[k];
                e6_small += NormSq(C0[k]-C2[k]);
            }
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= Norm(A1[0]);
            //std::cout<<"e6_small = "<<e6_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LU3[k] = A3[k];
            LAPNAME(getrf) (
                LAPCM N,N,LP(LU3[k].ptr()),N,&P3[k][0],&lapinfo);
            D3[k] = C3[k];
            LAPNAME(getrs) (
                LAPCM LAPNT, N,K,LP(LU3[k].ptr()),N,&P3[k][0],
                LP(D3[k].ptr()),N,&lapinfo);
#elif PART1 >= 2 && PART1 <= 5
            D3[k] = C3[k];
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
                N,K,one,BP(A3[k].cptr()),N,BP(D3[k].ptr()),N
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                C0[k] = MPART1(A3[k]) * D3[k];
                e6_blas += NormSq(C0[k]-C3[k]);
            }
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= Norm(A1[0]);
            //std::cout<<"e6_blas = "<<e6_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            D4[k] = A4[k].transpose().lu().solve(C4[k]);
#elif PART1 >= 2 && PART1 <= 5
            D4[k] = EPART1(A4[k].transpose()).solveTriangular(C4[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k].transpose()) * D4[k];
                e6_eigen += (C4[k]-temp1).squaredNorm();
            }
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(A1[0]);
            //std::cout<<"e6_eigen = "<<e6_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            D5[k] = A5[k].transpose().inverse() * C5[k];
#elif PART1 >= 2 && PART1 <= 5
            D5[k] = EPART1(A5[k].transpose()).solveTriangular(C5[k]);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k].transpose()) * D5[k];
                e6_smalleigen += (C5[k]-temp1).squaredNorm();
            }
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(A1[0]);
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
            B1[k] = MPART1(A1[k]).inverse();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A1[k]) * B1[k];
                e7_reg += NormSq(A0[k]-1);
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
            B2[k] = MPART1(A2[k]).inverse();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A2[k]) * B2[k];
                e7_small += NormSq(A0[k]-1);
            }
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= Norm(A1[0]);
            //std::cout<<"e7_small = "<<e7_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
#if (PART1 != 1) || !defined(TISCOMPLEX) || defined(MKL)
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            B3[k] = A3[k];
            LAPNAME(getrf) (
                LAPCM N,N,LP(B3[k].ptr()),N,&P3[k][0],&lapinfo);
            int lwork = N*LAP_BLOCKSIZE;
            std::vector<RT> work(lwork);
            LAPNAME(getri) (
                LAPCM N,LP(B3[k].ptr()),N,&P3[k][0],
                LP(&work[0]),lwork,&lapinfo);
#elif PART1 >= 2 && PART1 <= 5
            B3[k].setToIdentity();
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
                N,N,one,BP(A3[k].cptr()),N,BP(B3[k].ptr()),N
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A3[k]) * B3[k];
                e7_blas += NormSq(A0[k]-1);
            }
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= Norm(A1[0]);
            //std::cout<<"e7_blas = "<<e7_blas<<std::endl;
        }
#endif
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            B4[k] = A4[k].transpose().inverse();
#elif PART1 >= 2 && PART1 <= 5
            B4[k] = EPART1(A4[k].transpose()).solveTriangular(B4[k].Identity(N,N));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k].transpose()) * B4[k];
                e7_eigen += (temp1-temp1.Identity(N,N)).squaredNorm();
            }
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= Norm(A1[0]);
            //std::cout<<"e7_eigen = "<<e7_eigen<<std::endl;
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
            B5[k] = A5[k].transpose().inverse();
#elif PART1 >= 2 && PART1 <= 5
            B5[k] = EPART1(A5[k].transpose()).solveTriangular(B5[k].Identity());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k].transpose()) * B5[k];
                e7_smalleigen += (temp1-B5[k].Identity()).squaredNorm();
            }
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= Norm(A1[0]);
            //std::cout<<"e7_smalleigen = "<<e7_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // B = B.inverse()
        ClearCache();

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) { B1[k] = A1[k].transpose(); }
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 >= 2 && PART1 <= 5
            MPART1(B1[k].transpose()).invertSelf();
#else
            B1[k].transpose() = MPART1(B1[k].transpose()).inverse();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A1[k].transpose()) * MPART1(B1[k]);
                e8_reg += NormSq(A0[k]-1);
            }
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= Norm(A1[0]);
            //std::cout<<"e8_reg = "<<e8_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) { B2[k] = A2[k]; }
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 >= 2 && PART1 <= 5
            MPART1(B2[k].transpose()).invertSelf();
#else
            B2[k].transpose() = MPART1(B2[k].transpose()).inverse();
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A2[k]) * MPART1(B2[k]);
                e8_small += NormSq(A0[k]-1);
            }
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= Norm(A1[0]);
            //std::cout<<"e8_small = "<<e8_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
#if (PART1 != 1) || !defined(TISCOMPLEX) || defined(MKL)
        for (int k=0; k<nloops2; ++k) { B3[k] = A3[k]; }
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LAPNAME(getrf) (
                LAPCM N,N,LP(B3[k].ptr()),N,&P3[k][0],&lapinfo);
            int lwork = N*LAP_BLOCKSIZE;
            std::vector<RT> work(lwork);
            LAPNAME(getri) (
                LAPCM N,LP(B3[k].ptr()),N,&P3[k][0],
                LP(&work[0]),lwork,&lapinfo);
#elif PART1 >= 2 && PART1 <= 5
            B3[k].setToIdentity();
            BLASNAME(trsm) (
                BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
                N,N,one,BP(A3[k].cptr()),N,BP(B3[k].ptr()),N
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                A0[k] = MPART1(A3[k]) * MPART1(B3[k]);
                e8_blas += NormSq(A0[k]-1);
            }
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= Norm(A1[0]);
            //std::cout<<"e8_blas = "<<e8_blas<<std::endl;
        }
#endif
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) { B4[k] = A4[k]; }
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            B4[k].transpose() = B4[k].transpose().inverse();
#elif PART1 >= 2 && PART1 <= 5
            B4[k].transpose() = EPART1(A4[k].transpose()).solveTriangular(B4[k].transpose().Identity(N,N));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A4[k]) * EPART1(B4[k]);
                e8_eigen += (temp1-temp1.Identity(N,N)).squaredNorm();
            }
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= Norm(A1[0]);
            //std::cout<<"e8_eigen = "<<e8_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2; ++k) { B5[k] = A5[k]; }
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            B5[k].transpose() = B5[k].transpose().inverse().eval();
#elif PART1 >= 2 && PART1 <= 5
            B5[k].transpose() = EPART1(A5[k].transpose()).solveTriangular(B5[k].transpose().Identity());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = EPART1(A5[k]) * EPART1(B5[k]);
                e8_smalleigen += (temp1-B5[k].Identity()).squaredNorm();
            }
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= Norm(A1[0]);
            //std::cout<<"e8_smalleigen = "<<e8_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // F %= A
        ClearCache();

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) { F1[k] = E1[k]; }
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F1[k] %= MPART1(A1[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F1[k] * MPART1(A1[k]);
                e9_reg += NormSq(E0[k]-E1[k]);
            }
            e9_reg = sqrt(e9_reg/nloops2);
            e9_reg /= Norm(A1[0]);
            //std::cout<<"e9_reg = "<<e9_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) F2[k] = E2[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F2[k] %= MPART1(A2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F2[k] * MPART1(A2[k]);
                e9_small += NormSq(E0[k]-E2[k]);
            }
            e9_small = sqrt(e9_small/nloops2);
            e9_small /= Norm(A1[0]);
            //std::cout<<"e9_small = "<<e9_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        for (int k=0; k<nloops2; ++k) F3[k] = E3[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LU3[k] = A3[k];
            LAPNAME(getrf) (
                LAPCM N,N,LP(LU3[k].ptr()),N,&P3[k][0],&lapinfo);
            // This only works for RowMajor F3, so need a temporary.
            tmv::Matrix<T,tmv::RowMajor> temp = F3[k];
            LAPNAME(getrs) (
                LAPCM LAPT,N,K,LP(LU3[k].ptr()),N,&P3[k][0],
                LP(temp.ptr()),N,&lapinfo);
            F3[k] = temp;
#elif PART1 >= 2 && PART1 <= 5
            BLASNAME(trsm) (
                BLASCM BLASRight, BLASUPLO1x, BLAST, BLASDIAG1,
                K,N,one,BP(A3[k].cptr()),N,BP(F3[k].ptr()),K
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F3[k] * MPART1(A3[k]);
                e9_blas += NormSq(E0[k]-E3[k]);
            }
            e9_blas = sqrt(e9_blas/nloops2);
            e9_blas /= Norm(A1[0]);
            //std::cout<<"e9_blas = "<<e9_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) F4[k] = E4[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            F4[k].transpose() = A4[k].lu().solve(F4[k].transpose());
#elif PART1 >= 2 && PART1 <= 5
            F4[k].transpose() = 
                EPART1(A4[k].transpose()).transpose().solveTriangular(F4[k].transpose());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = F4[k] * EPART1(A4[k].transpose());
                e9_eigen += (E4[k]-temp1).squaredNorm();
            }
            e9_eigen = sqrt(e9_eigen/nloops2);
            e9_eigen /= Norm(A1[0]);
            //std::cout<<"e9_eigen = "<<e9_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) F5[k] = E5[k];
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            F5[k] *= A5[k].transpose().inverse();
#elif PART1 >= 2 && PART1 <= 5
            EPART1(A5[k].transpose()).transpose().solveTriangularInPlace(F5[k].transpose());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e9_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = F5[k] * EPART1(A5[k].transpose());
                e9_smalleigen += (E5[k]-temp1).squaredNorm();
            }
            e9_smalleigen = sqrt(e9_smalleigen/nloops2);
            e9_smalleigen /= Norm(A1[0]);
            //std::cout<<"e9_smalleigen = "<<e9_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // F = E % A
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F1[k] = E1[k] % MPART1(A1[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F1[k] * MPART1(A1[k]);
                e10_reg += NormSq(E0[k]-E1[k]);
            }
            e10_reg = sqrt(e10_reg/nloops2);
            e10_reg /= Norm(A1[0]);
            //std::cout<<"e10_reg = "<<e10_reg<<std::endl;
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            F2[k] = E2[k] % MPART1(A2[k]);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_small = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F2[k] * MPART1(A2[k]);
                e10_small += NormSq(E0[k]-E2[k]);
            }
            e10_small = sqrt(e10_small/nloops2);
            e10_small /= Norm(A1[0]);
            //std::cout<<"e10_small = "<<e10_small<<std::endl;
        }
#endif
#endif

#ifdef DOLAP
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            LU3[k] = A3[k];
            LAPNAME(getrf) (
                LAPCM N,N,LP(LU3[k].ptr()),N,&P3[k][0],&lapinfo);
            tmv::Matrix<T,tmv::RowMajor> temp = E3[k];
            LAPNAME(getrs) (
                LAPCM LAPT,N,K,LP(LU3[k].ptr()),N,&P3[k][0],
                LP(temp.ptr()),N,&lapinfo);
            F3[k] = temp;
#elif PART1 >= 2 && PART1 <= 5
            F3[k] = E3[k];
            BLASNAME(trsm) (
                BLASCM BLASRight, BLASUPLO1x, BLAST, BLASDIAG1,
                K,N,one,BP(A3[k].cptr()),N,BP(F3[k].ptr()),K
                BLAS1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_blas = 0.;
            for (int k=0; k<nloops2; ++k) {
                E0[k] = F3[k] * MPART1(A3[k]);
                e10_blas += NormSq(E0[k]-E3[k]);
            }
            e10_blas = sqrt(e10_blas/nloops2);
            e10_blas /= Norm(A1[0]);
            //std::cout<<"e10_blas = "<<e10_blas<<std::endl;
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if PART1 == 1
            F4[k].transpose() = A4[k].lu().solve(E4[k].transpose());
#elif PART1 >= 2 && PART1 <= 5
            F4[k].transpose() =
                EPART1(A4[k].transpose()).transpose().solveTriangular(E4[k].transpose());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = F4[k] * EPART1(A4[k].transpose());
                e10_eigen += (E4[k]-temp1).squaredNorm();
            }
            e10_eigen = sqrt(e10_eigen/nloops2);
            e10_eigen /= Norm(A1[0]);
            //std::cout<<"e10_eigen = "<<e10_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
#if PART1 == 1
            F5[k] = E5[k] * A5[k].transpose().inverse();
#elif PART1 >= 2 && PART1 <= 5
            F5[k].transpose() = 
                EPART1(A5[k].transpose()).transpose().solveTriangular(E5[k].transpose());
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e10_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                EIGENM temp1 = F5[k] * EPART1(A5[k].transpose());
                e10_smalleigen += (E5[k]-temp1).squaredNorm();
            }
            e10_smalleigen = sqrt(e10_smalleigen/nloops2);
            e10_smalleigen /= Norm(A1[0]);
            //std::cout<<"e10_smalleigen = "<<e10_smalleigen<<std::endl;
        }
#endif
#endif
#endif

#if 1 // d = A.det()
#ifdef ERRORCHECK
        for (int k=0; k<nloops2; ++k) {
            // I just use the TMV calculation for the error check.
            // So if there is an error there, it only turns up in comparison
            // with the other 4 answers.
            d0[k] = MPART1(A1[k]).det();
        }
#endif
        ClearCache();

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            d1[k] = MPART1(A1[k]).det();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_reg = 0.;
            for (int k=0; k<nloops2; ++k) {
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
            d2[k] = MPART1(A2[k]).det();
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
            d3[k] = T(1);
#if PART1 == 1
            LU3[k] = A3[k];
            LAPNAME(getrf) (
                LAPCM N,N,LP(LU3[k].ptr()),N,&P3[k][0],&lapinfo);
            for(int i=0;i<N;++i) d3[k] *= LU3[k](i,i);
            for(int i=0;i<N;++i) if (P3[k][i] != i+1) d3[k] = -d3[k];
#elif PART1 >= 2 && PART1 <= 5
            for(int i=0;i<N;++i) d3[k] *= A3[k](i,i);
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
            d4[k] = A4[k].transpose().determinant();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e12_eigen += tmv::TMV_NORM(d0[k]-d4[k]);
            }
            e12_eigen = sqrt(e12_eigen/nloops2);
            e12_eigen /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e12_eigen = "<<e12_eigen<<std::endl;
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            d5[k] = A5[k].transpose().determinant();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e12_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                e12_smalleigen += tmv::TMV_NORM(d0[k]-d5[k]);
            }
            e12_smalleigen = sqrt(e12_smalleigen/nloops2);
            e12_smalleigen /= tmv::TMV_ABS(d0[0]);
            //std::cout<<"e12_smalleigen = "<<e12_smalleigen<<std::endl;
        }
#endif
#endif
#endif
#ifdef SHOW_DOTS
        std::cout<<"."; std::cout.flush();
#endif
    }
    std::cout<<"\n";

#ifdef DO_SEPARATE_DECOMP
    std::cout<<"A -> PLU               "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"No decomposition:\n";
    std::cout<<"D /= A                 "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"D = C / A              "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"B = A.inverse()        "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"d = A.det()            "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
    std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
    std::cout<<"Including decomposition:\n";
#endif
    std::cout<<"D /= A                 "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
    std::cout<<"D = C / A              "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"B = A.inverse()        "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"B = B.inverse()        "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
    std::cout<<"F %= A                 "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
    std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
    std::cout<<"F = E % A              "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
    std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
    std::cout<<"d = A.det()            "<<t12_reg<<"  "<<t12_small<<"  "<<t12_blas;
    std::cout<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
#ifdef DO_SEPARATE_DECOMP
    std::cout<<"A -> PLU               "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"D /= A ND              "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"D = C / A ND           "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"B = A.inverse() ND     "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"d = A.det() ND         "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
    std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
#endif
    std::cout<<"D /= A                 "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
    std::cout<<"D = C / A              "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
    std::cout<<"B = A.inverse()        "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
    std::cout<<"B = B.inverse()        "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
    std::cout<<"F %= A                 "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
    std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
    std::cout<<"F = E % A              "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
    std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
    std::cout<<"d = A.det()            "<<e12_reg<<"  "<<e12_small<<"  "<<e12_blas;
    std::cout<<"  "<<e12_eigen<<"  "<<e12_smalleigen<<std::endl;
    std::cout<<"\n\n";
#endif
#ifdef DOREG
#ifdef TMV_NO_LIB
    for(int i=0;i<nloops2;++i) if (LU1[i]) delete LU1[i];
#endif
#endif
#ifdef DOSMALL
    for(int i=0;i<nloops2;++i) if (LU2[i]) delete LU2[i];
#endif
#ifdef DOEIGEN
    for(int i=0;i<nloops2;++i) if (LU4[i]) delete LU4[i];
#endif
#ifdef DOEIGENSMALL
    for(int i=0;i<nloops2x;++i) if (LU5[i]) delete LU5[i];
#endif
}
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

    std::vector<tmv::Matrix<T> > AC(nloops2,tmv::Matrix<T>(N,N));
    std::vector<tmv::Matrix<T,tmv::RowMajor> > AR(nloops2,tmv::Matrix<T,tmv::RowMajor>(N,N));
    std::vector<tmv::Matrix<T> > B(nloops2,tmv::Matrix<T>(N,N));
    std::vector<tmv::Matrix<T> > C(nloops2,tmv::Matrix<T>(N,K));
    std::vector<tmv::Matrix<T> > D(nloops2,tmv::Matrix<T>(N,K));
    std::vector<tmv::Matrix<T> > E(nloops2,tmv::Matrix<T>(K,N));
    std::vector<tmv::Matrix<T> > F(nloops2,tmv::Matrix<T>(K,N));
    for(int k=0;k<nloops2;++k) {
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) {
            AC[k](i,j) = RAND;
            B[k](j,i) = RAND;
#if 0
            if (i+2*j == N || 3*i+2*j == 2*N) {
                A[k](i,j) *= T(50);
                B[k](j,i) *= T(60);
            }
#endif
        }
#ifdef SMALL_OFFDIAG
        AC[k].diag().addToAll(1.);
        AC[k].upperTri().offDiag() /= 10000.;
        AC[k].lowerTri().offDiag() /= 10000.;
        B[k].diag().addToAll(1.);
        B[k].upperTri().offDiag() /= 10000.;
        B[k].lowerTri().offDiag() /= 10000.;
#endif
        for(int i=0;i<N;++i) for(int j=0;j<K;++j) {
            C[k](i,j) = RAND;
            E[k](j,i) = RAND;
        }
        D[k].setZero();
        F[k].setZero();
        AR[k] = AC[k];
    }

    std::cout<<"N,K = "<<N<<" , "<<K<<std::endl;
    std::cout<<"nloops = "<<nloops1<<" x "<<nloops2;
    std::cout<<" = "<<nloops1*nloops2<<std::endl;
#ifdef DOEIGENSMALL
    if (nloops2x == 0)
        std::cout<<"Matrix is too big for the stack, so no \"Eigen Known\" tests.\n";
#endif
    std::cout<<"\nTime for:               ";
    std::cout<<" TMV    TMV Small  LAPACK   Eigen Eigen Known"<<std::endl;

#ifdef DO_C
    LU_C(AC,B,C,D,E,F);
#endif

    for(int k=0;k<nloops2;++k) { D[k].setZero(); F[k].setZero(); }

#ifdef DO_R
    LU_R(AR,B,C,D,E,F);
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
