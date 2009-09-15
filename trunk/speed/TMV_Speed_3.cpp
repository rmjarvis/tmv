
// How big do you want the matrices to be?
#if 0
const int M = 32;
const int N = 32;
const int K = 32;
#elif 0
const int M = 979;
const int N = 949;
const int K = 999;
#else
const int M = 4;
const int N = 4;
const int K = 64;
#endif

// Define the type to use:
//#define TISFLOAT
#define TISCOMPLEX

// Define the parts of the matrices to use
// 1 = all (rectangle matrix)
// 2 = uppertri
// 3 = unituppertri
// 4 = lowertri
// 5 = unitlowertri
#define PART1 1
#define PART2 1
#define PART3 1

// Define whether you want to include the error checks.
//#define XDEBUG_PRODMM
#define ERRORCHECK

// Define which versions you want to test:
#define DOREG
//#define DOSMALL
//#define DOBLAS
//#define DOEIGEN
//#define DOEIGENSMALL

// Define which batches of functions you want to test:
//#define DOMULTMM_CCC
//#define DOMULTMM_RCC
//#define DOMULTMM_CRC
#define DOMULTMM_RRC

// Set this if you only want to do a single loop.
// Not so useful for timing, but useful for debugging.
// Also useful if M,N,K are very large to avoid overflow in product
//#define ONELOOP
//#define PRINTALGO_MM
//#define PRINTALGO_MM_BLOCK
//#define PRINTALGO_MM_OMP

// Normally I fill the matrices and vectors with random values, but 
// that can make it difficult to debug the arithmetic.
// This define uses simpler values that I can multiply in my head.
//#define SIMPLE_VALUES

// Set up the target number of operations and memory to use for testing
const int targetnflops = 100000000; // in real ops
const int targetmem = 10000000; // in bytes

// Include the BLAS library you want to test TMV against:
#if 0
#include "mkl.h"
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
const int nloops1 = 1;
const int nloops2 = 1;
#else
const int nloops2 = (
    (M*N*K*sizeof(T) > targetmem ? 1 : targetmem / (M*N*K) / (sizeof(T))));
const int nloops1 = (
    (M*N*K*nloops2*XFOUR > targetnflops ? 1 :
     targetnflops / (M*N*K*nloops2) / XFOUR ));
#endif

#if (PART1 == 1) // Full rectangle
#define INPART1(i,j) true
#define UNITPART1(i,j) false
#define MPART1(m) m
#define EPART1(m) m
#define COL_LEN1(j) M
#define COL_START1(j) 0
#define COL_END1(j) M
#define ROW_LEN1(i) N
#define ROW_START1(i) 0
#define ROW_END1(i) N

#elif (PART1 == 2) // UpperTri
#define INPART1(i,j) (i<=j)
#define UNITPART1(i,j) false
#define MPART1(m) m.UpperTri()
#define EPART1(m) m.part<Eigen::UpperTriangular>()
#define COL_LEN1(j) j+1
#define COL_START1(j) 0
#define COL_END1(j) j+1
#define ROW_LEN1(i) N-i
#define ROW_START1(i) i
#define ROW_END1(i) N

#elif (PART1 == 3) // UnitUpperTri
#define INPART1(i,j) (i<j)
#define UNITPART1(i,j) (i==j)
#define MPART1(m) m.UnitUpperTri()
#define EPART1(m) m.part<Eigen::UnitUpperTriangular>()
#define COL_LEN1(j) j
#define COL_START1(j) 0
#define COL_END1(j) j
#define ROW_LEN1(i) N-i-1
#define ROW_START1(i) i+1
#define ROW_END1(i) N

#elif (PART1 == 4) // LowerTri
#define INPART1(i,j) (i>=j)
#define UNITPART1(i,j) false
#define MPART1(m) m.LowerTri()
#define EPART1(m) m.part<Eigen::LowerTriangular>()
#define COL_LEN1(j) N-j
#define COL_START1(j) j
#define COL_END1(j) N
#define ROW_LEN1(i) i+1
#define ROW_START1(i) 0
#define ROW_END1(i) i+1

#elif (PART1 == 5) // UnitLowerTri
#define INPART1(i,j) (i>j)
#define UNITPART1(i,j) (i==j)
#define MPART1(m) m.UnitLowerTri()
#define EPART1(m) m.part<Eigen::UnitLowerTriangular>()
#define COL_LEN1(j) N-j-1
#define COL_START1(j) j+1
#define COL_END1(j) N
#define ROW_LEN1(i) i
#define ROW_START1(i) 0
#define ROW_END1(i) i

#endif

#if (PART2 == 1) // Full rectangle
#define INPART2(i,j) true
#define UNITPART2(i,j) false
#define MPART2(m) m
#define EPART2(m) m
#define COL_LEN2(j) M
#define COL_START2(j) 0
#define COL_END2(j) M
#define ROW_LEN2(i) N
#define ROW_START2(i) 0
#define ROW_END2(i) N

#elif (PART2 == 2) // UpperTri
#define INPART2(i,j) (i<=j)
#define UNITPART2(i,j) false
#define MPART2(m) m.UpperTri()
#define EPART2(m) m.part<Eigen::UpperTriangular>()
#define COL_LEN2(j) j+1
#define COL_START2(j) 0
#define COL_END2(j) j+1
#define ROW_LEN2(i) N-i
#define ROW_START2(i) i
#define ROW_END2(i) N

#elif (PART2 == 3) // UnitUpperTri
#define INPART2(i,j) (i<j)
#define UNITPART2(i,j) (i==j)
#define MPART2(m) m.UnitUpperTri()
#define EPART2(m) m.part<Eigen::UnitUpperTriangular>()
#define COL_LEN2(j) j
#define COL_START2(j) 0
#define COL_END2(j) j
#define ROW_LEN2(i) N-i-1
#define ROW_START2(i) i+1
#define ROW_END2(i) N

#elif (PART2 == 4) // LowerTri
#define INPART2(i,j) (i>=j)
#define UNITPART2(i,j) false
#define MPART2(m) m.LowerTri()
#define EPART2(m) m.part<Eigen::LowerTriangular>()
#define COL_LEN2(j) N-j
#define COL_START2(j) j
#define COL_END2(j) N
#define ROW_LEN2(i) i+1
#define ROW_START2(i) 0
#define ROW_END2(i) i+1

#elif (PART2 == 5) // UnitLowerTri
#define INPART2(i,j) (i>j)
#define UNITPART2(i,j) (i==j)
#define MPART2(m) m.UnitLowerTri()
#define EPART2(m) m.part<Eigen::UnitLowerTriangular>()
#define COL_LEN2(j) N-j-1
#define COL_START2(j) j+1
#define COL_END2(j) N
#define ROW_LEN2(i) i
#define ROW_START2(i) 0
#define ROW_END2(i) i

#endif

#if (PART3 == 1) // Full rectangle
#define INPART3(i,j) true
#define UNITPART3(i,j) false
#define MPART3(m) m
#define EPART3(m) m
#define COL_LEN3(j) M
#define COL_START3(j) 0
#define COL_END3(j) M
#define ROW_LEN3(i) N
#define ROW_START3(i) 0
#define ROW_END3(i) N

#elif (PART3 == 2) // UpperTri
#define INPART3(i,j) (i<=j)
#define UNITPART3(i,j) false
#define MPART3(m) m.UpperTri()
#define EPART3(m) m.part<Eigen::UpperTriangular>()
#define COL_LEN3(j) j+1
#define COL_START3(j) 0
#define COL_END3(j) j+1
#define ROW_LEN3(i) N-i
#define ROW_START3(i) i
#define ROW_END3(i) N

#elif (PART3 == 3) // UnitUpperTri
#define INPART3(i,j) (i<j)
#define UNITPART3(i,j) (i==j)
#define MPART3(m) m.UnitUpperTri()
#define EPART3(m) m.part<Eigen::UnitUpperTriangular>()
#define COL_LEN3(j) j
#define COL_START3(j) 0
#define COL_END3(j) j
#define ROW_LEN3(i) N-i-1
#define ROW_START3(i) i+1
#define ROW_END3(i) N

#elif (PART3 == 4) // LowerTri
#define INPART3(i,j) (i>=j)
#define UNITPART3(i,j) false
#define MPART3(m) m.LowerTri()
#define EPART3(m) m.part<Eigen::LowerTriangular>()
#define COL_LEN3(j) N-j
#define COL_START3(j) j
#define COL_END3(j) N
#define ROW_LEN3(i) i+1
#define ROW_START3(i) 0
#define ROW_END3(i) i+1

#elif (PART3 == 5) // UnitLowerTri
#define INPART3(i,j) (i>j)
#define UNITPART3(i,j) (i==j)
#define MPART3(m) m.UnitLowerTri()
#define EPART3(m) m.part<Eigen::UnitLowerTriangular>()
#define COL_LEN3(j) N-j-1
#define COL_START3(j) j+1
#define COL_END3(j) N
#define ROW_LEN3(i) i
#define ROW_START3(i) 0
#define ROW_END3(i) i

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

const int lap_worksize = (M > N ? M * 128 : N * 128);
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
#else
UNKNOWN LAP definition....  // Give compile error
#define LAPNAME1(c,x) c ## x
#endif

#define LAPNAME2(c,x) LAPNAME1(c,x)
#define LAPNAME(x) LAPNAME2(BLAS_LETTER,x)

#endif // DOBLAS

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

#define EIGENSMA Eigen::Matrix<T,M,K>
#define EIGENSMB Eigen::Matrix<T,K,M>
#define EIGENSMC Eigen::Matrix<T,K,N>
#define EIGENSMD Eigen::Matrix<T,N,K>
#define EIGENSME Eigen::Matrix<T,M,N>
#define ALLOC(X) Eigen::aligned_allocator< X >

#endif

#ifdef DOEIGENSMALL
const int nloops2x = (
    (M*N*sizeof(T) < 1536 * 1024 && 
     M*K*sizeof(T) < 1536 * 1024 && 
     K*N*sizeof(T) < 1536 * 1024 && 
     N!=10000 && M!=10000 && K!=10000) ? 
    nloops2 : 0 );
#endif

#include <sys/time.h>
#include <algorithm>
#include <numeric>
#include <iostream>
#include "TMV.h"


#if (defined(DOEIGEN) || defined(DOEIGENSMALL))

// Supposedly, these are required when using vector's of Eigen types,
// but they don't seem to work.  
// And leaving them out does seem to work.
//#define EIGEN_USE_NEW_STDVECTOR
//#include "Eigen/StdVector"
#include "Eigen/Core"
#include "Eigen/Array"

#endif

static void ClearCache()
{
  tmv::Vector<double> X(100000,8.);
  tmv::Vector<double> Y(100000,8.);
  tmv::Vector<double> Z(100000);
  Z = X + Y;
}

#ifdef DOMULTMM_CCC
static void MultMM_CCC(
    const std::vector<tmv::Matrix<T> >& A1,
    const std::vector<tmv::Matrix<T> >& C1,
    std::vector<tmv::Matrix<T> >& E1)
{
#ifdef ERRORCHECK
  std::vector<tmv::Matrix<T> > A0 = A1;
  std::vector<tmv::Matrix<T> > C0 = C1;
  std::vector<tmv::Matrix<T> > E0 = E1;
#endif

#ifdef DOSMALL
  std::vector<tmv::SmallMatrix<T,M,K> > A2(nloops2);
  std::vector<tmv::SmallMatrix<T,K,N> > C2(nloops2);
  std::vector<tmv::SmallMatrix<T,M,N> > E2(nloops2);

  for(int k=0; k<nloops2; ++k) {
    A2[k] = A1[k]; C2[k] = C1[k]; E2[k] = E1[k];
  }
#endif

#ifdef DOBLAS
  std::vector<tmv::Matrix<T> > A3 = A1;
  std::vector<tmv::Matrix<T> > C3 = C1;
  std::vector<tmv::Matrix<T> > E3 = E1;
#endif

#ifdef DOEIGEN
  std::vector<EIGENM,ALLOC(EIGENM) > A4(nloops2,EIGENM(M,K));
  std::vector<EIGENM,ALLOC(EIGENM) > C4(nloops2,EIGENM(K,N));
  std::vector<EIGENM,ALLOC(EIGENM) > E4(nloops2,EIGENM(M,N));
  for(int k=0;k<nloops2;++k) {
    for(int j=0;j<K;++j) for(int i=0;i<M;++i) A4[k](i,j) = A1[k](i,j);
    for(int j=0;j<N;++j) for(int i=0;i<K;++i) C4[k](i,j) = C1[k](i,j);
  }
#endif

#ifdef DOEIGENSMALL
  std::vector<EIGENSMA,ALLOC(EIGENSMA) > A5;
  std::vector<EIGENSMC,ALLOC(EIGENSMC) > C5;
  std::vector<EIGENSME,ALLOC(EIGENSME) > E5;
  if (nloops2x) { 
    A5.resize(nloops2x); C5.resize(nloops2x);
    E5.resize(nloops2x); 
  }
  for(int k=0;k<nloops2x;++k) {
    for(int j=0;j<K;++j) for(int i=0;i<M;++i) A5[k](i,j) = A1[k](i,j);
    for(int j=0;j<N;++j) for(int i=0;i<K;++i) C5[k](i,j) = C1[k](i,j);
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
#endif

  for (int n=0; n<nloops1; ++n) {

#if 1 // E = A * C
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = MPART1(A0[n2]) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = A0[n2].row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = MPART1(A1[k]) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_reg = 0.;
      for (int k=0; k<nloops2; ++k) e1_reg += NormSq(E1[k]-E0[k]);
      e1_reg = sqrt(e1_reg/nloops2);
      e1_reg /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = MPART1(A2[k]) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_small = 0.;
      for (int k=0; k<nloops2; ++k) e1_small += NormSq(E2[k]-E0[k]);
      e1_small = sqrt(e1_small/nloops2);
      e1_small /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLASNT,
	  M,N,K,one,BP(A3[k].cptr()),M,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,one,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,one,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_blas = 0.;
      for (int k=0; k<nloops2; ++k) e1_blas += NormSq(E3[k]-E0[k]);
      e1_blas = sqrt(e1_blas/nloops2);
      e1_blas /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = EPART1(A4[k]) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e1_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e1_eigen = sqrt(e1_eigen/nloops2);
      e1_eigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = EPART1(A5[k]) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e1_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e1_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e1_smalleigen = sqrt(e1_smalleigen/nloops2);
      e1_smalleigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif

#if 1 // E = -A * C
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = -MPART1(A0[n2]) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = -A0[n2].row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = -MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = -MPART1(A1[k]) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_reg = 0.;
      for (int k=0; k<nloops2; ++k) e2_reg += NormSq(E1[k]-E0[k]);
      e2_reg = sqrt(e2_reg/nloops2);
      e2_reg /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = -MPART1(A2[k]) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_small = 0.;
      for (int k=0; k<nloops2; ++k) e2_small += NormSq(E2[k]-E0[k]);
      e2_small = sqrt(e2_small/nloops2);
      e2_small /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLASNT,
	  M,N,K,mone,BP(A3[k].cptr()),M,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,mone,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,mone,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_blas = 0.;
      for (int k=0; k<nloops2; ++k) e2_blas += NormSq(E3[k]-E0[k]);
      e2_blas = sqrt(e2_blas/nloops2);
      e2_blas /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = -EPART1(A4[k]) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e2_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e2_eigen = sqrt(e2_eigen/nloops2);
      e2_eigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = -EPART1(A5[k]) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e2_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e2_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e2_smalleigen = sqrt(e2_smalleigen/nloops2);
      e2_smalleigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif

#if 1 // E = 7 * A * C
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = RT(7) * MPART1(A0[n2]) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = RT(7) * A0[n2].row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = RT(7) * MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = 7 * MPART1(A1[k]) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_reg = 0.;
      for (int k=0; k<nloops2; ++k) e3_reg += NormSq(E1[k]-E0[k]);
      e3_reg = sqrt(e3_reg/nloops2);
      e3_reg /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = 7 * MPART1(A2[k]) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_small = 0.;
      for (int k=0; k<nloops2; ++k) e3_small += NormSq(E2[k]-E0[k]);
      e3_small = sqrt(e3_small/nloops2);
      e3_small /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLASNT,
	  M,N,K,seven,BP(A3[k].cptr()),M,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,seven,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,seven,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_blas = 0.;
      for (int k=0; k<nloops2; ++k) e3_blas += NormSq(E3[k]-E0[k]);
      e3_blas = sqrt(e3_blas/nloops2);
      e3_blas /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = 7 * EPART1(A4[k]) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e3_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e3_eigen = sqrt(e3_eigen/nloops2);
      e3_eigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = 7 * EPART1(A5[k]) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e3_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e3_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e3_smalleigen = sqrt(e3_smalleigen/nloops2);
      e3_smalleigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif

#if 1 // E -= A * C
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) -= MPART1(A0[n2]) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) -= A0[n2].row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) -= MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) -= MPART1(A1[k]) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_reg = 0.;
      for (int k=0; k<nloops2; ++k) e4_reg += NormSq(E1[k]-E0[k]);
      e4_reg = sqrt(e4_reg/nloops2);
      e4_reg /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) -= MPART1(A2[k]) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_small = 0.;
      for (int k=0; k<nloops2; ++k) e4_small += NormSq(E2[k]-E0[k]);
      e4_small = sqrt(e4_small/nloops2);
      e4_small /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLASNT,
	  M,N,K,mone,BP(A3[k].cptr()),M,BP(C3[k].cptr()),K,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5)
      tmv::Matrix<T> temp = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,one,BP(A3[k].cptr()),M,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) -= temp;
#elif (PART2 >= 2 && PART2 <= 5)
      tmv::Matrix<T> temp = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,one,BP(C3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) -= temp;
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_blas = 0.;
      for (int k=0; k<nloops2; ++k) e4_blas += NormSq(E3[k]-E0[k]);
      e4_blas = sqrt(e4_blas/nloops2);
      e4_blas /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) -= EPART1(A4[k]) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e4_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e4_eigen = sqrt(e4_eigen/nloops2);
      e4_eigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) -= EPART1(A5[k]) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e4_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e4_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e4_smalleigen = sqrt(e4_smalleigen/nloops2);
      e4_smalleigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif

#if 1 // E += 8 * A * C
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += RT(8) * MPART1(A0[n2]) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += RT(8) * A0[n2].row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = RT(8) * MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += 8 * MPART1(A1[k]) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_reg = 0.;
      for (int k=0; k<nloops2; ++k) e5_reg += NormSq(E1[k]-E0[k]);
      e5_reg = sqrt(e5_reg/nloops2);
      e5_reg /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += 8 * MPART1(A2[k]) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_small = 0.;
      for (int k=0; k<nloops2; ++k) e5_small += NormSq(E2[k]-E0[k]);
      e5_small = sqrt(e5_small/nloops2);
      e5_small /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLASNT,
	  M,N,K,eight,BP(A3[k].cptr()),M,BP(C3[k].cptr()),K,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5)
      tmv::Matrix<T> temp = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,eight,BP(A3[k].cptr()),M,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += temp;
#elif (PART2 >= 2 && PART2 <= 5)
      tmv::Matrix<T> temp = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,eight,BP(C3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += temp;
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_blas = 0.;
      for (int k=0; k<nloops2; ++k) e5_blas += NormSq(E3[k]-E0[k]);
      e5_blas = sqrt(e5_blas/nloops2);
      e5_blas /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += 8 * EPART1(A4[k]) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e5_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e5_eigen = sqrt(e5_eigen/nloops2);
      e5_eigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += 8 * EPART1(A5[k]) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e5_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e5_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e5_smalleigen = sqrt(e5_smalleigen/nloops2);
      e5_smalleigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif

#if 1 // E = (7,1) * A * C
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(A0[n2]) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * A0[n2].row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(A1[k]) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_reg = 0.;
      for (int k=0; k<nloops2; ++k) e6_reg += NormSq(E1[k]-E0[k]);
      e6_reg = sqrt(e6_reg/nloops2);
      e6_reg /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(A2[k]) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_small = 0.;
      for (int k=0; k<nloops2; ++k) e6_small += NormSq(E2[k]-E0[k]);
      e6_small = sqrt(e6_small/nloops2);
      e6_small /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLASNT,
	  M,N,K,z71,BP(A3[k].cptr()),M,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z71,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z71,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_blas = 0.;
      for (int k=0; k<nloops2; ++k) e6_blas += NormSq(E3[k]-E0[k]);
      e6_blas = sqrt(e6_blas/nloops2);
      e6_blas /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(A4[k]) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e6_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e6_eigen = sqrt(e6_eigen/nloops2);
      e6_eigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(A5[k]) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e6_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e6_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e6_smalleigen = sqrt(e6_smalleigen/nloops2);
      e6_smalleigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * A * C
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(A0[n2]) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * A0[n2].row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(A1[k]) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_reg = 0.;
      for (int k=0; k<nloops2; ++k) e7_reg += NormSq(E1[k]-E0[k]);
      e7_reg = sqrt(e7_reg/nloops2);
      e7_reg /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(A2[k]) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_small = 0.;
      for (int k=0; k<nloops2; ++k) e7_small += NormSq(E2[k]-E0[k]);
      e7_small = sqrt(e7_small/nloops2);
      e7_small /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLASNT,
	  M,N,K,z89,BP(A3[k].cptr()),M,BP(C3[k].cptr()),K,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5)
      tmv::Matrix<T> temp = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z89,BP(A3[k].cptr()),M,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += temp;
#elif (PART2 >= 2 && PART2 <= 5)
      tmv::Matrix<T> temp = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z89,BP(C3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += temp;
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_blas = 0.;
      for (int k=0; k<nloops2; ++k) e7_blas += NormSq(E3[k]-E0[k]);
      e7_blas = sqrt(e7_blas/nloops2);
      e7_blas /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(A4[k]) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e7_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e7_eigen = sqrt(e7_eigen/nloops2);
      e7_eigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(A5[k]) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e7_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e7_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e7_smalleigen = sqrt(e7_smalleigen/nloops2);
      e7_smalleigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E = (7,1) * A* * C
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(A0[n2]).Conjugate() * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * A0[n2].row(i).Conjugate() * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(A0[n2]).Conjugate() * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(A1[k]).Conjugate() * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_reg = 0.;
      for (int k=0; k<nloops2; ++k) e8_reg += NormSq(E1[k]-E0[k]);
      e8_reg = sqrt(e8_reg/nloops2);
      e8_reg /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(A2[k]).Conjugate() * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_small = 0.;
      for (int k=0; k<nloops2; ++k) e8_small += NormSq(E2[k]-E0[k]);
      e8_small = sqrt(e8_small/nloops2);
      e8_small /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASDNAME(scal) (M*K,dmone,A3[k].Real().ptr()+1,2);
      BLASNAME(gemm) (BLASCM BLASNT, BLASNT,
	  M,N,K,z71,BP(A3[k].cptr()),M,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (M*K,dmone,A3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      BLASDNAME(scal) (M*K,dmone,A3[k].Real().ptr()+1,2);
      MPART3(E3[k]) = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z71,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      BLASDNAME(scal) (M*K,dmone,A3[k].Real().ptr()+1,2);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(A3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z71,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_blas = 0.;
      for (int k=0; k<nloops2; ++k) e8_blas += NormSq(E3[k]-E0[k]);
      e8_blas = sqrt(e8_blas/nloops2);
      e8_blas /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(A4[k]).conjugate() * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e8_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e8_eigen = sqrt(e8_eigen/nloops2);
      e8_eigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(A5[k]).conjugate() * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e8_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e8_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e8_smalleigen = sqrt(e8_smalleigen/nloops2);
      e8_smalleigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * A* * C
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(A0[n2]).Conjugate() * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * A0[n2].row(i).Conjugate() * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(A0[n2]).Conjugate() * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(A1[k]).Conjugate() * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_reg = 0.;
      for (int k=0; k<nloops2; ++k) e9_reg += NormSq(E1[k]-E0[k]);
      e9_reg = sqrt(e9_reg/nloops2);
      e9_reg /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(A2[k]).Conjugate() * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_small = 0.;
      for (int k=0; k<nloops2; ++k) e9_small += NormSq(E2[k]-E0[k]);
      e9_small = sqrt(e9_small/nloops2);
      e9_small /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASDNAME(scal) (M*K,dmone,A3[k].Real().ptr()+1,2);
      BLASNAME(gemm) (BLASCM BLASNT, BLASNT,
	  M,N,K,z89,BP(A3[k].cptr()),M,BP(C3[k].cptr()),K,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (M*K,dmone,A3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5)
      tmv::Matrix<T> temp = MPART2(C3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z89c,BP(A3[k].cptr()),M,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += temp.Conjugate();
#elif (PART2 >= 2 && PART2 <= 5)
      tmv::Matrix<T> temp = MPART1(A3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z89,BP(C3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += temp;
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_blas = 0.;
      for (int k=0; k<nloops2; ++k) e9_blas += NormSq(E3[k]-E0[k]);
      e9_blas = sqrt(e9_blas/nloops2);
      e9_blas /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(A4[k]).conjugate() * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e9_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e9_eigen = sqrt(e9_eigen/nloops2);
      e9_eigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(A5[k]).conjugate() * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e9_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e9_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e9_smalleigen = sqrt(e9_smalleigen/nloops2);
      e9_smalleigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E = (7,1) * A * C*
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
	for(int j=0;j<N;++j) {
	  E0[n2].col(j) = T(7,1) * A0[n2] * C0[n2].col(j).Conjugate();
	}
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(A0[n2]) * C0[n2].col(j).Conjugate();
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * A0[n2].row(i) * MPART2(C0[n2]).Conjugate();
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(A0[n2]) * temp2.col(j).Conjugate();
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(A1[k]) * MPART2(C1[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_reg = 0.;
      for (int k=0; k<nloops2; ++k) e10_reg += NormSq(E1[k]-E0[k]);
      e10_reg = sqrt(e10_reg/nloops2);
      e10_reg /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(A2[k]) * MPART2(C2[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_small = 0.;
      for (int k=0; k<nloops2; ++k) e10_small += NormSq(E2[k]-E0[k]);
      e10_small = sqrt(e10_small/nloops2);
      e10_small /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
      BLASNAME(gemm) (BLASCM BLASNT, BLASNT,
	  M,N,K,z71,BP(A3[k].cptr()),M,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(C3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z71,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
      MPART3(E3[k]) = MPART1(A3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z71,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_blas = 0.;
      for (int k=0; k<nloops2; ++k) e10_blas += NormSq(E3[k]-E0[k]);
      e10_blas = sqrt(e10_blas/nloops2);
      e10_blas /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(A4[k]) * EPART2(C4[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e10_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e10_eigen = sqrt(e10_eigen/nloops2);
      e10_eigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(A5[k]) * EPART2(C5[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e10_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e10_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e10_smalleigen = sqrt(e10_smalleigen/nloops2);
      e10_smalleigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * A * C*
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(A0[n2]) * C0[n2].col(j).Conjugate();
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * A0[n2].row(i) * MPART2(C0[n2]).Conjugate();
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(A0[n2]) * temp2.col(j).Conjugate();
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(A1[k]) * MPART2(C1[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_reg = 0.;
      for (int k=0; k<nloops2; ++k) e11_reg += NormSq(E1[k]-E0[k]);
      e11_reg = sqrt(e11_reg/nloops2);
      e11_reg /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(A2[k]) * MPART2(C2[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_small = 0.;
      for (int k=0; k<nloops2; ++k) e11_small += NormSq(E2[k]-E0[k]);
      e11_small = sqrt(e11_small/nloops2);
      e11_small /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
      BLASNAME(gemm) (BLASCM BLASNT, BLASNT,
	  M,N,K,z89,BP(A3[k].cptr()),M,BP(C3[k].cptr()),K,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5)
      tmv::Matrix<T> temp = MPART2(C3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z89,BP(A3[k].cptr()),M,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += temp;
#elif (PART2 >= 2 && PART2 <= 5)
      tmv::Matrix<T> temp = MPART1(A3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z89c,BP(C3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += temp.Conjugate();
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_blas = 0.;
      for (int k=0; k<nloops2; ++k) e11_blas += NormSq(E3[k]-E0[k]);
      e11_blas = sqrt(e11_blas/nloops2);
      e11_blas /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(A4[k]) * EPART2(C4[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e11_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e11_eigen = sqrt(e11_eigen/nloops2);
      e11_eigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(A5[k]) * EPART2(C5[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e11_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e11_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e11_smalleigen = sqrt(e11_smalleigen/nloops2);
      e11_smalleigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E = (7,1) * A* * C*
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(A0[n2]).Conjugate() * C0[n2].col(j).Conjugate();
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * A0[n2].row(i).Conjugate() * MPART2(C0[n2]).Conjugate();
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(A0[n2]).Conjugate() * temp2.col(j).Conjugate();
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(A1[k]).Conjugate() * MPART2(C1[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_reg = 0.;
      for (int k=0; k<nloops2; ++k) e12_reg += NormSq(E1[k]-E0[k]);
      e12_reg = sqrt(e12_reg/nloops2);
      e12_reg /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(A2[k]).Conjugate() * MPART2(C2[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_small = 0.;
      for (int k=0; k<nloops2; ++k) e12_small += NormSq(E2[k]-E0[k]);
      e12_small = sqrt(e12_small/nloops2);
      e12_small /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLASNT,
	  M,N,K,z71c,BP(A3[k].cptr()),M,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (M*N,dmone,E3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z71c,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      BLASDNAME(scal) (M*N,dmone,E3[k].Real().ptr()+1,2);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z71c,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      BLASDNAME(scal) (M*N,dmone,E3[k].Real().ptr()+1,2);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_blas = 0.;
      for (int k=0; k<nloops2; ++k) e12_blas += NormSq(E3[k]-E0[k]);
      e12_blas = sqrt(e12_blas/nloops2);
      e12_blas /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(A4[k]).conjugate() * EPART2(C4[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e12_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e12_eigen = sqrt(e12_eigen/nloops2);
      e12_eigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(A5[k]).conjugate() * EPART2(C5[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e12_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e12_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e12_smalleigen = sqrt(e12_smalleigen/nloops2);
      e12_smalleigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * A* * C*
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(A0[n2]).Conjugate() * C0[n2].col(j).Conjugate();
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * A0[n2].row(i).Conjugate() * MPART2(C0[n2]).Conjugate();
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(A0[n2]).Conjugate() * temp2.col(j).Conjugate();
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(A1[k]).Conjugate() * MPART2(C1[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_reg = 0.;
      for (int k=0; k<nloops2; ++k) e13_reg += NormSq(E1[k]-E0[k]);
      e13_reg = sqrt(e13_reg/nloops2);
      e13_reg /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(A2[k]).Conjugate() * MPART2(C2[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_small = 0.;
      for (int k=0; k<nloops2; ++k) e13_small += NormSq(E2[k]-E0[k]);
      e13_small = sqrt(e13_small/nloops2);
      e13_small /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASDNAME(scal) (M*N,dmone,E3[k].Real().ptr()+1,2);
      BLASNAME(gemm) (BLASCM BLASNT, BLASNT,
	  M,N,K,z89c,BP(A3[k].cptr()),M,BP(C3[k].cptr()),K,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (M*N,dmone,E3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5)
      tmv::Matrix<T> temp = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z89c,BP(A3[k].cptr()),M,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += temp.Conjugate();
#elif (PART2 >= 2 && PART2 <= 5)
      tmv::Matrix<T> temp = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z89c,BP(C3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += temp.Conjugate();
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_blas = 0.;
      for (int k=0; k<nloops2; ++k) e13_blas += NormSq(E3[k]-E0[k]);
      e13_blas = sqrt(e13_blas/nloops2);
      e13_blas /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(A4[k]).conjugate() * EPART2(C4[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e13_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e13_eigen = sqrt(e13_eigen/nloops2);
      e13_eigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(A5[k]).conjugate() * EPART2(C5[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e13_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e13_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e13_smalleigen = sqrt(e13_smalleigen/nloops2);
      e13_smalleigen /= Norm(A0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif
  }

  std::cout<<"E = A * C               "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
  std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
  std::cout<<"E = -A * C              "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
  std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
  std::cout<<"E = 7 * A * C           "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
  std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
  std::cout<<"E -= A * C              "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
  std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
  std::cout<<"E += 8 * A * C          "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
  std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
#ifdef TISCOMPLEX
  std::cout<<"E = (7,1) * A * C       "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
  std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A * C      "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
  std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * A* * C      "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
  std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A* * C     "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
  std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * A * C*      "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
  std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A * C*     "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
  std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * A* * C*     "<<t12_reg<<"  "<<t12_small<<"  "<<t12_blas;
  std::cout<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A* * C*    "<<t13_reg<<"  "<<t13_small<<"  "<<t13_blas;
  std::cout<<"  "<<t13_eigen<<"  "<<t13_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
  std::cout<<"errors:\n";
  std::cout<<"E = A * C               "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
  std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
  std::cout<<"E = -A * C              "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
  std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
  std::cout<<"E = 7 * A * C           "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
  std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
  std::cout<<"E -= A * C              "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
  std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
  std::cout<<"E += 8 * A * C          "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
  std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
#ifdef TISCOMPLEX
  std::cout<<"E = (7,1) * A * C       "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
  std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A * C      "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
  std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * A* * C      "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
  std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A* * C     "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
  std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * A * C*      "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
  std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A * C*     "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
  std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * A* * C*     "<<e12_reg<<"  "<<e12_small<<"  "<<e12_blas;
  std::cout<<"  "<<e12_eigen<<"  "<<e12_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A* * C*    "<<e13_reg<<"  "<<e13_small<<"  "<<e13_blas;
  std::cout<<"  "<<e13_eigen<<"  "<<e13_smalleigen<<std::endl;
#endif
  std::cout<<"\n\n";
#endif
}
#endif

#ifdef DOMULTMM_RCC
static void MultMM_RCC(
    const std::vector<tmv::Matrix<T> >& B1,
    const std::vector<tmv::Matrix<T> >& C1,
    std::vector<tmv::Matrix<T> >& E1)
{
#ifdef ERRORCHECK
  std::vector<tmv::Matrix<T> > B0 = B1;
  std::vector<tmv::Matrix<T> > C0 = C1;
  std::vector<tmv::Matrix<T> > E0 = E1;
#endif

#ifdef DOSMALL
  std::vector<tmv::SmallMatrix<T,K,M> > B2(nloops2);
  std::vector<tmv::SmallMatrix<T,K,N> > C2(nloops2);
  std::vector<tmv::SmallMatrix<T,M,N> > E2(nloops2);

  for(int k=0; k<nloops2; ++k) {
    B2[k] = B1[k]; C2[k] = C1[k]; E2[k] = E1[k];
  }
#endif

#ifdef DOBLAS
  std::vector<tmv::Matrix<T> > B3 = B1;
  std::vector<tmv::Matrix<T> > C3 = C1;
  std::vector<tmv::Matrix<T> > E3 = E1;
#endif

#ifdef DOEIGEN
  std::vector<EIGENM,ALLOC(EIGENM) > B4(nloops2,EIGENM(K,M));
  std::vector<EIGENM,ALLOC(EIGENM) > C4(nloops2,EIGENM(K,N));
  std::vector<EIGENM,ALLOC(EIGENM) > E4(nloops2,EIGENM(M,N));
  for(int k=0;k<nloops2;++k) {
    for(int j=0;j<K;++j) for(int i=0;i<M;++i) B4[k](j,i) = B1[k](j,i);
    for(int j=0;j<N;++j) for(int i=0;i<K;++i) C4[k](i,j) = C1[k](i,j);
  }
#endif

#ifdef DOEIGENSMALL
  std::vector<EIGENSMB,ALLOC(EIGENSMB) > B5;
  std::vector<EIGENSMC,ALLOC(EIGENSMC) > C5;
  std::vector<EIGENSME,ALLOC(EIGENSME) > E5;
  if (nloops2x) { 
    B5.resize(nloops2x); C5.resize(nloops2x);
    E5.resize(nloops2x); 
  }
  for(int k=0;k<nloops2x;++k) {
    for(int j=0;j<K;++j) for(int i=0;i<M;++i) B5[k](j,i) = B1[k](j,i);
    for(int j=0;j<N;++j) for(int i=0;i<K;++i) C5[k](i,j) = C1[k](i,j);
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
#endif

  for (int n=0; n<nloops1; ++n) {

#if 1 // E = B.Transpose() * C
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = MPART1(B0[n2].Transpose()) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = B0[n2].Transpose().row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = MPART1(B1[k].Transpose()) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_reg = 0.;
      for (int k=0; k<nloops2; ++k) e1_reg += NormSq(E1[k]-E0[k]);
      e1_reg = sqrt(e1_reg/nloops2);
      e1_reg /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = MPART1(B2[k].Transpose()) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_small = 0.;
      for (int k=0; k<nloops2; ++k) e1_small += NormSq(E2[k]-E0[k]);
      e1_small = sqrt(e1_small/nloops2);
      e1_small /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLASNT,
	  M,N,K,one,BP(B3[k].cptr()),K,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,one,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,one,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_blas = 0.;
      for (int k=0; k<nloops2; ++k) e1_blas += NormSq(E3[k]-E0[k]);
      e1_blas = sqrt(e1_blas/nloops2);
      e1_blas /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = EPART1(B4[k].transpose()) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e1_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e1_eigen = sqrt(e1_eigen/nloops2);
      e1_eigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = EPART1(B5[k].transpose()) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e1_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e1_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e1_smalleigen = sqrt(e1_smalleigen/nloops2);
      e1_smalleigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif

#if 1 // E = -B.Transpose() * C
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = -MPART1(B0[n2].Transpose()) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = -B0[n2].Transpose().row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = -MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = -MPART1(B1[k].Transpose()) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_reg = 0.;
      for (int k=0; k<nloops2; ++k) e2_reg += NormSq(E1[k]-E0[k]);
      e2_reg = sqrt(e2_reg/nloops2);
      e2_reg /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = -MPART1(B2[k].Transpose()) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_small = 0.;
      for (int k=0; k<nloops2; ++k) e2_small += NormSq(E2[k]-E0[k]);
      e2_small = sqrt(e2_small/nloops2);
      e2_small /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLASNT,
	  M,N,K,mone,BP(B3[k].cptr()),K,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,mone,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,mone,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_blas = 0.;
      for (int k=0; k<nloops2; ++k) e2_blas += NormSq(E3[k]-E0[k]);
      e2_blas = sqrt(e2_blas/nloops2);
      e2_blas /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = -EPART1(B4[k].transpose()) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e2_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e2_eigen = sqrt(e2_eigen/nloops2);
      e2_eigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = -EPART1(B5[k].transpose()) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e2_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e2_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e2_smalleigen = sqrt(e2_smalleigen/nloops2);
      e2_smalleigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif

#if 1 // E = 7 * B.Transpose() * C
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = RT(7) * MPART1(B0[n2].Transpose()) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = RT(7) * B0[n2].Transpose().row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = RT(7) * MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = 7 * MPART1(B1[k].Transpose()) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_reg = 0.;
      for (int k=0; k<nloops2; ++k) e3_reg += NormSq(E1[k]-E0[k]);
      e3_reg = sqrt(e3_reg/nloops2);
      e3_reg /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = 7 * MPART1(B2[k].Transpose()) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_small = 0.;
      for (int k=0; k<nloops2; ++k) e3_small += NormSq(E2[k]-E0[k]);
      e3_small = sqrt(e3_small/nloops2);
      e3_small /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLASNT,
	  M,N,K,seven,BP(B3[k].cptr()),K,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,seven,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,seven,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_blas = 0.;
      for (int k=0; k<nloops2; ++k) e3_blas += NormSq(E3[k]-E0[k]);
      e3_blas = sqrt(e3_blas/nloops2);
      e3_blas /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = 7 * EPART1(B4[k].transpose()) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e3_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e3_eigen = sqrt(e3_eigen/nloops2);
      e3_eigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = 7 * EPART1(B5[k].transpose()) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e3_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e3_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e3_smalleigen = sqrt(e3_smalleigen/nloops2);
      e3_smalleigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif

#if 1 // E -= B.Transpose() * C
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) -= MPART1(B0[n2].Transpose()) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) -= B0[n2].Transpose().row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) -= MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) -= MPART1(B1[k].Transpose()) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_reg = 0.;
      for (int k=0; k<nloops2; ++k) e4_reg += NormSq(E1[k]-E0[k]);
      e4_reg = sqrt(e4_reg/nloops2);
      e4_reg /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) -= MPART1(B2[k].Transpose()) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_small = 0.;
      for (int k=0; k<nloops2; ++k) e4_small += NormSq(E2[k]-E0[k]);
      e4_small = sqrt(e4_small/nloops2);
      e4_small /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLASNT,
	  M,N,K,mone,BP(B3[k].cptr()),K,BP(C3[k].cptr()),K,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,one,BP(B3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) -= MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,one,BP(C3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) -= MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_blas = 0.;
      for (int k=0; k<nloops2; ++k) e4_blas += NormSq(E3[k]-E0[k]);
      e4_blas = sqrt(e4_blas/nloops2);
      e4_blas /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) -= EPART1(B4[k].transpose()) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e4_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e4_eigen = sqrt(e4_eigen/nloops2);
      e4_eigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) -= EPART1(B5[k].transpose()) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e4_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e4_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e4_smalleigen = sqrt(e4_smalleigen/nloops2);
      e4_smalleigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif

#if 1 // E += 8 * B.Transpose() * C
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += RT(8) * MPART1(B0[n2].Transpose()) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += RT(8) * B0[n2].Transpose().row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = RT(8) * MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += 8 * MPART1(B1[k].Transpose()) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_reg = 0.;
      for (int k=0; k<nloops2; ++k) e5_reg += NormSq(E1[k]-E0[k]);
      e5_reg = sqrt(e5_reg/nloops2);
      e5_reg /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += 8 * MPART1(B2[k].Transpose()) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_small = 0.;
      for (int k=0; k<nloops2; ++k) e5_small += NormSq(E2[k]-E0[k]);
      e5_small = sqrt(e5_small/nloops2);
      e5_small /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLASNT,
	  M,N,K,eight,BP(B3[k].cptr()),K,BP(C3[k].cptr()),K,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,eight,BP(B3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,eight,BP(C3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_blas = 0.;
      for (int k=0; k<nloops2; ++k) e5_blas += NormSq(E3[k]-E0[k]);
      e5_blas = sqrt(e5_blas/nloops2);
      e5_blas /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += 8 * EPART1(B4[k].transpose()) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e5_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e5_eigen = sqrt(e5_eigen/nloops2);
      e5_eigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += 8 * EPART1(B5[k].transpose()) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e5_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e5_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e5_smalleigen = sqrt(e5_smalleigen/nloops2);
      e5_smalleigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif

#if 1 // E = (7,1) * B.Transpose() * C
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(B0[n2].Transpose()) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * B0[n2].Transpose().row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(B1[k].Transpose()) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_reg = 0.;
      for (int k=0; k<nloops2; ++k) e6_reg += NormSq(E1[k]-E0[k]);
      e6_reg = sqrt(e6_reg/nloops2);
      e6_reg /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(B2[k].Transpose()) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_small = 0.;
      for (int k=0; k<nloops2; ++k) e6_small += NormSq(E2[k]-E0[k]);
      e6_small = sqrt(e6_small/nloops2);
      e6_small /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLASNT,
	  M,N,K,z71,BP(B3[k].cptr()),K,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,z71,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z71,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_blas = 0.;
      for (int k=0; k<nloops2; ++k) e6_blas += NormSq(E3[k]-E0[k]);
      e6_blas = sqrt(e6_blas/nloops2);
      e6_blas /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(B4[k].transpose()) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e6_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e6_eigen = sqrt(e6_eigen/nloops2);
      e6_eigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(B5[k].transpose()) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e6_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e6_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e6_smalleigen = sqrt(e6_smalleigen/nloops2);
      e6_smalleigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * B.Transpose() * C
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(B0[n2].Transpose()) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * B0[n2].Transpose().row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(B1[k].Transpose()) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_reg = 0.;
      for (int k=0; k<nloops2; ++k) e7_reg += NormSq(E1[k]-E0[k]);
      e7_reg = sqrt(e7_reg/nloops2);
      e7_reg /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(B2[k].Transpose()) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_small = 0.;
      for (int k=0; k<nloops2; ++k) e7_small += NormSq(E2[k]-E0[k]);
      e7_small = sqrt(e7_small/nloops2);
      e7_small /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLASNT,
	  M,N,K,z89,BP(B3[k].cptr()),K,BP(C3[k].cptr()),K,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,z89,BP(B3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z89,BP(C3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_blas = 0.;
      for (int k=0; k<nloops2; ++k) e7_blas += NormSq(E3[k]-E0[k]);
      e7_blas = sqrt(e7_blas/nloops2);
      e7_blas /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(B4[k].transpose()) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e7_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e7_eigen = sqrt(e7_eigen/nloops2);
      e7_eigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(B5[k].transpose()) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e7_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e7_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e7_smalleigen = sqrt(e7_smalleigen/nloops2);
      e7_smalleigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E = (7,1) * B.Adjoint() * C
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(B0[n2].Adjoint()) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * B0[n2].Adjoint().row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(B0[n2].Adjoint()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(B1[k].Adjoint()) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_reg = 0.;
      for (int k=0; k<nloops2; ++k) e8_reg += NormSq(E1[k]-E0[k]);
      e8_reg = sqrt(e8_reg/nloops2);
      e8_reg /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(B2[k].Adjoint()) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_small = 0.;
      for (int k=0; k<nloops2; ++k) e8_small += NormSq(E2[k]-E0[k]);
      e8_small = sqrt(e8_small/nloops2);
      e8_small /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASCT, BLASNT,
	  M,N,K,z71,BP(B3[k].cptr()),K,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLASCT, BLASDIAG1,
	  M,N,z71,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(B3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z71,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_blas = 0.;
      for (int k=0; k<nloops2; ++k) e8_blas += NormSq(E3[k]-E0[k]);
      e8_blas = sqrt(e8_blas/nloops2);
      e8_blas /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(B4[k].adjoint()) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e8_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e8_eigen = sqrt(e8_eigen/nloops2);
      e8_eigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(B5[k].adjoint()) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e8_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e8_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e8_smalleigen = sqrt(e8_smalleigen/nloops2);
      e8_smalleigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * B.Adjoint() * C
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(B0[n2].Adjoint()) * C0[n2].col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * B0[n2].Adjoint().row(i) * MPART2(C0[n2]);
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]);
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(B0[n2].Adjoint()) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(B1[k].Adjoint()) * MPART2(C1[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_reg = 0.;
      for (int k=0; k<nloops2; ++k) e9_reg += NormSq(E1[k]-E0[k]);
      e9_reg = sqrt(e9_reg/nloops2);
      e9_reg /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(B2[k].Adjoint()) * MPART2(C2[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_small = 0.;
      for (int k=0; k<nloops2; ++k) e9_small += NormSq(E2[k]-E0[k]);
      e9_small = sqrt(e9_small/nloops2);
      e9_small /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASCT, BLASNT,
	  M,N,K,z89,BP(B3[k].cptr()),K,BP(C3[k].cptr()),K,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(C3[k]);
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLASCT, BLASDIAG1,
	  M,N,z89,BP(B3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(B3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z89,BP(C3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_blas = 0.;
      for (int k=0; k<nloops2; ++k) e9_blas += NormSq(E3[k]-E0[k]);
      e9_blas = sqrt(e9_blas/nloops2);
      e9_blas /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(B4[k].adjoint()) * EPART2(C4[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e9_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e9_eigen = sqrt(e9_eigen/nloops2);
      e9_eigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(B5[k].adjoint()) * EPART2(C5[k]);

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e9_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e9_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e9_smalleigen = sqrt(e9_smalleigen/nloops2);
      e9_smalleigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E = (7,1) * B * C*
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(B0[n2].Transpose()) * C0[n2].col(j).Conjugate();
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * B0[n2].Transpose().row(i) * MPART2(C0[n2]).Conjugate();
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]).Conjugate();
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(B1[k].Transpose()) * MPART2(C1[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_reg = 0.;
      for (int k=0; k<nloops2; ++k) e10_reg += NormSq(E1[k]-E0[k]);
      e10_reg = sqrt(e10_reg/nloops2);
      e10_reg /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(B2[k].Transpose()) * MPART2(C2[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_small = 0.;
      for (int k=0; k<nloops2; ++k) e10_small += NormSq(E2[k]-E0[k]);
      e10_small = sqrt(e10_small/nloops2);
      e10_small /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
      BLASNAME(gemm) (BLASCM BLAST, BLASNT,
	  M,N,K,z71,BP(B3[k].cptr()),K,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(C3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,z71,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
      MPART3(E3[k]) = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z71,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_blas = 0.;
      for (int k=0; k<nloops2; ++k) e10_blas += NormSq(E3[k]-E0[k]);
      e10_blas = sqrt(e10_blas/nloops2);
      e10_blas /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(B4[k].transpose()) * EPART2(C4[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e10_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e10_eigen = sqrt(e10_eigen/nloops2);
      e10_eigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(B5[k].transpose()) * EPART2(C5[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e10_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e10_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e10_smalleigen = sqrt(e10_smalleigen/nloops2);
      e10_smalleigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * B * C*
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(B0[n2].Transpose()) * C0[n2].col(j).Conjugate();
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * B0[n2].Transpose().row(i) * MPART2(C0[n2]).Conjugate();
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]).Conjugate();
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(B1[k].Transpose()) * MPART2(C1[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_reg = 0.;
      for (int k=0; k<nloops2; ++k) e11_reg += NormSq(E1[k]-E0[k]);
      e11_reg = sqrt(e11_reg/nloops2);
      e11_reg /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(B2[k].Transpose()) * MPART2(C2[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_small = 0.;
      for (int k=0; k<nloops2; ++k) e11_small += NormSq(E2[k]-E0[k]);
      e11_small = sqrt(e11_small/nloops2);
      e11_small /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
      BLASNAME(gemm) (BLASCM BLAST, BLASNT,
	  M,N,K,z89,BP(B3[k].cptr()),K,BP(C3[k].cptr()),K,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(C3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,z89,BP(B3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(B3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z89c,BP(C3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
      MPART3(E3[k]) += MPART3(temp.Conjugate());
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_blas = 0.;
      for (int k=0; k<nloops2; ++k) e11_blas += NormSq(E3[k]-E0[k]);
      e11_blas = sqrt(e11_blas/nloops2);
      e11_blas /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(B4[k].transpose()) * EPART2(C4[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e11_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e11_eigen = sqrt(e11_eigen/nloops2);
      e11_eigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(B5[k].transpose()) * EPART2(C5[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e11_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e11_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e11_smalleigen = sqrt(e11_smalleigen/nloops2);
      e11_smalleigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E = (7,1) * B.Adjoint() * C*
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(B0[n2].Adjoint()) * C0[n2].col(j).Conjugate();
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * B0[n2].Adjoint().row(i) * MPART2(C0[n2]).Conjugate();
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]).Conjugate();
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(B0[n2].Adjoint()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(B1[k].Adjoint()) * MPART2(C1[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_reg = 0.;
      for (int k=0; k<nloops2; ++k) e12_reg += NormSq(E1[k]-E0[k]);
      e12_reg = sqrt(e12_reg/nloops2);
      e12_reg /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(B2[k].Adjoint()) * MPART2(C2[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_small = 0.;
      for (int k=0; k<nloops2; ++k) e12_small += NormSq(E2[k]-E0[k]);
      e12_small = sqrt(e12_small/nloops2);
      e12_small /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLASNT,
	  M,N,K,z71c,BP(B3[k].cptr()),K,BP(C3[k].cptr()),K,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (M*N,dmone,E3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(C3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLASCT, BLASDIAG1,
	  M,N,z71,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
      MPART3(E3[k]) = MPART1(B3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z71,BP(C3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_blas = 0.;
      for (int k=0; k<nloops2; ++k) e12_blas += NormSq(E3[k]-E0[k]);
      e12_blas = sqrt(e12_blas/nloops2);
      e12_blas /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(B4[k].adjoint()) * EPART2(C4[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e12_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e12_eigen = sqrt(e12_eigen/nloops2);
      e12_eigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(B5[k].adjoint()) * EPART2(C5[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e12_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e12_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e12_smalleigen = sqrt(e12_smalleigen/nloops2);
      e12_smalleigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * B.Adjoint() * C*
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(B0[n2].Adjoint()) * C0[n2].col(j).Conjugate();
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * B0[n2].Adjoint().row(i) * MPART2(C0[n2]).Conjugate();
#else
	tmv::Matrix<T> temp2 = MPART2(C0[n2]).Conjugate();
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(B0[n2].Adjoint()) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(B1[k].Adjoint()) * MPART2(C1[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_reg = 0.;
      for (int k=0; k<nloops2; ++k) e13_reg += NormSq(E1[k]-E0[k]);
      e13_reg = sqrt(e13_reg/nloops2);
      e13_reg /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(B2[k].Adjoint()) * MPART2(C2[k]).Conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_small = 0.;
      for (int k=0; k<nloops2; ++k) e13_small += NormSq(E2[k]-E0[k]);
      e13_small = sqrt(e13_small/nloops2);
      e13_small /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASDNAME(scal) (M*N,dmone,E3[k].Real().ptr()+1,2);
      BLASNAME(gemm) (BLASCM BLAST, BLASNT,
	  M,N,K,z89c,BP(B3[k].cptr()),K,BP(C3[k].cptr()),K,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (M*N,dmone,E3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(C3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLASCT, BLASDIAG1,
	  M,N,z89,BP(B3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2, BLASNT, BLASDIAG2,
	  M,N,z89c,BP(C3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      BLASDNAME(scal) (K*N,dmone,C3[k].Real().ptr()+1,2);
      MPART3(E3[k]) += MPART3(temp.Conjugate());
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_blas = 0.;
      for (int k=0; k<nloops2; ++k) e13_blas += NormSq(E3[k]-E0[k]);
      e13_blas = sqrt(e13_blas/nloops2);
      e13_blas /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(B4[k].adjoint()) * EPART2(C4[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e13_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e13_eigen = sqrt(e13_eigen/nloops2);
      e13_eigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(B5[k].adjoint()) * EPART2(C5[k]).conjugate();

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e13_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e13_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e13_smalleigen = sqrt(e13_smalleigen/nloops2);
      e13_smalleigen /= Norm(B0[0])*Norm(C0[0]);
    }
#endif
#endif
#endif
#endif
  }

  std::cout<<"E = B * C               "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
  std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
  std::cout<<"E = -B * C              "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
  std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
  std::cout<<"E = 7 * B * C           "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
  std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
  std::cout<<"E -= B * C              "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
  std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
  std::cout<<"E += 8 * B * C          "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
  std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
#ifdef TISCOMPLEX
  std::cout<<"E = (7,1) * B * C       "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
  std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B * C      "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
  std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * B* * C      "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
  std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B* * C     "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
  std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * B * C*      "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
  std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B * C*     "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
  std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * B* * C*     "<<t12_reg<<"  "<<t12_small<<"  "<<t12_blas;
  std::cout<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B* * C*    "<<t13_reg<<"  "<<t13_small<<"  "<<t13_blas;
  std::cout<<"  "<<t13_eigen<<"  "<<t13_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
  std::cout<<"errors:\n";
  std::cout<<"E = B * C               "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
  std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
  std::cout<<"E = -B * C              "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
  std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
  std::cout<<"E = 7 * B * C           "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
  std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
  std::cout<<"E -= B * C              "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
  std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
  std::cout<<"E += 8 * B * C          "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
  std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
#ifdef TISCOMPLEX
  std::cout<<"E = (7,1) * B * C       "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
  std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B * C      "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
  std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * B* * C      "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
  std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B* * C     "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
  std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * B * C*      "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
  std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B * C*     "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
  std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * B* * C*     "<<e12_reg<<"  "<<e12_small<<"  "<<e12_blas;
  std::cout<<"  "<<e12_eigen<<"  "<<e12_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B* * C*    "<<e13_reg<<"  "<<e13_small<<"  "<<e13_blas;
  std::cout<<"  "<<e13_eigen<<"  "<<e13_smalleigen<<std::endl;
#endif
  std::cout<<"\n\n";
#endif
}
#endif

#ifdef DOMULTMM_CRC
static void MultMM_CRC(
    const std::vector<tmv::Matrix<T> >& A1,
    const std::vector<tmv::Matrix<T> >& D1,
    std::vector<tmv::Matrix<T> >& E1)
{
#ifdef ERRORCHECK
  std::vector<tmv::Matrix<T> > A0 = A1;
  std::vector<tmv::Matrix<T> > D0 = D1;
  std::vector<tmv::Matrix<T> > E0 = E1;
#endif

#ifdef DOSMALL
  std::vector<tmv::SmallMatrix<T,M,K> > A2(nloops2);
  std::vector<tmv::SmallMatrix<T,N,K> > D2(nloops2);
  std::vector<tmv::SmallMatrix<T,M,N> > E2(nloops2);

  for(int k=0; k<nloops2; ++k) {
    A2[k] = A1[k]; D2[k] = D1[k]; E2[k] = E1[k];
  }
#endif

#ifdef DOBLAS
  std::vector<tmv::Matrix<T> > A3 = A1;
  std::vector<tmv::Matrix<T> > D3 = D1;
  std::vector<tmv::Matrix<T> > E3 = E1;
#endif

#ifdef DOEIGEN
  std::vector<EIGENM,ALLOC(EIGENM) > A4(nloops2,EIGENM(M,K));
  std::vector<EIGENM,ALLOC(EIGENM) > D4(nloops2,EIGENM(N,K));
  std::vector<EIGENM,ALLOC(EIGENM) > E4(nloops2,EIGENM(M,N));
  for(int k=0;k<nloops2;++k) {
    for(int j=0;j<K;++j) for(int i=0;i<M;++i) A4[k](i,j) = A1[k](i,j);
    for(int j=0;j<N;++j) for(int i=0;i<K;++i) D4[k](j,i) = D1[k](j,i);
  }
#endif

#ifdef DOEIGENSMALL
  std::vector<EIGENSMA,ALLOC(EIGENSMA) > A5;
  std::vector<EIGENSMD,ALLOC(EIGENSMD) > D5;
  std::vector<EIGENSME,ALLOC(EIGENSME) > E5;
  if (nloops2x) { 
    A5.resize(nloops2x); D5.resize(nloops2x);
    E5.resize(nloops2x); 
  }
  for(int k=0;k<nloops2x;++k) {
    for(int j=0;j<K;++j) for(int i=0;i<M;++i) A5[k](i,j) = A1[k](i,j);
    for(int j=0;j<N;++j) for(int i=0;i<K;++i) D5[k](j,i) = D1[k](j,i);
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
#endif

  for (int n=0; n<nloops1; ++n) {

#if 1 // E = A * D.Transpose()
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = MPART1(A0[n2]) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = A0[n2].row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = MPART1(A1[k]) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_reg = 0.;
      for (int k=0; k<nloops2; ++k) e1_reg += NormSq(E1[k]-E0[k]);
      e1_reg = sqrt(e1_reg/nloops2);
      e1_reg /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = MPART1(A2[k]) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_small = 0.;
      for (int k=0; k<nloops2; ++k) e1_small += NormSq(E2[k]-E0[k]);
      e1_small = sqrt(e1_small/nloops2);
      e1_small /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLAST,
	  M,N,K,one,BP(A3[k].cptr()),M,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,one,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,one,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_blas = 0.;
      for (int k=0; k<nloops2; ++k) e1_blas += NormSq(E3[k]-E0[k]);
      e1_blas = sqrt(e1_blas/nloops2);
      e1_blas /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = EPART1(A4[k]) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e1_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e1_eigen = sqrt(e1_eigen/nloops2);
      e1_eigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = EPART1(A5[k]) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e1_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e1_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e1_smalleigen = sqrt(e1_smalleigen/nloops2);
      e1_smalleigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif

#if 1 // E = -A * D.Transpose()
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = -MPART1(A0[n2]) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = -A0[n2].row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = -MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = -MPART1(A1[k]) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_reg = 0.;
      for (int k=0; k<nloops2; ++k) e2_reg += NormSq(E1[k]-E0[k]);
      e2_reg = sqrt(e2_reg/nloops2);
      e2_reg /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = -MPART1(A2[k]) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_small = 0.;
      for (int k=0; k<nloops2; ++k) e2_small += NormSq(E2[k]-E0[k]);
      e2_small = sqrt(e2_small/nloops2);
      e2_small /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLAST,
	  M,N,K,mone,BP(A3[k].cptr()),M,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,mone,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,mone,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_blas = 0.;
      for (int k=0; k<nloops2; ++k) e2_blas += NormSq(E3[k]-E0[k]);
      e2_blas = sqrt(e2_blas/nloops2);
      e2_blas /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = -EPART1(A4[k]) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e2_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e2_eigen = sqrt(e2_eigen/nloops2);
      e2_eigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = -EPART1(A5[k]) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e2_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e2_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e2_smalleigen = sqrt(e2_smalleigen/nloops2);
      e2_smalleigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif

#if 1 // E = 7 * A * D.Transpose()
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = RT(7) * MPART1(A0[n2]) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = RT(7) * A0[n2].row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = RT(7) * MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = 7 * MPART1(A1[k]) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_reg = 0.;
      for (int k=0; k<nloops2; ++k) e3_reg += NormSq(E1[k]-E0[k]);
      e3_reg = sqrt(e3_reg/nloops2);
      e3_reg /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = 7 * MPART1(A2[k]) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_small = 0.;
      for (int k=0; k<nloops2; ++k) e3_small += NormSq(E2[k]-E0[k]);
      e3_small = sqrt(e3_small/nloops2);
      e3_small /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLAST,
	  M,N,K,seven,BP(A3[k].cptr()),M,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,seven,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,seven,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_blas = 0.;
      for (int k=0; k<nloops2; ++k) e3_blas += NormSq(E3[k]-E0[k]);
      e3_blas = sqrt(e3_blas/nloops2);
      e3_blas /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = 7 * EPART1(A4[k]) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e3_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e3_eigen = sqrt(e3_eigen/nloops2);
      e3_eigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = 7 * EPART1(A5[k]) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e3_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e3_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e3_smalleigen = sqrt(e3_smalleigen/nloops2);
      e3_smalleigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif

#if 1 // E -= A * D.Transpose()
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) -= MPART1(A0[n2]) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) -= A0[n2].row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) -= MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) -= MPART1(A1[k]) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_reg = 0.;
      for (int k=0; k<nloops2; ++k) e4_reg += NormSq(E1[k]-E0[k]);
      e4_reg = sqrt(e4_reg/nloops2);
      e4_reg /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) -= MPART1(A2[k]) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_small = 0.;
      for (int k=0; k<nloops2; ++k) e4_small += NormSq(E2[k]-E0[k]);
      e4_small = sqrt(e4_small/nloops2);
      e4_small /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLAST,
	  M,N,K,mone,BP(A3[k].cptr()),M,BP(D3[k].cptr()),N,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,one,BP(A3[k].cptr()),M,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) -= MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,one,BP(D3[k].cptr()),N,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) -= MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_blas = 0.;
      for (int k=0; k<nloops2; ++k) e4_blas += NormSq(E3[k]-E0[k]);
      e4_blas = sqrt(e4_blas/nloops2);
      e4_blas /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) -= EPART1(A4[k]) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e4_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e4_eigen = sqrt(e4_eigen/nloops2);
      e4_eigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) -= EPART1(A5[k]) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e4_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e4_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e4_smalleigen = sqrt(e4_smalleigen/nloops2);
      e4_smalleigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif

#if 1 // E += 8 * A * D.Transpose()
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += RT(8) * MPART1(A0[n2]) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += RT(8) * A0[n2].row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = RT(8) * MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += 8 * MPART1(A1[k]) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_reg = 0.;
      for (int k=0; k<nloops2; ++k) e5_reg += NormSq(E1[k]-E0[k]);
      e5_reg = sqrt(e5_reg/nloops2);
      e5_reg /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += 8 * MPART1(A2[k]) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_small = 0.;
      for (int k=0; k<nloops2; ++k) e5_small += NormSq(E2[k]-E0[k]);
      e5_small = sqrt(e5_small/nloops2);
      e5_small /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLAST,
	  M,N,K,eight,BP(A3[k].cptr()),M,BP(D3[k].cptr()),N,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,eight,BP(A3[k].cptr()),M,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,eight,BP(D3[k].cptr()),N,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_blas = 0.;
      for (int k=0; k<nloops2; ++k) e5_blas += NormSq(E3[k]-E0[k]);
      e5_blas = sqrt(e5_blas/nloops2);
      e5_blas /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += 8 * EPART1(A4[k]) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e5_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e5_eigen = sqrt(e5_eigen/nloops2);
      e5_eigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += 8 * EPART1(A5[k]) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e5_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e5_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e5_smalleigen = sqrt(e5_smalleigen/nloops2);
      e5_smalleigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif

#if 1 // E = (7,1) * A * D.Transpose()
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(A0[n2]) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * A0[n2].row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(A1[k]) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_reg = 0.;
      for (int k=0; k<nloops2; ++k) e6_reg += NormSq(E1[k]-E0[k]);
      e6_reg = sqrt(e6_reg/nloops2);
      e6_reg /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(A2[k]) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_small = 0.;
      for (int k=0; k<nloops2; ++k) e6_small += NormSq(E2[k]-E0[k]);
      e6_small = sqrt(e6_small/nloops2);
      e6_small /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLAST,
	  M,N,K,z71,BP(A3[k].cptr()),M,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z71,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,z71,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_blas = 0.;
      for (int k=0; k<nloops2; ++k) e6_blas += NormSq(E3[k]-E0[k]);
      e6_blas = sqrt(e6_blas/nloops2);
      e6_blas /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(A4[k]) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e6_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e6_eigen = sqrt(e6_eigen/nloops2);
      e6_eigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(A5[k]) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e6_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e6_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e6_smalleigen = sqrt(e6_smalleigen/nloops2);
      e6_smalleigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * A * D.Transpose()
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(A0[n2]) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * A0[n2].row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(A1[k]) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_reg = 0.;
      for (int k=0; k<nloops2; ++k) e7_reg += NormSq(E1[k]-E0[k]);
      e7_reg = sqrt(e7_reg/nloops2);
      e7_reg /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(A2[k]) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_small = 0.;
      for (int k=0; k<nloops2; ++k) e7_small += NormSq(E2[k]-E0[k]);
      e7_small = sqrt(e7_small/nloops2);
      e7_small /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLAST,
	  M,N,K,z89,BP(A3[k].cptr()),M,BP(D3[k].cptr()),N,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z89,BP(A3[k].cptr()),M,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,z89,BP(D3[k].cptr()),N,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_blas = 0.;
      for (int k=0; k<nloops2; ++k) e7_blas += NormSq(E3[k]-E0[k]);
      e7_blas = sqrt(e7_blas/nloops2);
      e7_blas /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(A4[k]) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e7_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e7_eigen = sqrt(e7_eigen/nloops2);
      e7_eigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(A5[k]) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e7_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e7_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e7_smalleigen = sqrt(e7_smalleigen/nloops2);
      e7_smalleigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E = (7,1) * A* * D.Transpose()
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(A0[n2]).Conjugate() * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * A0[n2].row(i).Conjugate() * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(A0[n2]).Conjugate() * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(A1[k]).Conjugate() * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_reg = 0.;
      for (int k=0; k<nloops2; ++k) e8_reg += NormSq(E1[k]-E0[k]);
      e8_reg = sqrt(e8_reg/nloops2);
      e8_reg /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(A2[k]).Conjugate() * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_small = 0.;
      for (int k=0; k<nloops2; ++k) e8_small += NormSq(E2[k]-E0[k]);
      e8_small = sqrt(e8_small/nloops2);
      e8_small /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASDNAME(scal) (M*K,dmone,A3[k].Real().ptr()+1,2);
      BLASNAME(gemm) (BLASCM BLASNT, BLAST,
	  M,N,K,z71,BP(A3[k].cptr()),M,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (M*K,dmone,A3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z71c,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      BLASDNAME(scal) (M*K,dmone,E3[k].Real().ptr()+1,2);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(A3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,z71,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_blas = 0.;
      for (int k=0; k<nloops2; ++k) e8_blas += NormSq(E3[k]-E0[k]);
      e8_blas = sqrt(e8_blas/nloops2);
      e8_blas /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(A4[k]).conjugate() * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e8_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e8_eigen = sqrt(e8_eigen/nloops2);
      e8_eigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(A5[k]).conjugate() * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e8_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e8_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e8_smalleigen = sqrt(e8_smalleigen/nloops2);
      e8_smalleigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * A* * D.Transpose()
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(A0[n2]).Conjugate() * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * A0[n2].row(i).Conjugate() * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(A0[n2]).Conjugate() * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(A1[k]).Conjugate() * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_reg = 0.;
      for (int k=0; k<nloops2; ++k) e9_reg += NormSq(E1[k]-E0[k]);
      e9_reg = sqrt(e9_reg/nloops2);
      e9_reg /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(A2[k]).Conjugate() * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_small = 0.;
      for (int k=0; k<nloops2; ++k) e9_small += NormSq(E2[k]-E0[k]);
      e9_small = sqrt(e9_small/nloops2);
      e9_small /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASDNAME(scal) (M*K,dmone,A3[k].Real().ptr()+1,2);
      BLASNAME(gemm) (BLASCM BLASNT, BLAST,
	  M,N,K,z89,BP(A3[k].cptr()),M,BP(D3[k].cptr()),N,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (M*K,dmone,A3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(D3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z89c,BP(A3[k].cptr()),M,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp.Conjugate());
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(A3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,z89,BP(D3[k].cptr()),N,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_blas = 0.;
      for (int k=0; k<nloops2; ++k) e9_blas += NormSq(E3[k]-E0[k]);
      e9_blas = sqrt(e9_blas/nloops2);
      e9_blas /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(A4[k]).conjugate() * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e9_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e9_eigen = sqrt(e9_eigen/nloops2);
      e9_eigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(A5[k]).conjugate() * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e9_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e9_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e9_smalleigen = sqrt(e9_smalleigen/nloops2);
      e9_smalleigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E = (7,1) * A * D.Adjoint()
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(A0[n2]) * D0[n2].Adjoint().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * A0[n2].row(i) * MPART2(D0[n2].Adjoint());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Adjoint());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(A1[k]) * MPART2(D1[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_reg = 0.;
      for (int k=0; k<nloops2; ++k) e10_reg += NormSq(E1[k]-E0[k]);
      e10_reg = sqrt(e10_reg/nloops2);
      e10_reg /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(A2[k]) * MPART2(D2[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_small = 0.;
      for (int k=0; k<nloops2; ++k) e10_small += NormSq(E2[k]-E0[k]);
      e10_small = sqrt(e10_small/nloops2);
      e10_small /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLASCT,
	  M,N,K,z71,BP(A3[k].cptr()),M,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z71,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLASCT, BLASDIAG2,
	  M,N,z71,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_blas = 0.;
      for (int k=0; k<nloops2; ++k) e10_blas += NormSq(E3[k]-E0[k]);
      e10_blas = sqrt(e10_blas/nloops2);
      e10_blas /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(A4[k]) * EPART2(D4[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e10_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e10_eigen = sqrt(e10_eigen/nloops2);
      e10_eigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(A5[k]) * EPART2(D5[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e10_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e10_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e10_smalleigen = sqrt(e10_smalleigen/nloops2);
      e10_smalleigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * A * D.Adjoint()
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(A0[n2]) * D0[n2].Adjoint().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * A0[n2].row(i) * MPART2(D0[n2].Adjoint());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Adjoint());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(A0[n2]) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(A1[k]) * MPART2(D1[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_reg = 0.;
      for (int k=0; k<nloops2; ++k) e11_reg += NormSq(E1[k]-E0[k]);
      e11_reg = sqrt(e11_reg/nloops2);
      e11_reg /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(A2[k]) * MPART2(D2[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_small = 0.;
      for (int k=0; k<nloops2; ++k) e11_small += NormSq(E2[k]-E0[k]);
      e11_small = sqrt(e11_small/nloops2);
      e11_small /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLASCT,
	  M,N,K,z89,BP(A3[k].cptr()),M,BP(D3[k].cptr()),N,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(D3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z89,BP(A3[k].cptr()),M,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(A3[k]);
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLASCT, BLASDIAG2,
	  M,N,z89,BP(D3[k].cptr()),N,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_blas = 0.;
      for (int k=0; k<nloops2; ++k) e11_blas += NormSq(E3[k]-E0[k]);
      e11_blas = sqrt(e11_blas/nloops2);
      e11_blas /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(A4[k]) * EPART2(D4[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e11_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e11_eigen = sqrt(e11_eigen/nloops2);
      e11_eigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(A5[k]) * EPART2(D5[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e11_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e11_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e11_smalleigen = sqrt(e11_smalleigen/nloops2);
      e11_smalleigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E = (7,1) * A* * D.Adjoint()
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(A0[n2]).Conjugate() * D0[n2].Adjoint().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * A0[n2].row(i).Conjugate() * MPART2(D0[n2].Adjoint());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Adjoint());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(A0[n2]).Conjugate() * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(A1[k]).Conjugate() * MPART2(D1[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_reg = 0.;
      for (int k=0; k<nloops2; ++k) e12_reg += NormSq(E1[k]-E0[k]);
      e12_reg = sqrt(e12_reg/nloops2);
      e12_reg /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(A2[k]).Conjugate() * MPART2(D2[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_small = 0.;
      for (int k=0; k<nloops2; ++k) e12_small += NormSq(E2[k]-E0[k]);
      e12_small = sqrt(e12_small/nloops2);
      e12_small /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASNT, BLAST,
	  M,N,K,z71c,BP(A3[k].cptr()),M,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (M*N,dmone,E3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z71c,BP(A3[k].cptr()),M,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      BLASDNAME(scal) (M*N,dmone,E3[k].Real().ptr()+1,2);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(A3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLASCT, BLASDIAG2,
	  M,N,z71,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_blas = 0.;
      for (int k=0; k<nloops2; ++k) e12_blas += NormSq(E3[k]-E0[k]);
      e12_blas = sqrt(e12_blas/nloops2);
      e12_blas /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(A4[k]).conjugate() * EPART2(D4[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e12_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e12_eigen = sqrt(e12_eigen/nloops2);
      e12_eigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(A5[k]).conjugate() * EPART2(D5[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e12_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e12_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e12_smalleigen = sqrt(e12_smalleigen/nloops2);
      e12_smalleigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * A* * D.Adjoint()
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(A0[n2]).Conjugate() * D0[n2].Adjoint().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * A0[n2].row(i).Conjugate() * MPART2(D0[n2].Adjoint());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Adjoint());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(A0[n2]).Conjugate() * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(A1[k]).Conjugate() * MPART2(D1[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_reg = 0.;
      for (int k=0; k<nloops2; ++k) e13_reg += NormSq(E1[k]-E0[k]);
      e13_reg = sqrt(e13_reg/nloops2);
      e13_reg /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(A2[k]).Conjugate() * MPART2(D2[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_small = 0.;
      for (int k=0; k<nloops2; ++k) e13_small += NormSq(E2[k]-E0[k]);
      e13_small = sqrt(e13_small/nloops2);
      e13_small /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASDNAME(scal) (M*N,dmone,E3[k].Real().ptr()+1,2);
      BLASNAME(gemm) (BLASCM BLASNT, BLAST,
	  M,N,K,z89c,BP(A3[k].cptr()),M,BP(D3[k].cptr()),N,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
      BLASDNAME(scal) (M*N,dmone,E3[k].Real().ptr()+1,2);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1, BLASNT, BLASDIAG1,
	  M,N,z89c,BP(A3[k].cptr()),M,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp.Conjugate());
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(A3[k].Conjugate());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLASCT, BLASDIAG2,
	  M,N,z89,BP(D3[k].cptr()),N,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_blas = 0.;
      for (int k=0; k<nloops2; ++k) e13_blas += NormSq(E3[k]-E0[k]);
      e13_blas = sqrt(e13_blas/nloops2);
      e13_blas /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(A4[k]).conjugate() * EPART2(D4[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e13_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e13_eigen = sqrt(e13_eigen/nloops2);
      e13_eigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(A5[k]).conjugate() * EPART2(D5[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e13_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e13_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e13_smalleigen = sqrt(e13_smalleigen/nloops2);
      e13_smalleigen /= Norm(A0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif
  }

  std::cout<<"E = A * D               "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
  std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
  std::cout<<"E = -A * D              "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
  std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
  std::cout<<"E = 7 * A * D           "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
  std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
  std::cout<<"E -= A * D              "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
  std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
  std::cout<<"E += 8 * A * D          "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
  std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
#ifdef TISCOMPLEX
  std::cout<<"E = (7,1) * A * D       "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
  std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A * D      "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
  std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * A* * D      "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
  std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A* * D     "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
  std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * A * D*      "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
  std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A * D*     "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
  std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * A* * D*     "<<t12_reg<<"  "<<t12_small<<"  "<<t12_blas;
  std::cout<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A* * D*    "<<t13_reg<<"  "<<t13_small<<"  "<<t13_blas;
  std::cout<<"  "<<t13_eigen<<"  "<<t13_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
  std::cout<<"errors:\n";
  std::cout<<"E = A * D               "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
  std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
  std::cout<<"E = -A * D              "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
  std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
  std::cout<<"E = 7 * A * D           "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
  std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
  std::cout<<"E -= A * D              "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
  std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
  std::cout<<"E += 8 * A * D          "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
  std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
#ifdef TISCOMPLEX
  std::cout<<"E = (7,1) * A * D       "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
  std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A * D      "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
  std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * A* * D      "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
  std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A* * D     "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
  std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * A * D*      "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
  std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A * D*     "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
  std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * A* * D*     "<<e12_reg<<"  "<<e12_small<<"  "<<e12_blas;
  std::cout<<"  "<<e12_eigen<<"  "<<e12_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * A* * D*    "<<e13_reg<<"  "<<e13_small<<"  "<<e13_blas;
  std::cout<<"  "<<e13_eigen<<"  "<<e13_smalleigen<<std::endl;
#endif
  std::cout<<"\n\n";
#endif
}
#endif


#ifdef DOMULTMM_RRC
static void MultMM_RRC(
    const std::vector<tmv::Matrix<T> >& B1,
    const std::vector<tmv::Matrix<T> >& D1,
    std::vector<tmv::Matrix<T> >& E1)
{
#ifdef ERRORCHECK
  std::vector<tmv::Matrix<T> > B0 = B1;
  std::vector<tmv::Matrix<T> > D0 = D1;
  std::vector<tmv::Matrix<T> > E0 = E1;
#endif

#ifdef DOSMALL
  std::vector<tmv::SmallMatrix<T,K,M> > B2(nloops2);
  std::vector<tmv::SmallMatrix<T,N,K> > D2(nloops2);
  std::vector<tmv::SmallMatrix<T,M,N> > E2(nloops2);

  for(int k=0; k<nloops2; ++k) {
    B2[k] = B1[k]; D2[k] = D1[k]; E2[k] = E1[k];
  }
#endif

#ifdef DOBLAS
  std::vector<tmv::Matrix<T> > B3 = B1;
  std::vector<tmv::Matrix<T> > D3 = D1;
  std::vector<tmv::Matrix<T> > E3 = E1;
#endif

#ifdef DOEIGEN
  std::vector<EIGENM,ALLOC(EIGENM) > B4(nloops2,EIGENM(K,M));
  std::vector<EIGENM,ALLOC(EIGENM) > D4(nloops2,EIGENM(N,K));
  std::vector<EIGENM,ALLOC(EIGENM) > E4(nloops2,EIGENM(M,N));
  for(int k=0;k<nloops2;++k) {
    for(int j=0;j<K;++j) for(int i=0;i<M;++i) B4[k](j,i) = B1[k](j,i);
    for(int j=0;j<N;++j) for(int i=0;i<K;++i) D4[k](j,i) = D1[k](j,i);
  }
#endif

#ifdef DOEIGENSMALL
  std::vector<EIGENSMB,ALLOC(EIGENSMB) > B5;
  std::vector<EIGENSMD,ALLOC(EIGENSMD) > D5;
  std::vector<EIGENSME,ALLOC(EIGENSME) > E5;
  if (nloops2x) { 
    B5.resize(nloops2x); D5.resize(nloops2x);
    E5.resize(nloops2x); 
  }
  for(int k=0;k<nloops2x;++k) {
    for(int j=0;j<K;++j) for(int i=0;i<M;++i) B5[k](j,i) = B1[k](j,i);
    for(int j=0;j<N;++j) for(int i=0;i<K;++i) D5[k](j,i) = D1[k](j,i);
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
#endif

  for (int n=0; n<nloops1; ++n) {

#if 1 // E = B * D
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = MPART1(B0[n2].Transpose()) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = B0[n2].Transpose().row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = MPART1(B1[k].Transpose()) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_reg = 0.;
      for (int k=0; k<nloops2; ++k) e1_reg += NormSq(E1[k]-E0[k]);
      e1_reg = sqrt(e1_reg/nloops2);
      e1_reg /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = MPART1(B2[k].Transpose()) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_small = 0.;
      for (int k=0; k<nloops2; ++k) e1_small += NormSq(E2[k]-E0[k]);
      e1_small = sqrt(e1_small/nloops2);
      e1_small /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLAST,
	  M,N,K,one,BP(B3[k].cptr()),K,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,one,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,one,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_blas = 0.;
      for (int k=0; k<nloops2; ++k) e1_blas += NormSq(E3[k]-E0[k]);
      e1_blas = sqrt(e1_blas/nloops2);
      e1_blas /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = EPART1(B4[k].transpose()) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e1_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e1_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e1_eigen = sqrt(e1_eigen/nloops2);
      e1_eigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = EPART1(B5[k].transpose()) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e1_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e1_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e1_smalleigen = sqrt(e1_smalleigen/nloops2);
      e1_smalleigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif

#if 1 // E = -B * D
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = -MPART1(B0[n2].Transpose()) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = -B0[n2].Transpose().row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = -MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = -MPART1(B1[k].Transpose()) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_reg = 0.;
      for (int k=0; k<nloops2; ++k) e2_reg += NormSq(E1[k]-E0[k]);
      e2_reg = sqrt(e2_reg/nloops2);
      e2_reg /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = -MPART1(B2[k].Transpose()) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_small = 0.;
      for (int k=0; k<nloops2; ++k) e2_small += NormSq(E2[k]-E0[k]);
      e2_small = sqrt(e2_small/nloops2);
      e2_small /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLAST,
	  M,N,K,mone,BP(B3[k].cptr()),K,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,mone,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,mone,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_blas = 0.;
      for (int k=0; k<nloops2; ++k) e2_blas += NormSq(E3[k]-E0[k]);
      e2_blas = sqrt(e2_blas/nloops2);
      e2_blas /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = -EPART1(B4[k].transpose()) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e2_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e2_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e2_eigen = sqrt(e2_eigen/nloops2);
      e2_eigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = -EPART1(B5[k].transpose()) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e2_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e2_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e2_smalleigen = sqrt(e2_smalleigen/nloops2);
      e2_smalleigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif

#if 1 // E = 7 * B * D
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = RT(7) * MPART1(B0[n2].Transpose()) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = RT(7) * B0[n2].Transpose().row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = RT(7) * MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = 7 * MPART1(B1[k].Transpose()) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_reg = 0.;
      for (int k=0; k<nloops2; ++k) e3_reg += NormSq(E1[k]-E0[k]);
      e3_reg = sqrt(e3_reg/nloops2);
      e3_reg /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = 7 * MPART1(B2[k].Transpose()) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_small = 0.;
      for (int k=0; k<nloops2; ++k) e3_small += NormSq(E2[k]-E0[k]);
      e3_small = sqrt(e3_small/nloops2);
      e3_small /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLAST,
	  M,N,K,seven,BP(B3[k].cptr()),K,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,seven,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,seven,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_blas = 0.;
      for (int k=0; k<nloops2; ++k) e3_blas += NormSq(E3[k]-E0[k]);
      e3_blas = sqrt(e3_blas/nloops2);
      e3_blas /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = 7 * EPART1(B4[k].transpose()) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e3_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e3_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e3_eigen = sqrt(e3_eigen/nloops2);
      e3_eigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = 7 * EPART1(B5[k].transpose()) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e3_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e3_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e3_smalleigen = sqrt(e3_smalleigen/nloops2);
      e3_smalleigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif

#if 1 // E -= B * D
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) -= MPART1(B0[n2].Transpose()) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) -= B0[n2].Transpose().row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) -= MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) -= MPART1(B1[k].Transpose()) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_reg = 0.;
      for (int k=0; k<nloops2; ++k) e4_reg += NormSq(E1[k]-E0[k]);
      e4_reg = sqrt(e4_reg/nloops2);
      e4_reg /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) -= MPART1(B2[k].Transpose()) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_small = 0.;
      for (int k=0; k<nloops2; ++k) e4_small += NormSq(E2[k]-E0[k]);
      e4_small = sqrt(e4_small/nloops2);
      e4_small /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLAST,
	  M,N,K,mone,BP(B3[k].cptr()),K,BP(D3[k].cptr()),N,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,one,BP(B3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) -= MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,one,BP(D3[k].cptr()),N,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) -= MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_blas = 0.;
      for (int k=0; k<nloops2; ++k) e4_blas += NormSq(E3[k]-E0[k]);
      e4_blas = sqrt(e4_blas/nloops2);
      e4_blas /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) -= EPART1(B4[k].transpose()) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e4_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e4_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e4_eigen = sqrt(e4_eigen/nloops2);
      e4_eigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) -= EPART1(B5[k].transpose()) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e4_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e4_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e4_smalleigen = sqrt(e4_smalleigen/nloops2);
      e4_smalleigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif

#if 1 // E += 8 * B * D
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += RT(8) * MPART1(B0[n2].Transpose()) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += RT(8) * B0[n2].Transpose().row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = RT(8) * MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += 8 * MPART1(B1[k].Transpose()) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_reg = 0.;
      for (int k=0; k<nloops2; ++k) e5_reg += NormSq(E1[k]-E0[k]);
      e5_reg = sqrt(e5_reg/nloops2);
      e5_reg /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += 8 * MPART1(B2[k].Transpose()) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_small = 0.;
      for (int k=0; k<nloops2; ++k) e5_small += NormSq(E2[k]-E0[k]);
      e5_small = sqrt(e5_small/nloops2);
      e5_small /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLAST,
	  M,N,K,eight,BP(B3[k].cptr()),K,BP(D3[k].cptr()),N,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,eight,BP(B3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,eight,BP(D3[k].cptr()),N,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_blas = 0.;
      for (int k=0; k<nloops2; ++k) e5_blas += NormSq(E3[k]-E0[k]);
      e5_blas = sqrt(e5_blas/nloops2);
      e5_blas /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += 8 * EPART1(B4[k].transpose()) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e5_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e5_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e5_eigen = sqrt(e5_eigen/nloops2);
      e5_eigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += 8 * EPART1(B5[k].transpose()) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e5_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e5_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e5_smalleigen = sqrt(e5_smalleigen/nloops2);
      e5_smalleigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif

#if 1 // E = (7,1) * B * D
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(B0[n2].Transpose()) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * B0[n2].Transpose().row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(B1[k].Transpose()) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_reg = 0.;
      for (int k=0; k<nloops2; ++k) e6_reg += NormSq(E1[k]-E0[k]);
      e6_reg = sqrt(e6_reg/nloops2);
      e6_reg /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(B2[k].Transpose()) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_small = 0.;
      for (int k=0; k<nloops2; ++k) e6_small += NormSq(E2[k]-E0[k]);
      e6_small = sqrt(e6_small/nloops2);
      e6_small /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLAST,
	  M,N,K,z71,BP(B3[k].cptr()),K,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,z71,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,z71,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_blas = 0.;
      for (int k=0; k<nloops2; ++k) e6_blas += NormSq(E3[k]-E0[k]);
      e6_blas = sqrt(e6_blas/nloops2);
      e6_blas /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(B4[k].transpose()) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e6_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e6_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e6_eigen = sqrt(e6_eigen/nloops2);
      e6_eigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(B5[k].transpose()) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e6_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e6_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e6_smalleigen = sqrt(e6_smalleigen/nloops2);
      e6_smalleigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * B * D
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(B0[n2].Transpose()) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * B0[n2].Transpose().row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(B1[k].Transpose()) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_reg = 0.;
      for (int k=0; k<nloops2; ++k) e7_reg += NormSq(E1[k]-E0[k]);
      e7_reg = sqrt(e7_reg/nloops2);
      e7_reg /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(B2[k].Transpose()) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_small = 0.;
      for (int k=0; k<nloops2; ++k) e7_small += NormSq(E2[k]-E0[k]);
      e7_small = sqrt(e7_small/nloops2);
      e7_small /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLAST,
	  M,N,K,z89,BP(B3[k].cptr()),K,BP(D3[k].cptr()),N,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,z89,BP(B3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,z89,BP(D3[k].cptr()),N,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_blas = 0.;
      for (int k=0; k<nloops2; ++k) e7_blas += NormSq(E3[k]-E0[k]);
      e7_blas = sqrt(e7_blas/nloops2);
      e7_blas /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(B4[k].transpose()) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e7_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e7_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e7_eigen = sqrt(e7_eigen/nloops2);
      e7_eigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(B5[k].transpose()) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e7_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e7_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e7_smalleigen = sqrt(e7_smalleigen/nloops2);
      e7_smalleigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E = (7,1) * B.Adjoint() * D
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(B0[n2].Adjoint()) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * B0[n2].Adjoint().row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(B0[n2].Adjoint()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(B1[k].Adjoint()) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_reg = 0.;
      for (int k=0; k<nloops2; ++k) e8_reg += NormSq(E1[k]-E0[k]);
      e8_reg = sqrt(e8_reg/nloops2);
      e8_reg /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(B2[k].Adjoint()) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_small = 0.;
      for (int k=0; k<nloops2; ++k) e8_small += NormSq(E2[k]-E0[k]);
      e8_small = sqrt(e8_small/nloops2);
      e8_small /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASCT, BLAST,
	  M,N,K,z71,BP(B3[k].cptr()),K,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLASCT, BLASDIAG1,
	  M,N,z71,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(B3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,z71,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_blas = 0.;
      for (int k=0; k<nloops2; ++k) e8_blas += NormSq(E3[k]-E0[k]);
      e8_blas = sqrt(e8_blas/nloops2);
      e8_blas /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(B4[k].adjoint()) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e8_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e8_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e8_eigen = sqrt(e8_eigen/nloops2);
      e8_eigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(B5[k].adjoint()) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e8_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e8_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e8_smalleigen = sqrt(e8_smalleigen/nloops2);
      e8_smalleigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * B.Adjoint() * D
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(B0[n2].Adjoint()) * D0[n2].Transpose().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * B0[n2].Adjoint().row(i) * MPART2(D0[n2].Transpose());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Transpose());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(B0[n2].Adjoint()) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(B1[k].Adjoint()) * MPART2(D1[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_reg = 0.;
      for (int k=0; k<nloops2; ++k) e9_reg += NormSq(E1[k]-E0[k]);
      e9_reg = sqrt(e9_reg/nloops2);
      e9_reg /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(B2[k].Adjoint()) * MPART2(D2[k].Transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_small = 0.;
      for (int k=0; k<nloops2; ++k) e9_small += NormSq(E2[k]-E0[k]);
      e9_small = sqrt(e9_small/nloops2);
      e9_small /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASCT, BLAST,
	  M,N,K,z89,BP(B3[k].cptr()),K,BP(D3[k].cptr()),N,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(D3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLASCT, BLASDIAG1,
	  M,N,z89,BP(B3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(B3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLAST, BLASDIAG2,
	  M,N,z89,BP(D3[k].cptr()),N,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_blas = 0.;
      for (int k=0; k<nloops2; ++k) e9_blas += NormSq(E3[k]-E0[k]);
      e9_blas = sqrt(e9_blas/nloops2);
      e9_blas /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(B4[k].adjoint()) * EPART2(D4[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e9_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e9_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e9_eigen = sqrt(e9_eigen/nloops2);
      e9_eigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(B5[k].adjoint()) * EPART2(D5[k].transpose());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e9_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e9_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e9_smalleigen = sqrt(e9_smalleigen/nloops2);
      e9_smalleigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E = (7,1) * B * D.Adjoint()
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(B0[n2].Transpose()) * D0[n2].Adjoint().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * B0[n2].Transpose().row(i) * MPART2(D0[n2].Adjoint());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Adjoint());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(B1[k].Transpose()) * MPART2(D1[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_reg = 0.;
      for (int k=0; k<nloops2; ++k) e10_reg += NormSq(E1[k]-E0[k]);
      e10_reg = sqrt(e10_reg/nloops2);
      e10_reg /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(B2[k].Transpose()) * MPART2(D2[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_small = 0.;
      for (int k=0; k<nloops2; ++k) e10_small += NormSq(E2[k]-E0[k]);
      e10_small = sqrt(e10_small/nloops2);
      e10_small /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLASCT,
	  M,N,K,z71,BP(B3[k].cptr()),K,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,z71,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLASCT, BLASDIAG2,
	  M,N,z71,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_blas = 0.;
      for (int k=0; k<nloops2; ++k) e10_blas += NormSq(E3[k]-E0[k]);
      e10_blas = sqrt(e10_blas/nloops2);
      e10_blas /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(B4[k].transpose()) * EPART2(D4[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e10_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e10_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e10_eigen = sqrt(e10_eigen/nloops2);
      e10_eigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(B5[k].transpose()) * EPART2(D5[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e10_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e10_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e10_smalleigen = sqrt(e10_smalleigen/nloops2);
      e10_smalleigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * B * D.Adjoint()
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(B0[n2].Transpose()) * D0[n2].Adjoint().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * B0[n2].Transpose().row(i) * MPART2(D0[n2].Adjoint());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Adjoint());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(B0[n2].Transpose()) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(B1[k].Transpose()) * MPART2(D1[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_reg = 0.;
      for (int k=0; k<nloops2; ++k) e11_reg += NormSq(E1[k]-E0[k]);
      e11_reg = sqrt(e11_reg/nloops2);
      e11_reg /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(B2[k].Transpose()) * MPART2(D2[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_small = 0.;
      for (int k=0; k<nloops2; ++k) e11_small += NormSq(E2[k]-E0[k]);
      e11_small = sqrt(e11_small/nloops2);
      e11_small /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLAST, BLASCT,
	  M,N,K,z89,BP(B3[k].cptr()),K,BP(D3[k].cptr()),N,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(D3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLAST, BLASDIAG1,
	  M,N,z89,BP(B3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(B3[k].Transpose());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLASCT, BLASDIAG2,
	  M,N,z89,BP(D3[k].cptr()),N,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_blas = 0.;
      for (int k=0; k<nloops2; ++k) e11_blas += NormSq(E3[k]-E0[k]);
      e11_blas = sqrt(e11_blas/nloops2);
      e11_blas /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(B4[k].transpose()) * EPART2(D4[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e11_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e11_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e11_eigen = sqrt(e11_eigen/nloops2);
      e11_eigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(B5[k].transpose()) * EPART2(D5[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e11_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e11_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e11_smalleigen = sqrt(e11_smalleigen/nloops2);
      e11_smalleigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E = (7,1) * B.Adjoint() * D.Adjoint()
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) = T(7,1) * MPART1(B0[n2].Adjoint()) * D0[n2].Adjoint().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) = T(7,1) * B0[n2].Adjoint().row(i) * MPART2(D0[n2].Adjoint());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Adjoint());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(7,1) * MPART1(B0[n2].Adjoint()) * temp2.col(j);
	MPART3(E0[n2]) = MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) = T(7,1) * MPART1(B1[k].Adjoint()) * MPART2(D1[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_reg = 0.;
      for (int k=0; k<nloops2; ++k) e12_reg += NormSq(E1[k]-E0[k]);
      e12_reg = sqrt(e12_reg/nloops2);
      e12_reg /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) = T(7,1) * MPART1(B2[k].Adjoint()) * MPART2(D2[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_small = 0.;
      for (int k=0; k<nloops2; ++k) e12_small += NormSq(E2[k]-E0[k]);
      e12_small = sqrt(e12_small/nloops2);
      e12_small /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASCT, BLASCT,
	  M,N,K,z71,BP(B3[k].cptr()),K,BP(D3[k].cptr()),N,
	  zero,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART2(D3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLASCT, BLASDIAG1,
	  M,N,z71,BP(B3[k].cptr()),K,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      MPART3(E3[k]) = MPART1(B3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLASCT, BLASDIAG2,
	  M,N,z71,BP(D3[k].cptr()),N,BP(E3[k].ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_blas = 0.;
      for (int k=0; k<nloops2; ++k) e12_blas += NormSq(E3[k]-E0[k]);
      e12_blas = sqrt(e12_blas/nloops2);
      e12_blas /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) = T(7,1) * EPART1(B4[k].adjoint()) * EPART2(D4[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e12_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e12_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e12_eigen = sqrt(e12_eigen/nloops2);
      e12_eigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) = T(7,1) * EPART1(B5[k].adjoint()) * EPART2(D5[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t12_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e12_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e12_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e12_smalleigen = sqrt(e12_smalleigen/nloops2);
      e12_smalleigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif

#if 1 // E += (8,9) * B.Adjoint() * D.Adjoint()
#ifdef TISCOMPLEX
    ClearCache();

#ifdef ERRORCHECK
    if (n == 0) {
      for (int n2=0; n2<nloops2; ++n2) {
#if (PART2 == 1)
	for(int j=0;j<N;++j) 
	  E0[n2].col(j) += T(8,9) * MPART1(B0[n2].Adjoint()) * D0[n2].Adjoint().col(j);
#elif (PART1 == 1) 
	for(int i=0;i<M;++i) 
	  E0[n2].row(i) += T(8,9) * B0[n2].Adjoint().row(i) * MPART2(D0[n2].Adjoint());
#else
	tmv::Matrix<T> temp2 = MPART2(D0[n2].Adjoint());
	tmv::Matrix<T> temp3(M,N);
	for(int j=0;j<N;++j) 
	  temp3.col(j) = T(8,9) * MPART1(B0[n2].Adjoint()) * temp2.col(j);
	MPART3(E0[n2]) += MPART3(temp3);
#endif
      }
    }
#endif

#ifdef DOREG
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E1[k]) += T(8,9) * MPART1(B1[k].Adjoint()) * MPART2(D1[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_reg += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_reg = 0.;
      for (int k=0; k<nloops2; ++k) e13_reg += NormSq(E1[k]-E0[k]);
      e13_reg = sqrt(e13_reg/nloops2);
      e13_reg /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      MPART3(E2[k]) += T(8,9) * MPART1(B2[k].Adjoint()) * MPART2(D2[k].Adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_small += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_small = 0.;
      for (int k=0; k<nloops2; ++k) e13_small += NormSq(E2[k]-E0[k]);
      e13_small = sqrt(e13_small/nloops2);
      e13_small /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOBLAS
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k) {
#if (PART1 == 1) && (PART2 == 1)
      BLASNAME(gemm) (BLASCM BLASCT, BLASCT,
	  M,N,K,z89,BP(B3[k].cptr()),K,BP(D3[k].cptr()),N,
	  one,BP(E3[k].ptr()),M BLAS1 BLAS1);
#elif (PART1 >= 2 && PART1 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART2(D3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASLeft, BLASUPLO1x, BLASCT, BLASDIAG1,
	  M,N,z89,BP(B3[k].cptr()),K,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#elif (PART2 >= 2 && PART2 <= 5) && (PART3 == 1)
      tmv::Matrix<T> temp = MPART1(B3[k].Adjoint());
      BLASNAME(trmm) (BLASCM BLASRight, BLASUPLO2x, BLASCT, BLASDIAG2,
	  M,N,z89,BP(D3[k].cptr()),N,BP(temp.ptr()),M 
	  BLAS1 BLAS1 BLAS1 BLAS1);
      MPART3(E3[k]) += MPART3(temp);
#endif
    }

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_blas += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_blas = 0.;
      for (int k=0; k<nloops2; ++k) e13_blas += NormSq(E3[k]-E0[k]);
      e13_blas = sqrt(e13_blas/nloops2);
      e13_blas /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGEN
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2; ++k)
      EPART3(E4[k]) += T(8,9) * EPART1(B4[k].adjoint()) * EPART2(D4[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_eigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0) {
      e13_eigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e13_eigen += tmv::TMV_NORM(E4[k](i,j)-E0[k](i,j));
      e13_eigen = sqrt(e13_eigen/nloops2);
      e13_eigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif

#ifdef DOEIGENSMALL
    gettimeofday(&tp,0);
    ta = tp.tv_sec + tp.tv_usec/1.e6;

    for (int k=0; k<nloops2x; ++k)
      EPART3(E5[k]) += T(8,9) * EPART1(B5[k].adjoint()) * EPART2(D5[k].adjoint());

    gettimeofday(&tp,0);
    tb = tp.tv_sec + tp.tv_usec/1.e6;
    t13_smalleigen += tb-ta;
#ifdef ERRORCHECK
    if (n == 0 && nloops2x) {
      e13_smalleigen = 0.;
      for (int k=0; k<nloops2; ++k) 
	for (int i=0;i<M;++i) for(int j=0;j<N;++j)
	  e13_smalleigen += tmv::TMV_NORM(E5[k](i,j)-E0[k](i,j));
      e13_smalleigen = sqrt(e13_smalleigen/nloops2);
      e13_smalleigen /= Norm(B0[0])*Norm(D0[0]);
    }
#endif
#endif
#endif
#endif
  }

  std::cout<<"E = B * D               "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
  std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
  std::cout<<"E = -B * D              "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
  std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
  std::cout<<"E = 7 * B * D           "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
  std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
  std::cout<<"E -= B * D              "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
  std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
  std::cout<<"E += 8 * B * D          "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
  std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
#ifdef TISCOMPLEX
  std::cout<<"E = (7,1) * B * D       "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
  std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B * D      "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
  std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * B* * D      "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
  std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B* * D     "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
  std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * B * D*      "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
  std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B * D*     "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
  std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * B* * D*     "<<t12_reg<<"  "<<t12_small<<"  "<<t12_blas;
  std::cout<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B* * D*    "<<t13_reg<<"  "<<t13_small<<"  "<<t13_blas;
  std::cout<<"  "<<t13_eigen<<"  "<<t13_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
  std::cout<<"errors:\n";
  std::cout<<"E = B * D               "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
  std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
  std::cout<<"E = -B * D              "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
  std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
  std::cout<<"E = 7 * B * D           "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
  std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
  std::cout<<"E -= B * D              "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
  std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
  std::cout<<"E += 8 * B * D          "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
  std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
#ifdef TISCOMPLEX
  std::cout<<"E = (7,1) * B * D       "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
  std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B * D      "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
  std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * B* * D      "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
  std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B* * D     "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
  std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * B * D*      "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
  std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B * D*     "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
  std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
  std::cout<<"E = (7,1) * B* * D*     "<<e12_reg<<"  "<<e12_small<<"  "<<e12_blas;
  std::cout<<"  "<<e12_eigen<<"  "<<e12_smalleigen<<std::endl;
  std::cout<<"E += (8,9) * B* * D*    "<<e13_reg<<"  "<<e13_small<<"  "<<e13_blas;
  std::cout<<"  "<<e13_eigen<<"  "<<e13_smalleigen<<std::endl;
#endif
  std::cout<<"\n\n";
#endif
}
#endif


#ifdef SIMPLE_VALUES
#ifdef TISCOMPLEX
//#define RAND ( T(1,2) )
#define RAND ( T(i,j) )
#else
//#define RAND ( T(1) )
#define RAND ( T(i) )
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

#ifdef MEMDEBUG
  atexit(&DumpUnfreed);
#endif

  std::vector<tmv::Matrix<T> > A(nloops2,tmv::Matrix<T>(M,K));
  std::vector<tmv::Matrix<T> > B(nloops2,tmv::Matrix<T>(K,M));
  std::vector<tmv::Matrix<T> > C(nloops2,tmv::Matrix<T>(K,N));
  std::vector<tmv::Matrix<T> > D(nloops2,tmv::Matrix<T>(N,K));
  std::vector<tmv::Matrix<T> > E(nloops2,tmv::Matrix<T>(M,N));
  for(int k=0;k<nloops2;++k) {
    for(int i=0;i<M;++i) for(int j=0;j<K;++j) {
      A[k](i,j) = RAND;
      B[k](j,i) = RAND;
    }
    for(int i=0;i<K;++i) for(int j=0;j<N;++j) {
      C[k](i,j) = RAND;
      D[k](j,i) = RAND;
    }
  }
  std::cout<<"M,N,K = "<<M<<" , "<<N<<" , "<<K<<std::endl;
  std::cout<<"nloops = "<<nloops1<<" x "<<nloops2;
  std::cout<<" = "<<nloops1*nloops2<<std::endl;
#ifdef DOEIGENSMALL
  if (nloops2x == 0)
    std::cout<<"Matrix is too big for the stack, so no \"Eigen Known\" tests.\n";
#endif
  std::cout<<"\nTime for:               ";
  std::cout<<"  TMV    "<<" TMV Small"<<"   BLAS   "<<"  Eigen  "<<"Eigen Known"<<std::endl;

#ifdef DOMULTMM_CCC
  MultMM_CCC(A,C,E);
#endif

#ifdef DOMULTMM_CRC
  MultMM_CRC(A,D,E);
#endif

#ifdef DOMULTMM_RCC
  MultMM_RCC(B,C,E);
#endif

#ifdef DOMULTMM_RRC
  MultMM_RRC(B,D,E);
#endif

  return 0;
}
#if 0
catch (int) {}
#else
catch (tmv::Error& e)
{
  std::cout<<"Caught error "<<e<<std::endl;
  return 1;
}
#endif
