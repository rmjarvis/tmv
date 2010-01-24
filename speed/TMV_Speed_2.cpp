
//#define PRINTALGO
//#define PRINTALGO_MV
//#define PRINTALGO_UV
//#define PRINTALGO_UD
//#define PRINTALGO_UM
//#define PRINTALGO_UU
//#define PRINTALGO_MD
//#define PRINTALGO_MM
//#define PRINTALGO_R1
//#define XDEBUG_PRODMV
//#define XDEBUG_PRODMM
//#define XDEBUG_OPRODVV
#undef NDEBUG

#include <iostream>
#include "TMV.h"

// How big do you want the matrices to be?
const int M = 4;
const int N = M;
#define AISSQUARE

// Define the type to use:
#define TISFLOAT
#define TISCOMPLEX

// Define the part of the matrix to use
// 1 = all (rectangle matrix)
// 2 = uppertri
// 3 = unituppertri
// 4 = lowertri
// 5 = unitlowertri
#define PART 1

// Define whether you want to include the error checks.
#define ERRORCHECK

// Define which versions you want to test:
#define DOREG
#define DOSMALL
#define DOBLAS
#define DOEIGEN
#define DOEIGENSMALL

// Define which batches of functions you want to test:
//#define DOMULTXM
//#define DOADDMM
//#define DONORM
//#define DOMULTMV
#define DORANK1
//#define DOMULTMD
//#define DOMULTDM
//#define DOSWAP

// Set this if you only want to do a single loop.
// Not so useful for timing, but useful for debugging.
//#define ONELOOP

// Normally I fill the matrices and vectors with random values, but 
// that can make it difficult to debug the arithmetic.
// This define uses simpler values that you can multiply in your head.
//#define SIMPLE_VALUES

// My algorithms for MultMV are designed to be fast when the vector has lots
// of zeros.  This option sets the first and last third of the vector to 0's.
#define ZERO_VALUES

// Set up the target number of operations and memory to use for testing
const int targetnflops = 100000000; // in real ops
const int targetmem = 1000000; // in bytes

// Include the BLAS library you want to test TMV against:
//#include "mkl.h"
extern "C" {
#include "util/fblas.h"
#include "util/flapack.h"
}

// Set these as appropriate given the BLAS and LAPACK libraries included above.
//#define CBLAS
//#define CLAPACK
#define FBLAS
#define FLAP
#define XLAP


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
    (M*N*sizeof(T) > targetmem ? 1 : targetmem / (M*N) / (sizeof(T))));
const int nloops1 = targetnflops / (M * N * nloops2) / XFOUR;
#endif

#if (PART == 1) // Full rectangle
#define INPART(i,j) true
#define UNITPART(i,j) false
#define MPART(m) m
#define MPART2(m) m
#define EPART(m) m
#define EPART2(m) m
#define EPART3(m) m
#define DO_PERMUTE
#define DO_TRANSPOSE_SELF
#define COL_LEN(j) M
#define COL_LEN2(j) M
#define COL_START(j) 0
#define COL_END(j) M
#define ROW_LEN(i) N
#define ROW_START(i) 0
#define ROW_END(i) N

#elif (PART == 2) // UpperTri
#define INPART(i,j) (i<=j)
#define UNITPART(i,j) false
#define MPART(m) m.upperTri().viewAsUnknownDiag()
#define MPART2(m) m.upperTri().viewAsUnknownDiag()
#define EPART(m) m.part<Eigen::UpperTriangular>()
#define EPART2(m) m.part<Eigen::UpperTriangular>()
#define EPART3(m) m.part<Eigen::UpperTriangular>()
#define COL_LEN(j) (j+1)
#define COL_LEN2(j) (j+1)
#define COL_START(j) 0
#define COL_END(j) (j+1)
#define ROW_LEN(i) (N-i)
#define ROW_START(i) i
#define ROW_END(i) N
#define UPLO_CHAR 'U'
#define DIAG_CHAR 'N'

#elif (PART == 3) // UnitUpperTri
#define INPART(i,j) (i<j)
#define UNITPART(i,j) (i==j)
#define MPART(m) m.unitUpperTri().viewAsUnknownDiag()
#define MPART2(m) m.upperTri().viewAsUnknownDiag()
#define EPART(m) m.part<Eigen::UnitUpperTriangular>()
#define EPART2(m) m.part<Eigen::UpperTriangular>()
#define EPART3(m) m.part<Eigen::StrictlyUpperTriangular>()
#define COL_LEN(j) j
#define COL_LEN2(j) j
#define COL_START(j) 0
#define COL_END(j) j
#define ROW_LEN(i) (N-i-1)
#define ROW_START(i) (i+1)
#define ROW_END(i) N
#define UPLO_CHAR 'U'
#define DIAG_CHAR 'U'

#elif (PART == 4) // LowerTri
#define INPART(i,j) (i>=j)
#define UNITPART(i,j) false
#define MPART(m) m.lowerTri().viewAsUnknownDiag()
#define MPART2(m) m.lowerTri().viewAsUnknownDiag()
#define EPART(m) m.part<Eigen::LowerTriangular>()
#define EPART2(m) m.part<Eigen::LowerTriangular>()
#define EPART3(m) m.part<Eigen::LowerTriangular>()
#define COL_LEN(j) (N-j)
#define COL_LEN2(j) (N-j)
#define COL_START(j) j
#define COL_END(j) N
#define ROW_LEN(i) (i+1)
#define ROW_START(i) 0
#define ROW_END(i) (i+1)
#define UPLO_CHAR 'L'
#define DIAG_CHAR 'N'

#elif (PART == 5) // UnitLowerTri
#define INPART(i,j) (i>j)
#define UNITPART(i,j) (i==j)
#define MPART(m) m.unitLowerTri().viewAsUnknownDiag()
#define MPART2(m) m.lowerTri().viewAsUnknownDiag()
#define EPART(m) m.part<Eigen::UnitLowerTriangular>()
#define EPART2(m) m.part<Eigen::LowerTriangular>()
#define EPART3(m) m.part<Eigen::StrictlyLowerTriangular>()
#define COL_LEN(j) (N-j-1)
#define COL_LEN2(j) (N-j-1)
#define COL_START(j) (j+1)
#define COL_END(j) N
#define ROW_LEN(i) i
#define ROW_START(i) 0
#define ROW_END(i) i
#define UPLO_CHAR 'L'
#define DIAG_CHAR 'U'

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
#if ((PART == 2) || (PART == 3))
#define BLASUPLO CblasUpper
#define BLASUPLO2 CblasLower
#elif ((PART == 4) || (PART == 5))
#define BLASUPLO CblasLower
#define BLASUPLO2 CblasUpper
#endif
#if ((PART == 2) || (PART == 4))
#define BLASDIAG CblasNonUnit
#elif ((PART == 3) || (PART == 5))
#define BLASDIAG CblasUnit
#endif
#define BLAS1
#elif defined FBLAS
const char BlasCh_N = 'N';
const char BlasCh_T = 'T';
const char BlasCh_C = 'C';
const char BlasCh_U = 'U';
const char BlasCh_L = 'L';
#define BLASNAME1(c,x) c ## x ## _
#define BLASINAME1(c,x) i ## c ## x ## _
#define BLASM1 -1
#define BLAS1 ,1
#define BLASCM 
#define BLASNT BlasCh_N
#define BLAST BlasCh_T
#define BLASCT BlasCh_C
#if ((PART == 2) || (PART == 3))
#define BLASUPLO BlasCh_U
#define BLASUPLO2 BlasCh_L
#elif ((PART == 4) || (PART == 5))
#define BLASUPLO BlasCh_L
#define BLASUPLO2 BlasCh_U
#endif
#if ((PART == 2) || (PART == 4))
#define BLASDIAG BlasCh_N
#elif ((PART == 3) || (PART == 5))
#define BLASDIAG BlasCh_U
#endif
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

#ifdef CLAP
#define LAPNAME1(c,x) clapack_ ## c ## x
#define LAPCM CblasColMajor,
RT lap_rwork[M*N];
T lap_work[M*N];
#define LAPWORK lap_work
#define LAPRWORK lap_rwork
#elif defined FLAP
const char LapCh_N = 'N';
const char LapCh_T = 'T';
const char LapCh_C = 'C';
#define LAPNAME1(c,x) c ## x ## _
#define LAPCM 
RT lap_rwork[M*N];
T lap_work[M*N];
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
#define EIGENV VectorXcf
#define EIGENM MatrixXcf
#else
#define EIGENV VectorXcd
#define EIGENM MatrixXcd
#endif
#else
#ifdef TISFLOAT
#define EIGENV VectorXf
#define EIGENM MatrixXf
#else
#define EIGENV VectorXd
#define EIGENM MatrixXd
#endif
#endif

#endif

#ifdef DOEIGENSMALL
const int nloops2x = (
    (M*N*sizeof(T) < 512 * 1024 && N!=10000 && M!=10000) ? nloops2 : 0 );
#endif

#include <sys/time.h>
#include <algorithm>
#include <numeric>
#include <iostream>

#if (defined(DOEIGEN) || defined(DOEIGENSMALL))
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

#ifdef DOMULTXM
static void MultXM(
    const std::vector<tmv::Matrix<T> >& A1,
    const std::vector<tmv::Matrix<T> >& B1)
{
    std::vector<tmv::Matrix<T> > C1 = A1;

#ifdef ERRORCHECK
    std::vector<tmv::Matrix<T> > A0 = A1;
    std::vector<tmv::Matrix<T> > B0 = B1;
    std::vector<tmv::Matrix<T> > C0 = A1;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallMatrix<T,M,N> > A2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,M> > B2(nloops2);
    std::vector<tmv::SmallMatrix<T,M,N> > C2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; B2[k] = B1[k]; C2[k] = C1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Matrix<T> > A3 = A1;
    std::vector<tmv::Matrix<T> > B3 = B1;
    std::vector<tmv::Matrix<T> > C3 = C1;
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENM> A4(nloops2,Eigen::EIGENM(M,N));
    std::vector<Eigen::EIGENM> B4(nloops2,Eigen::EIGENM(N,M));
    std::vector<Eigen::EIGENM> C4(nloops2,Eigen::EIGENM(M,N));
    for(int k=0;k<nloops2;++k) {
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) A4[k](i,j) = A1[k](i,j);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) B4[k](j,i) = B1[k](j,i);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) C4[k](j,i) = C1[k](j,i);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,M,N> > A5;
    std::vector<Eigen::Matrix<T,N,M> > B5;
    std::vector<Eigen::Matrix<T,M,N> > C5;
    if (nloops2x) {
        A5.resize(nloops2x); B5.resize(nloops2x); C5.resize(nloops2x);
    }
    for(int k=0;k<nloops2x;++k) {
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) A5[k](i,j) = A1[k](i,j);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) B5[k](j,i) = B1[k](j,i);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) C5[k](j,i) = C1[k](j,i);
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

#if 1 // C = A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = A0[k](i,j);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(C1[k]) = MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(C1[k]-C0[k]);
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(C2[k]) = MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(C2[k]-C0[k]);
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j)
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(C3[k]-C0[k]);
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            //EPART(C4[k]) = EPART(A4[k]);
            // Eigen won't compile the above statement for triangle parts.
            // They forgot to overwrite the defaulat assignment operator,
            // which doesn't work here.
            // However, the below statement does work, even though it looks like
            // it should give an error, since the shapes don't match.
            // Eigen uses the convention that the LHS determines the shape.
            // Also, it doesn't support using a UnitDiagonal triangle matrix
            // on the LHS, even when the RHS is also UnitDiagonal, so there is
            // no conflict, so we have to use StrictlyUpperTriangle or
            // StrictlyLowerTriangle instead.  That is done in EPART3.
            EPART3(C4[k]) = A4[k];
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e1_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            EPART3(C5[k]) = A5[k];
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e1_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // C = 8 * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = RT(8) * A0[k](i,j);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = RT(8);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = 8 * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) e2_reg += NormSq(C1[k]-C0[k]);
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= RT(8)*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = 8 * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) e2_small += NormSq(C2[k]-C0[k]);
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= RT(8)*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(scal) (M*N,eight,BP(C3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(scal) (
                    COL_LEN(j),eight,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(8));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) e2_blas += NormSq(C3[k]-C0[k]);
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= RT(8)*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = 8 * EPART(A4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e2_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= RT(8)*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = 8 * EPART(A5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e2_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= RT(8)*Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // C *= 7
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) *= RT(7);
                    else if (UNITPART(i,j))
                        C0[k](i,j) *= RT(7);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) *= 7;

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            MPART2(C2[k]) *= 7;

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
#if (PART == 1)
            BLASNAME(scal) (M*N,seven,BP(C3[k].ptr()),1);
#elif ((PART == 2) || (PART == 3))
            for(int j=0;j<N;++j) 
                BLASNAME(scal) (j+1,seven,BP(C3[k].col(j,0,j+1).ptr()),1);
#elif ((PART == 4) || (PART == 5))
            for(int j=0;j<N;++j)
                BLASNAME(scal) (N-j,seven,BP(C3[k].col(j,j,N).ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            EPART2(C4[k]) *= 7;

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e3_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) *= 7;

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e3_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // C = -A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = -A0[k](i,j);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = -RT(1);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = -MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) e4_reg += NormSq(C1[k]-C0[k]);
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = -MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) e4_small += NormSq(C2[k]-C0[k]);
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(scal) (M*N,mone,BP(C3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(scal) (
                    COL_LEN(j),mone,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(-1));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) e4_blas += NormSq(C3[k]-C0[k]);
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = -EPART(A4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e4_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = -EPART(A5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e4_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // C = -C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = -C0[k](i,j);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = -RT(1);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = -MPART(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) e5_reg += NormSq(C1[k]-C0[k]);
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = -MPART(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) e5_small += NormSq(C2[k]-C0[k]);
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(scal) (M*N,mone,BP(C3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(scal) (
                    COL_LEN(j),mone,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(-1));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) e5_blas += NormSq(C3[k]-C0[k]);
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = -EPART(C4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e5_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = -EPART(C5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e5_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#ifdef TISCOMPLEX
#if 1 // C *= complex(8,9)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) *= T(8,9);
                    else if (UNITPART(i,j))
                        C0[k](i,j) *= T(8,9);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) *= T(8,9);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) e6_reg += NormSq(C1[k]-C0[k]);
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) *= T(8,9);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) e6_small += NormSq(C2[k]-C0[k]);
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(scal) (M*N,z89,BP(C3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(scal) (
                    COL_LEN(j),z89,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag() *= T(8,9);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) e6_blas += NormSq(C3[k]-C0[k]);
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) *= T(8,9);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e6_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) *= T(8,9);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e6_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // C = complex(8,9) * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = T(8,9) * A0[k](i,j);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = T(8,9);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = T(8,9) * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) e7_reg += NormSq(C1[k]-C0[k]);
            e7_reg = sqrt(e7_reg/nloops2);
            e7_reg /= RT(12)*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = T(8,9) * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) e7_small += NormSq(C2[k]-C0[k]);
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= RT(12)*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(scal) (M*N,z89,BP(C3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(scal) (
                    COL_LEN(j),z89,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(8,9));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) e7_blas += NormSq(C3[k]-C0[k]);
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= RT(12)*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = T(8,9) * EPART(A4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e7_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= RT(12)*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = T(8,9) * EPART(A5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e7_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= RT(12)*Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // C = complex(8,9) * A*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = T(8,9) * std::conj(A0[k](i,j));
                    else if (UNITPART(i,j))
                        C0[k](i,j) = T(8,9);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = T(8,9) * MPART(A1[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) e8_reg += NormSq(C1[k]-C0[k]);
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= RT(12)*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = T(8,9) * MPART(A2[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) e8_small += NormSq(C2[k]-C0[k]);
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= RT(12)*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASDNAME(scal) (M*N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(scal) (M*N,z89,BP(C3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASDNAME(scal) (
                    COL_LEN(j),dmone,
                    C3[k].col(j,COL_START(j),COL_END(j)).realPart().ptr()+1,2);
                BLASNAME(scal) (
                    COL_LEN(j),z89,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(8,9));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) e8_blas += NormSq(C3[k]-C0[k]);
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= RT(12)*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = T(8,9) * EPART(A4[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e8_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= RT(12)*Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = T(8,9) * EPART(A5[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e8_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= RT(12)*Norm(A0[0]);
        }
#endif
#endif
#endif
#endif

#if 1 // C = B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = B0[k](j,i);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(C1[k]) = MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_reg = 0.;
            for (int k=0; k<nloops2; ++k) e9_reg += NormSq(C1[k]-C0[k]);
            e9_reg = sqrt(e9_reg/nloops2);
            e9_reg /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(C2[k]) = MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_small = 0.;
            for (int k=0; k<nloops2; ++k) e9_small += NormSq(C2[k]-C0[k]);
            e9_small = sqrt(e9_small/nloops2);
            e9_small /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_blas = 0.;
            for (int k=0; k<nloops2; ++k) e9_blas += NormSq(C3[k]-C0[k]);
            e9_blas = sqrt(e9_blas/nloops2);
            e9_blas /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART3(C4[k]) = EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e9_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e9_eigen = sqrt(e9_eigen/nloops2);
            e9_eigen /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART3(C5[k]) = EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e9_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e9_smalleigen = sqrt(e9_smalleigen/nloops2);
            e9_smalleigen /= Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = 8 * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = RT(8) * B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = RT(8);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = 8 * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_reg = 0.;
            for (int k=0; k<nloops2; ++k) e10_reg += NormSq(C1[k]-C0[k]);
            e10_reg = sqrt(e10_reg/nloops2);
            e10_reg /= RT(8)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = 8 * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_small = 0.;
            for (int k=0; k<nloops2; ++k) e10_small += NormSq(C2[k]-C0[k]);
            e10_small = sqrt(e10_small/nloops2);
            e10_small /= RT(8)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(scal) (
                    COL_LEN(j),eight,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
            }
#if (PART == 1)
            BLASNAME(scal) (M*N,eight,BP(C3[k].ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(8));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_blas = 0.;
            for (int k=0; k<nloops2; ++k) e10_blas += NormSq(C3[k]-C0[k]);
            e10_blas = sqrt(e10_blas/nloops2);
            e10_blas /= RT(8)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = 8 * EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e10_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e10_eigen = sqrt(e10_eigen/nloops2);
            e10_eigen /= RT(8)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = 8 * EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e10_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e10_smalleigen = sqrt(e10_smalleigen/nloops2);
            e10_smalleigen /= RT(8)*Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = -B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = -B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = -RT(1);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = -MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_reg = 0.;
            for (int k=0; k<nloops2; ++k) e11_reg += NormSq(C1[k]-C0[k]);
            e11_reg = sqrt(e11_reg/nloops2);
            e11_reg /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = -MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_small = 0.;
            for (int k=0; k<nloops2; ++k) e11_small += NormSq(C2[k]-C0[k]);
            e11_small = sqrt(e11_small/nloops2);
            e11_small /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(scal) (
                    COL_LEN(j),mone,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
            }
#if (PART == 1)
            BLASNAME(scal) (M*N,mone,BP(C3[k].ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(-1));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_blas = 0.;
            for (int k=0; k<nloops2; ++k) e11_blas += NormSq(C3[k]-C0[k]);
            e11_blas = sqrt(e11_blas/nloops2);
            e11_blas /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = -EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e11_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e11_eigen = sqrt(e11_eigen/nloops2);
            e11_eigen /= Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = -EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e11_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e11_smalleigen = sqrt(e11_smalleigen/nloops2);
            e11_smalleigen /= Norm(B0[0]);
        }
#endif
#endif
#endif

#ifdef TISCOMPLEX
#if 1 // C = complex(8,9) * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = T(8,9) * B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = T(8,9);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = T(8,9) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_reg = 0.;
            for (int k=0; k<nloops2; ++k) e12_reg += NormSq(C1[k]-C0[k]);
            e12_reg = sqrt(e12_reg/nloops2);
            e12_reg /= RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = T(8,9) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_small = 0.;
            for (int k=0; k<nloops2; ++k) e12_small += NormSq(C2[k]-C0[k]);
            e12_small = sqrt(e12_small/nloops2);
            e12_small /= RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(scal) (
                    COL_LEN(j),z89,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
            }
#if (PART == 1)
            BLASNAME(scal) (M*N,z89,BP(C3[k].ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(8,9));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_blas = 0.;
            for (int k=0; k<nloops2; ++k) e12_blas += NormSq(C3[k]-C0[k]);
            e12_blas = sqrt(e12_blas/nloops2);
            e12_blas /= RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = T(8,9) * EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e12_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e12_eigen = sqrt(e12_eigen/nloops2);
            e12_eigen /= RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = T(8,9) * EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e12_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e12_smalleigen = sqrt(e12_smalleigen/nloops2);
            e12_smalleigen /= RT(12)*Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = complex(8,9) * B*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = T(8,9) * std::conj(B0[k](j,i));
                    else if (UNITPART(i,j))
                        C0[k](i,j) = T(8,9);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = T(8,9) * MPART(B1[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_reg = 0.;
            for (int k=0; k<nloops2; ++k) e13_reg += NormSq(C1[k]-C0[k]);
            e13_reg = sqrt(e13_reg/nloops2);
            e13_reg /= RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = T(8,9) * MPART(B2[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_small = 0.;
            for (int k=0; k<nloops2; ++k) e13_small += NormSq(C2[k]-C0[k]);
            e13_small = sqrt(e13_small/nloops2);
            e13_small /= RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#if ((PART >= 2) && (PART <= 5))
                BLASDNAME(scal) (
                    COL_LEN(j),dmone,
                    C3[k].col(j,COL_START(j),COL_END(j)).realPart().ptr()+1,2);
                BLASNAME(scal) (
                    COL_LEN(j),z89,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
            }
#if (PART == 1)
            BLASDNAME(scal) (M*N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(scal) (M*N,z89,BP(C3[k].ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(8,9));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_blas = 0.;
            for (int k=0; k<nloops2; ++k) e13_blas += NormSq(C3[k]-C0[k]);
            e13_blas = sqrt(e13_blas/nloops2);
            e13_blas /= RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = T(8,9) * EPART(B4[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e13_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e13_eigen = sqrt(e13_eigen/nloops2);
            e13_eigen /= RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = T(8,9) * EPART(B5[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e13_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e13_smalleigen = sqrt(e13_smalleigen/nloops2);
            e13_smalleigen /= RT(12)*Norm(B0[0]);
        }
#endif
#endif
#endif
#endif
    }

    std::cout<<"C = A                   "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"C = 8*A                 "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"C *= 7                  "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"C = -A                  "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"C = -C                  "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C *= (8,9)              "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"C = (8,9) * A           "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"C = (8,9) * A*          "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
#endif
    std::cout<<"C = B                   "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
    std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
    std::cout<<"C = 8*B                 "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
    std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
    std::cout<<"C = -B                  "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
    std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C = (8,9) * B           "<<t12_reg<<"  "<<t12_small<<"  "<<t12_blas;
    std::cout<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;
    std::cout<<"C = (8,9) * B*          "<<t13_reg<<"  "<<t13_small<<"  "<<t13_blas;
    std::cout<<"  "<<t13_eigen<<"  "<<t13_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"C = A                   "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"C = 8 * A               "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"C *= 7                  "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"C = -A                  "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"C = -C                  "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C *= (8,9)              "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
    std::cout<<"C = (8,9) * A           "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
    std::cout<<"C = (8,9) * A*          "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
#endif
    std::cout<<"C = B                   "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
    std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
    std::cout<<"C = 8*B                 "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
    std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
    std::cout<<"C = -B                  "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
    std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C = (8,9) * B           "<<e12_reg<<"  "<<e12_small<<"  "<<e12_blas;
    std::cout<<"  "<<e12_eigen<<"  "<<e12_smalleigen<<std::endl;
    std::cout<<"C = (8,9) * B*          "<<e13_reg<<"  "<<e13_small<<"  "<<e13_blas;
    std::cout<<"  "<<e13_eigen<<"  "<<e13_smalleigen<<std::endl;
#endif
    std::cout<<"\n\n";
#endif

}
#endif

#ifdef DOADDMM
static void AddMM(
    const std::vector<tmv::Matrix<T> >& A1,
    const std::vector<tmv::Matrix<T> >& B1)
{
    std::vector<tmv::Matrix<T> > A1x = A1;
    std::vector<tmv::Matrix<T> > B1x = B1;
    std::vector<tmv::Matrix<T> > C1 = A1;

#ifdef ERRORCHECK
    std::vector<tmv::Matrix<T> > A0 = A1;
    std::vector<tmv::Matrix<T> > B0 = B1;
    std::vector<tmv::Matrix<T> > C0 = C1;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallMatrix<T,M,N> > A2(nloops2);
    std::vector<tmv::SmallMatrix<T,M,N> > A2x(nloops2);
    std::vector<tmv::SmallMatrix<T,N,M> > B2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,M> > B2x(nloops2);
    std::vector<tmv::SmallMatrix<T,M,N> > C2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; B2[k] = B1[k];
        A2x[k] = A1x[k]; B2x[k] = B1x[k];
        C2[k] = C1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Matrix<T> > A3 = A1;
    std::vector<tmv::Matrix<T> > A3x = A1x;
    std::vector<tmv::Matrix<T> > B3 = B1;
    std::vector<tmv::Matrix<T> > B3x = B1x;
    std::vector<tmv::Matrix<T> > C3 = C1;
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENM> A4(nloops2,Eigen::EIGENM(M,N));
    std::vector<Eigen::EIGENM> A4x(nloops2,Eigen::EIGENM(M,N));
    std::vector<Eigen::EIGENM> B4(nloops2,Eigen::EIGENM(N,M));
    std::vector<Eigen::EIGENM> B4x(nloops2,Eigen::EIGENM(N,M));
    std::vector<Eigen::EIGENM> C4(nloops2,Eigen::EIGENM(M,N));
    for(int k=0;k<nloops2;++k) {
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
            A4[k](i,j) = A1[k](i,j);
            A4x[k](i,j) = A1x[k](i,j);
            B4[k](j,i) = B1[k](j,i);
            B4x[k](j,i) = B1x[k](j,i);
            C4[k](i,j) = C1[k](i,j);
        }
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,M,N> > A5;
    std::vector<Eigen::Matrix<T,M,N> > A5x;
    std::vector<Eigen::Matrix<T,N,M> > B5;
    std::vector<Eigen::Matrix<T,N,M> > B5x;
    std::vector<Eigen::Matrix<T,M,N> > C5;
    if (nloops2x) 
    { 
        A5.resize(nloops2x); A5x.resize(nloops2x); 
        B5.resize(nloops2x); B5x.resize(nloops2x); C5.resize(nloops2x); 
    }
    for(int k=0;k<nloops2x;++k) {
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
            A5[k](i,j) = A1[k](i,j);
            A5x[k](i,j) = A1x[k](i,j);
            B5[k](j,i) = B1[k](j,i);
            B5x[k](j,i) = B1x[k](j,i);
            C5[k](i,j) = C1[k](i,j);
        }
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
#endif

    for (int n=0; n<nloops1; ++n) {

#if 1 // C += A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) += A0[k](i,j);
                    else if (UNITPART(i,j))
                        C0[k](i,j) += RT(1);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) += MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(C1[k]-C0[k]);
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) += MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(C2[k]-C0[k]);
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(axpy) (M*N,one,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(axpy) (
                    COL_LEN(j),one,
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().addToAll(T(1));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(C3[k]-C0[k]);
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) += EPART(A4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e1_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) += EPART(A5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e1_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // C += B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) += B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) += RT(1);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) += MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            MPART2(C2[k]) += MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(axpy) (
                    COL_LEN(j),one,
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().addToAll(T(1));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            EPART2(C4[k]) += EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e2_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) += EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e2_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // C -= B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) -= B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) -= RT(1);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) -= MPART(B1[k].transpose());


        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            MPART2(C2[k]) -= MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(axpy) (
                    COL_LEN(j),mone,
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().addToAll(T(-1));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            EPART2(C4[k]) -= EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e3_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) -= EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e3_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // C += 8 * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) += RT(8) * B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) += RT(8);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) += 8 * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            MPART2(C2[k]) += 8 * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(axpy) (
                    COL_LEN(j),eight,
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().addToAll(T(8));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            EPART2(C4[k]) += 8 * EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e4_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) += 8 * EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e4_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // C = A + Ax
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = A0[k](i,j) + A0[k](i,j);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = RT(2);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = MPART(A1[k]) + MPART(A1x[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            MPART2(C2[k]) = MPART(A2[k]) + MPART(A2x[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3x[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(axpy) (M*N,one,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(axpy) (
                    COL_LEN(j),one,
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(2));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            EPART2(C4[k]) = EPART(A4[k]) + EPART(A4x[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e5_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = EPART(A5[k]) + EPART(A5x[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e5_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = A + B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = A0[k](i,j) + B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = RT(2);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = MPART(A1[k]) + MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            MPART2(C2[k]) = MPART(A2[k]) + MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
#endif
            for(int j=0;j<N;++j) {
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
                BLASNAME(axpy) (
                    COL_LEN(j),one,
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(2));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            EPART2(C4[k]) = EPART(A4[k]) + EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e6_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = EPART(A5[k]) + EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e6_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = B + Bx
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = B0[k](j,i) + B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = RT(2);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = MPART(B1[k].transpose()) + MPART(B1x[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            MPART2(C2[k]) = MPART(B2[k].transpose()) + MPART(B2x[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(axpy) (
                    COL_LEN(j),one,
                    BP(B3x[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(2));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            EPART2(C4[k]) = EPART(B4[k].transpose()) + EPART(B4x[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e7_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = EPART(B5[k].transpose()) + EPART(B5x[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e7_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = A - B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = A0[k](i,j) - B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = RT(0);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = MPART(A1[k]) - MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            MPART2(C2[k]) = MPART(A2[k]) - MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
#endif
            for(int j=0;j<N;++j) {
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
                BLASNAME(axpy) (
                    COL_LEN(j),mone,
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(0));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            EPART2(C4[k]) = EPART(A4[k]) - EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e8_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = EPART(A5[k]) - EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e8_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = 7*A - 12*B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = RT(7) * A0[k](i,j) - RT(12) * B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = RT(-5);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = 7 * MPART(A1[k]) - 12 * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_reg = 0.;
            for (int k=0; k<nloops2; ++k) e9_reg += NormSq(C1[k]-C0[k]);
            e9_reg = sqrt(e9_reg/nloops2);
            e9_reg /= RT(7)*Norm(A0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = 7 * MPART(A2[k]) - 12 * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_small = 0.;
            for (int k=0; k<nloops2; ++k) e9_small += NormSq(C2[k]-C0[k]);
            e9_small = sqrt(e9_small/nloops2);
            e9_small /= RT(7)*Norm(A0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(scal) (M*N,seven,BP(C3[k].ptr()),1);
#endif
            for(int j=0;j<N;++j) {
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(scal) (
                    COL_LEN(j),seven,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
                BLASNAME(axpy) (
                    COL_LEN(j),mtwelve,
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(-5));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_blas = 0.;
            for (int k=0; k<nloops2; ++k) e9_blas += NormSq(C3[k]-C0[k]);
            e9_blas = sqrt(e9_blas/nloops2);
            e9_blas /= RT(7)*Norm(A0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = 7 * EPART(A4[k]) - 12 * EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e9_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e9_eigen = sqrt(e9_eigen/nloops2);
            e9_eigen /= RT(7)*Norm(A0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = 7 * EPART(A5[k]) - 12 * EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e9_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e9_smalleigen = sqrt(e9_smalleigen/nloops2);
            e9_smalleigen /= RT(7)*Norm(A0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = 7*C - 12*A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = RT(7) * C0[k](i,j) - RT(12) * A0[k](i,j);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = RT(-5);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = 7 * MPART(C1[k]) - 12 * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_reg = 0.;
            for (int k=0; k<nloops2; ++k) e10_reg += NormSq(C1[k]-C0[k]);
            e10_reg = sqrt(e10_reg/nloops2);
            e10_reg /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = 7 * MPART(C2[k]) - 12 * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_small = 0.;
            for (int k=0; k<nloops2; ++k) e10_small += NormSq(C2[k]-C0[k]);
            e10_small = sqrt(e10_small/nloops2);
            e10_small /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(scal) (M*N,seven,BP(C3[k].ptr()),1);
            BLASNAME(axpy) (M*N,mtwelve,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(scal) (
                    COL_LEN(j),seven,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(axpy) (
                    COL_LEN(j),mtwelve,
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(-5));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_blas = 0.;
            for (int k=0; k<nloops2; ++k) e10_blas += NormSq(C3[k]-C0[k]);
            e10_blas = sqrt(e10_blas/nloops2);
            e10_blas /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = 7 * EPART(C4[k]) - 12 * EPART(A4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e10_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e10_eigen = sqrt(e10_eigen/nloops2);
            e10_eigen /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = 7 * EPART(C5[k]) - 12 * EPART(A5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e10_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e10_smalleigen = sqrt(e10_smalleigen/nloops2);
            e10_smalleigen /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif
#endif

#ifdef TISCOMPLEX
#if 1 // C += complex(8,9) * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) += T(8,9) * A0[k](i,j);
                    else if (UNITPART(i,j))
                        C0[k](i,j) += T(8,9);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) += T(8,9) * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_reg = 0.;
            for (int k=0; k<nloops2; ++k) e11_reg += NormSq(C1[k]-C0[k]);
            e11_reg = sqrt(e11_reg/nloops2);
            e11_reg /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) += T(8,9) * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_small = 0.;
            for (int k=0; k<nloops2; ++k) e11_small += NormSq(C2[k]-C0[k]);
            e11_small = sqrt(e11_small/nloops2);
            e11_small /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(axpy) (M*N,z89,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(axpy) (
                    COL_LEN(j),z89,
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().addToAll(T(8,9));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_blas = 0.;
            for (int k=0; k<nloops2; ++k) e11_blas += NormSq(C3[k]-C0[k]);
            e11_blas = sqrt(e11_blas/nloops2);
            e11_blas /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) += T(8,9) * EPART(A4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e11_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e11_eigen = sqrt(e11_eigen/nloops2);
            e11_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) += T(8,9) * EPART(A5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e11_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e11_smalleigen = sqrt(e11_smalleigen/nloops2);
            e11_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // C += complex(8,9) * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) += T(8,9) * B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) += T(8,9);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) += T(8,9) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_reg = 0.;
            for (int k=0; k<nloops2; ++k) e12_reg += NormSq(C1[k]-C0[k]);
            e12_reg = sqrt(e12_reg/nloops2);
            e12_reg /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) += T(8,9) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_small = 0.;
            for (int k=0; k<nloops2; ++k) e12_small += NormSq(C2[k]-C0[k]);
            e12_small = sqrt(e12_small/nloops2);
            e12_small /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(axpy) (
                    COL_LEN(j),z89,
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            C3[k].diag().addToAll(T(8,9));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_blas = 0.;
            for (int k=0; k<nloops2; ++k) e12_blas += NormSq(C3[k]-C0[k]);
            e12_blas = sqrt(e12_blas/nloops2);
            e12_blas /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) += T(8,9) * EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e12_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e12_eigen = sqrt(e12_eigen/nloops2);
            e12_eigen /= Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) += T(8,9) * EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e12_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e12_smalleigen = sqrt(e12_smalleigen/nloops2);
            e12_smalleigen /= Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // C = 7*A + (8,9)*B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = RT(7)*A0[k](i,j) + T(8,9) * B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = T(15,9);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = 7 * MPART(A1[k]) + T(8,9) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_reg = 0.;
            for (int k=0; k<nloops2; ++k) e13_reg += NormSq(C1[k]-C0[k]);
            e13_reg = sqrt(e13_reg/nloops2);
            e13_reg /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = 7 * MPART(A2[k]) + T(8,9) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_small = 0.;
            for (int k=0; k<nloops2; ++k) e13_small += NormSq(C2[k]-C0[k]);
            e13_small = sqrt(e13_small/nloops2);
            e13_small /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(scal) (M*N,seven,BP(C3[k].ptr()),1);
#endif
            for(int j=0;j<N;++j) {
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(scal) (
                    COL_LEN(j),seven,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
                BLASNAME(axpy) (
                    COL_LEN(j),z89,
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(15,9));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_blas = 0.;
            for (int k=0; k<nloops2; ++k) e13_blas += NormSq(C3[k]-C0[k]);
            e13_blas = sqrt(e13_blas/nloops2);
            e13_blas /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = 7 * EPART(A4[k]) + T(8,9) * EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e13_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e13_eigen = sqrt(e13_eigen/nloops2);
            e13_eigen /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = 7 * EPART(A5[k]) + T(8,9) * EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e13_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e13_smalleigen = sqrt(e13_smalleigen/nloops2);
            e13_smalleigen /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = (7,1)*A + (8,9)*B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = T(7,1)*A0[k](i,j) + T(8,9) * B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = T(15,10);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = T(7,1) * MPART(A1[k]) + T(8,9) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_reg = 0.;
            for (int k=0; k<nloops2; ++k) e14_reg += NormSq(C1[k]-C0[k]);
            e14_reg = sqrt(e14_reg/nloops2);
            e14_reg /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = T(7,1) * MPART(A2[k]) + T(8,9) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_small = 0.;
            for (int k=0; k<nloops2; ++k) e14_small += NormSq(C2[k]-C0[k]);
            e14_small = sqrt(e14_small/nloops2);
            e14_small /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(scal) (M*N,z71,BP(C3[k].ptr()),1);
#endif
            for(int j=0;j<N;++j) {
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(scal) (
                    COL_LEN(j),z71,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
                BLASNAME(axpy) (
                    COL_LEN(j),z89,
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(15,10));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_blas = 0.;
            for (int k=0; k<nloops2; ++k) e14_blas += NormSq(C3[k]-C0[k]);
            e14_blas = sqrt(e14_blas/nloops2);
            e14_blas /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = T(7,1) * EPART(A4[k]) + T(8,9) * EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e14_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e14_eigen = sqrt(e14_eigen/nloops2);
            e14_eigen /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = T(7,1) * EPART(A5[k]) + T(8,9) * EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e14_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e14_smalleigen = sqrt(e14_smalleigen/nloops2);
            e14_smalleigen /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif
#endif

#if 1 // C = (7,1)*A - B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        C0[k](i,j) = T(7,1)*A0[k](i,j) - B0[k](j,i);
                    else if (UNITPART(i,j))
                        C0[k](i,j) = T(6,1);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C1[k]) = T(7,1) * MPART(A1[k]) - MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_reg = 0.;
            for (int k=0; k<nloops2; ++k) e15_reg += NormSq(C1[k]-C0[k]);
            e15_reg = sqrt(e15_reg/nloops2);
            e15_reg /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(C2[k]) = T(7,1) * MPART(A2[k]) - MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_small = 0.;
            for (int k=0; k<nloops2; ++k) e15_small += NormSq(C2[k]-C0[k]);
            e15_small = sqrt(e15_small/nloops2);
            e15_small /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(C3[k].ptr()),1);
            BLASNAME(scal) (M*N,z71,BP(C3[k].ptr()),1);
#endif
            for(int j=0;j<N;++j) {
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(scal) (
                    COL_LEN(j),z71,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
                BLASNAME(axpy) (
                    COL_LEN(j),mone,
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(C3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#if ((PART == 3) || (PART == 5))
            C3[k].diag().setAllTo(T(6,1));
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_blas = 0.;
            for (int k=0; k<nloops2; ++k) e15_blas += NormSq(C3[k]-C0[k]);
            e15_blas = sqrt(e15_blas/nloops2);
            e15_blas /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(C4[k]) = T(7,1) * EPART(A4[k]) - EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e15_eigen += tmv::TMV_NORM(C4[k](i,j)-C0[k](i,j));
            e15_eigen = sqrt(e15_eigen/nloops2);
            e15_eigen /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(C5[k]) = T(7,1) * EPART(A5[k]) - EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e15_smalleigen += tmv::TMV_NORM(C5[k](i,j)-C0[k](i,j));
            e15_smalleigen = sqrt(e15_smalleigen/nloops2);
            e15_smalleigen /= RT(7)*Norm(C0[0]) + RT(12)*Norm(B0[0]);
        }
#endif
#endif
#endif
#endif
    }

    std::cout<<"C += A                  "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"C += B                  "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"C -= B                  "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"C += 8*B                "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"C = A + Ax              "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
    std::cout<<"C = A + B               "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"C = B + Bx              "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"C = A - B               "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
    std::cout<<"C = 7*A - 12*B          "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
    std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
    std::cout<<"C = 7*C - 12*A          "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
    std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C += (8,9)*A            "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
    std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
    std::cout<<"C += (8,9)*B            "<<t12_reg<<"  "<<t12_small<<"  "<<t12_blas;
    std::cout<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;
    std::cout<<"C = 7*A + (8,9)*B       "<<t13_reg<<"  "<<t13_small<<"  "<<t13_blas;
    std::cout<<"  "<<t13_eigen<<"  "<<t13_smalleigen<<std::endl;
    std::cout<<"C = (7,1)*A + (8,9)*B   "<<t14_reg<<"  "<<t14_small<<"  "<<t14_blas;
    std::cout<<"  "<<t14_eigen<<"  "<<t14_smalleigen<<std::endl;
    std::cout<<"C = (7,1)*A - B         "<<t15_reg<<"  "<<t15_small<<"  "<<t15_blas;
    std::cout<<"  "<<t15_eigen<<"  "<<t15_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"C += A                  "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"C += B                  "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"C -= B                  "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"C += 8*B                "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"C = A + Ax              "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
    std::cout<<"C = A + B               "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
    std::cout<<"C = B + Bx              "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
    std::cout<<"C = A - B               "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
    std::cout<<"C = 7*A - 12*B          "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
    std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
    std::cout<<"C = 7*C - 12*A          "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
    std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"C += (8,9)*A            "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
    std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
    std::cout<<"C += (8,9)*B            "<<e12_reg<<"  "<<e12_small<<"  "<<e12_blas;
    std::cout<<"  "<<e12_eigen<<"  "<<e12_smalleigen<<std::endl;
    std::cout<<"C = 7*A + (8,9)*B       "<<e13_reg<<"  "<<e13_small<<"  "<<e13_blas;
    std::cout<<"  "<<e13_eigen<<"  "<<e13_smalleigen<<std::endl;
    std::cout<<"C = (7,1)*A + (8,9)*B   "<<e14_reg<<"  "<<e14_small<<"  "<<e14_blas;
    std::cout<<"  "<<e14_eigen<<"  "<<e14_smalleigen<<std::endl;
    std::cout<<"C = (7,1)*A - B         "<<e15_reg<<"  "<<e15_small<<"  "<<e15_blas;
    std::cout<<"  "<<e15_eigen<<"  "<<e15_smalleigen<<std::endl;
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
static void NormM(
    const std::vector<tmv::Matrix<T> >& A1)
{
#ifdef ERRORCHECK
    std::vector<tmv::Matrix<T> > A0 = A1;
    std::vector<T> D0(nloops2);
    std::vector<RT> F0(nloops2);
#endif

    std::vector<T> D1(nloops2);
    std::vector<RT> F1(nloops2);

#ifdef DOSMALL
    std::vector<tmv::SmallMatrix<T,M,N> > A2(nloops2);
    std::vector<T> D2(nloops2);
    std::vector<RT> F2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Matrix<T> > A3 = A1;
    std::vector<T> D3(nloops2);
    std::vector<RT> F3(nloops2);
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENM> A4(nloops2,Eigen::EIGENM(M,N));
    std::vector<T> D4(nloops2);
    std::vector<RT> F4(nloops2);
    for(int k=0;k<nloops2;++k) {
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) A4[k](i,j) = A1[k](i,j);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,M,N> > A5;
    if (nloops2x) { A5.resize(nloops2x); }
    std::vector<T> D5(nloops2x);
    std::vector<RT> F5(nloops2x);
    for(int k=0;k<nloops2x;++k) {
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) A5[k](i,j) = A1[k](i,j);
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

    for (int n=0; n<nloops1; ++n) {

#if 1 // F = Norm1(A)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                F0[k] = RT(0);
                for(int j=0;j<N;++j) {
                    RT sum(0);
                    for(int i=0;i<M;++i) {
                        if (INPART(i,j))
                            sum += std::abs(A0[k](i,j));
                        else if (UNITPART(i,j))
                            sum += RT(1);
                    }
                    if (sum > F0[k]) F0[k] = sum;
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = Norm1(MPART(A1[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= F0[0];
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = Norm1(MPART(A2[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= F0[0];
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;


        for (int k=0; k<nloops2; ++k) {
#ifdef XLAP
            char c = '1';
#if (PART == 1)
            F3[k] = LAPNAME(lange) (LAPCM c, M, N, LP(A3[k].cptr()), M, LAPRWORK);
#elif ((PART >= 2) && (PART <= 5))
            char u = UPLO_CHAR;
            char d = DIAG_CHAR;
            F3[k] = LAPNAME(lantr) (LAPCM c, u, d, M, N, LP(A3[k].cptr()), M, LAPRWORK);
#endif
#else
            RT max(0);
            for(int j=0;j<N;j++) {
                RT temp;
#ifdef TISCOMPLEX
                // BLAS doesn't have this capability
#if ((PART >= 1) && (PART <= 5))
                temp = real(std::accumulate(
                        A3[k].col(j,COL_START(j),COL_END(j)).begin(),
                        A3[k].col(j,COL_START(j),COL_END(j)).end(),
                        T(0.),Summer()));
#endif
#else
#if ((PART >= 1) && (PART <= 5))
                temp = std::abs(BLASNAME(asum) (
                        COL_LEN(j),
                        BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1));
#endif
#endif
#if ((PART == 3) || (PART == 5))
                temp += RT(1);
#endif
                if (temp > max) max = temp;
            }
            F3[k] = max;
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= F0[0];
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // Eigen doesn't have this capability
            RT max(0);
            for(int j=0;j<N;j++) {
                RT temp;
                temp = EPART(A4[k]).col(j).cwise().abs().sum();
                if (temp > max) max = temp;
            }
            F4[k] = max;
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e1_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= F0[0];
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            RT max(0);
            for(int j=0;j<N;j++) {
                RT temp;
                temp = EPART(A5[k]).col(j).cwise().abs().sum();
                if (temp > max) max = temp;
            }
            F5[k] = max;
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e1_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= F0[0];
        }
#endif
#endif
#endif

#if 1 // F = NormF(A)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                F0[k] = RT(0);
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        F0[k] += std::abs(A0[k](i,j)*A0[k](i,j));
                    else if (UNITPART(i,j))
                        F0[k] += RT(1);
                F0[k] = std::sqrt(F0[k]);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = NormF(MPART(A1[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= F0[0];
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = NormF(MPART(A2[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= F0[0];
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef XLAP
            char c = 'F';
#if (PART == 1)
            F3[k] = LAPNAME(lange) (LAPCM c, M, N, LP(A3[k].cptr()), M, LAPRWORK);
#elif ((PART >= 2) && (PART <= 5))
            char u = UPLO_CHAR;
            char d = DIAG_CHAR;
            F3[k] = LAPNAME(lantr) (LAPCM c, u, d, M, N, LP(A3[k].cptr()), M, LAPRWORK);
#endif
#else
            RT sum(0);
            for(int j=0;j<N;j++) {
#ifdef TISCOMPLEX
                RT temp;
#ifdef CBLAS
                BLASNAME(dotc)_sub (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(&temp));
#else
                BLASNAME(dotc) (
                    BP(&temp),COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1);
#endif
                sum += temp;
#else
                sum += BLASNAME(dot) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1);
#endif
            }
#if ((PART == 3) || (PART == 5))
            sum += RT(N);
#endif
            F3[k] = sqrt(sum);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= F0[0];
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = EPART(A4[k]).norm();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= F0[0];
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = EPART(A5[k]).norm();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e2_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= F0[0];
        }
#endif
#endif
#endif

#if 1 // F = NormInf(A)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                F0[k] = RT(0);
                for(int i=0;i<M;++i) {
                    RT sum(0);
                    for(int j=0;j<N;++j) 
                        if (INPART(i,j))
                            sum += std::abs(A0[k](i,j));
                        else if (UNITPART(i,j))
                            sum += RT(1);
                    if (sum > F0[k]) F0[k] = sum;
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = NormInf(MPART(A1[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= F0[0];
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = NormInf(MPART(A2[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= F0[0];
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef XLAP
            char c = 'I';
#if (PART == 1)
            F3[k] = LAPNAME(lange) (LAPCM c, M, N, LP(A3[k].cptr()), M, LAPRWORK);
#elif ((PART >= 2) && (PART <= 5))
            char u = UPLO_CHAR;
            char d = DIAG_CHAR;
            F3[k] = LAPNAME(lantr) (LAPCM c, u, d, M, N, LP(A3[k].cptr()), M, LAPRWORK);
#endif
#else
            RT max(0);
            for(int i=0;i<M;i++) {
                RT temp;
#ifdef TISCOMPLEX
                // BLAS doesn't have this capability
                temp = real(std::accumulate(
                        A3[k].row(i,ROW_START(i),ROW_END(i)).begin(),
                        A3[k].row(i,ROW_START(i),ROW_END(i)).end(),
                        T(0.),Summer()));
#else
                temp = std::abs(BLASNAME(asum) (
                        ROW_LEN(i),
                        A3[k].row(i,ROW_START(i),ROW_END(i)).cptr(),M));
#endif
#if ((PART == 3) || (PART == 5))
                temp += RT(1);
#endif
                if (temp > max) max = temp;
            }
            F3[k] = max;
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= F0[0];
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            RT max(0);
            for(int i=0;i<M;i++) {
                RT temp;
                temp = EPART(A4[k]).row(i).cwise().abs().sum();
                if (temp > max) max = temp;
            }
            F4[k] = max;
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e3_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= F0[0];
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            RT max(0);
            for(int i=0;i<M;i++) {
                RT temp;
                temp = EPART(A5[k]).row(i).cwise().abs().sum();
                if (temp > max) max = temp;
            }
            F5[k] = max;
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e3_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= F0[0];
        }
#endif
#endif
#endif

#if 1 // F = NormSq(A)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                F0[k] = RT(0);
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        F0[k] += std::abs(A0[k](i,j)*A0[k](i,j));
                    else if (UNITPART(i,j))
                        F0[k] += RT(1);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = NormSq(MPART(A1[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= F0[0];
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = NormSq(MPART(A2[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= F0[0];
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
#ifdef TISCOMPLEX
            T temp;
#ifdef CBLAS
            BLASNAME(dotc)_sub (N,BP(A3[k].cptr()),1,BP(A3[k].cptr()),1,BP(&temp));
#else
            BLASNAME(dotc) (BP(&temp),M*N,BP(A3[k].cptr()),1,BP(A3[k].cptr()),1);
#endif
            F3[k] = std::real(temp);
#else
            F3[k] = BLASNAME(dot) (M*N,BP(A3[k].cptr()),1,BP(A3[k].cptr()),1);
#endif
#elif ((PART >= 2) && (PART <= 5))
            F3[k] = RT(0);
            for(int j=0; j<N; ++j) {
#ifdef TISCOMPLEX
                T temp;
#ifdef CBLAS
                BLASNAME(dotc)_sub (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(&temp));
#else
                BLASNAME(dotc) (
                    BP(&temp),COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1);
#endif
                F3[k] += std::real(temp);
#else
                F3[k] += BLASNAME(dot) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1);
#endif
            }
#if ((PART == 3) || (PART == 5))
            F3[k] += RT(N);
#endif
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= F0[0];
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = EPART(A4[k]).squaredNorm();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e4_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= F0[0];
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = EPART(A5[k]).squaredNorm();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e4_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= F0[0];
        }
#endif
#endif
#endif

#if 1 // F = MaxAbsElement(A)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                F0[k] = RT(0);
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        if (std::abs(A0[k](i,j)) > F0[k]) F0[k] = std::abs(A0[k](i,j));
                        else if (UNITPART(i,j))
                            if (RT(1) > F0[k]) F0[k] = RT(1);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = MaxAbsElement(MPART(A1[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e5_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= F0[0];
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = MaxAbsElement(MPART(A2[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e5_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= F0[0];
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef XLAP
            char c = 'M';
#if (PART == 1)
            F3[k] = LAPNAME(lange) (LAPCM c, M, N, LP(A3[k].cptr()), M, LAPRWORK);
#elif ((PART >= 2) && (PART <= 5))
            char u = UPLO_CHAR;
            char d = DIAG_CHAR;
            F3[k] = LAPNAME(lantr) (LAPCM c, u, d, M, N, LP(A3[k].cptr()), M, LAPRWORK);
#endif
#else
#if (PART == 1)
#ifdef TISCOMPLEX
            // BLAS doesn't have this capability
            tmv::Vector<T>::iterator maxi = 
                std::max_element(
                    A3[k].linearView().begin(),A3[k].linearView().end(),
                    tmv::Compare<tmv::ASCEND,tmv::ABS_COMP,T>());
            F3[k] = std::abs(*maxi);
#else
            F3[k] = std::abs(A3[k](BLASINAME(amax) (M*N,BP(A3[k].cptr()),1) BLASM1));
#endif
#elif ((PART >= 2) && (PART <= 5))
            RT max(0);
            for(int j=0;j<N;j++) {
                RT temp;
#ifdef TISCOMPLEX
                tmv::Vector<T>::iterator maxi = 
                    std::max_element(
                        A3[k].col(j,COL_START(j),COL_END(j)).begin(),
                        A3[k].col(j,COL_START(j),COL_END(j)).end(),
                        tmv::Compare<tmv::ASCEND,tmv::ABS_COMP,T>());
                temp = std::abs(*maxi);
#else
                temp = std::abs(BLASNAME(asum) (
                        COL_LEN(j),
                        BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1));
#endif
#if ((PART == 3) || (PART == 5))
                if (temp < RT(1)) temp = RT(1);
#endif
                if (temp > max) max = temp;
            }
            F3[k] = max;
#endif
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e5_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= F0[0];
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = EPART(A4[k]).cwise().abs().maxCoeff();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e5_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= F0[0];
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = EPART(A5[k]).cwise().abs().maxCoeff();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e5_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= F0[0];
        }
#endif
#endif
#endif 

#if 1 // D = SumElements(A)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                D0[k] = T(0);
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        D0[k] += A0[k](i,j);
                    else if (UNITPART(i,j))
                        D0[k] += T(1);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = SumElements(MPART(A1[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e6_reg += std::abs((D1[k]-D0[k])*(D1[k]-D0[k]));
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= std::abs(D0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = SumElements(MPART(A2[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e6_small += std::abs((D2[k]-D0[k])*(D2[k]-D0[k]));
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= std::abs(D0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // No BLAS function for this.  Use STL.
#if (PART == 1)
            D3[k] = std::accumulate(
                A3[k].linearView().begin(),A3[k].linearView().end(),T(0.));
#elif ((PART >= 2) && (PART <= 5))
            D3[k] = T(0);
            for(int j=0;j<N;j++) {
                D3[k] += std::accumulate(
                    A3[k].col(j,COL_START(j),COL_END(j)).begin(),
                    A3[k].col(j,COL_START(j),COL_END(j)).end(),T(0.));
            }
#if ((PART == 3) || (PART == 5))
            D3[k] += T(N);
#endif
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e6_blas += std::abs((D3[k]-D0[k])*(D3[k]-D0[k]));
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= std::abs(D0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = EPART(A4[k]).sum();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e6_eigen += std::abs((D4[k]-D0[k])*(D4[k]-D0[k]));
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= std::abs(D0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = EPART(A5[k]).sum();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e6_smalleigen += std::abs((D5[k]-D0[k])*(D5[k]-D0[k]));
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= std::abs(D0[0]);
        }
#endif
#endif
#endif

#if 1 // F = SumAbsElements(A)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                F0[k] = RT(0);
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        F0[k] += std::abs(A0[k](i,j));
                    else if (UNITPART(i,j))
                        F0[k] += RT(1);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = SumAbsElements(MPART(A1[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e7_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e7_reg = sqrt(e7_reg/nloops2);
            e7_reg /= F0[0];
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = SumAbsElements(MPART(A2[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e7_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= F0[0];
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
#ifdef TISCOMPLEX
            // BLAS doesn't have this capability
            F3[k] = real(std::accumulate(
                    A3[k].linearView().begin(),A3[k].linearView().end(),
                    T(0.),Summer()));
#else
            F3[k] = std::abs(BLASNAME(asum) (M*N,BP(A3[k].cptr()),1));
#endif
#elif ((PART >= 2) && (PART <= 5))
            F3[k] = RT(0);
            for(int j=0;j<N;++j) {
#ifdef TISCOMPLEX
                // BLAS doesn't have this capability
                F3[k] += real(std::accumulate(
                        A3[k].col(j,COL_START(j),COL_END(j)).begin(),
                        A3[k].col(j,COL_START(j),COL_END(j)).end(),
                        T(0.),Summer()));
#else
                F3[k] += std::abs(BLASNAME(asum) (
                        COL_LEN(j),
                        BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1));
#endif
            }
#if ((PART == 3) || (PART == 5))
            F3[k] += RT(N);
#endif
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e7_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= F0[0];
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = EPART(A4[k]).cwise().abs().sum();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e7_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= F0[0];
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = EPART(A5[k]).cwise().abs().sum();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e7_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= F0[0];
        }
#endif
#endif
#endif

#if 1 // F = SumAbs2Elements(A)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                F0[k] = RT(0);
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                    {
#ifdef TISCOMPLEX
                        RT temp = std::abs(real(A0[k](i,j))) + std::abs(imag(A0[k](i,j)));
#else
                        RT temp = std::abs(A0[k](i,j));
#endif
                        F0[k] += temp;
                    }
                    else if (UNITPART(i,j))
                        F0[k] += RT(1);
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F1[k] = SumAbs2Elements(MPART(A1[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e8_reg += std::abs((F1[k]-F0[k])*(F1[k]-F0[k]));
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= F0[0];
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F2[k] = SumAbs2Elements(MPART(A2[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) 
                e8_small += std::abs((F2[k]-F0[k])*(F2[k]-F0[k]));
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= F0[0];
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            F3[k] = std::abs(BLASDZNAME(asum) (M*N,BP(A3[k].cptr()),1));
#elif ((PART >= 2) && (PART <= 5))
            F3[k] = RT(0);
            for(int j=0;j<N;++j) {
                F3[k] += std::abs(BLASDZNAME(asum) (
                        COL_LEN(j),
                        BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1));
            }
#if ((PART == 3) || (PART == 5))
            F3[k] += RT(N);
#endif
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) 
                e8_blas += std::abs((F3[k]-F0[k])*(F3[k]-F0[k]));
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= F0[0];
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            F4[k] = (
                EPART(A4[k]).real().cwise().abs().sum() + 
                EPART(A4[k]).imag().cwise().abs().sum());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                e8_eigen += std::abs((F4[k]-F0[k])*(F4[k]-F0[k]));
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= F0[0];
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            F5[k] = (
                EPART(A5[k]).real().cwise().abs().sum() + 
                EPART(A5[k]).imag().cwise().abs().sum());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                e8_smalleigen += std::abs((F5[k]-F0[k])*(F5[k]-F0[k]));
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= F0[0];
        }
#endif
#endif
#endif
    }

    std::cout<<"Norm1(A)                "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"NormF(A)                "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"NormInf(A)              "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"NormSq(A)               "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"MaxAbsElement(A)        "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
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
    std::cout<<"NormF(A)                "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"NormInf(A)              "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"NormSq(A)               "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"MaxAbsElement(A)        "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
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

#ifdef DOSWAP
static void SwapM(
    std::vector<tmv::Matrix<T> >& A1,
    std::vector<tmv::Matrix<T> >& B1)
{
    // Make plausible random permutation2:
    tmv::Vector<T> A_temp1 = A1[0].col(0);
    int P[M];
    A_temp1.sort(P);
    tmv::Vector<T> A_temp2 = A1[0].row(0);
    int Q[N];
    A_temp2.sort(Q);

    std::vector<tmv::Matrix<T> > A1x = A1;
    std::vector<tmv::Matrix<T> > B1x = B1;
    for(int k=0; k<nloops2; ++k) {
        A1x[k] += B1[k].transpose(); 
        B1x[k] -= A1[k].transpose();
    }

#ifdef ERRORCHECK
    std::vector<tmv::Matrix<T> > A0 = A1;
    std::vector<tmv::Matrix<T> > B0 = B1;
    std::vector<tmv::Matrix<T> > A0x = A1x;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallMatrix<T,M,N> > A2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,M> > B2(nloops2);
    std::vector<tmv::SmallMatrix<T,M,N> > A2x(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; B2[k] = B1[k]; A2x[k] = A1x[k]; 
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Matrix<T> > A3 = A1;
    std::vector<tmv::Matrix<T> > B3 = B1;
    std::vector<tmv::Matrix<T> > A3x = A1x;

    int P3[M];
    int Q3[N];
    for(int i=0;i<M;++i) P3[i] = P[i] + 1;
    for(int j=0;j<N;++j) Q3[j] = Q[j] + 1;
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENM> A4(nloops2,Eigen::EIGENM(M,N));
    std::vector<Eigen::EIGENM> B4(nloops2,Eigen::EIGENM(N,M));
    std::vector<Eigen::EIGENM> A4x(nloops2,Eigen::EIGENM(M,N));
    for(int k=0;k<nloops2;++k) {
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
            A4[k](i,j) = A1[k](i,j);
            B4[k](j,i) = B1[k](j,i);
            A4x[k](i,j) = A1x[k](i,j);
        }
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,M,N> > A5;
    std::vector<Eigen::Matrix<T,N,M> > B5;
    std::vector<Eigen::Matrix<T,M,N> > A5x;
    if (nloops2x) { 
        A5.resize(nloops2x); B5.resize(nloops2x); A5x.resize(nloops2x); 
    }
    for(int k=0;k<nloops2x;++k) {
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
            A5[k](i,j) = A1[k](i,j);
            B5[k](j,i) = B1[k](j,i);
            A5x[k](i,j) = A1x[k](i,j);
        }
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

    for (int n=0; n<nloops1; ++n) {

#if 1 // Swap(A,Ax);
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        std::swap(A0[k](i,j),A0x[k](i,j));
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) 
            Swap(MPART(A1[k]),MPART(A1x[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(A1x[k]-A0x[k]);
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(A1[k]-A0[k]);
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            Swap(MPART(A2[k]),MPART(A2x[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(A2x[k]-A0x[k]);
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(A2[k]-A0[k]);
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(swap) (M*N,BP(A3[k].ptr()),1,BP(A3x[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(swap) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).ptr()),1,
                    BP(A3x[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(A3x[k]-A0x[k]);
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(A3[k]-A0[k]);
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART3(A4[k]).swap(EPART3(A4x[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e1_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e1_eigen += tmv::TMV_NORM(A4x[k](i,j)-A0x[k](i,j));
            }
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART3(A5[k]).swap(EPART3(A5x[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e1_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e1_smalleigen += tmv::TMV_NORM(A5x[k](i,j)-A0x[k](i,j));
            }
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif
#endif

#if 1 // Swap(A[0::N/2,0::N/2],Ax[0::N/2,0::N/2])
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M/2;++i) for(int j=0;j<N/2;++j) 
                    if (INPART(i,j))
                        std::swap(A0[k](i,j),A0x[k](i,j));
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            Swap(
                MPART(A1[k].subMatrix(0,M/2,0,N/2)),
                MPART(A1x[k].subMatrix(0,M/2,0,N/2)));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) e2_reg += NormSq(A1[k]-A0[k]);
            for (int k=0; k<nloops2; ++k) e2_reg += NormSq(A1x[k]-A0x[k]);
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            Swap(
                MPART(A2[k].subMatrix(0,M/2,0,N/2)),
                MPART(A2x[k].subMatrix(0,M/2,0,N/2)));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) e2_small += NormSq(A2[k]-A0[k]);
            for (int k=0; k<nloops2; ++k) e2_small += NormSq(A2x[k]-A0x[k]);
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            for(int j=0;j<N/2;++j) 
                BLASNAME(swap) (M/2,BP(A3[k].ptr()+j*M),1,BP(A3x[k].ptr()+j*M),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N/2;++j) {
                int col_end = M/2;
                if (COL_END(j) < col_end) col_end = COL_END(j);
                int col_len = col_end - COL_START(j);
                BLASNAME(swap) (
                    col_len,
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).ptr()),1,
                    BP(A3x[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) e2_blas += NormSq(A3[k]-A0[k]);
            for (int k=0; k<nloops2; ++k) e2_blas += NormSq(A3x[k]-A0x[k]);
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART3(A4[k].block(0,0,M/2,N/2)).swap(EPART3(A4x[k].block(0,0,M/2,N/2)));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e2_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e2_eigen += tmv::TMV_NORM(A4x[k](i,j)-A0x[k](i,j));
            }
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        // ??? Eigen has a bug where it gets this wrong for N == 1
        if (N > 1) 
            for (int k=0; k<nloops2x; ++k)
                EPART3(A5[k].block(0,0,M/2,N/2)).swap(EPART3(A5x[k].block(0,0,M/2,N/2)));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e2_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e2_smalleigen += tmv::TMV_NORM(A5x[k](i,j)-A0x[k](i,j));
            }
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif
#endif

#if 1 // Swap(A,B);
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        std::swap(A0[k](i,j),B0[k](j,i));
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) 
            Swap(MPART(A1[k]),MPART(B1[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) e3_reg += NormSq(B1[k]-B0[k]);
            for (int k=0; k<nloops2; ++k) e3_reg += NormSq(A1[k]-A0[k]);
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            Swap(MPART(A2[k]),MPART(B2[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) e3_small += NormSq(B2[k]-B0[k]);
            for (int k=0; k<nloops2; ++k) e3_small += NormSq(A2[k]-A0[k]);
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(swap) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).ptr()),1,
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).ptr()),N);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) e3_blas += NormSq(B3[k]-B0[k]);
            for (int k=0; k<nloops2; ++k) e3_blas += NormSq(A3[k]-A0[k]);
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // Eigen requires doing this with block else it won't compile.
            //A4[k].swap(EPART(B4[k].transpose()));
            EPART3(A4[k].block(0,0,M,N)).swap(EPART3(B4[k].transpose()));
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e3_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e3_eigen += tmv::TMV_NORM(B4[k](j,i)-B0[k](j,i));
            }
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART3(A5[k].block(0,0,M,N)).swap(EPART3(B5[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e3_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e3_smalleigen += tmv::TMV_NORM(B5[k](j,i)-B0[k](j,i));
            }
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(A0[0]) + Norm(B0[0]);
        }
#endif
#endif
#endif

#ifdef DO_TRANSPOSE_SELF
#if 1 // A.transposeSelf()
#ifdef AISSQUARE
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<i;++j) 
                    if (INPART(i,j) && INPART(j,i))
                        std::swap(A0[k](i,j),A0[k](j,i));
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]).transposeSelf();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            MPART(A2[k]).transposeSelf();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
#if (PART == 1)
            for(int i=1;i<M;++i)
                BLASNAME(swap) (i,BP(A3[k].ptr()+i),M,BP(A3[k].ptr()+i*M),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            EPART(A4[k]).transposeInPlace();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e4_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) 
            EPART(A5[k]).transposeInPlace();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e4_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif
#endif
#endif

#ifdef TISCOMPLEX
#if 1 // A.conjugateSelf()
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                        A0[k](i,j) = std::conj(A0[k](i,j));
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]).conjugateSelf();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            MPART(A2[k]).conjugateSelf();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
#if (PART == 1)
            BLASDNAME(scal) (M*N,dmone,A3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASDNAME(scal) (
                    COL_LEN(j),dmone,
                    A3[k].col(j,COL_START(j),COL_END(j)).realPart().ptr()+1,2);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            EPART3(A4[k]) = EPART(A4[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e5_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            }
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART3(A5[k]) = EPART(A5[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e5_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            }
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // Swap(A,Ax.conjugate);
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) 
            for (int k=0; k<nloops2; ++k) 
                for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                    if (INPART(i,j))
                    {
                        T temp = conj(A0x[k](i,j));
                        A0x[k](i,j) = conj(A0[k](i,j));
                        A0[k](i,j) = temp;
                    }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) 
            Swap(MPART(A1[k]),MPART(A1x[k].conjugate()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) e6_reg += NormSq(A1x[k]-A0x[k]);
            for (int k=0; k<nloops2; ++k) e6_reg += NormSq(A1[k]-A0[k]);
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            Swap(MPART(A2[k]),MPART(A2x[k]).conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) e6_small += NormSq(A2x[k]-A0x[k]);
            for (int k=0; k<nloops2; ++k) e6_small += NormSq(A2[k]-A0[k]);
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(swap) (M*N,BP(A3[k].ptr()),1,BP(A3x[k].ptr()),1);
            BLASDNAME(scal) (M*N,dmone,A3[k].realPart().ptr()+1,2);
            BLASDNAME(scal) (M*N,dmone,A3x[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(swap) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).ptr()),1,
                    BP(A3x[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASDNAME(scal) (
                    COL_LEN(j),dmone,
                    A3[k].col(j,COL_START(j),COL_END(j)).realPart().ptr()+1,2);
                BLASDNAME(scal) (
                    COL_LEN(j),dmone,
                    A3x[k].col(j,COL_START(j),COL_END(j)).realPart().ptr()+1,2);
            }
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) e6_blas += NormSq(A3x[k]-A0x[k]);
            for (int k=0; k<nloops2; ++k) e6_blas += NormSq(A3[k]-A0[k]);
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            // Neither of these work.
            //A4[k].swap(A4x[k].conjugate());
            //A4[k].conjugate().swap(A4x[k]);
            EPART3(A4[k]).swap(EPART3(A4x[k]));
            EPART3(A4[k]) = EPART(A4[k]).conjugate();
            EPART3(A4x[k]) = EPART(A4x[k].conjugate());
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) {
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e6_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e6_eigen += tmv::TMV_NORM(A4x[k](i,j)-A0x[k](i,j));
            }
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) {
            EPART3(A5[k]).swap(EPART3(A5x[k]));
            EPART3(A5[k]) = EPART(A5[k]).conjugate();
            EPART3(A5x[k]) = EPART(A5x[k]).conjugate();
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) {
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e6_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e6_smalleigen += tmv::TMV_NORM(A5x[k](i,j)-A0x[k](i,j));
            }
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(A0[0]) + Norm(A0x[0]);
        }
#endif
#endif
#endif
#endif

#ifdef DO_PERMUTE
#if 1 // A.permuteRows(P)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) 
                    for(int j=0;j<N;++j) 
                        std::swap(A0[k](i,j),A0[k](P[i],j));
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]).permuteRows(P);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            MPART(A2[k]).permuteRows(P);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
#if (PART == 1)
#ifdef XLAP
            LAPNAME(laswp) (N,LP(A3[k].ptr()),M,1,M,P3,1);
#else
            for(int i=0;i<M;++i)
                if (P[i] != i)
                    BLASNAME(swap) (M,BP(A3[k].ptr()+i),M,BP(A3[k].ptr()+P[i]),M);
#endif
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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

        for (int k=0; k<nloops2; ++k) 
            for(int j=0;j<N;++j) 
                for(int i=0;i<M;++i) 
                    std::swap(A4[k](i,j),A4[k](P[i],j));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e7_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) 
            for(int j=0;j<N;++j) 
                for(int i=0;i<M;++i) 
                    std::swap(A5[k](i,j),A5[k](P[i],j));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e7_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // A.reversePermuteRows(P)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=M-1;i>=0;--i) 
                    for(int j=0;j<N;++j)
                        std::swap(A0[k](i,j),A0[k](P[i],j));
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]).reversePermuteRows(P);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
            MPART(A2[k]).reversePermuteRows(P);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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
#if (PART == 1)
#ifdef XLAP
            LAPNAME(laswp) (N,LP(A3[k].ptr()),M,1,M,P3,-1);
#else
            for(int i=M-1;i>=0;--i)
                if (P[i] != i)
                    BLASNAME(swap) (M,BP(A3[k].ptr()+i),M,BP(A3[k].ptr()+P[i]),M);
#endif
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
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

        for (int k=0; k<nloops2; ++k) 
            for(int j=0;j<N;++j)
                for(int i=M-1;i>=0;--i) 
                    std::swap(A4[k](i,j),A4[k](P[i],j));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e8_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) 
            for(int j=0;j<N;++j)
                for(int i=M-1;i>=0;--i) 
                    std::swap(A5[k](i,j),A5[k](P[i],j));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e8_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // A.permuteCols(Q)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) 
                    for(int j=0;j<N;++j) 
                        std::swap(A0[k](i,j),A0[k](i,Q[j]));
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]).permuteCols(Q);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_reg = 0.;
            for (int k=0; k<nloops2; ++k) e9_reg += NormSq(A1[k]-A0[k]);
            e9_reg = sqrt(e9_reg/nloops2);
            e9_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A2[k]).permuteCols(Q);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_small = 0.;
            for (int k=0; k<nloops2; ++k) e9_small += NormSq(A2[k]-A0[k]);
            e9_small = sqrt(e9_small/nloops2);
            e9_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            for(int j=0;j<N;++j)
                if (Q[j] != j)
                    BLASNAME(swap) (M,BP(A3[k].ptr()+j*M),1,BP(A3[k].ptr()+Q[j]*M),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_blas = 0.;
            for (int k=0; k<nloops2; ++k) e9_blas += NormSq(A3[k]-A0[k]);
            e9_blas = sqrt(e9_blas/nloops2);
            e9_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) 
            for(int j=0;j<N;++j) 
                for(int i=0;i<M;++i) 
                    std::swap(A4[k](i,j),A4[k](i,Q[j]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e9_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e9_eigen = sqrt(e9_eigen/nloops2);
            e9_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) 
            for(int j=0;j<N;++j) 
                for(int i=0;i<M;++i) 
                    std::swap(A5[k](i,j),A5[k](i,Q[j]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e9_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e9_smalleigen = sqrt(e9_smalleigen/nloops2);
            e9_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif

#if 1 // A.reversePermuteCols(Q)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) 
                    for(int j=N-1;j>=0;--j) 
                        std::swap(A0[k](i,j),A0[k](i,Q[j]));
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]).reversePermuteCols(Q);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_reg = 0.;
            for (int k=0; k<nloops2; ++k) e10_reg += NormSq(A1[k]-A0[k]);
            e10_reg = sqrt(e10_reg/nloops2);
            e10_reg /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A2[k]).reversePermuteCols(Q);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_small = 0.;
            for (int k=0; k<nloops2; ++k) e10_small += NormSq(A2[k]-A0[k]);
            e10_small = sqrt(e10_small/nloops2);
            e10_small /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            for(int j=N-1;j>=0;--j)
                if (Q[j] != j)
                    BLASNAME(swap) (M,BP(A3[k].ptr()+j*M),1,BP(A3[k].ptr()+Q[j]*M),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_blas = 0.;
            for (int k=0; k<nloops2; ++k) e10_blas += NormSq(A3[k]-A0[k]);
            e10_blas = sqrt(e10_blas/nloops2);
            e10_blas /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) 
            for(int j=N-1;j>=0;--j)
                for(int i=0;i<M;++i) 
                    std::swap(A4[k](i,j),A4[k](i,Q[j]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e10_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e10_eigen = sqrt(e10_eigen/nloops2);
            e10_eigen /= Norm(A0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k) 
            for(int j=N-1;j>=0;--j)
                for(int i=0;i<M;++i) 
                    std::swap(A5[k](i,j),A5[k](i,Q[j]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_smalleigen = 0.;
            for (int k=0; k<nloops2x; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e10_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e10_smalleigen = sqrt(e10_smalleigen/nloops2);
            e10_smalleigen /= Norm(A0[0]);
        }
#endif
#endif
#endif
#endif
    }

    std::cout<<"Swap(A,Ax)              "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"Swap(A',Ax')            "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"Swap(A,B)               "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
#ifdef DO_TRANSPOSE_SELF
#ifdef AISSQUARE
    std::cout<<"A.transposeSelf()       "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
#endif
#endif
#ifdef TISCOMPLEX
    std::cout<<"A.conjugateSelf()       "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
    std::cout<<"Swap(A,Ax.conjugate())  "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
#endif
#ifdef DO_PERMUTE
    std::cout<<"A.permuteRows(P)        "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"A.reversePermuteRows(P) "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
    std::cout<<"A.permuteCols(Q)        "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
    std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
    std::cout<<"A.reversePermuteCols(Q) "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
    std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"Swap(A,Ax)              "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"Swap(A',Ax')            "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"Swap(A,B)               "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
#ifdef DO_TRANSPOSE_SELF
#ifdef AISSQUARE
    std::cout<<"A.transposeSelf()       "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
#endif
#endif
#ifdef TISCOMPLEX
    std::cout<<"A.conjugateSelf()       "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
    std::cout<<"Swap(A,Ax.conjugate())  "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
#endif
#ifdef DOPERMUTE
    std::cout<<"A.permuteRows(P)        "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
    std::cout<<"A.reversePermuteRows(P) "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
    std::cout<<"A.permuteCols(Q)        "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
    std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
    std::cout<<"A.reversePermuteCols(Q) "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
    std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
#endif
    std::cout<<"\n\n";
#endif
}
#endif

#ifdef DOMULTMV
static void MultMV(
    const std::vector<tmv::Matrix<T> >& A1,
    const std::vector<tmv::Matrix<T> >& B1,
    const std::vector<tmv::Vector<T> >& C1,
    std::vector<tmv::Vector<T> >& D1)
{
#ifdef ERRORCHECK
    std::vector<tmv::Matrix<T> > A0 = A1;
    std::vector<tmv::Matrix<T> > B0 = B1;
    std::vector<tmv::Vector<T> > C0 = C1;
    std::vector<tmv::Vector<T> > D0 = D1;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallMatrix<T,M,N> > A2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,M> > B2(nloops2);
    std::vector<tmv::SmallVector<T,N> > C2(nloops2);
    std::vector<tmv::SmallVector<T,M> > D2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; B2[k] = B1[k]; C2[k] = C1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Matrix<T> > A3 = A1;
    std::vector<tmv::Matrix<T> > B3 = B1;
    std::vector<tmv::Vector<T> > C3 = C1;
    std::vector<tmv::Vector<T> > D3 = D1;
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENM> A4(nloops2,Eigen::EIGENM(M,N));
    std::vector<Eigen::EIGENM> B4(nloops2,Eigen::EIGENM(N,M));
    std::vector<Eigen::EIGENV> C4(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENV> D4(nloops2,Eigen::EIGENV(M));
    for(int k=0;k<nloops2;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) A4[k](i,j) = A1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) B4[k](j,i) = B1[k](j,i);
        for(int i=0;i<N;++i) C4[k](i) = C1[k](i);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,M,N> > A5;
    std::vector<Eigen::Matrix<T,N,M> > B5;
    std::vector<Eigen::Matrix<T,N,1> > C5;
    std::vector<Eigen::Matrix<T,M,1> > D5;
    if (nloops2x) { 
        A5.resize(nloops2x); B5.resize(nloops2x);
        C5.resize(nloops2x); D5.resize(nloops2x);
    }
    for(int k=0;k<nloops2x;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) A5[k](i,j) = A1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) B5[k](j,i) = B1[k](j,i);
        for(int i=0;i<N;++i) C5[k](i) = C1[k](i);
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
    double t17_reg=0., t17_small=0., t17_blas=0., t17_eigen=0., t17_smalleigen=0.;
    double t18_reg=0., t18_small=0., t18_blas=0., t18_eigen=0., t18_smalleigen=0.;
    double t19_reg=0., t19_small=0., t19_blas=0., t19_eigen=0., t19_smalleigen=0.;
    double t20_reg=0., t20_small=0., t20_blas=0., t20_eigen=0., t20_smalleigen=0.;
    double t21_reg=0., t21_small=0., t21_blas=0., t21_eigen=0., t21_smalleigen=0.;
    double t22_reg=0., t22_small=0., t22_blas=0., t22_eigen=0., t22_smalleigen=0.;
    double t23_reg=0., t23_small=0., t23_blas=0., t23_eigen=0., t23_smalleigen=0.;
    double t24_reg=0., t24_small=0., t24_blas=0., t24_eigen=0., t24_smalleigen=0.;
    double t25_reg=0., t25_small=0., t25_blas=0., t25_eigen=0., t25_smalleigen=0.;
    double t26_reg=0., t26_small=0., t26_blas=0., t26_eigen=0., t26_smalleigen=0.;
    double t27_reg=0., t27_small=0., t27_blas=0., t27_eigen=0., t27_smalleigen=0.;
    double t28_reg=0., t28_small=0., t28_blas=0., t28_eigen=0., t28_smalleigen=0.;
    double t29_reg=0., t29_small=0., t29_blas=0., t29_eigen=0., t29_smalleigen=0.;
    double t30_reg=0., t30_small=0., t30_blas=0., t30_eigen=0., t30_smalleigen=0.;
    double t31_reg=0., t31_small=0., t31_blas=0., t31_eigen=0., t31_smalleigen=0.;
    double t32_reg=0., t32_small=0., t32_blas=0., t32_eigen=0., t32_smalleigen=0.;
    double t33_reg=0., t33_small=0., t33_blas=0., t33_eigen=0., t33_smalleigen=0.;
    double t34_reg=0., t34_small=0., t34_blas=0., t34_eigen=0., t34_smalleigen=0.;
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
    double e17_reg=0., e17_small=0., e17_blas=0., e17_eigen=0., e17_smalleigen=0.;
    double e18_reg=0., e18_small=0., e18_blas=0., e18_eigen=0., e18_smalleigen=0.;
    double e19_reg=0., e19_small=0., e19_blas=0., e19_eigen=0., e19_smalleigen=0.;
    double e20_reg=0., e20_small=0., e20_blas=0., e20_eigen=0., e20_smalleigen=0.;
    double e21_reg=0., e21_small=0., e21_blas=0., e21_eigen=0., e21_smalleigen=0.;
    double e22_reg=0., e22_small=0., e22_blas=0., e22_eigen=0., e22_smalleigen=0.;
    double e23_reg=0., e23_small=0., e23_blas=0., e23_eigen=0., e23_smalleigen=0.;
    double e24_reg=0., e24_small=0., e24_blas=0., e24_eigen=0., e24_smalleigen=0.;
    double e25_reg=0., e25_small=0., e25_blas=0., e25_eigen=0., e25_smalleigen=0.;
    double e26_reg=0., e26_small=0., e26_blas=0., e26_eigen=0., e26_smalleigen=0.;
    double e27_reg=0., e27_small=0., e27_blas=0., e27_eigen=0., e27_smalleigen=0.;
    double e28_reg=0., e28_small=0., e28_blas=0., e28_eigen=0., e28_smalleigen=0.;
    double e29_reg=0., e29_small=0., e29_blas=0., e29_eigen=0., e29_smalleigen=0.;
    double e30_reg=0., e30_small=0., e30_blas=0., e30_eigen=0., e30_smalleigen=0.;
    double e31_reg=0., e31_small=0., e31_blas=0., e31_eigen=0., e31_smalleigen=0.;
    double e32_reg=0., e32_small=0., e32_blas=0., e32_eigen=0., e32_smalleigen=0.;
    double e33_reg=0., e33_small=0., e33_blas=0., e33_eigen=0., e33_smalleigen=0.;
    double e34_reg=0., e34_small=0., e34_blas=0., e34_eigen=0., e34_smalleigen=0.;
#endif

    for (int n=0; n<nloops1; ++n) {

#if 1 // D = A * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = MPART(A1[k]) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(D1[k]-D0[k]);
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = MPART(A2[k]) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(D2[k]-D0[k]);
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,one,BP(A3[k].cptr()),M,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASNT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(D3[k]-D0[k]);
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = EPART(A4[k]) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e1_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = EPART(A5[k]) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e1_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = B * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = MPART(B1[k].transpose()) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_reg += NormSq(D1[k]-D0[k]);
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = MPART(B2[k].transpose()) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) e2_small += NormSq(D2[k]-D0[k]);
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLAST,
                            N,M,one,BP(B3[k].cptr()),N,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLAST, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) e2_blas += NormSq(D3[k]-D0[k]);
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = EPART(B4[k].transpose()) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e2_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = EPART(B5[k].transpose()) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e2_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = -A * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) -= A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) -= C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = -MPART(A1[k]) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) e3_reg += NormSq(D1[k]-D0[k]);
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = -MPART(A2[k]) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) e3_small += NormSq(D2[k]-D0[k]);
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,mone,BP(A3[k].cptr()),M,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASNT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,mone,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) e3_blas += NormSq(D3[k]-D0[k]);
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = -EPART(A4[k]) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e3_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = -EPART(A5[k]) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e3_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = -B * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) -= B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) -= C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = -MPART(B1[k].transpose()) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) e4_reg += NormSq(D1[k]-D0[k]);
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = -MPART(B2[k].transpose()) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) e4_small += NormSq(D2[k]-D0[k]);
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLAST,
                            N,M,mone,BP(B3[k].cptr()),N,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLAST, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,mone,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) e4_blas += NormSq(D3[k]-D0[k]);
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = -EPART(B4[k].transpose()) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e4_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = -EPART(B5[k].transpose()) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e4_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = 7 * A * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += RT(7) * A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += RT(7) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = 7 * MPART(A1[k]) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) e5_reg += NormSq(D1[k]-D0[k]);
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = 7 * MPART(A2[k]) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) e5_small += NormSq(D2[k]-D0[k]);
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,seven,BP(A3[k].cptr()),M,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASNT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,seven,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) e5_blas += NormSq(D3[k]-D0[k]);
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = 7 * EPART(A4[k]) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e5_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = 7 * EPART(A5[k]) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e5_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = 7 * B * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += RT(7) * B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += RT(7) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = 7 * MPART(B1[k].transpose()) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) e6_reg += NormSq(D1[k]-D0[k]);
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = 7 * MPART(B2[k].transpose()) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) e6_small += NormSq(D2[k]-D0[k]);
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLAST,
                            N,M,seven,BP(B3[k].cptr()),N,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLAST, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,seven,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) e6_blas += NormSq(D3[k]-D0[k]);
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = 7 * EPART(B4[k].transpose()) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e6_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = 7 * EPART(B5[k].transpose()) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e6_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D -= A * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) -= A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) -= C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] -= MPART(A1[k]) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) e7_reg += NormSq(D1[k]-D0[k]);
            e7_reg = sqrt(e7_reg/nloops2);
            e7_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] -= MPART(A2[k]) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) e7_small += NormSq(D2[k]-D0[k]);
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,mone,BP(A3[k].cptr()),M,BP(C3[k].cptr()),1,
                            one,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASNT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(temp.ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(axpy) (N,mone,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) e7_blas += NormSq(D3[k]-D0[k]);
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] -= EPART(A4[k]) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e7_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] -= EPART(A5[k]) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e7_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D -= B * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) -= B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) -= C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] -= MPART(B1[k].transpose()) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) e8_reg += NormSq(D1[k]-D0[k]);
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] -= MPART(B2[k].transpose()) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) e8_small += NormSq(D2[k]-D0[k]);
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLAST,
                            N,M,mone,BP(B3[k].cptr()),N,BP(C3[k].cptr()),1,
                            one,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLAST, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(temp.ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(axpy) (N,mone,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) e8_blas += NormSq(D3[k]-D0[k]);
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] -= EPART(B4[k].transpose()) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e8_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] -= EPART(B5[k].transpose()) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e8_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += 8 * A * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += RT(8) * A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += RT(8) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] += 8 * MPART(A1[k]) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_reg = 0.;
            for (int k=0; k<nloops2; ++k) e9_reg += NormSq(D1[k]-D0[k]);
            e9_reg = sqrt(e9_reg/nloops2);
            e9_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] += 8 * MPART(A2[k]) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_small = 0.;
            for (int k=0; k<nloops2; ++k) e9_small += NormSq(D2[k]-D0[k]);
            e9_small = sqrt(e9_small/nloops2);
            e9_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,eight,BP(A3[k].cptr()),M,BP(C3[k].cptr()),1,
                            one,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASNT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(temp.ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(axpy) (N,eight,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_blas = 0.;
            for (int k=0; k<nloops2; ++k) e9_blas += NormSq(D3[k]-D0[k]);
            e9_blas = sqrt(e9_blas/nloops2);
            e9_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] += 8 * EPART(A4[k]) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e9_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e9_eigen = sqrt(e9_eigen/nloops2);
            e9_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] += 8 * EPART(A5[k]) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e9_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e9_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e9_smalleigen = sqrt(e9_smalleigen/nloops2);
            e9_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += 8 * B * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += RT(8) * B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += RT(8) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] += 8 * MPART(B1[k].transpose()) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_reg = 0.;
            for (int k=0; k<nloops2; ++k) e10_reg += NormSq(D1[k]-D0[k]);
            e10_reg = sqrt(e10_reg/nloops2);
            e10_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] += 8 * MPART(B2[k].transpose()) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_small = 0.;
            for (int k=0; k<nloops2; ++k) e10_small += NormSq(D2[k]-D0[k]);
            e10_small = sqrt(e10_small/nloops2);
            e10_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLAST,
                            N,M,eight,BP(B3[k].cptr()),N,BP(C3[k].cptr()),1,
                            one,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLAST, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(temp.ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(axpy) (N,eight,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_blas = 0.;
            for (int k=0; k<nloops2; ++k) e10_blas += NormSq(D3[k]-D0[k]);
            e10_blas = sqrt(e10_blas/nloops2);
            e10_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] += 8 * EPART(B4[k].transpose()) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e10_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e10_eigen = sqrt(e10_eigen/nloops2);
            e10_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] += 8 * EPART(B5[k].transpose()) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e10_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e10_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e10_smalleigen = sqrt(e10_smalleigen/nloops2);
            e10_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#ifdef TISCOMPLEX
#if 1 // D = (7,1) * A * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(7,1) * A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = T(7,1) * MPART(A1[k]) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_reg = 0.;
            for (int k=0; k<nloops2; ++k) e11_reg += NormSq(D1[k]-D0[k]);
            e11_reg = sqrt(e11_reg/nloops2);
            e11_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = T(7,1) * MPART(A2[k]) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_small = 0.;
            for (int k=0; k<nloops2; ++k) e11_small += NormSq(D2[k]-D0[k]);
            e11_small = sqrt(e11_small/nloops2);
            e11_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,z71,BP(A3[k].cptr()),M,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASNT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_blas = 0.;
            for (int k=0; k<nloops2; ++k) e11_blas += NormSq(D3[k]-D0[k]);
            e11_blas = sqrt(e11_blas/nloops2);
            e11_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = T(7,1) * EPART(A4[k]) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e11_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e11_eigen = sqrt(e11_eigen/nloops2);
            e11_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = T(7,1) * EPART(A5[k]) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e11_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e11_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e11_smalleigen = sqrt(e11_smalleigen/nloops2);
            e11_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * B * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(7,1) * B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = T(7,1) * MPART(B1[k].transpose()) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_reg = 0.;
            for (int k=0; k<nloops2; ++k) e12_reg += NormSq(D1[k]-D0[k]);
            e12_reg = sqrt(e12_reg/nloops2);
            e12_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = T(7,1) * MPART(B2[k].transpose()) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_small = 0.;
            for (int k=0; k<nloops2; ++k) e12_small += NormSq(D2[k]-D0[k]);
            e12_small = sqrt(e12_small/nloops2);
            e12_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLAST,
                            N,M,z71,BP(B3[k].cptr()),N,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLAST, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_blas = 0.;
            for (int k=0; k<nloops2; ++k) e12_blas += NormSq(D3[k]-D0[k]);
            e12_blas = sqrt(e12_blas/nloops2);
            e12_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = T(7,1) * EPART(B4[k].transpose()) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e12_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e12_eigen = sqrt(e12_eigen/nloops2);
            e12_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = T(7,1) * EPART(B5[k].transpose()) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e12_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e12_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e12_smalleigen = sqrt(e12_smalleigen/nloops2);
            e12_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * A * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(8,9) * A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += T(8,9) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] += T(8,9) * MPART(A1[k]) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_reg = 0.;
            for (int k=0; k<nloops2; ++k) e13_reg += NormSq(D1[k]-D0[k]);
            e13_reg = sqrt(e13_reg/nloops2);
            e13_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] += T(8,9) * MPART(A2[k]) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_small = 0.;
            for (int k=0; k<nloops2; ++k) e13_small += NormSq(D2[k]-D0[k]);
            e13_small = sqrt(e13_small/nloops2);
            e13_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,z89,BP(A3[k].cptr()),M,BP(C3[k].cptr()),1,
                            one,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASNT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(temp.ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(axpy) (N,z89,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_blas = 0.;
            for (int k=0; k<nloops2; ++k) e13_blas += NormSq(D3[k]-D0[k]);
            e13_blas = sqrt(e13_blas/nloops2);
            e13_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] += T(8,9) * EPART(A4[k]) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e13_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e13_eigen = sqrt(e13_eigen/nloops2);
            e13_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] += T(8,9) * EPART(A5[k]) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e13_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e13_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e13_smalleigen = sqrt(e13_smalleigen/nloops2);
            e13_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * B * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(8,9) * B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += T(8,9) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] += T(8,9) * MPART(B1[k].transpose()) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_reg = 0.;
            for (int k=0; k<nloops2; ++k) e14_reg += NormSq(D1[k]-D0[k]);
            e14_reg = sqrt(e14_reg/nloops2);
            e14_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] += T(8,9) * MPART(B2[k].transpose()) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_small = 0.;
            for (int k=0; k<nloops2; ++k) e14_small += NormSq(D2[k]-D0[k]);
            e14_small = sqrt(e14_small/nloops2);
            e14_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLAST,
                            N,M,z89,BP(B3[k].cptr()),N,BP(C3[k].cptr()),1,
                            one,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLAST, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(temp.ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(axpy) (N,z89,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_blas = 0.;
            for (int k=0; k<nloops2; ++k) e14_blas += NormSq(D3[k]-D0[k]);
            e14_blas = sqrt(e14_blas/nloops2);
            e14_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] += T(8,9) * EPART(B4[k].transpose()) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e14_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e14_eigen = sqrt(e14_eigen/nloops2);
            e14_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] += T(8,9) * EPART(B5[k].transpose()) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e14_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e14_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e14_smalleigen = sqrt(e14_smalleigen/nloops2);
            e14_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * A* * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(7,1) * std::conj(A0[k](i,j)) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = T(7,1) * MPART(A1[k]).conjugate() * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_reg = 0.;
            for (int k=0; k<nloops2; ++k) e15_reg += NormSq(D1[k]-D0[k]);
            e15_reg = sqrt(e15_reg/nloops2);
            e15_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = T(7,1) * MPART(A2[k]).conjugate() * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_small = 0.;
            for (int k=0; k<nloops2; ++k) e15_small += NormSq(D2[k]-D0[k]);
            e15_small = sqrt(e15_small/nloops2);
            e15_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,z71c,BP(A3[k].cptr()),M,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASDNAME(scal) (M,dmone,D3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2);
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASNT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_blas = 0.;
            for (int k=0; k<nloops2; ++k) e15_blas += NormSq(D3[k]-D0[k]);
            e15_blas = sqrt(e15_blas/nloops2);
            e15_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = T(7,1) * EPART(A4[k]).conjugate() * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e15_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e15_eigen = sqrt(e15_eigen/nloops2);
            e15_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = T(7,1) * EPART(A5[k]).conjugate() * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e15_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e15_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e15_smalleigen = sqrt(e15_smalleigen/nloops2);
            e15_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * B* * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(7,1) * std::conj(B0[k](j,i)) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = T(7,1) * MPART(B1[k].adjoint()) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_reg = 0.;
            for (int k=0; k<nloops2; ++k) e17_reg += NormSq(D1[k]-D0[k]);
            e17_reg = sqrt(e17_reg/nloops2);
            e17_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = T(7,1) * MPART(B2[k].adjoint()) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_small = 0.;
            for (int k=0; k<nloops2; ++k) e17_small += NormSq(D2[k]-D0[k]);
            e17_small = sqrt(e17_small/nloops2);
            e17_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLASCT,
                            N,M,z71,BP(B3[k].cptr()),N,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLASCT, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_blas = 0.;
            for (int k=0; k<nloops2; ++k) e17_blas += NormSq(D3[k]-D0[k]);
            e17_blas = sqrt(e17_blas/nloops2);
            e17_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = T(7,1) * EPART(B4[k].adjoint()) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e17_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e17_eigen = sqrt(e17_eigen/nloops2);
            e17_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = T(7,1) * EPART(B5[k].adjoint()) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e17_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e17_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e17_smalleigen = sqrt(e17_smalleigen/nloops2);
            e17_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * A* * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(8,9) * std::conj(A0[k](i,j)) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += T(8,9) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] += T(8,9) * MPART(A1[k]).conjugate() * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_reg = 0.;
            for (int k=0; k<nloops2; ++k) e16_reg += NormSq(D1[k]-D0[k]);
            e16_reg = sqrt(e16_reg/nloops2);
            e16_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] += T(8,9) * MPART(A2[k]).conjugate() * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_small = 0.;
            for (int k=0; k<nloops2; ++k) e16_small += NormSq(D2[k]-D0[k]);
            e16_small = sqrt(e16_small/nloops2);
            e16_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASDNAME(scal) (M,dmone,D3[k].realPart().ptr()+1,2);
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,z89c,BP(A3[k].cptr()),M,BP(C3[k].cptr()),1,
                            one,BP(D3[k].ptr()),1 BLAS1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASDNAME(scal) (M,dmone,D3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(temp.ptr()),1);
            BLASDNAME(scal) (N,dmone,temp.realPart().ptr()+1,2);
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASNT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(temp.ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASDNAME(scal) (N,dmone,temp.realPart().ptr()+1,2);
            BLASNAME(axpy) (N,z89,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_blas = 0.;
            for (int k=0; k<nloops2; ++k) e16_blas += NormSq(D3[k]-D0[k]);
            e16_blas = sqrt(e16_blas/nloops2);
            e16_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] += T(8,9) * EPART(A4[k]).conjugate() * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e16_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e16_eigen = sqrt(e16_eigen/nloops2);
            e16_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] += T(8,9) * EPART(A5[k]).conjugate() * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e16_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e16_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e16_smalleigen = sqrt(e16_smalleigen/nloops2);
            e16_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * B* * C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(8,9) * std::conj(B0[k](j,i)) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i) += T(8,9) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] += T(8,9) * MPART(B1[k].adjoint()) * C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_reg = 0.;
            for (int k=0; k<nloops2; ++k) e18_reg += NormSq(D1[k]-D0[k]);
            e18_reg = sqrt(e18_reg/nloops2);
            e18_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] += T(8,9) * MPART(B2[k].adjoint()) * C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_small = 0.;
            for (int k=0; k<nloops2; ++k) e18_small += NormSq(D2[k]-D0[k]);
            e18_small = sqrt(e18_small/nloops2);
            e18_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLASCT,
                            N,M,z89,BP(B3[k].cptr()),N,BP(C3[k].cptr()),1,
                            one,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLASCT, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(temp.ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(axpy) (N,z89,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_blas = 0.;
            for (int k=0; k<nloops2; ++k) e18_blas += NormSq(D3[k]-D0[k]);
            e18_blas = sqrt(e18_blas/nloops2);
            e18_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] += T(8,9) * EPART(B4[k].adjoint()) * C4[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e18_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e18_eigen = sqrt(e18_eigen/nloops2);
            e18_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] += T(8,9) * EPART(B5[k].adjoint()) * C5[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e18_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e18_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e18_smalleigen = sqrt(e18_smalleigen/nloops2);
            e18_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * A * C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(7,1) * A0[k](i,j) * std::conj(C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i) += T(7,1) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = T(7,1) * MPART(A1[k]) * C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_reg = 0.;
            for (int k=0; k<nloops2; ++k) e19_reg += NormSq(D1[k]-D0[k]);
            e19_reg = sqrt(e19_reg/nloops2);
            e19_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = T(7,1) * MPART(A2[k]) * C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_small = 0.;
            for (int k=0; k<nloops2; ++k) e19_small += NormSq(D2[k]-D0[k]);
            e19_small = sqrt(e19_small/nloops2);
            e19_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,z71,BP(A3[k].cptr()),M,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2);
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASNT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_blas = 0.;
            for (int k=0; k<nloops2; ++k) e19_blas += NormSq(D3[k]-D0[k]);
            e19_blas = sqrt(e19_blas/nloops2);
            e19_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = T(7,1) * EPART(A4[k]) * C4[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e19_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e19_eigen = sqrt(e19_eigen/nloops2);
            e19_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = T(7,1) * EPART(A5[k]) * C5[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e19_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e19_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e19_smalleigen = sqrt(e19_smalleigen/nloops2);
            e19_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * B * C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(7,1) * B0[k](j,i) * std::conj(C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i) += T(7,1) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = T(7,1) * MPART(B1[k].transpose()) * C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_reg = 0.;
            for (int k=0; k<nloops2; ++k) e20_reg += NormSq(D1[k]-D0[k]);
            e20_reg = sqrt(e20_reg/nloops2);
            e20_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = T(7,1) * MPART(B2[k].transpose()) * C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_small = 0.;
            for (int k=0; k<nloops2; ++k) e20_small += NormSq(D2[k]-D0[k]);
            e20_small = sqrt(e20_small/nloops2);
            e20_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLASCT,
                            N,M,z71c,BP(B3[k].cptr()),N,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
            BLASDNAME(scal) (M,dmone,D3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLAST, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_blas = 0.;
            for (int k=0; k<nloops2; ++k) e20_blas += NormSq(D3[k]-D0[k]);
            e20_blas = sqrt(e20_blas/nloops2);
            e20_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = T(7,1) * EPART(B4[k].transpose()) * C4[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e20_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e20_eigen = sqrt(e20_eigen/nloops2);
            e20_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = T(7,1) * EPART(B5[k].transpose()) * C5[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e20_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e20_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e20_smalleigen = sqrt(e20_smalleigen/nloops2);
            e20_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * A * C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(8,9) * A0[k](i,j) * std::conj(C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i) += T(8,9) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] += T(8,9) * MPART(A1[k]) * C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_reg = 0.;
            for (int k=0; k<nloops2; ++k) e21_reg += NormSq(D1[k]-D0[k]);
            e21_reg = sqrt(e21_reg/nloops2);
            e21_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] += T(8,9) * MPART(A2[k]) * C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_small = 0.;
            for (int k=0; k<nloops2; ++k) e21_small += NormSq(D2[k]-D0[k]);
            e21_small = sqrt(e21_small/nloops2);
            e21_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,z89,BP(A3[k].cptr()),M,BP(C3[k].cptr()),1,
                            one,BP(D3[k].ptr()),1 BLAS1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(temp.ptr()),1);
            BLASDNAME(scal) (N,dmone,temp.realPart().ptr()+1,2);
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASNT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(temp.ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(axpy) (N,z89,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_blas = 0.;
            for (int k=0; k<nloops2; ++k) e21_blas += NormSq(D3[k]-D0[k]);
            e21_blas = sqrt(e21_blas/nloops2);
            e21_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] += T(8,9) * EPART(A4[k]) * C4[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e21_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e21_eigen = sqrt(e21_eigen/nloops2);
            e21_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] += T(8,9) * EPART(A5[k]) * C5[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e21_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e21_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e21_smalleigen = sqrt(e21_smalleigen/nloops2);
            e21_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * B * C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(8,9) * B0[k](j,i) * std::conj(C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i) += T(8,9) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] += T(8,9) * MPART(B1[k].transpose()) * C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_reg = 0.;
            for (int k=0; k<nloops2; ++k) e22_reg += NormSq(D1[k]-D0[k]);
            e22_reg = sqrt(e22_reg/nloops2);
            e22_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] += T(8,9) * MPART(B2[k].transpose()) * C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_small = 0.;
            for (int k=0; k<nloops2; ++k) e22_small += NormSq(D2[k]-D0[k]);
            e22_small = sqrt(e22_small/nloops2);
            e22_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(gemv) (BLASCM BLAST,
                            N,M,z89,BP(B3[k].cptr()),N,BP(C3[k].cptr()),1,
                            one,BP(D3[k].ptr()),1 BLAS1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(temp.ptr()),1);
            BLASDNAME(scal) (N,dmone,temp.realPart().ptr()+1,2);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLAST, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(temp.ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(axpy) (N,z89,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_blas = 0.;
            for (int k=0; k<nloops2; ++k) e22_blas += NormSq(D3[k]-D0[k]);
            e22_blas = sqrt(e22_blas/nloops2);
            e22_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] += T(8,9) * EPART(B4[k].transpose()) * C4[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e22_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e22_eigen = sqrt(e22_eigen/nloops2);
            e22_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] += T(8,9) * EPART(B5[k].transpose()) * C5[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e22_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e22_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e22_smalleigen = sqrt(e22_smalleigen/nloops2);
            e22_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * A* * C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(7,1) * std::conj(A0[k](i,j) * C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i) += T(7,1) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = T(7,1) * MPART(A1[k]).conjugate() * C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e23_reg = 0.;
            for (int k=0; k<nloops2; ++k) e23_reg += NormSq(D1[k]-D0[k]);
            e23_reg = sqrt(e23_reg/nloops2);
            e23_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = T(7,1) * MPART(A2[k]).conjugate() * C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e23_small = 0.;
            for (int k=0; k<nloops2; ++k) e23_small += NormSq(D2[k]-D0[k]);
            e23_small = sqrt(e23_small/nloops2);
            e23_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,z71c,BP(A3[k].cptr()),M,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
            BLASDNAME(scal) (M,dmone,D3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASNT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e23_blas = 0.;
            for (int k=0; k<nloops2; ++k) e23_blas += NormSq(D3[k]-D0[k]);
            e23_blas = sqrt(e23_blas/nloops2);
            e23_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = T(7,1) * EPART(A4[k]).conjugate() * C4[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e23_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e23_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e23_eigen = sqrt(e23_eigen/nloops2);
            e23_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = T(7,1) * EPART(A5[k]).conjugate() * C5[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e23_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e23_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e23_smalleigen = sqrt(e23_smalleigen/nloops2);
            e23_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * B* * C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(7,1) * std::conj(B0[k](j,i) * C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i) += T(7,1) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] = T(7,1) * MPART(B1[k].adjoint()) * C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e24_reg = 0.;
            for (int k=0; k<nloops2; ++k) e24_reg += NormSq(D1[k]-D0[k]);
            e24_reg = sqrt(e24_reg/nloops2);
            e24_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] = T(7,1) * MPART(B2[k].adjoint()) * C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e24_small = 0.;
            for (int k=0; k<nloops2; ++k) e24_small += NormSq(D2[k]-D0[k]);
            e24_small = sqrt(e24_small/nloops2);
            e24_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(gemv) (BLASCM BLAST,
                            N,M,z71c,BP(B3[k].cptr()),N,BP(C3[k].cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
            BLASDNAME(scal) (M,dmone,D3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLASCT, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e24_blas = 0.;
            for (int k=0; k<nloops2; ++k) e24_blas += NormSq(D3[k]-D0[k]);
            e24_blas = sqrt(e24_blas/nloops2);
            e24_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] = T(7,1) * EPART(B4[k].adjoint()) * C4[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e24_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e24_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e24_eigen = sqrt(e24_eigen/nloops2);
            e24_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] = T(7,1) * EPART(B5[k].adjoint()) * C5[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e24_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e24_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e24_smalleigen = sqrt(e24_smalleigen/nloops2);
            e24_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * A* * C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(8,9) * std::conj(A0[k](i,j) * C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i) += T(8,9) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] += T(8,9) * MPART(A1[k]).conjugate() * C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e25_reg = 0.;
            for (int k=0; k<nloops2; ++k) e25_reg += NormSq(D1[k]-D0[k]);
            e25_reg = sqrt(e25_reg/nloops2);
            e25_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] += T(8,9) * MPART(A2[k]).conjugate() * C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e25_small = 0.;
            for (int k=0; k<nloops2; ++k) e25_small += NormSq(D2[k]-D0[k]);
            e25_small = sqrt(e25_small/nloops2);
            e25_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASDNAME(scal) (M,dmone,D3[k].realPart().ptr()+1,2);
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,z89c,BP(A3[k].cptr()),M,BP(C3[k].cptr()),1,
                            one,BP(D3[k].ptr()),1 BLAS1);
            BLASDNAME(scal) (M,dmone,D3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASNT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(temp.ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASDNAME(scal) (N,dmone,temp.realPart().ptr()+1,2);
            BLASNAME(axpy) (N,z89,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e25_blas = 0.;
            for (int k=0; k<nloops2; ++k) e25_blas += NormSq(D3[k]-D0[k]);
            e25_blas = sqrt(e25_blas/nloops2);
            e25_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] += T(8,9) * EPART(A4[k]).conjugate() * C4[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e25_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e25_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e25_eigen = sqrt(e25_eigen/nloops2);
            e25_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] += T(8,9) * EPART(A5[k]).conjugate() * C5[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e25_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e25_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e25_smalleigen = sqrt(e25_smalleigen/nloops2);
            e25_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * B* * C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i) += T(8,9) * std::conj(B0[k](j,i) * C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i) += T(8,9) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] += T(8,9) * MPART(B1[k].adjoint()) * C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e26_reg = 0.;
            for (int k=0; k<nloops2; ++k) e26_reg += NormSq(D1[k]-D0[k]);
            e26_reg = sqrt(e26_reg/nloops2);
            e26_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] += T(8,9) * MPART(B2[k].adjoint()) * C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e26_small = 0.;
            for (int k=0; k<nloops2; ++k) e26_small += NormSq(D2[k]-D0[k]);
            e26_small = sqrt(e26_small/nloops2);
            e26_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(gemv) (BLASCM BLASCT,
                            N,M,z89,BP(B3[k].cptr()),N,BP(C3[k].cptr()),1,
                            one,BP(D3[k].ptr()),1 BLAS1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(temp.ptr()),1);
            BLASDNAME(scal) (N,dmone,temp.realPart().ptr()+1,2);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLASCT, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(temp.ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(axpy) (N,z89,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e26_blas = 0.;
            for (int k=0; k<nloops2; ++k) e26_blas += NormSq(D3[k]-D0[k]);
            e26_blas = sqrt(e26_blas/nloops2);
            e26_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k] += T(8,9) * EPART(B4[k].adjoint()) * C4[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e26_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e26_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e26_eigen = sqrt(e26_eigen/nloops2);
            e26_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k] += T(8,9) * EPART(B5[k].adjoint()) * C5[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e26_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e26_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e26_smalleigen = sqrt(e26_smalleigen/nloops2);
            e26_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif
#endif

#ifdef AISSQUARE
#if 1 // D *= A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(j,i))
                            D0[k](i) += A0[k](j,i) * C0[k](j);
                        else if (UNITPART(j,i))
                            D0[k](i) += C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) D1[k] = C1[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] *= MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t27_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e27_reg = 0.;
            for (int k=0; k<nloops2; ++k) e27_reg += NormSq(D1[k]-D0[k]);
            e27_reg = sqrt(e27_reg/nloops2);
            e27_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] *= MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t27_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e27_small = 0.;
            for (int k=0; k<nloops2; ++k) e27_small += NormSq(D2[k]-D0[k]);
            e27_small = sqrt(e27_small/nloops2);
            e27_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        for (int k=0; k<nloops2; ++k) D3[k] = C3[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(D3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(gemv) (BLASCM BLAST,
                            N,M,one,BP(A3[k].cptr()),M,BP(temp.cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(trmv) (BLASCM BLASUPLO, BLAST, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t27_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e27_blas = 0.;
            for (int k=0; k<nloops2; ++k) e27_blas += NormSq(D3[k]-D0[k]);
            e27_blas = sqrt(e27_blas/nloops2);
            e27_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k].transpose() *= EPART(A4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t27_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e27_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e27_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e27_eigen = sqrt(e27_eigen/nloops2);
            e27_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k].transpose() *= EPART(A5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t27_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e27_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e27_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e27_smalleigen = sqrt(e27_smalleigen/nloops2);
            e27_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D *= B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(j,i))
                            D0[k](i) += B0[k](i,j) * C0[k](j);
                        else if (UNITPART(j,i))
                            D0[k](i) += C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) D1[k] = C1[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] *= MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t28_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e28_reg = 0.;
            for (int k=0; k<nloops2; ++k) e28_reg += NormSq(D1[k]-D0[k]);
            e28_reg = sqrt(e28_reg/nloops2);
            e28_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] *= MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t28_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e28_small = 0.;
            for (int k=0; k<nloops2; ++k) e28_small += NormSq(D2[k]-D0[k]);
            e28_small = sqrt(e28_small/nloops2);
            e28_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        for (int k=0; k<nloops2; ++k) D3[k] = C3[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(D3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,one,BP(B3[k].cptr()),N,BP(temp.cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLASNT, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t28_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e28_blas = 0.;
            for (int k=0; k<nloops2; ++k) e28_blas += NormSq(D3[k]-D0[k]);
            e28_blas = sqrt(e28_blas/nloops2);
            e28_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k].transpose() *= EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t28_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e28_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e28_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e28_eigen = sqrt(e28_eigen/nloops2);
            e28_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k].transpose() *= EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t28_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e28_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e28_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e28_smalleigen = sqrt(e28_smalleigen/nloops2);
            e28_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D *= 7 * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(j,i))
                            D0[k](i) += RT(7) * A0[k](j,i) * C0[k](j);
                        else if (UNITPART(j,i))
                            D0[k](i) += RT(7) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) D1[k] = C1[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] *= RT(7) * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t29_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e29_reg = 0.;
            for (int k=0; k<nloops2; ++k) e29_reg += NormSq(D1[k]-D0[k]);
            e29_reg = sqrt(e29_reg/nloops2);
            e29_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] *= RT(7) * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t29_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e29_small = 0.;
            for (int k=0; k<nloops2; ++k) e29_small += NormSq(D2[k]-D0[k]);
            e29_small = sqrt(e29_small/nloops2);
            e29_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        for (int k=0; k<nloops2; ++k) D3[k] = C3[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(D3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(gemv) (BLASCM BLAST,
                            N,M,seven,BP(A3[k].cptr()),M,BP(temp.cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(trmv) (BLASCM BLASUPLO, BLAST, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,seven,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t29_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e29_blas = 0.;
            for (int k=0; k<nloops2; ++k) e29_blas += NormSq(D3[k]-D0[k]);
            e29_blas = sqrt(e29_blas/nloops2);
            e29_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k].transpose() *= RT(7) * EPART(A4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t29_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e29_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e29_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e29_eigen = sqrt(e29_eigen/nloops2);
            e29_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k].transpose() *= RT(7) * EPART(A5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t29_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e29_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e29_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e29_smalleigen = sqrt(e29_smalleigen/nloops2);
            e29_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D *= 7 * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(j,i))
                            D0[k](i) += RT(7) * B0[k](i,j) * C0[k](j);
                        else if (UNITPART(j,i))
                            D0[k](i) += RT(7) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) D1[k] = C1[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] *= RT(7) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t30_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e30_reg = 0.;
            for (int k=0; k<nloops2; ++k) e30_reg += NormSq(D1[k]-D0[k]);
            e30_reg = sqrt(e30_reg/nloops2);
            e30_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] *= RT(7) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t30_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e30_small = 0.;
            for (int k=0; k<nloops2; ++k) e30_small += NormSq(D2[k]-D0[k]);
            e30_small = sqrt(e30_small/nloops2);
            e30_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        for (int k=0; k<nloops2; ++k) D3[k] = C3[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(D3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,seven,BP(B3[k].cptr()),N,BP(temp.cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLASNT, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,seven,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t30_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e30_blas = 0.;
            for (int k=0; k<nloops2; ++k) e30_blas += NormSq(D3[k]-D0[k]);
            e30_blas = sqrt(e30_blas/nloops2);
            e30_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k].transpose() *= RT(7) * EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t30_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e30_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e30_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e30_eigen = sqrt(e30_eigen/nloops2);
            e30_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k].transpose() *= RT(7) * EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t30_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e30_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e30_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e30_smalleigen = sqrt(e30_smalleigen/nloops2);
            e30_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#ifdef TISCOMPLEX
#if 1 // D *= (7,1) * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(j,i))
                            D0[k](i) += T(7,1) * A0[k](j,i) * C0[k](j);
                        else if (UNITPART(j,i))
                            D0[k](i) += T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) D1[k] = C1[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] *= T(7,1) * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t31_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e31_reg = 0.;
            for (int k=0; k<nloops2; ++k) e31_reg += NormSq(D1[k]-D0[k]);
            e31_reg = sqrt(e31_reg/nloops2);
            e31_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] *= T(7,1) * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t31_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e31_small = 0.;
            for (int k=0; k<nloops2; ++k) e31_small += NormSq(D2[k]-D0[k]);
            e31_small = sqrt(e31_small/nloops2);
            e31_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        for (int k=0; k<nloops2; ++k) D3[k] = C3[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(D3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(gemv) (BLASCM BLAST,
                            N,M,z71,BP(A3[k].cptr()),M,BP(temp.cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(trmv) (BLASCM BLASUPLO, BLAST, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t31_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e31_blas = 0.;
            for (int k=0; k<nloops2; ++k) e31_blas += NormSq(D3[k]-D0[k]);
            e31_blas = sqrt(e31_blas/nloops2);
            e31_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k].transpose() *= T(7,1) * EPART(A4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t31_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e31_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e31_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e31_eigen = sqrt(e31_eigen/nloops2);
            e31_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k].transpose() *= T(7,1) * EPART(A5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t31_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e31_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e31_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e31_smalleigen = sqrt(e31_smalleigen/nloops2);
            e31_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D *= (7,1) * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(j,i))
                            D0[k](i) += T(7,1) * B0[k](i,j) * C0[k](j);
                        else if (UNITPART(j,i))
                            D0[k](i) += T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) D1[k] = C1[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] *= T(7,1) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t32_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e32_reg = 0.;
            for (int k=0; k<nloops2; ++k) e32_reg += NormSq(D1[k]-D0[k]);
            e32_reg = sqrt(e32_reg/nloops2);
            e32_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] *= T(7,1) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t32_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e32_small = 0.;
            for (int k=0; k<nloops2; ++k) e32_small += NormSq(D2[k]-D0[k]);
            e32_small = sqrt(e32_small/nloops2);
            e32_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        for (int k=0; k<nloops2; ++k) D3[k] = C3[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(D3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,z71,BP(B3[k].cptr()),N,BP(temp.cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLASNT, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t32_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e32_blas = 0.;
            for (int k=0; k<nloops2; ++k) e32_blas += NormSq(D3[k]-D0[k]);
            e32_blas = sqrt(e32_blas/nloops2);
            e32_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k].transpose() *= T(7,1) * EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t32_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e32_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e32_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e32_eigen = sqrt(e32_eigen/nloops2);
            e32_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k].transpose() *= T(7,1) * EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t32_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e32_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e32_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e32_smalleigen = sqrt(e32_smalleigen/nloops2);
            e32_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D *= (7,1) * A*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(j,i))
                            D0[k](i) += T(7,1) * std::conj(A0[k](j,i)) * C0[k](j);
                        else if (UNITPART(j,i))
                            D0[k](i) += T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) D1[k] = C1[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] *= T(7,1) * MPART(A1[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t33_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e33_reg = 0.;
            for (int k=0; k<nloops2; ++k) e33_reg += NormSq(D1[k]-D0[k]);
            e33_reg = sqrt(e33_reg/nloops2);
            e33_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] *= T(7,1) * MPART(A2[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t33_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e33_small = 0.;
            for (int k=0; k<nloops2; ++k) e33_small += NormSq(D2[k]-D0[k]);
            e33_small = sqrt(e33_small/nloops2);
            e33_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        for (int k=0; k<nloops2; ++k) D3[k] = C3[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(D3[k].cptr()),1,BP(temp.ptr()),1);
            BLASNAME(gemv) (BLASCM BLASCT,
                            N,M,z71,BP(A3[k].cptr()),M,BP(temp.cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
#elif ((PART >= 2) && (PART <= 5))
            BLASNAME(trmv) (BLASCM BLASUPLO, BLASCT, BLASDIAG,
                            N,BP(A3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t33_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e33_blas = 0.;
            for (int k=0; k<nloops2; ++k) e33_blas += NormSq(D3[k]-D0[k]);
            e33_blas = sqrt(e33_blas/nloops2);
            e33_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k].transpose() *= T(7,1) * EPART(A4[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t33_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e33_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e33_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e33_eigen = sqrt(e33_eigen/nloops2);
            e33_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k].transpose() *= T(7,1) * EPART(A5[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t33_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e33_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e33_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e33_smalleigen = sqrt(e33_smalleigen/nloops2);
            e33_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D *= (7,1) * B*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    D0[k](i) = T(0);
                    for(int j=0;j<N;++j) {
                        if (INPART(j,i))
                            D0[k](i) += T(7,1) * std::conj(B0[k](i,j)) * C0[k](j);
                        else if (UNITPART(j,i))
                            D0[k](i) += T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        for (int k=0; k<nloops2; ++k) D1[k] = C1[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D1[k] *= T(7,1) * MPART(B1[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t34_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e34_reg = 0.;
            for (int k=0; k<nloops2; ++k) e34_reg += NormSq(D1[k]-D0[k]);
            e34_reg = sqrt(e34_reg/nloops2);
            e34_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        for (int k=0; k<nloops2; ++k) D2[k] = C2[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D2[k] *= T(7,1) * MPART(B2[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t34_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e34_small = 0.;
            for (int k=0; k<nloops2; ++k) e34_small += NormSq(D2[k]-D0[k]);
            e34_small = sqrt(e34_small/nloops2);
            e34_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        for (int k=0; k<nloops2; ++k) D3[k] = C3[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            tmv::Vector<T> temp(N);
            BLASNAME(copy) (N,BP(D3[k].cptr()),1,BP(temp.ptr()),1);
            BLASDNAME(scal) (N,dmone,temp.realPart().ptr()+1,2);
            BLASNAME(gemv) (BLASCM BLASNT,
                            M,N,z71c,BP(B3[k].cptr()),N,BP(temp.cptr()),1,
                            zero,BP(D3[k].ptr()),1 BLAS1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2);
            BLASNAME(trmv) (BLASCM BLASUPLO2, BLASNT, BLASDIAG,
                            N,BP(B3[k].cptr()),N,BP(D3[k].ptr()),1 BLAS1 BLAS1 BLAS1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t34_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e34_blas = 0.;
            for (int k=0; k<nloops2; ++k) e34_blas += NormSq(D3[k]-D0[k]);
            e34_blas = sqrt(e34_blas/nloops2);
            e34_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        for (int k=0; k<nloops2; ++k) D4[k] = C4[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            D4[k].transpose() *= T(7,1) * EPART(B4[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t34_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e34_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e34_eigen += tmv::TMV_NORM(D4[k](i)-D0[k](i));
            e34_eigen = sqrt(e34_eigen/nloops2);
            e34_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        for (int k=0; k<nloops2x; ++k) D5[k] = C5[k];

        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            D5[k].transpose() *= T(7,1) * EPART(B5[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t34_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e34_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) e34_smalleigen += tmv::TMV_NORM(D5[k](i)-D0[k](i));
            e34_smalleigen = sqrt(e34_smalleigen/nloops2);
            e34_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif
#endif
#endif

    }

    std::cout<<"D = A * C               "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"D = B * C               "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"D = -A * C              "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"D = -B * C              "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"D = 7 * A * C           "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
    std::cout<<"D = 7 * B * C           "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"D -= A * C              "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"D -= B * C              "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
    std::cout<<"D += 8 * A * C          "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
    std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
    std::cout<<"D += 8 * B * C          "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
    std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"D = (7,1) * A * C       "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
    std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * B * C       "<<t12_reg<<"  "<<t12_small<<"  "<<t12_blas;
    std::cout<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * A * C      "<<t13_reg<<"  "<<t13_small<<"  "<<t13_blas;
    std::cout<<"  "<<t13_eigen<<"  "<<t13_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * B * C      "<<t14_reg<<"  "<<t14_small<<"  "<<t14_blas;
    std::cout<<"  "<<t14_eigen<<"  "<<t14_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * A* * C      "<<t15_reg<<"  "<<t15_small<<"  "<<t15_blas;
    std::cout<<"  "<<t15_eigen<<"  "<<t15_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * B* * C      "<<t17_reg<<"  "<<t17_small<<"  "<<t17_blas;
    std::cout<<"  "<<t17_eigen<<"  "<<t17_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * A* * C     "<<t16_reg<<"  "<<t16_small<<"  "<<t16_blas;
    std::cout<<"  "<<t16_eigen<<"  "<<t16_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * B* * C     "<<t18_reg<<"  "<<t18_small<<"  "<<t18_blas;
    std::cout<<"  "<<t18_eigen<<"  "<<t18_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * A * C*      "<<t19_reg<<"  "<<t19_small<<"  "<<t19_blas;
    std::cout<<"  "<<t19_eigen<<"  "<<t19_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * B * C*      "<<t20_reg<<"  "<<t20_small<<"  "<<t20_blas;
    std::cout<<"  "<<t20_eigen<<"  "<<t20_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * A * C*     "<<t21_reg<<"  "<<t21_small<<"  "<<t21_blas;
    std::cout<<"  "<<t21_eigen<<"  "<<t21_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * B * C*     "<<t22_reg<<"  "<<t22_small<<"  "<<t22_blas;
    std::cout<<"  "<<t22_eigen<<"  "<<t22_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * A* * C*     "<<t23_reg<<"  "<<t23_small<<"  "<<t23_blas;
    std::cout<<"  "<<t23_eigen<<"  "<<t23_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * B* * C*     "<<t24_reg<<"  "<<t24_small<<"  "<<t24_blas;
    std::cout<<"  "<<t24_eigen<<"  "<<t24_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * A* * C*    "<<t25_reg<<"  "<<t25_small<<"  "<<t25_blas;
    std::cout<<"  "<<t25_eigen<<"  "<<t25_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * B* * C*    "<<t26_reg<<"  "<<t26_small<<"  "<<t26_blas;
    std::cout<<"  "<<t26_eigen<<"  "<<t26_smalleigen<<std::endl;
#endif
#ifdef AISSQUARE
    std::cout<<"D *= A                  "<<t27_reg<<"  "<<t27_small<<"  "<<t27_blas;
    std::cout<<"  "<<t27_eigen<<"  "<<t27_smalleigen<<std::endl;
    std::cout<<"D *= B                  "<<t28_reg<<"  "<<t28_small<<"  "<<t28_blas;
    std::cout<<"  "<<t28_eigen<<"  "<<t28_smalleigen<<std::endl;
    std::cout<<"D *= 7 * A              "<<t29_reg<<"  "<<t29_small<<"  "<<t29_blas;
    std::cout<<"  "<<t29_eigen<<"  "<<t29_smalleigen<<std::endl;
    std::cout<<"D *= 7 * B              "<<t30_reg<<"  "<<t30_small<<"  "<<t30_blas;
    std::cout<<"  "<<t30_eigen<<"  "<<t30_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"D *= (7,1) * A          "<<t31_reg<<"  "<<t31_small<<"  "<<t31_blas;
    std::cout<<"  "<<t31_eigen<<"  "<<t31_smalleigen<<std::endl;
    std::cout<<"D *= (7,1) * B          "<<t32_reg<<"  "<<t32_small<<"  "<<t32_blas;
    std::cout<<"  "<<t32_eigen<<"  "<<t32_smalleigen<<std::endl;
    std::cout<<"D *= (7,1) * A*         "<<t33_reg<<"  "<<t33_small<<"  "<<t33_blas;
    std::cout<<"  "<<t33_eigen<<"  "<<t33_smalleigen<<std::endl;
    std::cout<<"D *= (7,1) * B*         "<<t34_reg<<"  "<<t34_small<<"  "<<t34_blas;
    std::cout<<"  "<<t34_eigen<<"  "<<t34_smalleigen<<std::endl;
#endif
#endif

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"D = A * C               "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"D = B * C               "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"D = -A * C              "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"D = -B * C              "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"D = 7 * A * C           "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
    std::cout<<"D = 7 * B * C           "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
    std::cout<<"D -= A * C              "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
    std::cout<<"D -= B * C              "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
    std::cout<<"D += 8 * A * C          "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
    std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
    std::cout<<"D += 8 * B * C          "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
    std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"D = (7,1) * A * C       "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
    std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * B * C       "<<e12_reg<<"  "<<e12_small<<"  "<<e12_blas;
    std::cout<<"  "<<e12_eigen<<"  "<<e12_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * A * C      "<<e13_reg<<"  "<<e13_small<<"  "<<e13_blas;
    std::cout<<"  "<<e13_eigen<<"  "<<e13_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * B * C      "<<e14_reg<<"  "<<e14_small<<"  "<<e14_blas;
    std::cout<<"  "<<e14_eigen<<"  "<<e14_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * A* * C      "<<e15_reg<<"  "<<e15_small<<"  "<<e15_blas;
    std::cout<<"  "<<e15_eigen<<"  "<<e15_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * B* * C      "<<e17_reg<<"  "<<e17_small<<"  "<<e17_blas;
    std::cout<<"  "<<e17_eigen<<"  "<<e17_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * A* * C     "<<e16_reg<<"  "<<e16_small<<"  "<<e16_blas;
    std::cout<<"  "<<e16_eigen<<"  "<<e16_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * B* * C     "<<e18_reg<<"  "<<e18_small<<"  "<<e18_blas;
    std::cout<<"  "<<e18_eigen<<"  "<<e18_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * A * C*      "<<e19_reg<<"  "<<e19_small<<"  "<<e19_blas;
    std::cout<<"  "<<e19_eigen<<"  "<<e19_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * B * C*      "<<e20_reg<<"  "<<e20_small<<"  "<<e20_blas;
    std::cout<<"  "<<e20_eigen<<"  "<<e20_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * A * C*     "<<e21_reg<<"  "<<e21_small<<"  "<<e21_blas;
    std::cout<<"  "<<e21_eigen<<"  "<<e21_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * B * C*     "<<e22_reg<<"  "<<e22_small<<"  "<<e22_blas;
    std::cout<<"  "<<e22_eigen<<"  "<<e22_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * A* * C*     "<<e23_reg<<"  "<<e23_small<<"  "<<e23_blas;
    std::cout<<"  "<<e23_eigen<<"  "<<e23_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * B* * C*     "<<e24_reg<<"  "<<e24_small<<"  "<<e24_blas;
    std::cout<<"  "<<e24_eigen<<"  "<<e24_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * A* * C*    "<<e25_reg<<"  "<<e25_small<<"  "<<e25_blas;
    std::cout<<"  "<<e25_eigen<<"  "<<e25_smalleigen<<std::endl;
    std::cout<<"D += (8,9) * B* * C*    "<<e26_reg<<"  "<<e26_small<<"  "<<e26_blas;
    std::cout<<"  "<<e26_eigen<<"  "<<e26_smalleigen<<std::endl;
#endif
#ifdef AISSQUARE
    std::cout<<"D *= A                  "<<e27_reg<<"  "<<e27_small<<"  "<<e27_blas;
    std::cout<<"  "<<e27_eigen<<"  "<<e27_smalleigen<<std::endl;
    std::cout<<"D *= B                  "<<e28_reg<<"  "<<e28_small<<"  "<<e28_blas;
    std::cout<<"  "<<e28_eigen<<"  "<<e28_smalleigen<<std::endl;
    std::cout<<"D *= 7 * A              "<<e29_reg<<"  "<<e29_small<<"  "<<e29_blas;
    std::cout<<"  "<<e29_eigen<<"  "<<e29_smalleigen<<std::endl;
    std::cout<<"D *= 7 * B              "<<e30_reg<<"  "<<e30_small<<"  "<<e30_blas;
    std::cout<<"  "<<e30_eigen<<"  "<<e30_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"D *= (7,1) * A          "<<e31_reg<<"  "<<e31_small<<"  "<<e31_blas;
    std::cout<<"  "<<e31_eigen<<"  "<<e31_smalleigen<<std::endl;
    std::cout<<"D *= (7,1) * B          "<<e32_reg<<"  "<<e32_small<<"  "<<e32_blas;
    std::cout<<"  "<<e32_eigen<<"  "<<e32_smalleigen<<std::endl;
    std::cout<<"D *= (7,1) * A*         "<<e33_reg<<"  "<<e33_small<<"  "<<e33_blas;
    std::cout<<"  "<<e33_eigen<<"  "<<e33_smalleigen<<std::endl;
    std::cout<<"D *= (7,1) * B*         "<<e34_reg<<"  "<<e34_small<<"  "<<e34_blas;
    std::cout<<"  "<<e34_eigen<<"  "<<e34_smalleigen<<std::endl;
#endif
#endif
    std::cout<<"\n\n";
#endif
}
#endif

#if (PART != 1 && defined(DORANK1))
#undef DORANK1
#endif

#ifdef DORANK1
static void Rank1Update(
    std::vector<tmv::Matrix<T> >& A1,
    std::vector<tmv::Matrix<T> >& B1,
    const std::vector<tmv::Vector<T> >& C1,
    const std::vector<tmv::Vector<T> >& D1)
{
#ifdef ERRORCHECK
    std::vector<tmv::Matrix<T> > A0 = A1;
    std::vector<tmv::Matrix<T> > B0 = B1;
    std::vector<tmv::Vector<T> > C0 = C1;
    std::vector<tmv::Vector<T> > D0 = D1;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallMatrix<T,M,N> > A2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,M> > B2(nloops2);
    std::vector<tmv::SmallVector<T,N> > C2(nloops2);
    std::vector<tmv::SmallVector<T,M> > D2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        C2[k] = C1[k]; D2[k] = D1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Matrix<T> > A3 = A1;
    std::vector<tmv::Matrix<T> > B3 = B1;
    std::vector<tmv::Vector<T> > C3 = C1;
    std::vector<tmv::Vector<T> > D3 = D1;
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENM> A4(nloops2,Eigen::EIGENM(M,N));
    std::vector<Eigen::EIGENM> B4(nloops2,Eigen::EIGENM(N,M));
    std::vector<Eigen::EIGENV> C4(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENV> D4(nloops2,Eigen::EIGENV(M));
    for(int k=0;k<nloops2;++k) {
        for(int i=0;i<N;++i) C4[k](i) = C1[k](i);
        for(int i=0;i<M;++i) D4[k](i) = D1[k](i);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,M,N> > A5;
    std::vector<Eigen::Matrix<T,N,M> > B5;
    std::vector<Eigen::Matrix<T,N,1> > C5;
    std::vector<Eigen::Matrix<T,M,1> > D5;
    if (nloops2x) { 
        A5.resize(nloops2x); B5.resize(nloops2x);
        C5.resize(nloops2x); D5.resize(nloops2x);
    }
    for(int k=0;k<nloops2x;++k) {
        for(int i=0;i<N;++i) C5[k](i) = C1[k](i);
        for(int i=0;i<M;++i) D5[k](i) = D1[k](i);
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
    double t17_reg=0., t17_small=0., t17_blas=0., t17_eigen=0., t17_smalleigen=0.;
    double t18_reg=0., t18_small=0., t18_blas=0., t18_eigen=0., t18_smalleigen=0.;
    double t19_reg=0., t19_small=0., t19_blas=0., t19_eigen=0., t19_smalleigen=0.;
    double t20_reg=0., t20_small=0., t20_blas=0., t20_eigen=0., t20_smalleigen=0.;
    double t21_reg=0., t21_small=0., t21_blas=0., t21_eigen=0., t21_smalleigen=0.;
    double t22_reg=0., t22_small=0., t22_blas=0., t22_eigen=0., t22_smalleigen=0.;
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
    double e17_reg=0., e17_small=0., e17_blas=0., e17_eigen=0., e17_smalleigen=0.;
    double e18_reg=0., e18_small=0., e18_blas=0., e18_eigen=0., e18_smalleigen=0.;
    double e19_reg=0., e19_small=0., e19_blas=0., e19_eigen=0., e19_smalleigen=0.;
    double e20_reg=0., e20_small=0., e20_blas=0., e20_eigen=0., e20_smalleigen=0.;
    double e21_reg=0., e21_small=0., e21_blas=0., e21_eigen=0., e21_smalleigen=0.;
    double e22_reg=0., e22_small=0., e22_blas=0., e22_eigen=0., e22_smalleigen=0.;
#endif

    for (int n=0; n<nloops1; ++n) {

#if 1 // A = D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            A0[k](i,j) = D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]) = D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(A1[k]-A0[k]);
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A2[k]) = D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(A2[k]-A0[k]);
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            A3[k].setZero();
#ifdef TISCOMPLEX
            BLASNAME(geru) (BLASCM 
                            M,N,one,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                            BP(A3[k].ptr()),M);
#else
            BLASNAME(ger) (BLASCM 
                           M,N,one,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                           BP(A3[k].ptr()),M);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(A3[k]-A0[k]);
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(A4[k]) = D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e1_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(A5[k]) = D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e1_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // B = D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            B0[k](j,i) = D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B1[k].transpose()) = D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) e2_reg += NormSq(B1[k]-B0[k]);
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B2[k].transpose()) = D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) e2_small += NormSq(B2[k]-B0[k]);
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            B3[k].setZero();
#ifdef TISCOMPLEX
            BLASNAME(geru) (BLASCM 
                            N,M,one,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                            BP(B3[k].ptr()),N);
#else
            BLASNAME(ger) (BLASCM 
                           N,M,one,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                           BP(B3[k].ptr()),N);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) e2_blas += NormSq(B3[k]-B0[k]);
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(B4[k].transpose()) = D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e2_eigen += tmv::TMV_NORM(B4[k](j,i)-B0[k](j,i));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(B5[k].transpose()) = D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e2_smalleigen += tmv::TMV_NORM(B5[k](j,i)-B0[k](j,i));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // A = -D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            A0[k](i,j) = -D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]) = -D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) e3_reg += NormSq(A1[k]-A0[k]);
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A2[k]) = -D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) e3_small += NormSq(A2[k]-A0[k]);
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            A3[k].setZero();
#ifdef TISCOMPLEX
            BLASNAME(geru) (BLASCM 
                            M,N,mone,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                            BP(A3[k].ptr()),M);
#else
            BLASNAME(ger) (BLASCM 
                           M,N,mone,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                           BP(A3[k].ptr()),M);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) e3_blas += NormSq(A3[k]-A0[k]);
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(A4[k]) = -D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e3_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(A5[k]) = -D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e3_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // B = -D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            B0[k](j,i) = -D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B1[k].transpose()) = -D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) e4_reg += NormSq(B1[k]-B0[k]);
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B2[k].transpose()) = -D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) e4_small += NormSq(B2[k]-B0[k]);
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            B3[k].setZero();
#ifdef TISCOMPLEX
            BLASNAME(geru) (BLASCM 
                            N,M,mone,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                            BP(B3[k].ptr()),N);
#else
            BLASNAME(ger) (BLASCM 
                           N,M,mone,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                           BP(B3[k].ptr()),N);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) e4_blas += NormSq(B3[k]-B0[k]);
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(B4[k].transpose()) = -D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e4_eigen += tmv::TMV_NORM(B4[k](j,i)-B0[k](j,i));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(B5[k].transpose()) = -D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e4_smalleigen += tmv::TMV_NORM(B5[k](j,i)-B0[k](j,i));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // A = 7 * D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            A0[k](i,j) = RT(7) * D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]) = 7 * D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) e5_reg += NormSq(A1[k]-A0[k]);
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A2[k]) = 7 * D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) e5_small += NormSq(A2[k]-A0[k]);
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            A3[k].setZero();
#ifdef TISCOMPLEX
            BLASNAME(geru) (BLASCM 
                            M,N,seven,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                            BP(A3[k].ptr()),M);
#else
            BLASNAME(ger) (BLASCM 
                           M,N,seven,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                           BP(A3[k].ptr()),M);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) e5_blas += NormSq(A3[k]-A0[k]);
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(A4[k]) = 7 * D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e5_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(A5[k]) = 7 * D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e5_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // B = 7 * D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            B0[k](j,i) = RT(7) * D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B1[k].transpose()) = 7 * D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) e6_reg += NormSq(B1[k]-B0[k]);
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B2[k].transpose()) = 7 * D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) e6_small += NormSq(B2[k]-B0[k]);
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            B3[k].setZero();
#ifdef TISCOMPLEX
            BLASNAME(geru) (BLASCM 
                            N,M,seven,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                            BP(B3[k].ptr()),N);
#else
            BLASNAME(ger) (BLASCM 
                           N,M,seven,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                           BP(B3[k].ptr()),N);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) e6_blas += NormSq(B3[k]-B0[k]);
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(B4[k].transpose()) = 7 * D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e6_eigen += tmv::TMV_NORM(B4[k](j,i)-B0[k](j,i));
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(B5[k].transpose()) = 7 * D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e6_smalleigen += tmv::TMV_NORM(B5[k](j,i)-B0[k](j,i));
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // A -= D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            A0[k](i,j) -= D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]) -= D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) e7_reg += NormSq(A1[k]-A0[k]);
            e7_reg = sqrt(e7_reg/nloops2);
            e7_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A2[k]) -= D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) e7_small += NormSq(A2[k]-A0[k]);
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
            BLASNAME(geru) (BLASCM 
                            M,N,mone,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                            BP(A3[k].ptr()),M);
#else
            BLASNAME(ger) (BLASCM 
                           M,N,mone,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                           BP(A3[k].ptr()),M);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) e7_blas += NormSq(A3[k]-A0[k]);
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(A4[k]) -= D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e7_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(A5[k]) -= D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e7_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // B -= D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            B0[k](j,i) -= D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B1[k].transpose()) -= D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) e8_reg += NormSq(B1[k]-B0[k]);
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B2[k].transpose()) -= D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) e8_small += NormSq(B2[k]-B0[k]);
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
            BLASNAME(geru) (BLASCM 
                            N,M,mone,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                            BP(B3[k].ptr()),N);
#else
            BLASNAME(ger) (BLASCM 
                           N,M,mone,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                           BP(B3[k].ptr()),N);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) e8_blas += NormSq(B3[k]-B0[k]);
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(B4[k].transpose()) -= D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e8_eigen += tmv::TMV_NORM(B4[k](j,i)-B0[k](j,i));
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(B5[k].transpose()) -= D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e8_smalleigen += tmv::TMV_NORM(B5[k](j,i)-B0[k](j,i));
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // A += 8 * D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            A0[k](i,j) += RT(8) * D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]) += 8 * D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_reg = 0.;
            for (int k=0; k<nloops2; ++k) e9_reg += NormSq(A1[k]-A0[k]);
            e9_reg = sqrt(e9_reg/nloops2);
            e9_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A2[k]) += 8 * D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_small = 0.;
            for (int k=0; k<nloops2; ++k) e9_small += NormSq(A2[k]-A0[k]);
            e9_small = sqrt(e9_small/nloops2);
            e9_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
            BLASNAME(geru) (BLASCM 
                            M,N,eight,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                            BP(A3[k].ptr()),M);
#else
            BLASNAME(ger) (BLASCM 
                           M,N,eight,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                           BP(A3[k].ptr()),M);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_blas = 0.;
            for (int k=0; k<nloops2; ++k) e9_blas += NormSq(A3[k]-A0[k]);
            e9_blas = sqrt(e9_blas/nloops2);
            e9_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(A4[k]) += 8 * D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e9_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e9_eigen = sqrt(e9_eigen/nloops2);
            e9_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(A5[k]) += 8 * D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e9_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e9_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e9_smalleigen = sqrt(e9_smalleigen/nloops2);
            e9_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // B += 8 * D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            B0[k](j,i) += RT(8) * D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B1[k].transpose()) += 8 * D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_reg = 0.;
            for (int k=0; k<nloops2; ++k) e10_reg += NormSq(B1[k]-B0[k]);
            e10_reg = sqrt(e10_reg/nloops2);
            e10_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B2[k].transpose()) += 8 * D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_small = 0.;
            for (int k=0; k<nloops2; ++k) e10_small += NormSq(B2[k]-B0[k]);
            e10_small = sqrt(e10_small/nloops2);
            e10_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#ifdef TISCOMPLEX
            BLASNAME(geru) (BLASCM 
                            N,M,eight,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                            BP(B3[k].ptr()),N);
#else
            BLASNAME(ger) (BLASCM 
                           N,M,eight,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                           BP(B3[k].ptr()),N);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_blas = 0.;
            for (int k=0; k<nloops2; ++k) e10_blas += NormSq(B3[k]-B0[k]);
            e10_blas = sqrt(e10_blas/nloops2);
            e10_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(B4[k].transpose()) += 8 * D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e10_eigen += tmv::TMV_NORM(B4[k](j,i)-B0[k](j,i));
            e10_eigen = sqrt(e10_eigen/nloops2);
            e10_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(B5[k].transpose()) += 8 * D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e10_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e10_smalleigen += tmv::TMV_NORM(B5[k](j,i)-B0[k](j,i));
            e10_smalleigen = sqrt(e10_smalleigen/nloops2);
            e10_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#ifdef TISCOMPLEX
#if 1 // A = (7,1) * D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            A0[k](i,j) = T(7,1) * D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]) = T(7,1) * D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_reg = 0.;
            for (int k=0; k<nloops2; ++k) e11_reg += NormSq(A1[k]-A0[k]);
            e11_reg = sqrt(e11_reg/nloops2);
            e11_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A2[k]) = T(7,1) * D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_small = 0.;
            for (int k=0; k<nloops2; ++k) e11_small += NormSq(A2[k]-A0[k]);
            e11_small = sqrt(e11_small/nloops2);
            e11_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            A3[k].setZero();
            BLASNAME(geru) (BLASCM 
                            M,N,z71,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                            BP(A3[k].ptr()),M);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_blas = 0.;
            for (int k=0; k<nloops2; ++k) e11_blas += NormSq(A3[k]-A0[k]);
            e11_blas = sqrt(e11_blas/nloops2);
            e11_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(A4[k]) = T(7,1) * D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e11_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e11_eigen = sqrt(e11_eigen/nloops2);
            e11_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(A5[k]) = T(7,1) * D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e11_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e11_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e11_smalleigen = sqrt(e11_smalleigen/nloops2);
            e11_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // B = (7,1) * D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            B0[k](j,i) = T(7,1) * D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B1[k].transpose()) = T(7,1) * D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_reg = 0.;
            for (int k=0; k<nloops2; ++k) e12_reg += NormSq(B1[k]-B0[k]);
            e12_reg = sqrt(e12_reg/nloops2);
            e12_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B2[k].transpose()) = T(7,1) * D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_small = 0.;
            for (int k=0; k<nloops2; ++k) e12_small += NormSq(B2[k]-B0[k]);
            e12_small = sqrt(e12_small/nloops2);
            e12_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            B3[k].setZero();
            BLASNAME(geru) (BLASCM 
                            N,M,z71,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                            BP(B3[k].ptr()),N);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_blas = 0.;
            for (int k=0; k<nloops2; ++k) e12_blas += NormSq(B3[k]-B0[k]);
            e12_blas = sqrt(e12_blas/nloops2);
            e12_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(B4[k].transpose()) = T(7,1) * D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e12_eigen += tmv::TMV_NORM(B4[k](j,i)-B0[k](j,i));
            e12_eigen = sqrt(e12_eigen/nloops2);
            e12_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(B5[k].transpose()) = T(7,1) * D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e12_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e12_smalleigen += tmv::TMV_NORM(B5[k](j,i)-B0[k](j,i));
            e12_smalleigen = sqrt(e12_smalleigen/nloops2);
            e12_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // A += (8,9) * D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            A0[k](i,j) += T(8,9) * D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]) += T(8,9) * D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_reg = 0.;
            for (int k=0; k<nloops2; ++k) e13_reg += NormSq(A1[k]-A0[k]);
            e13_reg = sqrt(e13_reg/nloops2);
            e13_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A2[k]) += T(8,9) * D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_small = 0.;
            for (int k=0; k<nloops2; ++k) e13_small += NormSq(A2[k]-A0[k]);
            e13_small = sqrt(e13_small/nloops2);
            e13_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(geru) (BLASCM 
                            M,N,z89,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                            BP(A3[k].ptr()),M);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_blas = 0.;
            for (int k=0; k<nloops2; ++k) e13_blas += NormSq(A3[k]-A0[k]);
            e13_blas = sqrt(e13_blas/nloops2);
            e13_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(A4[k]) += T(8,9) * D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e13_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e13_eigen = sqrt(e13_eigen/nloops2);
            e13_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(A5[k]) += T(8,9) * D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e13_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e13_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e13_smalleigen = sqrt(e13_smalleigen/nloops2);
            e13_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // B += (8,9) * D ^ C
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            B0[k](j,i) += T(8,9) * D0[k](i) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B1[k].transpose()) += T(8,9) * D1[k] ^ C1[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_reg = 0.;
            for (int k=0; k<nloops2; ++k) e14_reg += NormSq(B1[k]-B0[k]);
            e14_reg = sqrt(e14_reg/nloops2);
            e14_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B2[k].transpose()) += T(8,9) * D2[k] ^ C2[k];

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_small = 0.;
            for (int k=0; k<nloops2; ++k) e14_small += NormSq(B2[k]-B0[k]);
            e14_small = sqrt(e14_small/nloops2);
            e14_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(geru) (BLASCM 
                            N,M,z89,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                            BP(B3[k].ptr()),N);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_blas = 0.;
            for (int k=0; k<nloops2; ++k) e14_blas += NormSq(B3[k]-B0[k]);
            e14_blas = sqrt(e14_blas/nloops2);
            e14_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(B4[k].transpose()) += T(8,9) * D4[k] * C4[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e14_eigen += tmv::TMV_NORM(B4[k](j,i)-B0[k](j,i));
            e14_eigen = sqrt(e14_eigen/nloops2);
            e14_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(B5[k].transpose()) += T(8,9) * D5[k] * C5[k].transpose();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e14_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e14_smalleigen += tmv::TMV_NORM(B5[k](j,i)-B0[k](j,i));
            e14_smalleigen = sqrt(e14_smalleigen/nloops2);
            e14_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // A = (7,1) * D* ^ C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            A0[k](i,j) = T(7,1) * std::conj(D0[k](i) * C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]) = T(7,1) * D1[k].conjugate() ^ C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_reg = 0.;
            for (int k=0; k<nloops2; ++k) e15_reg += NormSq(A1[k]-A0[k]);
            e15_reg = sqrt(e15_reg/nloops2);
            e15_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A2[k]) = T(7,1) * D2[k].conjugate() ^ C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_small = 0.;
            for (int k=0; k<nloops2; ++k) e15_small += NormSq(A2[k]-A0[k]);
            e15_small = sqrt(e15_small/nloops2);
            e15_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            A3[k].setZero();
            BLASNAME(geru) (BLASCM 
                            M,N,z71c,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                            BP(A3[k].ptr()),M);
            BLASDNAME(scal) (M*N,dmone,A3[k].realPart().ptr()+1,2);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_blas = 0.;
            for (int k=0; k<nloops2; ++k) e15_blas += NormSq(A3[k]-A0[k]);
            e15_blas = sqrt(e15_blas/nloops2);
            e15_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(A4[k]) = T(7,1) * D4[k].conjugate() * C4[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e15_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e15_eigen = sqrt(e15_eigen/nloops2);
            e15_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(A5[k]) = T(7,1) * D5[k].conjugate() * C5[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e15_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e15_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e15_smalleigen = sqrt(e15_smalleigen/nloops2);
            e15_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // B = (7,1) * D* ^ C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            B0[k](j,i) = T(7,1) * std::conj(D0[k](i) * C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B1[k].transpose()) = T(7,1) * D1[k].conjugate() ^ C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_reg = 0.;
            for (int k=0; k<nloops2; ++k) e16_reg += NormSq(B1[k]-B0[k]);
            e16_reg = sqrt(e16_reg/nloops2);
            e16_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B2[k].transpose()) = T(7,1) * D2[k].conjugate() ^ C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_small = 0.;
            for (int k=0; k<nloops2; ++k) e16_small += NormSq(B2[k]-B0[k]);
            e16_small = sqrt(e16_small/nloops2);
            e16_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            B3[k].setZero();
            BLASNAME(geru) (BLASCM 
                            N,M,z71c,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                            BP(B3[k].ptr()),N);
            BLASDNAME(scal) (M*N,dmone,B3[k].realPart().ptr()+1,2);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_blas = 0.;
            for (int k=0; k<nloops2; ++k) e16_blas += NormSq(B3[k]-B0[k]);
            e16_blas = sqrt(e16_blas/nloops2);
            e16_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(B4[k].transpose()) = T(7,1) * D4[k].conjugate() * C4[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e16_eigen += tmv::TMV_NORM(B4[k](j,i)-B0[k](j,i));
            e16_eigen = sqrt(e16_eigen/nloops2);
            e16_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(B5[k].transpose()) = T(7,1) * D5[k].conjugate() * C5[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e16_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e16_smalleigen += tmv::TMV_NORM(B5[k](j,i)-B0[k](j,i));
            e16_smalleigen = sqrt(e16_smalleigen/nloops2);
            e16_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // A += (8,9) * D* ^ C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            A0[k](i,j) += T(8,9) * std::conj(D0[k](i) * C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]) += T(8,9) * D1[k].conjugate() ^ C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_reg = 0.;
            for (int k=0; k<nloops2; ++k) e17_reg += NormSq(A1[k]-A0[k]);
            e17_reg = sqrt(e17_reg/nloops2);
            e17_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A2[k]) += T(8,9) * D2[k].conjugate() ^ C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_small = 0.;
            for (int k=0; k<nloops2; ++k) e17_small += NormSq(A2[k]-A0[k]);
            e17_small = sqrt(e17_small/nloops2);
            e17_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASDNAME(scal) (M*N,dmone,A3[k].realPart().ptr()+1,2);
            BLASNAME(geru) (BLASCM 
                            M,N,z89c,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                            BP(A3[k].ptr()),M);
            BLASDNAME(scal) (M*N,dmone,A3[k].realPart().ptr()+1,2);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_blas = 0.;
            for (int k=0; k<nloops2; ++k) e17_blas += NormSq(A3[k]-A0[k]);
            e17_blas = sqrt(e17_blas/nloops2);
            e17_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(A4[k]) += T(8,9) * D4[k].conjugate() * C4[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e17_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e17_eigen = sqrt(e17_eigen/nloops2);
            e17_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(A5[k]) += T(8,9) * D5[k].conjugate() * C5[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e17_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e17_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e17_smalleigen = sqrt(e17_smalleigen/nloops2);
            e17_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // B += (8,9) * D* ^ C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            B0[k](j,i) += T(8,9) * std::conj(D0[k](i) * C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B1[k].transpose()) += T(8,9) * D1[k].conjugate() ^ C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_reg = 0.;
            for (int k=0; k<nloops2; ++k) e18_reg += NormSq(B1[k]-B0[k]);
            e18_reg = sqrt(e18_reg/nloops2);
            e18_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B2[k].transpose()) += T(8,9) * D2[k].conjugate() ^ C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_small = 0.;
            for (int k=0; k<nloops2; ++k) e18_small += NormSq(B2[k]-B0[k]);
            e18_small = sqrt(e18_small/nloops2);
            e18_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASDNAME(scal) (M*N,dmone,B3[k].realPart().ptr()+1,2);
            BLASNAME(geru) (BLASCM 
                            N,M,z89c,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                            BP(B3[k].ptr()),N);
            BLASDNAME(scal) (M*N,dmone,B3[k].realPart().ptr()+1,2);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_blas = 0.;
            for (int k=0; k<nloops2; ++k) e18_blas += NormSq(B3[k]-B0[k]);
            e18_blas = sqrt(e18_blas/nloops2);
            e18_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(B4[k].transpose()) += T(8,9) * D4[k].conjugate() * C4[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e18_eigen += tmv::TMV_NORM(B4[k](j,i)-B0[k](j,i));
            e18_eigen = sqrt(e18_eigen/nloops2);
            e18_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(B5[k].transpose()) += T(8,9) * D5[k].conjugate() * C5[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e18_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e18_smalleigen += tmv::TMV_NORM(B5[k](j,i)-B0[k](j,i));
            e18_smalleigen = sqrt(e18_smalleigen/nloops2);
            e18_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // A = (7,1) * D ^ C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            A0[k](i,j) = T(7,1) * D0[k](i) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]) = T(7,1) * D1[k] ^ C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_reg = 0.;
            for (int k=0; k<nloops2; ++k) e19_reg += NormSq(A1[k]-A0[k]);
            e19_reg = sqrt(e19_reg/nloops2);
            e19_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A2[k]) = T(7,1) * D2[k] ^ C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_small = 0.;
            for (int k=0; k<nloops2; ++k) e19_small += NormSq(A2[k]-A0[k]);
            e19_small = sqrt(e19_small/nloops2);
            e19_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            A3[k].setZero();
            BLASNAME(gerc) (BLASCM 
                            M,N,z71,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                            BP(A3[k].ptr()),M);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_blas = 0.;
            for (int k=0; k<nloops2; ++k) e19_blas += NormSq(A3[k]-A0[k]);
            e19_blas = sqrt(e19_blas/nloops2);
            e19_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(A4[k]) = T(7,1) * D4[k] * C4[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e19_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e19_eigen = sqrt(e19_eigen/nloops2);
            e19_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(A5[k]) = T(7,1) * D5[k] * C5[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e19_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e19_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e19_smalleigen = sqrt(e19_smalleigen/nloops2);
            e19_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // B = (7,1) * D ^ C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            B0[k](j,i) = T(7,1) * D0[k](i) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B1[k].transpose()) = T(7,1) * D1[k] ^ C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_reg = 0.;
            for (int k=0; k<nloops2; ++k) e20_reg += NormSq(B1[k]-B0[k]);
            e20_reg = sqrt(e20_reg/nloops2);
            e20_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B2[k].transpose()) = T(7,1) * D2[k] ^ C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_small = 0.;
            for (int k=0; k<nloops2; ++k) e20_small += NormSq(B2[k]-B0[k]);
            e20_small = sqrt(e20_small/nloops2);
            e20_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            B3[k].setZero();
            BLASNAME(gerc) (BLASCM 
                            N,M,z71c,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                            BP(B3[k].ptr()),N);
            BLASDNAME(scal) (M*N,dmone,B3[k].realPart().ptr()+1,2);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_blas = 0.;
            for (int k=0; k<nloops2; ++k) e20_blas += NormSq(B3[k]-B0[k]);
            e20_blas = sqrt(e20_blas/nloops2);
            e20_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(B4[k].transpose()) = T(7,1) * D4[k] * C4[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e20_eigen += tmv::TMV_NORM(B4[k](j,i)-B0[k](j,i));
            e20_eigen = sqrt(e20_eigen/nloops2);
            e20_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(B5[k].transpose()) = T(7,1) * D5[k] * C5[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e20_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e20_smalleigen += tmv::TMV_NORM(B5[k](j,i)-B0[k](j,i));
            e20_smalleigen = sqrt(e20_smalleigen/nloops2);
            e20_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // A += (8,9) * D ^ C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            A0[k](i,j) += T(8,9) * D0[k](i) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A1[k]) += T(8,9) * D1[k] ^ C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_reg = 0.;
            for (int k=0; k<nloops2; ++k) e21_reg += NormSq(A1[k]-A0[k]);
            e21_reg = sqrt(e21_reg/nloops2);
            e21_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(A2[k]) += T(8,9) * D2[k] ^ C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_small = 0.;
            for (int k=0; k<nloops2; ++k) e21_small += NormSq(A2[k]-A0[k]);
            e21_small = sqrt(e21_small/nloops2);
            e21_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASNAME(gerc) (BLASCM 
                            M,N,z89,BP(D3[k].cptr()),1,BP(C3[k].cptr()),1,
                            BP(A3[k].ptr()),M);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_blas = 0.;
            for (int k=0; k<nloops2; ++k) e21_blas += NormSq(A3[k]-A0[k]);
            e21_blas = sqrt(e21_blas/nloops2);
            e21_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(A4[k]) += T(8,9) * D4[k] * C4[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e21_eigen += tmv::TMV_NORM(A4[k](i,j)-A0[k](i,j));
            e21_eigen = sqrt(e21_eigen/nloops2);
            e21_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(A5[k]) += T(8,9) * D5[k] * C5[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e21_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e21_smalleigen += tmv::TMV_NORM(A5[k](i,j)-A0[k](i,j));
            e21_smalleigen = sqrt(e21_smalleigen/nloops2);
            e21_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // B += (8,9) * D ^ C*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            B0[k](j,i) += T(8,9) * D0[k](i) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B1[k].transpose()) += T(8,9) * D1[k] ^ C1[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_reg = 0.;
            for (int k=0; k<nloops2; ++k) e22_reg += NormSq(B1[k]-B0[k]);
            e22_reg = sqrt(e22_reg/nloops2);
            e22_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART(B2[k].transpose()) += T(8,9) * D2[k] ^ C2[k].conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_small = 0.;
            for (int k=0; k<nloops2; ++k) e22_small += NormSq(B2[k]-B0[k]);
            e22_small = sqrt(e22_small/nloops2);
            e22_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            BLASDNAME(scal) (M*N,dmone,B3[k].realPart().ptr()+1,2);
            BLASNAME(gerc) (BLASCM 
                            N,M,z89c,BP(C3[k].cptr()),1,BP(D3[k].cptr()),1,
                            BP(B3[k].ptr()),N);
            BLASDNAME(scal) (M*N,dmone,B3[k].realPart().ptr()+1,2);
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_blas = 0.;
            for (int k=0; k<nloops2; ++k) e22_blas += NormSq(B3[k]-B0[k]);
            e22_blas = sqrt(e22_blas/nloops2);
            e22_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART(B4[k].transpose()) += T(8,9) * D4[k] * C4[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e22_eigen += tmv::TMV_NORM(B4[k](j,i)-B0[k](j,i));
            e22_eigen = sqrt(e22_eigen/nloops2);
            e22_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART(B5[k].transpose()) += T(8,9) * D5[k] * C5[k].adjoint();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e22_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for (int j=0;j<N;++j)
                    e22_smalleigen += tmv::TMV_NORM(B5[k](j,i)-B0[k](j,i));
            e22_smalleigen = sqrt(e22_smalleigen/nloops2);
            e22_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif
#endif

    }

    std::cout<<"A = D ^ C               "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"B = D ^ C               "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"A = -A ^ C              "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"B = -D ^ C              "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"A = 7 * D ^ C           "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
    std::cout<<"B = 7 * D ^ C           "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"A -= D ^ C              "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"B -= D ^ C              "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
    std::cout<<"A += 8 * D ^ C          "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
    std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
    std::cout<<"B += 8 * D ^ C          "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
    std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"A = (7,1) * D ^ C       "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
    std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
    std::cout<<"B = (7,1) * D ^ C       "<<t12_reg<<"  "<<t12_small<<"  "<<t12_blas;
    std::cout<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;
    std::cout<<"A += (8,9) * D ^ C      "<<t13_reg<<"  "<<t13_small<<"  "<<t13_blas;
    std::cout<<"  "<<t13_eigen<<"  "<<t13_smalleigen<<std::endl;
    std::cout<<"B += (8,9) * D ^ C      "<<t14_reg<<"  "<<t14_small<<"  "<<t14_blas;
    std::cout<<"  "<<t14_eigen<<"  "<<t14_smalleigen<<std::endl;
    std::cout<<"A = (7,1) * D ^ C*      "<<t19_reg<<"  "<<t19_small<<"  "<<t19_blas;
    std::cout<<"  "<<t19_eigen<<"  "<<t19_smalleigen<<std::endl;
    std::cout<<"B = (7,1) * D ^ C*      "<<t20_reg<<"  "<<t20_small<<"  "<<t20_blas;
    std::cout<<"  "<<t20_eigen<<"  "<<t20_smalleigen<<std::endl;
    std::cout<<"A += (8,9) * D ^ C*     "<<t21_reg<<"  "<<t21_small<<"  "<<t21_blas;
    std::cout<<"  "<<t21_eigen<<"  "<<t21_smalleigen<<std::endl;
    std::cout<<"B += (8,9) * D ^ C*     "<<t22_reg<<"  "<<t22_small<<"  "<<t22_blas;
    std::cout<<"  "<<t22_eigen<<"  "<<t22_smalleigen<<std::endl;
    std::cout<<"A = (7,1) * D* ^ C*     "<<t15_reg<<"  "<<t15_small<<"  "<<t15_blas;
    std::cout<<"  "<<t15_eigen<<"  "<<t15_smalleigen<<std::endl;
    std::cout<<"B = (7,1) * D* ^ C*     "<<t16_reg<<"  "<<t16_small<<"  "<<t16_blas;
    std::cout<<"  "<<t16_eigen<<"  "<<t16_smalleigen<<std::endl;
    std::cout<<"A += (8,9) * D* ^ C*    "<<t17_reg<<"  "<<t17_small<<"  "<<t17_blas;
    std::cout<<"  "<<t17_eigen<<"  "<<t17_smalleigen<<std::endl;
    std::cout<<"B += (8,9) * D* ^ C*    "<<t18_reg<<"  "<<t18_small<<"  "<<t18_blas;
    std::cout<<"  "<<t18_eigen<<"  "<<t18_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"A = D ^ C               "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"B = D ^ C               "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"A = -D ^ C              "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"B = -D ^ C              "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"A = 7 * D ^ C           "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
    std::cout<<"B = 7 * D ^ C           "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
    std::cout<<"A -= D ^ C              "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
    std::cout<<"B -= D ^ C              "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
    std::cout<<"A += 8 * D ^ C          "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
    std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
    std::cout<<"B += 8 * D ^ C          "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
    std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"A = (7,1) * D ^ C       "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
    std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
    std::cout<<"B = (7,1) * D ^ C       "<<e12_reg<<"  "<<e12_small<<"  "<<e12_blas;
    std::cout<<"  "<<e12_eigen<<"  "<<e12_smalleigen<<std::endl;
    std::cout<<"A += (8,9) * D ^ C      "<<e13_reg<<"  "<<e13_small<<"  "<<e13_blas;
    std::cout<<"  "<<e13_eigen<<"  "<<e13_smalleigen<<std::endl;
    std::cout<<"B += (8,9) * D ^ C      "<<e14_reg<<"  "<<e14_small<<"  "<<e14_blas;
    std::cout<<"  "<<e14_eigen<<"  "<<e14_smalleigen<<std::endl;
    std::cout<<"A = (7,1) * D ^ C*      "<<e19_reg<<"  "<<e19_small<<"  "<<e19_blas;
    std::cout<<"  "<<e19_eigen<<"  "<<e19_smalleigen<<std::endl;
    std::cout<<"B = (7,1) * D ^ C*      "<<e20_reg<<"  "<<e20_small<<"  "<<e20_blas;
    std::cout<<"  "<<e20_eigen<<"  "<<e20_smalleigen<<std::endl;
    std::cout<<"A += (8,9) * D ^ C*     "<<e21_reg<<"  "<<e21_small<<"  "<<e21_blas;
    std::cout<<"  "<<e21_eigen<<"  "<<e21_smalleigen<<std::endl;
    std::cout<<"B += (8,9) * D ^ C*     "<<e22_reg<<"  "<<e22_small<<"  "<<e22_blas;
    std::cout<<"  "<<e22_eigen<<"  "<<e22_smalleigen<<std::endl;
    std::cout<<"A = (7,1) * D* ^ C*     "<<e15_reg<<"  "<<e15_small<<"  "<<e15_blas;
    std::cout<<"  "<<e15_eigen<<"  "<<e15_smalleigen<<std::endl;
    std::cout<<"B = (7,1) * D* ^ C*     "<<e16_reg<<"  "<<e16_small<<"  "<<e16_blas;
    std::cout<<"  "<<e16_eigen<<"  "<<e16_smalleigen<<std::endl;
    std::cout<<"A += (8,9) * D* ^ C*    "<<e17_reg<<"  "<<e17_small<<"  "<<e17_blas;
    std::cout<<"  "<<e17_eigen<<"  "<<e17_smalleigen<<std::endl;
    std::cout<<"B += (8,9) * D* ^ C*    "<<e18_reg<<"  "<<e18_small<<"  "<<e18_blas;
    std::cout<<"  "<<e18_eigen<<"  "<<e18_smalleigen<<std::endl;
#endif
    std::cout<<"\n\n";
#endif
}
#endif

#ifdef DOMULTMD
static void MultMD(
    const std::vector<tmv::Matrix<T> >& A1,
    const std::vector<tmv::Matrix<T> >& B1,
    const std::vector<tmv::Vector<T> >& C1)
{
    std::vector<tmv::Matrix<T> > D1 = A1;

#ifdef ERRORCHECK
    std::vector<tmv::Matrix<T> > A0 = A1;
    std::vector<tmv::Matrix<T> > B0 = B1;
    std::vector<tmv::Vector<T> > C0 = C1;
    std::vector<tmv::Matrix<T> > D0 = D1;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallMatrix<T,M,N> > A2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,M> > B2(nloops2);
    std::vector<tmv::SmallVector<T,N> > C2(nloops2);
    std::vector<tmv::SmallMatrix<T,M,N> > D2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; B2[k] = B1[k]; C2[k] = C1[k]; D2[k] = D1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Matrix<T> > A3 = A1;
    std::vector<tmv::Matrix<T> > B3 = B1;
    std::vector<tmv::Vector<T> > C3 = C1;
    std::vector<tmv::Matrix<T> > D3 = D1;
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENM> A4(nloops2,Eigen::EIGENM(M,N));
    std::vector<Eigen::EIGENM> B4(nloops2,Eigen::EIGENM(N,M));
    std::vector<Eigen::EIGENV> C4(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENM> D4(nloops2,Eigen::EIGENM(M,N));
    for(int k=0;k<nloops2;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) A4[k](i,j) = A1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) B4[k](j,i) = B1[k](j,i);
        for(int i=0;i<N;++i) C4[k](i) = C1[k](i);
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) D4[k](i,j) = D1[k](i,j);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,M,N> > A5;
    std::vector<Eigen::Matrix<T,N,M> > B5;
    std::vector<Eigen::Matrix<T,N,1> > C5;
    std::vector<Eigen::Matrix<T,M,N> > D5;
    if (nloops2x) { 
        A5.resize(nloops2x); B5.resize(nloops2x);
        C5.resize(nloops2x); D5.resize(nloops2x);
    }
    for(int k=0;k<nloops2x;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) A5[k](i,j) = A1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) B5[k](j,i) = B1[k](j,i);
        for(int i=0;i<N;++i) C5[k](i) = C1[k](i);
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) D5[k](i,j) = D1[k](i,j);
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
    double t17_reg=0., t17_small=0., t17_blas=0., t17_eigen=0., t17_smalleigen=0.;
    double t18_reg=0., t18_small=0., t18_blas=0., t18_eigen=0., t18_smalleigen=0.;
    double t19_reg=0., t19_small=0., t19_blas=0., t19_eigen=0., t19_smalleigen=0.;
    double t20_reg=0., t20_small=0., t20_blas=0., t20_eigen=0., t20_smalleigen=0.;
    double t21_reg=0., t21_small=0., t21_blas=0., t21_eigen=0., t21_smalleigen=0.;
    double t22_reg=0., t22_small=0., t22_blas=0., t22_eigen=0., t22_smalleigen=0.;
    double t23_reg=0., t23_small=0., t23_blas=0., t23_eigen=0., t23_smalleigen=0.;
    double t24_reg=0., t24_small=0., t24_blas=0., t24_eigen=0., t24_smalleigen=0.;
    double t25_reg=0., t25_small=0., t25_blas=0., t25_eigen=0., t25_smalleigen=0.;
    double t26_reg=0., t26_small=0., t26_blas=0., t26_eigen=0., t26_smalleigen=0.;
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
    double e17_reg=0., e17_small=0., e17_blas=0., e17_eigen=0., e17_smalleigen=0.;
    double e18_reg=0., e18_small=0., e18_blas=0., e18_eigen=0., e18_smalleigen=0.;
    double e19_reg=0., e19_small=0., e19_blas=0., e19_eigen=0., e19_smalleigen=0.;
    double e20_reg=0., e20_small=0., e20_blas=0., e20_eigen=0., e20_smalleigen=0.;
    double e21_reg=0., e21_small=0., e21_blas=0., e21_eigen=0., e21_smalleigen=0.;
    double e22_reg=0., e22_small=0., e22_blas=0., e22_eigen=0., e22_smalleigen=0.;
    double e23_reg=0., e23_small=0., e23_blas=0., e23_eigen=0., e23_smalleigen=0.;
    double e24_reg=0., e24_small=0., e24_blas=0., e24_eigen=0., e24_smalleigen=0.;
    double e25_reg=0., e25_small=0., e25_blas=0., e25_eigen=0., e25_smalleigen=0.;
    double e26_reg=0., e26_small=0., e26_blas=0., e26_eigen=0., e26_smalleigen=0.;
#endif

    for (int n=0; n<nloops1; ++n) {

#if 1 // D = A * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = MPART(A1[k]) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(D1[k]-D0[k]);
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = MPART(A2[k]) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(D2[k]-D0[k]);
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
#endif
            for(int j=0;j<N;++j) {
                T x = C3[k](j);
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(D3[k]-D0[k]);
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = EPART(A4[k]) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e1_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = EPART(A5[k]) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e1_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = B * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) 
            MPART2(D1[k]) = MPART(B1[k].transpose()) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_reg += NormSq(D1[k]-D0[k]);
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = MPART(B2[k].transpose()) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) e2_small += NormSq(D2[k]-D0[k]);
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                T x = C3[k](j);
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) e2_blas += NormSq(D3[k]-D0[k]);
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = EPART(B4[k].transpose()) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e2_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = EPART(B5[k].transpose()) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e2_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = -A * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = -A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = -C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = -MPART(A1[k]) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) e3_reg += NormSq(D1[k]-D0[k]);
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = -MPART(A2[k]) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) e3_small += NormSq(D2[k]-D0[k]);
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
#endif
            for(int j=0;j<N;++j) {
                T x = -C3[k](j);
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,mone,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) e3_blas += NormSq(D3[k]-D0[k]);
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = -EPART(A4[k]) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e3_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = -EPART(A5[k]) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e3_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = -B * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = -B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = -C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = -MPART(B1[k].transpose()) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) e4_reg += NormSq(D1[k]-D0[k]);
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = -MPART(B2[k].transpose()) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) e4_small += NormSq(D2[k]-D0[k]);
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                T x = -C3[k](j);
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,mone,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) e4_blas += NormSq(D3[k]-D0[k]);
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = -EPART(B4[k].transpose()) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e4_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = -EPART(B5[k].transpose()) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e4_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = 7 * A * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = RT(7) * A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = RT(7) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = 7 * MPART(A1[k]) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) e5_reg += NormSq(D1[k]-D0[k]);
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = 7 * MPART(A2[k]) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) e5_small += NormSq(D2[k]-D0[k]);
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
#endif
            for(int j=0;j<N;++j) {
                T x = T(7)*C3[k](j);
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,seven,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) e5_blas += NormSq(D3[k]-D0[k]);
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = 7 * EPART(A4[k]) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e5_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = 7 * EPART(A5[k]) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e5_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = 7 * B * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = RT(7) * B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = RT(7) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = 7 * MPART(B1[k].transpose()) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) e6_reg += NormSq(D1[k]-D0[k]);
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = 7 * MPART(B2[k].transpose()) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) e6_small += NormSq(D2[k]-D0[k]);
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                T x = T(7) * C3[k](j);
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,seven,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) e6_blas += NormSq(D3[k]-D0[k]);
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = 7 * EPART(B4[k].transpose()) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e6_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = 7 * EPART(B5[k].transpose()) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e6_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D -= A * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) -= A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) -= C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) -= MPART(A1[k]) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) e7_reg += NormSq(D1[k]-D0[k]);
            e7_reg = sqrt(e7_reg/nloops2);
            e7_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) -= MPART(A2[k]) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) e7_small += NormSq(D2[k]-D0[k]);
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                tmv::Vector<T> temp(COL_LEN(j));
                T x = -C3[k](j);
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    COL_LEN(j),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,mone,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) e7_blas += NormSq(D3[k]-D0[k]);
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) -= EPART(A4[k]) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e7_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) -= EPART(A5[k]) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e7_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D -= B * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) -= B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) -= C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) -= MPART(B1[k].transpose()) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) e8_reg += NormSq(D1[k]-D0[k]);
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) -= MPART(B2[k].transpose()) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) e8_small += NormSq(D2[k]-D0[k]);
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                tmv::Vector<T> temp(COL_LEN(j));
                T x = -C3[k](j);
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    COL_LEN(j),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,mone,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) e8_blas += NormSq(D3[k]-D0[k]);
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) -= EPART(B4[k].transpose()) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e8_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) -= EPART(B5[k].transpose()) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e8_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += 8 * A * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += RT(8) * A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) += RT(8) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += 8 * MPART(A1[k]) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_reg = 0.;
            for (int k=0; k<nloops2; ++k) e9_reg += NormSq(D1[k]-D0[k]);
            e9_reg = sqrt(e9_reg/nloops2);
            e9_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += 8 * MPART(A2[k]) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_small = 0.;
            for (int k=0; k<nloops2; ++k) e9_small += NormSq(D2[k]-D0[k]);
            e9_small = sqrt(e9_small/nloops2);
            e9_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                tmv::Vector<T> temp(COL_LEN(j));
                T x = T(8) * C3[k](j);
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    COL_LEN(j),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,eight,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_blas = 0.;
            for (int k=0; k<nloops2; ++k) e9_blas += NormSq(D3[k]-D0[k]);
            e9_blas = sqrt(e9_blas/nloops2);
            e9_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += 8 * EPART(A4[k]) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e9_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e9_eigen = sqrt(e9_eigen/nloops2);
            e9_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += 8 * EPART(A5[k]) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e9_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e9_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e9_smalleigen = sqrt(e9_smalleigen/nloops2);
            e9_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += 8 * B * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += RT(8) * B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) += RT(8) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += 8 * MPART(B1[k].transpose()) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_reg = 0.;
            for (int k=0; k<nloops2; ++k) e10_reg += NormSq(D1[k]-D0[k]);
            e10_reg = sqrt(e10_reg/nloops2);
            e10_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += 8 * MPART(B2[k].transpose()) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_small = 0.;
            for (int k=0; k<nloops2; ++k) e10_small += NormSq(D2[k]-D0[k]);
            e10_small = sqrt(e10_small/nloops2);
            e10_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                tmv::Vector<T> temp(COL_LEN(j));
                T x = T(8) * C3[k](j);
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    COL_LEN(j),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,eight,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_blas = 0.;
            for (int k=0; k<nloops2; ++k) e10_blas += NormSq(D3[k]-D0[k]);
            e10_blas = sqrt(e10_blas/nloops2);
            e10_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += 8 * EPART(B4[k].transpose()) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e10_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e10_eigen = sqrt(e10_eigen/nloops2);
            e10_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += 8 * EPART(B5[k].transpose()) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e10_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e10_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e10_smalleigen = sqrt(e10_smalleigen/nloops2);
            e10_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#ifdef TISCOMPLEX
#if 1 // D = (7,1) * A * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * MPART(A1[k]) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_reg = 0.;
            for (int k=0; k<nloops2; ++k) e11_reg += NormSq(D1[k]-D0[k]);
            e11_reg = sqrt(e11_reg/nloops2);
            e11_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * MPART(A2[k]) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_small = 0.;
            for (int k=0; k<nloops2; ++k) e11_small += NormSq(D2[k]-D0[k]);
            e11_small = sqrt(e11_small/nloops2);
            e11_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
#endif
            for(int j=0;j<N;++j) {
                T x = T(7,1)*C3[k](j);
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_blas = 0.;
            for (int k=0; k<nloops2; ++k) e11_blas += NormSq(D3[k]-D0[k]);
            e11_blas = sqrt(e11_blas/nloops2);
            e11_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * EPART(A4[k]) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e11_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e11_eigen = sqrt(e11_eigen/nloops2);
            e11_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * EPART(A5[k]) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e11_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e11_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e11_smalleigen = sqrt(e11_smalleigen/nloops2);
            e11_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * B * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * MPART(B1[k].transpose()) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_reg = 0.;
            for (int k=0; k<nloops2; ++k) e12_reg += NormSq(D1[k]-D0[k]);
            e12_reg = sqrt(e12_reg/nloops2);
            e12_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * MPART(B2[k].transpose()) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_small = 0.;
            for (int k=0; k<nloops2; ++k) e12_small += NormSq(D2[k]-D0[k]);
            e12_small = sqrt(e12_small/nloops2);
            e12_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                T x = T(7,1) * C3[k](j);
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_blas = 0.;
            for (int k=0; k<nloops2; ++k) e12_blas += NormSq(D3[k]-D0[k]);
            e12_blas = sqrt(e12_blas/nloops2);
            e12_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * EPART(B4[k].transpose()) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e12_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e12_eigen = sqrt(e12_eigen/nloops2);
            e12_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * EPART(B5[k].transpose()) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e12_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e12_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e12_smalleigen = sqrt(e12_smalleigen/nloops2);
            e12_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * A * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * A0[k](i,j) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * MPART(A1[k]) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_reg = 0.;
            for (int k=0; k<nloops2; ++k) e13_reg += NormSq(D1[k]-D0[k]);
            e13_reg = sqrt(e13_reg/nloops2);
            e13_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * MPART(A2[k]) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_small = 0.;
            for (int k=0; k<nloops2; ++k) e13_small += NormSq(D2[k]-D0[k]);
            e13_small = sqrt(e13_small/nloops2);
            e13_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                tmv::Vector<T> temp(COL_LEN(j));
                T x = T(8,9) * C3[k](j);
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    COL_LEN(j),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_blas = 0.;
            for (int k=0; k<nloops2; ++k) e13_blas += NormSq(D3[k]-D0[k]);
            e13_blas = sqrt(e13_blas/nloops2);
            e13_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * EPART(A4[k]) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e13_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e13_eigen = sqrt(e13_eigen/nloops2);
            e13_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * EPART(A5[k]) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e13_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e13_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e13_smalleigen = sqrt(e13_smalleigen/nloops2);
            e13_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * B * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * B0[k](j,i) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * MPART(B1[k].transpose()) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_reg = 0.;
            for (int k=0; k<nloops2; ++k) e14_reg += NormSq(D1[k]-D0[k]);
            e14_reg = sqrt(e14_reg/nloops2);
            e14_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * MPART(B2[k].transpose()) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_small = 0.;
            for (int k=0; k<nloops2; ++k) e14_small += NormSq(D2[k]-D0[k]);
            e14_small = sqrt(e14_small/nloops2);
            e14_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                tmv::Vector<T> temp(COL_LEN(j));
                T x = T(8,9) * C3[k](j);
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    COL_LEN(j),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_blas = 0.;
            for (int k=0; k<nloops2; ++k) e14_blas += NormSq(D3[k]-D0[k]);
            e14_blas = sqrt(e14_blas/nloops2);
            e14_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * EPART(B4[k].transpose()) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e14_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e14_eigen = sqrt(e14_eigen/nloops2);
            e14_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * EPART(B5[k].transpose()) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e14_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e14_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e14_smalleigen = sqrt(e14_smalleigen/nloops2);
            e14_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * A* * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(A0[k](i,j)) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * MPART(A1[k]).conjugate() * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_reg = 0.;
            for (int k=0; k<nloops2; ++k) e15_reg += NormSq(D1[k]-D0[k]);
            e15_reg = sqrt(e15_reg/nloops2);
            e15_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * MPART(A2[k]).conjugate() * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_small = 0.;
            for (int k=0; k<nloops2; ++k) e15_small += NormSq(D2[k]-D0[k]);
            e15_small = sqrt(e15_small/nloops2);
            e15_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASDNAME(scal) (M*N,dmone,D3[k].realPart().ptr()+1,2);
#endif
            for(int j=0;j<N;++j) {
                T x = T(7,1)*C3[k](j);
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASDNAME(scal) (
                    COL_LEN(j),dmone,
                    D3[k].col(j,COL_START(j),COL_END(j)).realPart().ptr()+1,2);
#endif
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_blas = 0.;
            for (int k=0; k<nloops2; ++k) e15_blas += NormSq(D3[k]-D0[k]);
            e15_blas = sqrt(e15_blas/nloops2);
            e15_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * EPART(A4[k]).conjugate() * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e15_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e15_eigen = sqrt(e15_eigen/nloops2);
            e15_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * EPART(A5[k]).conjugate() * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e15_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e15_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e15_smalleigen = sqrt(e15_smalleigen/nloops2);
            e15_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * B* * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(B0[k](j,i)) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * MPART(B1[k].adjoint()) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_reg = 0.;
            for (int k=0; k<nloops2; ++k) e17_reg += NormSq(D1[k]-D0[k]);
            e17_reg = sqrt(e17_reg/nloops2);
            e17_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * MPART(B2[k].adjoint()) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_small = 0.;
            for (int k=0; k<nloops2; ++k) e17_small += NormSq(D2[k]-D0[k]);
            e17_small = sqrt(e17_small/nloops2);
            e17_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                T x = T(7,1) * C3[k](j);
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASDNAME(scal) (
                    COL_LEN(j),dmone,
                    D3[k].col(j,COL_START(j),COL_END(j)).realPart().ptr()+1,2);
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_blas = 0.;
            for (int k=0; k<nloops2; ++k) e17_blas += NormSq(D3[k]-D0[k]);
            e17_blas = sqrt(e17_blas/nloops2);
            e17_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * EPART(B4[k].adjoint()) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e17_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e17_eigen = sqrt(e17_eigen/nloops2);
            e17_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * EPART(B5[k].adjoint()) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e17_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e17_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e17_smalleigen = sqrt(e17_smalleigen/nloops2);
            e17_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * A* * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(A0[k](i,j)) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * MPART(A1[k]).conjugate() * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_reg = 0.;
            for (int k=0; k<nloops2; ++k) e16_reg += NormSq(D1[k]-D0[k]);
            e16_reg = sqrt(e16_reg/nloops2);
            e16_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * MPART(A2[k]).conjugate() * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_small = 0.;
            for (int k=0; k<nloops2; ++k) e16_small += NormSq(D2[k]-D0[k]);
            e16_small = sqrt(e16_small/nloops2);
            e16_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                tmv::Vector<T> temp(COL_LEN(j));
                T x = T(8,9) * C3[k](j);
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(temp.ptr()),1);
                BLASDNAME(scal) (COL_LEN(j),dmone,temp.realPart().ptr()+1,2);
                BLASNAME(axpy) (
                    COL_LEN(j),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_blas = 0.;
            for (int k=0; k<nloops2; ++k) e16_blas += NormSq(D3[k]-D0[k]);
            e16_blas = sqrt(e16_blas/nloops2);
            e16_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * EPART(A4[k]).conjugate() * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e16_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e16_eigen = sqrt(e16_eigen/nloops2);
            e16_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * EPART(A5[k]).conjugate() * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e16_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e16_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e16_smalleigen = sqrt(e16_smalleigen/nloops2);
            e16_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * B* * diag(C)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(B0[k](j,i)) * C0[k](j);
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * MPART(B1[k].adjoint()) * DiagMatrixViewOf(C1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_reg = 0.;
            for (int k=0; k<nloops2; ++k) e18_reg += NormSq(D1[k]-D0[k]);
            e18_reg = sqrt(e18_reg/nloops2);
            e18_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * MPART(B2[k].adjoint()) * DiagMatrixViewOf(C2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_small = 0.;
            for (int k=0; k<nloops2; ++k) e18_small += NormSq(D2[k]-D0[k]);
            e18_small = sqrt(e18_small/nloops2);
            e18_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                tmv::Vector<T> temp(COL_LEN(j));
                T x = T(8,9) * C3[k](j);
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(temp.ptr()),1);
                BLASDNAME(scal) (COL_LEN(j),dmone,temp.realPart().ptr()+1,2);
                BLASNAME(axpy) (
                    COL_LEN(j),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_blas = 0.;
            for (int k=0; k<nloops2; ++k) e18_blas += NormSq(D3[k]-D0[k]);
            e18_blas = sqrt(e18_blas/nloops2);
            e18_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * EPART(B4[k].adjoint()) * C4[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e18_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e18_eigen = sqrt(e18_eigen/nloops2);
            e18_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * EPART(B5[k].adjoint()) * C5[k].asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e18_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e18_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e18_smalleigen = sqrt(e18_smalleigen/nloops2);
            e18_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * A * diag(C*)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * A0[k](i,j) * std::conj(C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * MPART(A1[k]) * DiagMatrixViewOf(C1[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_reg = 0.;
            for (int k=0; k<nloops2; ++k) e19_reg += NormSq(D1[k]-D0[k]);
            e19_reg = sqrt(e19_reg/nloops2);
            e19_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * MPART(A2[k]) * DiagMatrixViewOf(C2[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_small = 0.;
            for (int k=0; k<nloops2; ++k) e19_small += NormSq(D2[k]-D0[k]);
            e19_small = sqrt(e19_small/nloops2);
            e19_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
#endif
            for(int j=0;j<N;++j) {
                T x = T(7,1)*std::conj(C3[k](j));
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2*(N+1));
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_blas = 0.;
            for (int k=0; k<nloops2; ++k) e19_blas += NormSq(D3[k]-D0[k]);
            e19_blas = sqrt(e19_blas/nloops2);
            e19_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * EPART(A4[k]) * C4[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e19_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e19_eigen = sqrt(e19_eigen/nloops2);
            e19_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * EPART(A5[k]) * C5[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e19_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e19_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e19_smalleigen = sqrt(e19_smalleigen/nloops2);
            e19_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * B * diag(C*)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * B0[k](j,i) * std::conj(C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * MPART(B1[k].transpose()) * DiagMatrixViewOf(C1[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_reg = 0.;
            for (int k=0; k<nloops2; ++k) e20_reg += NormSq(D1[k]-D0[k]);
            e20_reg = sqrt(e20_reg/nloops2);
            e20_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * MPART(B2[k].transpose()) * DiagMatrixViewOf(C2[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_small = 0.;
            for (int k=0; k<nloops2; ++k) e20_small += NormSq(D2[k]-D0[k]);
            e20_small = sqrt(e20_small/nloops2);
            e20_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                T x = T(7,1) * std::conj(C3[k](j));
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2*(N+1));
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_blas = 0.;
            for (int k=0; k<nloops2; ++k) e20_blas += NormSq(D3[k]-D0[k]);
            e20_blas = sqrt(e20_blas/nloops2);
            e20_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * EPART(B4[k].transpose()) * C4[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e20_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e20_eigen = sqrt(e20_eigen/nloops2);
            e20_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * EPART(B5[k].transpose()) * C5[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e20_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e20_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e20_smalleigen = sqrt(e20_smalleigen/nloops2);
            e20_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * A * diag(C*)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * A0[k](i,j) * std::conj(C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * MPART(A1[k]) * DiagMatrixViewOf(C1[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_reg = 0.;
            for (int k=0; k<nloops2; ++k) e21_reg += NormSq(D1[k]-D0[k]);
            e21_reg = sqrt(e21_reg/nloops2);
            e21_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * MPART(A2[k]) * DiagMatrixViewOf(C2[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_small = 0.;
            for (int k=0; k<nloops2; ++k) e21_small += NormSq(D2[k]-D0[k]);
            e21_small = sqrt(e21_small/nloops2);
            e21_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                tmv::Vector<T> temp(COL_LEN(j));
                T x = T(8,9) * std::conj(C3[k](j));
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    COL_LEN(j),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_blas = 0.;
            for (int k=0; k<nloops2; ++k) e21_blas += NormSq(D3[k]-D0[k]);
            e21_blas = sqrt(e21_blas/nloops2);
            e21_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * EPART(A4[k]) * C4[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e21_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e21_eigen = sqrt(e21_eigen/nloops2);
            e21_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * EPART(A5[k]) * C5[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e21_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e21_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e21_smalleigen = sqrt(e21_smalleigen/nloops2);
            e21_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * B * diag(C*)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * B0[k](j,i) * std::conj(C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * MPART(B1[k].transpose()) * DiagMatrixViewOf(C1[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_reg = 0.;
            for (int k=0; k<nloops2; ++k) e22_reg += NormSq(D1[k]-D0[k]);
            e22_reg = sqrt(e22_reg/nloops2);
            e22_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * MPART(B2[k].transpose()) * DiagMatrixViewOf(C2[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_small = 0.;
            for (int k=0; k<nloops2; ++k) e22_small += NormSq(D2[k]-D0[k]);
            e22_small = sqrt(e22_small/nloops2);
            e22_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                tmv::Vector<T> temp(COL_LEN(j));
                T x = T(8,9) * std::conj(C3[k](j));
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    COL_LEN(j),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_blas = 0.;
            for (int k=0; k<nloops2; ++k) e22_blas += NormSq(D3[k]-D0[k]);
            e22_blas = sqrt(e22_blas/nloops2);
            e22_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * EPART(B4[k].transpose()) * C4[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e22_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e22_eigen = sqrt(e22_eigen/nloops2);
            e22_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * EPART(B5[k].transpose()) * C5[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e22_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e22_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e22_smalleigen = sqrt(e22_smalleigen/nloops2);
            e22_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * A* * diag(C*)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(A0[k](i,j) * C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * MPART(A1[k]).conjugate() * DiagMatrixViewOf(C1[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e23_reg = 0.;
            for (int k=0; k<nloops2; ++k) e23_reg += NormSq(D1[k]-D0[k]);
            e23_reg = sqrt(e23_reg/nloops2);
            e23_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * MPART(A2[k]).conjugate() * DiagMatrixViewOf(C2[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e23_small = 0.;
            for (int k=0; k<nloops2; ++k) e23_small += NormSq(D2[k]-D0[k]);
            e23_small = sqrt(e23_small/nloops2);
            e23_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
#endif
            for(int j=0;j<N;++j) {
                T x = T(7,-1)*C3[k](j);
#if ((PART >= 2) && (PART <= 5))
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#if ((PART >= 2) && (PART <= 5))
                BLASDNAME(scal) (
                    COL_LEN(j),dmone,
                    D3[k].col(j,COL_START(j),COL_END(j)).realPart().ptr()+1,2);
#endif
            }
#if (PART == 1)
            BLASDNAME(scal) (M*N,dmone,D3[k].realPart().ptr()+1,2);
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2*(N+1));
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e23_blas = 0.;
            for (int k=0; k<nloops2; ++k) e23_blas += NormSq(D3[k]-D0[k]);
            e23_blas = sqrt(e23_blas/nloops2);
            e23_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * EPART(A4[k]).conjugate() * C4[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e23_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e23_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e23_eigen = sqrt(e23_eigen/nloops2);
            e23_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * EPART(A5[k]).conjugate() * C5[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e23_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e23_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e23_smalleigen = sqrt(e23_smalleigen/nloops2);
            e23_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * B* * diag(C*)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(B0[k](j,i) * C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * MPART(B1[k].adjoint()) * DiagMatrixViewOf(C1[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e24_reg = 0.;
            for (int k=0; k<nloops2; ++k) e24_reg += NormSq(D1[k]-D0[k]);
            e24_reg = sqrt(e24_reg/nloops2);
            e24_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * MPART(B2[k].adjoint()) * DiagMatrixViewOf(C2[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e24_small = 0.;
            for (int k=0; k<nloops2; ++k) e24_small += NormSq(D2[k]-D0[k]);
            e24_small = sqrt(e24_small/nloops2);
            e24_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                T x = T(7,1) * std::conj(C3[k](j));
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASDNAME(scal) (
                    COL_LEN(j),dmone,
                    D3[k].col(j,COL_START(j),COL_END(j)).realPart().ptr()+1,2);
                BLASNAME(scal) (
                    COL_LEN(j),BX(x),
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2*(N+1));
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e24_blas = 0.;
            for (int k=0; k<nloops2; ++k) e24_blas += NormSq(D3[k]-D0[k]);
            e24_blas = sqrt(e24_blas/nloops2);
            e24_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * EPART(B4[k].adjoint()) * C4[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e24_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e24_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e24_eigen = sqrt(e24_eigen/nloops2);
            e24_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * EPART(B5[k].adjoint()) * C5[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e24_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e24_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e24_smalleigen = sqrt(e24_smalleigen/nloops2);
            e24_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * A* * diag(C*)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(A0[k](i,j) * C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * MPART(A1[k]).conjugate() * DiagMatrixViewOf(C1[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e25_reg = 0.;
            for (int k=0; k<nloops2; ++k) e25_reg += NormSq(D1[k]-D0[k]);
            e25_reg = sqrt(e25_reg/nloops2);
            e25_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * MPART(A2[k]).conjugate() * DiagMatrixViewOf(C2[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e25_small = 0.;
            for (int k=0; k<nloops2; ++k) e25_small += NormSq(D2[k]-D0[k]);
            e25_small = sqrt(e25_small/nloops2);
            e25_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                tmv::Vector<T> temp(COL_LEN(j));
                T x = T(8,9) * std::conj(C3[k](j));
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(temp.ptr()),1);
                BLASDNAME(scal) (COL_LEN(j),dmone,temp.realPart().ptr()+1,2);
                BLASNAME(axpy) (
                    COL_LEN(j),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e25_blas = 0.;
            for (int k=0; k<nloops2; ++k) e25_blas += NormSq(D3[k]-D0[k]);
            e25_blas = sqrt(e25_blas/nloops2);
            e25_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * EPART(A4[k]).conjugate() * C4[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e25_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e25_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e25_eigen = sqrt(e25_eigen/nloops2);
            e25_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * EPART(A5[k]).conjugate() * C5[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e25_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e25_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e25_smalleigen = sqrt(e25_smalleigen/nloops2);
            e25_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * B* * diag(C*)
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(B0[k](j,i) * C0[k](j));
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * MPART(B1[k].adjoint()) * DiagMatrixViewOf(C1[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e26_reg = 0.;
            for (int k=0; k<nloops2; ++k) e26_reg += NormSq(D1[k]-D0[k]);
            e26_reg = sqrt(e26_reg/nloops2);
            e26_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * MPART(B2[k].adjoint()) * DiagMatrixViewOf(C2[k].conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e26_small = 0.;
            for (int k=0; k<nloops2; ++k) e26_small += NormSq(D2[k]-D0[k]);
            e26_small = sqrt(e26_small/nloops2);
            e26_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int j=0;j<N;++j) {
                tmv::Vector<T> temp(COL_LEN(j));
                T x = T(8,9) * std::conj(C3[k](j));
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(B3[k].row(j,COL_START(j),COL_END(j)).cptr()),N,
                    BP(temp.ptr()),1);
                BLASDNAME(scal) (COL_LEN(j),dmone,temp.realPart().ptr()+1,2);
                BLASNAME(axpy) (
                    COL_LEN(j),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e26_blas = 0.;
            for (int k=0; k<nloops2; ++k) e26_blas += NormSq(D3[k]-D0[k]);
            e26_blas = sqrt(e26_blas/nloops2);
            e26_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * EPART(B4[k].adjoint()) * C4[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e26_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e26_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e26_eigen = sqrt(e26_eigen/nloops2);
            e26_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * EPART(B5[k].adjoint()) * C5[k].conjugate().asDiagonal();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e26_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e26_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e26_smalleigen = sqrt(e26_smalleigen/nloops2);
            e26_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif
#endif

    }

    std::cout<<"D = A * diag(C)         "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas;
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"D = B * diag(C)         "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"D = -A * diag(C)        "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"D = -B * diag(C)        "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"D = 7 * A * diag(C)     "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
    std::cout<<"D = 7 * B * diag(C)     "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"D -= A * diag(C)        "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"D -= B * diag(C)        "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
    std::cout<<"D += 8 * A * diag(C)    "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
    std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
    std::cout<<"D += 8 * B * diag(C)    "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
    std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"D = (7,1)*A*diag(C)     "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
    std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*B*diag(C)     "<<t12_reg<<"  "<<t12_small<<"  "<<t12_blas;
    std::cout<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*A*diag(C)    "<<t13_reg<<"  "<<t13_small<<"  "<<t13_blas;
    std::cout<<"  "<<t13_eigen<<"  "<<t13_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*B*diag(C)    "<<t14_reg<<"  "<<t14_small<<"  "<<t14_blas;
    std::cout<<"  "<<t14_eigen<<"  "<<t14_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*A* *diag(C)   "<<t15_reg<<"  "<<t15_small<<"  "<<t15_blas;
    std::cout<<"  "<<t15_eigen<<"  "<<t15_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*B* *diag(C)   "<<t17_reg<<"  "<<t17_small<<"  "<<t17_blas;
    std::cout<<"  "<<t17_eigen<<"  "<<t17_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*A* *diag(C)  "<<t16_reg<<"  "<<t16_small<<"  "<<t16_blas;
    std::cout<<"  "<<t16_eigen<<"  "<<t16_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*B* *diag(C)  "<<t18_reg<<"  "<<t18_small<<"  "<<t18_blas;
    std::cout<<"  "<<t18_eigen<<"  "<<t18_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*A*diag(C*)    "<<t19_reg<<"  "<<t19_small<<"  "<<t19_blas;
    std::cout<<"  "<<t19_eigen<<"  "<<t19_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*B*diag(C*)    "<<t20_reg<<"  "<<t20_small<<"  "<<t20_blas;
    std::cout<<"  "<<t20_eigen<<"  "<<t20_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*A*diag(C*)   "<<t21_reg<<"  "<<t21_small<<"  "<<t21_blas;
    std::cout<<"  "<<t21_eigen<<"  "<<t21_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*B*diag(C*)   "<<t22_reg<<"  "<<t22_small<<"  "<<t22_blas;
    std::cout<<"  "<<t22_eigen<<"  "<<t22_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*A* *diag(C*)  "<<t23_reg<<"  "<<t23_small<<"  "<<t23_blas;
    std::cout<<"  "<<t23_eigen<<"  "<<t23_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*B* *diag(C*)  "<<t24_reg<<"  "<<t24_small<<"  "<<t24_blas;
    std::cout<<"  "<<t24_eigen<<"  "<<t24_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*A* *diag(C*) "<<t25_reg<<"  "<<t25_small<<"  "<<t25_blas;
    std::cout<<"  "<<t25_eigen<<"  "<<t25_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*B* *diag(C*) "<<t26_reg<<"  "<<t26_small<<"  "<<t26_blas;
    std::cout<<"  "<<t26_eigen<<"  "<<t26_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"D = A * diag(C)         "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"D = B * diag(C)         "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"D = -A * diag(C)        "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"D = -B * diag(C)        "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"D = 7 * A * diag(C)     "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
    std::cout<<"D = 7 * B * diag(C)     "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
    std::cout<<"D -= A * diag(C)        "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
    std::cout<<"D -= B * diag(C)        "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
    std::cout<<"D += 8 * A * diag(C)    "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
    std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
    std::cout<<"D += 8 * B * diag(C)    "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
    std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"D = (7,1) * A * diag(C) "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
    std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * B * diag(C) "<<e12_reg<<"  "<<e12_small<<"  "<<e12_blas;
    std::cout<<"  "<<e12_eigen<<"  "<<e12_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*A*diag(C)    "<<e13_reg<<"  "<<e13_small<<"  "<<e13_blas;
    std::cout<<"  "<<e13_eigen<<"  "<<e13_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*B*diag(C)    "<<e14_reg<<"  "<<e14_small<<"  "<<e14_blas;
    std::cout<<"  "<<e14_eigen<<"  "<<e14_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*A* *diag(C)   "<<e15_reg<<"  "<<e15_small<<"  "<<e15_blas;
    std::cout<<"  "<<e15_eigen<<"  "<<e15_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*B* *diag(C)   "<<e17_reg<<"  "<<e17_small<<"  "<<e17_blas;
    std::cout<<"  "<<e17_eigen<<"  "<<e17_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*A* *diag(C)  "<<e16_reg<<"  "<<e16_small<<"  "<<e16_blas;
    std::cout<<"  "<<e16_eigen<<"  "<<e16_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*B* *diag(C)  "<<e18_reg<<"  "<<e18_small<<"  "<<e18_blas;
    std::cout<<"  "<<e18_eigen<<"  "<<e18_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*A*diag(C*)    "<<e19_reg<<"  "<<e19_small<<"  "<<e19_blas;
    std::cout<<"  "<<e19_eigen<<"  "<<e19_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*B*diag(C*)    "<<e20_reg<<"  "<<e20_small<<"  "<<e20_blas;
    std::cout<<"  "<<e20_eigen<<"  "<<e20_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*A*diag(C*)   "<<e21_reg<<"  "<<e21_small<<"  "<<e21_blas;
    std::cout<<"  "<<e21_eigen<<"  "<<e21_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*B*diag(C*)   "<<e22_reg<<"  "<<e22_small<<"  "<<e22_blas;
    std::cout<<"  "<<e22_eigen<<"  "<<e22_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*A* *diag(C*)  "<<e23_reg<<"  "<<e23_small<<"  "<<e23_blas;
    std::cout<<"  "<<e23_eigen<<"  "<<e23_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*B* *diag(C*)  "<<e24_reg<<"  "<<e24_small<<"  "<<e24_blas;
    std::cout<<"  "<<e24_eigen<<"  "<<e24_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*A* *diag(C*) "<<e25_reg<<"  "<<e25_small<<"  "<<e25_blas;
    std::cout<<"  "<<e25_eigen<<"  "<<e25_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*B* *diag(C*) "<<e26_reg<<"  "<<e26_small<<"  "<<e26_blas;
    std::cout<<"  "<<e26_eigen<<"  "<<e26_smalleigen<<std::endl;
#endif
    std::cout<<"\n\n";
#endif
}
#endif

#ifdef DOMULTDM
static void MultDM(
    const std::vector<tmv::Matrix<T> >& A1,
    const std::vector<tmv::Matrix<T> >& B1,
    const std::vector<tmv::Vector<T> >& C1)
{
    std::vector<tmv::Matrix<T> > D1 = A1;

#ifdef ERRORCHECK
    std::vector<tmv::Matrix<T> > A0 = A1;
    std::vector<tmv::Matrix<T> > B0 = B1;
    std::vector<tmv::Vector<T> > C0 = C1;
    std::vector<tmv::Matrix<T> > D0 = D1;
#endif

#ifdef DOSMALL
    std::vector<tmv::SmallMatrix<T,M,N> > A2(nloops2);
    std::vector<tmv::SmallMatrix<T,N,M> > B2(nloops2);
    std::vector<tmv::SmallVector<T,M> > C2(nloops2);
    std::vector<tmv::SmallMatrix<T,M,N> > D2(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2[k] = A1[k]; B2[k] = B1[k]; C2[k] = C1[k]; D2[k] = D1[k];
    }
#endif

#ifdef DOBLAS
    std::vector<tmv::Matrix<T> > A3 = A1;
    std::vector<tmv::Matrix<T> > B3 = B1;
    std::vector<tmv::Vector<T> > C3 = C1;
    std::vector<tmv::Matrix<T> > D3 = D1;
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENM> A4(nloops2,Eigen::EIGENM(M,N));
    std::vector<Eigen::EIGENM> B4(nloops2,Eigen::EIGENM(N,M));
    std::vector<Eigen::EIGENV> C4(nloops2,Eigen::EIGENV(M));
    std::vector<Eigen::EIGENM> D4(nloops2,Eigen::EIGENM(M,N));
    for(int k=0;k<nloops2;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) A4[k](i,j) = A1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) B4[k](j,i) = B1[k](j,i);
        for(int i=0;i<M;++i) C4[k](i) = C1[k](i);
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) D4[k](i,j) = D1[k](i,j);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<T,M,N> > A5;
    std::vector<Eigen::Matrix<T,N,M> > B5;
    std::vector<Eigen::Matrix<T,M,1> > C5;
    std::vector<Eigen::Matrix<T,M,N> > D5;
    if (nloops2x) { 
        A5.resize(nloops2x); B5.resize(nloops2x);
        C5.resize(nloops2x); D5.resize(nloops2x);
    }
    for(int k=0;k<nloops2x;++k) {
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) A5[k](i,j) = A1[k](i,j);
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) B5[k](j,i) = B1[k](j,i);
        for(int i=0;i<M;++i) C5[k](i) = C1[k](i);
        for(int j=0;j<N;++j) for(int i=0;i<M;++i) D5[k](i,j) = D1[k](i,j);
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
    double t17_reg=0., t17_small=0., t17_blas=0., t17_eigen=0., t17_smalleigen=0.;
    double t18_reg=0., t18_small=0., t18_blas=0., t18_eigen=0., t18_smalleigen=0.;
    double t19_reg=0., t19_small=0., t19_blas=0., t19_eigen=0., t19_smalleigen=0.;
    double t20_reg=0., t20_small=0., t20_blas=0., t20_eigen=0., t20_smalleigen=0.;
    double t21_reg=0., t21_small=0., t21_blas=0., t21_eigen=0., t21_smalleigen=0.;
    double t22_reg=0., t22_small=0., t22_blas=0., t22_eigen=0., t22_smalleigen=0.;
    double t23_reg=0., t23_small=0., t23_blas=0., t23_eigen=0., t23_smalleigen=0.;
    double t24_reg=0., t24_small=0., t24_blas=0., t24_eigen=0., t24_smalleigen=0.;
    double t25_reg=0., t25_small=0., t25_blas=0., t25_eigen=0., t25_smalleigen=0.;
    double t26_reg=0., t26_small=0., t26_blas=0., t26_eigen=0., t26_smalleigen=0.;
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
    double e17_reg=0., e17_small=0., e17_blas=0., e17_eigen=0., e17_smalleigen=0.;
    double e18_reg=0., e18_small=0., e18_blas=0., e18_eigen=0., e18_smalleigen=0.;
    double e19_reg=0., e19_small=0., e19_blas=0., e19_eigen=0., e19_smalleigen=0.;
    double e20_reg=0., e20_small=0., e20_blas=0., e20_eigen=0., e20_smalleigen=0.;
    double e21_reg=0., e21_small=0., e21_blas=0., e21_eigen=0., e21_smalleigen=0.;
    double e22_reg=0., e22_small=0., e22_blas=0., e22_eigen=0., e22_smalleigen=0.;
    double e23_reg=0., e23_small=0., e23_blas=0., e23_eigen=0., e23_smalleigen=0.;
    double e24_reg=0., e24_small=0., e24_blas=0., e24_eigen=0., e24_smalleigen=0.;
    double e25_reg=0., e25_small=0., e25_blas=0., e25_eigen=0., e25_smalleigen=0.;
    double e26_reg=0., e26_small=0., e26_blas=0., e26_eigen=0., e26_smalleigen=0.;
#endif

    for (int n=0; n<nloops1; ++n) {

#if 1 // D = diag(C) * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = A0[k](i,j) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = DiagMatrixViewOf(C1[k]) * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_reg = 0.;
            for (int k=0; k<nloops2; ++k) e1_reg += NormSq(D1[k]-D0[k]);
            e1_reg = sqrt(e1_reg/nloops2);
            e1_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = DiagMatrixViewOf(C2[k]) * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_small = 0.;
            for (int k=0; k<nloops2; ++k) e1_small += NormSq(D2[k]-D0[k]);
            e1_small = sqrt(e1_small/nloops2);
            e1_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
            for(int i=0;i<M;++i) {
                T x = C3[k](i);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_blas = 0.;
            for (int k=0; k<nloops2; ++k) e1_blas += NormSq(D3[k]-D0[k]);
            e1_blas = sqrt(e1_blas/nloops2);
            e1_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = C4[k].asDiagonal() * EPART(A4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e1_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e1_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e1_eigen = sqrt(e1_eigen/nloops2);
            e1_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = C5[k].asDiagonal() * EPART(A5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e1_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e1_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e1_smalleigen = sqrt(e1_smalleigen/nloops2);
            e1_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = diag(C) * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = B0[k](j,i) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = DiagMatrixViewOf(C1[k]) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_reg = 0.;
            for (int k=0; k<nloops2; ++k) 
                e2_reg += NormSq(D1[k]-D0[k]);
            e2_reg = sqrt(e2_reg/nloops2);
            e2_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = DiagMatrixViewOf(C2[k]) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_small = 0.;
            for (int k=0; k<nloops2; ++k) e2_small += NormSq(D2[k]-D0[k]);
            e2_small = sqrt(e2_small/nloops2);
            e2_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                T x = C3[k](i);
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(B3[k].col(i,ROW_START(i),ROW_END(i)).cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_blas = 0.;
            for (int k=0; k<nloops2; ++k) e2_blas += NormSq(D3[k]-D0[k]);
            e2_blas = sqrt(e2_blas/nloops2);
            e2_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = C4[k].asDiagonal() * EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e2_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e2_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e2_eigen = sqrt(e2_eigen/nloops2);
            e2_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = C5[k].asDiagonal() * EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e2_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e2_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e2_smalleigen = sqrt(e2_smalleigen/nloops2);
            e2_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = -diag(C) * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = -A0[k](i,j) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = -C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = -DiagMatrixViewOf(C1[k]) * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_reg = 0.;
            for (int k=0; k<nloops2; ++k) e3_reg += NormSq(D1[k]-D0[k]);
            e3_reg = sqrt(e3_reg/nloops2);
            e3_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = -DiagMatrixViewOf(C2[k]) * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_small = 0.;
            for (int k=0; k<nloops2; ++k) e3_small += NormSq(D2[k]-D0[k]);
            e3_small = sqrt(e3_small/nloops2);
            e3_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
            for(int i=0;i<M;++i) {
                T x = -C3[k](i);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,mone,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_blas = 0.;
            for (int k=0; k<nloops2; ++k) e3_blas += NormSq(D3[k]-D0[k]);
            e3_blas = sqrt(e3_blas/nloops2);
            e3_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            // Note: Without the (), this is _really_ slow for largish matrices.
            EPART2(D4[k]) = -(C4[k].asDiagonal() * EPART(A4[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e3_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e3_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e3_eigen = sqrt(e3_eigen/nloops2);
            e3_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = -(C5[k].asDiagonal() * EPART(A5[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e3_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e3_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e3_smalleigen = sqrt(e3_smalleigen/nloops2);
            e3_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = -diag(C) * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = -B0[k](j,i) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = -C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = -DiagMatrixViewOf(C1[k]) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_reg = 0.;
            for (int k=0; k<nloops2; ++k) e4_reg += NormSq(D1[k]-D0[k]);
            e4_reg = sqrt(e4_reg/nloops2);
            e4_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = -DiagMatrixViewOf(C2[k]) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_small = 0.;
            for (int k=0; k<nloops2; ++k) e4_small += NormSq(D2[k]-D0[k]);
            e4_small = sqrt(e4_small/nloops2);
            e4_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                T x = -C3[k](i);
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(B3[k].col(i,ROW_START(i),ROW_END(i)).cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,mone,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_blas = 0.;
            for (int k=0; k<nloops2; ++k) e4_blas += NormSq(D3[k]-D0[k]);
            e4_blas = sqrt(e4_blas/nloops2);
            e4_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = -(C4[k].asDiagonal() * EPART(B4[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e4_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e4_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e4_eigen = sqrt(e4_eigen/nloops2);
            e4_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = -(C5[k].asDiagonal() * EPART(B5[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e4_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e4_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e4_smalleigen = sqrt(e4_smalleigen/nloops2);
            e4_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = 7 * diag(C) * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = RT(7) * A0[k](i,j) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = RT(7) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = 7 * DiagMatrixViewOf(C1[k]) * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_reg = 0.;
            for (int k=0; k<nloops2; ++k) e5_reg += NormSq(D1[k]-D0[k]);
            e5_reg = sqrt(e5_reg/nloops2);
            e5_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = 7 * DiagMatrixViewOf(C2[k]) * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_small = 0.;
            for (int k=0; k<nloops2; ++k) e5_small += NormSq(D2[k]-D0[k]);
            e5_small = sqrt(e5_small/nloops2);
            e5_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
            for(int i=0;i<M;++i) {
                T x = T(7)*C3[k](i);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,seven,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_blas = 0.;
            for (int k=0; k<nloops2; ++k) e5_blas += NormSq(D3[k]-D0[k]);
            e5_blas = sqrt(e5_blas/nloops2);
            e5_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = 7 * (C4[k].asDiagonal() * EPART(A4[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e5_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e5_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e5_eigen = sqrt(e5_eigen/nloops2);
            e5_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = 7 * (C5[k].asDiagonal() * EPART(A5[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e5_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e5_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e5_smalleigen = sqrt(e5_smalleigen/nloops2);
            e5_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = 7 * diag(C) * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = RT(7) * B0[k](j,i) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = RT(7) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = 7 * DiagMatrixViewOf(C1[k]) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_reg = 0.;
            for (int k=0; k<nloops2; ++k) e6_reg += NormSq(D1[k]-D0[k]);
            e6_reg = sqrt(e6_reg/nloops2);
            e6_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = 7 * DiagMatrixViewOf(C2[k]) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_small = 0.;
            for (int k=0; k<nloops2; ++k) e6_small += NormSq(D2[k]-D0[k]);
            e6_small = sqrt(e6_small/nloops2);
            e6_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                T x = T(7) * C3[k](i);
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(B3[k].col(i,ROW_START(i),ROW_END(i)).cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,seven,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_blas = 0.;
            for (int k=0; k<nloops2; ++k) e6_blas += NormSq(D3[k]-D0[k]);
            e6_blas = sqrt(e6_blas/nloops2);
            e6_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = 7 * (C4[k].asDiagonal() * EPART(B4[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e6_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e6_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e6_eigen = sqrt(e6_eigen/nloops2);
            e6_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = 7 * (C5[k].asDiagonal() * EPART(B5[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e6_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e6_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e6_smalleigen = sqrt(e6_smalleigen/nloops2);
            e6_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D -= diag(C) * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) -= A0[k](i,j) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) -= C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) -= DiagMatrixViewOf(C1[k]) * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_reg = 0.;
            for (int k=0; k<nloops2; ++k) e7_reg += NormSq(D1[k]-D0[k]);
            e7_reg = sqrt(e7_reg/nloops2);
            e7_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) -= DiagMatrixViewOf(C2[k]) * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_small = 0.;
            for (int k=0; k<nloops2; ++k) e7_small += NormSq(D2[k]-D0[k]);
            e7_small = sqrt(e7_small/nloops2);
            e7_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                tmv::Vector<T> temp(ROW_LEN(i));
                T x = -C3[k](i);
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(A3[k].row(i,ROW_START(i),ROW_END(i)).cptr()),M,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    ROW_LEN(i),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,mone,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_blas = 0.;
            for (int k=0; k<nloops2; ++k) e7_blas += NormSq(D3[k]-D0[k]);
            e7_blas = sqrt(e7_blas/nloops2);
            e7_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) -= C4[k].asDiagonal() * EPART(A4[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e7_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e7_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e7_eigen = sqrt(e7_eigen/nloops2);
            e7_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) -= C5[k].asDiagonal() * EPART(A5[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e7_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e7_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e7_smalleigen = sqrt(e7_smalleigen/nloops2);
            e7_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D -= diag(C) * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) -= B0[k](j,i) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) -= C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) -= DiagMatrixViewOf(C1[k]) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_reg = 0.;
            for (int k=0; k<nloops2; ++k) e8_reg += NormSq(D1[k]-D0[k]);
            e8_reg = sqrt(e8_reg/nloops2);
            e8_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) -= DiagMatrixViewOf(C2[k]) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_small = 0.;
            for (int k=0; k<nloops2; ++k) e8_small += NormSq(D2[k]-D0[k]);
            e8_small = sqrt(e8_small/nloops2);
            e8_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                tmv::Vector<T> temp(ROW_LEN(i));
                T x = -C3[k](i);
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(B3[k].col(i,ROW_START(i),ROW_END(i)).cptr()),1,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    ROW_LEN(i),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,mone,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_blas = 0.;
            for (int k=0; k<nloops2; ++k) e8_blas += NormSq(D3[k]-D0[k]);
            e8_blas = sqrt(e8_blas/nloops2);
            e8_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) -= C4[k].asDiagonal() * EPART(B4[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e8_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e8_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e8_eigen = sqrt(e8_eigen/nloops2);
            e8_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) -= C5[k].asDiagonal() * EPART(B5[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e8_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e8_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e8_smalleigen = sqrt(e8_smalleigen/nloops2);
            e8_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += 8 * diag(C) * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += RT(8) * A0[k](i,j) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) += RT(8) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += 8 * DiagMatrixViewOf(C1[k]) * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_reg = 0.;
            for (int k=0; k<nloops2; ++k) e9_reg += NormSq(D1[k]-D0[k]);
            e9_reg = sqrt(e9_reg/nloops2);
            e9_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += 8 * DiagMatrixViewOf(C2[k]) * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_small = 0.;
            for (int k=0; k<nloops2; ++k) e9_small += NormSq(D2[k]-D0[k]);
            e9_small = sqrt(e9_small/nloops2);
            e9_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                tmv::Vector<T> temp(ROW_LEN(i));
                T x = T(8) * C3[k](i);
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(A3[k].row(i,ROW_START(i),ROW_END(i)).cptr()),M,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    ROW_LEN(i),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,eight,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_blas = 0.;
            for (int k=0; k<nloops2; ++k) e9_blas += NormSq(D3[k]-D0[k]);
            e9_blas = sqrt(e9_blas/nloops2);
            e9_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += 8 * (C4[k].asDiagonal() * EPART(A4[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e9_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e9_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e9_eigen = sqrt(e9_eigen/nloops2);
            e9_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += 8 * (C5[k].asDiagonal() * EPART(A5[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e9_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e9_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e9_smalleigen = sqrt(e9_smalleigen/nloops2);
            e9_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += 8 * diag(C) * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += RT(8) * B0[k](j,i) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) += RT(8) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += 8 * DiagMatrixViewOf(C1[k]) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_reg = 0.;
            for (int k=0; k<nloops2; ++k) e10_reg += NormSq(D1[k]-D0[k]);
            e10_reg = sqrt(e10_reg/nloops2);
            e10_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += 8 * DiagMatrixViewOf(C2[k]) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_small = 0.;
            for (int k=0; k<nloops2; ++k) e10_small += NormSq(D2[k]-D0[k]);
            e10_small = sqrt(e10_small/nloops2);
            e10_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                tmv::Vector<T> temp(ROW_LEN(i));
                T x = T(8) * C3[k](i);
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(B3[k].col(i,ROW_START(i),ROW_END(i)).cptr()),1,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    ROW_LEN(i),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,eight,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_blas = 0.;
            for (int k=0; k<nloops2; ++k) e10_blas += NormSq(D3[k]-D0[k]);
            e10_blas = sqrt(e10_blas/nloops2);
            e10_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += 8 * (C4[k].asDiagonal() * EPART(B4[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e10_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e10_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e10_eigen = sqrt(e10_eigen/nloops2);
            e10_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += 8 * (C5[k].asDiagonal() * EPART(B5[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e10_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e10_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e10_smalleigen = sqrt(e10_smalleigen/nloops2);
            e10_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#ifdef TISCOMPLEX
#if 1 // D = (7,1) * diag(C) * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * A0[k](i,j) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * DiagMatrixViewOf(C1[k]) * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_reg = 0.;
            for (int k=0; k<nloops2; ++k) e11_reg += NormSq(D1[k]-D0[k]);
            e11_reg = sqrt(e11_reg/nloops2);
            e11_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * DiagMatrixViewOf(C2[k]) * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_small = 0.;
            for (int k=0; k<nloops2; ++k) e11_small += NormSq(D2[k]-D0[k]);
            e11_small = sqrt(e11_small/nloops2);
            e11_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
            for(int i=0;i<M;++i) {
                T x = T(7,1)*C3[k](i);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_blas = 0.;
            for (int k=0; k<nloops2; ++k) e11_blas += NormSq(D3[k]-D0[k]);
            e11_blas = sqrt(e11_blas/nloops2);
            e11_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * (C4[k].asDiagonal() * EPART(A4[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e11_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e11_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e11_eigen = sqrt(e11_eigen/nloops2);
            e11_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * (C5[k].asDiagonal() * EPART(A5[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e11_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e11_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e11_smalleigen = sqrt(e11_smalleigen/nloops2);
            e11_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * diag(C) * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * B0[k](j,i) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * DiagMatrixViewOf(C1[k]) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_reg = 0.;
            for (int k=0; k<nloops2; ++k) e12_reg += NormSq(D1[k]-D0[k]);
            e12_reg = sqrt(e12_reg/nloops2);
            e12_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * DiagMatrixViewOf(C2[k]) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_small = 0.;
            for (int k=0; k<nloops2; ++k) e12_small += NormSq(D2[k]-D0[k]);
            e12_small = sqrt(e12_small/nloops2);
            e12_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                T x = T(7,1) * C3[k](i);
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(B3[k].col(i,ROW_START(i),ROW_END(i)).cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_blas = 0.;
            for (int k=0; k<nloops2; ++k) e12_blas += NormSq(D3[k]-D0[k]);
            e12_blas = sqrt(e12_blas/nloops2);
            e12_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * (C4[k].asDiagonal() * EPART(B4[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e12_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e12_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e12_eigen = sqrt(e12_eigen/nloops2);
            e12_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * (C5[k].asDiagonal() * EPART(B5[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e12_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e12_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e12_smalleigen = sqrt(e12_smalleigen/nloops2);
            e12_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * diag(C) * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * A0[k](i,j) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * DiagMatrixViewOf(C1[k]) * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_reg = 0.;
            for (int k=0; k<nloops2; ++k) e13_reg += NormSq(D1[k]-D0[k]);
            e13_reg = sqrt(e13_reg/nloops2);
            e13_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * DiagMatrixViewOf(C2[k]) * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_small = 0.;
            for (int k=0; k<nloops2; ++k) e13_small += NormSq(D2[k]-D0[k]);
            e13_small = sqrt(e13_small/nloops2);
            e13_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            tmv::Matrix<T> temp(M,N);
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(temp.ptr()),1);
            for(int i=0;i<M;++i) {
                T x = T(8,9) * C3[k](i);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(temp.row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#if (PART == 1)
            BLASNAME(axpy) (M*N,one,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(axpy) (
                    COL_LEN(j),one,
                    BP(temp.col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_blas = 0.;
            for (int k=0; k<nloops2; ++k) e13_blas += NormSq(D3[k]-D0[k]);
            e13_blas = sqrt(e13_blas/nloops2);
            e13_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * (C4[k].asDiagonal() * EPART(A4[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e13_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e13_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e13_eigen = sqrt(e13_eigen/nloops2);
            e13_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * (C5[k].asDiagonal() * EPART(A5[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e13_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e13_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e13_smalleigen = sqrt(e13_smalleigen/nloops2);
            e13_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * diag(C) * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * B0[k](j,i) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * DiagMatrixViewOf(C1[k]) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_reg = 0.;
            for (int k=0; k<nloops2; ++k) e14_reg += NormSq(D1[k]-D0[k]);
            e14_reg = sqrt(e14_reg/nloops2);
            e14_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * DiagMatrixViewOf(C2[k]) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_small = 0.;
            for (int k=0; k<nloops2; ++k) e14_small += NormSq(D2[k]-D0[k]);
            e14_small = sqrt(e14_small/nloops2);
            e14_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                tmv::Vector<T> temp(ROW_LEN(i));
                T x = T(8,9) * C3[k](i);
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(B3[k].col(i,ROW_START(i),ROW_END(i)).cptr()),1,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    ROW_LEN(i),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_blas = 0.;
            for (int k=0; k<nloops2; ++k) e14_blas += NormSq(D3[k]-D0[k]);
            e14_blas = sqrt(e14_blas/nloops2);
            e14_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * (C4[k].asDiagonal() * EPART(B4[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e14_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e14_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e14_eigen = sqrt(e14_eigen/nloops2);
            e14_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * (C5[k].asDiagonal() * EPART(B5[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e14_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e14_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e14_smalleigen = sqrt(e14_smalleigen/nloops2);
            e14_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * diag(C) * A*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(A0[k](i,j)) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * DiagMatrixViewOf(C1[k]) * MPART(A1[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_reg = 0.;
            for (int k=0; k<nloops2; ++k) e15_reg += NormSq(D1[k]-D0[k]);
            e15_reg = sqrt(e15_reg/nloops2);
            e15_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * DiagMatrixViewOf(C2[k]) * MPART(A2[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_small = 0.;
            for (int k=0; k<nloops2; ++k) e15_small += NormSq(D2[k]-D0[k]);
            e15_small = sqrt(e15_small/nloops2);
            e15_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
            BLASDNAME(scal) (M*N,dmone,D3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
                BLASDNAME(scal) (
                    COL_LEN(j),dmone,
                    D3[k].col(j,COL_START(j),COL_END(j)).realPart().ptr()+1,2);
            }
#endif
            for(int i=0;i<M;++i) {
                T x = T(7,1)*C3[k](i);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_blas = 0.;
            for (int k=0; k<nloops2; ++k) e15_blas += NormSq(D3[k]-D0[k]);
            e15_blas = sqrt(e15_blas/nloops2);
            e15_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * (C4[k].asDiagonal() * EPART(A4[k]).conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e15_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e15_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e15_eigen = sqrt(e15_eigen/nloops2);
            e15_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * (C5[k].asDiagonal() * EPART(A5[k]).conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e15_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e15_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e15_smalleigen = sqrt(e15_smalleigen/nloops2);
            e15_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * diag(C) * B*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(B0[k](j,i)) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * DiagMatrixViewOf(C1[k]) * MPART(B1[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_reg = 0.;
            for (int k=0; k<nloops2; ++k) e17_reg += NormSq(D1[k]-D0[k]);
            e17_reg = sqrt(e17_reg/nloops2);
            e17_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * DiagMatrixViewOf(C2[k]) * MPART(B2[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_small = 0.;
            for (int k=0; k<nloops2; ++k) e17_small += NormSq(D2[k]-D0[k]);
            e17_small = sqrt(e17_small/nloops2);
            e17_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                T x = T(7,1) * C3[k](i);
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(B3[k].col(i,ROW_START(i),ROW_END(i)).cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
                BLASDNAME(scal) (
                    ROW_LEN(i),dmone,
                    D3[k].row(i,ROW_START(i),ROW_END(i)).realPart().ptr()+1,2*M);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_blas = 0.;
            for (int k=0; k<nloops2; ++k) e17_blas += NormSq(D3[k]-D0[k]);
            e17_blas = sqrt(e17_blas/nloops2);
            e17_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * (C4[k].asDiagonal() * EPART(B4[k].adjoint()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e17_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e17_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e17_eigen = sqrt(e17_eigen/nloops2);
            e17_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * (C5[k].asDiagonal() * EPART(B5[k].adjoint()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t17_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e17_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e17_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e17_smalleigen = sqrt(e17_smalleigen/nloops2);
            e17_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * diag(C) * A*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(A0[k](i,j)) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * DiagMatrixViewOf(C1[k]) * MPART(A1[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_reg = 0.;
            for (int k=0; k<nloops2; ++k) e16_reg += NormSq(D1[k]-D0[k]);
            e16_reg = sqrt(e16_reg/nloops2);
            e16_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * DiagMatrixViewOf(C2[k]) * MPART(A2[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_small = 0.;
            for (int k=0; k<nloops2; ++k) e16_small += NormSq(D2[k]-D0[k]);
            e16_small = sqrt(e16_small/nloops2);
            e16_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            tmv::Matrix<T> temp(M,N);
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(temp.ptr()),1);
            BLASDNAME(scal) (M*N,dmone,temp.realPart().ptr()+1,2);
            for(int i=0;i<M;++i) {
                T x = T(8,9) * C3[k](i);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(temp.row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#if (PART == 1)
            BLASNAME(axpy) (M*N,one,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(axpy) (
                    COL_LEN(j),one,
                    BP(temp.col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_blas = 0.;
            for (int k=0; k<nloops2; ++k) e16_blas += NormSq(D3[k]-D0[k]);
            e16_blas = sqrt(e16_blas/nloops2);
            e16_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * (C4[k].asDiagonal() * EPART(A4[k]).conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e16_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e16_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e16_eigen = sqrt(e16_eigen/nloops2);
            e16_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * (C5[k].asDiagonal() * EPART(A5[k]).conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e16_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e16_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e16_smalleigen = sqrt(e16_smalleigen/nloops2);
            e16_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * diag(C) * B*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(B0[k](j,i)) * C0[k](i);
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * C0[k](j);
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * DiagMatrixViewOf(C1[k]) * MPART(B1[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_reg = 0.;
            for (int k=0; k<nloops2; ++k) e18_reg += NormSq(D1[k]-D0[k]);
            e18_reg = sqrt(e18_reg/nloops2);
            e18_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * DiagMatrixViewOf(C2[k]) * MPART(B2[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_small = 0.;
            for (int k=0; k<nloops2; ++k) e18_small += NormSq(D2[k]-D0[k]);
            e18_small = sqrt(e18_small/nloops2);
            e18_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                tmv::Vector<T> temp(ROW_LEN(i));
                T x = T(8,9) * C3[k](i);
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(B3[k].col(i,ROW_START(i),ROW_END(i)).cptr()),1,
                    BP(temp.ptr()),1);
                BLASDNAME(scal) (ROW_LEN(i),dmone,temp.realPart().ptr()+1,2);
                BLASNAME(axpy) (
                    ROW_LEN(i),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_blas = 0.;
            for (int k=0; k<nloops2; ++k) e18_blas += NormSq(D3[k]-D0[k]);
            e18_blas = sqrt(e18_blas/nloops2);
            e18_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * (C4[k].asDiagonal() * EPART(B4[k].adjoint()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e18_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e18_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e18_eigen = sqrt(e18_eigen/nloops2);
            e18_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * (C5[k].asDiagonal() * EPART(B5[k].adjoint()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t18_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e18_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e18_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e18_smalleigen = sqrt(e18_smalleigen/nloops2);
            e18_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * diag(C*) * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * A0[k](i,j) * std::conj(C0[k](i));
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * DiagMatrixViewOf(C1[k].conjugate()) * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_reg = 0.;
            for (int k=0; k<nloops2; ++k) e19_reg += NormSq(D1[k]-D0[k]);
            e19_reg = sqrt(e19_reg/nloops2);
            e19_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * DiagMatrixViewOf(C2[k].conjugate()) * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_small = 0.;
            for (int k=0; k<nloops2; ++k) e19_small += NormSq(D2[k]-D0[k]);
            e19_small = sqrt(e19_small/nloops2);
            e19_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
            for(int i=0;i<M;++i) {
                T x = T(7,1)*std::conj(C3[k](i));
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2*(N+1));
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_blas = 0.;
            for (int k=0; k<nloops2; ++k) e19_blas += NormSq(D3[k]-D0[k]);
            e19_blas = sqrt(e19_blas/nloops2);
            e19_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * (C4[k].conjugate().asDiagonal() * EPART(A4[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e19_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e19_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e19_eigen = sqrt(e19_eigen/nloops2);
            e19_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * (C5[k].conjugate().asDiagonal() * EPART(A5[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t19_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e19_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e19_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e19_smalleigen = sqrt(e19_smalleigen/nloops2);
            e19_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * diag(C*) * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * B0[k](j,i) * std::conj(C0[k](i));
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * DiagMatrixViewOf(C1[k].conjugate()) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_reg = 0.;
            for (int k=0; k<nloops2; ++k) e20_reg += NormSq(D1[k]-D0[k]);
            e20_reg = sqrt(e20_reg/nloops2);
            e20_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * DiagMatrixViewOf(C2[k].conjugate()) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_small = 0.;
            for (int k=0; k<nloops2; ++k) e20_small += NormSq(D2[k]-D0[k]);
            e20_small = sqrt(e20_small/nloops2);
            e20_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                T x = T(7,1) * std::conj(C3[k](i));
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(B3[k].col(i,ROW_START(i),ROW_END(i)).cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2*(N+1));
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_blas = 0.;
            for (int k=0; k<nloops2; ++k) e20_blas += NormSq(D3[k]-D0[k]);
            e20_blas = sqrt(e20_blas/nloops2);
            e20_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * (C4[k].conjugate().asDiagonal() * EPART(B4[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e20_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e20_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e20_eigen = sqrt(e20_eigen/nloops2);
            e20_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * (C5[k].conjugate().asDiagonal() * EPART(B5[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t20_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e20_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e20_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e20_smalleigen = sqrt(e20_smalleigen/nloops2);
            e20_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * diag(C*) * A
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * A0[k](i,j) * std::conj(C0[k](i));
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * DiagMatrixViewOf(C1[k].conjugate()) * MPART(A1[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_reg = 0.;
            for (int k=0; k<nloops2; ++k) e21_reg += NormSq(D1[k]-D0[k]);
            e21_reg = sqrt(e21_reg/nloops2);
            e21_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * DiagMatrixViewOf(C2[k].conjugate()) * MPART(A2[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_small = 0.;
            for (int k=0; k<nloops2; ++k) e21_small += NormSq(D2[k]-D0[k]);
            e21_small = sqrt(e21_small/nloops2);
            e21_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            tmv::Matrix<T> temp(M,N);
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(temp.ptr()),1);
            for(int i=0;i<M;++i) {
                T x = T(8,9) * std::conj(C3[k](i));
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(temp.row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#if (PART == 1)
            BLASNAME(axpy) (M*N,one,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(axpy) (
                    COL_LEN(j),one,
                    BP(temp.col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_blas = 0.;
            for (int k=0; k<nloops2; ++k) e21_blas += NormSq(D3[k]-D0[k]);
            e21_blas = sqrt(e21_blas/nloops2);
            e21_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * (C4[k].conjugate().asDiagonal() * EPART(A4[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e21_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e21_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e21_eigen = sqrt(e21_eigen/nloops2);
            e21_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * (C5[k].conjugate().asDiagonal() * EPART(A5[k]));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t21_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e21_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e21_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e21_smalleigen = sqrt(e21_smalleigen/nloops2);
            e21_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * diag(C*) * B
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * B0[k](j,i) * std::conj(C0[k](i));
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * DiagMatrixViewOf(C1[k].conjugate()) * MPART(B1[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_reg = 0.;
            for (int k=0; k<nloops2; ++k) e22_reg += NormSq(D1[k]-D0[k]);
            e22_reg = sqrt(e22_reg/nloops2);
            e22_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * DiagMatrixViewOf(C2[k].conjugate()) * MPART(B2[k].transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_small = 0.;
            for (int k=0; k<nloops2; ++k) e22_small += NormSq(D2[k]-D0[k]);
            e22_small = sqrt(e22_small/nloops2);
            e22_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                tmv::Vector<T> temp(ROW_LEN(i));
                T x = T(8,9) * std::conj(C3[k](i));
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(B3[k].col(i,ROW_START(i),ROW_END(i)).cptr()),1,
                    BP(temp.ptr()),1);
                BLASNAME(axpy) (
                    ROW_LEN(i),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_blas = 0.;
            for (int k=0; k<nloops2; ++k) e22_blas += NormSq(D3[k]-D0[k]);
            e22_blas = sqrt(e22_blas/nloops2);
            e22_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * (C4[k].conjugate().asDiagonal() * EPART(B4[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e22_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e22_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e22_eigen = sqrt(e22_eigen/nloops2);
            e22_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * (C5[k].conjugate().asDiagonal() * EPART(B5[k].transpose()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t22_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e22_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e22_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e22_smalleigen = sqrt(e22_smalleigen/nloops2);
            e22_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * diag(C*) * A*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(A0[k](i,j) * C0[k](i));
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * DiagMatrixViewOf(C1[k].conjugate()) * MPART(A1[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e23_reg = 0.;
            for (int k=0; k<nloops2; ++k) e23_reg += NormSq(D1[k]-D0[k]);
            e23_reg = sqrt(e23_reg/nloops2);
            e23_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * DiagMatrixViewOf(C2[k].conjugate()) * MPART(A2[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e23_small = 0.;
            for (int k=0; k<nloops2; ++k) e23_small += NormSq(D2[k]-D0[k]);
            e23_small = sqrt(e23_small/nloops2);
            e23_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if (PART == 1)
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(D3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASNAME(copy) (
                    COL_LEN(j),
                    BP(A3[k].col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
            }
#endif
            for(int i=0;i<M;++i) {
                T x = T(7,-1)*C3[k](i);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#if (PART == 1)
            BLASDNAME(scal) (M*N,dmone,D3[k].realPart().ptr()+1,2);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) {
                BLASDNAME(scal) (
                    COL_LEN(j),dmone,
                    D3[k].col(j,COL_START(j),COL_END(j)).realPart().ptr()+1,2);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2*(N+1));
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e23_blas = 0.;
            for (int k=0; k<nloops2; ++k) e23_blas += NormSq(D3[k]-D0[k]);
            e23_blas = sqrt(e23_blas/nloops2);
            e23_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * (C4[k].conjugate().asDiagonal() * EPART(A4[k]).conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e23_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e23_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e23_eigen = sqrt(e23_eigen/nloops2);
            e23_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * (C5[k].conjugate().asDiagonal() * EPART(A5[k]).conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t23_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e23_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e23_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e23_smalleigen = sqrt(e23_smalleigen/nloops2);
            e23_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D = (7,1) * diag(C*) * B*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(B0[k](j,i) * C0[k](i));
                        else if (UNITPART(i,j))
                            D0[k](i,j) = T(7,1) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) = T(7,1) * DiagMatrixViewOf(C1[k].conjugate()) * MPART(B1[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e24_reg = 0.;
            for (int k=0; k<nloops2; ++k) e24_reg += NormSq(D1[k]-D0[k]);
            e24_reg = sqrt(e24_reg/nloops2);
            e24_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) = T(7,1) * DiagMatrixViewOf(C2[k].conjugate()) * MPART(B2[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e24_small = 0.;
            for (int k=0; k<nloops2; ++k) e24_small += NormSq(D2[k]-D0[k]);
            e24_small = sqrt(e24_small/nloops2);
            e24_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                T x = T(7,1) * std::conj(C3[k](i));
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(B3[k].col(i,ROW_START(i),ROW_END(i)).cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
                BLASDNAME(scal) (
                    ROW_LEN(i),dmone,
                    D3[k].row(i,ROW_START(i),ROW_END(i)).realPart().ptr()+1,2*M);
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASNAME(copy) (N,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,D3[k].realPart().ptr()+1,2*(N+1));
            BLASNAME(scal) (N,z71,BP(D3[k].ptr()),N+1);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e24_blas = 0.;
            for (int k=0; k<nloops2; ++k) e24_blas += NormSq(D3[k]-D0[k]);
            e24_blas = sqrt(e24_blas/nloops2);
            e24_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) = T(7,1) * (C4[k].conjugate().asDiagonal() * EPART(B4[k].adjoint()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e24_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e24_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e24_eigen = sqrt(e24_eigen/nloops2);
            e24_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) = T(7,1) * (C5[k].conjugate().asDiagonal() * EPART(B5[k].adjoint()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t24_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e24_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e24_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e24_smalleigen = sqrt(e24_smalleigen/nloops2);
            e24_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * diag(C*) * A*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(A0[k](i,j) * C0[k](i));
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * DiagMatrixViewOf(C1[k].conjugate()) * MPART(A1[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e25_reg = 0.;
            for (int k=0; k<nloops2; ++k) e25_reg += NormSq(D1[k]-D0[k]);
            e25_reg = sqrt(e25_reg/nloops2);
            e25_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * DiagMatrixViewOf(C2[k].conjugate()) * MPART(A2[k]).conjugate();

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e25_small = 0.;
            for (int k=0; k<nloops2; ++k) e25_small += NormSq(D2[k]-D0[k]);
            e25_small = sqrt(e25_small/nloops2);
            e25_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
            tmv::Matrix<T> temp(M,N);
            BLASNAME(copy) (M*N,BP(A3[k].cptr()),1,BP(temp.ptr()),1);
            BLASDNAME(scal) (M*N,dmone,temp.realPart().ptr()+1,2);
            for(int i=0;i<M;++i) {
                T x = T(8,9) * std::conj(C3[k](i));
                BLASNAME(scal) (
                    ROW_LEN(i),BX(x),
                    BP(temp.row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#if (PART == 1)
            BLASNAME(axpy) (M*N,one,BP(temp.cptr()),1,BP(D3[k].ptr()),1);
#elif ((PART >= 2) && (PART <= 5))
            for(int j=0;j<N;++j) 
                BLASNAME(axpy) (
                    COL_LEN(j),one,
                    BP(temp.col(j,COL_START(j),COL_END(j)).cptr()),1,
                    BP(D3[k].col(j,COL_START(j),COL_END(j)).ptr()),1);
#endif
#if ((PART == 3) || (PART == 5))
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e25_blas = 0.;
            for (int k=0; k<nloops2; ++k) e25_blas += NormSq(D3[k]-D0[k]);
            e25_blas = sqrt(e25_blas/nloops2);
            e25_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * (C4[k].conjugate().asDiagonal() * EPART(A4[k]).conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e25_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e25_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e25_eigen = sqrt(e25_eigen/nloops2);
            e25_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * (C5[k].conjugate().asDiagonal() * EPART(A5[k]).conjugate());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t25_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e25_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e25_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e25_smalleigen = sqrt(e25_smalleigen/nloops2);
            e25_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif

#if 1 // D += (8,9) * diag(C*) * B*
        ClearCache();

#ifdef ERRORCHECK
        if (n == 0) {
            for (int k=0; k<nloops2; ++k) {
                for(int i=0;i<M;++i) {
                    for(int j=0;j<N;++j) {
                        if (INPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(B0[k](j,i) * C0[k](i));
                        else if (UNITPART(i,j))
                            D0[k](i,j) += T(8,9) * std::conj(C0[k](j));
                    }
                }
            }
        }
#endif

#ifdef DOREG
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D1[k]) += T(8,9) * DiagMatrixViewOf(C1[k].conjugate()) * MPART(B1[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_reg += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e26_reg = 0.;
            for (int k=0; k<nloops2; ++k) e26_reg += NormSq(D1[k]-D0[k]);
            e26_reg = sqrt(e26_reg/nloops2);
            e26_reg /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            MPART2(D2[k]) += T(8,9) * DiagMatrixViewOf(C2[k].conjugate()) * MPART(B2[k].adjoint());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_small += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e26_small = 0.;
            for (int k=0; k<nloops2; ++k) e26_small += NormSq(D2[k]-D0[k]);
            e26_small = sqrt(e26_small/nloops2);
            e26_small /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOBLAS
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k) {
#if ((PART >= 1) && (PART <= 5))
            for(int i=0;i<M;++i) {
                tmv::Vector<T> temp(ROW_LEN(i));
                T x = T(8,9) * std::conj(C3[k](i));
                BLASNAME(copy) (
                    ROW_LEN(i),
                    BP(B3[k].col(i,ROW_START(i),ROW_END(i)).cptr()),1,
                    BP(temp.ptr()),1);
                BLASDNAME(scal) (ROW_LEN(i),dmone,temp.realPart().ptr()+1,2);
                BLASNAME(axpy) (
                    ROW_LEN(i),BX(x),
                    BP(temp.cptr()),1,
                    BP(D3[k].row(i,ROW_START(i),ROW_END(i)).ptr()),M);
            }
#endif
#if ((PART == 3) || (PART == 5))
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
            BLASNAME(axpy) (N,z89,BP(C3[k].cptr()),1,BP(D3[k].ptr()),N+1);
            BLASDNAME(scal) (N,dmone,C3[k].realPart().ptr()+1,2);
#endif
        }

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_blas += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e26_blas = 0.;
            for (int k=0; k<nloops2; ++k) e26_blas += NormSq(D3[k]-D0[k]);
            e26_blas = sqrt(e26_blas/nloops2);
            e26_blas /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EPART2(D4[k]) += T(8,9) * (C4[k].conjugate().asDiagonal() * EPART(B4[k].adjoint()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_eigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0) {
            e26_eigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e26_eigen += tmv::TMV_NORM(D4[k](i,j)-D0[k](i,j));
            e26_eigen = sqrt(e26_eigen/nloops2);
            e26_eigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EPART2(D5[k]) += T(8,9) * (C5[k].conjugate().asDiagonal() * EPART(B5[k].adjoint()));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t26_smalleigen += tb-ta;
#ifdef ERRORCHECK
        if (n == 0 && nloops2x) {
            e26_smalleigen = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<M;++i) for(int j=0;j<N;++j)
                    e26_smalleigen += tmv::TMV_NORM(D5[k](i,j)-D0[k](i,j));
            e26_smalleigen = sqrt(e26_smalleigen/nloops2);
            e26_smalleigen /= Norm(A0[0])*Norm(C0[0]);
        }
#endif
#endif
#endif
#endif

    }

    std::cout<<"D = diag(C) * A         "<<t1_reg<<"  "<<t1_small<<"  "<<t1_blas; 
    std::cout<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"D = diag(C) * B         "<<t2_reg<<"  "<<t2_small<<"  "<<t2_blas;
    std::cout<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"D = -diag(C) * A        "<<t3_reg<<"  "<<t3_small<<"  "<<t3_blas;
    std::cout<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"D = -diag(C) * B        "<<t4_reg<<"  "<<t4_small<<"  "<<t4_blas;
    std::cout<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"D = 7 * diag(C) * A     "<<t5_reg<<"  "<<t5_small<<"  "<<t5_blas;
    std::cout<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
    std::cout<<"D = 7 * diag(C) * B     "<<t6_reg<<"  "<<t6_small<<"  "<<t6_blas;
    std::cout<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"D -= diag(C) * A        "<<t7_reg<<"  "<<t7_small<<"  "<<t7_blas;
    std::cout<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"D -= diag(C) * B        "<<t8_reg<<"  "<<t8_small<<"  "<<t8_blas;
    std::cout<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
    std::cout<<"D += 8 * diag(C) * A    "<<t9_reg<<"  "<<t9_small<<"  "<<t9_blas;
    std::cout<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
    std::cout<<"D += 8 * diag(C) * B    "<<t10_reg<<"  "<<t10_small<<"  "<<t10_blas;
    std::cout<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"D = (7,1)*diag(C)*A     "<<t11_reg<<"  "<<t11_small<<"  "<<t11_blas;
    std::cout<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*diag(C)*B     "<<t12_reg<<"  "<<t12_small<<"  "<<t12_blas;
    std::cout<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C)*A    "<<t13_reg<<"  "<<t13_small<<"  "<<t13_blas;
    std::cout<<"  "<<t13_eigen<<"  "<<t13_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C)*B    "<<t14_reg<<"  "<<t14_small<<"  "<<t14_blas;
    std::cout<<"  "<<t14_eigen<<"  "<<t14_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*diag(C)*A*    "<<t15_reg<<"  "<<t15_small<<"  "<<t15_blas;
    std::cout<<"  "<<t15_eigen<<"  "<<t15_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*diag(C)*B*    "<<t17_reg<<"  "<<t17_small<<"  "<<t17_blas;
    std::cout<<"  "<<t17_eigen<<"  "<<t17_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C)*A*   "<<t16_reg<<"  "<<t16_small<<"  "<<t16_blas;
    std::cout<<"  "<<t16_eigen<<"  "<<t16_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C)*B*   "<<t18_reg<<"  "<<t18_small<<"  "<<t18_blas;
    std::cout<<"  "<<t18_eigen<<"  "<<t18_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*diag(C*)*A    "<<t19_reg<<"  "<<t19_small<<"  "<<t19_blas;
    std::cout<<"  "<<t19_eigen<<"  "<<t19_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*diag(C*)*B    "<<t20_reg<<"  "<<t20_small<<"  "<<t20_blas;
    std::cout<<"  "<<t20_eigen<<"  "<<t20_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C*)*A   "<<t21_reg<<"  "<<t21_small<<"  "<<t21_blas;
    std::cout<<"  "<<t21_eigen<<"  "<<t21_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C*)*B   "<<t22_reg<<"  "<<t22_small<<"  "<<t22_blas;
    std::cout<<"  "<<t22_eigen<<"  "<<t22_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*diag(C*)*A*   "<<t23_reg<<"  "<<t23_small<<"  "<<t23_blas;
    std::cout<<"  "<<t23_eigen<<"  "<<t23_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*diag(C*)*B*   "<<t24_reg<<"  "<<t24_small<<"  "<<t24_blas;
    std::cout<<"  "<<t24_eigen<<"  "<<t24_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C*)*A*  "<<t25_reg<<"  "<<t25_small<<"  "<<t25_blas;
    std::cout<<"  "<<t25_eigen<<"  "<<t25_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C*)*B*  "<<t26_reg<<"  "<<t26_small<<"  "<<t26_blas;
    std::cout<<"  "<<t26_eigen<<"  "<<t26_smalleigen<<std::endl;
#endif

#ifdef ERRORCHECK
    std::cout<<"errors:\n";
    std::cout<<"D = diag(C) * A         "<<e1_reg<<"  "<<e1_small<<"  "<<e1_blas;
    std::cout<<"  "<<e1_eigen<<"  "<<e1_smalleigen<<std::endl;
    std::cout<<"D = diag(C) * B         "<<e2_reg<<"  "<<e2_small<<"  "<<e2_blas;
    std::cout<<"  "<<e2_eigen<<"  "<<e2_smalleigen<<std::endl;
    std::cout<<"D = -diag(C) * A        "<<e3_reg<<"  "<<e3_small<<"  "<<e3_blas;
    std::cout<<"  "<<e3_eigen<<"  "<<e3_smalleigen<<std::endl;
    std::cout<<"D = -diag(C) * B        "<<e4_reg<<"  "<<e4_small<<"  "<<e4_blas;
    std::cout<<"  "<<e4_eigen<<"  "<<e4_smalleigen<<std::endl;
    std::cout<<"D = 7 * diag(C) * A     "<<e5_reg<<"  "<<e5_small<<"  "<<e5_blas;
    std::cout<<"  "<<e5_eigen<<"  "<<e5_smalleigen<<std::endl;
    std::cout<<"D = 7 * diag(C) * B     "<<e6_reg<<"  "<<e6_small<<"  "<<e6_blas;
    std::cout<<"  "<<e6_eigen<<"  "<<e6_smalleigen<<std::endl;
    std::cout<<"D -= diag(C) * A        "<<e7_reg<<"  "<<e7_small<<"  "<<e7_blas;
    std::cout<<"  "<<e7_eigen<<"  "<<e7_smalleigen<<std::endl;
    std::cout<<"D -= diag(C) * B        "<<e8_reg<<"  "<<e8_small<<"  "<<e8_blas;
    std::cout<<"  "<<e8_eigen<<"  "<<e8_smalleigen<<std::endl;
    std::cout<<"D += 8 * diag(C) * A    "<<e9_reg<<"  "<<e9_small<<"  "<<e9_blas;
    std::cout<<"  "<<e9_eigen<<"  "<<e9_smalleigen<<std::endl;
    std::cout<<"D += 8 * diag(C) * B    "<<e10_reg<<"  "<<e10_small<<"  "<<e10_blas;
    std::cout<<"  "<<e10_eigen<<"  "<<e10_smalleigen<<std::endl;
#ifdef TISCOMPLEX
    std::cout<<"D = (7,1) * diag(C) * A "<<e11_reg<<"  "<<e11_small<<"  "<<e11_blas;
    std::cout<<"  "<<e11_eigen<<"  "<<e11_smalleigen<<std::endl;
    std::cout<<"D = (7,1) * diag(C) * B "<<e12_reg<<"  "<<e12_small<<"  "<<e12_blas;
    std::cout<<"  "<<e12_eigen<<"  "<<e12_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C)*A    "<<e13_reg<<"  "<<e13_small<<"  "<<e13_blas;
    std::cout<<"  "<<e13_eigen<<"  "<<e13_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C)*B    "<<e14_reg<<"  "<<e14_small<<"  "<<e14_blas;
    std::cout<<"  "<<e14_eigen<<"  "<<e14_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*diag(C)*A*    "<<e15_reg<<"  "<<e15_small<<"  "<<e15_blas;
    std::cout<<"  "<<e15_eigen<<"  "<<e15_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*diag(C)*B*    "<<e17_reg<<"  "<<e17_small<<"  "<<e17_blas;
    std::cout<<"  "<<e17_eigen<<"  "<<e17_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C)*A*   "<<e16_reg<<"  "<<e16_small<<"  "<<e16_blas;
    std::cout<<"  "<<e16_eigen<<"  "<<e16_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C)*B*   "<<e18_reg<<"  "<<e18_small<<"  "<<e18_blas;
    std::cout<<"  "<<e18_eigen<<"  "<<e18_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*diag(C*)*A    "<<e19_reg<<"  "<<e19_small<<"  "<<e19_blas;
    std::cout<<"  "<<e19_eigen<<"  "<<e19_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*diag(C*)*B    "<<e20_reg<<"  "<<e20_small<<"  "<<e20_blas;
    std::cout<<"  "<<e20_eigen<<"  "<<e20_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C*)*A   "<<e21_reg<<"  "<<e21_small<<"  "<<e21_blas;
    std::cout<<"  "<<e21_eigen<<"  "<<e21_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C*)*B   "<<e22_reg<<"  "<<e22_small<<"  "<<e22_blas;
    std::cout<<"  "<<e22_eigen<<"  "<<e22_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*diag(C*)*A*   "<<e23_reg<<"  "<<e23_small<<"  "<<e23_blas;
    std::cout<<"  "<<e23_eigen<<"  "<<e23_smalleigen<<std::endl;
    std::cout<<"D = (7,1)*diag(C*)*B*   "<<e24_reg<<"  "<<e24_small<<"  "<<e24_blas;
    std::cout<<"  "<<e24_eigen<<"  "<<e24_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C*)*A*  "<<e25_reg<<"  "<<e25_small<<"  "<<e25_blas;
    std::cout<<"  "<<e25_eigen<<"  "<<e25_smalleigen<<std::endl;
    std::cout<<"D += (8,9)*diag(C*)*B*  "<<e26_reg<<"  "<<e26_small<<"  "<<e26_blas;
    std::cout<<"  "<<e26_eigen<<"  "<<e26_smalleigen<<std::endl;
#endif
    std::cout<<"\n\n";
#endif
}
#endif


#ifdef SIMPLE_VALUES
#ifdef TISCOMPLEX
#define RAND ( T(1+i,1+j) )
#else
#define RAND ( T(1+i+j) )
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
    srand(417191879);

#ifdef MEMDEBUG
    atexit(&DumpUnfreed);
#endif

    std::vector<tmv::Matrix<T> > A(nloops2,tmv::Matrix<T>(M,N));
    std::vector<tmv::Matrix<T> > B(nloops2,tmv::Matrix<T>(N,M));
    std::vector<tmv::Vector<T> > C(nloops2,tmv::Vector<T>(N));
    std::vector<tmv::Vector<T> > D(nloops2,tmv::Vector<T>(M));
    for(int k=0;k<nloops2;++k) {
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
            A[k](i,j) = RAND;
            B[k](j,i) = RAND;
        }
        for(int j=0;j<N;++j) {
            int i=0;
            C[k](j) = RAND;
        }
        for(int i=0;i<M;++i) {
            int j=0;
            D[k](i) = RAND;
        }
#ifdef ZERO_VALUES
        C[k].subVector(0,N/3).setZero();
        C[k].subVector(2*N/3,N).setZero();
#endif
    }
    std::cout<<"M,N = "<<M<<" , "<<N<<std::endl;
    std::cout<<"nloops = "<<nloops1<<" x "<<nloops2;
    std::cout<<" = "<<nloops1*nloops2<<std::endl;
#ifdef DOEIGENSMALL
    if (nloops2x == 0)
        std::cout<<"Matrix is too big for the stack, so no \"Eigen Known\" tests.\n";
#endif
    std::cout<<"\nTime for:               ";
    std::cout<<"  TMV    "<<" TMV Small"<<"   BLAS   "<<"  Eigen  "<<"Eigen Known"<<std::endl;

#ifdef DOMULTXM
    MultXM(A,B);
#endif

#ifdef DOADDMM
    AddMM(A,B);
#endif

#ifdef DONORM
    NormM(A);
#endif

#ifdef DOMULTMV
    MultMV(A,B,C,D);
#endif

#ifdef DORANK1
    Rank1Update(A,B,C,D);
#endif

#ifdef DOMULTMD
    MultMD(A,B,C);
#endif

#ifdef DOMULTDM
    MultDM(A,B,D);
#endif

#ifdef DOSWAP
    SwapM(A,B);
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

