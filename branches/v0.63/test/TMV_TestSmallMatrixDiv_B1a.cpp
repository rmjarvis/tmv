
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::StorageType stor, int N> 
static void TestSmallNonSquareDiv_B1a()
{
    tmv::SmallMatrix<T,6,N,stor> m;

    for(int i=0;i<6;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
    m.diag() *= T(11);
    m /= T(7);
    if (N > 1) m(1,0) = -2;
    if (N > 2) m(2,0) = 7;
    if (N > 3) m(3,0) = -10;
    if (N > 2) m(2,2) = 30;

    tmv::SmallMatrix<T,6,N,stor> a1 = m;
    tmv::SmallMatrix<std::complex<T>,6,N,stor> c1 = a1 * std::complex<T>(1,2);

    tmv::SmallMatrix<T,N,6,stor> a0;
    tmv::SmallMatrix<std::complex<T>,N,6,stor> c0;

    TestMatrixDivArith3a<T>(tmv::QR,a0,c0,a1,c1,"NonSquare"); 
#ifdef XTEST
    tmv::SmallMatrix<T,6,N,stor,tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,6,N,stor,tmv::FortranStyle> c1f = c1;

    TestMatrixDivArith3a<T>(tmv::QR,a0,c0,a1f,c1f,"NonSquare"); 
#endif
}

template <class T> 
void TestSmallMatrixDiv_B1a()
{
    TestSmallNonSquareDiv_B1a<T,tmv::ColMajor,2>();
    TestSmallNonSquareDiv_B1a<T,tmv::ColMajor,5>();
#ifdef XTEST
    TestSmallNonSquareDiv_B1a<T,tmv::ColMajor,1>();
    TestSmallNonSquareDiv_B1a<T,tmv::ColMajor,3>();
    TestSmallNonSquareDiv_B1a<T,tmv::ColMajor,4>();
    TestSmallNonSquareDiv_B1a<T,tmv::RowMajor,1>();
    TestSmallNonSquareDiv_B1a<T,tmv::RowMajor,2>();
    TestSmallNonSquareDiv_B1a<T,tmv::RowMajor,3>();
    TestSmallNonSquareDiv_B1a<T,tmv::RowMajor,4>();
    TestSmallNonSquareDiv_B1a<T,tmv::RowMajor,5>();
#endif
}

#ifdef INST_DOUBLE
template void TestSmallMatrixDiv_B1a<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrixDiv_B1a<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrixDiv_B1a<long double>();
#endif
