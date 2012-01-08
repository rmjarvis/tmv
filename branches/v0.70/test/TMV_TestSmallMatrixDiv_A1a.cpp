
#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::StorageType stor, int N> 
static void TestSmallSquareDiv_A1a()
{
    tmv::SmallMatrix<T,N,N,stor> m;

    for(int i=0;i<N;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
    m.diag() *= T(11);
    m /= T(7);
    if (N > 1) m(1,0) = -2;
    if (N > 2) m(2,0) = 7;
    if (N > 3) m(3,0) = -10;
    if (N > 2) m(2,2) = 30;

    tmv::SmallMatrix<std::complex<T>,N,N,stor> c = m;

    tmv::SmallMatrix<T,N,N,stor> a1 = m;
    tmv::SmallMatrix<std::complex<T>,N,N,stor> c1 = a1 * std::complex<T>(1,2);
    c1.diag().addToAll(std::complex<T>(3,1));

    TestMatrixDivArith3a(tmv::LU,a1,c1,"Square"); 
#if (XTEST & 2)
    tmv::SmallMatrix<T,N,N,stor|tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,stor|tmv::FortranStyle> c1f = c1;

    TestMatrixDivArith3a(tmv::LU,a1f,c1f,"Square"); 
#endif
}

template <class T> 
void TestSmallMatrixDiv_A1a()
{
    TestSmallSquareDiv_A1a<T,tmv::ColMajor,2>();
    TestSmallSquareDiv_A1a<T,tmv::ColMajor,5>();
#if (XTEST & 2)
    TestSmallSquareDiv_A1a<T,tmv::ColMajor,1>();
    TestSmallSquareDiv_A1a<T,tmv::ColMajor,3>();
    TestSmallSquareDiv_A1a<T,tmv::ColMajor,4>();
    TestSmallSquareDiv_A1a<T,tmv::RowMajor,1>();
    TestSmallSquareDiv_A1a<T,tmv::RowMajor,2>();
    TestSmallSquareDiv_A1a<T,tmv::RowMajor,3>();
    TestSmallSquareDiv_A1a<T,tmv::RowMajor,4>();
    TestSmallSquareDiv_A1a<T,tmv::RowMajor,5>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixDiv_A1a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixDiv_A1a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixDiv_A1a<long double>();
#endif
