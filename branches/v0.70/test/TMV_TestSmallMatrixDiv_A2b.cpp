
#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::StorageType stor, int N> 
static void TestSmallSquareDiv_A2b()
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
    tmv::SmallMatrix<T,N,N,tmv::RowMajor> a2b = m.transpose();

    tmv::SmallMatrix<std::complex<T>,N,N,stor> c1 = a1 * std::complex<T>(1,2);
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor> c2b = a2b * std::complex<T>(-3,4);
    c1.diag().addToAll(std::complex<T>(3,1));
    c2b.diag().addToAll(std::complex<T>(-5,8));
    c2b.row(0).addToAll(std::complex<T>(-2,-11));

    tmv::SmallMatrix<T,N,N,stor> a3;
    tmv::SmallMatrix<std::complex<T>,N,N,stor> c3;

    TestMatrixDivArith3b(tmv::LU,a1,a2b,a3,c1,c2b,c3,"Square/Square"); 
#if (XTEST & 32)
    tmv::SmallMatrix<T,N,N,stor|tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,stor|tmv::FortranStyle> c1f = c1;

    tmv::SmallMatrix<T,N,N,tmv::RowMajor|tmv::FortranStyle> a2fb = a2b;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor|tmv::FortranStyle> c2fb = c2b;
    tmv::SmallMatrix<T,N,N,stor|tmv::FortranStyle> a3f = a3;
    tmv::SmallMatrix<std::complex<T>,N,N,stor|tmv::FortranStyle> c3f = c3;

    TestMatrixDivArith3b(tmv::LU,a1f,a2b,a3,c1f,c2b,c3,"Square/Square"); 
    TestMatrixDivArith3b(tmv::LU,a1f,a2fb,a3,c1f,c2fb,c3,"Square/Square"); 
    TestMatrixDivArith3b(tmv::LU,a1f,a2fb,a3f,c1f,c2fb,c3f,"Square/Square"); 
#endif
}

template <class T> 
void TestSmallMatrixDiv_A2b()
{
    TestSmallSquareDiv_A2b<T,tmv::ColMajor,2>();
    TestSmallSquareDiv_A2b<T,tmv::ColMajor,5>();
#if (XTEST & 2)
    TestSmallSquareDiv_A2b<T,tmv::ColMajor,1>();
    TestSmallSquareDiv_A2b<T,tmv::ColMajor,3>();
    TestSmallSquareDiv_A2b<T,tmv::ColMajor,4>();
    TestSmallSquareDiv_A2b<T,tmv::RowMajor,1>();
    TestSmallSquareDiv_A2b<T,tmv::RowMajor,2>();
    TestSmallSquareDiv_A2b<T,tmv::RowMajor,3>();
    TestSmallSquareDiv_A2b<T,tmv::RowMajor,4>();
    TestSmallSquareDiv_A2b<T,tmv::RowMajor,5>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixDiv_A2b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixDiv_A2b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixDiv_A2b<long double>();
#endif
