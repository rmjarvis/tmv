
#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::StorageType stor, int N> 
static void TestSmallSquareDiv_A1b()
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

    tmv::SmallVector<T,N> b = m.row(0);
    tmv::SmallVector<std::complex<T>,N> e = c.row(0);
    tmv::SmallVector<T,N> x = m.col(0);;
    tmv::SmallVector<std::complex<T>,N> y = c.col(0);

    TestMatrixDivArith3d(tmv::LU,a1,b,x,c1,e,y,"V/Square"); 
#if (XTEST & 32)
    tmv::SmallMatrix<T,N,N,stor|tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,stor|tmv::FortranStyle> c1f = c1;

    tmv::SmallVector<T,N,tmv::FortranStyle> bf = b;
    tmv::SmallVector<std::complex<T>,N,tmv::FortranStyle> ef = e;
    tmv::SmallVector<T,N,tmv::FortranStyle> xf = x;
    tmv::SmallVector<std::complex<T>,N,tmv::FortranStyle> yf = y;

    TestMatrixDivArith3d(tmv::LU,a1f,b,x,c1f,e,y,"V/Square"); 
    TestMatrixDivArith3d(tmv::LU,a1f,bf,x,c1f,ef,y,"V/Square"); 
    TestMatrixDivArith3d(tmv::LU,a1f,bf,xf,c1f,ef,yf,"V/Square"); 
#endif
}

template <class T> 
void TestSmallMatrixDiv_A1b()
{
    TestSmallSquareDiv_A1b<T,tmv::ColMajor,2>();
    TestSmallSquareDiv_A1b<T,tmv::ColMajor,5>();
#if (XTEST & 2)
    TestSmallSquareDiv_A1b<T,tmv::ColMajor,1>();
    TestSmallSquareDiv_A1b<T,tmv::ColMajor,3>();
    TestSmallSquareDiv_A1b<T,tmv::ColMajor,4>();
    TestSmallSquareDiv_A1b<T,tmv::RowMajor,1>();
    TestSmallSquareDiv_A1b<T,tmv::RowMajor,2>();
    TestSmallSquareDiv_A1b<T,tmv::RowMajor,3>();
    TestSmallSquareDiv_A1b<T,tmv::RowMajor,4>();
    TestSmallSquareDiv_A1b<T,tmv::RowMajor,5>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixDiv_A1b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixDiv_A1b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixDiv_A1b<long double>();
#endif
