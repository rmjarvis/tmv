
#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::StorageType stor, int N> 
static void TestSmallNonSquareDiv_B1b()
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

    tmv::SmallVector<T,N> b = m.row(0);
    tmv::SmallVector<std::complex<T>,N> e = c1.row(0);
    tmv::SmallVector<T,6> x = m.col(0);
    tmv::SmallVector<std::complex<T>,6> y = c1.col(0);

    TestMatrixDivArith3d(tmv::QR,a1,x,b,c1,y,e,"V/NonSquare"); 
#if (XTEST & 32)
    tmv::SmallMatrix<T,6,N,stor|tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,6,N,stor|tmv::FortranStyle> c1f = c1;

    tmv::SmallVector<T,N,tmv::FortranStyle> bf = b;
    tmv::SmallVector<std::complex<T>,N,tmv::FortranStyle> ef = e;
    tmv::SmallVector<T,6,tmv::FortranStyle> xf = x;
    tmv::SmallVector<std::complex<T>,6,tmv::FortranStyle> yf = y;

    TestMatrixDivArith3d(tmv::QR,a1f,x,b,c1f,y,e,"V/NonSquare"); 
    TestMatrixDivArith3d(tmv::QR,a1f,xf,b,c1f,yf,e,"V/NonSquare"); 
    TestMatrixDivArith3d(tmv::QR,a1f,xf,bf,c1f,yf,ef,"V/NonSquare"); 
#endif
}

template <class T> 
void TestSmallMatrixDiv_B1b()
{
    TestSmallNonSquareDiv_B1b<T,tmv::ColMajor,2>();
    TestSmallNonSquareDiv_B1b<T,tmv::ColMajor,5>();
#if (XTEST & 2)
    TestSmallNonSquareDiv_B1b<T,tmv::ColMajor,1>();
    TestSmallNonSquareDiv_B1b<T,tmv::ColMajor,3>();
    TestSmallNonSquareDiv_B1b<T,tmv::ColMajor,4>();
    TestSmallNonSquareDiv_B1b<T,tmv::RowMajor,1>();
    TestSmallNonSquareDiv_B1b<T,tmv::RowMajor,2>();
    TestSmallNonSquareDiv_B1b<T,tmv::RowMajor,3>();
    TestSmallNonSquareDiv_B1b<T,tmv::RowMajor,4>();
    TestSmallNonSquareDiv_B1b<T,tmv::RowMajor,5>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixDiv_B1b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixDiv_B1b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixDiv_B1b<long double>();
#endif
