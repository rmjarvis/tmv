
#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::StorageType stor, int N> 
static void TestSmallNonSquareDiv_B3b()
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
    tmv::SmallMatrix<T,N,N,tmv::RowMajor> a2b = m.transpose() * m;
    if (N > 1) a2b.row(1) *= T(3);
    if (N > 2) a2b.col(2).addToAll(-4);
    tmv::SmallMatrix<std::complex<T>,6,N,stor> c1 = a1 * std::complex<T>(1,2);
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor> c2b =
        a2b * std::complex<T>(-3,4);

    tmv::SmallMatrix<T,N,6,stor> a4;
    tmv::SmallMatrix<std::complex<T>,N,6,stor> c4;

    TestMatrixDivArith3c(tmv::QR,a1,a2b,a4,c1,c2b,c4,"Square/NonSquare"); 
#if (XTEST & 32)
    tmv::SmallMatrix<T,6,N,stor|tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,6,N,stor|tmv::FortranStyle> c1f = c1;

    tmv::SmallMatrix<T,N,N,tmv::RowMajor|tmv::FortranStyle> a2fb = a2b;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor|tmv::FortranStyle> c2fb = c2b;
    tmv::SmallMatrix<T,N,6,stor|tmv::FortranStyle> a4f = a4;
    tmv::SmallMatrix<std::complex<T>,N,6,stor|tmv::FortranStyle> c4f = c4;

    TestMatrixDivArith3c(tmv::QR,a1f,a2b,a4,c1f,c2b,c4,"Square/NonSquare"); 
    TestMatrixDivArith3c(tmv::QR,a1f,a2fb,a4,c1f,c2fb,c4,"Square/NonSquare"); 
    TestMatrixDivArith3c(tmv::QR,a1f,a2fb,a4f,c1f,c2fb,c4f,"Square/NonSquare"); 
#endif
}

template <class T> 
void TestSmallMatrixDiv_B3b()
{
    TestSmallNonSquareDiv_B3b<T,tmv::ColMajor,2>();
    TestSmallNonSquareDiv_B3b<T,tmv::ColMajor,5>();
#if (XTEST & 2)
    TestSmallNonSquareDiv_B3b<T,tmv::ColMajor,1>();
    TestSmallNonSquareDiv_B3b<T,tmv::ColMajor,3>();
    TestSmallNonSquareDiv_B3b<T,tmv::ColMajor,4>();
    TestSmallNonSquareDiv_B3b<T,tmv::RowMajor,1>();
    TestSmallNonSquareDiv_B3b<T,tmv::RowMajor,2>();
    TestSmallNonSquareDiv_B3b<T,tmv::RowMajor,3>();
    TestSmallNonSquareDiv_B3b<T,tmv::RowMajor,4>();
    TestSmallNonSquareDiv_B3b<T,tmv::RowMajor,5>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixDiv_B3b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixDiv_B3b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixDiv_B3b<long double>();
#endif
