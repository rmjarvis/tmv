
#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#define NOLDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::StorageType stor, int N> 
static void TestSmallNonSquareDiv_B2a()
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
    tmv::SmallMatrix<T,6,6,tmv::ColMajor> a3a = m * m.transpose();
    tmv::SmallMatrix<T,N,6,stor> a4;
    a3a.row(5) *= T(7);
    if (N > 3) a3a.col(3).addToAll(7);
    tmv::SmallMatrix<std::complex<T>,6,N,stor> c1 = a1 * std::complex<T>(1,2);
    tmv::SmallMatrix<std::complex<T>,6,6,tmv::ColMajor> c3a =
        a3a * std::complex<T>(-4,8);
    tmv::SmallMatrix<std::complex<T>,N,6,stor> c4;

    TestMatrixDivArith3b(tmv::QR,a1,a3a,a4,c1,c3a,c4,"Square/NonSquare"); 
#if (XTEST & 32)
    tmv::SmallMatrix<T,6,N,stor|tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,6,N,stor|tmv::FortranStyle> c1f = c1;

    tmv::SmallMatrix<T,6,6,tmv::ColMajor|tmv::FortranStyle> a3fa = a3a;
    tmv::SmallMatrix<std::complex<T>,6,6,tmv::ColMajor|tmv::FortranStyle> c3fa = c3a;
    tmv::SmallMatrix<T,N,6,stor|tmv::FortranStyle> a4f = a4;
    tmv::SmallMatrix<std::complex<T>,N,6,stor|tmv::FortranStyle> c4f = c4;

    TestMatrixDivArith3b(tmv::QR,a1f,a3a,a4,c1f,c3a,c4,"Square/NonSquare"); 
    TestMatrixDivArith3b(tmv::QR,a1f,a3fa,a4,c1f,c3fa,c4,"Square/NonSquare"); 
    TestMatrixDivArith3b(tmv::QR,a1f,a3fa,a4f,c1f,c3fa,c4f,"Square/NonSquare"); 
#endif
}

template <class T> 
void TestSmallMatrixDiv_B2a()
{
    TestSmallNonSquareDiv_B2a<T,tmv::ColMajor,2>();
    TestSmallNonSquareDiv_B2a<T,tmv::ColMajor,5>();
#if (XTEST & 2)
    TestSmallNonSquareDiv_B2a<T,tmv::ColMajor,1>();
    TestSmallNonSquareDiv_B2a<T,tmv::ColMajor,3>();
    TestSmallNonSquareDiv_B2a<T,tmv::ColMajor,4>();
    TestSmallNonSquareDiv_B2a<T,tmv::RowMajor,1>();
    TestSmallNonSquareDiv_B2a<T,tmv::RowMajor,2>();
    TestSmallNonSquareDiv_B2a<T,tmv::RowMajor,3>();
    TestSmallNonSquareDiv_B2a<T,tmv::RowMajor,4>();
    TestSmallNonSquareDiv_B2a<T,tmv::RowMajor,5>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixDiv_B2a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixDiv_B2a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixDiv_B2a<long double>();
#endif
