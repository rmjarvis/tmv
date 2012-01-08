
#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#define NOLDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::StorageType stor, int N> 
static void TestSmallNonSquareDiv_B4a()
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

    tmv::SmallMatrix<T,N,7,stor> a5;
    for(int i=0;i<N;++i) for(int j=0;j<7;++j) a5(i,j) = T(1-3*i+2*j);
    a5.subMatrix(0,N,1,7) -= a1.transpose();
    tmv::SmallMatrix<std::complex<T>,N,7,stor> c5 = a5 * std::complex<T>(2,-2);
    c5.subMatrix(0,N,1,7) -= c1.transpose();
    c5.col(3) *= std::complex<T>(-1,3);
    c5.row(0).addToAll(std::complex<T>(1,9));

    tmv::SmallMatrix<T,6,7,stor> a7;
    for(int i=0;i<6;++i) for(int j=0;j<7;++j) a7(i,j) = T(5+2*i-2*j);
    a7.subMatrix(0,6,1,1+N) -= T(2)*a1;
    tmv::SmallMatrix<std::complex<T>,6,7,stor> c7 = a7 * std::complex<T>(-1,-3);
    c7.subMatrix(0,6,1,1+N) -= T(2)*c1;
    c7.col(3) *= std::complex<T>(-1,3);
    c7.row(4).addToAll(std::complex<T>(1,9));

    TestMatrixDivArith3b(tmv::QR,a1,a7,a5,c1,c7,c5,"NonSquare/NonSquare"); 
#if (XTEST & 32)
    tmv::SmallMatrix<T,6,N,stor|tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,6,N,stor|tmv::FortranStyle> c1f = c1;
    tmv::SmallMatrix<T,N,7,stor|tmv::FortranStyle> a5f = a5;
    tmv::SmallMatrix<std::complex<T>,N,7,stor|tmv::FortranStyle> c5f = c5;
    tmv::SmallMatrix<T,6,7,stor|tmv::FortranStyle> a7f = a7;
    tmv::SmallMatrix<std::complex<T>,6,7,stor|tmv::FortranStyle> c7f = c7;

    TestMatrixDivArith3b(tmv::QR,a1f,a7,a5,c1f,c7,c5,"NonSquare/NonSquare"); 
    TestMatrixDivArith3b(tmv::QR,a1f,a7f,a5,c1f,c7f,c5,"NonSquare/NonSquare"); 
    TestMatrixDivArith3b(tmv::QR,a1f,a7f,a5f,c1f,c7f,c5f,"NonSquare/NonSquare"); 
#endif
}

template <class T> 
void TestSmallMatrixDiv_B4a()
{
    TestSmallNonSquareDiv_B4a<T,tmv::ColMajor,2>();
    TestSmallNonSquareDiv_B4a<T,tmv::ColMajor,5>();
#if (XTEST & 2)
    TestSmallNonSquareDiv_B4a<T,tmv::ColMajor,1>();
    TestSmallNonSquareDiv_B4a<T,tmv::ColMajor,3>();
    TestSmallNonSquareDiv_B4a<T,tmv::ColMajor,4>();
    TestSmallNonSquareDiv_B4a<T,tmv::RowMajor,1>();
    TestSmallNonSquareDiv_B4a<T,tmv::RowMajor,2>();
    TestSmallNonSquareDiv_B4a<T,tmv::RowMajor,3>();
    TestSmallNonSquareDiv_B4a<T,tmv::RowMajor,4>();
    TestSmallNonSquareDiv_B4a<T,tmv::RowMajor,5>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixDiv_B4a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixDiv_B4a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixDiv_B4a<long double>();
#endif
