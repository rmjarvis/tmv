
#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::StorageType stor, int N> 
static void TestSmallSquareDiv_A5()
{
    tmv::SmallMatrix<T,N,N,stor> m;

    for(int i=0;i<N;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
    m.diag() *= T(11);
    m /= T(7);
    if (N > 1) m(1,0) = -2;
    if (N > 2) m(2,0) = 7;
    if (N > 3) m(3,0) = -10;
    if (N > 2) m(2,2) = 30;

    tmv::SmallMatrix<std::complex<T>,N,N,stor> c = m*std::complex<T>(-1,4);

    tmv::SmallMatrix<T,N,N,stor> a1 = m;
    tmv::SmallMatrix<std::complex<T>,N,N,stor> c1 = a1 * std::complex<T>(1,2);
    c1.diag().addToAll(std::complex<T>(3,1));

    tmv::SmallMatrix<T,N,7,stor> a3;
    for(int i=0;i<N;++i) for(int j=0;j<7;++j) a3(i,j) = T(1-3*i+2*j);
    tmv::SmallMatrix<T,7,N,stor> a4 = a3.transpose();
    a4.subMatrix(1,1+N,0,N) -= a1;

    tmv::SmallMatrix<std::complex<T>,N,7,stor> c3 = a3*std::complex<T>(1,2);
    tmv::SmallMatrix<std::complex<T>,7,N,stor> c4 = c3.adjoint();
    c4.subMatrix(1,1+N,0,N) -= c1;

    tmv::SmallMatrix<T,7,N,stor> a4b;
    tmv::SmallMatrix<std::complex<T>,7,N,stor> c4b;

    TestMatrixDivArith3c(tmv::LU,a1,a4,a4b,c1,c4,c4b,"NonSquare/Square");
#if (XTEST & 32)
    tmv::SmallMatrix<T,N,N,stor|tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,stor|tmv::FortranStyle> c1f = c1;

    tmv::SmallMatrix<T,7,N,stor|tmv::FortranStyle> a4f = a4;
    tmv::SmallMatrix<std::complex<T>,7,N,stor|tmv::FortranStyle> c4f = c4;
    tmv::SmallMatrix<T,7,N,stor|tmv::FortranStyle> a4fb = a4;
    tmv::SmallMatrix<std::complex<T>,7,N,stor|tmv::FortranStyle> c4fb = c4;

    TestMatrixDivArith3c(tmv::LU,a1f,a4,a4b,c1f,c4,c4b,"NonSquare/Square");
    TestMatrixDivArith3c(tmv::LU,a1f,a4f,a4b,c1f,c4f,c4b,"NonSquare/Square");
    TestMatrixDivArith3c(tmv::LU,a1f,a4f,a4fb,c1f,c4f,c4fb,"NonSquare/Square");
#endif
}

template <class T> 
void TestSmallMatrixDiv_A5()
{
    TestSmallSquareDiv_A5<T,tmv::ColMajor,2>();
    TestSmallSquareDiv_A5<T,tmv::ColMajor,5>();
#if (XTEST & 2)
    TestSmallSquareDiv_A5<T,tmv::ColMajor,1>();
    TestSmallSquareDiv_A5<T,tmv::ColMajor,3>();
    TestSmallSquareDiv_A5<T,tmv::ColMajor,4>();
    TestSmallSquareDiv_A5<T,tmv::RowMajor,1>();
    TestSmallSquareDiv_A5<T,tmv::RowMajor,2>();
    TestSmallSquareDiv_A5<T,tmv::RowMajor,3>();
    TestSmallSquareDiv_A5<T,tmv::RowMajor,4>();
    TestSmallSquareDiv_A5<T,tmv::RowMajor,5>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixDiv_A5<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixDiv_A5<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixDiv_A5<long double>();
#endif
