
#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::StorageType stor, int N> 
static void TestSmallNonSquareDiv_B5a()
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

    tmv::SmallMatrix<T,7,N,stor> a4;
    for(int i=0;i<7;++i) for(int j=0;j<N;++j) a4(i,j) = T(1-3*i+2*j);
    a4.subMatrix(0,6,0,N) += a1;
    tmv::SmallMatrix<std::complex<T>,7,N,stor> c4 = a4*std::complex<T>(1,2);
    c4.subMatrix(0,6,0,N) += c1;
    if (N > 1) c4.col(1) *= std::complex<T>(2,1);
    c4.row(2).addToAll(std::complex<T>(-7,2));

    tmv::SmallMatrix<T,7,6,stor> a6;
    for(int i=0;i<7;++i) for(int j=0;j<6;++j) a6(i,j) = T(5+2*i-2*j);
    a6.subMatrix(1,7,1,1+N) += a1;
    tmv::SmallMatrix<std::complex<T>,7,6,stor> c6 = a6*std::complex<T>(1,2);
    c6.subMatrix(1,7,1,1+N) += c1;
    c6.col(1) *= std::complex<T>(2,1);
    c6.row(5).addToAll(std::complex<T>(-7,2));

    TestMatrixDivArith3c(tmv::QR,a1,a4,a6,c1,c4,c6,"NonSquare/NonSquare"); 
#if (XTEST & 32)
    tmv::SmallMatrix<T,6,N,stor|tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,6,N,stor|tmv::FortranStyle> c1f = c1;
    tmv::SmallMatrix<T,7,N,stor|tmv::FortranStyle> a4f = a4;
    tmv::SmallMatrix<std::complex<T>,7,N,stor|tmv::FortranStyle> c4f = c4;
    tmv::SmallMatrix<T,7,6,stor|tmv::FortranStyle> a6f = a6;
    tmv::SmallMatrix<std::complex<T>,7,6,stor|tmv::FortranStyle> c6f = c6;

    TestMatrixDivArith3c(tmv::QR,a1f,a4,a6,c1f,c4,c6,"NonSquare/NonSquare"); 
    TestMatrixDivArith3c(tmv::QR,a1f,a4f,a6,c1f,c4f,c6,"NonSquare/NonSquare"); 
    TestMatrixDivArith3c(tmv::QR,a1f,a4f,a6f,c1f,c4f,c6f,"NonSquare/NonSquare"); 

#endif
}

template <class T> 
void TestSmallMatrixDiv_B5a()
{
    TestSmallNonSquareDiv_B5a<T,tmv::ColMajor,2>();
    TestSmallNonSquareDiv_B5a<T,tmv::ColMajor,5>();
#if (XTEST & 2)
    TestSmallNonSquareDiv_B5a<T,tmv::ColMajor,1>();
    TestSmallNonSquareDiv_B5a<T,tmv::ColMajor,3>();
    TestSmallNonSquareDiv_B5a<T,tmv::ColMajor,4>();
    TestSmallNonSquareDiv_B5a<T,tmv::RowMajor,1>();
    TestSmallNonSquareDiv_B5a<T,tmv::RowMajor,2>();
    TestSmallNonSquareDiv_B5a<T,tmv::RowMajor,3>();
    TestSmallNonSquareDiv_B5a<T,tmv::RowMajor,4>();
    TestSmallNonSquareDiv_B5a<T,tmv::RowMajor,5>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixDiv_B5a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixDiv_B5a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixDiv_B5a<long double>();
#endif
