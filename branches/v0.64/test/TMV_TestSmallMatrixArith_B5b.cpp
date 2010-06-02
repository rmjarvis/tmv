#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test3.h"
#define INORDER
#define NOMULTEQ
#include "TMV_TestMatrixArith.h"

template <class T, int N> 
static void DoTestSmallMatrixArith_B5b()
{
    tmv::SmallMatrix<T,N,N,tmv::RowMajor> a1;
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) {
        a1(i,j) = T(2.9+4.3*i-5.1*j);
    }
    a1(0,0) = 14;
    if (N > 1) a1(1,0) = -2;
    if (N > 2) a1(2,0) = 7;
    if (N > 3) a1(3,0) = -10;
    if (N > 2) a1(2,2) = 30;

    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor> ca1 = 
        std::complex<T>(3,2)*a1;
    if (N > 3) ca1(2,3) += std::complex<T>(2.4,3.7);
    if (N > 1) ca1(1,0) *= std::complex<T>(0.8,2.8);
    if (N > 1) ca1.col(1) *= std::complex<T>(-1.1,3.6);
    if (N > 3) ca1.row(3) += 
        tmv::SmallVector<std::complex<T>,N>(std::complex<T>(1.8,9.2));

    tmv::SmallMatrix<T,7,N,tmv::RowMajor> a3;
    for(int i=0;i<7;++i) for(int j=0;j<N;++j) a3(i,j) = T(1-3*i+2*j);
    a3.subMatrix(2,N+2,0,N) += a1;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::RowMajor> ca3 = 
        a3*std::complex<T>(1,2);
    ca3.subMatrix(2,N+2,0,N) += ca1;
    if (N > 1) ca3.col(1) *= std::complex<T>(2,1);
    ca3.row(6).addToAll(std::complex<T>(-7,2));

    tmv::SmallMatrix<T,N,7,tmv::RowMajor> a5 = a3.transpose();
    a5.subMatrix(0,N,1,N+1) -= a1;
    tmv::SmallMatrix<std::complex<T>,N,7,tmv::RowMajor> ca5 = ca3.adjoint();
    ca5.subMatrix(0,N,1,N+1) -= ca1;
    ca5.col(3) *= std::complex<T>(-1,3);
    ca5.row(0).addToAll(std::complex<T>(1,9));

    if (showstartdone) {
        std::cout<<"B5b\n";
    }
    TestMatrixArith5<T>(a1,ca1,a5,ca5,"NonSquare");
#ifdef XTEST
    tmv::SmallMatrix<T,N,N,tmv::ColMajor> a2 = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::ColMajor> ca2 = ca1;
    tmv::SmallMatrix<T,N,7,tmv::ColMajor> a6 = a5;
    tmv::SmallMatrix<std::complex<T>,N,7,tmv::ColMajor> ca6 = ca5;

    tmv::SmallMatrix<T,N,N,tmv::ColMajor,tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::ColMajor,tmv::FortranStyle> ca1f = ca1;
    tmv::SmallMatrix<T,N,7,tmv::RowMajor,tmv::FortranStyle> a5f = a5;
    tmv::SmallMatrix<std::complex<T>,N,7,tmv::RowMajor,tmv::FortranStyle> ca5f = ca5;

    TestMatrixArith5<T>(a1,ca1,a6,ca6,"NonSquare");
    TestMatrixArith5<T>(a2,ca2,a5,ca5,"NonSquare");
    TestMatrixArith5<T>(a2,ca2,a6,ca6,"NonSquare");
    TestMatrixArith5<T>(a1f,ca1f,a5,ca5,"NonSquare");
    TestMatrixArith5<T>(a1f,ca1f,a5f,ca5f,"NonSquare");
#endif
}

template <class T> 
void TestSmallMatrixArith_B5b()
{
    DoTestSmallMatrixArith_B5b<T,2>();
    DoTestSmallMatrixArith_B5b<T,5>();
#ifdef XTEST
    DoTestSmallMatrixArith_B5b<T,1>();
    DoTestSmallMatrixArith_B5b<T,3>();
    DoTestSmallMatrixArith_B5b<T,4>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_B5b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_B5b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_B5b<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_B5b<int>();
#endif
