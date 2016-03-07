#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#define INORDER
#define NOMULTEQ
#include "TMV_TestMatrixArith.h"

template <class T, int N> 
static void DoTestSmallMatrixArith_B5c()
{
    tmv::SmallMatrix<T,N,N,tmv::RowMajor> a1;
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) {
        a1(i,j) = T(3+4*i-5*j);
    }
    a1(0,0) = 14;
    if (N > 1) a1(1,0) = -2;
    if (N > 2) a1(2,0) = 7;
    if (N > 3) a1(3,0) = -10;
    if (N > 2) a1(2,2) = 30;

    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor> ca1 = 
        std::complex<T>(3,2)*a1;
    if (N > 3) ca1(2,3) += std::complex<T>(2,4);
    if (N > 1) ca1(1,0) *= std::complex<T>(1,3);
    if (N > 1) ca1.col(1) *= std::complex<T>(-1,4);
    if (N > 3) ca1.row(3) += 
        tmv::SmallVector<std::complex<T>,N>(std::complex<T>(2,9));

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
        std::cout<<"B5c"<<std::endl;
    }
    TestMatrixArith5(a3,ca3,a5,ca5,"NonSquare");
#if (XTEST & 2)
    tmv::SmallMatrix<T,7,N,tmv::ColMajor> a4 = a3;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::ColMajor> ca4 = ca3;
    tmv::SmallMatrix<T,N,7,tmv::ColMajor> a6 = a5;
    tmv::SmallMatrix<std::complex<T>,N,7,tmv::ColMajor> ca6 = ca5;
    TestMatrixArith5(a3,ca3,a6,ca6,"NonSquare");
    TestMatrixArith5(a4,ca4,a5,ca5,"NonSquare");
    TestMatrixArith5(a4,ca4,a6,ca6,"NonSquare");
#endif

#if (XTEST & 32)
    tmv::SmallMatrix<T,7,N,tmv::RowMajor|tmv::FortranStyle> a3f = a3;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::RowMajor|tmv::FortranStyle> ca3f = ca3;
    tmv::SmallMatrix<T,N,7,tmv::RowMajor|tmv::FortranStyle> a5f = a5;
    tmv::SmallMatrix<std::complex<T>,N,7,tmv::RowMajor|tmv::FortranStyle> ca5f = ca5;

    TestMatrixArith5(a3f,ca3f,a5,ca5,"NonSquare");
    TestMatrixArith5(a3f,ca3f,a5f,ca5f,"NonSquare");
#endif
}

template <class T> 
void TestSmallMatrixArith_B5c()
{
    DoTestSmallMatrixArith_B5c<T,2>();
    DoTestSmallMatrixArith_B5c<T,5>();
#if (XTEST & 2)
    DoTestSmallMatrixArith_B5c<T,1>();
    DoTestSmallMatrixArith_B5c<T,3>();
    DoTestSmallMatrixArith_B5c<T,4>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_B5c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_B5c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_B5c<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_B5c<int>();
#endif
