#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"

#define INORDER
#define NOELEMMULT

#include "TMV_TestMatrixArith.h"

template <class T, int N> 
static void DoTestSmallMatrixArith_B6a()
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
    tmv::SmallMatrix<T,7,N,tmv::ColMajor> a4 = a3;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::ColMajor> ca4 = ca3;

    if (showstartdone) {
        std::cout<<"B6a"<<std::endl;
    }
    TestMatrixArith6(a3,ca3,a1,ca1,a4,ca4,"NonSquare");
#if (XTEST & 2)
    tmv::SmallMatrix<T,7,N,tmv::RowMajor> a3b = a3;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::RowMajor> ca3b = ca3;
    tmv::SmallMatrix<T,7,N,tmv::ColMajor> a4b = a4;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::ColMajor> ca4b = ca4;

    tmv::SmallMatrix<T,N,N,tmv::ColMajor> a2 = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::ColMajor> ca2 = ca1;
    TestMatrixArith6(a3,ca3,a1,ca1,a3b,ca3b,"NonSquare");
    TestMatrixArith6(a3,ca3,a2,ca2,a3b,ca3b,"NonSquare");
    TestMatrixArith6(a3,ca3,a2,ca2,a4,ca4,"NonSquare");
    TestMatrixArith6(a4,ca4,a1,ca1,a3,ca3,"NonSquare");
    TestMatrixArith6(a4,ca4,a1,ca1,a4b,ca4b,"NonSquare");
    TestMatrixArith6(a4,ca4,a2,ca2,a3,ca3,"NonSquare");
    TestMatrixArith6(a4,ca4,a2,ca2,a4b,ca4b,"NonSquare");
#endif

#if (XTEST & 32)
    tmv::SmallMatrix<T,N,N,tmv::RowMajor|tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor|tmv::FortranStyle> ca1f = ca1;
    tmv::SmallMatrix<T,7,N,tmv::RowMajor|tmv::FortranStyle> a3f = a3;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::RowMajor|tmv::FortranStyle> ca3f = ca3;
    tmv::SmallMatrix<T,7,N,tmv::ColMajor|tmv::FortranStyle> a4f = a4;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::ColMajor|tmv::FortranStyle> ca4f = ca4;

    TestMatrixArith6(a3f,ca3f,a1,ca1,a4,ca4,"NonSquare");
    TestMatrixArith6(a3f,ca3f,a1f,ca1f,a4,ca4,"NonSquare");
    TestMatrixArith6(a3f,ca3f,a1f,ca1f,a4f,ca4f,"NonSquare");
#endif
}

template <class T> 
void TestSmallMatrixArith_B6a()
{
    DoTestSmallMatrixArith_B6a<T,2>();
    DoTestSmallMatrixArith_B6a<T,5>();
#if (XTEST & 2)
    DoTestSmallMatrixArith_B6a<T,1>();
    DoTestSmallMatrixArith_B6a<T,3>();
    DoTestSmallMatrixArith_B6a<T,4>();
#endif
}


#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_B6a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_B6a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_B6a<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_B6a<int>();
#endif
