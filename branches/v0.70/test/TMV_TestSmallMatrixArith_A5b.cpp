#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV_TestMatrixArith.h"

template <class T, int N> 
static void DoTestSmallMatrixArith_A5b()
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

    tmv::SmallMatrix<T,N,N,tmv::ColMajor> a2 = a1.transpose();
    if (N > 1) a2.row(1) *= T(3);
    if (N > 2) a2.col(2) -= tmv::SmallVector<T,N>(4);
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::ColMajor> ca2 = ca1;
    ca2 -= T(0.13)*a2;
    ca2 *= std::complex<T>(1,-2);

    if (showstartdone) {
        std::cout<<"A5b"<<std::endl;
    }
    TestMatrixArith5(a2,ca2,a1,ca1,"Square");
#if (XTEST & 2)
    tmv::SmallMatrix<T,N,N,tmv::ColMajor> a2b = a2;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::ColMajor> ca2b = ca2;
    TestMatrixArith5(a2,ca2,a2b,ca2b,"Square");
#endif

#if (XTEST & 32)
    tmv::SmallMatrix<T,N,N,tmv::RowMajor|tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor|tmv::FortranStyle> ca1f = ca1;
    tmv::SmallMatrix<T,N,N,tmv::ColMajor|tmv::FortranStyle> a2f = a2;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::ColMajor|tmv::FortranStyle> ca2f = ca2;

    TestMatrixArith5(a2f,ca2f,a1,ca1,"Square");
    TestMatrixArith5(a2f,ca2f,a1f,ca1f,"Square");
#endif
}

template <class T> 
void TestSmallMatrixArith_A5b()
{
    DoTestSmallMatrixArith_A5b<T,2>();
    DoTestSmallMatrixArith_A5b<T,5>();
#if (XTEST & 2)
    DoTestSmallMatrixArith_A5b<T,1>();
    DoTestSmallMatrixArith_A5b<T,3>();
    DoTestSmallMatrixArith_A5b<T,4>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_A5b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_A5b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_A5b<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_A5b<int>();
#endif
