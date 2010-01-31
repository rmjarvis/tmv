#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV_TestMatrixArith.h"

template <class T, int N> 
static void DoTestSmallMatrixArith_A5a()
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

    tmv::SmallMatrix<T,N,N,tmv::ColMajor> a2 = a1.transpose();
    if (N > 1) a2.row(1) *= T(3.1);
    if (N > 2) a2.col(2) -= tmv::SmallVector<T,N>(4.9);
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::ColMajor> ca2 = ca1;
    ca2 -= T(0.13)*a2;
    ca2 *= std::complex<T>(1.1,-2.5);

    tmv::SmallMatrix<T,N,N> a1x;
    tmv::SmallMatrix<std::complex<T>,N,N> ca1x;

    if (showstartdone) {
        std::cout<<"A5a\n";
    }
    TestMatrixArith5<T>(a1x,ca1x,a1x,ca1x,a1,ca1,a2,ca2,"Square");

#ifdef XTEST
    tmv::SmallMatrix<T,N,N,tmv::RowMajor> a1b = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor> ca1b = ca1;
    tmv::SmallMatrix<T,N,N,tmv::RowMajor,tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor,tmv::FortranStyle> ca1f = ca1;
    tmv::SmallMatrix<T,N,N,tmv::ColMajor,tmv::FortranStyle> a2f = a2;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::ColMajor,tmv::FortranStyle> ca2f = ca2;

    TestMatrixArith5<T>(a1x,ca1x,a1x,ca1x,a1,ca1,a1b,ca1b,"Square");
    TestMatrixArith5<T>(a1x,ca1x,a1x,ca1x,a1f,ca1f,a2,ca2,"Square");
    TestMatrixArith5<T>(a1x,ca1x,a1x,ca1x,a1f,ca1f,a2f,ca2f,"Square");
#endif
}

template <class T> 
void TestSmallMatrixArith_A5a()
{
    DoTestSmallMatrixArith_A5a<T,2>();
    DoTestSmallMatrixArith_A5a<T,5>();
#ifdef XTEST
    DoTestSmallMatrixArith_A5a<T,1>();
    DoTestSmallMatrixArith_A5a<T,3>();
    DoTestSmallMatrixArith_A5a<T,4>();
#endif
}


#ifdef INST_DOUBLE
template void TestSmallMatrixArith_A5a<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrixArith_A5a<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrixArith_A5a<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrixArith_A5a<int>();
#endif
