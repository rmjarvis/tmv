#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test3.h"
#define NONSQUARE
#include "TMV_TestMatrixArith.h"

template <class T, int N> 
static void DoTestSmallMatrixArith_B1()
{
    tmv::SmallMatrix<T,N,N,tmv::RowMajor> a1;
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) {
        a1(i,j) = T(2.9+4.3*i-5.1*j);
    }
    a1.diag().addToAll(T(6));
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

    tmv::SmallMatrix<T,7,N> a3x;
    tmv::SmallMatrix<std::complex<T>,7,N> ca3x;
    tmv::SmallMatrix<T,N,7> a5x;
    tmv::SmallMatrix<std::complex<T>,N,7> ca5x;

    if (showstartdone) {
        std::cout<<"B1\n";
    }
    TestMatrixArith1<T>(a3x,ca3x,a3,ca3,"NonSquare");
    TestMatrixArith1<T>(a5x,ca5x,a5,ca5,"NonSquare");
#ifdef XTEST
    tmv::SmallMatrix<T,7,N,tmv::ColMajor> a4 = a3;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::ColMajor> ca4 = ca3;
    tmv::SmallMatrix<T,N,7,tmv::ColMajor> a6 = a5;
    tmv::SmallMatrix<std::complex<T>,N,7,tmv::ColMajor> ca6 = ca5;

    tmv::SmallMatrix<T,7,N,tmv::RowMajor,tmv::FortranStyle> a3f = a3;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::RowMajor,tmv::FortranStyle> ca3f = ca3;
    tmv::SmallMatrix<T,7,N,tmv::ColMajor,tmv::FortranStyle> a4f = a4;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::ColMajor,tmv::FortranStyle> ca4f = ca4;
    tmv::SmallMatrix<T,N,7,tmv::RowMajor,tmv::FortranStyle> a5f = a5;
    tmv::SmallMatrix<std::complex<T>,N,7,tmv::RowMajor,tmv::FortranStyle> ca5f = ca5;
    tmv::SmallMatrix<T,N,7,tmv::ColMajor,tmv::FortranStyle> a6f = a6;
    tmv::SmallMatrix<std::complex<T>,N,7,tmv::ColMajor,tmv::FortranStyle> ca6f = ca6;

    TestMatrixArith1<T>(a3x,ca3x,a4,ca4,"NonSquare");
    TestMatrixArith1<T>(a5x,ca5x,a6,ca6,"NonSquare");
    TestMatrixArith1<T>(a3x,ca3x,a3f,ca3f,"NonSquare");
    TestMatrixArith1<T>(a3x,ca3x,a4f,ca4f,"NonSquare");
    TestMatrixArith1<T>(a5x,ca5x,a5f,ca5f,"NonSquare");
    TestMatrixArith1<T>(a5x,ca5x,a6f,ca6f,"NonSquare");
#endif
}

template <class T> 
void TestSmallMatrixArith_B1()
{
    DoTestSmallMatrixArith_B1<T,2>();
    DoTestSmallMatrixArith_B1<T,5>();
#ifdef XTEST
    DoTestSmallMatrixArith_B1<T,1>();
    DoTestSmallMatrixArith_B1<T,3>();
    DoTestSmallMatrixArith_B1<T,4>();
#endif
}


#ifdef INST_DOUBLE
template void TestSmallMatrixArith_B1<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrixArith_B1<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrixArith_B1<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrixArith_B1<int>();
#endif
