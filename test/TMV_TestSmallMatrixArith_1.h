
#include "TMV_Test.h"
#include "TMV_Test_3.h"

#include "TMV.h"
#include <fstream>

#include "TMV_TestMatrixArith.h"
#define CT std::complex<T>

template <class T> void TestSmallMatrixArith_1a();
template <class T> void TestSmallMatrixArith_1b();
template <class T> void TestSmallMatrixArith_1c();
template <class T> void TestSmallMatrixArith_1d();
template <class T> void TestSmallMatrixArith_1e();
template <class T> void TestSmallMatrixArith_1f();
template <class T> void TestSmallMatrixArith_1g();

template <class T, int M, int N> void TestSmallMatrixArith_1(std::string label)
{
    if (showstartdone) {
        std::cout<<"Start Arith_1: M,N = "<<M<<','<<N<<std::endl;
    }

    tmv::SmallMatrix<T,M,N,tmv::RowMajor> a1x;
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
        a1x(i,j) = T(2.9+4.3*i-5.1*j);
    }
    a1x.diag(0,0,std::min(M,N)).addToAll(T(6));
    a1x(0,0) = 14; 
    if (M > 1) a1x(1,0) = -2; 
    if (M > 2) a1x(2,0) = 7; 
    if (M > 3) a1x(3,0) = -10;
    if (M > 2 && N > 2) a1x(2,2) = 30;

    tmv::SmallMatrix<CT,M,N,tmv::RowMajor> ca1x = CT(3,2)*a1x;
    if (M > 2 && N > 3) ca1x(2,3) += CT(2.4,3.7);
    if (M > 1) ca1x(1,0) *= CT(0.8,2.8);
    if (N > 1) ca1x.col(1) *= CT(-1.1,3.6);
    if (M > 3) ca1x.row(3) += tmv::SmallVector<CT,N>(CT(1.8,9.2));

#ifndef NONSQUARE
    // These next two is to make sure Det is calculable without overflow.
    if (N > 10) {
        a1x /= T(N*N); a1x += T(1);
        ca1x /= T(N*N); ca1x += T(1);
    }
#endif

    tmv::SmallMatrix<T,N,M,tmv::ColMajor> a2x = a1x.transpose();
    if (N > 1) a2x.row(1) *= T(3.1);
    if (M > 2) a2x.col(2) -= tmv::SmallVector<T,N>(4.9);
    tmv::SmallMatrix<CT,N,M,tmv::ColMajor> ca2x = ca1x.transpose();
    ca2x -= T(1.3)*a2x;
    ca2x *= CT(1.1,-2.5);

    tmv::SmallMatrixView<T,M,N,N,1> a1 = a1x.view();
    tmv::SmallMatrixView<T,N,M,1,N> a2 = a2x.view();
    tmv::SmallMatrixView<CT,M,N,N,1> ca1 = ca1x.view();
    tmv::SmallMatrixView<CT,N,M,1,N> ca2 = ca2x.view();

    TestMatrixArith1<T>(a1,ca1,label+" ColMajor");
    TestMatrixArith1<T>(a2,ca2,label+" RowMajor");

#if (XTEST & 1)
    tmv::SmallMatrix<T,3*M,4*N,tmv::ColMajor> a3x;
    tmv::SmallMatrixView<T,M,N> a3 = a3x.subMatrix(0,3*M,0,4*N,3,4);
    a3 = a1;
    tmv::SmallMatrix<CT,3*M,4*N,tmv::ColMajor> ca3x;
    tmv::SmallMatrixView<CT,M,N> ca3 = ca3x.subMatrix(0,3*M,0,4*N,3,4);
    ca3 = ca1;

    TestMatrixArith1<T>(a3,ca3,label+" NonMajor");
#endif
}
