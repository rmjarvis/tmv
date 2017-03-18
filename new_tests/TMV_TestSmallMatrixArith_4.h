
#include "TMV_Test.h"
#include "TMV_Test_3.h"

#include "TMV.h"
#include <fstream>

#include "TMV_TestMatrixArith.h"
#define CT std::complex<T>

template <class T> void TestSmallMatrixArith_4a();
template <class T> void TestSmallMatrixArith_4b();
template <class T> void TestSmallMatrixArith_4c();
template <class T> void TestSmallMatrixArith_4d();
template <class T> void TestSmallMatrixArith_4e();
template <class T> void TestSmallMatrixArith_4f();
template <class T> void TestSmallMatrixArith_4g();

template <class T, int M, int N> void TestSmallMatrixArith_4(std::string label)
{
    if (showstartdone) {
        std::cout<<"Start Arith_4: M,N = "<<M<<','<<N<<std::endl;
    }

    tmv::SmallMatrix<T,M,N,tmv::RowMajor> a1r;
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
        a1r(i,j) = T(2.9+4.3*i-5.1*j);
    }
    a1r.diag(0,0,std::min(M,N)).addToAll(T(6));
    a1r(0,0) = 14; 
    if (M > 1) a1r(1,0) = -2; 
    if (M > 2) a1r(2,0) = 7; 
    if (M > 3) a1r(3,0) = -10;
    if (M > 2 && N > 2) a1r(2,2) = 30;

    tmv::SmallMatrix<CT,M,N,tmv::RowMajor> ca1r = CT(3,2)*a1r;
    if (M > 2 && N > 3) ca1r(2,3) += CT(2.4,3.7);
    if (M > 1) ca1r(1,0) *= CT(0.8,2.8);
    if (N > 1) ca1r.col(1) *= CT(-1.1,3.6);
    if (M > 3) ca1r.row(3) += tmv::SmallVector<CT,N>(CT(1.8,9.2));

    tmv::SmallMatrix<T,M,N,tmv::RowMajor> a2r = a1r;
    if (M > 1) a2r.row(1) *= T(3.1);
    if (N > 2) a2r.col(2) -= tmv::SmallVector<T,M>(4.9);
    tmv::SmallMatrix<CT,M,N,tmv::RowMajor> ca2r = ca1r;
    ca2r -= T(1.3)*a2r;
    ca2r *= CT(1.1,-2.5);

    tmv::SmallMatrix<T,M,N,tmv::ColMajor> a1c = a1r;
    tmv::SmallMatrix<T,M,N,tmv::ColMajor> a2c = a2r;
    tmv::SmallMatrix<CT,M,N,tmv::ColMajor> ca1c = ca1r;
    tmv::SmallMatrix<CT,M,N,tmv::ColMajor> ca2c = ca2r;

    TestMatrixArith4(a1r,ca1r,a2r,ca2r,label+" RM/RM");
    TestMatrixArith4(a1c,ca1c,a2c,ca2c,label+" CM/CM");
    TestMatrixArith4(a1r,ca1r,a2c,ca2c,label+" RM/CM");
    TestMatrixArith4(a1c,ca1c,a2r,ca2r,label+" CM/RM");

#if (XTEST & 1)
    tmv::SmallMatrix<T,3*M,4*N,tmv::ColMajor> a1nx;
    tmv::SmallMatrixView<T,M,N> a1n = a1nx.subMatrix(0,3*M,0,4*N,3,4);
    a1n = a1r;
    tmv::SmallMatrix<CT,3*M,4*N,tmv::ColMajor> ca1nx;
    tmv::SmallMatrixView<CT,M,N> ca1n = ca1nx.subMatrix(0,3*M,0,4*N,3,4);
    ca1n = ca1r;

    TestMatrixArith4(a1n,ca1n,a2r,ca2r,label+" NM/RM");
    TestMatrixArith4(a1r,ca1r,a1n,ca1n,label+" RM/NM");
#if (XTEST & 2)
    TestMatrixArith4(a1n,ca1n,a2c,ca2c,label+" NM/CM");
    TestMatrixArith4(a1c,ca1c,a1n,ca1n,label+" CM/NM");
#endif
#endif
}
