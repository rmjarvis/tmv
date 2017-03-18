
#include "TMV_Test.h"
#include "TMV_Test_3.h"

#include "TMV.h"
#include <fstream>

#ifdef NONSQUARE
#define NOMULTEQ
#endif

#include "TMV_TestMatrixArith.h"
#define CT std::complex<T>

template <class T> void TestSmallMatrixArith_6a();
template <class T> void TestSmallMatrixArith_6b();
template <class T> void TestSmallMatrixArith_6c();
template <class T> void TestSmallMatrixArith_6d();
template <class T> void TestSmallMatrixArith_6e();
template <class T> void TestSmallMatrixArith_6f();
template <class T> void TestSmallMatrixArith_6g();
template <class T> void TestSmallMatrixArith_6h();
template <class T> void TestSmallMatrixArith_6i();
template <class T> void TestSmallMatrixArith_6j();

template <class T, int M, int N, int K> 
void TestSmallMatrixArith_6(std::string label)
{
    if (showstartdone) {
        std::cout<<"Start Arith_6A: M,N,K = "<<M<<','<<N<<','<<K<<std::endl;
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

    tmv::SmallMatrix<T,N,K,tmv::RowMajor> a2r;
    for(int i=0;i<N;++i) for(int j=0;j<K;++j) {
        a2r(i,j) = T(2.9+4.3*i-5.1*j);
    }
    a2r.diag(0,0,std::min(N,K)).addToAll(T(6));
    a2r(0,0) = 14; 
    if (N > 1) a2r(1,0) = -2; 
    if (N > 2) a2r(2,0) = 7; 
    if (N > 3) a2r(3,0) = -10;
    if (N > 2 && K > 2) a2r(2,2) = 30;
    if (N > 1) a2r.row(1) *= T(3.1);
    if (K > 2) a2r.col(2) -= tmv::SmallVector<T,N>(4.9);
    tmv::SmallMatrix<CT,N,K,tmv::RowMajor> ca2r = CT(3,2)*a2r;
    if (N > 2 && K > 3) ca2r(2,3) += CT(2.4,3.7);
    if (N > 1) ca2r(1,0) *= CT(0.8,2.8);
    if (K > 1) ca2r.col(1) *= CT(-1.1,3.6);
    if (N > 3) ca2r.row(3) += tmv::SmallVector<CT,K>(CT(1.8,9.2));
    ca2r -= T(1.3)*a2r;
    ca2r *= CT(1.1,-2.5);

    tmv::SmallMatrix<T,M,N,tmv::ColMajor> a1c = a1r;
    tmv::SmallMatrix<T,N,K,tmv::ColMajor> a2c = a2r;
    tmv::SmallMatrix<CT,M,N,tmv::ColMajor> ca1c = ca1r;
    tmv::SmallMatrix<CT,N,K,tmv::ColMajor> ca2c = ca2r;

    TestMatrixArith6x(a1c,ca1c,a2c,ca2c,label+" CM/CM");
#if (XTEST & 2)
    TestMatrixArith6x(a1r,ca1r,a2r,ca2r,label+" RM/RM");
    TestMatrixArith6x(a1r,ca1r,a2c,ca2c,label+" RM/CM");
    TestMatrixArith6x(a1c,ca1c,a2r,ca2r,label+" CM/RM");
#endif

#if (XTEST & 1)
    tmv::SmallMatrix<T,3*M,4*N,tmv::ColMajor> a1nx;
    tmv::SmallMatrixView<T,M,N> a1n = a1nx.subMatrix(0,3*M,0,4*N,3,4);
    a1n = a1r;
    tmv::SmallMatrix<CT,3*M,4*N,tmv::ColMajor> ca1nx;
    tmv::SmallMatrixView<CT,M,N> ca1n = ca1nx.subMatrix(0,3*M,0,4*N,3,4);
    ca1n = ca1r;
    tmv::SmallMatrix<T,3*N,4*K,tmv::ColMajor> a2nx;
    tmv::SmallMatrixView<T,N,K> a2n = a2nx.subMatrix(0,3*N,0,4*K,3,4);
    a2n = a2r;
    tmv::SmallMatrix<CT,3*N,4*K,tmv::ColMajor> ca2nx;
    tmv::SmallMatrixView<CT,N,K> ca2n = ca2nx.subMatrix(0,3*N,0,4*K,3,4);
    ca2n = ca2r;

    TestMatrixArith6x(a1n,ca1n,a2r,ca2r,label+" NM/RM");
    TestMatrixArith6x(a1r,ca1r,a2n,ca2n,label+" RM/NM");
#if (XTEST & 2)
    TestMatrixArith6x(a1n,ca1n,a2c,ca2c,label+" NM/CM");
    TestMatrixArith6x(a1c,ca1c,a2n,ca2n,label+" CM/NM");
#endif
#endif

    tmv::SmallMatrix<T,M,K,tmv::ColMajor> a3c(T(0));
    tmv::SmallMatrix<CT,M,K,tmv::ColMajor> ca3c(CT(0));

    TestMatrixArith6(a1c,ca1c,a2c,ca2c,a3c,ca3c,label+" CM/CM/CM");
    TestMatrixArith6(a1r,ca1r,a2r,ca2r,a3c,ca3c,label+" RM/RM/CM");
    TestMatrixArith6(a1r,ca1r,a2c,ca2c,a3c,ca3c,label+" RM/CM/CM");
    TestMatrixArith6(a1c,ca1c,a2r,ca2r,a3c,ca3c,label+" CM/RM/CM");

#if (XTEST & 2)
    tmv::SmallMatrix<T,M,K,tmv::RowMajor> a3r(T(0));
    tmv::SmallMatrix<CT,M,K,tmv::RowMajor> ca3r(CT(0));

    TestMatrixArith6(a1c,ca1c,a2c,ca2c,a3r,ca3r,label+" CM/CM/RM");
    TestMatrixArith6(a1r,ca1r,a2r,ca2r,a3r,ca3r,label+" RM/RM/RM");
    TestMatrixArith6(a1r,ca1r,a2c,ca2c,a3r,ca3r,label+" RM/CM/RM");
    TestMatrixArith6(a1c,ca1c,a2r,ca2r,a3r,ca3r,label+" CM/RM/RM");
#endif
}
