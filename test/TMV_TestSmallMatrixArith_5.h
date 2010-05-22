
#include "TMV_Test.h"
#include "TMV_Test3.h"

#include "TMV.h"
#include <fstream>

#ifdef NONSQUARE
#define NOMULTEQ
#endif

#include "TMV_TestMatrixArith.h"
#define CT std::complex<T>

template <class T> void TestSmallMatrixArith_5a();
template <class T> void TestSmallMatrixArith_5b();
template <class T> void TestSmallMatrixArith_5c();
template <class T> void TestSmallMatrixArith_5d();
template <class T> void TestSmallMatrixArith_5e();
template <class T> void TestSmallMatrixArith_5f();
template <class T> void TestSmallMatrixArith_5g();
template <class T> void TestSmallMatrixArith_5h();
template <class T> void TestSmallMatrixArith_5i();
template <class T> void TestSmallMatrixArith_5j();

template <class T, int M, int N, int K> 
void TestSmallMatrixArith_5(std::string label)
{
    if (showstartdone) {
        std::cout<<"Start Arith_5: M,N,K = "<<M<<','<<N<<','<<K<<std::endl;
    }

    tmv::SmallMatrix<T,M,N,tmv::RowMajor> a1rx;
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
        a1rx(i,j) = T(2.9+4.3*i-5.1*j);
    }
    a1rx.diag(0,0,std::min(M,N)).addToAll(T(6));
    a1rx(0,0) = 14; 
    if (M > 1) a1rx(1,0) = -2; 
    if (M > 2) a1rx(2,0) = 7; 
    if (M > 3) a1rx(3,0) = -10;
    if (M > 2 && N > 2) a1rx(2,2) = 30;

    tmv::SmallMatrix<CT,M,N,tmv::RowMajor> ca1rx = CT(3,2)*a1rx;
    if (M > 2 && N > 3) ca1rx(2,3) += CT(2.4,3.7);
    if (M > 1) ca1rx(1,0) *= CT(0.8,2.8);
    if (N > 1) ca1rx.col(1) *= CT(-1.1,3.6);
    if (M > 3) ca1rx.row(3) += tmv::SmallVector<CT,N>(CT(1.8,9.2));

    tmv::SmallMatrix<T,N,K,tmv::RowMajor> a2rx;
    for(int i=0;i<N;++i) for(int j=0;j<K;++j) {
        a2rx(i,j) = T(2.9+4.3*i-5.1*j);
    }
    a2rx.diag(0,0,std::min(N,K)).addToAll(T(6));
    a2rx(0,0) = 14; 
    if (N > 1) a2rx(1,0) = -2; 
    if (N > 2) a2rx(2,0) = 7; 
    if (N > 3) a2rx(3,0) = -10;
    if (N > 2 && K > 2) a2rx(2,2) = 30;
    if (N > 1) a2rx.row(1) *= T(3.1);
    if (K > 2) a2rx.col(2) -= tmv::SmallVector<T,N>(4.9);
    tmv::SmallMatrix<CT,N,K,tmv::RowMajor> ca2rx = CT(3,2)*a2rx;
    if (N > 2 && K > 3) ca2rx(2,3) += CT(2.4,3.7);
    if (N > 1) ca2rx(1,0) *= CT(0.8,2.8);
    if (K > 1) ca2rx.col(1) *= CT(-1.1,3.6);
    if (N > 3) ca2rx.row(3) += tmv::SmallVector<CT,K>(CT(1.8,9.2));
    ca2rx -= T(1.3)*a2rx;
    ca2rx *= CT(1.1,-2.5);

    tmv::SmallMatrix<T,M,N,tmv::ColMajor> a1cx = a1rx;
    tmv::SmallMatrix<T,N,K,tmv::ColMajor> a2cx = a2rx;
    tmv::SmallMatrix<CT,M,N,tmv::ColMajor> ca1cx = ca1rx;
    tmv::SmallMatrix<CT,N,K,tmv::ColMajor> ca2cx = ca2rx;

    tmv::SmallMatrixView<T,M,N,N,1> a1r = a1rx.view();
    tmv::SmallMatrixView<T,M,N,1,M> a1c = a1cx.view();
    tmv::SmallMatrixView<CT,M,N,N,1> ca1r = ca1rx.view();
    tmv::SmallMatrixView<CT,M,N,1,M> ca1c = ca1cx.view();
    tmv::SmallMatrixView<T,N,K,K,1> a2r = a2rx.view();
    tmv::SmallMatrixView<T,N,K,1,N> a2c = a2cx.view();
    tmv::SmallMatrixView<CT,N,K,K,1> ca2r = ca2rx.view();
    tmv::SmallMatrixView<CT,N,K,1,N> ca2c = ca2cx.view();

    TestMatrixArith5<T>(a1r,ca1r,a2r,ca2r,label+" RM/RM");
    TestMatrixArith5<T>(a1c,ca1c,a2c,ca2c,label+" CM/CM");
    TestMatrixArith5<T>(a1r,ca1r,a2c,ca2c,label+" RM/CM");
    TestMatrixArith5<T>(a1c,ca1c,a2r,ca2r,label+" CM/RM");

#if (XTEST & 1)
    tmv::SmallMatrix<T,3*M,4*N,tmv::ColMajor> a1nx;
    tmv::SmallMatrixView<T,M,N,3,12*M> a1n = a1nx.subMatrix(0,3*M,0,4*N,3,4);
    a1n = a1rx;
    tmv::SmallMatrix<CT,3*M,4*N,tmv::ColMajor> ca1nx;
    tmv::SmallMatrixView<CT,M,N,3,12*M> ca1n = ca1nx.subMatrix(0,3*M,0,4*N,3,4);
    ca1n = ca1rx;
    tmv::SmallMatrix<T,3*N,4*K,tmv::ColMajor> a2nx;
    tmv::SmallMatrixView<T,N,K,3,12*N> a2n = a2nx.subMatrix(0,3*N,0,4*K,3,4);
    a2n = a2rx;
    tmv::SmallMatrix<CT,3*N,4*K,tmv::ColMajor> ca2nx;
    tmv::SmallMatrixView<CT,N,K,3,12*N> ca2n = ca2nx.subMatrix(0,3*N,0,4*K,3,4);
    ca2n = ca2rx;

    TestMatrixArith5<T>(a1n,ca1n,a2r,ca2r,label+" NM/RM");
    TestMatrixArith5<T>(a1r,ca1r,a2n,ca2n,label+" RM/NM");
#if (XTEST & 2)
    TestMatrixArith5<T>(a1n,ca1n,a2c,ca2c,label+" NM/CM");
    TestMatrixArith5<T>(a1c,ca1c,a2n,ca2n,label+" CM/NM");
#endif
#endif
}

