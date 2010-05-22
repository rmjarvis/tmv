
#include "TMV_Test.h"
#include "TMV_Test3.h"

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

template <class T, int M, int N> void TestSmallMatrixArith_4(std::string label)
{
    if (showstartdone) {
        std::cout<<"Start Arith_4: M,N = "<<M<<','<<N<<std::endl;
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

    tmv::SmallMatrix<T,M,N,tmv::RowMajor> a2rx = a1rx;
    if (M > 1) a2rx.row(1) *= T(3.1);
    if (N > 2) a2rx.col(2) -= tmv::SmallVector<T,M>(4.9);
    tmv::SmallMatrix<CT,M,N,tmv::RowMajor> ca2rx = ca1rx;
    ca2rx -= T(1.3)*a2rx;
    ca2rx *= CT(1.1,-2.5);

    tmv::SmallMatrix<T,M,N,tmv::ColMajor> a1cx = a1rx;
    tmv::SmallMatrix<T,M,N,tmv::ColMajor> a2cx = a2rx;
    tmv::SmallMatrix<CT,M,N,tmv::ColMajor> ca1cx = ca1rx;
    tmv::SmallMatrix<CT,M,N,tmv::ColMajor> ca2cx = ca2rx;

    tmv::SmallMatrixView<T,M,N,N,1> a1r = a1rx.view();
    tmv::SmallMatrixView<T,M,N,1,M> a1c = a1cx.view();
    tmv::SmallMatrixView<CT,M,N,N,1> ca1r = ca1rx.view();
    tmv::SmallMatrixView<CT,M,N,1,M> ca1c = ca1cx.view();
    tmv::SmallMatrixView<T,M,N,N,1> a2r = a2rx.view();
    tmv::SmallMatrixView<T,M,N,1,M> a2c = a2cx.view();
    tmv::SmallMatrixView<CT,M,N,N,1> ca2r = ca2rx.view();
    tmv::SmallMatrixView<CT,M,N,1,M> ca2c = ca2cx.view();

    TestMatrixArith4<T>(a1r,ca1r,a2r,ca2r,label+" RM/RM");
    TestMatrixArith4<T>(a1c,ca1c,a2c,ca2c,label+" CM/CM");
    TestMatrixArith4<T>(a1r,ca1r,a2c,ca2c,label+" RM/CM");
    TestMatrixArith4<T>(a1c,ca1c,a2r,ca2r,label+" CM/RM");

#if (XTEST & 1)
    tmv::SmallMatrix<T,3*M,4*N,tmv::ColMajor> anx;
    tmv::SmallMatrixView<T,M,N,3,12*M> an = anx.subMatrix(0,3*M,0,4*N,3,4);
    an = a1rx;
    tmv::SmallMatrix<CT,3*M,4*N,tmv::ColMajor> canx;
    tmv::SmallMatrixView<CT,M,N,3,12*M> can = canx.subMatrix(0,3*M,0,4*N,3,4);
    can = ca1rx;

    TestMatrixArith4<T>(an,can,a2r,ca2r,label+" NM/RM");
    TestMatrixArith4<T>(a1r,ca1r,an,can,label+" RM/NM");
#if (XTEST & 2)
    TestMatrixArith4<T>(an,can,a2c,ca2c,label+" NM/CM");
    TestMatrixArith4<T>(a1c,ca1c,an,can,label+" CM/NM");
#endif
#endif
}
