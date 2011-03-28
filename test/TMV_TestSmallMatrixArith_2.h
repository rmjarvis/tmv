
#include "TMV_Test.h"
#include "TMV_Test_3.h"

#include "TMV.h"
#include <fstream>

#include "TMV_TestMatrixArith.h"
#define CT std::complex<T>

template <class T> void TestSmallMatrixArith_2a();
template <class T> void TestSmallMatrixArith_2b();
template <class T> void TestSmallMatrixArith_2c();
template <class T> void TestSmallMatrixArith_2d();
template <class T> void TestSmallMatrixArith_2e();
template <class T> void TestSmallMatrixArith_2f();
template <class T> void TestSmallMatrixArith_2g();

template <class T, int M, int N> void TestSmallMatrixArith_2(std::string label)
{
    if (showstartdone) {
        std::cout<<"Start Arith_2: M,N = "<<M<<','<<N<<std::endl;
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

    tmv::SmallMatrix<T,N,M,tmv::ColMajor> a2x = a1x.transpose();
    if (N > 1) a2x.row(1) *= T(3.1);
    if (M > 2) a2x.col(2) -= tmv::SmallVector<T,N>(4.9);
    tmv::SmallMatrix<CT,N,M,tmv::ColMajor> ca2x = ca1x.transpose();
    ca2x -= T(1.3)*a2x;
    ca2x *= CT(1.1,-2.5);

    tmv::SmallVector<T,N> v1x = a1x.row(0);
    tmv::SmallVector<T,M> v2x = a1x.col(0);
    tmv::SmallVector<CT,N> cv1x = ca1x.row(0);
    tmv::SmallVector<CT,M> cv2x = ca1x.col(0);
    tmv::SmallVectorView<T,N,1> v1 = v1x.view();
    tmv::SmallVectorView<T,M,1> v2 = v2x.view();
    tmv::SmallVectorView<CT,N,1> cv1 = cv1x.view();
    tmv::SmallVectorView<CT,M,1> cv2 = cv2x.view();

    tmv::SmallMatrixView<T,M,N,N,1> a1 = a1x.view();
    tmv::SmallMatrixView<T,N,M,1,N> a2 = a2x.view();
    tmv::SmallMatrixView<CT,M,N,N,1> ca1 = ca1x.view();
    tmv::SmallMatrixView<CT,N,M,1,N> ca2 = ca2x.view();

    TestMatrixArith2a<T>(a1,ca1,v1,cv1,v2,cv2,label+" ColMajor");
    TestMatrixArith2a<T>(a2,ca2,v2,cv2,v1,cv1,label+" RowMajor");

    TestMatrixArith2b<T>(a1,ca1,v1,cv1,v2,cv2,label+" ColMajor");
    TestMatrixArith2b<T>(a2,ca2,v2,cv2,v1,cv1,label+" RowMajor");

#if (XTEST & 1)
    tmv::SmallMatrix<T,3*M,4*N,tmv::ColMajor> a3x;
    tmv::SmallMatrixView<T,M,N,3,12*M> a3 = a3x.subMatrix(0,3*M,0,4*N,3,4);
    a3 = a1;
    tmv::SmallMatrix<CT,3*M,4*N,tmv::ColMajor> ca3x;
    tmv::SmallMatrixView<CT,M,N,3,12*M> ca3 = ca3x.subMatrix(0,3*M,0,4*N,3,4);
    ca3 = ca1;

    tmv::SmallVector<T,5*N> v15;
    tmv::SmallVectorView<T,N,5> v1s = v15.subVector(0,5*N,5);
    v1s = v1;
    tmv::SmallVector<T,5*M> v25;
    tmv::SmallVectorView<T,M,5> v2s = v25.subVector(0,5*M,5);
    v2s = v2;
    tmv::SmallVector<CT,5*N> cv15;
    tmv::SmallVectorView<CT,N,5> cv1s = cv15.subVector(0,5*N,5);
    cv1s = cv1;
    tmv::SmallVector<CT,5*M> cv25;
    tmv::SmallVectorView<CT,M,5> cv2s = cv25.subVector(0,5*M,5);
    cv2s = cv2;

    TestMatrixArith2a<T>(a1,ca1,v1s,cv1s,v2,cv2,label+" ColMajor Step1");
    TestMatrixArith2a<T>(a1,ca1,v1,cv1,v2s,cv2s,label+" ColMajor Step2");
    TestMatrixArith2a<T>(a1,ca1,v1s,cv1s,v2s,cv2s,label+" ColMajor Step12");
    TestMatrixArith2a<T>(a2,ca2,v2s,cv2s,v1,cv1,label+" RowMajor Step1");
    TestMatrixArith2a<T>(a2,ca2,v2,cv2,v1s,cv1s,label+" RowMajor Step2");
    TestMatrixArith2a<T>(a2,ca2,v2s,cv2s,v1s,cv1s,label+" RowMajor Step12");
    TestMatrixArith2a<T>(a3,ca3,v1,cv1,v2,cv2,label+" NonMajor");
    TestMatrixArith2a<T>(a3,ca3,v1s,cv1s,v2,cv2,label+" NonMajor Step1");
    TestMatrixArith2a<T>(a3,ca3,v1,cv1,v2s,cv2s,label+" NonMajor Step2");
    TestMatrixArith2a<T>(a3,ca3,v1s,cv1s,v2s,cv2s,label+" NonMajor Step12");

    TestMatrixArith2b<T>(a1,ca1,v1s,cv1s,v2,cv2,label+" ColMajor Step1");
    TestMatrixArith2b<T>(a1,ca1,v1,cv1,v2s,cv2s,label+" ColMajor Step2");
    TestMatrixArith2b<T>(a1,ca1,v1s,cv1s,v2s,cv2s,label+" ColMajor Step12");
    TestMatrixArith2b<T>(a2,ca2,v2s,cv2s,v1,cv1,label+" RowMajor Step1");
    TestMatrixArith2b<T>(a2,ca2,v2,cv2,v1s,cv1s,label+" RowMajor Step2");
    TestMatrixArith2b<T>(a2,ca2,v2s,cv2s,v1s,cv1s,label+" RowMajor Step12");
    TestMatrixArith2b<T>(a3,ca3,v1,cv1,v2,cv2,label+" NonMajor");
    TestMatrixArith2b<T>(a3,ca3,v1s,cv1s,v2,cv2,label+" NonMajor Step1");
    TestMatrixArith2b<T>(a3,ca3,v1,cv1,v2s,cv2s,label+" NonMajor Step2");
    TestMatrixArith2b<T>(a3,ca3,v1s,cv1s,v2s,cv2s,label+" NonMajor Step12");
#endif
#if (XTEST & 2)
    tmv::MatrixView<T> a1m = a1;
    tmv::MatrixView<T> a2m = a2;
    tmv::MatrixView<CT> ca1m = ca1;
    tmv::MatrixView<CT> ca2m = ca2;
    tmv::VectorView<T> v1v = v1;
    tmv::VectorView<T> v2v = v2;
    tmv::VectorView<CT> cv1v = cv1;
    tmv::VectorView<CT> cv2v = cv2;

    TestMatrixArith2a<T>(a1m,ca1m,v1,cv1,v2,cv2,label+" ColMajor");
    TestMatrixArith2a<T>(a2m,ca2m,v2,cv2,v1,cv1,label+" RowMajor");
    TestMatrixArith2a<T>(a1,ca1,v1v,cv1v,v2v,cv2v,label+" ColMajor");
    TestMatrixArith2a<T>(a2,ca2,v2v,cv2v,v1v,cv1v,label+" RowMajor");

    TestMatrixArith2b<T>(a1m,ca1m,v1,cv1,v2,cv2,label+" ColMajor");
    TestMatrixArith2b<T>(a2m,ca2m,v2,cv2,v1,cv1,label+" RowMajor");
    TestMatrixArith2b<T>(a1,ca1,v1v,cv1v,v2v,cv2v,label+" ColMajor");
    TestMatrixArith2b<T>(a2,ca2,v2v,cv2v,v1v,cv1v,label+" RowMajor");
#endif
}
