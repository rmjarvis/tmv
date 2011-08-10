
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV.h"
#include <fstream>

// Remove this once it's ok to test Det()
#define NODIV

#include "TMV_TestMatrixArith.h"
#define CT std::complex<T>

template <class T> void TestSmallMatrixArith_7a();
template <class T> void TestSmallMatrixArith_7b();
template <class T> void TestSmallMatrixArith_7c();
template <class T> void TestSmallMatrixArith_7d();
template <class T> void TestSmallMatrixArith_7e();
template <class T> void TestSmallMatrixArith_7f();
template <class T> void TestSmallMatrixArith_7g();

template <class T, int M, int N> void TestSmallMatrixArith_7(std::string label)
{
    if (showstartdone) {
        std::cout<<"Start Arith_7: M,N = "<<M<<','<<N<<std::endl;
    }

    tmv::SmallMatrix<T,M,N,tmv::RowMajor> a1;
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
        a1(i,j) = T(2.9+4.3*i-5.1*j);
    }
    a1.diag(0,0,std::min(M,N)).addToAll(T(6));
    a1(0,0) = 14; 
    if (M > 1) a1(1,0) = -2; 
    if (M > 2) a1(2,0) = 7; 
    if (M > 3) a1(3,0) = -10;
    if (M > 2 && N > 2) a1(2,2) = 30;

    tmv::SmallMatrix<CT,M,N,tmv::RowMajor> ca1 = CT(3,2)*a1;
    if (M > 2 && N > 3) ca1(2,3) += CT(2.4,3.7);
    if (M > 1) ca1(1,0) *= CT(0.8,2.8);
    if (N > 1) ca1.col(1) *= CT(-1.1,3.6);
    if (M > 3) ca1.row(3) += tmv::SmallVector<CT,N>(CT(1.8,9.2));

    tmv::SmallMatrix<T,N,M,tmv::ColMajor> a2 = a1.transpose();
    if (N > 1) a2.row(1) *= T(3.1);
    if (M > 2) a2.col(2) -= tmv::SmallVector<T,N>(4.9);
    tmv::SmallMatrix<CT,N,M,tmv::ColMajor> ca2 = ca1.transpose();
    ca2 -= T(1.3)*a2;
    ca2 *= CT(1.1,-2.5);

    tmv::SmallVector<T,M> v1 = a1.col(0);
    tmv::SmallVector<T,N> v2 = a1.row(0);
    tmv::SmallVector<CT,M> cv1 = ca1.col(0);
    tmv::SmallVector<CT,N> cv2 = ca1.row(0);

    TestMatrixArith7<T>(a1,ca1,v1,cv1,v2,cv2,label+" ColMajor");
    TestMatrixArith7<T>(a2,ca2,v2,cv2,v1,cv1,label+" RowMajor");

#if (XTEST & 1)
    tmv::SmallMatrix<T,3*M,4*N,tmv::ColMajor> a3x;
    tmv::SmallMatrixView<T,M,N> a3 = a3x.subMatrix(0,3*M,0,4*N,3,4);
    a3 = a1;
    tmv::SmallMatrix<CT,3*M,4*N,tmv::ColMajor> ca3x;
    tmv::SmallMatrixView<CT,M,N> ca3 = ca3x.subMatrix(0,3*M,0,4*N,3,4);
    ca3 = ca1;

    tmv::SmallVector<T,5*M> v3x;
    tmv::SmallVectorView<T,M> v3 = v3x.subVector(0,5*M,5);
    v3 = v1;
    tmv::SmallVector<T,5*N> v4x;
    tmv::SmallVectorView<T,N> v4 = v4x.subVector(0,5*N,5);
    v4 = v2;
    tmv::SmallVector<CT,5*M> cv3x;
    tmv::SmallVectorView<CT,M> cv3 = cv3x.subVector(0,5*M,5);
    cv3 = cv1;
    tmv::SmallVector<CT,5*N> cv4x;
    tmv::SmallVectorView<CT,N> cv4 = cv4x.subVector(0,5*N,5);
    cv4 = cv2;

    TestMatrixArith7<T>(a1,ca1,v3,cv3,v2,cv2,label+" ColMajor Step1");
    TestMatrixArith7<T>(a1,ca1,v1,cv1,v4,cv4,label+" ColMajor Step2");
    TestMatrixArith7<T>(a1,ca1,v3,cv3,v4,cv4,label+" ColMajor Step12");
    TestMatrixArith7<T>(a2,ca2,v4,cv4,v1,cv1,label+" RowMajor Step1");
    TestMatrixArith7<T>(a2,ca2,v2,cv2,v3,cv3,label+" RowMajor Step2");
    TestMatrixArith7<T>(a2,ca2,v4,cv4,v3,cv3,label+" RowMajor Step12");
    TestMatrixArith7<T>(a3,ca3,v1,cv1,v2,cv2,label+" NonMajor");
    TestMatrixArith7<T>(a3,ca3,v3,cv3,v2,cv2,label+" NonMajor Step1");
    TestMatrixArith7<T>(a3,ca3,v1,cv1,v4,cv4,label+" NonMajor Step2");
    TestMatrixArith7<T>(a3,ca3,v3,cv3,v4,cv4,label+" NonMajor Step12");
#endif
#if (XTEST & 2)
    tmv::MatrixView<T> a1m = a1.view();
    tmv::MatrixView<CT> ca1m = ca1.view();
    tmv::VectorView<T> v1v = v1.view();
    tmv::VectorView<T> v2v = v2.view();
    tmv::VectorView<CT> cv1v = cv1.view();
    tmv::VectorView<CT> cv2v = cv2.view();

    TestMatrixArith7<T>(a1m,ca1m,v1,cv1,v2,cv2,label+" ColMajor");
    TestMatrixArith7<T>(a1,ca1,v1v,cv1v,v2v,cv2v,label+" ColMajor");
#endif
}
