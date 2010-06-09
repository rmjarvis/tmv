#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV_TestMatrixArith.h"

template <class T, int N> 
static void DoTestSmallMatrixArith_A3a()
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

    tmv::SmallVector<T,N> v1 = a1.row(0);
    tmv::SmallVector<std::complex<T>,N> cv1 = ca1.row(0);
    tmv::SmallVector<T,N> v2 = v1.col(0);
    tmv::SmallVector<std::complex<T>,N> cv2 = cv1.col(0);

    if (showstartdone) {
        std::cout<<"A3a\n";
    }
    TestMatrixArith3a<T>(a1,ca1,v1,cv1,v2,cv2,"Square");
    TestMatrixArith3b<T>(a1,ca1,v1,cv1,v2,cv2,"Square");

#ifdef XTEST
    tmv::SmallMatrix<T,N,N,tmv::RowMajor,tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor,tmv::FortranStyle> ca1f = ca1;

    tmv::SmallVector<T,N,tmv::FortranStyle> v1f = v1;
    tmv::SmallVector<std::complex<T>,N,tmv::FortranStyle> cv1f = cv1;
    tmv::SmallVector<T,N,tmv::FortranStyle> v2f = v1;
    tmv::SmallVector<std::complex<T>,N,tmv::FortranStyle> cv2f = cv1;

    TestMatrixArith3a<T>(a1f,ca1f,v1,cv1,v2,cv2,"Square");
    TestMatrixArith3a<T>(a1f,ca1f,v1f,cv1f,v2,cv2,"Square");
    TestMatrixArith3a<T>(a1f,ca1f,v1f,cv1f,v2f,cv2f,"Square");
    TestMatrixArith3b<T>(a1f,ca1f,v1,cv1,v2,cv2,"Square");
    TestMatrixArith3b<T>(a1f,ca1f,v1f,cv1f,v2,cv2,"Square");
    TestMatrixArith3b<T>(a1f,ca1f,v1f,cv1f,v2f,cv2f,"Square");
#endif
}

template <class T> 
void TestSmallMatrixArith_A3a()
{
    DoTestSmallMatrixArith_A3a<T,2>();
    DoTestSmallMatrixArith_A3a<T,5>();
#ifdef XTEST
    DoTestSmallMatrixArith_A3a<T,1>();
    DoTestSmallMatrixArith_A3a<T,3>();
    DoTestSmallMatrixArith_A3a<T,4>();
#endif
}


#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_A3a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_A3a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_A3a<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_A3a<int>();
#endif
