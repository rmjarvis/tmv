#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV_TestMatrixArith.h"

template <class T, int N> 
static void DoTestSmallMatrixArith_B3a()
{
    std::cout<<"Start B3a "<<std::endl;
    tmv::SmallMatrix<T,N,N,tmv::RowMajor> a1;
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) {
        a1(i,j) = T(3+4*i-5*j);
    }
    a1(0,0) = 14;
    if (N > 1) a1(1,0) = -2;
    if (N > 2) a1(2,0) = 7;
    if (N > 3) a1(3,0) = -10;
    if (N > 2) a1(2,2) = 30;
    std::cout<<"a1 = "<<a1<<std::endl;

    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor> ca1 = 
        std::complex<T>(3,2)*a1;
    std::cout<<"ca1 = "<<ca1<<std::endl;
    if (N > 3) ca1(2,3) += std::complex<T>(2,4);
    std::cout<<"ca1 => "<<ca1<<std::endl;
    if (N > 1) ca1(1,0) *= std::complex<T>(1,3);
    std::cout<<"ca1 => "<<ca1<<std::endl;
    if (N > 1) ca1.col(1) *= std::complex<T>(-1,4);
    std::cout<<"ca1 => "<<ca1<<std::endl;
    if (N > 3) ca1.row(3) += 
        tmv::SmallVector<std::complex<T>,N>(std::complex<T>(2,9));
    std::cout<<"ca1 => "<<ca1<<std::endl;

    tmv::SmallMatrix<T,7,N,tmv::RowMajor> a3;
    for(int i=0;i<7;++i) for(int j=0;j<N;++j) a3(i,j) = T(1-3*i+2*j);
    a3.subMatrix(2,N+2,0,N) += a1;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::RowMajor> ca3 = 
        a3*std::complex<T>(1,2);
    ca3.subMatrix(2,N+2,0,N) += ca1;
    if (N > 1) ca3.col(1) *= std::complex<T>(2,1);
    ca3.row(6).addToAll(std::complex<T>(-7,2));
    std::cout<<"a3 = "<<a3<<std::endl;
    std::cout<<"ca3 = "<<ca3<<std::endl;

    tmv::SmallVector<T,N> v1 = a1.row(0);
    tmv::SmallVector<std::complex<T>,N> cv1 = ca1.row(0);
    tmv::SmallVector<T,7> v2 = a3.col(0);
    tmv::SmallVector<std::complex<T>,7> cv2 = ca3.col(0);
    std::cout<<"v1 = "<<v1<<std::endl;
    std::cout<<"cv1 = "<<cv1<<std::endl;
    std::cout<<"v2 = "<<v2<<std::endl;
    std::cout<<"cv2 = "<<cv2<<std::endl;

#if 0
    if (showstartdone) {
        std::cout<<"B3a"<<std::endl;
    }
    TestMatrixArith3a<T>(a3,ca3,v1,cv1,v2,cv2,"NonSquare");
    TestMatrixArith3b<T>(a3,ca3,v1,cv1,v2,cv2,"NonSquare");

#if (XTEST & 2)
    tmv::SmallMatrix<T,7,N,tmv::ColMajor> a4 = a3;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::ColMajor> ca4 = ca3;
    TestMatrixArith3a<T>(a4,ca4,v1,cv1,v2,cv2,"NonSquare");
    TestMatrixArith3b<T>(a4,ca4,v1,cv1,v2,cv2,"NonSquare");
#endif

#if (XTEST & 32)
    tmv::SmallMatrix<T,7,N,tmv::RowMajor,tmv::FortranStyle> a3f = a3;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::RowMajor,tmv::FortranStyle> ca3f = ca3;
    tmv::SmallMatrix<T,7,N,tmv::ColMajor,tmv::FortranStyle> a4f = a4;
    tmv::SmallMatrix<std::complex<T>,7,N,tmv::ColMajor,tmv::FortranStyle> ca4f = ca4;

    tmv::SmallVector<T,N,tmv::FortranStyle> v1f = v1;
    tmv::SmallVector<std::complex<T>,N,tmv::FortranStyle> cv1f = cv1;
    tmv::SmallVector<T,7,tmv::FortranStyle> v2f = v2;
    tmv::SmallVector<std::complex<T>,7,tmv::FortranStyle> cv2f = cv2;

    TestMatrixArith3a<T>(a3f,ca3f,v1,cv1,v2,cv2,"NonSquare");
    TestMatrixArith3a<T>(a4f,ca4f,v1,cv1,v2,cv2,"NonSquare");
    TestMatrixArith3a<T>(a3f,ca3f,v1f,cv1f,v2,cv2,"NonSquare");
    TestMatrixArith3a<T>(a3f,ca3f,v1f,cv1f,v2f,cv2f,"NonSquare");

    TestMatrixArith3b<T>(a3f,ca3f,v1,cv1,v2,cv2,"NonSquare");
    TestMatrixArith3b<T>(a4f,ca4f,v1,cv1,v2,cv2,"NonSquare");
    TestMatrixArith3b<T>(a3f,ca3f,v1f,cv1f,v2,cv2,"NonSquare");
    TestMatrixArith3b<T>(a3f,ca3f,v1f,cv1f,v2f,cv2f,"NonSquare");
#endif
#endif
}

template <class T> 
void TestSmallMatrixArith_B3a()
{
    DoTestSmallMatrixArith_B3a<T,2>();
    DoTestSmallMatrixArith_B3a<T,5>();
#if (XTEST & 2)
    DoTestSmallMatrixArith_B3a<T,1>();
    DoTestSmallMatrixArith_B3a<T,3>();
    DoTestSmallMatrixArith_B3a<T,4>();
#endif
}


#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_B3a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_B3a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_B3a<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_B3a<int>();
#endif
