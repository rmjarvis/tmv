#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV_TestMatrixArith.h"

template <class T, int N> 
static void DoTestSmallMatrixArith_B3d()
{
    tmv::SmallMatrix<T,N,N,tmv::RowMajor> a1;
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) {
        a1(i,j) = T(3+4*i-5*j);
    }
    a1(0,0) = 14;
    if (N > 1) a1(1,0) = -2;
    if (N > 2) a1(2,0) = 7;
    if (N > 3) a1(3,0) = -10;
    if (N > 2) a1(2,2) = 30;

    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor> ca1 = 
        std::complex<T>(3,2)*a1;
    if (N > 3) ca1(2,3) += std::complex<T>(2,4);
    if (N > 1) ca1(1,0) *= std::complex<T>(1,3);
    if (N > 1) ca1.col(1) *= std::complex<T>(-1,4);
    if (N > 3) ca1.row(3) += 
        tmv::SmallVector<std::complex<T>,N>(std::complex<T>(2,9));

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

    tmv::SmallVector<T,N> v1 = a1.row(0);
    tmv::SmallVector<std::complex<T>,N> cv1 = ca1.row(0);
    tmv::SmallVector<T,7> v2 = a3.col(0);
    tmv::SmallVector<std::complex<T>,7> cv2 = ca3.col(0);

    if (showstartdone) {
        std::cout<<"B3d"<<std::endl;
    }
    TestMatrixArith3b(a5,ca5,v2,cv2,v1,cv1,"NonSquare");

#if (XTEST & 2)
    tmv::SmallMatrix<T,N,7,tmv::ColMajor> a6 = a5;
    tmv::SmallMatrix<std::complex<T>,N,7,tmv::ColMajor> ca6 = ca5;
#endif

#if (XTEST & 32)
    tmv::SmallMatrix<T,N,7,tmv::RowMajor|tmv::FortranStyle> a5f = a5;
    tmv::SmallMatrix<std::complex<T>,N,7,tmv::RowMajor|tmv::FortranStyle> ca5f = ca5;
    tmv::SmallMatrix<T,N,7,tmv::ColMajor|tmv::FortranStyle> a6f = a6;
    tmv::SmallMatrix<std::complex<T>,N,7,tmv::ColMajor|tmv::FortranStyle> ca6f = ca6;

    tmv::SmallVector<T,N,tmv::FortranStyle> v1f = v1;
    tmv::SmallVector<std::complex<T>,N,tmv::FortranStyle> cv1f = cv1;
    tmv::SmallVector<T,7,tmv::FortranStyle> v2f = v2;
    tmv::SmallVector<std::complex<T>,7,tmv::FortranStyle> cv2f = cv2;

    TestMatrixArith3b(a5f,ca5f,v2,cv2,v1,cv1,"NonSquare");
    TestMatrixArith3b(a6f,ca6f,v2,cv2,v1,cv1,"NonSquare");
    TestMatrixArith3b(a5f,ca5f,v2f,cv2f,v1,cv1,"NonSquare");
    TestMatrixArith3b(a5f,ca5f,v2f,cv2f,v1f,cv1f,"NonSquare");
#endif
}

template <class T> 
void TestSmallMatrixArith_B3d()
{
    DoTestSmallMatrixArith_B3d<T,2>();
    DoTestSmallMatrixArith_B3d<T,5>();
#if (XTEST & 2)
    DoTestSmallMatrixArith_B3d<T,1>();
    DoTestSmallMatrixArith_B3d<T,3>();
    DoTestSmallMatrixArith_B3d<T,4>();
#endif
}


#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_B3d<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_B3d<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_B3d<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_B3d<int>();
#endif
