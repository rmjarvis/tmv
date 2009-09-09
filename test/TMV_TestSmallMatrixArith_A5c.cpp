// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV_Mat.h"

#include "TMV_TestMatrixArith.h"

template <class T, int N> static void DoTestSmallMatrixArith_A5c()
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

  if (N > 10) {
    a1 /= T(N*N); a1 += T(1); 
    ca1 /= T(N*N); ca1 += T(1); 
  }


  tmv::SmallMatrix<T,N,N,tmv::RowMajor> a1b = a1;
  tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor> ca1b = ca1;

  if (showstartdone) {
    std::cout<<"A5c\n";
  }
  TestMatrixArith5<T>(a1,ca1,a1b,ca1b,"Square");

#ifdef XTEST
  tmv::SmallMatrix<T,N,N,tmv::RowMajor,tmv::FortranStyle> a1f = a1;
  tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor,tmv::FortranStyle> ca1f = ca1;
  tmv::SmallMatrix<T,N,N,tmv::ColMajor,tmv::FortranStyle> a2f = a1;
  tmv::SmallMatrix<std::complex<T>,N,N,tmv::ColMajor,tmv::FortranStyle> ca2f = ca1;

  TestMatrixArith5<T>(a1f,ca1f,a2f,ca2f,"Square");
#endif
}

template <class T> void TestSmallMatrixArith_A5c()
{
  DoTestSmallMatrixArith_A5c<T,2>();
  DoTestSmallMatrixArith_A5c<T,12>();
#ifdef XTEST
  DoTestSmallMatrixArith_A5c<T,1>();
  DoTestSmallMatrixArith_A5c<T,3>();
  DoTestSmallMatrixArith_A5c<T,4>();
  DoTestSmallMatrixArith_A5c<T,22>();
#endif
}


#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_A5c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_A5c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_A5c<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_A5c<int>();
#endif
