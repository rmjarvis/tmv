// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"

#include "TMV_TestMatrixArith.h"

template <class T, int N> static void DoTestSmallMatrixArith_B4c()
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

  tmv::SmallMatrix<T,7,N,tmv::RowMajor> a3;
  for(int i=0;i<7;++i) for(int j=0;j<N;++j) a3(i,j) = T(1-3*i+2*j);
  a3.SubMatrix(2,N+2,0,N) += a1;
  tmv::SmallMatrix<std::complex<T>,7,N,tmv::RowMajor> ca3 = 
  a3*std::complex<T>(1,2);
  ca3.SubMatrix(2,N+2,0,N) += ca1;
  if (N > 1) ca3.col(1) *= std::complex<T>(2,1);
  ca3.row(6).AddToAll(std::complex<T>(-7,2));

  tmv::SmallMatrix<T,N,7,tmv::RowMajor> a5 = a3.Transpose();
  tmv::SmallMatrix<std::complex<T>,N,7,tmv::RowMajor> ca5 = ca3.Adjoint();
  tmv::SmallMatrix<T,N,7,tmv::ColMajor> a6 = a5;
  tmv::SmallMatrix<std::complex<T>,N,7,tmv::ColMajor> ca6 = ca5;

  tmv::SmallMatrix<T,N,7> a5x;
  tmv::SmallMatrix<std::complex<T>,N,7> ca5x;

  if (showstartdone) {
    std::cout<<"B4c\n";
  }
  TestMatrixArith4<T>(a5x,ca5x,a5,ca5,a6,ca6,"NonSquare");

#ifdef XTEST
  tmv::SmallMatrix<T,N,7,tmv::RowMajor> a5b = a5;
  tmv::SmallMatrix<std::complex<T>,N,7,tmv::RowMajor> ca5b = ca5;

  tmv::SmallMatrix<T,N,7,tmv::RowMajor,tmv::FortranStyle> a5f = a5;
  tmv::SmallMatrix<std::complex<T>,N,7,tmv::RowMajor,tmv::FortranStyle> ca5f = ca5;
  tmv::SmallMatrix<T,N,7,tmv::ColMajor,tmv::FortranStyle> a6f = a6;
  tmv::SmallMatrix<std::complex<T>,N,7,tmv::ColMajor,tmv::FortranStyle> ca6f = ca6;

  TestMatrixArith4<T>(a5x,ca5x,a5,ca5,a5b,ca5b,"NonSquare");
  TestMatrixArith4<T>(a5x,ca5x,a5f,ca5f,a6,ca6,"NonSquare");
  TestMatrixArith4<T>(a5x,ca5x,a5f,ca5f,a6f,ca6f,"NonSquare");
#endif
}

template <class T> void TestSmallMatrixArith_B4c()
{
  DoTestSmallMatrixArith_B4c<T,2>();
  DoTestSmallMatrixArith_B4c<T,5>();
#ifdef XTEST
  DoTestSmallMatrixArith_B4c<T,1>();
  DoTestSmallMatrixArith_B4c<T,3>();
  DoTestSmallMatrixArith_B4c<T,4>();
#endif
}

#ifdef INST_DOUBLE
template void TestSmallMatrixArith_B4c<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrixArith_B4c<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrixArith_B4c<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrixArith_B4c<int>();
#endif
