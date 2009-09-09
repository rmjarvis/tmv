// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV_Mat.h"

#define INORDER
#define NOMULTEQ
#include "TMV_TestMatrixArith.h"

template <class T, int N> static void DoTestSmallMatrixArith_B6c()
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

  tmv::SmallMatrix<T,7,N,tmv::RowMajor> a3;
  for(int i=0;i<7;++i) for(int j=0;j<N;++j) a3(i,j) = T(1-3*i+2*j);
  if (N > 7) a3 += a1.Rows(1,8);
  else a3.SubMatrix(2,N+2,0,N) += a1;
  tmv::SmallMatrix<std::complex<T>,7,N,tmv::RowMajor> ca3 = 
  a3*std::complex<T>(1,2);
  if (N > 7) ca3 += ca1.Rows(1,8);
  else ca3.SubMatrix(2,N+2,0,N) += ca1;
  if (N > 1) ca3.col(1) *= std::complex<T>(2,1);
  ca3.row(6).AddToAll(std::complex<T>(-7,2));

  tmv::SmallMatrix<T,N,7,tmv::RowMajor> a5 = a3.Transpose();
  if (N > 7) a5 -= a1.Cols(0,7);
  else a5.SubMatrix(0,N,1,N+1) -= a1;
  tmv::SmallMatrix<std::complex<T>,N,7,tmv::RowMajor> ca5 = ca3.Adjoint();
  if (N > 7) ca5 -= ca1.Cols(0,7);
  else ca5.SubMatrix(0,N,1,N+1) -= ca1;
  ca5.col(3) *= std::complex<T>(-1,3);
  ca5.row(0).AddToAll(std::complex<T>(1,9));

  tmv::SmallMatrix<T,7,7,tmv::RowMajor> a7;
  tmv::SmallMatrix<std::complex<T>,7,7,tmv::RowMajor> ca7;

  if (showstartdone) {
    std::cout<<"B6c\n";
  }
  TestMatrixArith6<T>(a3,ca3,a5,ca5,a7,ca7,"NonSquare");
#ifdef XTEST
  tmv::SmallMatrix<T,7,N,tmv::ColMajor> a4 = a3;
  tmv::SmallMatrix<std::complex<T>,7,N,tmv::ColMajor> ca4 = ca3;
  tmv::SmallMatrix<T,N,7,tmv::ColMajor> a6 = a5;
  tmv::SmallMatrix<std::complex<T>,N,7,tmv::ColMajor> ca6 = ca5;
  tmv::SmallMatrix<T,7,7,tmv::ColMajor> a8 = a7;
  tmv::SmallMatrix<std::complex<T>,7,7,tmv::ColMajor> ca8 = ca7;

  tmv::SmallMatrix<T,7,N,tmv::RowMajor,tmv::FortranStyle> a3f = a3;
  tmv::SmallMatrix<std::complex<T>,7,N,tmv::RowMajor,tmv::FortranStyle> ca3f = ca3;
  tmv::SmallMatrix<T,N,7,tmv::RowMajor,tmv::FortranStyle> a5f = a5;
  tmv::SmallMatrix<std::complex<T>,N,7,tmv::RowMajor,tmv::FortranStyle> ca5f = ca5;
  tmv::SmallMatrix<T,7,7,tmv::ColMajor,tmv::FortranStyle> a7f = a7;
  tmv::SmallMatrix<std::complex<T>,7,7,tmv::ColMajor,tmv::FortranStyle> ca7f = ca7;

  TestMatrixArith6<T>(a3,ca3,a5,ca5,a8,ca8,"NonSquare");
  TestMatrixArith6<T>(a3,ca3,a6,ca6,a7,ca7,"NonSquare");
  TestMatrixArith6<T>(a3,ca3,a6,ca6,a8,ca8,"NonSquare");
  TestMatrixArith6<T>(a4,ca4,a5,ca5,a7,ca7,"NonSquare");
  TestMatrixArith6<T>(a4,ca4,a5,ca5,a8,ca8,"NonSquare");
  TestMatrixArith6<T>(a4,ca4,a6,ca6,a7,ca7,"NonSquare");
  TestMatrixArith6<T>(a4,ca4,a6,ca6,a8,ca8,"NonSquare");
  TestMatrixArith6<T>(a3f,ca3f,a5,ca5,a7,ca7,"NonSquare");
  TestMatrixArith6<T>(a3f,ca3f,a5f,ca5f,a7,ca7,"NonSquare");
  TestMatrixArith6<T>(a3f,ca3f,a5f,ca5f,a7f,ca7f,"NonSquare");
#endif
}

template <class T> void TestSmallMatrixArith_B6c()
{
  DoTestSmallMatrixArith_B6c<T,2>();
  DoTestSmallMatrixArith_B6c<T,12>();
#ifdef XTEST
  DoTestSmallMatrixArith_B6c<T,1>();
  DoTestSmallMatrixArith_B6c<T,3>();
  DoTestSmallMatrixArith_B6c<T,4>();
  DoTestSmallMatrixArith_B6c<T,22>();
  DoTestSmallMatrixArith_B6c<T,555>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_B6c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_B6c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_B6c<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_B6c<int>();
#endif
