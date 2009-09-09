// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_Test.h"
#include "TMV_Test3.h"

#include "TMV.h"
#include <fstream>

// Remove this once it's ok to test Det()
#define NODIV

#include "TMV_TestMatrixArith.h"
#define CT std::complex<T>

template <class T> void TestSmallMatrixArith_1a();
template <class T> void TestSmallMatrixArith_1b();
template <class T> void TestSmallMatrixArith_1c();
template <class T> void TestSmallMatrixArith_1d();
template <class T> void TestSmallMatrixArith_1e();
template <class T> void TestSmallMatrixArith_1f();

template <class T, int M, int N> void TestSmallMatrixArith_1(std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Arith_1: M,N = "<<M<<','<<N<<std::endl;
  }

  tmv::SmallMatrix<T,M,N,tmv::RowMajor> a1x;
  for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
    a1x(i,j) = T(2.9+4.3*i-5.1*j);
  }
  a1x.diag().AddToAll(T(6));
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

#ifndef NONSQUARE
  // These next two is to make sure Det is calculable without overflow.
  if (N > 10) {
    a1x /= T(N*N); a1x += T(1);
    ca1x /= T(N*N); ca1x += T(1);
  }
#endif

  tmv::SmallMatrix<T,N,M,tmv::ColMajor> a2x = a1x.Transpose();
  if (N > 1) a2x.row(1) *= T(3.1);
  if (M > 2) a2x.col(2) -= tmv::SmallVector<T,N>(4.9);
  tmv::SmallMatrix<CT,N,M,tmv::ColMajor> ca2x = ca1x.Transpose();
  ca2x -= T(1.3)*a2x;
  ca2x *= CT(1.1,-2.5);

  tmv::SmallMatrixView<T,M,N,1,M> a1 = a1x.View();
  tmv::SmallMatrixView<T,N,M,M,1> a2 = a2x.View();
  tmv::SmallMatrixView<CT,M,N,1,M> ca1 = ca1x.View();
  tmv::SmallMatrixView<CT,N,M,M,1> ca2 = ca2x.View();

  TestMatrixArith1<T>(a1,ca1,label+" ColMajor");
  TestMatrixArith1<T>(a2,ca2,label+" RowMajor");

#ifdef XTEST
  tmv::SmallMatrix<T,3*M,4*N,tmv::ColMajor> a3x;
  tmv::SmallMatrixView<T,M,N,3,12*M> a3 = a3x.SubMatrix(0,3*M,0,4*N,3,4);
  a3 = a1;
  tmv::SmallMatrix<CT,3*M,4*N,tmv::ColMajor> ca3x;
  tmv::SmallMatrixView<CT,M,N,3,12*M> ca3 = ca3x.SubMatrix(0,3*M,0,4*N,3,4);
  ca3 = ca1;

  TestMatrixArith1<T>(a3,ca3,label+" NonMajor");
#endif
}
