// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV_Mat.h"
#include <fstream>

template <class T> void TestAllSmallMatrixB()
{
  TestSmallMatrixArith_B1<T>();
  TestSmallMatrixArith_B2a<T>();
  TestSmallMatrixArith_B2b<T>();
  TestSmallMatrixArith_B3a<T>();
  TestSmallMatrixArith_B3b<T>();
  TestSmallMatrixArith_B4a<T>();
  TestSmallMatrixArith_B4b<T>();
  TestSmallMatrixArith_B4c<T>();
  TestSmallMatrixArith_B4d<T>();
  TestSmallMatrixArith_B5a<T>();
  TestSmallMatrixArith_B5b<T>();
  TestSmallMatrixArith_B5c<T>();
  TestSmallMatrixArith_B5d<T>();
  TestSmallMatrixArith_B6a<T>();
  TestSmallMatrixArith_B6b<T>();
  TestSmallMatrixArith_B6c<T>();
  TestSmallMatrixArith_B6d<T>();
  TestSmallMatrixArith_B7a<T>();
  TestSmallMatrixArith_B7b<T>();
  std::cout<<"SmallMatrix<"<<tmv::TypeText(T())<<"> NonSquare Arithmetic passed all tests\n";
}

#ifdef TEST_DOUBLE
template void TestAllSmallMatrixB<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllSmallMatrixB<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllSmallMatrixB<long double>();
#endif
#ifdef TEST_INT
template void TestAllSmallMatrixB<int>();
#endif
