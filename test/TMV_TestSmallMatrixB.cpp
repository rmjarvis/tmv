#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"
#include <fstream>

template <class T> void TestAllSmallMatrixB()
{
  TestSmallMatrixArith_B1<T>();
  TestSmallMatrixArith_B3a<T>();
  TestSmallMatrixArith_B3b<T>();
  TestSmallMatrixArith_B4a<T>();
  TestSmallMatrixArith_B4b<T>();
  TestSmallMatrixArith_B4c<T>();
  TestSmallMatrixArith_B5a<T>();
  TestSmallMatrixArith_B5b<T>();
  TestSmallMatrixArith_B5c<T>();
  TestSmallMatrixArith_B6a<T>();
  TestSmallMatrixArith_B6b<T>();
  std::cout<<"SmallMatrix<"<<tmv::Type(T())<<"> NonSquare Arithmetic passed all tests\n";
}

#ifdef INST_DOUBLE
template void TestAllSmallMatrixB<double>();
#endif
#ifdef INST_FLOAT
template void TestAllSmallMatrixB<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllSmallMatrixB<long double>();
#endif
#ifdef INST_INT
template void TestAllSmallMatrixB<int>();
#endif
