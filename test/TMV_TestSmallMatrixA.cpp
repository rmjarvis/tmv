#ifdef NDEBUG
#undef NDEBUG
#endif
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"
#include <fstream>

template <class T> void TestAllSmallMatrixA()
{
  TestSmallMatrixArith_A1<T>();
  TestSmallMatrixArith_A2a<T>();
  TestSmallMatrixArith_A2b<T>();
  TestSmallMatrixArith_A2c<T>();
  TestSmallMatrixArith_A3a<T>();
  TestSmallMatrixArith_A3b<T>();
  TestSmallMatrixArith_A3c<T>();
  TestSmallMatrixArith_A4a<T>();
  TestSmallMatrixArith_A4b<T>();
  TestSmallMatrixArith_A4c<T>();
  TestSmallMatrixArith_A5a<T>();
  TestSmallMatrixArith_A5b<T>();
  TestSmallMatrixArith_A5c<T>();
  TestSmallMatrixArith_A6a<T>();
  TestSmallMatrixArith_A6b<T>();
  TestSmallMatrixArith_A6c<T>();
  std::cout<<"SmallMatrix<"<<tmv::Type(T())<<"> Square Arithmetic passed all tests\n";
}

#ifdef INST_DOUBLE
template void TestAllSmallMatrixA<double>();
#endif
#ifdef INST_FLOAT
template void TestAllSmallMatrixA<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllSmallMatrixA<long double>();
#endif
#ifdef INST_INT
template void TestAllSmallMatrixA<int>();
#endif
