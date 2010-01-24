#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include <fstream>

template <class T> void TestAllSmallMatrixA()
{
    TestSmallMatrixArith_A1a<T>();
    TestSmallMatrixArith_A1b<T>();
    TestSmallMatrixArith_A2a<T>();
    TestSmallMatrixArith_A2b<T>();
    TestSmallMatrixArith_A2c<T>();
    TestSmallMatrixArith_A3a<T>();
    TestSmallMatrixArith_A3b<T>();
    TestSmallMatrixArith_A3c<T>();
    TestSmallMatrixArith_A3d<T>();
    TestSmallMatrixArith_A3e<T>();
    TestSmallMatrixArith_A4a<T>();
    TestSmallMatrixArith_A4b<T>();
    TestSmallMatrixArith_A4c<T>();
    TestSmallMatrixArith_A4d<T>();
    TestSmallMatrixArith_A4e<T>();
    TestSmallMatrixArith_A5a<T>();
    TestSmallMatrixArith_A5b<T>();
    TestSmallMatrixArith_A5c<T>();
    TestSmallMatrixArith_A5d<T>();
    TestSmallMatrixArith_A5e<T>();
    TestSmallMatrixArith_A6a<T>();
    TestSmallMatrixArith_A6b<T>();
    TestSmallMatrixArith_A6c<T>();
    TestSmallMatrixArith_A6d<T>();
    TestSmallMatrixArith_A6e<T>();
    TestSmallMatrixArith_A7a<T>();
    TestSmallMatrixArith_A7b<T>();
    TestSmallMatrixArith_A7c<T>();
    std::cout<<"SmallMatrix<"<<tmv::TMV_Text(T())<<"> Square Arithmetic passed all tests\n";
}

#ifdef TEST_DOUBLE
template void TestAllSmallMatrixA<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllSmallMatrixA<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllSmallMatrixA<long double>();
#endif
#ifdef TEST_INT
template void TestAllSmallMatrixA<int>();
#endif
