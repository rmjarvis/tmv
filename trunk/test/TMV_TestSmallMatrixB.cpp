
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV.h"
#include <fstream>

template <class T> void TestAllSmallMatrixB()
{
    TestSmallMatrixArith_4<T>();
    TestSmallMatrixArith_5<T>();
    TestSmallMatrixArith_6<T>();
    std::cout<<"SmallMatrix<"<<Text(T())<<"> Arithmetic passed all MM tests\n";
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
