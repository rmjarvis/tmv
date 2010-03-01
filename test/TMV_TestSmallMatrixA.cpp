#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include <fstream>

template <class T> void TestSmallMatrixA()
{
    TestSmallMatrixArith_1<T>();
    TestSmallMatrixArith_2<T>();
    TestSmallMatrixArith_3<T>();
    TestSmallMatrixArith_7<T>();
    std::cout<<"SmallMatrix<"<<tmv::TMV_Text(T())<<
        "> Square Arithmetic passed all tests\n";
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixA<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixA<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixA<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixA<int>();
#endif
