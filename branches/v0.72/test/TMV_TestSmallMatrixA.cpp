#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include <fstream>

template <class T> 
void TestAllSmallMatrixA()
{
    TestSmallMatrixArith_A1<T>();
    TestSmallMatrixArith_A2a<T>();
    TestSmallMatrixArith_A2b<T>();
    TestSmallMatrixArith_A3a<T>();
    TestSmallMatrixArith_A3b<T>();
    TestSmallMatrixArith_A3c<T>();
    TestSmallMatrixArith_A3d<T>();
    TestSmallMatrixArith_A4a<T>();
    TestSmallMatrixArith_A4b<T>();
    TestSmallMatrixArith_A5a<T>();
    TestSmallMatrixArith_A5b<T>();
    TestSmallMatrixArith_A6a<T>();
    TestSmallMatrixArith_A6b<T>();
    TestSmallMatrixArith_A7<T>();
    std::cout<<"SmallMatrix<"<<tmv::TMV_Text(T())<<
        "> Square Arithmetic passed all tests\n";
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
