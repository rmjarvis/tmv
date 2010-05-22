#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include <fstream>

template <class T> void TestAllSmallMatrixA()
{
    TestSmallMatrixArith_1<T>();
    TestSmallMatrixArith_2<T>();
    TestSmallMatrixArith_3<T>();
    TestSmallMatrixArith_7<T>();
    std::cout<<"SmallMatrix<"<<tmv::TMV_Text(T())<<
        "> Arithmetic passed all MV tests\n";
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
