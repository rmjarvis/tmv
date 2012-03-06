
#include "TMV_Test.h"
#include "TMV_Test_3.h"

template <class T> 
void TestAllSmallMatrixDivA()
{
    TestSmallMatrixDiv_A1a<T>();
    TestSmallMatrixDiv_A1b<T>();
    TestSmallMatrixDiv_A1c<T>();
    TestSmallMatrixDiv_A2a<T>();
    TestSmallMatrixDiv_A2b<T>();
    TestSmallMatrixDiv_A3a<T>();
    TestSmallMatrixDiv_A3b<T>();
    TestSmallMatrixDiv_A4<T>();
    TestSmallMatrixDiv_A5<T>();

    std::cout<<"Square SmallMatrix<"<<tmv::TMV_Text(T())<<"> Division ";
    std::cout<<" passed all arithemtic tests\n";
}

#ifdef TEST_DOUBLE
template void TestAllSmallMatrixDivA<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllSmallMatrixDivA<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllSmallMatrixDivA<long double>();
#endif
