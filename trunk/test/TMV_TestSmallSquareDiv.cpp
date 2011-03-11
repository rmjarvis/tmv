
#include "TMV_Test.h"
#include "TMV_Test_3.h"

template <class T> 
void TestAllSmallSquareDiv()
{
    TestSmallSquareDiv_a<T>();
    TestSmallSquareDiv_b<T>();
    TestSmallSquareDiv_c<T>();
    TestSmallSquareDiv_d<T>();
    TestSmallSquareDiv_e<T>();
    TestSmallSquareDiv_f<T>();
    TestSmallSquareDiv_g<T>();

    std::cout<<"Square SmallMatrix<"<<tmv::TMV_Text(T())<<"> Division ";
    std::cout<<"passed all tests\n";
}

#ifdef TEST_DOUBLE
template void TestAllSmallSquareDiv<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllSmallSquareDiv<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllSmallSquareDiv<long double>();
#endif
