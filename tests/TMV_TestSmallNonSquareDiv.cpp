
#include "TMV_Test.h"
#include "TMV_Test_3.h"

template <class T> 
void TestAllSmallNonSquareDiv()
{
    TestSmallNonSquareDiv_a<T>();
    TestSmallNonSquareDiv_b<T>();
    TestSmallNonSquareDiv_c<T>();
    TestSmallNonSquareDiv_d<T>();
    TestSmallNonSquareDiv_e<T>();
    TestSmallNonSquareDiv_f<T>();
    TestSmallNonSquareDiv_g<T>();

    std::cout<<"NonSquare SmallMatrix<"<<Text(T())<<"> Division ";
    std::cout<<"passed all tests\n";
}

#ifdef TEST_DOUBLE
template void TestAllSmallNonSquareDiv<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllSmallNonSquareDiv<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllSmallNonSquareDiv<long double>();
#endif
