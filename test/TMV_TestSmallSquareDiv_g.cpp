
#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_g()
{
    TestSmallSquareDiv<T,tmv::ColMajor,39>();
}

#ifdef TEST_DOUBLE
template void TestSmallSquareDiv_g<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallSquareDiv_g<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallSquareDiv_g<long double>();
#endif
