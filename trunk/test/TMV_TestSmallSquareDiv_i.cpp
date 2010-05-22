
#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_i()
{
#if XTEST & 2
    TestSmallSquareDiv<T,tmv::ColMajor,189>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallSquareDiv_i<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallSquareDiv_i<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallSquareDiv_i<long double>();
#endif
