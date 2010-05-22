
#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_c()
{
    TestSmallSquareDiv<T,tmv::ColMajor,4>();
}

#ifdef TEST_DOUBLE
template void TestSmallSquareDiv_c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallSquareDiv_c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallSquareDiv_c<long double>();
#endif
