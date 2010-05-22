
#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_d()
{
    TestSmallSquareDiv<T,tmv::RowMajor,5>();
}

#ifdef TEST_DOUBLE
template void TestSmallSquareDiv_d<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallSquareDiv_d<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallSquareDiv_d<long double>();
#endif
