
#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_b()
{
    TestSmallSquareDiv<T,tmv::RowMajor,2>("RM 2 2");
#if XTEST & 2
    TestSmallSquareDiv<T,tmv::ColMajor,2>("CM 2 2");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallSquareDiv_b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallSquareDiv_b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallSquareDiv_b<long double>();
#endif
