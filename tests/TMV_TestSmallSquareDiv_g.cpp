
#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_g()
{
    TestSmallSquareDiv<T,tmv::ColMajor,39>("CM 39 39");
#if XTEST & 2
    TestSmallSquareDiv<T,tmv::RowMajor,328>("RM 328 328");
#endif
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
