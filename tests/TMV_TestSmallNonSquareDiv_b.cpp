
#include "TMV_TestSmallNonSquareDiv.h"

template <class T> 
void TestSmallNonSquareDiv_b()
{
    TestSmallNonSquareDiv<T,tmv::RowMajor,3,2>("RM 3 2");
#if XTEST & 2
    TestSmallNonSquareDiv<T,tmv::ColMajor,3,2>("CM 3 2");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallNonSquareDiv_b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallNonSquareDiv_b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallNonSquareDiv_b<long double>();
#endif
