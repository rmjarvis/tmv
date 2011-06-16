
#include "TMV_TestSmallNonSquareDiv.h"

template <class T> 
void TestSmallNonSquareDiv_a()
{
    TestSmallNonSquareDiv<T,tmv::ColMajor,7,1>("CM 7 1");
#if XTEST & 2
    TestSmallNonSquareDiv<T,tmv::RowMajor,6,1>("RM 6 1");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallNonSquareDiv_a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallNonSquareDiv_a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallNonSquareDiv_a<long double>();
#endif
