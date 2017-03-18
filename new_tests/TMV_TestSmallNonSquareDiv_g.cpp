
#include "TMV_TestSmallNonSquareDiv.h"

template <class T> 
void TestSmallNonSquareDiv_g()
{
    TestSmallNonSquareDiv<T,tmv::ColMajor,41,39>("CM 41 39");
#if XTEST & 2
    TestSmallNonSquareDiv<T,tmv::RowMajor,328,229>("RM 328 229");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallNonSquareDiv_g<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallNonSquareDiv_g<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallNonSquareDiv_g<long double>();
#endif
