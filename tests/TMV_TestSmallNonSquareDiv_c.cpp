
#include "TMV_TestSmallNonSquareDiv.h"

template <class T> 
void TestSmallNonSquareDiv_c()
{
    TestSmallNonSquareDiv<T,tmv::ColMajor,4,3>("CM 4 3");
#if XTEST & 2
    TestSmallNonSquareDiv<T,tmv::RowMajor,5,3>("RM 5 3");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallNonSquareDiv_c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallNonSquareDiv_c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallNonSquareDiv_c<long double>();
#endif
