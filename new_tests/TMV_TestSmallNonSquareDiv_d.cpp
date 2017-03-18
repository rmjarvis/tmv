
#include "TMV_TestSmallNonSquareDiv.h"

template <class T> 
void TestSmallNonSquareDiv_d()
{
    TestSmallNonSquareDiv<T,tmv::RowMajor,6,4>("RM 6 4");
#if XTEST & 2
    TestSmallNonSquareDiv<T,tmv::ColMajor,5,4>("CM 5 4");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallNonSquareDiv_d<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallNonSquareDiv_d<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallNonSquareDiv_d<long double>();
#endif
