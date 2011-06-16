
#include "TMV_TestSmallNonSquareDiv.h"

template <class T> 
void TestSmallNonSquareDiv_e()
{
    TestSmallNonSquareDiv<T,tmv::ColMajor,10,5>("CM 10 5");
#if XTEST & 2
    TestSmallNonSquareDiv<T,tmv::RowMajor,20,5>("RM 20 5");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallNonSquareDiv_e<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallNonSquareDiv_e<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallNonSquareDiv_e<long double>();
#endif
