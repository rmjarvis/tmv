
#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_e()
{
    TestSmallSquareDiv<T,tmv::ColMajor,5>("CM 5 5");
#if XTEST & 2
    TestSmallSquareDiv<T,tmv::RowMajor,5>("RM 5 5");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallSquareDiv_e<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallSquareDiv_e<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallSquareDiv_e<long double>();
#endif
