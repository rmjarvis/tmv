
#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_f()
{
    TestSmallSquareDiv<T,tmv::RowMajor,12>("RM 12 12");
#if XTEST & 2
    TestSmallSquareDiv<T,tmv::ColMajor,12>("CM 12 12");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallSquareDiv_f<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallSquareDiv_f<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallSquareDiv_f<long double>();
#endif
