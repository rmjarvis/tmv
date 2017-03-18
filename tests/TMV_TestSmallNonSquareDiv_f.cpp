
#include "TMV_TestSmallNonSquareDiv.h"

template <class T> 
void TestSmallNonSquareDiv_f()
{
    TestSmallNonSquareDiv<T,tmv::RowMajor,35,12>("RM 35 12");
#if XTEST & 2
    TestSmallNonSquareDiv<T,tmv::ColMajor,65,12>("CM 65 12");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallNonSquareDiv_f<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallNonSquareDiv_f<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallNonSquareDiv_f<long double>();
#endif
