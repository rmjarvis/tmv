
#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_d()
{
    TestSmallSquareDiv<T,tmv::ColMajor,4>("CM 4 4");
#if XTEST & 2
    TestSmallSquareDiv<T,tmv::RowMajor,4>("RM 4 4");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallSquareDiv_d<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallSquareDiv_d<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallSquareDiv_d<long double>();
#endif
