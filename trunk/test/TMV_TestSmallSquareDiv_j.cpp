
#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_j()
{
#if XTEST & 2
    TestSmallSquareDiv<T,tmv::RowMajor,373>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallSquareDiv_j<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallSquareDiv_j<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallSquareDiv_j<long double>();
#endif
