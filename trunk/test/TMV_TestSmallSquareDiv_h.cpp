
#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_h()
{
#if XTEST & 2
    TestSmallSquareDiv<T,tmv::RowMajor,1>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallSquareDiv_h<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallSquareDiv_h<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallSquareDiv_h<long double>();
#endif
