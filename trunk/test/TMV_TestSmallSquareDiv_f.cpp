
#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_f()
{
#if XTEST & 2
    TestSmallSquareDiv<T,tmv::RowMajor,2>();
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
