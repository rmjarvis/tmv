#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_b()
{
    TestSmallSquareDiv<T,tmv::RowMajor,3>();
}

#ifdef TEST_DOUBLE
template void TestSmallSquareDiv_b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallSquareDiv_b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallSquareDiv_b<long double>();
#endif