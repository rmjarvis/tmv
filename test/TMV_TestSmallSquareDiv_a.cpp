
#include "TMV_TestSmallSquareDiv.h"

template <class T> 
void TestSmallSquareDiv_a()
{
    TestSmallSquareDiv<T,tmv::ColMajor,1>();
    TestSmallSquareDiv<T,tmv::RowMajor,1>();
}

#ifdef TEST_DOUBLE
template void TestSmallSquareDiv_a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallSquareDiv_a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallSquareDiv_a<long double>();
#endif
