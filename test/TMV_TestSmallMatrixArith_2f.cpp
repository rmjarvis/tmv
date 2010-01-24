
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_2.h"

template <class T> void TestSmallMatrixArith_2f()
{
    TestSmallMatrixArith_2<T,339,607>("339 607");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_2f<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_2f<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_2f<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_2f<int>();
#endif
