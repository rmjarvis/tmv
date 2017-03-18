
#include "TMV_TestSmallMatrixArith_3.h"

template <class T> void TestSmallMatrixArith_3d()
{
    TestSmallMatrixArith_3<T,4,4>("4 4");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_3d<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_3d<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_3d<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_3d<int>();
#endif
