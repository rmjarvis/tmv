
#include "TMV_TestSmallMatrixArith_3.h"

template <class T> void TestSmallMatrixArith_3()
{
    TestSmallMatrixArith_3a<T>();
    TestSmallMatrixArith_3b<T>();
    TestSmallMatrixArith_3c<T>();
    TestSmallMatrixArith_3d<T>();
    TestSmallMatrixArith_3e<T>();
    TestSmallMatrixArith_3f<T>();
    TestSmallMatrixArith_3g<T>();
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_3<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_3<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_3<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_3<int>();
#endif
