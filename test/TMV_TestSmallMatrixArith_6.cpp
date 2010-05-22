
#include "TMV_TestSmallMatrixArith_6.h"

template <class T> void TestSmallMatrixArith_6()
{
    TestSmallMatrixArith_6a<T>();
    TestSmallMatrixArith_6b<T>();
    TestSmallMatrixArith_6c<T>();
    TestSmallMatrixArith_6d<T>();
    TestSmallMatrixArith_6e<T>();
    TestSmallMatrixArith_6f<T>();
    TestSmallMatrixArith_6g<T>();
    TestSmallMatrixArith_6h<T>();
    TestSmallMatrixArith_6i<T>();
    TestSmallMatrixArith_6j<T>();
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_6<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_6<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_6<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_6<int>();
#endif
