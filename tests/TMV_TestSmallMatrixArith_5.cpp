
#include "TMV_TestSmallMatrixArith_5.h"

template <class T> void TestSmallMatrixArith_5()
{
    TestSmallMatrixArith_5a<T>();
    TestSmallMatrixArith_5b<T>();
    TestSmallMatrixArith_5c<T>();
    TestSmallMatrixArith_5d<T>();
    TestSmallMatrixArith_5e<T>();
    TestSmallMatrixArith_5f<T>();
    TestSmallMatrixArith_5g<T>();
    TestSmallMatrixArith_5h<T>();
    TestSmallMatrixArith_5i<T>();
    TestSmallMatrixArith_5j<T>();
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_5<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_5<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_5<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_5<int>();
#endif
