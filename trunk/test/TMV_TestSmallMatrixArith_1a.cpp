
//#define PRINTALGO_NormM
//#define PRINTALGO_SVD
//#define XDEBUG_SVD
//#undef NDEBUG

#include "TMV_TestSmallMatrixArith_1.h"

template <class T> void TestSmallMatrixArith_1a()
{
    TestSmallMatrixArith_1<T,1,1>("1 1");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_1a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_1a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_1a<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_1a<int>();
#endif
