#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include <fstream>

template <class T> void TestSmallMatrixB()
{
    TestSmallMatrixArith_4<T>();
    TestSmallMatrixArith_5<T>();
    TestSmallMatrixArith_6<T>();
    std::cout<<"SmallMatrix<"<<tmv::TMV_Text(T())<<
        "> NonSquare Arithmetic passed all tests\n";
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixB<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixB<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixB<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixB<int>();
#endif
