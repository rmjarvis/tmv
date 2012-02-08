#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include <fstream>

template <class T> void TestAllSmallMatrixB()
{
    TestSmallMatrixArith_B1<T>();
    TestSmallMatrixArith_B2a<T>();
    TestSmallMatrixArith_B2b<T>();
    TestSmallMatrixArith_B3a<T>();
    TestSmallMatrixArith_B3b<T>();
    TestSmallMatrixArith_B3c<T>();
    TestSmallMatrixArith_B3d<T>();
    TestSmallMatrixArith_B4a<T>();
    TestSmallMatrixArith_B4b<T>();
    TestSmallMatrixArith_B4c<T>();
    TestSmallMatrixArith_B4d<T>();
    TestSmallMatrixArith_B5a<T>();
    TestSmallMatrixArith_B5b<T>();
    TestSmallMatrixArith_B5c<T>();
    TestSmallMatrixArith_B5d<T>();
    TestSmallMatrixArith_B6a<T>();
    TestSmallMatrixArith_B6b<T>();
    TestSmallMatrixArith_B6c<T>();
    TestSmallMatrixArith_B6d<T>();
    TestSmallMatrixArith_B7a<T>();
    TestSmallMatrixArith_B7b<T>();
    std::cout<<"SmallMatrix<"<<tmv::TMV_Text(T())<<
        "> NonSquare Arithmetic passed all tests\n";
}

#ifdef TEST_DOUBLE
template void TestAllSmallMatrixB<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllSmallMatrixB<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllSmallMatrixB<long double>();
#endif
#ifdef TEST_INT
template void TestAllSmallMatrixB<int>();
#endif
