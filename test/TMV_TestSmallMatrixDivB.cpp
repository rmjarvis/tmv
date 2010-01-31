
#include "TMV_Test.h"
#include "TMV_Test3.h"

template <class T> 
void TestAllSmallMatrixDivB()
{
    TestSmallMatrixDiv_B1a<T>();
    TestSmallMatrixDiv_B1b<T>();
    TestSmallMatrixDiv_B1c<T>();
    TestSmallMatrixDiv_B2a<T>();
    TestSmallMatrixDiv_B2b<T>();
    TestSmallMatrixDiv_B3a<T>();
    TestSmallMatrixDiv_B3b<T>();
    TestSmallMatrixDiv_B4a<T>();
    TestSmallMatrixDiv_B4b<T>();
    TestSmallMatrixDiv_B5a<T>();
    TestSmallMatrixDiv_B5b<T>();

    std::cout<<"NonSquare SmallMatrix<"<<tmv::TMV_Text(T())<<"> Division ";
    std::cout<<" passed all arithemtic tests\n";
}

#ifdef INST_DOUBLE
template void TestAllSmallMatrixDivB<double>();
#endif
#ifdef INST_FLOAT
template void TestAllSmallMatrixDivB<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllSmallMatrixDivB<long double>();
#endif
