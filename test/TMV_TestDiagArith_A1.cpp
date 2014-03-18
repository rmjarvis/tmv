
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"
#include <fstream>

#define NOSV
#include "TMV_TestMatrixArith.h"

template <class T> void TestDiagMatrixArith_A1()
{
    const int N = 6;

    tmv::DiagMatrix<T> a(N);
    tmv::DiagMatrix<T> b(N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        if (i == j) {
            a(i,j) = T(-7+2*i);
            b(i,j) = T(-5+3*i);
        }

    tmv::DiagMatrix<std::complex<T> > ca = a*std::complex<T>(1,2);
    tmv::DiagMatrix<std::complex<T> > cb = b*std::complex<T>(-5,-1);

    tmv::DiagMatrixView<T> av = a.view();
    tmv::DiagMatrixView<std::complex<T> > cav = ca.view();

    TestMatrixArith1(av,cav, "Diag");
}

#ifdef TEST_DOUBLE
template void TestDiagMatrixArith_A1<double>();
#endif
#ifdef TEST_FLOAT
template void TestDiagMatrixArith_A1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestDiagMatrixArith_A1<long double>();
#endif
#ifdef TEST_INT
template void TestDiagMatrixArith_A1<int>();
#endif
