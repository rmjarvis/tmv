
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"
#include <fstream>

#include "TMV_TestMatrixArith.h"

template <class T> void TestDiagMatrixArith_A5()
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

    tmv::DiagMatrixView<T> bv = b.view();
    tmv::DiagMatrixView<std::complex<T> > cbv = cb.view();

    TestMatrixArith5(av,cav,bv,cbv, "Diag/Diag");
}

#ifdef TEST_DOUBLE
template void TestDiagMatrixArith_A5<double>();
#endif
#ifdef TEST_FLOAT
template void TestDiagMatrixArith_A5<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestDiagMatrixArith_A5<long double>();
#endif
#ifdef TEST_INT
template void TestDiagMatrixArith_A5<int>();
#endif
