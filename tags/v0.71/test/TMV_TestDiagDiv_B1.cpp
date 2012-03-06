
#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestDiagDiv_B1()
{
    const int N = 10;

    tmv::DiagMatrix<T> a(N);
    tmv::DiagMatrix<T> b(N);
    for (int i=0; i<N; ++i) {
        a(i,i) = T(3+5*i);
        b(i,i) = T(5-2*i);
    }

    tmv::DiagMatrix<std::complex<T> > ca = a*std::complex<T>(1,2);
    tmv::DiagMatrix<std::complex<T> > cb = b*std::complex<T>(-5,-1);

    tmv::Matrix<T> p(N,N,5);
    p.diag().addToAll(N*10);
    tmv::Matrix<std::complex<T> > cp = p*std::complex<T>(2,3);
    cp += ca;

    tmv::MatrixView<T> pv = p.view();
    tmv::DiagMatrixView<T> bv = b.view();
    tmv::MatrixView<std::complex<T> > cpv = cp.view();
    tmv::DiagMatrixView<std::complex<T> > cbv = cb.view();

    TestMatrixDivArith1(tmv::LU,pv,bv,cpv,cbv,"Diag/SquareM");
}

#ifdef TEST_DOUBLE
template void TestDiagDiv_B1<double>();
#endif
#ifdef TEST_FLOAT
template void TestDiagDiv_B1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestDiagDiv_B1<long double>();
#endif
