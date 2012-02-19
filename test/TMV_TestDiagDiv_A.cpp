#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestDiagDiv_A()
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

    tmv::Matrix<T> r(N,0,T(1));
    tmv::Matrix<std::complex<T> > cr(N,0,std::complex<T>(1));

    tmv::DiagMatrixView<T> av = a.view();
    tmv::DiagMatrixView<T> bv = b.view();
    tmv::DiagMatrixView<std::complex<T> > cav = ca.view();
    tmv::DiagMatrixView<std::complex<T> > cbv = cb.view();

    TestMatrixDivArith2(tmv::LU,av,bv,cav,cbv,"Diag/Diag");
}

#ifdef TEST_DOUBLE
template void TestDiagDiv_A<double>();
#endif
#ifdef TEST_FLOAT
template void TestDiagDiv_A<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestDiagDiv_A<long double>();
#endif
