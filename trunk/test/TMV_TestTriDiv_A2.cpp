#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> void TestTriDiv_A2() 
{
    const int N = 10;

    tmv::Matrix<T> m(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        m(i,j) = T(0.4+0.02*i-0.05*j);
    m.diag().addToAll(5);
    m.diag(1).addToAll(T(0.32));
    m.diag(-1).addToAll(T(0.91));

    tmv::Matrix<std::complex<T> > cm(m);
    cm += std::complex<T>(10,2);
    cm.diag(1) *= std::complex<T>(T(-0.5),T(-0.8));
    cm.diag(-1) *= std::complex<T>(T(-0.7),T(0.1));

    tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a1(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca1(cm);
    tmv::UpperTriMatrix<T,tmv::UnitDiag> a2(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2(cm);
    a1.saveDiv();
    a2.saveDiv();
    ca1.saveDiv();
    ca2.saveDiv();
    a1.setDiv();
    a2.setDiv();
    ca1.setDiv();
    ca2.setDiv();

    tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a1x(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca1x(cm);
    tmv::UpperTriMatrix<T,tmv::UnitDiag> a2x(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2x(cm);
    tmv::LowerTriMatrix<T,tmv::NonUnitDiag> b1x(m);
    tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag> cb1x(cm);
    tmv::LowerTriMatrix<T,tmv::UnitDiag> b2x(m);
    tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cb2x(cm);

    TestMatrixDivArith1<T>(
        tmv::LU,b2x,cb2x,a1.view(),a2.transpose(),ca1.view(),ca2.transpose(),
        "L/U");
    TestMatrixDivArith1<T>(
        tmv::LU,a2x,ca2x,a1.transpose(),a2.view(),ca1.transpose(),ca2.view(),
        "U/L");
    TestMatrixDivArith1<T>(
        tmv::LU,b1x,cb1x,a2.view(),a1.transpose(),ca2.view(),ca1.transpose(),
        "L/U");
    TestMatrixDivArith1<T>(
        tmv::LU,a1x,ca1x,a2.transpose(),a1.view(),ca2.transpose(),ca1.view(),
        "U/L");

#ifdef XTEST
#if (XTEST & 2)
    tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a1b(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca1b(cm);
    tmv::UpperTriMatrix<T,tmv::UnitDiag> a2b(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2b(cm);

    TestMatrixDivArith1<T>(
        tmv::LU,b1x,cb1x,a1.view(),a1b.transpose(),ca1.view(),ca1b.transpose(),
        "L/U");
    TestMatrixDivArith1<T>(
        tmv::LU,a1x,ca1x,a1.transpose(),a1b.view(),ca1.transpose(),ca1b.view(),
        "U/L");
    TestMatrixDivArith1<T>(
        tmv::LU,b2x,cb2x,a2.view(),a2b.transpose(),ca2.view(),ca2b.transpose(),
        "L/U");
    TestMatrixDivArith1<T>(
        tmv::LU,a2x,ca2x,a2.transpose(),a2b.view(),ca2.transpose(),ca2b.view(),
        "U/L");
#endif
#endif
}

#ifdef TEST_DOUBLE
template void TestTriDiv_A2<double>();
#endif
#ifdef TEST_FLOAT
template void TestTriDiv_A2<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestTriDiv_A2<long double>();
#endif
