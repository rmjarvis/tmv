#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestTriDiv_B1() 
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

    tmv::Matrix<T> mx(m);
    tmv::Matrix<std::complex<T> > cmx(cm);
    m.divideUsing(tmv::LU);
    m.saveDiv();
    m.setDiv();

    TestMatrixDivArith1<T>(
        tmv::LU,mx,cmx,a1.view(),m.view(),ca1.view(),cm.transpose(),"M/U");
    TestMatrixDivArith1<T>(
        tmv::LU,mx,cmx,a1.transpose(),m.view(),ca1.transpose(),cm.view(),
        "M/L");

#if (XTEST & 2)
    tmv::UpperTriMatrix<T,tmv::UnitDiag> a2(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2(cm);
    TestMatrixDivArith1<T>(
        tmv::LU,mx,cmx,a2.view(),m.view(),ca2.view(),cm.transpose(),"M/U");
    TestMatrixDivArith1<T>(
        tmv::LU,mx,cmx,a2.transpose(),m.view(),ca2.transpose(),cm.view(),
        "M/L");
#endif
}

#ifdef TEST_DOUBLE
template void TestTriDiv_B1<double>();
#endif
#ifdef TEST_FLOAT
template void TestTriDiv_B1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestTriDiv_B1<long double>();
#endif
