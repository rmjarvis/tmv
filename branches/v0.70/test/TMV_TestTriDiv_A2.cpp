#include "TMV_Test.h"
#include "TMV_Test_1.h"
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

    tmv::UpperTriMatrixView<T> a1v = a1.view();
    tmv::UpperTriMatrixView<T> a2v = a2.view();
    tmv::UpperTriMatrixView<std::complex<T> > ca1v = ca1.view();
    tmv::UpperTriMatrixView<std::complex<T> > ca2v = ca2.view();
    tmv::LowerTriMatrixView<T> a1t = a1.transpose();
    tmv::LowerTriMatrixView<T> a2t = a2.transpose();
    tmv::LowerTriMatrixView<std::complex<T> > ca1t = ca1.transpose();
    tmv::LowerTriMatrixView<std::complex<T> > ca2t = ca2.transpose();

    TestMatrixDivArith1(tmv::LU,a1v,a2t,ca1v,ca2t,"L/U");
    TestMatrixDivArith1(tmv::LU,a1t,a2v,ca1t,ca2v,"U/L");
    TestMatrixDivArith1(tmv::LU,a2v,a1t,ca2v,ca1t,"L/U");
    TestMatrixDivArith1(tmv::LU,a2t,a1v,ca2t,ca1v,"U/L");

#if (XTEST & 2)
    tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a1b(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca1b(cm);
    tmv::UpperTriMatrix<T,tmv::UnitDiag> a2b(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2b(cm);

    tmv::UpperTriMatrixView<T> a1bv = a1b.view();
    tmv::UpperTriMatrixView<T> a2bv = a2b.view();
    tmv::UpperTriMatrixView<std::complex<T> > ca1bv = ca1b.view();
    tmv::UpperTriMatrixView<std::complex<T> > ca2bv = ca2b.view();
    tmv::LowerTriMatrixView<T> a1bt = a1b.transpose();
    tmv::LowerTriMatrixView<T> a2bt = a2b.transpose();
    tmv::LowerTriMatrixView<std::complex<T> > ca1bt = ca1b.transpose();
    tmv::LowerTriMatrixView<std::complex<T> > ca2bt = ca2b.transpose();

    TestMatrixDivArith1(tmv::LU,a1v,a1bt,ca1v,ca1bt,"L/U");
    TestMatrixDivArith1(tmv::LU,a1t,a1bv,ca1t,ca1bv,"U/L");
    TestMatrixDivArith1(tmv::LU,a2v,a2bt,ca2v,ca2bt,"L/U");
    TestMatrixDivArith1(tmv::LU,a2t,a2bv,ca2t,ca2bv,"U/L");
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
