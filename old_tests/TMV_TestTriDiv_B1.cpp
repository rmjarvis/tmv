#include "TMV_Test.h"
#include "TMV_Test_1.h"
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

    tmv::UpperTriMatrixView<T> a1v = a1.view();
    tmv::UpperTriMatrixView<std::complex<T> > ca1v = ca1.view();
    tmv::LowerTriMatrixView<T> a1t = a1.transpose();
    tmv::LowerTriMatrixView<std::complex<T> > ca1t = ca1.transpose();
    tmv::MatrixView<T> mv = m.view();
    tmv::MatrixView<std::complex<T> > cmv = cm.view();
    mv.divideUsing(tmv::LU);
    mv.saveDiv();
    cmv.divideUsing(tmv::LU);
    cmv.saveDiv();

    TestMatrixDivArith1(tmv::LU,a1v,mv,ca1v,cmv,"M/U");
    TestMatrixDivArith1(tmv::LU,a1t,mv,ca1t,cmv,"M/L");

#if (XTEST & 2)
    tmv::UpperTriMatrix<T,tmv::UnitDiag> a2(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2(cm);

    tmv::UpperTriMatrixView<T> a2v = a2.view();
    tmv::UpperTriMatrixView<std::complex<T> > ca2v = ca2.view();
    tmv::LowerTriMatrixView<T> a2t = a1.transpose();
    tmv::LowerTriMatrixView<std::complex<T> > ca2t = ca2.transpose();

    tmv::MatrixView<T> mt = m.transpose();
    tmv::MatrixView<std::complex<T> > cmt = cm.transpose();
    mt.divideUsing(tmv::LU);
    mt.saveDiv();
    cmt.divideUsing(tmv::LU);
    cmt.saveDiv();

    TestMatrixDivArith1(tmv::LU,a2v,mv,ca2v,cmv,"M/U");
    TestMatrixDivArith1(tmv::LU,a2t,mv,ca2t,cmv,"M/L");
    TestMatrixDivArith1(tmv::LU,a2v,mt,ca2v,cmt,"M/U");
    TestMatrixDivArith1(tmv::LU,a2t,mt,ca2t,cmt,"M/L");
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
