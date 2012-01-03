#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"

template <class T1, class T2> inline bool CanLDivEq(
    const tmv::UpperTriMatrixView<T1>& a, const tmv::UpperTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanRDivEq(
    const tmv::UpperTriMatrixView<T1>& a, const tmv::UpperTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanLDivEq(
    const tmv::LowerTriMatrixView<T1>& a, const tmv::LowerTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanRDivEq(
    const tmv::LowerTriMatrixView<T1>& a, const tmv::LowerTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestTriDiv_A1() 
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

    TestMatrixDivArith2(tmv::LU,a1v,a2v,ca1v,ca2v,"U/U 1");
    TestMatrixDivArith2(tmv::LU,a1t,a2t,ca1t,ca2t,"L/L 1");
    TestMatrixDivArith2(tmv::LU,a2v,a1v,ca2v,ca1v,"U/U 2");
    TestMatrixDivArith2(tmv::LU,a2t,a1t,ca2t,ca1t,"L/L 2");

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

    TestMatrixDivArith1(tmv::LU,a1v,a1bv,ca1v,ca1bv,"U/U 3");
    TestMatrixDivArith1(tmv::LU,a1t,a1bt,ca1t,ca1bt,"L/L 3");
    TestMatrixDivArith1(tmv::LU,a2v,a2bv,ca2v,ca2bv,"U/U 4");
    TestMatrixDivArith1(tmv::LU,a2t,a2bt,ca2t,ca2bt,"L/L 4");
#endif
}

#ifdef TEST_DOUBLE
template void TestTriDiv_A1<double>();
#endif
#ifdef TEST_FLOAT
template void TestTriDiv_A1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestTriDiv_A1<long double>();
#endif
