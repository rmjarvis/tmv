#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"

template <class T1, class T2> inline bool CanLDivEq(
    const tmv::UpperTriMatrixView<T1>& a, const tmv::DiagMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanRDivEq(
    const tmv::UpperTriMatrixView<T1>& a, const tmv::DiagMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanLDivEq(
    const tmv::LowerTriMatrixView<T1>& a, const tmv::DiagMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanRDivEq(
    const tmv::LowerTriMatrixView<T1>& a, const tmv::DiagMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestTriDiv_C1() 
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

    tmv::DiagMatrix<T> d(m);
    tmv::DiagMatrix<std::complex<T> > cd(cm);

    tmv::UpperTriMatrixView<T> a1v = a1.view();
    tmv::UpperTriMatrixView<std::complex<T> > ca1v = ca1.view();
    tmv::LowerTriMatrixView<T> a1t = a1.transpose();
    tmv::LowerTriMatrixView<std::complex<T> > ca1t = ca1.transpose();
    tmv::DiagMatrixView<T> dv = d.view();
    tmv::DiagMatrixView<std::complex<T> > cdv = cd.view();

    TestMatrixDivArith1(tmv::LU,dv,a1v,cdv,ca1v,"U/D");
    TestMatrixDivArith1(tmv::LU,dv,a1t,cdv,ca1t,"L/D");

#if (XTEST & 2)
    tmv::UpperTriMatrix<T,tmv::UnitDiag> a2(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2(cm);

    tmv::UpperTriMatrixView<T> a2v = a2.view();
    tmv::UpperTriMatrixView<std::complex<T> > ca2v = ca2.view();
    tmv::LowerTriMatrixView<T> a2t = a2.transpose();
    tmv::LowerTriMatrixView<std::complex<T> > ca2t = ca2.transpose();

    TestMatrixDivArith1(tmv::LU,dv,a2v,cdv,ca2v,"U/D");
    TestMatrixDivArith1(tmv::LU,dv,a2t,cdv,ca2t,"L/D");
#endif
}

#ifdef TEST_DOUBLE
template void TestTriDiv_C1<double>();
#endif
#ifdef TEST_FLOAT
template void TestTriDiv_C1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestTriDiv_C1<long double>();
#endif
