#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestTriDiv_C2() 
{
    const int N = 10;

    tmv::Matrix<T> m(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        m(i,j) = T(0.4+0.02*i-0.05*j);
    m.diag().AddToAll(5);
    m.diag(1).AddToAll(T(0.32));
    m.diag(-1).AddToAll(T(0.91));

    tmv::Matrix<std::complex<T> > cm(m);
    cm += std::complex<T>(10,2);
    cm.diag(1) *= std::complex<T>(T(-0.5),T(-0.8));
    cm.diag(-1) *= std::complex<T>(T(-0.7),T(0.1));

    tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a1(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca1(cm);
    a1.SaveDiv();
    ca1.SaveDiv();
    a1.SetDiv();
    ca1.SetDiv();

    tmv::DiagMatrix<T> d(m);
    tmv::DiagMatrix<std::complex<T> > cd(cm);
    tmv::DiagMatrix<T> dx(m);
    tmv::DiagMatrix<std::complex<T> > cdx(cm);

    TestMatrixDivArith1<T>(tmv::LU,dx,cdx,a1.View(),d.View(),
                           ca1.View(),cd.View(),"D/U");
    TestMatrixDivArith1<T>(tmv::LU,dx,cdx,a1.Transpose(),d.View(),
                           ca1.Transpose(),cd.View(),"D/L");

#ifdef XTEST
    tmv::UpperTriMatrix<T,tmv::UnitDiag> a2(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2(cm);
    a2.SaveDiv();
    ca2.SaveDiv();
    a2.SetDiv();
    ca2.SetDiv();

    TestMatrixDivArith1<T>(tmv::LU,dx,cdx,a2.View(),d.View(),
                           ca2.View(),cd.View(),"D/U");
    TestMatrixDivArith1<T>(tmv::LU,dx,cdx,a2.Transpose(),d.View(),
                           ca2.Transpose(),cd.View(),"D/L");
#endif
}

#ifdef INST_DOUBLE
template void TestTriDiv_C2<double>();
#endif
#ifdef INST_FLOAT
template void TestTriDiv_C2<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestTriDiv_C2<long double>();
#endif
