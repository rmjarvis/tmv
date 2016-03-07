
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"
#include <fstream>

#define NOADDEQ
#define NOMULTEQ
#include "TMV_TestMatrixArith.h"

template <class T> void TestDiagMatrixArith_B5a()
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

    tmv::Matrix<T> p(N,N,5);
    p.diag().addToAll(N*10);
    tmv::Matrix<std::complex<T> > cp = p*std::complex<T>(2,3);
    cp += ca;

    tmv::Matrix<T> q(2*N,N,-2);
    q.rowRange(0,N) += p;
    q.rowRange(N,2*N) -= p;
    q.rowRange(N/2,3*N/2) += T(4)*p;
    tmv::Matrix<std::complex<T> > cq = q*std::complex<T>(-1,4);
    cq.rowRange(0,N) -= cp;
    cq.rowRange(N,2*N) += cp;
    cq.rowRange(N/2,3*N/2) -= T(4)*cp;

    tmv::Matrix<T> r(N,0,T(1));
    tmv::Matrix<std::complex<T> > cr(N,0,T(1));

    tmv::DiagMatrixView<T> av = a.view();
    tmv::DiagMatrixView<std::complex<T> > cav = ca.view();
    tmv::MatrixView<T> pv = p.view();
    tmv::MatrixView<std::complex<T> > cpv = cp.view();

    TestMatrixArith5(av,cav,pv,cpv, "Diag/SquareM");

    tmv::MatrixView<T> qv = q.view();
    tmv::MatrixView<std::complex<T> > cqv = cq.view();
    tmv::MatrixView<T> rv = r.view();
    tmv::MatrixView<std::complex<T> > crv = cr.view();

    TestMatrixArith5(av,cav,qv,cqv, "Diag/NonSquareM");
    TestMatrixArith5(av,cav,rv,crv, "Diag/DegenM");
}

#ifdef TEST_DOUBLE
template void TestDiagMatrixArith_B5a<double>();
#endif
#ifdef TEST_FLOAT
template void TestDiagMatrixArith_B5a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestDiagMatrixArith_B5a<long double>();
#endif
#ifdef TEST_INT
template void TestDiagMatrixArith_B5a<int>();
#endif
