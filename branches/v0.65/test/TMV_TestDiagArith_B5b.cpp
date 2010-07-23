
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"
#include <fstream>

#define NOSV
#include "TMV_TestMatrixArith.h"

template <class T> void TestDiagMatrixArith_B5b()
{
#if (XTEST & 2)
    const int N = 10;

    tmv::DiagMatrix<T> a(N);
    tmv::DiagMatrix<T> b(N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        if (i == j) {
            a(i,j) = T(3+i+5*j);
            b(i,j) = T(5+2*i+4*j);
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

    tmv::DiagMatrix<T> ax = a;
    tmv::DiagMatrix<std::complex<T> > cax = ca;
    tmv::Matrix<T> px = p;
    tmv::Matrix<std::complex<T> > cpx = cp;
    tmv::Matrix<T> qx = q;
    tmv::Matrix<std::complex<T> > cqx = cq;
    tmv::Matrix<T> rx = r;
    tmv::Matrix<std::complex<T> > crx = cr;

    tmv::DiagMatrixView<T> av = a.view();
    tmv::DiagMatrixView<std::complex<T> > cav = ca.view();
    tmv::MatrixView<T> pv = p.view();
    tmv::MatrixView<std::complex<T> > cpv = cp.view();

    tmv::MatrixView<T> qv = q.view();
    tmv::MatrixView<std::complex<T> > cqv = cq.view();
    tmv::MatrixView<T> rv = r.view();
    tmv::MatrixView<std::complex<T> > crv = cr.view();

    TestMatrixArith5<T>(pv,cpv,av,cav, "SquareM/Diag");
    TestMatrixArith5<T>(qv,cqv,av,cav, "NonSquareM/Diag");
#if (XTEST & 8)
    TestMatrixArith5<T>(rv,crv,av,cav, "DegenM/Diag");
#endif
#endif
}

#ifdef TEST_DOUBLE
template void TestDiagMatrixArith_B5b<double>();
#endif
#ifdef TEST_FLOAT
template void TestDiagMatrixArith_B5b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestDiagMatrixArith_B5b<long double>();
#endif
#ifdef TEST_INT
template void TestDiagMatrixArith_B5b<int>();
#endif
