
#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestDiagDiv_B2()
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
    cq.rowRange(N,2*N) -= cp;
    cq.rowRange(N/2,3*N/2) -= cp;

    tmv::Matrix<T> r(N,0,T(1));
    tmv::Matrix<std::complex<T> > cr(N,0,std::complex<T>(1));

    tmv::MatrixView<T> pv = p.view();
    tmv::MatrixView<T> qv = q.view();
    tmv::MatrixView<T> rv = r.view();
    tmv::DiagMatrixView<T> av = a.view();
    tmv::MatrixView<std::complex<T> > cpv = cp.view();
    tmv::MatrixView<std::complex<T> > cqv = cq.view();
    tmv::MatrixView<std::complex<T> > crv = cr.view();
    tmv::DiagMatrixView<std::complex<T> > cav = ca.view();

    tmv::Matrix<T> px(N,N);
    tmv::DiagMatrix<T> bx(N);
    tmv::Matrix<T> qx(2*N,N);
    tmv::Matrix<T> rx(N,0);
    tmv::Matrix<std::complex<T> > cpx(N,N);
    tmv::DiagMatrix<std::complex<T> > cbx(N);
    tmv::Matrix<std::complex<T> > cqx(2*N,N);
    tmv::Matrix<std::complex<T> > crx(N,0);

    TestMatrixDivArith1<T>(tmv::LU,px,cpx,av,pv,cav,cpv,"SquareM/Diag");
    TestMatrixDivArith1<T>(tmv::LU,qx,cqx,av,qv,cav,cqv,"NonSqaureM/Diag");
#if (XTEST & 8)
    TestMatrixDivArith1<T>(tmv::LU,rx,crx,av,rv,cav,crv,"DegenM/Diag");
#endif
}

#ifdef TMV_INST_DOUBLE
template void TestDiagDiv_B2<double>();
#endif
#ifdef TMV_INST_FLOAT
template void TestDiagDiv_B2<float>();
#endif
#ifdef TMV_INST_LONGDOUBLE
template void TestDiagDiv_B2<long double>();
#endif