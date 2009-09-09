// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestDiagDiv_B2()
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
  p.diag().AddToAll(N*10);
  tmv::Matrix<std::complex<T> > cp = p*std::complex<T>(2,3);
  cp += ca;

  tmv::Matrix<T> q(2*N,N,-2);
  q.Rows(0,N) += p;
  q.Rows(N,2*N) -= p;
  q.Rows(N/2,3*N/2) += T(4)*p;
  tmv::Matrix<std::complex<T> > cq = q*std::complex<T>(-1,4);
  cq.Rows(0,N) -= cp;
  cq.Rows(N,2*N) -= cp;
  cq.Rows(N/2,3*N/2) -= cp;

  tmv::Matrix<T> r(N,0,T(1));
  tmv::Matrix<std::complex<T> > cr(N,0,std::complex<T>(1));

  tmv::Matrix<T> px(N,N);
  tmv::DiagMatrix<T> bx(N);
  tmv::Matrix<T> qx(2*N,N);
  tmv::Matrix<T> rx(N,0);
  tmv::Matrix<std::complex<T> > cpx(N,N);
  tmv::DiagMatrix<std::complex<T> > cbx(N);
  tmv::Matrix<std::complex<T> > cqx(2*N,N);
  tmv::Matrix<std::complex<T> > crx(N,0);

  TestMatrixDivArith1<T>(tmv::LU,px,cpx,a.View(),p.View(),ca.View(),cp.View(),
      "SquareM/Diag");
#ifdef XTEST
  TestMatrixDivArith1<T>(tmv::LU,qx,cqx,a.View(),q.View(),ca.View(),cq.View(),
      "NonSqaureM/Diag");
  TestMatrixDivArith1<T>(tmv::LU,rx,crx,a.View(),r.View(),ca.View(),cr.View(),
      "DegenM/Diag");
#endif
}

#ifdef INST_DOUBLE
template void TestDiagDiv_B2<double>();
#endif
#ifdef INST_FLOAT
template void TestDiagDiv_B2<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestDiagDiv_B2<long double>();
#endif
