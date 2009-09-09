// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV_Diag.h"
#include "TMV_Mat.h"
#include <fstream>

#define NOADDEQ
#define NOMULTEQ
#define NOSV
#include "TMV_TestMatrixArith.h"

template <class T> void TestDiagMatrixArith_B6a()
{
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
  p.diag().AddToAll(N*10);
  tmv::Matrix<std::complex<T> > cp = p*std::complex<T>(2,3);
  cp += ca;

  tmv::Matrix<T> q(2*N,N,-2);
  q.Rows(0,N) += p;
  q.Rows(N,2*N) -= p;
  q.Rows(N/2,3*N/2) += T(4)*p;
  tmv::Matrix<std::complex<T> > cq = q*std::complex<T>(-1,4);
  cq.Rows(0,N) -= cp;
  cq.Rows(N,2*N) += cp;
  cq.Rows(N/2,3*N/2) -= T(4)*cp;

  tmv::Matrix<T> r(N,0,T(1));
  tmv::Matrix<std::complex<T> > cr(N,0,T(1));

  tmv::DiagMatrixView<T> av = a.View();
  tmv::DiagMatrixView<std::complex<T> > cav = ca.View();
  tmv::MatrixView<T> pv = p.View();
  tmv::MatrixView<std::complex<T> > cpv = cp.View();

  TestMatrixArith6x<T>(av,cav,pv,cpv, "Diag/SquareM");

  tmv::MatrixView<T> qv = q.View();
  tmv::MatrixView<std::complex<T> > cqv = cq.View();
  tmv::MatrixView<T> rv = r.View();
  tmv::MatrixView<std::complex<T> > crv = cr.View();

  TestMatrixArith6x<T>(av,cav,qv,cqv, "Diag/NonSquareM");
  TestMatrixArith6x<T>(av,cav,rv,crv, "Diag/DegenM");
}

#ifdef TEST_DOUBLE
template void TestDiagMatrixArith_B6a<double>();
#endif
#ifdef TEST_FLOAT
template void TestDiagMatrixArith_B6a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestDiagMatrixArith_B6a<long double>();
#endif
#ifdef TEST_INT
template void TestDiagMatrixArith_B6a<int>();
#endif
