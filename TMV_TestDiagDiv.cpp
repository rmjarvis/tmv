
#include "TMV.h"
#include "TMV_Diag.h"
#include "TMV_Test.h"
#include "TMV_TestMatrixDivArith.h"

template <class T> void TestDiagDiv()
{
  const int N = 10;

  tmv::DiagMatrix<T> a(N);
  tmv::DiagMatrix<T> b(N);
  for (int i=0; i<N; ++i) {
    a(i,i) = 3.+5*i;
    b(i,i) = 5.-2*i;
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

  TestMatrixDivArith<T>(tmv::LU,a.View(),p.View(),ca.View(),cp.View(),
      "SquareM/Diag");
#ifdef XTEST
  TestMatrixDivArith<T>(tmv::LU,a.View(),q.View(),ca.View(),cq.View(),
      "NonSqaureM/Diag");
  TestMatrixDivArith<T>(tmv::LU,a.View(),r.View(),ca.View(),cr.View(),
      "DegenM/Diag");
  TestMatrixDivArith<T>(tmv::LU,a.View(),b.View(),ca.View(),cb.View(),
      "Diag/Diag");
  TestMatrixDivArith<T>(tmv::LU,p.View(),b.View(),cp.View(),cb.View(),
      "Diag/SquareM");
#endif

  std::cout<<"DiagMatrix<"<<tmv::Type(T())<<"> Division passed all tests\n";
}

template void TestDiagDiv<double>();
#ifndef NOFLOAT
template void TestDiagDiv<float>();
#endif
#ifdef LONGDOUBLE
template void TestDiagDiv<long double>();
#endif
