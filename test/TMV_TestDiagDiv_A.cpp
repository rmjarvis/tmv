// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestDiagDiv_A()
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

  tmv::Matrix<T> r(N,0,T(1));
  tmv::Matrix<std::complex<T> > cr(N,0,std::complex<T>(1));

  tmv::DiagMatrix<T> bx(N);
  tmv::DiagMatrix<std::complex<T> > cbx(N);

  TestMatrixDivArith2<T>(tmv::LU,bx,cbx,a.View(),b.View(),ca.View(),cb.View(),
      "Diag/Diag");
}

template <class T> void TestDiagDiv()
{
  TestDiagDiv_A<T>();
  TestDiagDiv_B1<T>();
  TestDiagDiv_B2<T>();
  std::cout<<"DiagMatrix<"<<tmv::TypeText(T())<<"> Division passed all tests\n";
}

#ifdef INST_DOUBLE
template void TestDiagDiv<double>();
template void TestDiagDiv_A<double>();
#endif
#ifdef INST_FLOAT
template void TestDiagDiv<float>();
template void TestDiagDiv_A<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestDiagDiv<long double>();
template void TestDiagDiv_A<long double>();
#endif
