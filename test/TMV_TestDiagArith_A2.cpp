// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV_Diag.h"
#include "TMV_Mat.h"
#include <fstream>

#define NODIV
#define NOSV
#include "TMV_TestMatrixArith.h"

template <class T> void TestDiagMatrixArith_A2()
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

  tmv::DiagMatrixView<T> av = a.View();
  tmv::DiagMatrixView<std::complex<T> > cav = ca.View();

  TestMatrixArith2<T>(av,cav, "Diag");
}

#ifdef TEST_DOUBLE
template void TestDiagMatrixArith_A2<double>();
#endif
#ifdef TEST_FLOAT
template void TestDiagMatrixArith_A2<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestDiagMatrixArith_A2<long double>();
#endif
#ifdef TEST_INT
template void TestDiagMatrixArith_A2<int>();
#endif
