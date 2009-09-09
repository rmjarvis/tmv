// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV_Diag.h"
#include "TMV_Mat.h"
#include <fstream>

#define NODIV
#define NOSV
#include "TMV_TestMatrixArith.h"

template <class T> void TestDiagMatrixArith_A6()
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

  tmv::DiagMatrixView<T> bv = b.View();
  tmv::DiagMatrixView<std::complex<T> > cbv = cb.View();

  tmv::DiagMatrix<T> c = a;
  tmv::DiagMatrix<std::complex<T> > cc = ca;
  tmv::DiagMatrixView<T> cv = c.View();
  tmv::DiagMatrixView<std::complex<T> > ccv = cc.View();

  tmv::Matrix<T> m = a;
  tmv::Matrix<std::complex<T> > cm = ca;
  tmv::MatrixView<T> mv = m.View();
  tmv::MatrixView<std::complex<T> > cmv = cm.View();

  TestMatrixArith6<T>(av,cav,bv,cbv,cv,ccv, "Diag/Diag");
  TestMatrixArith6<T>(av,cav,bv,cbv,mv,cmv, "Diag/Diag/Matrix");
}

#ifdef TEST_DOUBLE
template void TestDiagMatrixArith_A6<double>();
#endif
#ifdef TEST_FLOAT
template void TestDiagMatrixArith_A6<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestDiagMatrixArith_A6<long double>();
#endif
#ifdef TEST_INT
template void TestDiagMatrixArith_A6<int>();
#endif
