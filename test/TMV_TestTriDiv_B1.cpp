// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestTriDiv_B1() 
{
  const int N = 10;

  tmv::Matrix<T> m(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) m(i,j) = T(0.4+0.02*i-0.05*j);
  m.diag().AddToAll(5);
  m.diag(1).AddToAll(T(0.32));
  m.diag(-1).AddToAll(T(0.91));

  tmv::Matrix<std::complex<T> > cm(m);
  cm += std::complex<T>(10,2);
  cm.diag(1) *= std::complex<T>(T(-0.5),T(-0.8));
  cm.diag(-1) *= std::complex<T>(T(-0.7),T(0.1));

  tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a1(m);
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca1(cm);

  tmv::Matrix<T> mx(m);
  tmv::Matrix<std::complex<T> > cmx(cm);
  m.DivideUsing(tmv::LU);
  m.SaveDiv();
  m.SetDiv();

  TestMatrixDivArith1<T>(tmv::LU,mx,cmx,a1.View(),m.View(),
      ca1.View(),cm.Transpose(),"M/U");
  TestMatrixDivArith1<T>(tmv::LU,mx,cmx,a1.Transpose(),m.View(),
      ca1.Transpose(),cm.View(),"M/L");

#ifdef XTEST
  tmv::UpperTriMatrix<T,tmv::UnitDiag> a2(m);
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2(cm);
  TestMatrixDivArith1<T>(tmv::LU,mx,cmx,a2.View(),m.View(),
      ca2.View(),cm.Transpose(),"M/U");
  TestMatrixDivArith1<T>(tmv::LU,mx,cmx,a2.Transpose(),m.View(),
      ca2.Transpose(),cm.View(),"M/L");
#endif
}

#ifdef TEST_DOUBLE
template void TestTriDiv_B1<double>();
#endif
#ifdef TEST_FLOAT
template void TestTriDiv_B1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestTriDiv_B1<long double>();
#endif
