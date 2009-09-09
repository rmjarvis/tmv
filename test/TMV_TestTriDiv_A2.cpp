// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> void TestTriDiv_A2() 
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
  tmv::UpperTriMatrix<T,tmv::UnitDiag> a2(m);
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2(cm);
  a1.SaveDiv();
  a2.SaveDiv();
  ca1.SaveDiv();
  ca2.SaveDiv();
  a1.SetDiv();
  a2.SetDiv();
  ca1.SetDiv();
  ca2.SetDiv();

  tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a1x(m);
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca1x(cm);
  tmv::UpperTriMatrix<T,tmv::UnitDiag> a2x(m);
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2x(cm);
  tmv::LowerTriMatrix<T,tmv::NonUnitDiag> b1x(m);
  tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag> cb1x(cm);
  tmv::LowerTriMatrix<T,tmv::UnitDiag> b2x(m);
  tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cb2x(cm);

  TestMatrixDivArith1<T>(tmv::LU,b2x,cb2x,a1.View(),a2.Transpose(),
      ca1.View(),ca2.Transpose(),"L/U");
  TestMatrixDivArith1<T>(tmv::LU,a2x,ca2x,a1.Transpose(),a2.View(),
      ca1.Transpose(),ca2.View(),"U/L");
  TestMatrixDivArith1<T>(tmv::LU,b1x,cb1x,a2.View(),a1.Transpose(),
      ca2.View(),ca1.Transpose(),"L/U");
  TestMatrixDivArith1<T>(tmv::LU,a1x,ca1x,a2.Transpose(),a1.View(),
      ca2.Transpose(),ca1.View(),"U/L");

#ifdef XTEST
  tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a1b(m);
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca1b(cm);
  tmv::UpperTriMatrix<T,tmv::UnitDiag> a2b(m);
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2b(cm);

  TestMatrixDivArith1<T>(tmv::LU,b1x,cb1x,a1.View(),a1b.Transpose(),
      ca1.View(),ca1b.Transpose(),"L/U");
  TestMatrixDivArith1<T>(tmv::LU,a1x,ca1x,a1.Transpose(),a1b.View(),
      ca1.Transpose(),ca1b.View(),"U/L");
  TestMatrixDivArith1<T>(tmv::LU,b2x,cb2x,a2.View(),a2b.Transpose(),
      ca2.View(),ca2b.Transpose(),"L/U");
  TestMatrixDivArith1<T>(tmv::LU,a2x,ca2x,a2.Transpose(),a2b.View(),
      ca2.Transpose(),ca2b.View(),"U/L");
#endif
}

#ifdef INST_DOUBLE
template void TestTriDiv_A2<double>();
#endif
#ifdef INST_FLOAT
template void TestTriDiv_A2<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestTriDiv_A2<long double>();
#endif
