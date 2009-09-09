// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

template <class T1, class T2> inline bool CanLDivEq(
    const tmv::UpperTriMatrixView<T1>& a, const tmv::DiagMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanRDivEq(
    const tmv::UpperTriMatrixView<T1>& a, const tmv::DiagMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanLDivEq(
    const tmv::LowerTriMatrixView<T1>& a, const tmv::DiagMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanRDivEq(
    const tmv::LowerTriMatrixView<T1>& a, const tmv::DiagMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestTriDiv_C1() 
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
  tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a1x(m);
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca1x(cm);
  tmv::LowerTriMatrix<T,tmv::NonUnitDiag> b1x(m);
  tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag> cb1x(cm);

  tmv::DiagMatrix<T> d(m);
  tmv::DiagMatrix<std::complex<T> > cd(cm);

  TestMatrixDivArith1<T>(tmv::LU,a1x,ca1x,d.View(),a1.View(),
      cd.View(),ca1.View(),"U/D");
  TestMatrixDivArith1<T>(tmv::LU,b1x,cb1x,d.View(),a1.Transpose(),
      cd.View(),ca1.Transpose(),"L/D");

#ifdef XTEST
  tmv::UpperTriMatrix<T,tmv::UnitDiag> a2(m);
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2(cm);
  tmv::UpperTriMatrix<T,tmv::UnitDiag> a2x(m);
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2x(cm);
  tmv::LowerTriMatrix<T,tmv::UnitDiag> b2x(m);
  tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cb2x(cm);

  TestMatrixDivArith1<T>(tmv::LU,a2x,ca2x,d.View(),a2.View(),
      cd.View(),ca2.View(),"U/D");
  TestMatrixDivArith1<T>(tmv::LU,b2x,cb2x,d.View(),a2.Transpose(),
      cd.View(),ca2.Transpose(),"L/D");
#endif
}

#ifdef INST_DOUBLE
template void TestTriDiv_C1<double>();
#endif
#ifdef INST_FLOAT
template void TestTriDiv_C1<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestTriDiv_C1<long double>();
#endif
