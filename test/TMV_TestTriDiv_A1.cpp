// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

template <class T1, class T2> inline bool CanLDivEq(
    const tmv::UpperTriMatrixView<T1>& a, const tmv::UpperTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanRDivEq(
    const tmv::UpperTriMatrixView<T1>& a, const tmv::UpperTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanLDivEq(
    const tmv::LowerTriMatrixView<T1>& a, const tmv::LowerTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanRDivEq(
    const tmv::LowerTriMatrixView<T1>& a, const tmv::LowerTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::DiagType D> void TestTriDiv() 
{
  const int N = 10;

  tmv::Matrix<T> m(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) m(i,j) = T(0.4+0.02*i-0.05*j);
  m.diag().AddToAll(5);
  m.diag(1).AddToAll(T(0.32));
  m.diag(-1).AddToAll(T(0.91));

  tmv::UpperTriMatrix<T,D> a(m);
  m = a;
  a.SaveDiv();
  m.SaveDiv();

  tmv::Vector<T> b(N);
  for (int i=0;i<N;++i) b(i) = T(i+7);

  a.SetDiv();
  m.DivideUsing(tmv::LU);
  m.SetDiv();

  if (showacc) {
    std::cout<<"b = "<<b<<std::endl;
    std::cout<<"a = "<<a<<std::endl;
    std::cout<<"m = "<<m<<std::endl;
  }

  T eps = EPS * Norm(m) * Norm(m.Inverse());

  tmv::Vector<T> x1 = b/a;
  tmv::Vector<T> x2 = b/m;
  if (showacc) {
    std::cout<<"x1 = "<<x1<<std::endl;
    std::cout<<"x2 = "<<x2<<std::endl;
    std::cout<<"a*x1-b = "<<a*x1-b<<std::endl;
    std::cout<<"m*x2-b = "<<m*x2-b<<std::endl;
    std::cout<<"Norm(x1-x2) = "<<Norm(x1-x2)<<std::endl;
    std::cout<<"EPS*Norm(x1) = "<<eps*Norm(x1)<<std::endl;
  }
  Assert(Norm(x1-x2) < eps*Norm(x1),"Tri b/a");

  x1 = b%a;
  x2 = b%m;
  if (showacc) {
    std::cout<<"x1 = "<<x1<<std::endl;
    std::cout<<"x2 = "<<x2<<std::endl;
    std::cout<<"x1*a-b = "<<x1*a-b<<std::endl;
    std::cout<<"x2*m-b = "<<x2*m-b<<std::endl;
    std::cout<<"Norm(x1-x2) = "<<Norm(x1-x2)<<std::endl;
    std::cout<<"EPS*Norm(b) = "<<eps*Norm(x1)<<std::endl;
  }
  Assert(Norm(x1-x2) < eps*Norm(x1),"Tri b%a");

  tmv::UpperTriMatrix<T,D> ainv = a.Inverse();
  tmv::Matrix<T> minv = m.Inverse();
  if (showacc) {
    std::cout<<"ainv = "<<ainv<<std::endl;
    std::cout<<"minv = "<<minv<<std::endl;
    std::cout<<"Norm(ainv-minv) = "<<Norm(ainv-minv)<<std::endl;
    std::cout<<"EPS*Norm(ainv) = "<<eps*Norm(ainv)<<std::endl;
  }
  Assert(Norm(ainv-minv) < eps*Norm(ainv),"Tri Inverse");

  if (showacc) {
    std::cout<<"a.Det = "<<a.Det()<<", m.Det = "<<m.Det()<<std::endl;
    std::cout<<"abs(adet-mdet) = "<<std::abs(a.Det()-m.Det());
    std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(m.Det())<<std::endl;
    std::cout<<"a.LogDet = "<<a.LogDet();
    std::cout<<", m.LogDet = "<<m.LogDet()<<std::endl;
  }
  Assert(std::abs(m.Det()-a.Det()) < eps*std::abs(m.Det()),"Tri Det");
  T asign, msign;
  Assert(std::abs(m.LogDet(&msign)-a.LogDet(&asign)) < N*eps,"Tri LogDet");
  Assert(std::abs(asign-msign) < N*eps,"Tri LogDet - sign");
  Assert(std::abs(a.Det() - asign*std::exp(a.LogDet())) < eps*std::abs(m.Det()),
      "Tri Det--LogDet");

  tmv::Matrix<std::complex<T> > cm(m);
  cm += std::complex<T>(10,2);
  cm.diag(1) *= std::complex<T>(T(-0.5),T(-0.8));
  cm.diag(-1) *= std::complex<T>(T(-0.7),T(0.1));

  tmv::UpperTriMatrix<std::complex<T>,D> ca(cm);
  cm = ca;
  ca.SaveDiv();
  cm.SaveDiv();

  cm.DivideUsing(tmv::LU);
  cm.SetDiv();
  ca.SetDiv();

  T ceps = EPS * Norm(cm) * Norm(cm.Inverse());

  if (showacc) {
    std::cout<<"ca.Det = "<<ca.Det()<<", cm.Det = "<<cm.Det()<<std::endl;
    std::cout<<"abs(cadet-cmdet) = "<<std::abs(ca.Det()-cm.Det());
    std::cout<<"  EPS*abs(cmdet) = "<<ceps*std::abs(cm.Det())<<std::endl;
  }
  Assert(std::abs(ca.Det()-cm.Det()) < ceps*std::abs(cm.Det()),"Tri CDet");
  std::complex<T> casign, cmsign;
  Assert(std::abs(cm.LogDet(&cmsign)-ca.LogDet(&casign)) < N*eps,"Tri CLogDet");
  Assert(std::abs(casign-cmsign) < N*eps,"Tri CLogDet - sign");
  Assert(std::abs(ca.Det() - casign*std::exp(ca.LogDet())) < 
      eps*std::abs(cm.Det()),"Tri CDet--LogDet");

  tmv::Vector<std::complex<T> > e(b);
  e(1) += std::complex<T>(-1,5);
  e(2) -= std::complex<T>(-1,5);

  tmv::Vector<std::complex<T> > y1 = b/ca;
  tmv::Vector<std::complex<T> > y2 = b/cm;
  if (showacc) {
    std::cout<<"y1 = "<<y1<<std::endl;
    std::cout<<"y2 = "<<y2<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm(y1-y2)<<std::endl;
    std::cout<<"EPS*Norm(y1) = "<<ceps*Norm(y1)<<std::endl;
  }
  Assert(Norm(y1-y2) < ceps*Norm(y1),"Tri b/c");
  y1 = b%ca;
  y2 = b%cm;
  if (showacc) {
    std::cout<<"y1 = "<<y1<<std::endl;
    std::cout<<"y2 = "<<y1<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm(y1-y2)<<std::endl;
    std::cout<<"EPS*Norm(y1) = "<<ceps*Norm(y1)<<std::endl;
  }
  Assert(Norm(y1-y2) < ceps*Norm(y1),"Tri b%c");

  // test complex / real
  y1 = e/a;
  y2 = e/m;
  Assert(Norm(y1-y2) < eps*Norm(y1),"Tri e/m");
  y1 = e%a;
  y2 = e%m;
  Assert(Norm(y1-y2) < eps*Norm(y1),"Tri e%m");

  // test complex / complex
  y1 = e/ca;
  y2 = e/cm;
  Assert(Norm(y1-y2) < ceps*Norm(y1),"Tri e/c");
  y1 = e%ca;
  y2 = e%cm;
  Assert(Norm(y1-y2) < ceps*Norm(y1),"Tri e%c");
}

template <class T> void TestTriDiv_A1() 
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

  TestMatrixDivArith2<T>(tmv::LU,a2x,ca2x,a1.View(),a2.View(),
      ca1.View(),ca2.View(),"U/U 1");
  TestMatrixDivArith2<T>(tmv::LU,b2x,cb2x,a1.Transpose(),a2.Transpose(),
      ca1.Transpose(),ca2.Transpose(),"L/L 1");
  TestMatrixDivArith2<T>(tmv::LU,a1x,ca1x,a2.View(),a1.View(),
      ca2.View(),ca1.View(),"U/U 2");
  TestMatrixDivArith2<T>(tmv::LU,b1x,cb1x,a2.Transpose(),a1.Transpose(),
      ca2.Transpose(),ca1.Transpose(),"L/L 2");

#ifdef XTEST
  tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a1b(m);
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca1b(cm);
  tmv::UpperTriMatrix<T,tmv::UnitDiag> a2b(m);
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2b(cm);

  TestMatrixDivArith1<T>(tmv::LU,a1x,ca1x,a1.View(),a1b.View(),
      ca1.View(),ca1b.View(),"U/U 3");
  TestMatrixDivArith1<T>(tmv::LU,b1x,cb1x,a1.Transpose(),a1b.Transpose(),
      ca1.Transpose(),ca1b.Transpose(),"L/L 3");
  TestMatrixDivArith1<T>(tmv::LU,a2x,ca2x,a2.View(),a2b.View(),
      ca2.View(),ca2b.View(),"U/U 4");
  TestMatrixDivArith1<T>(tmv::LU,b2x,cb2x,a2.Transpose(),a2b.Transpose(),
      ca2.Transpose(),ca2b.Transpose(),"L/L 4");
#endif
}

template <class T> void TestAllTriDiv()
{
  TestTriDiv<T,tmv::NonUnitDiag>();
  TestTriDiv<T,tmv::UnitDiag>();
  TestTriDiv_A1<T>();
  TestTriDiv_A2<T>();
  TestTriDiv_B1<T>();
  TestTriDiv_B2<T>();
  TestTriDiv_C1<T>();
  TestTriDiv_C2<T>();
  std::cout<<"TriMatrix<"<<tmv::TypeText(T())<<"> Division passed all tests\n";
}

#ifdef INST_DOUBLE
template void TestAllTriDiv<double>();
#endif
#ifdef INST_FLOAT
template void TestAllTriDiv<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllTriDiv<long double>();
#endif
