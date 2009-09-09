#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::DiagType D> inline void TestTriDiv() 
{
  const int N = 10;

  tmv::Matrix<T> m(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) m(i,j) = 0.4+0.02*i-0.05*j;
  m.diag().AddToAll(5);
  m.diag(1).AddToAll(0.32);
  m.diag(-1).AddToAll(0.91);

  tmv::UpperTriMatrix<T,D> a(m);
  m = a;
  a.SaveDiv();
  m.SaveDiv();

  tmv::Vector<T> b(N);
  for (int i=0;i<N;++i) b(i) = i+7.;

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
  }
  Assert(std::abs(m.Det()-a.Det()) < eps*std::abs(m.Det()),"Tri Det");

  tmv::Matrix<std::complex<T> > cm(m);
  cm += std::complex<T>(10,2);
  cm.diag(1) *= std::complex<T>(-0.5,-0.8);
  cm.diag(-1) *= std::complex<T>(-0.7,0.1);

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

  tmv::Matrix<T> m2(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j)  m2(i,j) = 5.-2*i+4*j;
  m2.diag().AddToAll(100.);
  tmv::Matrix<std::complex<T> > cm2 = m2*std::complex<T>(1,4);
  tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a3(m2);
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca3(m2);

  TestMatrixDivArith2<T>(tmv::LU,a.View(),m2.View(),ca.View(),cm2.View(),
      "SquareM/Tri");
  TestMatrixDivArith1<T>(tmv::LU,a.View(),a3.View(),ca.View(),ca3.View(),
      "Tri/Tri");
  TestMatrixDivArith1<T>(tmv::LU,m2.View(),a.View(),cm2.View(),ca.View(),
      "Tri/SquareM");

#ifdef XTEST
  tmv::UpperTriMatrix<T,tmv::UnitDiag> a2(m2);
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2(m2);

  tmv::Matrix<T> m3(2*N,N);
  m3.Rows(0,N) = m2;
  m3.Rows(N,2*N) = -m2;
  tmv::Matrix<std::complex<T> > cm3 = m3*std::complex<T>(2,-3);

  tmv::Matrix<T> m4(0,N);
  tmv::Matrix<std::complex<T> > cm4(0,N);

  tmv::DiagMatrix<T> d(m2);
  tmv::DiagMatrix<std::complex<T> > cd(cm2);

  TestMatrixDivArith1<T>(tmv::LU,a.View(),m3.View(),ca.View(),cm3.View(),
      "NonSquareM/Tri");
  TestMatrixDivArith1<T>(tmv::LU,a.View(),m4.View(),ca.View(),cm4.View(),
      "DegenM/Tri");
  TestMatrixDivArith1<T>(tmv::LU,a.View(),a2.View(),ca.View(),ca2.View(),
      "Tri/Tri");
  TestMatrixDivArith1<T>(tmv::LU,a.View(),d.View(),ca.View(),cd.View(),
      "Diag/Tri");
  TestMatrixDivArith1<T>(tmv::LU,d.View(),a.View(),cd.View(),ca.View(),
      "Tri/Diag");
#endif

}

template <class T> void TestAllTriDiv()
{
  TestTriDiv<T,tmv::NonUnitDiag>();
  TestTriDiv<T,tmv::UnitDiag>();
  std::cout<<"TriMatrix<"<<tmv::Type(T())<<"> Division passed all tests\n";
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
#ifdef INST_INT
template void TestAllTriDiv<int>();
#endif
