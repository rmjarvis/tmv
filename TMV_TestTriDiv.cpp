#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Tri.h"
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
  if (showdiv) {
    Assert(a.CheckDecomp(&std::cout),"CheckDecomp");
  } else {
    Assert(a.CheckDecomp(),"CheckDecomp");
  }

  if (showacc) {
    std::cout<<"b = "<<b<<std::endl;
    std::cout<<"a = "<<a<<std::endl;
    std::cout<<"m = "<<m<<std::endl;
  }

  tmv::Vector<T> x1 = b/a;
  tmv::Vector<T> x2 = b/m;
  if (showacc) {
    std::cout<<"x1 = "<<x1<<std::endl;
    std::cout<<"x2 = "<<x2<<std::endl;
    std::cout<<"a*x1-b = "<<a*x1-b<<std::endl;
    std::cout<<"m*x2-b = "<<m*x2-b<<std::endl;
    std::cout<<"Norm(x1-x2) = "<<Norm(x1-x2)<<std::endl;
    std::cout<<"EPS*Norm(m)*Norm(minv)*Norm(b) = "<<EPS*Norm(m)*Norm(m.Inverse())*Norm(b)<<std::endl;
  }
  Assert(Norm(x1-x2) < EPS*Norm(m)*Norm(m.Inverse())*Norm(b),"Tri b/a");

  x1 = b%a;
  x2 = b%m;
  if (showacc) {
    std::cout<<"x1 = "<<x1<<std::endl;
    std::cout<<"x2 = "<<x2<<std::endl;
    std::cout<<"x1*a-b = "<<x1*a-b<<std::endl;
    std::cout<<"x2*m-b = "<<x2*m-b<<std::endl;
    std::cout<<"Norm(x1-x2) = "<<Norm(x1-x2)<<std::endl;
    std::cout<<"EPS*Norm(m)*Norm(minv)*Norm(b) = "<<EPS*Norm(m)*Norm(m.Inverse())*Norm(b)<<std::endl;
  }
  Assert(Norm(x1-x2) < EPS*Norm(m)*Norm(m.Inverse())*Norm(b),"Tri b%a");

  tmv::UpperTriMatrix<T,D> ainv = a.Inverse();
  tmv::Matrix<T> minv = m.Inverse();
  if (showacc) {
    std::cout<<"ainv = "<<ainv<<std::endl;
    std::cout<<"minv = "<<minv<<std::endl;
    std::cout<<"Norm(ainv-minv) = "<<Norm(ainv-minv)<<std::endl;
    std::cout<<"EPS*Norm(ainv)^2*Norm(a) = "<<EPS*Norm(ainv)*Norm(ainv)*Norm(a)<<std::endl;
  }
  Assert(Norm(ainv-minv) < EPS*Norm(ainv)*Norm(ainv)*Norm(a),"Tri Inverse");

  if (showacc) {
    std::cout<<"a.Det = "<<a.Det()<<", m.Det = "<<m.Det()<<std::endl;
    std::cout<<"abs(adet-mdet) = "<<std::abs(a.Det()-m.Det());
    std::cout<<"  EPS*abs(mdet) = "<<EPS*std::abs(m.Det())<<std::endl;
  }
  Assert(std::abs(m.Det()-a.Det()) < EPS*m.colsize()*std::abs(m.Det()),"Tri Det");

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
  if (showdiv) {
    Assert(ca.CheckDecomp(&std::cout),"CheckDecomp");
  } else {
    Assert(ca.CheckDecomp(),"CheckDecomp");
  }

  if (showacc) {
    std::cout<<"ca.Det = "<<ca.Det()<<", cm.Det = "<<cm.Det()<<std::endl;
    std::cout<<"abs(cadet-cmdet) = "<<std::abs(ca.Det()-cm.Det());
    std::cout<<"  EPS*abs(cmdet) = "<<EPS*std::abs(cm.Det())<<std::endl;
  }
  Assert(std::abs(ca.Det()-cm.Det()) < EPS*Norm(cm)*std::abs(cm.Det()),"Tri CDet");

  tmv::Vector<std::complex<T> > e(b);
  e(1) += std::complex<T>(-1,5);
  e(2) -= std::complex<T>(-1,5);

  tmv::Vector<std::complex<T> > y1 = b/ca;
  tmv::Vector<std::complex<T> > y2 = b/cm;
  if (showacc) {
    std::cout<<"y1 = "<<y1<<std::endl;
    std::cout<<"y2 = "<<y2<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm(y1-y2)<<std::endl;
    std::cout<<"EPS*Norm(cm)*Norm(cminv)*Norm(b) = "<<EPS*Norm(cm)*Norm(cm.Inverse())*Norm(b)<<std::endl;
  }
  Assert(Norm(y1-y2) < EPS*Norm(cm)*Norm(cm.Inverse())*Norm(b),"Tri b/c");
  y1 = b%ca;
  y2 = b%cm;
  if (showacc) {
    std::cout<<"y1 = "<<y1<<std::endl;
    std::cout<<"y2 = "<<y1<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm(y1-y2)<<std::endl;
    std::cout<<"EPS*Norm(cm)*Norm(cminv)*Norm(b) = "<<EPS*Norm(cm)*Norm(cm.Inverse())*Norm(b)<<std::endl;
  }
  Assert(Norm(y1-y2) < EPS*Norm(cm)*Norm(cm.Inverse())*Norm(b),"Tri b%c");

  // test complex / real
  y1 = e/a;
  y2 = e/m;
  Assert(Norm(y1-y2) < EPS*Norm(m)*Norm(m.Inverse())*Norm(e),"Tri e/m");
  y1 = e%a;
  y2 = e%m;
  Assert(Norm(y1-y2) < EPS*Norm(m)*Norm(m.Inverse())*Norm(e),"Tri e%m");

  // test complex / complex
  y1 = e/ca;
  y2 = e/cm;
  Assert(Norm(y1-y2) < EPS*Norm(cm)*Norm(cm.Inverse())*Norm(e),"Tri e/c");
  y1 = e%ca;
  y2 = e%cm;
  Assert(Norm(y1-y2) < EPS*Norm(cm)*Norm(cm.Inverse())*Norm(e),"Tri e%c");

  tmv::Matrix<T> m2(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j)  m2(i,j) = 5.-2*i+4*j;
  m2.diag().AddToAll(100.);
  tmv::Matrix<std::complex<T> > cm2 = m2*std::complex<T>(1,4);

  tmv::UpperTriMatrix<T,tmv::UnitDiag> a2(m2);
  tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a3(m2);
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2(m2);
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca3(m2);

  tmv::Matrix<T> m3(2*N,N);
  m3.Rows(0,N) = m2;
  m3.Rows(N,2*N) = -m2;
  tmv::Matrix<std::complex<T> > cm3 = m3*std::complex<T>(2,-3);

  tmv::Matrix<T> m4(0,N);
  tmv::Matrix<std::complex<T> > cm4(0,N);

  TestMatrixDivArith<T>(tmv::LU,a.View(),m2.View(),ca.View(),cm2.View(),
      "SquareM/Tri");
#ifdef XTEST
  TestMatrixDivArith<T>(tmv::LU,a.View(),m3.View(),ca.View(),cm3.View(),
      "NonSquareM/Tri");
  TestMatrixDivArith<T>(tmv::LU,a.View(),a2.View(),ca.View(),ca2.View(),
      "Tri/Tri");
  TestMatrixDivArith<T>(tmv::LU,a.View(),a3.View(),ca.View(),ca3.View(),
      "Tri/Tri");
  TestMatrixDivArith<T>(tmv::LU,m2.View(),a.View(),cm2.View(),ca.View(),
      "Tri/SquareM");
#endif

}

template <class T> void TestAllTriDiv()
{
  TestTriDiv<T,tmv::NonUnitDiag>();
  TestTriDiv<T,tmv::UnitDiag>();
  std::cout<<"TriMatrix<"<<tmv::Type(T())<<"> Division passed all tests\n";
}

template void TestAllTriDiv<double>();
#ifndef NOFLOAT
template void TestAllTriDiv<float>();
#endif
#ifdef LONGDOUBLE
template void TestAllTriDiv<long double>();
#endif
