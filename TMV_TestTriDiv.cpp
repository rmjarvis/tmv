//#define SHOWACC
//#define SHOWTESTS
#define TESTDIV

#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Tri.h"
#include "TMV_TestMatrixDiv.h"

using tmv::Matrix;
using tmv::Vector;
using tmv::UpperTriMatrix;
using tmv::UnitDiag;
using tmv::NonUnitDiag;

template <class T, tmv::DiagType D> void TestTriDiv() 
{
  const int N = 10;

  Matrix<T> m(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) m(i,j) = 0.4+0.02*i-0.05*j;
  m.diag().AddToAll(5);
  m.diag(1).AddToAll(0.32);
  m.diag(-1).AddToAll(0.91);

  UpperTriMatrix<T,D> a(m);
  m = a;

  Vector<T> b(N);
  for (int i=0;i<N;++i) b(i) = i+7.;

  a.SetDiv();
  m.DivideUsing(tmv::LU);
  m.SetDiv();
  CheckDecomposition<T>(a,tmv::LU);

#ifdef SHOWACC
  cout<<"b = "<<b<<endl;
  cout<<"a = "<<a<<endl;
  cout<<"m = "<<m<<endl;
#endif

  Vector<T> x1 = b/a;
  Vector<T> x2 = b/m;
#ifdef SHOWACC
  cout<<"x1 = "<<x1<<endl;
  cout<<"x2 = "<<x2<<endl;
  cout<<"a*x1-b = "<<a*x1-b<<endl;
  cout<<"m*x2-b = "<<m*x2-b<<endl;
  cout<<"Norm(x1-x2) = "<<Norm(x1-x2)<<endl;
  cout<<"EPS*Norm(m)*Norm(minv)*Norm(b) = "<<EPS*Norm(m)*Norm(m.Inverse())*Norm(b)<<endl;
#endif
  Assert(Norm(x1-x2) < EPS*Norm(m)*Norm(m.Inverse()*Norm(b)),"Tri b/a");

  x1 = b%a;
  x2 = b%m;
#ifdef SHOWACC
  cout<<"x1 = "<<x1<<endl;
  cout<<"x2 = "<<x2<<endl;
  cout<<"x1*a-b = "<<x1*a-b<<endl;
  cout<<"x2*m-b = "<<x2*m-b<<endl;
  cout<<"Norm(x1-x2) = "<<Norm(x1-x2)<<endl;
  cout<<"EPS*Norm(m)*Norm(minv)*Norm(b) = "<<EPS*Norm(m)*Norm(m.Inverse())*Norm(b)<<endl;
#endif
  Assert(Norm(x1-x2) < EPS*Norm(m)*Norm(m.Inverse()*Norm(b)),"Tri b%a");

  UpperTriMatrix<T,D> ainv = a.TInverse();
  Matrix<T> minv = m.Inverse();
#ifdef SHOWACC
  cout<<"ainv = "<<ainv<<endl;
  cout<<"minv = "<<minv<<endl;
  cout<<"Norm(ainv-minv) = "<<Norm(ainv-minv)<<endl;
  cout<<"EPS*Norm(ainv)^2*Norm(a) = "<<EPS*Norm(ainv)*Norm(ainv)*Norm(a)<<endl;
#endif
  Assert(Norm(ainv-minv) < EPS*Norm(ainv)*Norm(ainv)*Norm(a),"Tri Inverse");

#ifdef SHOWACC
  cout<<"a.Det = "<<a.Det()<<", m.Det = "<<m.Det()<<endl;
  cout<<"abs(adet-mdet) = "<<abs(a.Det()-m.Det());
  cout<<"  EPS*abs(mdet) = "<<EPS*abs(m.Det())<<endl;
#endif
  Assert(abs(m.Det()-a.Det()) < EPS*m.colsize()*abs(m.Det()),"Tri Det");

  Matrix<complex<T> > cm(m);
  cm += complex<T>(10,2);
  cm.diag(1) *= complex<T>(-0.5,-0.8);
  cm.diag(-1) *= complex<T>(-0.7,0.1);

  UpperTriMatrix<complex<T>,D> ca(cm);
  cm = ca;

  cm.DivideUsing(tmv::LU);
  cm.SetDiv();
  ca.SetDiv();
  CheckDecomposition<complex<T> >(ca,tmv::LU);

#ifdef SHOWACC
  cout<<"ca.Det = "<<ca.Det()<<", cm.Det = "<<cm.Det()<<endl;
  cout<<"abs(cadet-cmdet) = "<<abs(ca.Det()-cm.Det());
  cout<<"  EPS*abs(cmdet) = "<<EPS*abs(cm.Det())<<endl;
#endif
  Assert(abs(ca.Det()-cm.Det()) < EPS*Norm(cm)*abs(cm.Det()),"Tri CDet");

  Vector<complex<T> > e(b);
  e(1) += complex<T>(-1,5);
  e(2) -= complex<T>(-1,5);

  Vector<complex<T> > y1 = b/ca;
  Vector<complex<T> > y2 = b/cm;
#ifdef SHOWACC
  cout<<"y1 = "<<y1<<endl;
  cout<<"y2 = "<<y2<<endl;
  cout<<"Norm(diff) = "<<Norm(y1-y2)<<endl;
  cout<<"EPS*Norm(cm)*Norm(cminv)*Norm(b) = "<<EPS*Norm(cm)*Norm(cm.Inverse())*Norm(b)<<endl;
#endif
  Assert(Norm(y1-y2) < EPS*Norm(cm)*Norm(cm.Inverse())*Norm(b),"Tri b/c");
  y1 = b%ca;
  y2 = b%cm;
#ifdef SHOWACC
  cout<<"y1 = "<<y1<<endl;
  cout<<"y2 = "<<y1<<endl;
  cout<<"Norm(diff) = "<<Norm(y1-y2)<<endl;
  cout<<"EPS*Norm(cm)*Norm(cminv)*Norm(b) = "<<EPS*Norm(cm)*Norm(cm.Inverse())*Norm(b)<<endl;
#endif
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

  Matrix<T> m2(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j)  m2(i,j) = 5.-2*i+4*j;
  m2.diag().AddToAll(100.);
  Matrix<complex<T> > cm2 = m2*complex<T>(1,4);

  UpperTriMatrix<T,UnitDiag> a2(m2);
  UpperTriMatrix<T,NonUnitDiag> a3(m2);
  UpperTriMatrix<complex<T>,UnitDiag> ca2(m2);
  UpperTriMatrix<complex<T>,NonUnitDiag> ca3(m2);

  Matrix<T> m3(2*N,N);
  m3.Rows(0,N) = m2;
  m3.Rows(N,2*N) = -m2;
  Matrix<complex<T> > cm3 = m3*complex<T>(2,-3);

  Matrix<T> m4(0,N);
  Matrix<complex<T> > cm4(0,N);

  TestMatrixDivArith<T>(tmv::LU,a,a2,ca,ca2,"Tri/Tri");
  TestMatrixDivArith<T>(tmv::LU,a,a3,ca,ca3,"Tri/Tri");
  TestMatrixDivArith<T>(tmv::LU,a,m2,ca,cm2,"SquareM/Tri");
  TestMatrixDivArith<T>(tmv::LU,a,m3,ca,cm3,"NonSquareM/Tri");
  TestMatrixDivArith<T>(tmv::LU,m2,a,cm2,ca,"Tri/SquareM");

  cout<<"TriMatrix<"<<tmv::Type(T())<<",";
  cout<<tmv::Text(D)<<"> Division passed all tests\n";
}

template <class T> void TestAllTriDiv()
{
  TestTriDiv<T,NonUnitDiag>();
  TestTriDiv<T,UnitDiag>();
}

template void TestAllTriDiv<double>();
#ifndef NOFLOAT
template void TestAllTriDiv<float>();
#endif
#ifdef LONGDOUBLE
template void TestAllTriDiv<long double>();
#endif
