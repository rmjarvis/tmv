//#define SHOWACC
//#define SHOWTESTS
//#define SHOWCHECK
#define TESTDIV

#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_TestMatrixDiv.h"

using tmv::Matrix;
using tmv::Vector;

template <class T> void TestSquareDiv(tmv::DivType dt)
{
  Matrix<T> m(4,4);
  m.DivideUsing(dt);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) m(i,j) = 2.+4*i-5*j;
  m(0,0) = 14.;
  m(1,0) = -2.;
  m(2,0) = 7.;
  m(3,0) = -10.;
  m(2,2) = 30.;

  Vector<T> b(4);
  b(0) = 2;
  b(1) = -10;
  b(2) = 5;
  b(3) = -5;

  m.SetDiv();
  CheckDecomposition<T>(m,dt);
  Vector<T> x = b/m;
  Vector<T> b2 = m*x;
#ifdef SHOWACC
  cout<<"b = "<<b<<endl;
  cout<<"x = b/m = "<<x<<endl;
  cout<<"b2 = "<<b2<<endl;
  cout<<"Norm(b-b2) = "<<Norm(b-b2)<<endl;
#endif
  Assert(Norm(b2-b) < EPS*(Norm(m)*Norm(b)),"Square b/m");

  x = b%m;
  b2 = x*m;
#ifdef SHOWACC
  cout<<"b = "<<b<<endl;
  cout<<"x = b%m = "<<x<<endl;
  cout<<"b2 = "<<b2<<endl;
  cout<<"Norm(b-b2) = "<<Norm(b-b2)<<endl;
#endif
  Assert(Norm(b2-b) < EPS*(Norm(m)*Norm(b)),"Square b%m");

  Matrix<T> minv = m.Inverse();
  Matrix<T> id = m*minv;
  Assert(Norm(id-tmv::Eye<T>(4)) < EPS*Norm(minv)*Norm(m),
      "Square Inverse");

  Matrix<T> mtm = m.Adjoint() * m;
  Assert(Norm(m.InverseATA()-mtm.Inverse()) < EPS*Norm(mtm)*Norm(mtm.Inverse()),
	"Square InverseATA");

  T mdet = 28800.;
#ifdef SHOWACC
  cout<<"abs(det-mdet) = "<<abs(m.Det()-mdet);
  cout<<"  EPS*abs(mdet) = "<<EPS*abs(mdet)<<endl;
#endif
  Assert(abs(m.Det()-mdet) < EPS*m.colsize()*abs(mdet),"Square Det");

  Matrix<complex<T> > c(4,4);
  c.DivideUsing(dt);
  c = m;
  c(2,3) += complex<T>(2,3);
  c(1,0) *= complex<T>(0,2);
  c.col(1) *= complex<T>(-1,3);
  c.row(3) += Vector<complex<T> >(4,complex<T>(1,9));

  c.SetDiv();
  CheckDecomposition<complex<T> >(c,dt);
  complex<T> cdet(-103604,101272);
#ifdef SHOWACC
  cout<<"cdet = "<<cdet<<endl;
  cout<<"C.Det = "<<c.Det()<<endl;
  cout<<"abs(det-cdet) = "<<abs(c.Det()-cdet);
  cout<<"  EPS*abs(cdet) = "<<EPS*abs(cdet)<<endl;
#endif
  Assert(abs(c.Det()-cdet) < EPS*c.colsize()*abs(cdet),"Square CDet");

  Matrix<complex<T> > cinv = c.Inverse();
  Matrix<complex<T> > cid = c*cinv;
  Assert(Norm(cid-tmv::Eye<T>(4)) < EPS*Norm(c)*Norm(cinv),
      "Square CInverse");

  Matrix<complex<T> > ctc = c.Adjoint() * c;
  Assert(Norm(c.InverseATA()-ctc.Inverse()) < EPS*Norm(ctc)*Norm(ctc.Inverse()),
      "Square CInverseATA");

  Vector<complex<T> > e(4);
  e = b*complex<T>(1,2);
  e(1) += complex<T>(-1,5);
  e(2) -= complex<T>(-1,5);

  // test real / complex
  Vector<complex<T> > y = b/c;
  Vector<complex<T> > b3 = c*y;
  Assert(Norm(b3-b) < EPS*Norm(c)*Norm(b),"Square b/c");
  y = b%c;
  b3 = y*c;
  Assert(Norm(b3-b) < EPS*Norm(c)*Norm(b),"Square b%c");

  // test complex / real
  y = e/m;
  b3 = m*y;
  Assert(Norm(b3-e) < EPS*Norm(m)*Norm(e),"Square e/m");
  y = e%m;
  b3 = y*m;
  Assert(Norm(b3-e) < EPS*Norm(m)*Norm(e),"Square e%m");

  // test complex / complex
  y = e/c;
  b3 = c*y;
  Assert(Norm(b3-e) < EPS*Norm(c)*Norm(e),"Square e/c");
  y = e%c;
  b3 = y*c;
  Assert(Norm(b3-e) < EPS*Norm(c)*Norm(e),"Square e%c");

  // test really big one
  complex<T> p[30] = { 3,6,1,6,8,3,3,3,34,25, 76,4,67,52,3,2,1,2,4,6, 57,4,24,7,2,1,33,64,23,9};
  complex<T> q[30] = { 12,3,5,34,52,4,234,243,42,648, 71,5,4,35,5,5,45,3,52,5, 36,32,2,53,63,5,2,43,12,1};
  complex<T> r[30] = { 142,3,51,2,27,42,4,23,42,14, 14,24,2,82,4,24,6,1,6,72, 62,42,13,265,325,37,3,52,32,13};
  Vector<complex<T> > P(30,p);
  Vector<complex<T> > Q(30,q);
  Vector<complex<T> > R(30,r);
  Q *= complex<T>(-4,10);
  Matrix<complex<T> > M = tmv::OuterProduct(P,Q);
  M.DivideUsing(dt);
  M.diag() += Vector<complex<T> >(30,complex<T>(100,-200));
  M(23,10) *= 12.;
  M(12,1) -= 1000.;
  M(6,2) += complex<T>(23,891);
  M(15,0) *= complex<T>(615,12);
  Vector<complex<T> > S = R/M;
  Vector<complex<T> > R2 = M*S;
  //CheckDecomposition<complex<T> >(M,dt);
  Assert(Norm(R2-R) < EPS*Norm(M)*Norm(R),"Square R/M");
  S = R%M;
  R2 = S*M;
  Assert(Norm(R2-R) < EPS*Norm(M)*Norm(R),"Square R%M");

  Matrix<T,RowMajor> a1 = m;
  Matrix<T,ColMajor> a2 = m.Transpose();
  a2.row(1) *= T(3);
  a2.col(2) -= Vector<T>(4,4.);
  Matrix<complex<T>,RowMajor> c1 = a1 * complex<T>(1,2);
  Matrix<complex<T>,ColMajor> c2 = a2 * complex<T>(-3,4);
  c1.diag().AddToAll(complex<T>(3,1));
  c2.diag().AddToAll(complex<T>(-5,8));
  c1.row(3).AddToAll(complex<T>(1,-6));
  c2.row(0).AddToAll(complex<T>(-2,-11));

  TestMatrixDivArith<T>(dt,a1,a2,c1,c2,"Square"); 

  Matrix<T,RowMajor> a3(7,4);
  for(int i=0;i<7;++i) for(int j=0;j<4;++j) a3(i,j) = 1-3*i+2*j;
  Matrix<T,ColMajor> a4 = a3.Transpose();
  a3.SubMatrix(2,6,0,4) += a1;
  a4.SubMatrix(0,4,1,5) -= a2;

  Matrix<complex<T>,RowMajor> c3 = a3*complex<T>(1,2);
  Matrix<complex<T>,ColMajor> c4 = c3.Adjoint();
  c3.SubMatrix(2,6,0,4) += c1;
  c4.SubMatrix(0,4,1,5) -= c2;
  c3.col(1) *= complex<T>(2,1);
  c3.row(2).AddToAll(complex<T>(-7,2));
  c4.col(3) *= complex<T>(-1,3);
  c4.row(0).AddToAll(complex<T>(1,9));

  TestMatrixDivArith<T>(dt,a1,a3,c1,c3,"Square/NonSquare");
  TestMatrixDivArith<T>(dt,a1,a4,c1,c4,"Square/NonSquare");

  Matrix<T> a5(4,0,1);
  Matrix<T> a6(0,4,1);
  Matrix<complex<T> > c5 = a5;
  Matrix<complex<T> > c6 = a6;

  TestMatrixDivArith<T>(dt,a1,a5,c1,c5,"Square/Degenerate");
  TestMatrixDivArith<T>(dt,a1,a6,c1,c6,"Square/Degenerate");

  cout<<"Square Matrix<"<<tmv::Type(T())<<"> Division using ";
  cout<<tmv::Text(dt)<<" passed all tests\n";
}

template <class T> void TestNonSquareDiv(tmv::DivType dt)
{
  Matrix<T> m(6,4);
  m.DivideUsing(dt);
  for(int i=0;i<6;++i) for(int j=0;j<4;++j) m(i,j) = 2.+4*i-5*j;
  m(0,0) = 14.;
  m(1,0) = -2.;
  m(2,0) = 7.;
  m(3,0) = -10.;
  m(2,2) = 30.;

  Vector<T> x(4);
  x(0) = 2;
  x(1) = -10;
  x(2) = 5;
  x(3) = -5;

  Vector<T> b = m * x;
  m.SetDiv();
  Vector<T> x2 = b/m;
  CheckDecomposition<T>(m,dt);
  Assert(Norm(x2-x) < EPS*Norm(m)*Norm(b),"NonSquare exact b/m");

  Vector<T> b2 = x%m;
  x2 = b2*m;
  Assert(Norm(x2-x) < EPS*Norm(m)*Norm(x),"NonSquare x%m");

  Matrix<T> mtm = m.Transpose()*m;
#ifdef SHOWACC
  cout<<"mtm.det = "<<mtm.Det()<<endl;
  cout<<"m.det^2 = "<<m.Det()*m.Det()<<endl;
#endif
  Assert(abs(m.Det()*m.Det()-mtm.Det()) < EPS*abs(mtm.Det()),"NonSquare Det");

  b(0) += 100.;
  x = b/m;
  b2 = m*x;
  T refnorm = Norm(b2-b);
  Vector<T> dx = T(sqrt(EPS))*Norm(x)*tmv::BasisVector<T>(4,0);
  b2 = m*(x+dx);
  Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (1)");
  b2 = m*(x-dx);
  Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (2)");
  dx = T(sqrt(EPS))*Norm(x)*tmv::BasisVector<T>(4,1);
  b2 = m*(x+dx);
  Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (3)");
  b2 = m*(x-dx);
  Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (4)");
  dx = T(sqrt(EPS))*Norm(x)*tmv::BasisVector<T>(4,2);
  b2 = m*(x+dx);
  Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (5)");
  b2 = m*(x-dx);
  Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (6)");
  dx = T(sqrt(EPS))*Norm(x)*tmv::BasisVector<T>(4,3);
  b2 = m*(x+dx);
  Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (7)");
  b2 = m*(x-dx);
  Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (8)");

  Matrix<T> minv = m.Inverse();
  Matrix<T> id = minv*m;
  Assert(Norm(id-tmv::Eye<T>(4)) < EPS*Norm(m),"NonSquare Inverse");
  Matrix<T> nonid = m*minv;
  Assert(Norm(nonid-nonid.Transpose()) < EPS*Norm(m),
      "NonSquare Pseudo-Inverse");

  Assert(Norm(m.InverseATA()-mtm.Inverse()) < EPS*Norm(mtm),"NonSquare InverseATA");

  Matrix<complex<T> > c(6,4);
  c.DivideUsing(dt);
  c = m;
  c(2,3) += complex<T>(2,3);
  c(1,0) *= complex<T>(0,2);
  c.col(1) *= complex<T>(-1,3);
  c.row(3) += Vector<complex<T> >(4,complex<T>(1,9));

  Vector<complex<T> > y(4);
  y(0) = complex<T>(2,9);
  y(1) = complex<T>(-10,4);
  y(2) = complex<T>(5,-1);
  y(3) = complex<T>(-5,-2);

  Vector<complex<T> > e = c * y;
  Vector<complex<T> > y2 = e/c;
  Assert(Norm(y2-y) < EPS*Norm(c)*Norm(y),"NonSquare exact e/c");

  Vector<complex<T> > e2 = y%c;
  y2 = e2*c;
  Assert(Norm(y2-y) < EPS*Norm(c)*Norm(y),"NonSquare e%c");

  Matrix<complex<T> > ctc = c.Adjoint()*c;
#ifdef SHOWACC
  cout<<"|c.det|^2 = "<<norm(c.Det())<<endl;
  cout<<"ctc.det = "<<ctc.Det()<<endl;
#endif
  Assert(abs(norm(c.Det())-ctc.Det()) < EPS*abs(ctc.Det()),"NonSquare CDet");

  Matrix<complex<T> > cinv = c.Inverse();
  Matrix<complex<T> > cid = cinv*c;
  Assert(Norm(cid-tmv::Eye<T>(4)) < EPS*Norm(c),"NonSquare CInverse");
  Matrix<complex<T> > cnonid = c*cinv;
  Assert(Norm(cnonid-cnonid.Adjoint()) < EPS*Norm(c),
      "NonSquare CPseudo-Inverse");

  Assert(Norm(c.InverseATA()-ctc.Inverse()) < EPS*Norm(ctc),"NonSquare CInverseATA");

  // Test short matrix (M < N)
  Matrix<T,ColMajor> ms = m.Transpose();
  ms.DivideUsing(dt);

  b = x * ms;
  x2 = b%ms;
  CheckDecomposition<T>(ms,dt);
  Assert(Norm(x2-x) < EPS*Norm(ms)*Norm(x),"NonSquare exact b%ms");

  b2 = x/ms;
  x2 = ms*b2;
  Assert(Norm(x2-x) < EPS*Norm(ms)*Norm(x),"NonSquare x/ms");

  // Test really long matrix
  Matrix<complex<T> > a(30,10);
  a.DivideUsing(dt);
  for(int i=0;i<30;++i) for(int j=0;j<10;++j) a(i,j) = 7.-13.*i+11.*j;
  a.SubMatrix(0,10,0,10) += complex<T>(30.,20.);
  a.SubMatrix(10,20,0,10) -= complex<T>(50.,-123.);
  a.SubMatrix(20,30,0,10) += complex<T>(10.,-75.);
  a.SubMatrix(1,10,1,10) += complex<T>(99.,100.);
  a.SubMatrix(2,10,2,10) -= complex<T>(51.,37.);

  Vector<complex<T> > s(10);
  for(int i=0;i<10;++i) s(i) = i+2;

  Vector<complex<T> > t = a * s;
  Vector<complex<T> > s2 = t/a;
  Assert(Norm(s2-s) < EPS*Norm(a)*Norm(s),"NonSquare t/a");

  Vector<complex<T> > t2 = s%a;
  s2 = t2*a;
  Assert(Norm(s2-s) < EPS*Norm(a)*Norm(s),"NonSquare t%a");

  // Test QR Downdate:
  Matrix<complex<T> > q = a;
  Matrix<complex<T> > r(10,10);
  QR_Decompose(q.View(),r.View());
  QR_DownDate(r.View(),a.row(29));
  tmv::ConstMatrixView<complex<T> > ax1 = a.SubMatrix(0,29,0,10);
  Matrix<complex<T> > rtr = r.Adjoint()*r;
  Matrix<complex<T> > ata1 = ax1.Adjoint()*ax1;
  Assert(Norm(rtr-ata1) < EPS*Norm(ata1),"QR_DownDate (R1)");
  QR_DownDate(r.View(),a.row(28));
  tmv::ConstMatrixView<complex<T> > ax2 = a.SubMatrix(0,28,0,10);
  rtr = r.Adjoint()*r;
  Matrix<complex<T> > ata2 = ax2.Adjoint()*ax2;
  Assert(Norm(rtr-ata2) < EPS*Norm(ata2),"QR_DownDate (R2)");
  q = a;
  QR_Decompose(q.View(),r.View());
  QR_DownDate(q.View(),r.View(),a.row(29));
  Assert(Norm(r.Adjoint()*r - ata1) < EPS*Norm(ata1),"QR_DownDate (QR) rtr");
  Assert(Norm(q.SubMatrix(0,29,0,10)*r-ax1) < EPS*Norm(ax1),"QR_DownDate (QR) qr");
  QR_DownDate(q.SubMatrix(0,29,0,10),r.View(),a.row(28));
  Assert(Norm(r.Adjoint()*r - ata2) < EPS*Norm(ata2),"QR_DownDate (QR) rtr");
  Assert(Norm(q.SubMatrix(0,28,0,10)*r-ax2) < EPS*Norm(ax2),"QR_DownDate (QR) qr");

  // Test with some identical eigenvalues.
  // First make an arbitrary unitary matrix:
  q = a;
  q.DivideUsing(dt);
  QR_Decompose(q.View(),r.View());
  r.Zero();
  r(0,0) = 1;
  r(1,1) = 5;
  r(2,2) = 1;
  r(3,3) = 1;
  r(4,4) = -3;
  r(5,5) = 7;
  r(6,6) = 1;
  r(7,7) = 1;
  r(8,8) = -2;
  r(9,9) = 1;
  q = q*r;

  t = q * s;
  s2 = t/q;
  //CheckDecomposition<complex<T> >(q,dt);
  Assert(Norm(s2-s) < EPS*Norm(q)*Norm(s),"NonSquare t/q");

  t2 = s%q;
  s2 = t2*q;
  Assert(Norm(s2-s) < EPS*Norm(q)*Norm(s),"NonSquare t%q");

  Matrix<T> a1 = m;
  Matrix<T> a2 = m.Transpose() * m;
  Matrix<T> a3 = m * m.Transpose();
  a2.row(1) *= T(3);
  a2.col(2).AddToAll(-4);
  a3.row(5) *= T(7);
  a3.col(3).AddToAll(7);
  Matrix<complex<T> > c1 = a1 * complex<T>(1,2);
  Matrix<complex<T> > c2 = a2 * complex<T>(-3,4);
  Matrix<complex<T> > c3 = a3 * complex<T>(-4,8);

  TestMatrixDivArith<T>(dt,a1,a2,c1,c2,"NonSquare/Square"); 
  TestMatrixDivArith<T>(dt,a1,a3,c1,c3,"NonSquare/Square"); 

  Matrix<T,RowMajor> a4(7,4);
  for(int i=0;i<7;++i) for(int j=0;j<4;++j) a4(i,j) = 1-3*i+2*j;
  Matrix<T,ColMajor> a5 = a4.Transpose();
  a4.SubMatrix(0,6,0,4) += a1;
  a5.SubMatrix(0,4,1,5) -= a2;
  Matrix<complex<T>,RowMajor> c4 = a4*complex<T>(1,2);
  Matrix<complex<T>,ColMajor> c5 = c4.Adjoint();
  c4.SubMatrix(0,6,0,4) += c1;
  c5.SubMatrix(0,4,1,5) -= c2;
  c4.col(1) *= complex<T>(2,1);
  c4.row(2).AddToAll(complex<T>(-7,2));
  c5.col(3) *= complex<T>(-1,3);
  c5.row(0).AddToAll(complex<T>(1,9));

  Matrix<T,RowMajor> a6(9,6);
  for(int i=0;i<9;++i) for(int j=0;j<6;++j) a6(i,j) = 5+2*i-2*j;
  Matrix<T,ColMajor> a7 = a6.Transpose();
  a6.SubMatrix(2,8,1,5) += a1;
  a7.SubMatrix(0,6,4,8) -= T(2)*a1;
  Matrix<complex<T>,RowMajor> c6 = a6*complex<T>(1,2);
  Matrix<complex<T>,ColMajor> c7 = c6.Adjoint();
  c6.SubMatrix(2,8,1,5) += c1;
  c7.SubMatrix(0,6,4,8) -= T(2)*c1;
  c6.col(1) *= complex<T>(2,1);
  c6.row(5).AddToAll(complex<T>(-7,2));
  c7.col(7) *= complex<T>(-1,3);
  c7.row(4).AddToAll(complex<T>(1,9));

  TestMatrixDivArith<T>(dt,a1,a4,c1,c4,"NonSquare/NonSquare");
  TestMatrixDivArith<T>(dt,a1,a5,c1,c5,"NonSquare/NonSquare");
  TestMatrixDivArith<T>(dt,a1,a6,c1,c6,"NonSquare/NonSquare");
  TestMatrixDivArith<T>(dt,a1,a7,c1,c7,"NonSquare/NonSquare");

  Matrix<T> a8(4,0,1);
  Matrix<T> a9(0,4,1);
  Matrix<T> a10(6,0,1);
  Matrix<T> a11(0,6,1);
  Matrix<complex<T> > c8 = a8;
  Matrix<complex<T> > c9 = a9;
  Matrix<complex<T> > c10 = a10;
  Matrix<complex<T> > c11 = a11;

  TestMatrixDivArith<T>(dt,a1,a8,c1,c8,"NonSquare/Degenerate");
  TestMatrixDivArith<T>(dt,a1,a9,c1,c9,"NonSquare/Degenerate");
  TestMatrixDivArith<T>(dt,a1,a10,c1,c10,"NonSquare/Degenerate");
  TestMatrixDivArith<T>(dt,a1,a11,c1,c11,"NonSquare/Degenerate");

  cout<<"NonSquare Matrix<"<<tmv::Type(T())<<"> Division using ";
  cout<<tmv::Text(dt)<<" passed all tests\n";
}

template <class T> void TestSingularDiv(tmv::DivType dt)
{
  Matrix<T> m(4,4);
  m.DivideUsing(dt);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) m(i,j) = 2.+4*i-5*j;
  m(2,2) = 30.;

  Vector<T> x(4);
  x(0) = 2;
  x(1) = -10;
  x(2) = 5;
  x(3) = -5;

  Vector<T> b = m * x;
  Vector<T> x2 = b/m;
  CheckDecomposition<T>(m,dt);

  Vector<T> b2 = m*x2;
  Assert(Norm(b2-b) < EPS*Norm(m)*Norm(x),"Singular exact b/m");

  b = x * m;
  x2 = b%m;
  b2 = x2*m;
  Assert(Norm(b2-b) < EPS*Norm(m)*Norm(x),"Singular exact b%m");

  b(0) += 100.;
  x = b/m;
  b2 = m*x;
  T norm = Norm(b2-b);
  // Shouldn't need the 3 here, but float and g++ fail without it.
  // icc with T=double doesn't need the 3.
  Vector<T> dx = T(3*sqrt(EPS))*tmv::BasisVector<T>(4,0);
  b2 = m*(x+dx);
  Assert(Norm(b2-b) >= norm,"Singular Least Squares b/m (1)");
  b2 = m*(x-dx);
  Assert(Norm(b2-b) >= norm,"Singular Least Squares b/m (2)");
  dx = T(3*sqrt(EPS))*tmv::BasisVector<T>(4,1);
  b2 = m*(x+dx);
  Assert(Norm(b2-b) >= norm,"Singular Least Squares b/m (3)");
  b2 = m*(x-dx);
  Assert(Norm(b2-b) >= norm,"Singular Least Squares b/m (4)");
  dx = T(3*sqrt(EPS))*tmv::BasisVector<T>(4,2);
  b2 = m*(x+dx);
  Assert(Norm(b2-b) >= norm,"Singular Least Squares b/m (5)");
  b2 = m*(x-dx);
  Assert(Norm(b2-b) >= norm,"Singular Least Squares b/m (6)");
  dx = T(3*sqrt(EPS))*tmv::BasisVector<T>(4,3);
  b2 = m*(x+dx);
  Assert(Norm(b2-b) >= norm,"Singular Least Squares b/m (7)");
  b2 = m*(x-dx);
  Assert(Norm(b2-b) >= norm,"Singular Least Squares b/m (8)");

  Matrix<T> minv = m.Inverse();

#ifdef SHOWACC
  cout<<"m = "<<m<<endl;
  cout<<"minv = "<<minv<<endl;
  cout<<"minv*m = "<<minv*m<<endl;
  cout<<"m*minv = "<<m*minv<<endl;
  cout<<"m*minv*m = "<<m*minv*m<<endl;
  cout<<"minv*m*minv = "<<minv*m*minv<<endl;
  cout<<"m*minv-(m*minv)T = "<<m*minv-Transpose(m*minv)<<endl;
  cout<<"minv*m-(minv*m)T = "<<minv*m-Transpose(minv*m)<<endl;
#endif

  Assert(Norm(m*minv*m - m) < EPS*Norm(m),"Singular Inverse M*X*M != M");
  Assert(Norm(minv*m*minv - minv) < EPS*Norm(m),"Singular Inverse X*M*X != X");
  Assert(Norm((m*minv)-(m*minv).Transpose()) < EPS*Norm(m),
      "Singular Inverse M*X != (M*X)T");
  if (dt != tmv::QRP) { // QRP doesn't get this right.
    Assert(Norm((minv*m)-(minv*m).Transpose()) < EPS*Norm(m),
	"Singular Inverse X*M != (X*M)T");
  }

  // Try big one with many singular values.
  Matrix<complex<T> > mm(30,30);
  mm.DivideUsing(dt);
  for(int i=0;i<30;++i) for(int j=0;j<30;++j) mm(i,j) = 4.-17.*i+23.*j;
  mm(20,20) += complex<T>(200.,-999.);
  mm(12,12) += complex<T>(500.,-104.);
  mm(7,7) += complex<T>(300.,123.);
  mm(28,28) += complex<T>(700.,231.);
  mm(24,24) += complex<T>(400.,-120.);

  Vector<complex<T> > xx(30);
  for(int i=0;i<30;++i) xx(i) = 10.+i;
  Vector<complex<T> > bb = mm*xx;
  Vector<complex<T> > xx2 = bb/mm;
  Vector<complex<T> > bb2 = mm*xx2;
#ifdef SHOWACC
  CheckDecomposition<complex<T> >(mm,dt);
  cout<<"Norm(bb2-bb) = "<<Norm(bb2-bb);
  cout<<", EPS*Norm(mm)*Norm(xx) = "<<EPS*Norm(mm)*Norm(xx)<<endl;
#endif
  Assert(Norm(bb2-bb) < 2*EPS*Norm(mm)*Norm(xx),"Singular exact bb/mm");

  bb = xx * mm;
  xx2 = bb%mm;
  bb2 = xx2*mm;
#ifdef SHOWACC
  cout<<"Norm(bb2-bb) = "<<Norm(bb2-bb);
  cout<<", EPS*Norm(mm)*Norm(xx) = "<<EPS*Norm(mm)*Norm(xx)<<endl;
#endif
  Assert(Norm(bb2-bb) < 2*EPS*Norm(mm)*Norm(xx),"Singular exact bb%mm");

  cout<<"Singular Matrix<"<<tmv::Type(T())<<"> Division using ";
  cout<<tmv::Text(dt)<<" passed all tests\n";
}

template <class T> void TestAllMatrixDiv()
{
  TestNonSquareDiv<T>(tmv::QR);
  TestSquareDiv<T>(tmv::LU);
  TestSquareDiv<T>(tmv::QR);
  TestSquareDiv<T>(tmv::QRP);
  TestSquareDiv<T>(tmv::SV);
  TestNonSquareDiv<T>(tmv::QR);
  TestNonSquareDiv<T>(tmv::QRP);
  TestNonSquareDiv<T>(tmv::SV);
  TestSingularDiv<T>(tmv::QRP);
  TestSingularDiv<T>(tmv::SV);
}

template void TestAllMatrixDiv<double>();
#ifndef NOFLOAT
template void TestAllMatrixDiv<float>();
#endif
#ifdef LONGDOUBLE
template void TestAllMatrixDiv<long double>();
#endif
