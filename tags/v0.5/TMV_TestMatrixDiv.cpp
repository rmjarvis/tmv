//#define SHOWACC
//#define SHOWTESTS
//#define SHOWCHECK

#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_TestMatrixDiv.h"

using tmv::Matrix;
using tmv::Vector;

template <class T> void TestSquareDiv(tmv::DivType dt)
{
  Matrix<T> m(4,4);
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

  m.DivideUsing(dt);
  m.SetDiv();
#ifdef SHOWCHECK
  Assert(m.CheckDecomp(&cout),"CheckDecomp");
#else
  Assert(m.CheckDecomp(),"CheckDecomp");
#endif
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
#ifdef SHOWACC
  cout<<"minv = "<<minv<<endl;
  cout<<"m*minv = "<<id<<endl;
  cout<<"minv*m = "<<minv*m<<endl;
  cout<<"Norm(id-I) = "<<Norm(id-tmv::Eye<T>(4))<<endl;
#endif
  Assert(Norm(id-tmv::Eye<T>(4)) < EPS*Norm(minv)*Norm(m),
      "Square Inverse");

  Matrix<T> mtm = m.Adjoint() * m;
#ifdef SHOWACC
  cout<<"mtm = "<<mtm<<endl;
  cout<<"mtm.inv = "<<mtm.Inverse()<<endl;
  cout<<"m.invata = "<<m.InverseATA()<<endl;
  cout<<"minv*minvt = "<<minv*minv.Adjoint()<<endl;
  cout<<"Norm(diff) = "<<Norm(m.InverseATA()-mtm.Inverse())<<endl;
#endif
  Assert(Norm(m.InverseATA()-mtm.Inverse()) < EPS*Norm(mtm)*Norm(mtm.Inverse()),
	"Square InverseATA");

  T mdet = 28800.;
#ifdef SHOWACC
  cout<<"abs(det-mdet) = "<<abs(m.Det()-mdet);
  cout<<"  EPS*abs(mdet) = "<<EPS*m.colsize()*abs(mdet)<<endl;
#endif
  Assert(abs(m.Det()-mdet) < EPS*m.colsize()*abs(mdet),"Square Det");

  Matrix<complex<T> > c(4,4);
  c = m;
  c(2,3) += complex<T>(2,3);
  c(1,0) *= complex<T>(0,2);
  c.col(1) *= complex<T>(-1,3);
  c.row(3) += Vector<complex<T> >(4,complex<T>(1,9));

  c.DivideUsing(dt);
  c.SetDiv();
#ifdef SHOWCHECK
  Assert(c.CheckDecomp(&cout),"CheckDecomp");
#else
  Assert(c.CheckDecomp(),"CheckDecomp");
#endif
  complex<T> cdet(-103604,101272);
#ifdef SHOWACC
  cout<<"cdet = "<<cdet<<endl;
  cout<<"C.Det = "<<c.Det()<<endl;
  cout<<"abs(det-cdet) = "<<abs(c.Det()-cdet);
  cout<<"  EPS*abs(cdet) = "<<EPS*c.colsize()*abs(cdet)<<endl;
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
  T p[100] = { 
    3,6,1,6,8,3,3,3,34,25, 76,4,67,52,3,2,1,2,4,6, 57,4,24,7,2,1,33,64,23,9,
    3,6,1,6,8,3,3,3,34,25, 76,4,67,52,3,2,1,2,4,6, 57,4,24,7,2,1,33,64,23,9,
    3,6,1,6,8,3,3,3,34,25, 76,4,67,52,3,2,1,2,4,6, 57,4,24,7,2,1,33,64,23,9,
    3,6,1,6,8,3,3,3,34,25 
  };
  T q[100] = { 
    12,3,5,34,52,4,4,23,42,68, 71,5,4,5,5,5,45,3,52,5, 36,32,2,53,6,5,2,43,1,1,
    12,3,5,34,52,4,4,23,42,68, 71,5,4,5,5,5,45,3,52,5, 36,32,2,53,6,5,2,43,1,1,
    12,3,5,34,52,4,4,23,42,68, 71,5,4,5,5,5,45,3,52,5, 36,32,2,53,6,5,2,43,1,1,
    12,3,5,34,52,4,4,23,42,68
  };
  T r[100] = {
    42,3,51,2,7,42,4,3,42,14, 14,24,2,82,4,24,6,1,6,7, 2,42,1,2,35,7,3,5,32,13,
    42,3,51,2,7,42,4,3,42,14, 14,24,2,82,4,24,6,1,6,7, 2,42,1,2,35,7,3,5,32,13,
    42,3,51,2,7,42,4,3,42,14, 14,24,2,82,4,24,6,1,6,7, 2,42,1,2,35,7,3,5,32,13,
    42,3,51,2,7,42,4,3,42,14
  };
  Vector<T> P(100,p);
  Vector<T> Q(100,q);
  Vector<T> R(100,r);
  Matrix<T,ColMajor> M = P ^ Q;
  Matrix<complex<T>,ColMajor> CM = P ^ (complex<T>(-4,10)*Q);
  M.diag().AddToAll(T(215));
  CM.diag().AddToAll(complex<T>(103,-53));
  M.row(23) *= T(12);
  M(12,1) -= T(142);
  CM(6,2) += complex<T>(23,89);
  CM.col(15) *= complex<T>(61,12);
  M.SubVector(65,5,1,3,29) *= T(2);
  M.SubVector(98,12,-1,2,18).AddToAll(T(197));
  CM.SubVector(53,0,1,3,31) *= complex<T>(2,-1);
  CM.SubVector(88,18,-1,2,23).AddToAll(complex<T>(197,174));
  M.DivideUsing(dt);
  CM.DivideUsing(dt);
  Vector<T> S = R/M;
  Vector<T> R2 = M*S;
#ifdef SHOWCHECK
  Assert(M.CheckDecomp(&cout),"CheckDecomp");
#else
  Assert(M.CheckDecomp(),"CheckDecomp");
#endif
#ifdef SHOWACC
  cout<<"R/M Norm(R2-R) = "<<Norm(R2-R)<<endl;
  cout<<"EPS*Norm(M)*Norm(R) = "<<EPS*Norm(M)*Norm(R)<<endl;
#endif
  Assert(Norm(R2-R) < 10*EPS*Norm(M)*Norm(R),"Square R/M");
  S = R%M;
  R2 = S*M;
#ifdef SHOWACC
  cout<<"R%M Norm(R2-R) = "<<Norm(R2-R)<<endl;
  cout<<"EPS*Norm(M)*Norm(R) = "<<EPS*Norm(M)*Norm(R)<<endl;
#endif
  Assert(Norm(R2-R) < 10*EPS*Norm(M)*Norm(R),"Square R%M");
  Vector<complex<T> > CS = R/CM;
  Vector<complex<T> > CR2 = CM*CS;
#ifdef SHOWCHECK
  Assert(CM.CheckDecomp(&cout),"CheckDecomp");
#else
  Assert(CM.CheckDecomp(),"CheckDecomp");
#endif
#ifdef SHOWACC
  cerr<<"R = "<<R<<endl;
  cerr<<"CS = "<<CS<<endl;
  cerr<<"CR2 = "<<CR2<<endl;
  cout<<"R/CM Norm(CR2-R) = "<<Norm(CR2-R)<<endl;
  cout<<"EPS*Norm(CM)*Norm(R) = "<<EPS*Norm(CM)*Norm(R)<<endl;
#endif
  Assert(Norm(CR2-R) < 10*EPS*Norm(CM)*Norm(R),"Square R/CM");
  CS = R%CM;
  CR2 = CS*CM;
#ifdef SHOWACC
  cout<<"R%CM Norm(CR2-R) = "<<Norm(CR2-R)<<endl;
  cout<<"EPS*Norm(CM)*Norm(R) = "<<EPS*Norm(CM)*Norm(R)<<endl;
#endif
  Assert(Norm(CR2-R) < 10*EPS*Norm(CM)*Norm(R),"Square R%CM");
  Vector<complex<T> > CR = complex<T>(3,-4)*R;
  CS = CR/CM;
  CR2 = CM*CS;
#ifdef SHOWACC
  cout<<"CR/CM Norm(CR2-CR) = "<<Norm(CR2-CR)<<endl;
  cout<<"EPS*Norm(CM)*Norm(CR) = "<<EPS*Norm(CM)*Norm(CR)<<endl;
#endif
  Assert(Norm(CR2-CR) < 10*EPS*Norm(CM)*Norm(CR),"Square CR/CM");
  CS = CR%CM;
  CR2 = CS*CM;
#ifdef SHOWACC
  cout<<"CR%CM Norm(CR2-CR) = "<<Norm(CR2-CR)<<endl;
  cout<<"EPS*Norm(CM)*Norm(CR) = "<<EPS*Norm(CM)*Norm(CR)<<endl;
#endif
  Assert(Norm(CR2-CR) < 10*EPS*Norm(CM)*Norm(CR),"Square CR%CM");
  CS = CR/M;
  CR2 = M*CS;
#ifdef SHOWACC
  cout<<"CR/M Norm(CR2-CR) = "<<Norm(CR2-CR)<<endl;
  cout<<"EPS*Norm(M)*Norm(CR) = "<<EPS*Norm(M)*Norm(CR)<<endl;
#endif
  Assert(Norm(CR2-CR) < 10*EPS*Norm(M)*Norm(CR),"Square CR/M");
  CS = CR%M;
  CR2 = CS*M;
#ifdef SHOWACC
  cout<<"CR%M Norm(CR2-CR) = "<<Norm(CR2-CR)<<endl;
  cout<<"EPS*Norm(M)*Norm(CR) = "<<EPS*Norm(M)*Norm(CR)<<endl;
#endif
  Assert(Norm(CR2-CR) < 10*EPS*Norm(M)*Norm(CR),"Square CR%M");

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

template <class T> void TestSquareSV()
{
  Matrix<T> m(4,4);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) m(i,j) = T(2.+4*i-5*j);
  m(0,0) = T(14);
  m(1,0) = T(-2);
  m(2,0) = T(7);
  m(3,0) = T(-10);
  m.diag() *= T(30);

  m.DivideUsing(tmv::SV);
  m.SetDiv();
  Matrix<T> U = m.SVD().GetU();
  Vector<RealType(T)> S = m.SVD().GetS();
  Matrix<T> V = m.SVD().GetV();

  Matrix<T> m1 = m;
  m1.DivideUsing(tmv::SVF);
  m1.SetDiv();
  Matrix<T> m2 = m;
  m2.DivideUsing(tmv::SVS);
  m2.SetDiv();
  Matrix<T> m3 = m;
  m3.DivideUsing(tmv::SVU);
  m3.SetDiv();
  Matrix<T> m4 = m;
  m4.DivideUsing(tmv::SVV);
  m4.SetDiv();

  Assert(Norm(m1.SVFD().GetS()-S) < EPS*(Norm(m)),"Square SVF S");
  Assert(Norm(m2.SVFD().GetS()-S) < EPS*(Norm(m)),"Square SVS S");
  Assert(Norm(m3.SVFD().GetS()-S) < EPS*(Norm(m)),"Square SVU S");
  Assert(Norm(m4.SVFD().GetS()-S) < EPS*(Norm(m)),"Square SVV S");

  Assert(Norm(m1.SVFD().GetU()-U) < EPS*(Norm(m)),"Square SVF U");
  Assert(Norm(m3.SVFD().GetU()-U) < EPS*(Norm(m)),"Square SVU U");

  Assert(Norm(m1.SVFD().GetV()-V) < EPS*(Norm(m)),"Square SVF V");
  Assert(Norm(m4.SVFD().GetV()-V) < EPS*(Norm(m)),"Square SVV V");

  Matrix<complex<T> > c(4,4);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) 
    c(i,j) = complex<T>(2.+4*i-5*j,3.-i);
  c(0,0) *= T(14);
  c(1,0) *= T(-2);
  c(2,0) *= T(7);
  c(3,0) *= T(-10);
  c.diag() *= T(30);

  c.DivideUsing(tmv::SV);
  c.SetDiv();
  Matrix<complex<T> > cU = c.SVD().GetU();
  Vector<T> cS = c.SVD().GetS();
  Matrix<complex<T> > cV = c.SVD().GetV();

  Matrix<complex<T> > c1 = c;
  c1.DivideUsing(tmv::SVF);
  c1.SetDiv();
  Matrix<complex<T> > c2 = c;
  c2.DivideUsing(tmv::SVS);
  c2.SetDiv();
  Matrix<complex<T> > c3 = c;
  c3.DivideUsing(tmv::SVU);
  c3.SetDiv();
  Matrix<complex<T> > c4 = c;
  c4.DivideUsing(tmv::SVV);
  c4.SetDiv();

  Assert(Norm(c1.SVFD().GetS()-cS) < EPS*(Norm(c)),"Square C SVF S");
  Assert(Norm(c2.SVFD().GetS()-cS) < EPS*(Norm(c)),"Square C SVS S");
  Assert(Norm(c3.SVFD().GetS()-cS) < EPS*(Norm(c)),"Square C SVU S");
  Assert(Norm(c4.SVFD().GetS()-cS) < EPS*(Norm(c)),"Square C SVV S");

  Assert(Norm(c1.SVFD().GetU()-cU) < EPS*(Norm(c)),"Square C SVF U");
  Assert(Norm(c3.SVFD().GetU()-cU) < EPS*(Norm(c)),"Square C SVU U");

  Assert(Norm(c1.SVFD().GetV()-cV) < EPS*(Norm(c)),"Square C SVF V");
  Assert(Norm(c4.SVFD().GetV()-cV) < EPS*(Norm(c)),"Square C SVV V");
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

  m.SetDiv();
#ifdef SHOWCHECK
  Assert(m.CheckDecomp(&cout),"CheckDecomp");
#else
  Assert(m.CheckDecomp(),"CheckDecomp");
#endif

  Vector<T> b = m * x;
  Vector<T> x2 = b/m;
  Assert(Norm(x2-x) < EPS*Norm(m)*Norm(b),"NonSquare exact b/m");

  Vector<T> b2 = x%m;
  x2 = b2*m;
  Assert(Norm(x2-x) < EPS*Norm(m)*Norm(x),"NonSquare x%m");

  Matrix<T> mtm = m.Transpose()*m;
#ifdef SHOWACC
  cout<<"mtm.det = "<<mtm.Det()<<endl;
  cout<<"m.det^2 = "<<m.Det()*m.Det()<<endl;
#endif
  Assert(abs(m.Det()*m.Det()-mtm.Det()) < EPS*m.colsize()*abs(mtm.Det()),"NonSquare Det");

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
  Matrix<T> nonid = m*minv;
#ifdef SHOWACC
  cout<<"minv = "<<minv<<endl;
  cout<<"minv*m = "<<id<<endl;
  cout<<"m*minv = "<<nonid<<endl;
  cout<<"Norm(id-I) = "<<Norm(id-tmv::Eye<T>(4))<<endl;
#endif
  Assert(Norm(id-tmv::Eye<T>(4)) < EPS*Norm(m)*Norm(minv),"NonSquare Inverse");
  Assert(Norm(nonid-nonid.Transpose()) < EPS*Norm(m)*Norm(minv),
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

#ifdef SHOWCHECK
  Assert(c.CheckDecomp(&cout),"CheckDecomp");
#else
  Assert(c.CheckDecomp(),"CheckDecomp");
#endif

  Assert(Norm(y2-y) < EPS*Norm(c)*Norm(y),"NonSquare exact e/c");

  Vector<complex<T> > e2 = y%c;
  y2 = e2*c;
  Assert(Norm(y2-y) < EPS*Norm(c)*Norm(y),"NonSquare e%c");

  Matrix<complex<T> > ctc = c.Adjoint()*c;
#ifdef SHOWACC
  cout<<"|c.det|^2 = "<<norm(c.Det())<<endl;
  cout<<"ctc.det = "<<ctc.Det()<<endl;
  cout<<"abs(|c.det|^2-ctc.det) = "<<abs(norm(c.Det())-ctc.Det())<<endl;
  cerr<<"eps*ctc.det = "<<EPS*c.colsize()*abs(ctc.Det())<<endl;
#endif
  Assert(abs(norm(c.Det())-ctc.Det()) < EPS*c.colsize()*abs(ctc.Det()),"NonSquare CDet");

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
#ifdef SHOWCHECK
  Assert(ms.CheckDecomp(&cout),"CheckDecomp");
#else
  Assert(ms.CheckDecomp(),"CheckDecomp");
#endif
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
#ifdef SHOWCHECK
  Assert(q.CheckDecomp(&cout),"CheckDecomp");
#else
  Assert(q.CheckDecomp(),"CheckDecomp");
#endif
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

template <class T> void TestNonSquareSV()
{
  Matrix<T> m(6,4);
  for(int i=0;i<6;++i) for(int j=0;j<4;++j) m(i,j) = T(2.+4*i-5*j);
  m(0,0) = T(14);
  m(1,0) = T(-2);
  m(2,0) = T(7);
  m(3,0) = T(-10);
  m.diag() *= T(30);

  m.DivideUsing(tmv::SV);
  m.SetDiv();
  Matrix<T> U = m.SVD().GetU();
  Vector<RealType(T)> S = m.SVD().GetS();
  Matrix<T> V = m.SVD().GetV();

  Matrix<T> m1 = m;
  m1.DivideUsing(tmv::SVF);
  m1.SetDiv();
  Matrix<T> m2 = m;
  m2.DivideUsing(tmv::SVS);
  m2.SetDiv();
  Matrix<T> m3 = m;
  m3.DivideUsing(tmv::SVU);
  m3.SetDiv();
  Matrix<T> m4 = m;
  m4.DivideUsing(tmv::SVV);
  m4.SetDiv();

  Assert(Norm(m1.SVFD().GetS()-S) < EPS*(Norm(m)),"NonSquare SVF S");
  Assert(Norm(m2.SVFD().GetS()-S) < EPS*(Norm(m)),"NonSquare SVS S");
  Assert(Norm(m3.SVFD().GetS()-S) < EPS*(Norm(m)),"NonSquare SVU S");
  Assert(Norm(m4.SVFD().GetS()-S) < EPS*(Norm(m)),"NonSquare SVV S");

  Assert(Norm(m1.SVFD().GetU()-U) < EPS*(Norm(m)),"NonSquare SVF U");
  Assert(Norm(m3.SVFD().GetU()-U) < EPS*(Norm(m)),"NonSquare SVU U");

  Assert(Norm(m1.SVFD().GetV()-V) < EPS*(Norm(m)),"NonSquare SVF V");
  Assert(Norm(m4.SVFD().GetV()-V) < EPS*(Norm(m)),"NonSquare SVV V");

  Matrix<T> a(15,4);
  for(int i=0;i<15;++i) for(int j=0;j<4;++j) a(i,j) = T(2.+4*i-5*j);
  a(0,0) = T(14);
  a(1,0) = T(-2);
  a(2,0) = T(7);
  a(3,0) = T(-10);
  a.diag() *= T(30);

  a.DivideUsing(tmv::SV);
  a.SetDiv();
  Matrix<T> aU = a.SVD().GetU();
  Vector<RealType(T)> aS = a.SVD().GetS();
  Matrix<T> aV = a.SVD().GetV();

  Matrix<T> a1 = a;
  a1.DivideUsing(tmv::SVF);
  a1.SetDiv();
  Matrix<T> a2 = a;
  a2.DivideUsing(tmv::SVS);
  a2.SetDiv();
  Matrix<T> a3 = a;
  a3.DivideUsing(tmv::SVU);
  a3.SetDiv();
  Matrix<T> a4 = a;
  a4.DivideUsing(tmv::SVV);
  a4.SetDiv();

  Assert(Norm(a1.SVFD().GetS()-aS) < EPS*(Norm(a)),"Tall NonSquare SVF S");
  Assert(Norm(a2.SVFD().GetS()-aS) < EPS*(Norm(a)),"Tall NonSquare SVS S");
  Assert(Norm(a3.SVFD().GetS()-aS) < EPS*(Norm(a)),"Tall NonSquare SVU S");
  Assert(Norm(a4.SVFD().GetS()-aS) < EPS*(Norm(a)),"Tall NonSquare SVV S");

  Assert(Norm(a1.SVFD().GetU()-aU) < EPS*(Norm(a)),"Tall NonSquare SVF U");
  Assert(Norm(a3.SVFD().GetU()-aU) < EPS*(Norm(a)),"Tall NonSquare SVU U");

  Assert(Norm(a1.SVFD().GetV()-aV) < EPS*(Norm(a)),"Tall NonSquare SVF V");
  Assert(Norm(a4.SVFD().GetV()-aV) < EPS*(Norm(a)),"Tall NonSquare SVV V");

  Matrix<complex<T> > c(6,4);
  for(int i=0;i<6;++i) for(int j=0;j<4;++j)
    c(i,j) = complex<T>(2.+4*i-5*j,3.-i);
  c(0,0) = T(14);
  c(1,0) = T(-2);
  c(2,0) = T(7);
  c(3,0) = T(-10);
  c.diag() *= T(30);

  c.DivideUsing(tmv::SV);
  c.SetDiv();
  Matrix<complex<T> > cU = c.SVD().GetU();
  Vector<T> cS = c.SVD().GetS();
  Matrix<complex<T> > cV = c.SVD().GetV();

  Matrix<complex<T> > c1 = c;
  c1.DivideUsing(tmv::SVF);
  c1.SetDiv();
  Matrix<complex<T> > c2 = c;
  c2.DivideUsing(tmv::SVS);
  c2.SetDiv();
  Matrix<complex<T> > c3 = c;
  c3.DivideUsing(tmv::SVU);
  c3.SetDiv();
  Matrix<complex<T> > c4 = c;
  c4.DivideUsing(tmv::SVV);
  c4.SetDiv();

  Assert(Norm(c1.SVFD().GetS()-cS) < EPS*(Norm(c)),"NonSquare C SVF S");
  Assert(Norm(c2.SVFD().GetS()-cS) < EPS*(Norm(c)),"NonSquare C SVS S");
  Assert(Norm(c3.SVFD().GetS()-cS) < EPS*(Norm(c)),"NonSquare C SVU S");
  Assert(Norm(c4.SVFD().GetS()-cS) < EPS*(Norm(c)),"NonSquare C SVV S");

  Assert(Norm(c1.SVFD().GetU()-cU) < EPS*(Norm(c)),"NonSquare C SVF U");
  Assert(Norm(c3.SVFD().GetU()-cU) < EPS*(Norm(c)),"NonSquare C SVU U");

  Assert(Norm(c1.SVFD().GetV()-cV) < EPS*(Norm(c)),"NonSquare C SVF V");
  Assert(Norm(c4.SVFD().GetV()-cV) < EPS*(Norm(c)),"NonSquare C SVV V");

  Matrix<complex<T> > d(15,4);
  for(int i=0;i<15;++i) for(int j=0;j<4;++j)
    d(i,j) = complex<T>(2.+4*i-5*j,3.-i);
  d(0,0) = T(14);
  d(1,0) = T(-2);
  d(2,0) = T(7);
  d(3,0) = T(-10);
  d.diag() *= T(30);

  d.DivideUsing(tmv::SV);
  d.SetDiv();
  Matrix<complex<T> > dU = d.SVD().GetU();
  Vector<T> dS = d.SVD().GetS();
  Matrix<complex<T> > dV = d.SVD().GetV();

  Matrix<complex<T> > d1 = d;
  d1.DivideUsing(tmv::SVF);
  d1.SetDiv();
  Matrix<complex<T> > d2 = d;
  d2.DivideUsing(tmv::SVS);
  d2.SetDiv();
  Matrix<complex<T> > d3 = d;
  d3.DivideUsing(tmv::SVU);
  d3.SetDiv();
  Matrix<complex<T> > d4 = d;
  d4.DivideUsing(tmv::SVV);
  d4.SetDiv();

  Assert(Norm(d1.SVFD().GetS()-dS) < EPS*(Norm(d)),"Tall NonSquare C SVF S");
  Assert(Norm(d2.SVFD().GetS()-dS) < EPS*(Norm(d)),"Tall NonSquare C SVS S");
  Assert(Norm(d3.SVFD().GetS()-dS) < EPS*(Norm(d)),"Tall NonSquare C SVU S");
  Assert(Norm(d4.SVFD().GetS()-dS) < EPS*(Norm(d)),"Tall NonSquare C SVV S");

  Assert(Norm(d1.SVFD().GetU()-dU) < EPS*(Norm(d)),"Tall NonSquare C SVF U");
  Assert(Norm(d3.SVFD().GetU()-dU) < EPS*(Norm(d)),"Tall NonSquare C SVU U");

  Assert(Norm(d1.SVFD().GetV()-dV) < EPS*(Norm(d)),"Tall NonSquare C SVF V");
  Assert(Norm(d4.SVFD().GetV()-dV) < EPS*(Norm(d)),"Tall NonSquare C SVV V");
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

  m.SetDiv();
#ifdef SHOWCHECK
  Assert(m.CheckDecomp(&cout),"CheckDecomp");
#else
  Assert(m.CheckDecomp(),"CheckDecomp");
#endif

  Vector<T> b = m * x;
  Vector<T> x2 = b/m;
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
#ifdef SHOWCHECK
  Assert(mm.CheckDecomp(&cout),"CheckDecomp");
#else
  Assert(mm.CheckDecomp(),"CheckDecomp");
#endif
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

template <class T> void TestSingularSV()
{
  Matrix<T> m(4,4);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) m(i,j) = 2.+4*i-5*j;
  m(2,2) = 30.;

  m.DivideUsing(tmv::SV);
  m.SetDiv();
  Matrix<T> U = m.SVD().GetU();
  Vector<RealType(T)> S = m.SVD().GetS();
  Matrix<T> V = m.SVD().GetV();

  Matrix<T> m1 = m;
  m1.DivideUsing(tmv::SVF);
  m1.SetDiv();
  Matrix<T> m2 = m;
  m2.DivideUsing(tmv::SVS);
  m2.SetDiv();
  Matrix<T> m3 = m;
  m3.DivideUsing(tmv::SVU);
  m3.SetDiv();
  Matrix<T> m4 = m;
  m4.DivideUsing(tmv::SVV);
  m4.SetDiv();

  Assert(Norm(m1.SVFD().GetS()-S) < EPS*(Norm(m)),"Singular SVF S");
  Assert(Norm(m2.SVFD().GetS()-S) < EPS*(Norm(m)),"Singular SVS S");
  Assert(Norm(m3.SVFD().GetS()-S) < EPS*(Norm(m)),"Singular SVU S");
  Assert(Norm(m4.SVFD().GetS()-S) < EPS*(Norm(m)),"Singular SVV S");

  Matrix<complex<T> > c(4,4);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) 
    c(i,j) = complex<T>(2.+4*i-5*j,3.-i);
  m(2,2) += 30.;

  c.DivideUsing(tmv::SV);
  c.SetDiv();
  Matrix<complex<T> > cU = c.SVD().GetU();
  Vector<T> cS = c.SVD().GetS();
  Matrix<complex<T> > cV = c.SVD().GetV();

  Matrix<complex<T> > c1 = c;
  c1.DivideUsing(tmv::SVF);
  c1.SetDiv();
  Matrix<complex<T> > c2 = c;
  c2.DivideUsing(tmv::SVS);
  c2.SetDiv();
  Matrix<complex<T> > c3 = c;
  c3.DivideUsing(tmv::SVU);
  c3.SetDiv();
  Matrix<complex<T> > c4 = c;
  c4.DivideUsing(tmv::SVV);
  c4.SetDiv();

  Assert(Norm(c1.SVFD().GetS()-cS) < EPS*(Norm(c)),"Singular C SVF S");
  Assert(Norm(c2.SVFD().GetS()-cS) < EPS*(Norm(c)),"Singular C SVS S");
  Assert(Norm(c3.SVFD().GetS()-cS) < EPS*(Norm(c)),"Singular C SVU S");
  Assert(Norm(c4.SVFD().GetS()-cS) < EPS*(Norm(c)),"Singular C SVV S");
}

template <class T> void TestAllMatrixDiv()
{
  TestSquareDiv<T>(tmv::LU);
  TestSquareDiv<T>(tmv::QR);
  TestSquareDiv<T>(tmv::QRP);
  TestSquareDiv<T>(tmv::SVF);
  TestSquareSV<T>();
  TestSquareDiv<T>(tmv::SV);
  TestNonSquareDiv<T>(tmv::QR);
  TestNonSquareDiv<T>(tmv::QRP);
  TestNonSquareDiv<T>(tmv::SVF);
  TestNonSquareSV<T>();
  TestNonSquareDiv<T>(tmv::SV);
  TestSingularDiv<T>(tmv::QRP);
  TestSingularDiv<T>(tmv::SVF);
  TestSingularSV<T>();
  TestSingularDiv<T>(tmv::SV);
}

template void TestAllMatrixDiv<double>();
#ifndef NOFLOAT
template void TestAllMatrixDiv<float>();
#endif
#ifdef LONGDOUBLE
template void TestAllMatrixDiv<long double>();
#endif
