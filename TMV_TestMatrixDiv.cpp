#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_TestMatrixDiv.h"
#include "TMV_Tri.h"

using tmv::Matrix;
using tmv::Vector;
using tmv::UpperTriMatrix;
using tmv::StorageType;
using tmv::RowMajor;
using tmv::ColMajor;
using tmv::DivType;

template <class T, StorageType stor> void TestSquareDiv(DivType dt)
{
  Matrix<T,stor> m(4,4);
  m.SaveDiv();

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
  if (showdiv) {
    Assert(m.CheckDecomp(&cout),"CheckDecomp");
  } else {
    Assert(m.CheckDecomp(),"CheckDecomp");
  }
  Vector<T> x = b/m;
  Vector<T> b2 = m*x;
  if (showacc) {
    cout<<"b = "<<b<<endl;
    cout<<"x = b/m = "<<x<<endl;
    cout<<"b2 = "<<b2<<endl;
    cout<<"Norm(b-b2) = "<<Norm(b-b2)<<endl;
  }
  Assert(Norm(b2-b) < EPS*(Norm(m)*Norm(b)),"Square b/m");

  x = b%m;
  b2 = x*m;
  if (showacc) {
    cout<<"b = "<<b<<endl;
    cout<<"x = b%m = "<<x<<endl;
    cout<<"b2 = "<<b2<<endl;
    cout<<"Norm(b-b2) = "<<Norm(b-b2)<<endl;
  }
  Assert(Norm(b2-b) < EPS*(Norm(m)*Norm(b)),"Square b%m");

  Matrix<T> minv = m.Inverse();
  Matrix<T> id = m*minv;
  if (showacc) {
    cout<<"minv = "<<minv<<endl;
    cout<<"m*minv = "<<id<<endl;
    cout<<"minv*m = "<<minv*m<<endl;
    cout<<"Norm(id-I) = "<<Norm(id-T(1))<<endl;
  }
  Assert(Norm(id-T(1)) < EPS*Norm(minv)*Norm(m),"Square Inverse");

  Matrix<T> mtm = m.Adjoint() * m;
  Matrix<T> mata(4,4);
  m.InverseATA(mata);
  if (showacc) {
    cout<<"mtm = "<<mtm<<endl;
    cout<<"mtm.inv = "<<mtm.Inverse()<<endl;
    cout<<"m.invata = "<<mata<<endl;
    cout<<"minv*minvt = "<<minv*minv.Adjoint()<<endl;
    cout<<"Norm(diff) = "<<Norm(mata-mtm.Inverse())<<endl;
  }
  Assert(Norm(mata-mtm.Inverse()) < EPS*Norm(mtm)*Norm(mtm.Inverse()),
	"Square InverseATA");

  T mdet = 28800.;
  if (showacc) {
    cout<<"abs(det-mdet) = "<<abs(m.Det()-mdet);
    cout<<"  EPS*abs(mdet) = "<<EPS*m.colsize()*abs(mdet)<<endl;
  }
  Assert(abs(m.Det()-mdet) < EPS*m.colsize()*abs(mdet),"Square Det");

  Matrix<complex<T>,stor> c(4,4);
  c.SaveDiv();
  c = m;
  c(2,3) += complex<T>(2,3);
  c(1,0) *= complex<T>(0,2);
  c.col(1) *= complex<T>(-1,3);
  c.row(3) += Vector<complex<T> >(4,complex<T>(1,9));

  c.DivideUsing(dt);
  c.SetDiv();
  if (showdiv) {
    Assert(c.CheckDecomp(&cout),"CheckDecomp");
  } else {
    Assert(c.CheckDecomp(),"CheckDecomp");
  }
  complex<T> cdet(-103604,101272);
  if (showacc) {
    cout<<"cdet = "<<cdet<<endl;
    cout<<"C.Det = "<<c.Det()<<endl;
    cout<<"abs(det-cdet) = "<<abs(c.Det()-cdet);
    cout<<"  EPS*abs(cdet) = "<<EPS*c.colsize()*abs(cdet)<<endl;
  }
  Assert(abs(c.Det()-cdet) < EPS*c.colsize()*abs(cdet),"Square CDet");

  Matrix<complex<T> > cinv = c.Inverse();
  Matrix<complex<T> > cid = c*cinv;
  Assert(Norm(cid-T(1)) < EPS*Norm(c)*Norm(cinv),"Square CInverse");

  Matrix<complex<T> > ctc = c.Adjoint() * c;
  Matrix<complex<T> > cata(4,4);
  c.InverseATA(cata);
  Assert(Norm(cata-ctc.Inverse()) < EPS*Norm(ctc)*Norm(ctc.Inverse()),
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
  Matrix<T,stor> M = P ^ Q;
  Matrix<complex<T>,stor> CM = P ^ (complex<T>(-4,10)*Q);
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
  M.SaveDiv();
  CM.SaveDiv();
  Vector<T> S = R/M;
  Vector<T> R2 = M*S;
  Assert(M.CheckDecomp(),"CheckDecomp");
  if (showacc) {
    cout<<"R/M Norm(R2-R) = "<<Norm(R2-R)<<endl;
    cout<<"EPS*Norm(M)*Norm(R) = "<<EPS*Norm(M)*Norm(R)<<endl;
  }
  Assert(Norm(R2-R) < 10*EPS*Norm(M)*Norm(R),"Square R/M");
  S = R%M;
  R2 = S*M;
  if (showacc) {
    cout<<"R%M Norm(R2-R) = "<<Norm(R2-R)<<endl;
    cout<<"EPS*Norm(M)*Norm(R) = "<<EPS*Norm(M)*Norm(R)<<endl;
  }
  Assert(Norm(R2-R) < 10*EPS*Norm(M)*Norm(R),"Square R%M");
  //cerr<<"Passed test\n";
  Vector<complex<T> > CS = R/CM;
  //cerr<<"CS = "<<CS<<endl;
  Vector<complex<T> > CR2 = CM*CS;
  //cerr<<"CR2 = "<<CR2<<endl;
  Assert(CM.CheckDecomp(),"CheckDecomp");
  if (showacc) {
    cout<<"R/CM Norm(CR2-R) = "<<Norm(CR2-R)<<endl;
    cout<<"EPS*Norm(CM)*Norm(R) = "<<EPS*Norm(CM)*Norm(R)<<endl;
  }
  Assert(Norm(CR2-R) < 10*EPS*Norm(CM)*Norm(R),"Square R/CM");
  CS = R%CM;
  CR2 = CS*CM;
  if (showacc) {
    cout<<"R%CM Norm(CR2-R) = "<<Norm(CR2-R)<<endl;
    cout<<"EPS*Norm(CM)*Norm(R) = "<<EPS*Norm(CM)*Norm(R)<<endl;
  }
  Assert(Norm(CR2-R) < 10*EPS*Norm(CM)*Norm(R),"Square R%CM");
  Vector<complex<T> > CR = complex<T>(3,-4)*R;
  CS = CR/CM;
  CR2 = CM*CS;
  if (showacc) {
    cout<<"CR/CM Norm(CR2-CR) = "<<Norm(CR2-CR)<<endl;
    cout<<"EPS*Norm(CM)*Norm(CR) = "<<EPS*Norm(CM)*Norm(CR)<<endl;
  }
  Assert(Norm(CR2-CR) < 10*EPS*Norm(CM)*Norm(CR),"Square CR/CM");
  CS = CR%CM;
  CR2 = CS*CM;
  if (showacc) {
    cout<<"CR%CM Norm(CR2-CR) = "<<Norm(CR2-CR)<<endl;
    cout<<"EPS*Norm(CM)*Norm(CR) = "<<EPS*Norm(CM)*Norm(CR)<<endl;
  }
  Assert(Norm(CR2-CR) < 10*EPS*Norm(CM)*Norm(CR),"Square CR%CM");
  CS = CR/M;
  CR2 = M*CS;
  if (showacc) {
    cout<<"CR/M Norm(CR2-CR) = "<<Norm(CR2-CR)<<endl;
    cout<<"EPS*Norm(M)*Norm(CR) = "<<EPS*Norm(M)*Norm(CR)<<endl;
  }
  Assert(Norm(CR2-CR) < 10*EPS*Norm(M)*Norm(CR),"Square CR/M");
  CS = CR%M;
  CR2 = CS*M;
  if (showacc) {
    cout<<"CR%M Norm(CR2-CR) = "<<Norm(CR2-CR)<<endl;
    cout<<"EPS*Norm(M)*Norm(CR) = "<<EPS*Norm(M)*Norm(CR)<<endl;
  }
  Assert(Norm(CR2-CR) < 10*EPS*Norm(M)*Norm(CR),"Square CR%M");

  Matrix<T,stor> a1 = m;
  Matrix<T,stor> a2 = m.Transpose();
  a2.row(1) *= T(3);
  a2.col(2) -= Vector<T>(4,4.);

  Matrix<complex<T>,stor> c1 = a1 * complex<T>(1,2);
  Matrix<complex<T>,stor> c2 = a2 * complex<T>(-3,4);
  c1.diag().AddToAll(complex<T>(3,1));
  c2.diag().AddToAll(complex<T>(-5,8));
  c1.row(3).AddToAll(complex<T>(1,-6));
  c2.row(0).AddToAll(complex<T>(-2,-11));

  TestMatrixDivArith<T>(dt,a1.View(),a2.View(),c1.View(),c2.View(),"Square"); 
#ifdef XTEST
  Matrix<T,stor,tmv::FortranStyle> a1f = a1;
  Matrix<T,stor,tmv::FortranStyle> a2f = a2;
  Matrix<complex<T>,stor,tmv::FortranStyle> c1f = c1;
  Matrix<complex<T>,stor,tmv::FortranStyle> c2f = c2;
  TestMatrixDivArith<T>(dt,a1f.View(),a2.View(),c1f.View(),c2.View(),
      "Square"); 
  TestMatrixDivArith<T>(dt,a1.View(),a2f.View(),c1.View(),c2f.View(),
      "Square"); 
  TestMatrixDivArith<T>(dt,a1f.View(),a2f.View(),c1f.View(),c2f.View(),
      "Square"); 
#endif

  Matrix<T,stor> a3(7,4);
  for(int i=0;i<7;++i) for(int j=0;j<4;++j) a3(i,j) = 1-3*i+2*j;
  Matrix<T,stor> a4 = a3.Transpose();
  a3.SubMatrix(2,6,0,4) += a1;
  a4.SubMatrix(0,4,1,5) -= a2;

  Matrix<complex<T>,stor> c3 = a3*complex<T>(1,2);
  Matrix<complex<T>,stor> c4 = c3.Adjoint();
  c3.SubMatrix(2,6,0,4) += c1;
  c4.SubMatrix(0,4,1,5) -= c2;
  c3.col(1) *= complex<T>(2,1);
  c3.row(2).AddToAll(complex<T>(-7,2));
  c4.col(3) *= complex<T>(-1,3);
  c4.row(0).AddToAll(complex<T>(1,9));

  TestMatrixDivArith<T>(dt,a1.View(),a3.View(),c1.View(),c3.View(),
      "Square/NonSquare");
  TestMatrixDivArith<T>(dt,a1.View(),a4.View(),c1.View(),c4.View(),
      "Square/NonSquare");

#ifdef XTEST
  Matrix<T,stor> a5(4,0,1);
  Matrix<T,stor> a6(0,4,1);
  Matrix<complex<T>,stor> c5 = a5;
  Matrix<complex<T>,stor> c6 = a6;

  TestMatrixDivArith<T>(dt,a1.View(),a5.View(),c1.View(),c5.View(),
      "Square/Degenerate");
  TestMatrixDivArith<T>(dt,a1.View(),a6.View(),c1.View(),c6.View(),
      "Square/Degenerate");
#endif

  if (stor == ColMajor) {
#ifdef XTEST
    TestSquareDiv<T,RowMajor>(dt);
  } else {
#endif
    cout<<"Square Matrix<"<<tmv::Type(T())<<"> Division using ";
    cout<<tmv::Text(dt)<<" passed all tests\n";
  }
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

  Matrix<T> m2 = m;
  m2.DivideUsing(tmv::SVS);
  m2.SetDiv();
  Matrix<T> m3 = m;
  m3.DivideUsing(tmv::SVU);
  m3.SetDiv();
  Matrix<T> m4 = m;
  m4.DivideUsing(tmv::SVV);
  m4.SetDiv();

  Assert(Norm(m2.SVD().GetS()-S) < EPS*(Norm(m)),"Square SVS S");
  Assert(Norm(m3.SVD().GetS()-S) < EPS*(Norm(m)),"Square SVU S");
  Assert(Norm(m4.SVD().GetS()-S) < EPS*(Norm(m)),"Square SVV S");

  Assert(Norm(m3.SVD().GetU()-U) < EPS*(Norm(m)),"Square SVU U");

  Assert(Norm(m4.SVD().GetV()-V) < EPS*(Norm(m)),"Square SVV V");

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

  Matrix<complex<T> > c2 = c;
  c2.DivideUsing(tmv::SVS);
  c2.SetDiv();
  Matrix<complex<T> > c3 = c;
  c3.DivideUsing(tmv::SVU);
  c3.SetDiv();
  Matrix<complex<T> > c4 = c;
  c4.DivideUsing(tmv::SVV);
  c4.SetDiv();

  Assert(Norm(c2.SVD().GetS()-cS) < EPS*(Norm(c)),"Square C SVS S");
  Assert(Norm(c3.SVD().GetS()-cS) < EPS*(Norm(c)),"Square C SVU S");
  Assert(Norm(c4.SVD().GetS()-cS) < EPS*(Norm(c)),"Square C SVV S");

  Assert(Norm(c3.SVD().GetU()-cU) < EPS*(Norm(c)),"Square C SVU U");

  Assert(Norm(c4.SVD().GetV()-cV) < EPS*(Norm(c)),"Square C SVV V");
}

template <class T, StorageType stor> void TestNonSquareDiv(DivType dt)
{
  Matrix<T,stor> m(6,4);
  m.SaveDiv();
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
  if (showdiv) {
    Assert(m.CheckDecomp(&cout),"CheckDecomp");
  } else {
    Assert(m.CheckDecomp(),"CheckDecomp");
  }

  Vector<T> b = m * x;
  Vector<T> x2 = b/m;
  Assert(Norm(x2-x) < EPS*Norm(m)*Norm(b),"NonSquare exact b/m");

  Vector<T> b2 = x%m;
  x2 = b2*m;
  Assert(Norm(x2-x) < EPS*Norm(m)*Norm(x),"NonSquare x%m");

  Matrix<T> mtm = m.Transpose()*m;
  if (showacc) {
    cout<<"mtm.det = "<<mtm.Det()<<endl;
    cout<<"m.det^2 = "<<m.Det()*m.Det()<<endl;
  }
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
  if (showacc) {
    cout<<"minv = "<<minv<<endl;
    cout<<"minv*m = "<<id<<endl;
    cout<<"m*minv = "<<nonid<<endl;
    cout<<"Norm(id-I) = "<<Norm(id-T(1))<<endl;
  }
  Assert(Norm(id-T(1)) < EPS*Norm(m)*Norm(minv),"NonSquare Inverse");
  Assert(Norm(nonid-nonid.Transpose()) < EPS*Norm(m)*Norm(minv),
      "NonSquare Pseudo-Inverse");

  Matrix<T> mata(4,4);
  m.InverseATA(mata);
  Assert(Norm(mata-mtm.Inverse()) < EPS*Norm(mtm),"NonSquare InverseATA");

  Matrix<complex<T>,stor> c(6,4);
  c.SaveDiv();
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

  if (showdiv) {
    Assert(c.CheckDecomp(&cout),"CheckDecomp");
  } else {
    Assert(c.CheckDecomp(),"CheckDecomp");
  }

  Assert(Norm(y2-y) < EPS*Norm(c)*Norm(y),"NonSquare exact e/c");

  Vector<complex<T> > e2 = y%c;
  y2 = e2*c;
  Assert(Norm(y2-y) < EPS*Norm(c)*Norm(y),"NonSquare e%c");

  Matrix<complex<T> > ctc = c.Adjoint()*c;
  if (showacc) {
    cout<<"|c.det|^2 = "<<norm(c.Det())<<endl;
    cout<<"ctc.det = "<<ctc.Det()<<endl;
    cout<<"abs(|c.det|^2-ctc.det) = "<<abs(norm(c.Det())-ctc.Det())<<endl;
    cout<<"eps*ctc.det = "<<EPS*c.colsize()*abs(ctc.Det())<<endl;
  }
  Assert(abs(norm(c.Det())-ctc.Det()) < EPS*c.colsize()*abs(ctc.Det()),"NonSquare CDet");

  Matrix<complex<T> > cinv = c.Inverse();
  Matrix<complex<T> > cid = cinv*c;
  Assert(Norm(cid-T(1)) < EPS*Norm(c),"NonSquare CInverse");
  Matrix<complex<T> > cnonid = c*cinv;
  Assert(Norm(cnonid-cnonid.Adjoint()) < EPS*Norm(c),
      "NonSquare CPseudo-Inverse");

  Matrix<complex<T> > cata(4,4);
  c.InverseATA(cata);
  Assert(Norm(cata-ctc.Inverse()) < EPS*Norm(ctc),"NonSquare CInverseATA");

  // Test short matrix (M < N)
  Matrix<T,stor> ms = m.Transpose();
  ms.DivideUsing(dt);
  ms.SaveDiv();

  b = x * ms;
  x2 = b%ms;
  if (showdiv) {
    Assert(ms.CheckDecomp(&cout),"CheckDecomp");
  } else {
    Assert(ms.CheckDecomp(),"CheckDecomp");
  }
  Assert(Norm(x2-x) < EPS*Norm(ms)*Norm(x),"NonSquare exact b%ms");

  b2 = x/ms;
  x2 = ms*b2;
  Assert(Norm(x2-x) < EPS*Norm(ms)*Norm(x),"NonSquare x/ms");

  // Test really long matrix
  Matrix<complex<T>,stor> a(30,10);
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

  // Test QR Update/Downdate:

  Matrix<complex<T>,stor> q30 = a;
  UpperTriMatrix<complex<T>,tmv::NonUnitDiag,stor> r30(10,10);
  QR_Decompose(q30.View(),r30.View());
  Assert(Norm(q30*r30-a) < EPS*a.NormSq(),"QR_Decompose");
  Assert(Norm(r30.Adjoint()*r30-a.Adjoint()*a) < EPS*r30.NormSq(),
      "QR_Decompose (RtR)");

  Matrix<complex<T>,stor> q10 = a.Rows(0,10);
  UpperTriMatrix<complex<T>,tmv::NonUnitDiag,stor> r10(10,10);
  QR_Decompose(q10.View(),r10.View());
  UpperTriMatrix<complex<T>,tmv::NonUnitDiag,stor> r = r10;
  Matrix<complex<T>,stor> a1030 = a.Rows(10,30);
  QR_Update(r.View(),a1030.View());
  Assert(Norm(r.Adjoint()*r-r30.Adjoint()*r30) < EPS*a.NormSq(),
      "QR_Update");
  r = r30;

  a1030 = a.Rows(10,30);
  QR_Downdate(r.View(),a1030.View());
  Assert(Norm(r.Adjoint()*r-r10.Adjoint()*r10) < EPS*a.NormSq(),
      "QR_Downdate");
  r = r10;

  Matrix<complex<T>,stor> a1020 = a.Rows(10,20);
  Matrix<complex<T>,stor> a2030 = a.Rows(20,30);
  QR_Update(r.View(),a1020.View());
  QR_Update(r.View(),a2030.View());
  Assert(Norm(r.Adjoint()*r-r30.Adjoint()*r30) < EPS*a.NormSq(),
      "QR_Update (double)");
  r = r30;

  a1020 = a.Rows(10,20);
  a2030 = a.Rows(20,30);
  QR_Downdate(r.View(),a1020.View());
  QR_Downdate(r.View(),a2030.View());
  Assert(Norm(r.Adjoint()*r-r10.Adjoint()*r10) < EPS*a.NormSq(),
      "QR_Downdate (double)");
  r = r10;

  Matrix<complex<T>,stor> q29 = a.Rows(0,29);
  UpperTriMatrix<complex<T>,tmv::NonUnitDiag,stor> r29(10,10);
  QR_Decompose(q29.View(),r29.View());
  r = r30;
  Vector<complex<T> > a29 = a.row(29);
  QR_Downdate(r.View(),a29.View());
  Assert(Norm(r.Adjoint()*r-r29.Adjoint()*r29) < EPS*a.NormSq(),
      "QR_Downdate (single row)");
  r = r29;

  a29 = a.row(29);
  QR_Update(r.View(),a29.View());
  Assert(Norm(r.Adjoint()*r-r30.Adjoint()*r30) < EPS*a.NormSq(),
      "QR_Downdate (single row)");

  // Test with some identical eigenvalues.
  // First make an arbitrary unitary matrix:
  Matrix<complex<T>,stor> q = a;
  q.DivideUsing(dt);
  q.SaveDiv();
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
  Assert(q.CheckDecomp(),"CheckDecomp");
  Assert(Norm(s2-s) < EPS*Norm(q)*Norm(s),"NonSquare t/q");

  t2 = s%q;
  s2 = t2*q;
  Assert(Norm(s2-s) < EPS*Norm(q)*Norm(s),"NonSquare t%q");

  Matrix<T,stor> a1 = m;
  Matrix<T,stor> a2 = m.Transpose() * m;
  Matrix<T,stor> a3 = m * m.Transpose();
  a2.row(1) *= T(3);
  a2.col(2).AddToAll(-4);
  a3.row(5) *= T(7);
  a3.col(3).AddToAll(7);
  Matrix<complex<T>,stor> c1 = a1 * complex<T>(1,2);
  Matrix<complex<T>,stor> c2 = a2 * complex<T>(-3,4);
  Matrix<complex<T>,stor> c3 = a3 * complex<T>(-4,8);

  TestMatrixDivArith<T>(dt,a1.View(),a2.View(),c1.View(),c2.View(),
      "NonSquare/Square"); 
  TestMatrixDivArith<T>(dt,a1.View(),a3.View(),c1.View(),c3.View(),
      "NonSquare/Square"); 

  Matrix<T,stor> a4(7,4);
  for(int i=0;i<7;++i) for(int j=0;j<4;++j) a4(i,j) = 1-3*i+2*j;
  Matrix<T,stor> a5 = a4.Transpose();
  a4.SubMatrix(0,6,0,4) += a1;
  a5.SubMatrix(0,4,1,5) -= a2;
  Matrix<complex<T>,stor> c4 = a4*complex<T>(1,2);
  Matrix<complex<T>,stor> c5 = c4.Adjoint();
  c4.SubMatrix(0,6,0,4) += c1;
  c5.SubMatrix(0,4,1,5) -= c2;
  c4.col(1) *= complex<T>(2,1);
  c4.row(2).AddToAll(complex<T>(-7,2));
  c5.col(3) *= complex<T>(-1,3);
  c5.row(0).AddToAll(complex<T>(1,9));

  Matrix<T,stor> a6(9,6);
  for(int i=0;i<9;++i) for(int j=0;j<6;++j) a6(i,j) = 5+2*i-2*j;
  Matrix<T,stor> a7 = a6.Transpose();
  a6.SubMatrix(2,8,1,5) += a1;
  a7.SubMatrix(0,6,4,8) -= T(2)*a1;
  Matrix<complex<T>,stor> c6 = a6*complex<T>(1,2);
  Matrix<complex<T>,stor> c7 = c6.Adjoint();
  c6.SubMatrix(2,8,1,5) += c1;
  c7.SubMatrix(0,6,4,8) -= T(2)*c1;
  c6.col(1) *= complex<T>(2,1);
  c6.row(5).AddToAll(complex<T>(-7,2));
  c7.col(7) *= complex<T>(-1,3);
  c7.row(4).AddToAll(complex<T>(1,9));

  TestMatrixDivArith<T>(dt,a1.View(),a4.View(),c1.View(),c4.View(),
      "NonSquare/NonSquare");
  TestMatrixDivArith<T>(dt,a1.View(),a5.View(),c1.View(),c5.View(),
      "NonSquare/NonSquare");
  TestMatrixDivArith<T>(dt,a1.View(),a6.View(),c1.View(),c6.View(),
      "NonSquare/NonSquare");
  TestMatrixDivArith<T>(dt,a1.View(),a7.View(),c1.View(),c7.View(),
      "NonSquare/NonSquare");

#ifdef XTEST
  Matrix<T,stor> a8(4,0,1);
  Matrix<T,stor> a9(0,4,1);
  Matrix<T,stor> a10(6,0,1);
  Matrix<T,stor> a11(0,6,1);
  Matrix<complex<T>,stor> c8 = a8;
  Matrix<complex<T>,stor> c9 = a9;
  Matrix<complex<T>,stor> c10 = a10;
  Matrix<complex<T>,stor> c11 = a11;

  TestMatrixDivArith<T>(dt,a1.View(),a8.View(),c1.View(),c8.View(),
      "NonSquare/Degenerate");
  TestMatrixDivArith<T>(dt,a1.View(),a9.View(),c1.View(),c9.View(),
      "NonSquare/Degenerate");
  TestMatrixDivArith<T>(dt,a1.View(),a10.View(),c1.View(),c10.View(),
      "NonSquare/Degenerate");
  TestMatrixDivArith<T>(dt,a1.View(),a11.View(),c1.View(),c11.View(),
      "NonSquare/Degenerate");
#endif

  if (stor == ColMajor) {
#ifdef XTEST
    TestNonSquareDiv<T,RowMajor>(dt);
  } else {
#endif
    cout<<"NonSquare Matrix<"<<tmv::Type(T())<<"> Division using ";
    cout<<tmv::Text(dt)<<" passed all tests\n";
  }
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

  Matrix<T> m2 = m;
  m2.DivideUsing(tmv::SVS);
  m2.SetDiv();
  Matrix<T> m3 = m;
  m3.DivideUsing(tmv::SVU);
  m3.SetDiv();
  Matrix<T> m4 = m;
  m4.DivideUsing(tmv::SVV);
  m4.SetDiv();

  Assert(Norm(m2.SVD().GetS()-S) < EPS*(Norm(m)),"NonSquare SVS S");
  Assert(Norm(m3.SVD().GetS()-S) < EPS*(Norm(m)),"NonSquare SVU S");
  Assert(Norm(m4.SVD().GetS()-S) < EPS*(Norm(m)),"NonSquare SVV S");

  Assert(Norm(m3.SVD().GetU()-U) < EPS*(Norm(m)),"NonSquare SVU U");

  Assert(Norm(m4.SVD().GetV()-V) < EPS*(Norm(m)),"NonSquare SVV V");

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

  Matrix<T> a2 = a;
  a2.DivideUsing(tmv::SVS);
  a2.SetDiv();
  Matrix<T> a3 = a;
  a3.DivideUsing(tmv::SVU);
  a3.SetDiv();
  Matrix<T> a4 = a;
  a4.DivideUsing(tmv::SVV);
  a4.SetDiv();

  Assert(Norm(a2.SVD().GetS()-aS) < EPS*(Norm(a)),"Tall NonSquare SVS S");
  Assert(Norm(a3.SVD().GetS()-aS) < EPS*(Norm(a)),"Tall NonSquare SVU S");
  Assert(Norm(a4.SVD().GetS()-aS) < EPS*(Norm(a)),"Tall NonSquare SVV S");

  Assert(Norm(a3.SVD().GetU()-aU) < EPS*(Norm(a)),"Tall NonSquare SVU U");

  Assert(Norm(a4.SVD().GetV()-aV) < EPS*(Norm(a)),"Tall NonSquare SVV V");

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

  Matrix<complex<T> > c2 = c;
  c2.DivideUsing(tmv::SVS);
  c2.SetDiv();
  Matrix<complex<T> > c3 = c;
  c3.DivideUsing(tmv::SVU);
  c3.SetDiv();
  Matrix<complex<T> > c4 = c;
  c4.DivideUsing(tmv::SVV);
  c4.SetDiv();

  Assert(Norm(c2.SVD().GetS()-cS) < EPS*(Norm(c)),"NonSquare C SVS S");
  Assert(Norm(c3.SVD().GetS()-cS) < EPS*(Norm(c)),"NonSquare C SVU S");
  Assert(Norm(c4.SVD().GetS()-cS) < EPS*(Norm(c)),"NonSquare C SVV S");

  Assert(Norm(c3.SVD().GetU()-cU) < EPS*(Norm(c)),"NonSquare C SVU U");

  Assert(Norm(c4.SVD().GetV()-cV) < EPS*(Norm(c)),"NonSquare C SVV V");

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

  Matrix<complex<T> > d2 = d;
  d2.DivideUsing(tmv::SVS);
  d2.SetDiv();
  Matrix<complex<T> > d3 = d;
  d3.DivideUsing(tmv::SVU);
  d3.SetDiv();
  Matrix<complex<T> > d4 = d;
  d4.DivideUsing(tmv::SVV);
  d4.SetDiv();

  Assert(Norm(d2.SVD().GetS()-dS) < EPS*(Norm(d)),"Tall NonSquare C SVS S");
  Assert(Norm(d3.SVD().GetS()-dS) < EPS*(Norm(d)),"Tall NonSquare C SVU S");
  Assert(Norm(d4.SVD().GetS()-dS) < EPS*(Norm(d)),"Tall NonSquare C SVV S");

  Assert(Norm(d3.SVD().GetU()-dU) < EPS*(Norm(d)),"Tall NonSquare C SVU U");

  Assert(Norm(d4.SVD().GetV()-dV) < EPS*(Norm(d)),"Tall NonSquare C SVV V");
}

template <class T, StorageType stor> void TestSingularDiv(DivType dt)
{
  Matrix<T,stor> m(4,4);
  m.SaveDiv();
  m.DivideUsing(dt);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) m(i,j) = 2.+4*i-5*j;
  m(2,2) = 30.;

  Vector<T> x(4);
  x(0) = 2;
  x(1) = -10;
  x(2) = 5;
  x(3) = -5;

  m.SetDiv();
  if (showdiv) {
    Assert(m.CheckDecomp(&cout),"CheckDecomp");
  } else {
    Assert(m.CheckDecomp(),"CheckDecomp");
  }

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

  if (showacc) {
    cout<<"m = "<<m<<endl;
    cout<<"minv = "<<minv<<endl;
    cout<<"minv*m = "<<minv*m<<endl;
    cout<<"m*minv = "<<m*minv<<endl;
    cout<<"m*minv*m = "<<m*minv*m<<endl;
    cout<<"minv*m*minv = "<<minv*m*minv<<endl;
    cout<<"m*minv-(m*minv)T = "<<m*minv-Transpose(m*minv)<<endl;
    cout<<"minv*m-(minv*m)T = "<<minv*m-Transpose(minv*m)<<endl;
  }

  Assert(Norm(m*minv*m - m) < EPS*Norm(m),"Singular Inverse M*X*M != M");
  Assert(Norm(minv*m*minv - minv) < EPS*Norm(m),"Singular Inverse X*M*X != X");
  Assert(Norm((m*minv)-(m*minv).Transpose()) < EPS*Norm(m),
      "Singular Inverse M*X != (M*X)T");
  if (dt != tmv::QRP) { // QRP doesn't get this right.
    Assert(Norm((minv*m)-(minv*m).Transpose()) < EPS*Norm(m),
	"Singular Inverse X*M != (X*M)T");
  }

  // Try big one with many singular values.
  Matrix<complex<T>,stor> mm(30,30);
  mm.SaveDiv();
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
  Assert(mm.CheckDecomp(),"CheckDecomp");
  if (showacc) {
    cout<<"Norm(bb2-bb) = "<<Norm(bb2-bb);
    cout<<", EPS*Norm(mm)*Norm(xx) = "<<EPS*Norm(mm)*Norm(xx)<<endl;
  }
  Assert(Norm(bb2-bb) < 30*EPS*Norm(mm)*Norm(xx),"Singular exact bb/mm");

  bb = xx * mm;
  xx2 = bb%mm;
  bb2 = xx2*mm;
  if (showacc) {
    cout<<"Norm(bb2-bb) = "<<Norm(bb2-bb);
    cout<<", EPS*Norm(mm)*Norm(xx) = "<<EPS*Norm(mm)*Norm(xx)<<endl;
  }
  Assert(Norm(bb2-bb) < 30*EPS*Norm(mm)*Norm(xx),"Singular exact bb%mm");

  if (stor == ColMajor) {
#ifdef XTEST
    TestSingularDiv<T,RowMajor>(dt);
  } else {
#endif
    cout<<"Singular Matrix<"<<tmv::Type(T())<<"> Division using ";
    cout<<tmv::Text(dt)<<" passed all tests\n";
  }
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

  Matrix<T> m2 = m;
  m2.DivideUsing(tmv::SVS);
  m2.SetDiv();
  Matrix<T> m3 = m;
  m3.DivideUsing(tmv::SVU);
  m3.SetDiv();
  Matrix<T> m4 = m;
  m4.DivideUsing(tmv::SVV);
  m4.SetDiv();

  Assert(Norm(m2.SVD().GetS()-S) < EPS*(Norm(m)),"Singular SVS S");
  Assert(Norm(m3.SVD().GetS()-S) < EPS*(Norm(m)),"Singular SVU S");
  Assert(Norm(m4.SVD().GetS()-S) < EPS*(Norm(m)),"Singular SVV S");

  Matrix<complex<T> > c(4,4);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) 
    c(i,j) = complex<T>(2.+4*i-5*j,3.-i);
  m(2,2) += 30.;

  c.DivideUsing(tmv::SV);
  c.SetDiv();
  Matrix<complex<T> > cU = c.SVD().GetU();
  Vector<T> cS = c.SVD().GetS();
  Matrix<complex<T> > cV = c.SVD().GetV();

  Matrix<complex<T> > c2 = c;
  c2.DivideUsing(tmv::SVS);
  c2.SetDiv();
  Matrix<complex<T> > c3 = c;
  c3.DivideUsing(tmv::SVU);
  c3.SetDiv();
  Matrix<complex<T> > c4 = c;
  c4.DivideUsing(tmv::SVV);
  c4.SetDiv();

  Assert(Norm(c2.SVD().GetS()-cS) < EPS*(Norm(c)),"Singular C SVS S");
  Assert(Norm(c3.SVD().GetS()-cS) < EPS*(Norm(c)),"Singular C SVU S");
  Assert(Norm(c4.SVD().GetS()-cS) < EPS*(Norm(c)),"Singular C SVV S");
}

template <class T> void TestAllMatrixDiv()
{
  TestSquareDiv<T,ColMajor>(tmv::LU);
  TestSquareDiv<T,ColMajor>(tmv::QR);
  TestSquareDiv<T,ColMajor>(tmv::QRP);
  TestSquareDiv<T,ColMajor>(tmv::SV);
  TestSquareSV<T>();
  TestNonSquareDiv<T,ColMajor>(tmv::QR);
  TestNonSquareDiv<T,ColMajor>(tmv::QRP);
  TestNonSquareDiv<T,ColMajor>(tmv::SV);
  TestNonSquareSV<T>();
  TestSingularDiv<T,ColMajor>(tmv::QRP);
  TestSingularDiv<T,ColMajor>(tmv::SV);
  TestSingularSV<T>();
}

template void TestAllMatrixDiv<double>();
#ifndef NOFLOAT
template void TestAllMatrixDiv<float>();
#endif
#ifdef LONGDOUBLE
template void TestAllMatrixDiv<long double>();
#endif
