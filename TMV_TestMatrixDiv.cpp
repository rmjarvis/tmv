#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_TestMatrixDivArith.h"
#include "TMV_Tri.h"

template <class T, tmv::StorageType stor> inline void TestSquareDiv(tmv::DivType dt)
{
  tmv::Matrix<T,stor> m(4,4);
  m.SaveDiv();

  for(int i=0;i<4;++i) for(int j=0;j<4;++j) m(i,j) = 2.+4*i-5*j;
  m(0,0) = 14.;
  m(1,0) = -2.;
  m(2,0) = 7.;
  m(3,0) = -10.;
  m(2,2) = 30.;

  tmv::Vector<T> b(4);
  b(0) = 2;
  b(1) = -10;
  b(2) = 5;
  b(3) = -5;

  m.DivideUsing(dt);
  m.SetDiv();
  if (showdiv) {
    Assert(m.CheckDecomp(&std::cout),"CheckDecomp");
  } else {
    Assert(m.CheckDecomp(),"CheckDecomp");
  }
  tmv::Vector<T> x = b/m;
  tmv::Vector<T> b2 = m*x;
  if (showacc) {
    std::cout<<"b = "<<b<<std::endl;
    std::cout<<"x = b/m = "<<x<<std::endl;
    std::cout<<"b2 = "<<b2<<std::endl;
    std::cout<<"Norm(b-b2) = "<<Norm(b-b2)<<std::endl;
  }
  Assert(Norm(b2-b) < EPS*(Norm(m)*Norm(b)),"Square b/m");

  x = b%m;
  b2 = x*m;
  if (showacc) {
    std::cout<<"b = "<<b<<std::endl;
    std::cout<<"x = b%m = "<<x<<std::endl;
    std::cout<<"b2 = "<<b2<<std::endl;
    std::cout<<"Norm(b-b2) = "<<Norm(b-b2)<<std::endl;
  }
  Assert(Norm(b2-b) < EPS*(Norm(m)*Norm(b)),"Square b%m");

  tmv::Matrix<T> minv = m.Inverse();
  tmv::Matrix<T> id = m*minv;
  if (showacc) {
    std::cout<<"minv = "<<minv<<std::endl;
    std::cout<<"m*minv = "<<id<<std::endl;
    std::cout<<"minv*m = "<<minv*m<<std::endl;
    std::cout<<"Norm(id-I) = "<<Norm(id-T(1))<<std::endl;
  }
  Assert(Norm(id-T(1)) < EPS*Norm(minv)*Norm(m),"Square Inverse");

  tmv::Matrix<T> mtm = m.Adjoint() * m;
  tmv::Matrix<T> mata(4,4);
  m.InverseATA(mata);
  if (showacc) {
    std::cout<<"mtm = "<<mtm<<std::endl;
    std::cout<<"mtm.inv = "<<mtm.Inverse()<<std::endl;
    std::cout<<"m.invata = "<<mata<<std::endl;
    std::cout<<"minv*minvt = "<<minv*minv.Adjoint()<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm(mata-mtm.Inverse())<<std::endl;
  }
  Assert(Norm(mata-mtm.Inverse()) < EPS*Norm(mtm)*Norm(mtm.Inverse()),
	"Square InverseATA");

  T mdet = 28800.;
  if (showacc) {
    std::cout<<"abs(det-mdet) = "<<std::abs(m.Det()-mdet);
    std::cout<<"  EPS*abs(mdet) = "<<EPS*m.colsize()*std::abs(mdet)<<std::endl;
  }
  Assert(std::abs(m.Det()-mdet) < EPS*m.colsize()*std::abs(mdet),"Square Det");

  tmv::Matrix<std::complex<T>,stor> c(4,4);
  c.SaveDiv();
  c = m;
  c(2,3) += std::complex<T>(2,3);
  c(1,0) *= std::complex<T>(0,2);
  c.col(1) *= std::complex<T>(-1,3);
  c.row(3) += tmv::Vector<std::complex<T> >(4,std::complex<T>(1,9));

  c.DivideUsing(dt);
  c.SetDiv();
  if (showdiv) {
    Assert(c.CheckDecomp(&std::cout),"CheckDecomp");
  } else {
    Assert(c.CheckDecomp(),"CheckDecomp");
  }
  std::complex<T> cdet(-103604,101272);
  if (showacc) {
    std::cout<<"cdet = "<<cdet<<std::endl;
    std::cout<<"C.Det = "<<c.Det()<<std::endl;
    std::cout<<"abs(det-cdet) = "<<std::abs(c.Det()-cdet);
    std::cout<<"  EPS*abs(cdet) = "<<EPS*c.colsize()*std::abs(cdet)<<std::endl;
  }
  Assert(std::abs(c.Det()-cdet) < EPS*c.colsize()*std::abs(cdet),"Square CDet");

  tmv::Matrix<std::complex<T> > cinv = c.Inverse();
  tmv::Matrix<std::complex<T> > cid = c*cinv;
  Assert(Norm(cid-T(1)) < EPS*Norm(c)*Norm(cinv),"Square CInverse");

  tmv::Matrix<std::complex<T> > ctc = c.Adjoint() * c;
  tmv::Matrix<std::complex<T> > cata(4,4);
  c.InverseATA(cata);
  Assert(Norm(cata-ctc.Inverse()) < EPS*Norm(ctc)*Norm(ctc.Inverse()),
      "Square CInverseATA");

  tmv::Vector<std::complex<T> > e(4);
  e = b*std::complex<T>(1,2);
  e(1) += std::complex<T>(-1,5);
  e(2) -= std::complex<T>(-1,5);

  // test real / complex
  tmv::Vector<std::complex<T> > y = b/c;
  tmv::Vector<std::complex<T> > b3 = c*y;
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
  tmv::Vector<T> P(100,p);
  tmv::Vector<T> Q(100,q);
  tmv::Vector<T> R(100,r);
  tmv::Matrix<T,stor> M = P ^ Q;
  tmv::Matrix<std::complex<T>,stor> CM = P ^ (std::complex<T>(-4,10)*Q);
  M.diag().AddToAll(T(215));
  CM.diag().AddToAll(std::complex<T>(103,-53));
  M.row(23) *= T(12);
  M(12,1) -= T(142);
  CM(6,2) += std::complex<T>(23,89);
  CM.col(15) *= std::complex<T>(61,12);
  M.SubVector(65,5,1,3,29) *= T(2);
  M.SubVector(98,12,-1,2,18).AddToAll(T(197));
  CM.SubVector(53,0,1,3,31) *= std::complex<T>(2,-1);
  CM.SubVector(88,18,-1,2,23).AddToAll(std::complex<T>(197,174));
  M.DivideUsing(dt);
  CM.DivideUsing(dt);
  M.SaveDiv();
  CM.SaveDiv();
  tmv::Vector<T> S = R/M;
  tmv::Vector<T> R2 = M*S;
  Assert(M.CheckDecomp(),"CheckDecomp");
  if (showacc) {
    std::cout<<"R/M Norm(R2-R) = "<<Norm(R2-R)<<std::endl;
    std::cout<<"EPS*Norm(M)*Norm(R) = "<<EPS*Norm(M)*Norm(R)<<std::endl;
  }
  Assert(Norm(R2-R) < 10*EPS*Norm(M)*Norm(R),"Square R/M");
  S = R%M;
  R2 = S*M;
  if (showacc) {
    std::cout<<"R%M Norm(R2-R) = "<<Norm(R2-R)<<std::endl;
    std::cout<<"EPS*Norm(M)*Norm(R) = "<<EPS*Norm(M)*Norm(R)<<std::endl;
  }
  Assert(Norm(R2-R) < 10*EPS*Norm(M)*Norm(R),"Square R%M");
  tmv::Vector<std::complex<T> > CS = R/CM;
  tmv::Vector<std::complex<T> > CR2 = CM*CS;
  Assert(CM.CheckDecomp(),"CheckDecomp");
  if (showacc) {
    std::cout<<"R/CM Norm(CR2-R) = "<<Norm(CR2-R)<<std::endl;
    std::cout<<"EPS*Norm(CM)*Norm(R) = "<<EPS*Norm(CM)*Norm(R)<<std::endl;
  }
  Assert(Norm(CR2-R) < 10*EPS*Norm(CM)*Norm(R),"Square R/CM");
  CS = R%CM;
  CR2 = CS*CM;
  if (showacc) {
    std::cout<<"R%CM Norm(CR2-R) = "<<Norm(CR2-R)<<std::endl;
    std::cout<<"EPS*Norm(CM)*Norm(R) = "<<EPS*Norm(CM)*Norm(R)<<std::endl;
  }
  Assert(Norm(CR2-R) < 10*EPS*Norm(CM)*Norm(R),"Square R%CM");
  tmv::Vector<std::complex<T> > CR = std::complex<T>(3,-4)*R;
  CS = CR/CM;
  CR2 = CM*CS;
  if (showacc) {
    std::cout<<"CR/CM Norm(CR2-CR) = "<<Norm(CR2-CR)<<std::endl;
    std::cout<<"EPS*Norm(CM)*Norm(CR) = "<<EPS*Norm(CM)*Norm(CR)<<std::endl;
  }
  Assert(Norm(CR2-CR) < 10*EPS*Norm(CM)*Norm(CR),"Square CR/CM");
  CS = CR%CM;
  CR2 = CS*CM;
  if (showacc) {
    std::cout<<"CR%CM Norm(CR2-CR) = "<<Norm(CR2-CR)<<std::endl;
    std::cout<<"EPS*Norm(CM)*Norm(CR) = "<<EPS*Norm(CM)*Norm(CR)<<std::endl;
  }
  Assert(Norm(CR2-CR) < 10*EPS*Norm(CM)*Norm(CR),"Square CR%CM");
  CS = CR/M;
  CR2 = M*CS;
  if (showacc) {
    std::cout<<"CR/M Norm(CR2-CR) = "<<Norm(CR2-CR)<<std::endl;
    std::cout<<"EPS*Norm(M)*Norm(CR) = "<<EPS*Norm(M)*Norm(CR)<<std::endl;
  }
  Assert(Norm(CR2-CR) < 10*EPS*Norm(M)*Norm(CR),"Square CR/M");
  CS = CR%M;
  CR2 = CS*M;
  if (showacc) {
    std::cout<<"CR%M Norm(CR2-CR) = "<<Norm(CR2-CR)<<std::endl;
    std::cout<<"EPS*Norm(M)*Norm(CR) = "<<EPS*Norm(M)*Norm(CR)<<std::endl;
  }
  Assert(Norm(CR2-CR) < 10*EPS*Norm(M)*Norm(CR),"Square CR%M");

  tmv::Matrix<T,stor> a1 = m;
  tmv::Matrix<T,stor> a2 = m.Transpose();
  a2.row(1) *= T(3);
  a2.col(2) -= tmv::Vector<T>(4,4.);

  tmv::Matrix<std::complex<T>,stor> c1 = a1 * std::complex<T>(1,2);
  tmv::Matrix<std::complex<T>,stor> c2 = a2 * std::complex<T>(-3,4);
  c1.diag().AddToAll(std::complex<T>(3,1));
  c2.diag().AddToAll(std::complex<T>(-5,8));
  c1.row(3).AddToAll(std::complex<T>(1,-6));
  c2.row(0).AddToAll(std::complex<T>(-2,-11));

  TestMatrixDivArith<T>(dt,a1.View(),a2.View(),c1.View(),c2.View(),"Square"); 
#ifdef XTEST
  tmv::Matrix<T,stor,tmv::FortranStyle> a1f = a1;
  tmv::Matrix<T,stor,tmv::FortranStyle> a2f = a2;
  tmv::Matrix<std::complex<T>,stor,tmv::FortranStyle> c1f = c1;
  tmv::Matrix<std::complex<T>,stor,tmv::FortranStyle> c2f = c2;
  TestMatrixDivArith<T>(dt,a1f.View(),a2.View(),c1f.View(),c2.View(),
      "Square"); 
  TestMatrixDivArith<T>(dt,a1.View(),a2f.View(),c1.View(),c2f.View(),
      "Square"); 
  TestMatrixDivArith<T>(dt,a1f.View(),a2f.View(),c1f.View(),c2f.View(),
      "Square"); 
#endif

  tmv::Matrix<T,stor> a3(7,4);
  for(int i=0;i<7;++i) for(int j=0;j<4;++j) a3(i,j) = 1-3*i+2*j;
  tmv::Matrix<T,stor> a4 = a3.Transpose();
  a3.SubMatrix(2,6,0,4) += a1;
  a4.SubMatrix(0,4,1,5) -= a2;

  tmv::Matrix<std::complex<T>,stor> c3 = a3*std::complex<T>(1,2);
  tmv::Matrix<std::complex<T>,stor> c4 = c3.Adjoint();
  c3.SubMatrix(2,6,0,4) += c1;
  c4.SubMatrix(0,4,1,5) -= c2;
  c3.col(1) *= std::complex<T>(2,1);
  c3.row(2).AddToAll(std::complex<T>(-7,2));
  c4.col(3) *= std::complex<T>(-1,3);
  c4.row(0).AddToAll(std::complex<T>(1,9));

  TestMatrixDivArith<T>(dt,a1.View(),a3.View(),c1.View(),c3.View(),
      "Square/NonSquare");
  TestMatrixDivArith<T>(dt,a1.View(),a4.View(),c1.View(),c4.View(),
      "Square/NonSquare");

#ifdef XTEST
  tmv::Matrix<T,stor> a5(4,0,1);
  tmv::Matrix<T,stor> a6(0,4,1);
  tmv::Matrix<std::complex<T>,stor> c5 = a5;
  tmv::Matrix<std::complex<T>,stor> c6 = a6;

  TestMatrixDivArith<T>(dt,a1.View(),a5.View(),c1.View(),c5.View(),
      "Square/Degenerate");
  TestMatrixDivArith<T>(dt,a1.View(),a6.View(),c1.View(),c6.View(),
      "Square/Degenerate");
#endif

  if (stor == tmv::ColMajor) {
#ifdef XTEST
    TestSquareDiv<T,tmv::RowMajor>(dt);
  } else {
#endif
    std::cout<<"Square Matrix<"<<tmv::Type(T())<<"> Division using ";
    std::cout<<tmv::Text(dt)<<" passed all tests\n";
  }
}

template <class T> inline void TestSquareSV()
{
  tmv::Matrix<T> m(4,4);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) m(i,j) = T(2.+4*i-5*j);
  m(0,0) = T(14);
  m(1,0) = T(-2);
  m(2,0) = T(7);
  m(3,0) = T(-10);
  m.diag() *= T(30);

  m.DivideUsing(tmv::SV);
  m.SetDiv();
  tmv::Matrix<T> U = m.SVD().GetU();
  tmv::Vector<RealType(T)> S = m.SVD().GetS();
  tmv::Matrix<T> V = m.SVD().GetV();

  tmv::Matrix<T> m2 = m;
  m2.DivideUsing(tmv::SVS);
  m2.SetDiv();
  tmv::Matrix<T> m3 = m;
  m3.DivideUsing(tmv::SVU);
  m3.SetDiv();
  tmv::Matrix<T> m4 = m;
  m4.DivideUsing(tmv::SVV);
  m4.SetDiv();

  Assert(Norm(m2.SVD().GetS()-S) < EPS*(Norm(m)),"Square SVS S");
  Assert(Norm(m3.SVD().GetS()-S) < EPS*(Norm(m)),"Square SVU S");
  Assert(Norm(m4.SVD().GetS()-S) < EPS*(Norm(m)),"Square SVV S");

  Assert(Norm(m3.SVD().GetU()-U) < EPS*(Norm(m)),"Square SVU U");

  Assert(Norm(m4.SVD().GetV()-V) < EPS*(Norm(m)),"Square SVV V");

  tmv::Matrix<std::complex<T> > c(4,4);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) 
    c(i,j) = std::complex<T>(2.+4*i-5*j,3.-i);
  c(0,0) *= T(14);
  c(1,0) *= T(-2);
  c(2,0) *= T(7);
  c(3,0) *= T(-10);
  c.diag() *= T(30);

  c.DivideUsing(tmv::SV);
  c.SetDiv();
  tmv::Matrix<std::complex<T> > cU = c.SVD().GetU();
  tmv::Vector<T> cS = c.SVD().GetS();
  tmv::Matrix<std::complex<T> > cV = c.SVD().GetV();

  tmv::Matrix<std::complex<T> > c2 = c;
  c2.DivideUsing(tmv::SVS);
  c2.SetDiv();
  tmv::Matrix<std::complex<T> > c3 = c;
  c3.DivideUsing(tmv::SVU);
  c3.SetDiv();
  tmv::Matrix<std::complex<T> > c4 = c;
  c4.DivideUsing(tmv::SVV);
  c4.SetDiv();

  Assert(Norm(c2.SVD().GetS()-cS) < EPS*(Norm(c)),"Square C SVS S");
  Assert(Norm(c3.SVD().GetS()-cS) < EPS*(Norm(c)),"Square C SVU S");
  Assert(Norm(c4.SVD().GetS()-cS) < EPS*(Norm(c)),"Square C SVV S");

  Assert(Norm(c3.SVD().GetU()-cU) < EPS*(Norm(c)),"Square C SVU U");

  Assert(Norm(c4.SVD().GetV()-cV) < EPS*(Norm(c)),"Square C SVV V");
}

template <class T, tmv::StorageType stor> inline void TestNonSquareDiv(tmv::DivType dt)
{
  tmv::Matrix<T,stor> m(6,4);
  m.SaveDiv();
  m.DivideUsing(dt);
  for(int i=0;i<6;++i) for(int j=0;j<4;++j) m(i,j) = 2.+4*i-5*j;
  m(0,0) = 14.;
  m(1,0) = -2.;
  m(2,0) = 7.;
  m(3,0) = -10.;
  m(2,2) = 30.;

  tmv::Vector<T> x(4);
  x(0) = 2;
  x(1) = -10;
  x(2) = 5;
  x(3) = -5;

  m.SetDiv();
  if (showdiv) {
    Assert(m.CheckDecomp(&std::cout),"CheckDecomp");
  } else {
    Assert(m.CheckDecomp(),"CheckDecomp");
  }

  tmv::Vector<T> b = m * x;
  tmv::Vector<T> x2 = b/m;
  Assert(Norm(x2-x) < EPS*Norm(m)*Norm(b),"NonSquare exact b/m");

  tmv::Vector<T> b2 = x%m;
  x2 = b2*m;
  Assert(Norm(x2-x) < EPS*Norm(m)*Norm(x),"NonSquare x%m");

  b(0) += 100.;
  x = b/m;
  b2 = m*x;
  T refnorm = Norm(b2-b);
  tmv::Vector<T> dx = T(sqrt(EPS))*Norm(x)*tmv::BasisVector<T>(4,0);
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

  tmv::Matrix<T> minv = m.Inverse();
  tmv::Matrix<T> id = minv*m;
  tmv::Matrix<T> nonid = m*minv;
  if (showacc) {
    std::cout<<"minv = "<<minv<<std::endl;
    std::cout<<"minv*m = "<<id<<std::endl;
    std::cout<<"m*minv = "<<nonid<<std::endl;
    std::cout<<"Norm(id-I) = "<<Norm(id-T(1))<<std::endl;
  }
  Assert(Norm(id-T(1)) < EPS*Norm(m)*Norm(minv),"NonSquare Inverse");
  Assert(Norm(nonid-nonid.Transpose()) < EPS*Norm(m)*Norm(minv),
      "NonSquare Pseudo-Inverse");

  tmv::Matrix<T> mata(4,4);
  m.InverseATA(mata);
  tmv::Matrix<T> mtm = m.Transpose()*m;
  Assert(Norm(mata-mtm.Inverse()) < EPS*Norm(mtm),"NonSquare InverseATA");

  tmv::Matrix<std::complex<T>,stor> c(6,4);
  c.SaveDiv();
  c.DivideUsing(dt);
  c = m;
  c(2,3) += std::complex<T>(2,3);
  c(1,0) *= std::complex<T>(0,2);
  c.col(1) *= std::complex<T>(-1,3);
  c.row(3) += tmv::Vector<std::complex<T> >(4,std::complex<T>(1,9));

  tmv::Vector<std::complex<T> > y(4);
  y(0) = std::complex<T>(2,9);
  y(1) = std::complex<T>(-10,4);
  y(2) = std::complex<T>(5,-1);
  y(3) = std::complex<T>(-5,-2);

  tmv::Vector<std::complex<T> > e = c * y;
  tmv::Vector<std::complex<T> > y2 = e/c;

  if (showdiv) {
    Assert(c.CheckDecomp(&std::cout),"CheckDecomp");
  } else {
    Assert(c.CheckDecomp(),"CheckDecomp");
  }

  Assert(Norm(y2-y) < EPS*Norm(c)*Norm(y),"NonSquare exact e/c");

  tmv::Vector<std::complex<T> > e2 = y%c;
  y2 = e2*c;
  Assert(Norm(y2-y) < EPS*Norm(c)*Norm(y),"NonSquare e%c");

  tmv::Matrix<std::complex<T> > cinv = c.Inverse();
  tmv::Matrix<std::complex<T> > cid = cinv*c;
  Assert(Norm(cid-T(1)) < EPS*Norm(c),"NonSquare CInverse");
  tmv::Matrix<std::complex<T> > cnonid = c*cinv;
  Assert(Norm(cnonid-cnonid.Adjoint()) < EPS*Norm(c),
      "NonSquare CPseudo-Inverse");

  tmv::Matrix<std::complex<T> > cata(4,4);
  c.InverseATA(cata);
  tmv::Matrix<std::complex<T> > ctc = c.Adjoint()*c;
  Assert(Norm(cata-ctc.Inverse()) < EPS*Norm(ctc),"NonSquare CInverseATA");

  // Test short matrix (M < N)
  tmv::Matrix<T,stor> ms = m.Transpose();
  ms.DivideUsing(dt);
  ms.SaveDiv();

  b = x * ms;
  x2 = b%ms;
  if (showdiv) {
    Assert(ms.CheckDecomp(&std::cout),"CheckDecomp");
  } else {
    Assert(ms.CheckDecomp(),"CheckDecomp");
  }
  Assert(Norm(x2-x) < EPS*Norm(ms)*Norm(x),"NonSquare exact b%ms");

  b2 = x/ms;
  x2 = ms*b2;
  Assert(Norm(x2-x) < EPS*Norm(ms)*Norm(x),"NonSquare x/ms");

  // Test really long matrix
  tmv::Matrix<std::complex<T>,stor> a(30,10);
  a.DivideUsing(dt);
  for(int i=0;i<30;++i) for(int j=0;j<10;++j) a(i,j) = 7.-13.*i+11.*j;
  a.SubMatrix(0,10,0,10) += std::complex<T>(30.,20.);
  a.SubMatrix(10,20,0,10) -= std::complex<T>(50.,-123.);
  a.SubMatrix(20,30,0,10) += std::complex<T>(10.,-75.);
  a.SubMatrix(1,10,1,10) += std::complex<T>(99.,100.);
  a.SubMatrix(2,10,2,10) -= std::complex<T>(51.,37.);

  tmv::Vector<std::complex<T> > s(10);
  for(int i=0;i<10;++i) s(i) = i+2;

  tmv::Vector<std::complex<T> > t = a * s;
  tmv::Vector<std::complex<T> > s2 = t/a;
  Assert(Norm(s2-s) < EPS*Norm(a)*Norm(s),"NonSquare t/a");

  tmv::Vector<std::complex<T> > t2 = s%a;
  s2 = t2*a;
  Assert(Norm(s2-s) < EPS*Norm(a)*Norm(s),"NonSquare t%a");

  // Test QR Update/Downdate:

  tmv::Matrix<std::complex<T>,stor> q30 = a;
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,stor> r30(10,10);
  QR_Decompose(q30.View(),r30.View());
  Assert(Norm(q30*r30-a) < EPS*a.NormSq(),"QR_Decompose");
  Assert(Norm(r30.Adjoint()*r30-a.Adjoint()*a) < EPS*r30.NormSq(),
      "QR_Decompose (RtR)");

  tmv::Matrix<std::complex<T>,stor> q10 = a.Rows(0,10);
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,stor> r10(10,10);
  QR_Decompose(q10.View(),r10.View());
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,stor> r = r10;
  tmv::Matrix<std::complex<T>,stor> a1030 = a.Rows(10,30);
  QR_Update(r.View(),a1030.View());
  Assert(Norm(r.Adjoint()*r-r30.Adjoint()*r30) < EPS*a.NormSq(),
      "QR_Update");
  r = r30;

  a1030 = a.Rows(10,30);
  QR_Downdate(r.View(),a1030.View());
  Assert(Norm(r.Adjoint()*r-r10.Adjoint()*r10) < EPS*a.NormSq(),
      "QR_Downdate");
  r = r10;

  tmv::Matrix<std::complex<T>,stor> a1020 = a.Rows(10,20);
  tmv::Matrix<std::complex<T>,stor> a2030 = a.Rows(20,30);
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

  tmv::Matrix<std::complex<T>,stor> q29 = a.Rows(0,29);
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,stor> r29(10,10);
  QR_Decompose(q29.View(),r29.View());
  r = r30;
  tmv::Vector<std::complex<T> > a29 = a.row(29);
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
  tmv::Matrix<std::complex<T>,stor> q = a;
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

  tmv::Matrix<T,stor> a1 = m;
  tmv::Matrix<T,stor> a2 = m.Transpose() * m;
  tmv::Matrix<T,stor> a3 = m * m.Transpose();
  a2.row(1) *= T(3);
  a2.col(2).AddToAll(-4);
  a3.row(5) *= T(7);
  a3.col(3).AddToAll(7);
  tmv::Matrix<std::complex<T>,stor> c1 = a1 * std::complex<T>(1,2);
  tmv::Matrix<std::complex<T>,stor> c2 = a2 * std::complex<T>(-3,4);
  tmv::Matrix<std::complex<T>,stor> c3 = a3 * std::complex<T>(-4,8);

  TestMatrixDivArith<T>(dt,a1.View(),a2.View(),c1.View(),c2.View(),
      "NonSquare/Square"); 
  TestMatrixDivArith<T>(dt,a1.View(),a3.View(),c1.View(),c3.View(),
      "NonSquare/Square"); 

  tmv::Matrix<T,stor> a4(7,4);
  for(int i=0;i<7;++i) for(int j=0;j<4;++j) a4(i,j) = 1-3*i+2*j;
  tmv::Matrix<T,stor> a5 = a4.Transpose();
  a4.SubMatrix(0,6,0,4) += a1;
  a5.SubMatrix(0,4,1,5) -= a2;
  tmv::Matrix<std::complex<T>,stor> c4 = a4*std::complex<T>(1,2);
  tmv::Matrix<std::complex<T>,stor> c5 = c4.Adjoint();
  c4.SubMatrix(0,6,0,4) += c1;
  c5.SubMatrix(0,4,1,5) -= c2;
  c4.col(1) *= std::complex<T>(2,1);
  c4.row(2).AddToAll(std::complex<T>(-7,2));
  c5.col(3) *= std::complex<T>(-1,3);
  c5.row(0).AddToAll(std::complex<T>(1,9));

  tmv::Matrix<T,stor> a6(9,6);
  for(int i=0;i<9;++i) for(int j=0;j<6;++j) a6(i,j) = 5+2*i-2*j;
  tmv::Matrix<T,stor> a7 = a6.Transpose();
  a6.SubMatrix(2,8,1,5) += a1;
  a7.SubMatrix(0,6,4,8) -= T(2)*a1;
  tmv::Matrix<std::complex<T>,stor> c6 = a6*std::complex<T>(1,2);
  tmv::Matrix<std::complex<T>,stor> c7 = c6.Adjoint();
  c6.SubMatrix(2,8,1,5) += c1;
  c7.SubMatrix(0,6,4,8) -= T(2)*c1;
  c6.col(1) *= std::complex<T>(2,1);
  c6.row(5).AddToAll(std::complex<T>(-7,2));
  c7.col(7) *= std::complex<T>(-1,3);
  c7.row(4).AddToAll(std::complex<T>(1,9));

  TestMatrixDivArith<T>(dt,a1.View(),a4.View(),c1.View(),c4.View(),
      "NonSquare/NonSquare");
  TestMatrixDivArith<T>(dt,a1.View(),a5.View(),c1.View(),c5.View(),
      "NonSquare/NonSquare");
  TestMatrixDivArith<T>(dt,a1.View(),a6.View(),c1.View(),c6.View(),
      "NonSquare/NonSquare");
  TestMatrixDivArith<T>(dt,a1.View(),a7.View(),c1.View(),c7.View(),
      "NonSquare/NonSquare");

#ifdef XTEST
  tmv::Matrix<T,stor> a8(4,0,1);
  tmv::Matrix<T,stor> a9(0,4,1);
  tmv::Matrix<T,stor> a10(6,0,1);
  tmv::Matrix<T,stor> a11(0,6,1);
  tmv::Matrix<std::complex<T>,stor> c8 = a8;
  tmv::Matrix<std::complex<T>,stor> c9 = a9;
  tmv::Matrix<std::complex<T>,stor> c10 = a10;
  tmv::Matrix<std::complex<T>,stor> c11 = a11;

  TestMatrixDivArith<T>(dt,a1.View(),a8.View(),c1.View(),c8.View(),
      "NonSquare/Degenerate");
  TestMatrixDivArith<T>(dt,a1.View(),a9.View(),c1.View(),c9.View(),
      "NonSquare/Degenerate");
  TestMatrixDivArith<T>(dt,a1.View(),a10.View(),c1.View(),c10.View(),
      "NonSquare/Degenerate");
  TestMatrixDivArith<T>(dt,a1.View(),a11.View(),c1.View(),c11.View(),
      "NonSquare/Degenerate");
#endif

  if (stor == tmv::ColMajor) {
#ifdef XTEST
    TestNonSquareDiv<T,tmv::RowMajor>(dt);
  } else {
#endif
    std::cout<<"NonSquare Matrix<"<<tmv::Type(T())<<"> Division using ";
    std::cout<<tmv::Text(dt)<<" passed all tests\n";
  }
}

template <class T> inline void TestNonSquareSV()
{
  tmv::Matrix<T> m(6,4);
  for(int i=0;i<6;++i) for(int j=0;j<4;++j) m(i,j) = T(2.+4*i-5*j);
  m(0,0) = T(14);
  m(1,0) = T(-2);
  m(2,0) = T(7);
  m(3,0) = T(-10);
  m.diag() *= T(30);

  m.DivideUsing(tmv::SV);
  m.SetDiv();
  tmv::Matrix<T> U = m.SVD().GetU();
  tmv::Vector<RealType(T)> S = m.SVD().GetS();
  tmv::Matrix<T> V = m.SVD().GetV();

  tmv::Matrix<T> m2 = m;
  m2.DivideUsing(tmv::SVS);
  m2.SetDiv();
  tmv::Matrix<T> m3 = m;
  m3.DivideUsing(tmv::SVU);
  m3.SetDiv();
  tmv::Matrix<T> m4 = m;
  m4.DivideUsing(tmv::SVV);
  m4.SetDiv();

  Assert(Norm(m2.SVD().GetS()-S) < EPS*(Norm(m)),"NonSquare SVS S");
  Assert(Norm(m3.SVD().GetS()-S) < EPS*(Norm(m)),"NonSquare SVU S");
  Assert(Norm(m4.SVD().GetS()-S) < EPS*(Norm(m)),"NonSquare SVV S");

  Assert(Norm(m3.SVD().GetU()-U) < EPS*(Norm(m)),"NonSquare SVU U");

  Assert(Norm(m4.SVD().GetV()-V) < EPS*(Norm(m)),"NonSquare SVV V");

  tmv::Matrix<T> a(15,4);
  for(int i=0;i<15;++i) for(int j=0;j<4;++j) a(i,j) = T(2.+4*i-5*j);
  a(0,0) = T(14);
  a(1,0) = T(-2);
  a(2,0) = T(7);
  a(3,0) = T(-10);
  a.diag() *= T(30);

  a.DivideUsing(tmv::SV);
  a.SetDiv();
  tmv::Matrix<T> aU = a.SVD().GetU();
  tmv::Vector<RealType(T)> aS = a.SVD().GetS();
  tmv::Matrix<T> aV = a.SVD().GetV();

  tmv::Matrix<T> a2 = a;
  a2.DivideUsing(tmv::SVS);
  a2.SetDiv();
  tmv::Matrix<T> a3 = a;
  a3.DivideUsing(tmv::SVU);
  a3.SetDiv();
  tmv::Matrix<T> a4 = a;
  a4.DivideUsing(tmv::SVV);
  a4.SetDiv();

  Assert(Norm(a2.SVD().GetS()-aS) < EPS*(Norm(a)),"Tall NonSquare SVS S");
  Assert(Norm(a3.SVD().GetS()-aS) < EPS*(Norm(a)),"Tall NonSquare SVU S");
  Assert(Norm(a4.SVD().GetS()-aS) < EPS*(Norm(a)),"Tall NonSquare SVV S");

  Assert(Norm(a3.SVD().GetU()-aU) < EPS*(Norm(a)),"Tall NonSquare SVU U");

  Assert(Norm(a4.SVD().GetV()-aV) < EPS*(Norm(a)),"Tall NonSquare SVV V");

  tmv::Matrix<std::complex<T> > c(6,4);
  for(int i=0;i<6;++i) for(int j=0;j<4;++j)
    c(i,j) = std::complex<T>(2.+4*i-5*j,3.-i);
  c(0,0) = T(14);
  c(1,0) = T(-2);
  c(2,0) = T(7);
  c(3,0) = T(-10);
  c.diag() *= T(30);

  c.DivideUsing(tmv::SV);
  c.SetDiv();
  tmv::Matrix<std::complex<T> > cU = c.SVD().GetU();
  tmv::Vector<T> cS = c.SVD().GetS();
  tmv::Matrix<std::complex<T> > cV = c.SVD().GetV();

  tmv::Matrix<std::complex<T> > c2 = c;
  c2.DivideUsing(tmv::SVS);
  c2.SetDiv();
  tmv::Matrix<std::complex<T> > c3 = c;
  c3.DivideUsing(tmv::SVU);
  c3.SetDiv();
  tmv::Matrix<std::complex<T> > c4 = c;
  c4.DivideUsing(tmv::SVV);
  c4.SetDiv();

  Assert(Norm(c2.SVD().GetS()-cS) < EPS*(Norm(c)),"NonSquare C SVS S");
  Assert(Norm(c3.SVD().GetS()-cS) < EPS*(Norm(c)),"NonSquare C SVU S");
  Assert(Norm(c4.SVD().GetS()-cS) < EPS*(Norm(c)),"NonSquare C SVV S");

  Assert(Norm(c3.SVD().GetU()-cU) < EPS*(Norm(c)),"NonSquare C SVU U");

  Assert(Norm(c4.SVD().GetV()-cV) < EPS*(Norm(c)),"NonSquare C SVV V");

  tmv::Matrix<std::complex<T> > d(15,4);
  for(int i=0;i<15;++i) for(int j=0;j<4;++j)
    d(i,j) = std::complex<T>(2.+4*i-5*j,3.-i);
  d(0,0) = T(14);
  d(1,0) = T(-2);
  d(2,0) = T(7);
  d(3,0) = T(-10);
  d.diag() *= T(30);

  d.DivideUsing(tmv::SV);
  d.SetDiv();
  tmv::Matrix<std::complex<T> > dU = d.SVD().GetU();
  tmv::Vector<T> dS = d.SVD().GetS();
  tmv::Matrix<std::complex<T> > dV = d.SVD().GetV();

  tmv::Matrix<std::complex<T> > d2 = d;
  d2.DivideUsing(tmv::SVS);
  d2.SetDiv();
  tmv::Matrix<std::complex<T> > d3 = d;
  d3.DivideUsing(tmv::SVU);
  d3.SetDiv();
  tmv::Matrix<std::complex<T> > d4 = d;
  d4.DivideUsing(tmv::SVV);
  d4.SetDiv();

  Assert(Norm(d2.SVD().GetS()-dS) < EPS*(Norm(d)),"Tall NonSquare C SVS S");
  Assert(Norm(d3.SVD().GetS()-dS) < EPS*(Norm(d)),"Tall NonSquare C SVU S");
  Assert(Norm(d4.SVD().GetS()-dS) < EPS*(Norm(d)),"Tall NonSquare C SVV S");

  Assert(Norm(d3.SVD().GetU()-dU) < EPS*(Norm(d)),"Tall NonSquare C SVU U");

  Assert(Norm(d4.SVD().GetV()-dV) < EPS*(Norm(d)),"Tall NonSquare C SVV V");
}

template <class T, tmv::StorageType stor> inline void TestSingularDiv(tmv::DivType dt)
{
  tmv::Matrix<T,stor> m(4,4);
  m.SaveDiv();
  m.DivideUsing(dt);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) m(i,j) = 2.+4*i-5*j;
  m(2,2) = 30.;

  tmv::Vector<T> x(4);
  x(0) = 2;
  x(1) = -10;
  x(2) = 5;
  x(3) = -5;

  m.SetDiv();
  if (showdiv) {
    Assert(m.CheckDecomp(&std::cout),"CheckDecomp");
  } else {
    Assert(m.CheckDecomp(),"CheckDecomp");
  }

  tmv::Vector<T> b = m * x;
  tmv::Vector<T> x2 = b/m;
  tmv::Vector<T> b2 = m*x2;
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
  tmv::Vector<T> dx = T(3*sqrt(EPS))*tmv::BasisVector<T>(4,0);
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

  tmv::Matrix<T> minv = m.Inverse();

  if (showacc) {
    std::cout<<"m = "<<m<<std::endl;
    std::cout<<"minv = "<<minv<<std::endl;
    std::cout<<"minv*m = "<<minv*m<<std::endl;
    std::cout<<"m*minv = "<<m*minv<<std::endl;
    std::cout<<"m*minv*m = "<<m*minv*m<<std::endl;
    std::cout<<"minv*m*minv = "<<minv*m*minv<<std::endl;
    std::cout<<"m*minv-(m*minv)T = "<<m*minv-Transpose(m*minv)<<std::endl;
    std::cout<<"minv*m-(minv*m)T = "<<minv*m-Transpose(minv*m)<<std::endl;
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
  tmv::Matrix<std::complex<T>,stor> mm(30,30);
  mm.SaveDiv();
  mm.DivideUsing(dt);
  for(int i=0;i<30;++i) for(int j=0;j<30;++j) mm(i,j) = 4.-17.*i+23.*j;
  mm(20,20) += std::complex<T>(200.,-999.);
  mm(12,12) += std::complex<T>(500.,-104.);
  mm(7,7) += std::complex<T>(300.,123.);
  mm(28,28) += std::complex<T>(700.,231.);
  mm(24,24) += std::complex<T>(400.,-120.);

  tmv::Vector<std::complex<T> > xx(30);
  for(int i=0;i<30;++i) xx(i) = 10.+i;
  tmv::Vector<std::complex<T> > bb = mm*xx;
  tmv::Vector<std::complex<T> > xx2 = bb/mm;
  tmv::Vector<std::complex<T> > bb2 = mm*xx2;
  Assert(mm.CheckDecomp(),"CheckDecomp");
  if (showacc) {
    std::cout<<"Norm(bb2-bb) = "<<Norm(bb2-bb);
    std::cout<<", EPS*Norm(mm)*Norm(xx) = "<<EPS*Norm(mm)*Norm(xx)<<std::endl;
  }
  Assert(Norm(bb2-bb) < 30*EPS*Norm(mm)*Norm(xx),"Singular exact bb/mm");

  bb = xx * mm;
  xx2 = bb%mm;
  bb2 = xx2*mm;
  if (showacc) {
    std::cout<<"Norm(bb2-bb) = "<<Norm(bb2-bb);
    std::cout<<", EPS*Norm(mm)*Norm(xx) = "<<EPS*Norm(mm)*Norm(xx)<<std::endl;
  }
  Assert(Norm(bb2-bb) < 30*EPS*Norm(mm)*Norm(xx),"Singular exact bb%mm");

  if (stor == tmv::ColMajor) {
#ifdef XTEST
    TestSingularDiv<T,tmv::RowMajor>(dt);
  } else {
#endif
    std::cout<<"Singular Matrix<"<<tmv::Type(T())<<"> Division using ";
    std::cout<<tmv::Text(dt)<<" passed all tests\n";
  }
}

template <class T> inline void TestSingularSV()
{
  tmv::Matrix<T> m(4,4);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) m(i,j) = 2.+4*i-5*j;
  m(2,2) = 30.;

  m.DivideUsing(tmv::SV);
  m.SetDiv();
  tmv::Matrix<T> U = m.SVD().GetU();
  tmv::Vector<RealType(T)> S = m.SVD().GetS();
  tmv::Matrix<T> V = m.SVD().GetV();

  tmv::Matrix<T> m2 = m;
  m2.DivideUsing(tmv::SVS);
  m2.SetDiv();
  tmv::Matrix<T> m3 = m;
  m3.DivideUsing(tmv::SVU);
  m3.SetDiv();
  tmv::Matrix<T> m4 = m;
  m4.DivideUsing(tmv::SVV);
  m4.SetDiv();

  Assert(Norm(m2.SVD().GetS()-S) < EPS*(Norm(m)),"Singular SVS S");
  Assert(Norm(m3.SVD().GetS()-S) < EPS*(Norm(m)),"Singular SVU S");
  Assert(Norm(m4.SVD().GetS()-S) < EPS*(Norm(m)),"Singular SVV S");

  tmv::Matrix<std::complex<T> > c(4,4);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) 
    c(i,j) = std::complex<T>(2.+4*i-5*j,3.-i);
  m(2,2) += 30.;

  c.DivideUsing(tmv::SV);
  c.SetDiv();
  tmv::Matrix<std::complex<T> > cU = c.SVD().GetU();
  tmv::Vector<T> cS = c.SVD().GetS();
  tmv::Matrix<std::complex<T> > cV = c.SVD().GetV();

  tmv::Matrix<std::complex<T> > c2 = c;
  c2.DivideUsing(tmv::SVS);
  c2.SetDiv();
  tmv::Matrix<std::complex<T> > c3 = c;
  c3.DivideUsing(tmv::SVU);
  c3.SetDiv();
  tmv::Matrix<std::complex<T> > c4 = c;
  c4.DivideUsing(tmv::SVV);
  c4.SetDiv();

  Assert(Norm(c2.SVD().GetS()-cS) < EPS*(Norm(c)),"Singular C SVS S");
  Assert(Norm(c3.SVD().GetS()-cS) < EPS*(Norm(c)),"Singular C SVU S");
  Assert(Norm(c4.SVD().GetS()-cS) < EPS*(Norm(c)),"Singular C SVV S");
}

template <class T> void TestAllMatrixDiv()
{
  TestSquareDiv<T,tmv::ColMajor>(tmv::LU);
  TestSquareDiv<T,tmv::ColMajor>(tmv::QR);
  TestSquareDiv<T,tmv::ColMajor>(tmv::QRP);
  TestSquareDiv<T,tmv::ColMajor>(tmv::SV);
  TestSquareSV<T>();
  TestNonSquareDiv<T,tmv::ColMajor>(tmv::QR);
  TestNonSquareDiv<T,tmv::ColMajor>(tmv::QRP);
  TestNonSquareDiv<T,tmv::ColMajor>(tmv::SV);
  TestNonSquareSV<T>();
  TestSingularDiv<T,tmv::ColMajor>(tmv::QRP);
  TestSingularDiv<T,tmv::ColMajor>(tmv::SV);
  TestSingularSV<T>();
}

template void TestAllMatrixDiv<double>();
#ifndef NOFLOAT
template void TestAllMatrixDiv<float>();
#endif
#ifdef LONGDOUBLE
template void TestAllMatrixDiv<long double>();
#endif
