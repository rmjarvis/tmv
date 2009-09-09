
#define STARTAT 0

#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_TestMatrixDiv.h"

using tmv::Matrix;
using tmv::BandMatrix;
using tmv::BandMatrixView;
using tmv::RowMajor;
using tmv::ColMajor;
using tmv::DiagMajor;
using tmv::UnitDiag;
using tmv::NonUnitDiag;

template <class T> void TestBandDiv(tmv::DivType dt)
{
  const int N = 10;

  vector<BandMatrixView<T> > b;
  vector<BandMatrixView<complex<T> > > cb;

  Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-2*j;
  a1.diag().AddToAll(T(10)*N);
  Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+j;
  a2.diag().AddToAll(T(10)*N);

  Matrix<complex<T> > ca1 = a1 * complex<T>(3,-4);
  Matrix<complex<T> > ca2 = a2 * complex<T>(3,-4);

  Vector<T> v1(N);
  Vector<T> v2(N-1);
  for (int i=0; i<N; ++i) v1(i) = 16.-3*i; 
  for (int i=0; i<N-1; ++i) v2(i) = -7.+2*i; 
  Vector<complex<T> > cv1 = v1 * complex<T>(-1,2);
  Vector<complex<T> > cv2 = v2 * complex<T>(-1,2);

  BandMatrix<T,RowMajor> B1(a1,3,3);
  BandMatrix<complex<T>,RowMajor> CB1(ca1,3,3);
  b.push_back(B1.SubBandMatrix(0,N,0,N,1,3));
  cb.push_back(CB1.SubBandMatrix(0,N,0,N,1,3));
  b.push_back(B1.SubBandMatrix(0,N,0,N,3,1));
  cb.push_back(CB1.SubBandMatrix(0,N,0,N,3,1));
  b.push_back(B1.SubBandMatrix(0,N,0,N,0,3));
  cb.push_back(CB1.SubBandMatrix(0,N,0,N,0,3));
  b.push_back(B1.SubBandMatrix(0,N,0,N,3,0));
  cb.push_back(CB1.SubBandMatrix(0,N,0,N,3,0));
  BandMatrix<T,ColMajor> B2(a1,3,3);
  BandMatrix<complex<T>,ColMajor> CB2(ca1,3,3);
  b.push_back(B2.SubBandMatrix(0,N,0,N,1,3));
  cb.push_back(CB2.SubBandMatrix(0,N,0,N,1,3));
  BandMatrix<T,DiagMajor> B3(a1,3,3);
  BandMatrix<complex<T>,DiagMajor> CB3(ca1,3,3);
  b.push_back(B3.SubBandMatrix(0,N,0,N,1,3));
  cb.push_back(CB3.SubBandMatrix(0,N,0,N,1,3));
  b.push_back(B3.SubBandMatrix(0,N,0,N,1,1));
  cb.push_back(CB3.SubBandMatrix(0,N,0,N,1,1));

  BandMatrix<T,ColMajor> B8(a2,6,6);
  BandMatrix<complex<T>,ColMajor> CB8(ca2,6,6);

#ifdef XTEST
  BandMatrix<T,DiagMajor> B4 = tmv::UpperBiDiagMatrix(v1,v2);
  BandMatrix<complex<T>,DiagMajor> CB4 = tmv::UpperBiDiagMatrix(cv1,cv2);
  BandMatrix<T,DiagMajor> B5 = tmv::LowerBiDiagMatrix(v2,v1);
  BandMatrix<complex<T>,DiagMajor> CB5 = tmv::LowerBiDiagMatrix(cv2,cv1);
  BandMatrix<T,DiagMajor> B6 = tmv::TriDiagMatrix(v2,v1,v2);
  BandMatrix<complex<T>,DiagMajor> CB6 = tmv::TriDiagMatrix(cv2,cv1,cv2);
  BandMatrix<T,RowMajor> B7(a1,3,N-1);
  BandMatrix<complex<T>,RowMajor> CB7(ca1,3,N-1);
  BandMatrix<T,ColMajor> B7b(a1,3,N-1);
  BandMatrix<complex<T>,ColMajor> CB7b(ca1,3,N-1);
  BandMatrix<T,RowMajor> B7c(a1,3,N-2);
  BandMatrix<complex<T>,RowMajor> CB7c(ca1,3,N-2);
  BandMatrix<T,ColMajor> B7d(a1,3,N-2);
  BandMatrix<complex<T>,ColMajor> CB7d(ca1,3,N-2);
  if (doallarith) {
    b.push_back(B1.View());
    cb.push_back(CB1.View());
    b.push_back(B2.View());
    cb.push_back(CB2.View());
    b.push_back(B2.SubBandMatrix(0,N,0,N,0,3));
    cb.push_back(CB2.SubBandMatrix(0,N,0,N,0,3));
    b.push_back(B2.SubBandMatrix(0,N,0,N,3,0));
    cb.push_back(CB2.SubBandMatrix(0,N,0,N,3,0));
    b.push_back(B3.View());
    cb.push_back(CB3.View());
    b.push_back(B3.SubBandMatrix(0,N,0,N,0,3));
    cb.push_back(CB3.SubBandMatrix(0,N,0,N,0,3));
    b.push_back(B3.SubBandMatrix(0,N,0,N,3,0));
    cb.push_back(CB3.SubBandMatrix(0,N,0,N,3,0));
    b.push_back(B4.View());
    cb.push_back(CB4.View());
    b.push_back(B5.View());
    cb.push_back(CB5.View());
    b.push_back(B6.View());
    cb.push_back(CB6.View());
    b.push_back(B7.View());
    cb.push_back(CB7.View());
    b.push_back(B7b.View());
    cb.push_back(CB7b.View());
    b.push_back(B7c.View());
    cb.push_back(CB7c.View());
    b.push_back(B7d.View());
    cb.push_back(CB7d.View());
  }

  b.push_back(B8.SubBandMatrix(0,2*N,0,2*N,3,3,2,2));
  cb.push_back(CB8.SubBandMatrix(0,2*N,0,2*N,3,3,2,2));

  BandMatrix<T,RowMajor> B8c(a2,6,6);
  BandMatrix<complex<T>,RowMajor> CB8c(ca2,6,6);
  BandMatrix<T,DiagMajor> B8b(a2,6,6);
  BandMatrix<complex<T>,DiagMajor> CB8b(ca2,6,6);
  if (doallarith) {
    b.push_back(B8.SubBandMatrix(0,2*N,0,2*N,0,3,2,2));
    cb.push_back(CB8.SubBandMatrix(0,2*N,0,2*N,0,3,2,2));
    b.push_back(B8.SubBandMatrix(0,2*N,0,2*N,3,0,2,2));
    cb.push_back(CB8.SubBandMatrix(0,2*N,0,2*N,3,0,2,2));
    b.push_back(B8c.SubBandMatrix(0,2*N,0,2*N,0,3,2,2));
    cb.push_back(CB8c.SubBandMatrix(0,2*N,0,2*N,0,3,2,2));
    b.push_back(B8c.SubBandMatrix(0,2*N,0,2*N,3,0,2,2));
    cb.push_back(CB8c.SubBandMatrix(0,2*N,0,2*N,3,0,2,2));
    b.push_back(B8b.SubBandMatrix(0,2*N,0,2*N,0,3,2,2));
    cb.push_back(CB8b.SubBandMatrix(0,2*N,0,2*N,0,3,2,2));
    b.push_back(B8b.SubBandMatrix(0,2*N,0,2*N,3,0,2,2));
    cb.push_back(CB8b.SubBandMatrix(0,2*N,0,2*N,3,0,2,2));
  }
#endif

  size_t nsquare = b.size();

  b.push_back(B8.SubBandMatrix(0,N,0,2*N,3,3));
  cb.push_back(CB8.SubBandMatrix(0,N,0,2*N,3,3));
  b.push_back(B8.SubBandMatrix(0,N,0,N+2,4,4));
  cb.push_back(CB8.SubBandMatrix(0,N,0,N+2,4,4));

#ifdef XTEST
  BandMatrix<T,DiagMajor> B9 = tmv::UpperBiDiagMatrix(v1,v1);
  BandMatrix<complex<T>,DiagMajor> CB9 = tmv::UpperBiDiagMatrix(cv1,cv1);
  BandMatrix<T,DiagMajor> B10 = tmv::LowerBiDiagMatrix(v2,v2);
  BandMatrix<complex<T>,DiagMajor> CB10 = tmv::LowerBiDiagMatrix(cv2,cv2);
  BandMatrix<T,DiagMajor> B11 = tmv::TriDiagMatrix(v2,v1,v1);
  BandMatrix<complex<T>,DiagMajor> CB11 = tmv::TriDiagMatrix(cv2,cv1,cv1);
  if (doallarith) {
    b.push_back(B9.View());
    cb.push_back(CB9.View());
    b.push_back(B10.View());
    cb.push_back(CB10.View());
    b.push_back(B11.View());
    cb.push_back(CB11.View());
  }
#endif

  size_t ntot = b.size();

  Matrix<T> a3 = a1.Cols(0,N/2);
  Matrix<complex<T> > ca3 = ca1.Cols(0,N/2);
  Matrix<T> a4 = a1.Rows(0,N/2);
  Matrix<complex<T> > ca4 = ca1.Rows(0,N/2);
  Matrix<T> a5 = a1.Cols(0,0);
  Matrix<complex<T> > ca5 = ca1.Cols(0,0);
  Matrix<T> a6 = a1.Rows(0,0);
  Matrix<complex<T> > ca6 = ca1.Rows(0,0);
  BandMatrix<T> b1 = B1;
  BandMatrix<complex<T> > cb1 = b1 * complex<T>(2.,-3.);
  BandMatrix<T> b2 = b1.SubBandMatrix(0,N/2,0,N,b1.nlo(),b1.nhi());
  BandMatrix<complex<T> > cb2 = cb1.SubBandMatrix(
      0,N/2,0,N,cb1.nlo(),cb1.nhi());
  BandMatrix<T> b3 = b1.SubBandMatrix(0,N,0,N/2,b1.nlo(),b1.nhi());
  BandMatrix<complex<T> > cb3 = cb1.SubBandMatrix(
      0,N,0,N/2,cb1.nlo(),cb1.nhi());
  Matrix<T> a7 = a1;
  Matrix<complex<T> > ca7 = ca1;
  a7.diag().AddToAll(T(10)*N);
  ca7.diag().AddToAll(T(10)*N);

  size_t ntest = dt == tmv::LU ? nsquare : ntot;

  for(size_t i=STARTAT;i<ntest;i++) {
    if (showstartdone) 
      cout<<"Start loop: i = "<<i<<"\nbi = "<<tmv::Type(b[i])<<"  "<<b[i]<<endl;
    const BandMatrixView<T>& bi = b[i];
    const BandMatrixView<complex<T> >& cbi = cb[i];
    bi.SaveDiv();
    cbi.SaveDiv();
    Assert(bi.colsize() == size_t(N),"Band colsize == N");

    Matrix<T> m = bi;
    m.SaveDiv();
    bi.DivideUsing(dt);
    bi.SetDiv();
    m.DivideUsing(dt);
    m.SetDiv();

    ostream* divout = showdiv ? &cout : 0;
    Assert(bi.CheckDecomp(divout),"CheckDecomp");
    T eps = EPS*Norm(m)*Norm(m.Inverse());

    Vector<T> x1 = v1/bi;
    Vector<T> x2 = v1/m;
    if (showacc) {
      cout<<"v/b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(v1)<<endl;
    }
    Assert(Norm(x1-x2) < eps*Norm(v1),"Band v/b");

    if (bi.IsSquare()) {
      x1 = v1%bi;
      x2 = v1%m;
      if (showacc) {
	cout<<"v%b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(v1)<<endl;
      }
      Assert(Norm(x1-x2) < eps*Norm(v1),"Band v%b");
    }

    Matrix<T,ColMajor> binv = bi.Inverse();
    Matrix<T,ColMajor> minv = m.Inverse();
    if (showacc) {
      cout<<"minv = "<<minv<<endl;
      cout<<"binv = "<<binv<<endl;
      cout<<"Norm(minv-binv) = "<<Norm(minv-binv)<<"  "<<eps*Norm(binv)<<endl;
    }
    Assert(Norm(binv-minv) < eps*Norm(binv),"Band Inverse");

    if (showacc) {
      cout<<"b.Det = "<<bi.Det()<<", m.Det = "<<m.Det()<<endl;
      cout<<"abs(bdet-mdet) = "<<abs(bi.Det()-m.Det());
      cout<<"  EPS*abs(mdet) = "<<eps*abs(m.Det())<<endl;
      cout<<"abs(abs(bdet)-abs(mdet)) = "<<abs(abs(bi.Det())-abs(m.Det()));
      cout<<"  EPS*abs(mdet) = "<<eps*abs(m.Det())<<endl;
    }
    if (m.IsSquare())
      Assert(abs(m.Det()-bi.Det()) < eps*abs(m.Det()),"Band Det");
    else
      Assert(abs(abs(m.Det())-abs(bi.Det())) < eps*abs(m.Det()),
	  "Band Det");

    cbi.DivideUsing(dt);
    cbi.SetDiv();
    Assert(cbi.CheckDecomp(divout),"CheckDecomp");

    Matrix<complex<T> > cm(cbi);
    cm.SaveDiv();
    cm.DivideUsing(dt);
    cm.SetDiv();
    T ceps = EPS*Norm(cm)*Norm(cm.Inverse());

    if (showacc) {
      cout<<"cbi.Det = "<<cbi.Det()<<", cm.Det = "<<cm.Det()<<endl;
      cout<<"abs(cbidet-cmdet) = "<<abs(cbi.Det()-cm.Det());
      cout<<"  cbidet/cmdet = "<<cbi.Det()/cm.Det();
      cout<<"  EPS*abs(cmdet) = "<<ceps*abs(cm.Det())<<endl;
      cout<<"abs(abs(bdet)-abs(mdet)) = "<<abs(abs(bi.Det())-abs(m.Det()));
      cout<<"  EPS*abs(mdet) = "<<ceps*abs(m.Det())<<endl;
    }
    if (cm.IsSquare())
      Assert(abs(cbi.Det()-cm.Det()) < ceps*abs(cm.Det()),
	  "Band CDet");
    else
      Assert(abs(abs(cbi.Det())-abs(cm.Det())) < ceps*abs(cm.Det()),
	  "Band CDet");

    Vector<complex<T> > cv(v1 * complex<T>(1,1));
    cv(1) += complex<T>(-1,5);
    cv(2) -= complex<T>(-1,5);

    // test real / complex
    Vector<complex<T> > y1 = v1/cbi;
    Vector<complex<T> > y2 = v1/cm;
    if (showacc) {
      cout<<"v/cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(v1)<<endl;
    }
    Assert(Norm(y1-y2) < ceps*Norm(v1),"Band v/cb");

    // test complex / real
    y1 = cv/bi;
    y2 = cv/m;
    if (showacc) {
      cout<<"cv/b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(cv)<<endl;
    }
    Assert(Norm(y1-y2) < eps*Norm(cv),"Band cv/b");

    // test complex / complex
    y1 = cv/cbi;
    y2 = cv/cm;
    if (showacc) {
      cout<<"cv/cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(cv)<<endl;
    }
    Assert(Norm(y1-y2) < ceps*Norm(cv),"Band cv/cb");

    if (bi.IsSquare()) {
      y1 = v1%cbi;
      y2 = v1%cm;
      if (showacc) {
	cout<<"v%cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(v1)<<endl;
      }
      Assert(Norm(y1-y2) < ceps*Norm(v1),"Band v%cb");

      y1 = cv%bi;
      y2 = cv%m;
      if (showacc) {
	cout<<"cv%b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(cv)<<endl;
      }
      Assert(Norm(y1-y2) < eps*Norm(cv),"Band cv%b");
      y1 = cv%cbi;
      y2 = cv%cm;
      if (showacc) {
	cout<<"cv%cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(cv)<<endl;
      }
      Assert(Norm(y1-y2) < ceps*Norm(cv),"Band cv%cb");
    }

    TestMatrixDivArith<T>(dt,b[i],a1.View(),cbi,ca1.View(),
	"SquareMatrix/Band");
#ifdef XTEST
    TestMatrixDivArith<T>(dt,b[i],a3.View(),cbi,ca3.View(),
	"NonSquareMatrix/Band");
    TestMatrixDivArith<T>(dt,b[i],a4.View(),cbi,ca4.View(),
	"NonSquareMatrix/Band");
    TestMatrixDivArith<T>(dt,b[i],a5.View(),cbi,ca5.View(),
	"DegenerateMatrix/Band");
    TestMatrixDivArith<T>(dt,b[i],a6.View(),cbi,ca6.View(),
	"DegenerateMatrix/Band");
    TestMatrixDivArith<T>(dt,b[i],b1.View(),cbi,cb1.View(),
	"SquareBand/Band");
    TestMatrixDivArith<T>(dt,b[i],b2.View(),cbi,cb2.View(),
	"NonSquareBand/Band");
    TestMatrixDivArith<T>(dt,b[i],b3.View(),cbi,cb3.View(),
	"NonSquareBand/Band");
    TestMatrixDivArith<T>(dt,a7.View(),b[i],ca7.View(),cbi,
	"Band/SquareMatrix");
#endif
  }

  cout<<"BandMatrix<"<<tmv::Type(T())<<"> Division using ";
  cout<<tmv::Text(dt)<<" passed all tests\n";
}

template <class T> void TestAllBandDiv()
{
  TestBandDiv<T>(tmv::LU);
  TestBandDiv<T>(tmv::QR);
  TestBandDiv<T>(tmv::SV);
}

template void TestAllBandDiv<double>();
#ifndef NOFLOAT
template void TestAllBandDiv<float>();
#endif
#ifdef LONGDOUBLE
template void TestAllBandDiv<long double>();
#endif
