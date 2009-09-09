//#define SHOWACC
//#define SHOWTESTS
//#define SHOWCHECK
//#define DOALLARITH

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

#ifdef DOALLARITH
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
  BandMatrix<T,DiagMajor> B4 = tmv::UpperBiDiagMatrix(v1,v2);
  BandMatrix<complex<T>,DiagMajor> CB4 = tmv::UpperBiDiagMatrix(cv1,cv2);
  b.push_back(B4.View());
  cb.push_back(CB4.View());
  BandMatrix<T,DiagMajor> B5 = tmv::LowerBiDiagMatrix(v2,v1);
  BandMatrix<complex<T>,DiagMajor> CB5 = tmv::LowerBiDiagMatrix(cv2,cv1);
  b.push_back(B5.View());
  cb.push_back(CB5.View());
  BandMatrix<T,DiagMajor> B6 = tmv::TriDiagMatrix(v2,v1,v2);
  BandMatrix<complex<T>,DiagMajor> CB6 = tmv::TriDiagMatrix(cv2,cv1,cv2);
  b.push_back(B6.View());
  cb.push_back(CB6.View());
#endif

  BandMatrix<T,ColMajor> B7(a2,6,6);
  BandMatrix<complex<T>,ColMajor> CB7(ca2,6,6);
  b.push_back(B7.SubBandMatrix(0,2*N,0,2*N,3,3,2,2));
  cb.push_back(CB7.SubBandMatrix(0,2*N,0,2*N,3,3,2,2));

#ifdef DOALLARITH
  b.push_back(B7.SubBandMatrix(0,2*N,0,2*N,0,3,2,2));
  cb.push_back(CB7.SubBandMatrix(0,2*N,0,2*N,0,3,2,2));
  b.push_back(B7.SubBandMatrix(0,2*N,0,2*N,3,0,2,2));
  cb.push_back(CB7.SubBandMatrix(0,2*N,0,2*N,3,0,2,2));
  BandMatrix<T,RowMajor> B7c(a2,6,6);
  BandMatrix<complex<T>,RowMajor> CB7c(ca2,6,6);
  b.push_back(B7c.SubBandMatrix(0,2*N,0,2*N,0,3,2,2));
  cb.push_back(CB7c.SubBandMatrix(0,2*N,0,2*N,0,3,2,2));
  b.push_back(B7c.SubBandMatrix(0,2*N,0,2*N,3,0,2,2));
  cb.push_back(CB7c.SubBandMatrix(0,2*N,0,2*N,3,0,2,2));
  BandMatrix<T,DiagMajor> B7b(a2,6,6);
  BandMatrix<complex<T>,DiagMajor> CB7b(ca2,6,6);
  b.push_back(B7b.SubBandMatrix(0,2*N,0,2*N,0,3,2,2));
  cb.push_back(CB7b.SubBandMatrix(0,2*N,0,2*N,0,3,2,2));
  b.push_back(B7b.SubBandMatrix(0,2*N,0,2*N,3,0,2,2));
  cb.push_back(CB7b.SubBandMatrix(0,2*N,0,2*N,3,0,2,2));
#endif

  size_t nsquare = b.size();

  b.push_back(B7.SubBandMatrix(0,N,0,2*N,3,3));
  cb.push_back(CB7.SubBandMatrix(0,N,0,2*N,3,3));
  b.push_back(B7.SubBandMatrix(0,N,0,N+2,4,4));
  cb.push_back(CB7.SubBandMatrix(0,N,0,N+2,4,4));

#ifdef DOALLARITH
  BandMatrix<T,DiagMajor> B8 = tmv::UpperBiDiagMatrix(v1,v1);
  BandMatrix<complex<T>,DiagMajor> CB8 = tmv::UpperBiDiagMatrix(cv1,cv1);
  b.push_back(B8.View());
  cb.push_back(CB8.View());
  BandMatrix<T,DiagMajor> B9 = tmv::LowerBiDiagMatrix(v2,v2);
  BandMatrix<complex<T>,DiagMajor> CB9 = tmv::LowerBiDiagMatrix(cv2,cv2);
  b.push_back(B9.View());
  cb.push_back(CB9.View());
  BandMatrix<T,DiagMajor> B10 = tmv::TriDiagMatrix(v2,v1,v1);
  BandMatrix<complex<T>,DiagMajor> CB10 = tmv::TriDiagMatrix(cv2,cv1,cv1);
  b.push_back(B10.View());
  cb.push_back(CB10.View());
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
    const BandMatrixView<T>& bi = b[i];
    const BandMatrixView<complex<T> >& cbi = cb[i];
    bi.SaveDiv();
    cbi.SaveDiv();
#ifdef SHOWACC
    cout<<"Start loop: i = "<<i<<"\nbi = "<<tmv::Type(bi)<<"  "<<bi<<endl;
#endif
    Assert(bi.colsize() == size_t(N),"Band colsize == N");

    Matrix<T> m = bi;
    m.SaveDiv();
    bi.DivideUsing(dt);
    bi.SetDiv();
    m.DivideUsing(dt);
    m.SetDiv();
#ifdef SHOWCHECK
    Assert(bi.CheckDecomp(&cout),"CheckDecomp");
#else
    Assert(bi.CheckDecomp(),"CheckDecomp");
#endif
    T eps = EPS*Norm(m)*Norm(m.Inverse());

    Vector<T> x1 = v1/bi;
    Vector<T> x2 = v1/m;
#ifdef SHOWACC
    cout<<"v/b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(v1)<<endl;
#endif
    Assert(Norm(x1-x2) < eps*Norm(v1),"Band v/b");

    if (bi.IsSquare()) {
      x1 = v1%bi;
      x2 = v1%m;
#ifdef SHOWACC
      cout<<"v%b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(v1)<<endl;
#endif
      Assert(Norm(x1-x2) < eps*Norm(v1),"Band v%b");
    }

    Matrix<T,ColMajor> binv = bi.Inverse();
    Matrix<T,ColMajor> minv = m.Inverse();
#ifdef SHOWACC
    cout<<"minv = "<<minv<<endl;
    cout<<"binv = "<<binv<<endl;
    cout<<"Norm(minv-binv) = "<<Norm(minv-binv)<<"  "<<eps*Norm(binv)<<endl;
#endif
    Assert(Norm(binv-minv) < eps*Norm(binv),"Band Inverse");

#ifdef SHOWACC
    cout<<"b.Det = "<<bi.Det()<<", m.Det = "<<m.Det()<<endl;
    cout<<"abs(bdet-mdet) = "<<abs(bi.Det()-m.Det());
    cout<<"  EPS*abs(mdet) = "<<eps*abs(m.Det())<<endl;
    cout<<"abs(abs(bdet)-abs(mdet)) = "<<abs(abs(bi.Det())-abs(m.Det()));
    cout<<"  EPS*abs(mdet) = "<<eps*abs(m.Det())<<endl;
#endif
    if (m.IsSquare())
      Assert(abs(m.Det()-bi.Det()) < eps*abs(m.Det()),"Band Det");
    else
      Assert(abs(abs(m.Det())-abs(bi.Det())) < eps*abs(m.Det()),
	  "Band Det");

    cbi.DivideUsing(dt);
    cbi.SetDiv();
#ifdef SHOWCHECK
    Assert(cbi.CheckDecomp(&cout),"CheckDecomp");
#else
    Assert(cbi.CheckDecomp(),"CheckDecomp");
#endif

    Matrix<complex<T> > cm(cbi);
    cm.SaveDiv();
    cm.DivideUsing(dt);
    cm.SetDiv();
    T ceps = EPS*Norm(cm)*Norm(cm.Inverse());

#ifdef SHOWACC
    cout<<"cbi.Det = "<<cbi.Det()<<", cm.Det = "<<cm.Det()<<endl;
    cout<<"abs(cbidet-cmdet) = "<<abs(cbi.Det()-cm.Det());
    cout<<"  cbidet/cmdet = "<<cbi.Det()/cm.Det();
    cout<<"  EPS*abs(cmdet) = "<<ceps*abs(cm.Det())<<endl;
    cout<<"abs(abs(bdet)-abs(mdet)) = "<<abs(abs(bi.Det())-abs(m.Det()));
    cout<<"  EPS*abs(mdet) = "<<ceps*abs(m.Det())<<endl;
#endif
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
#ifdef SHOWACC
    cout<<"v/cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(v1)<<endl;
#endif
    Assert(Norm(y1-y2) < ceps*Norm(v1),"Band v/cb");

    // test complex / real
    y1 = cv/bi;
    y2 = cv/m;
#ifdef SHOWACC
    cout<<"cv/b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(cv)<<endl;
#endif
    Assert(Norm(y1-y2) < eps*Norm(cv),"Band cv/b");

    // test complex / complex
    y1 = cv/cbi;
    y2 = cv/cm;
#ifdef SHOWACC
    cout<<"cv/cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(cv)<<endl;
#endif
    Assert(Norm(y1-y2) < ceps*Norm(cv),"Band cv/cb");

    if (bi.IsSquare()) {
      y1 = v1%cbi;
      y2 = v1%cm;
#ifdef SHOWACC
      cout<<"v%cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(v1)<<endl;
#endif
      Assert(Norm(y1-y2) < ceps*Norm(v1),"Band v%cb");

      y1 = cv%bi;
      y2 = cv%m;
#ifdef SHOWACC
      cout<<"cv%b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(cv)<<endl;
#endif
      Assert(Norm(y1-y2) < eps*Norm(cv),"Band cv%b");
      y1 = cv%cbi;
      y2 = cv%cm;
#ifdef SHOWACC
      cout<<"cv%cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(cv)<<endl;
#endif
      Assert(Norm(y1-y2) < ceps*Norm(cv),"Band cv%cb");
    }

    TestMatrixDivArith<T>(dt,b[i],a1,cbi,ca1,"SquareMatrix/Band");
    TestMatrixDivArith<T>(dt,b[i],a3,cbi,ca3,"NonSquareMatrix/Band");
    TestMatrixDivArith<T>(dt,b[i],a4,cbi,ca4,"NonSquareMatrix/Band");
    TestMatrixDivArith<T>(dt,b[i],a5,cbi,ca5,"DegenerateMatrix/Band");
    TestMatrixDivArith<T>(dt,b[i],a6,cbi,ca6,"DegenerateMatrix/Band");
    TestMatrixDivArith<T>(dt,b[i],b1,cbi,cb1,"SquareBand/Band");
    TestMatrixDivArith<T>(dt,b[i],b2,cbi,cb2,"NonSquareBand/Band");
    TestMatrixDivArith<T>(dt,b[i],b3,cbi,cb3,"NonSquareBand/Band");
    TestMatrixDivArith<T>(dt,a7,b[i],ca7,cbi,"Band/SquareMatrix");
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
