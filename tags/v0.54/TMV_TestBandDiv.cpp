
#define START 0

#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_TestMatrixDivArith.h"

template <class T> inline void TestBandDiv(tmv::DivType dt)
{
  const int N = 10;

  std::vector<tmv::BandMatrixView<T> > b;
  std::vector<tmv::BandMatrixView<std::complex<T> > > cb;

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-2*j;
  a1.diag().AddToAll(T(10)*N);
  tmv::Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+j;
  a2.diag().AddToAll(T(10)*N);

  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);
  tmv::Matrix<std::complex<T> > ca2 = a2 * std::complex<T>(3,-4);

  tmv::Vector<T> v1(N);
  tmv::Vector<T> v2(N-1);
  for (int i=0; i<N; ++i) v1(i) = 16.-3*i; 
  for (int i=0; i<N-1; ++i) v2(i) = -7.+2*i; 
  tmv::Vector<std::complex<T> > cv1 = v1 * std::complex<T>(-1,2);
  tmv::Vector<std::complex<T> > cv2 = v2 * std::complex<T>(-1,2);

  tmv::BandMatrix<T,tmv::RowMajor> B1(a1,3,3);
  tmv::BandMatrix<std::complex<T>,tmv::RowMajor> CB1(ca1,3,3);
  b.push_back(B1.SubBandMatrix(0,N,0,N,1,3));
  cb.push_back(CB1.SubBandMatrix(0,N,0,N,1,3));
  b.push_back(B1.SubBandMatrix(0,N,0,N,3,1));
  cb.push_back(CB1.SubBandMatrix(0,N,0,N,3,1));
  b.push_back(B1.SubBandMatrix(0,N,0,N,0,3));
  cb.push_back(CB1.SubBandMatrix(0,N,0,N,0,3));
  b.push_back(B1.SubBandMatrix(0,N,0,N,3,0));
  cb.push_back(CB1.SubBandMatrix(0,N,0,N,3,0));
  tmv::BandMatrix<T,tmv::ColMajor> B2(a1,3,3);
  tmv::BandMatrix<std::complex<T>,tmv::ColMajor> CB2(ca1,3,3);
  b.push_back(B2.SubBandMatrix(0,N,0,N,1,3));
  cb.push_back(CB2.SubBandMatrix(0,N,0,N,1,3));
  tmv::BandMatrix<T,tmv::DiagMajor> B3(a1,3,3);
  tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> CB3(ca1,3,3);
  b.push_back(B3.SubBandMatrix(0,N,0,N,1,3));
  cb.push_back(CB3.SubBandMatrix(0,N,0,N,1,3));
  b.push_back(B3.SubBandMatrix(0,N,0,N,1,1));
  cb.push_back(CB3.SubBandMatrix(0,N,0,N,1,1));

  tmv::BandMatrix<T,tmv::ColMajor> B8(a2,6,6);
  tmv::BandMatrix<std::complex<T>,tmv::ColMajor> CB8(ca2,6,6);

#ifdef XTEST
  tmv::BandMatrix<T,tmv::DiagMajor> B4 = tmv::UpperBiDiagMatrix(v1,v2);
  tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> CB4 = tmv::UpperBiDiagMatrix(cv1,cv2);
  tmv::BandMatrix<T,tmv::DiagMajor> B5 = tmv::LowerBiDiagMatrix(v2,v1);
  tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> CB5 = tmv::LowerBiDiagMatrix(cv2,cv1);
  tmv::BandMatrix<T,tmv::DiagMajor> B6 = tmv::TriDiagMatrix(v2,v1,v2);
  tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> CB6 = tmv::TriDiagMatrix(cv2,cv1,cv2);
  tmv::BandMatrix<T,tmv::RowMajor> B7(a1,3,N-1);
  tmv::BandMatrix<std::complex<T>,tmv::RowMajor> CB7(ca1,3,N-1);
  tmv::BandMatrix<T,tmv::ColMajor> B7b(a1,3,N-1);
  tmv::BandMatrix<std::complex<T>,tmv::ColMajor> CB7b(ca1,3,N-1);
  tmv::BandMatrix<T,tmv::RowMajor> B7c(a1,3,N-2);
  tmv::BandMatrix<std::complex<T>,tmv::RowMajor> CB7c(ca1,3,N-2);
  tmv::BandMatrix<T,tmv::ColMajor> B7d(a1,3,N-2);
  tmv::BandMatrix<std::complex<T>,tmv::ColMajor> CB7d(ca1,3,N-2);
  tmv::BandMatrix<T,tmv::DiagMajor> B7e(a1,3,N-2);
  tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> CB7e(ca1,3,N-2);
  tmv::BandMatrix<T,tmv::RowMajor> B7f(a1,N-3,N-3);
  tmv::BandMatrix<std::complex<T>,tmv::RowMajor> CB7f(ca1,N-3,N-3);
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
    b.push_back(B7e.View());
    cb.push_back(CB7e.View());
    b.push_back(B7f.View());
    cb.push_back(CB7f.View());
  }

  b.push_back(B8.SubBandMatrix(0,2*N,0,2*N,3,3,2,2));
  cb.push_back(CB8.SubBandMatrix(0,2*N,0,2*N,3,3,2,2));

  tmv::BandMatrix<T,tmv::RowMajor> B8c(a2,6,6);
  tmv::BandMatrix<std::complex<T>,tmv::RowMajor> CB8c(ca2,6,6);
  tmv::BandMatrix<T,tmv::DiagMajor> B8b(a2,6,6);
  tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> CB8b(ca2,6,6);
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
  tmv::BandMatrix<T,tmv::DiagMajor> B9 = tmv::UpperBiDiagMatrix(v1,v1);
  tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> CB9 = tmv::UpperBiDiagMatrix(cv1,cv1);
  tmv::BandMatrix<T,tmv::DiagMajor> B10 = tmv::LowerBiDiagMatrix(v2,v2);
  tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> CB10 = tmv::LowerBiDiagMatrix(cv2,cv2);
  tmv::BandMatrix<T,tmv::DiagMajor> B11 = tmv::TriDiagMatrix(v2,v1,v1);
  tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> CB11 = tmv::TriDiagMatrix(cv2,cv1,cv1);
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

  tmv::Matrix<T> a3 = a1.Cols(0,N/2);
  tmv::Matrix<std::complex<T> > ca3 = ca1.Cols(0,N/2);
  tmv::Matrix<T> a4 = a1.Rows(0,N/2);
  tmv::Matrix<std::complex<T> > ca4 = ca1.Rows(0,N/2);
  tmv::Matrix<T> a5 = a1.Cols(0,0);
  tmv::Matrix<std::complex<T> > ca5 = ca1.Cols(0,0);
  tmv::Matrix<T> a6 = a1.Rows(0,0);
  tmv::Matrix<std::complex<T> > ca6 = ca1.Rows(0,0);
  tmv::BandMatrix<T> b1 = B1;
  tmv::BandMatrix<std::complex<T> > cb1 = b1 * std::complex<T>(2.,-3.);
  tmv::BandMatrix<T> b2 = b1.SubBandMatrix(0,N/2,0,N,b1.nlo(),b1.nhi());
  tmv::BandMatrix<std::complex<T> > cb2 = cb1.SubBandMatrix(
      0,N/2,0,N,cb1.nlo(),cb1.nhi());
  tmv::BandMatrix<T> b3 = b1.SubBandMatrix(0,N,0,N/2,b1.nlo(),b1.nhi());
  tmv::BandMatrix<std::complex<T> > cb3 = cb1.SubBandMatrix(
      0,N,0,N/2,cb1.nlo(),cb1.nhi());
  tmv::Matrix<T> a7 = a1;
  tmv::Matrix<std::complex<T> > ca7 = ca1;
  a7.diag().AddToAll(T(10)*N);
  ca7.diag().AddToAll(T(10)*N);

  size_t ntest = dt == tmv::LU ? nsquare : ntot;

  for(size_t i=START;i<ntest;i++) {
    if (showstartdone) 
      std::cout<<"Start loop: i = "<<i<<"\nbi = "<<tmv::Type(b[i])<<"  "<<b[i]<<std::endl;
    const tmv::BandMatrixView<T>& bi = b[i];
    const tmv::BandMatrixView<std::complex<T> >& cbi = cb[i];
    bi.SaveDiv();
    cbi.SaveDiv();
    Assert(bi.colsize() == size_t(N),"Band colsize == N");

    tmv::Matrix<T> m = bi;
    m.SaveDiv();
    bi.DivideUsing(dt);
    bi.SetDiv();
    m.DivideUsing(dt);
    m.SetDiv();

    std::ostream* divout = showdiv ? &std::cout : 0;
    Assert(bi.CheckDecomp(divout),"CheckDecomp");
    T eps = EPS*Norm(m)*Norm(m.Inverse());

    tmv::Vector<T> x1 = v1/bi;
    tmv::Vector<T> x2 = v1/m;
    if (showacc) {
      std::cout<<"v/b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(v1)<<std::endl;
    }
    Assert(Norm(x1-x2) < eps*Norm(v1),"Band v/b");

    if (bi.IsSquare()) {
      x1 = v1%bi;
      x2 = v1%m;
      if (showacc) {
	std::cout<<"v%b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(v1)<<std::endl;
      }
      Assert(Norm(x1-x2) < eps*Norm(v1),"Band v%b");
    }

    tmv::Matrix<T,tmv::ColMajor> binv = bi.Inverse();
    tmv::Matrix<T,tmv::ColMajor> minv = m.Inverse();
    if (showacc) {
      std::cout<<"minv = "<<minv<<std::endl;
      std::cout<<"binv = "<<binv<<std::endl;
      std::cout<<"Norm(minv-binv) = "<<Norm(minv-binv)<<"  "<<eps*Norm(binv)<<std::endl;
    }
    Assert(Norm(binv-minv) < eps*Norm(binv),"Band Inverse");

    if (m.IsSquare()) {
      if (showacc) {
	std::cout<<"b.Det = "<<bi.Det()<<", m.Det = "<<m.Det()<<std::endl;
	std::cout<<"abs(bdet-mdet) = "<<std::abs(bi.Det()-m.Det());
	std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(m.Det())<<std::endl;
	std::cout<<"abs(abs(bdet)-abs(mdet)) = "<<std::abs(std::abs(bi.Det())-std::abs(m.Det()));
	std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(m.Det())<<std::endl;
      }
      Assert(std::abs(m.Det()-bi.Det()) < eps*std::abs(m.Det()),"Band Det");
    }

    cbi.DivideUsing(dt);
    cbi.SetDiv();
    Assert(cbi.CheckDecomp(divout),"CheckDecomp");

    tmv::Matrix<std::complex<T> > cm(cbi);
    cm.SaveDiv();
    cm.DivideUsing(dt);
    cm.SetDiv();
    T ceps = EPS*Norm(cm)*Norm(cm.Inverse());

    if (cm.IsSquare()) {
      if (showacc) {
	std::cout<<"cbi.Det = "<<cbi.Det()<<", cm.Det = "<<cm.Det()<<std::endl;
	std::cout<<"abs(cbidet-cmdet) = "<<std::abs(cbi.Det()-cm.Det());
	std::cout<<"  cbidet/cmdet = "<<cbi.Det()/cm.Det();
	std::cout<<"  EPS*abs(cmdet) = "<<ceps*std::abs(cm.Det())<<std::endl;
	std::cout<<"abs(abs(bdet)-abs(mdet)) = "<<std::abs(std::abs(bi.Det())-std::abs(m.Det()));
	std::cout<<"  EPS*abs(mdet) = "<<ceps*std::abs(m.Det())<<std::endl;
      }
      Assert(abs(cbi.Det()-cm.Det()) < ceps*std::abs(cm.Det()),
	  "Band CDet");
    }

    tmv::Vector<std::complex<T> > cv(v1 * std::complex<T>(1,1));
    cv(1) += std::complex<T>(-1,5);
    cv(2) -= std::complex<T>(-1,5);

    // test real / complex
    tmv::Vector<std::complex<T> > y1 = v1/cbi;
    tmv::Vector<std::complex<T> > y2 = v1/cm;
    if (showacc) {
      std::cout<<"v/cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(v1)<<std::endl;
    }
    Assert(Norm(y1-y2) < ceps*Norm(v1),"Band v/cb");

    // test complex / real
    y1 = cv/bi;
    y2 = cv/m;
    if (showacc) {
      std::cout<<"cv/b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(cv)<<std::endl;
    }
    Assert(Norm(y1-y2) < eps*Norm(cv),"Band cv/b");

    // test complex / complex
    y1 = cv/cbi;
    y2 = cv/cm;
    if (showacc) {
      std::cout<<"cv/cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(cv)<<std::endl;
    }
    Assert(Norm(y1-y2) < ceps*Norm(cv),"Band cv/cb");

    if (bi.IsSquare()) {
      y1 = v1%cbi;
      y2 = v1%cm;
      if (showacc) {
	std::cout<<"v%cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(v1)<<std::endl;
      }
      Assert(Norm(y1-y2) < ceps*Norm(v1),"Band v%cb");

      y1 = cv%bi;
      y2 = cv%m;
      if (showacc) {
	std::cout<<"cv%b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(cv)<<std::endl;
      }
      Assert(Norm(y1-y2) < eps*Norm(cv),"Band cv%b");
      y1 = cv%cbi;
      y2 = cv%cm;
      if (showacc) {
	std::cout<<"cv%cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(cv)<<std::endl;
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

  std::cout<<"BandMatrix<"<<tmv::Type(T())<<"> Division using ";
  std::cout<<tmv::Text(dt)<<" passed all tests\n";
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
