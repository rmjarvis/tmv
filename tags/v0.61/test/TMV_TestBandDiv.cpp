#ifdef NDEBUG
#undef NDEBUG
#endif

#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"
#include "TMV_TestBandArith.h"

template <class T> void TestBandDiv(tmv::DivType dt)
{
  const int N = 10;

  std::vector<tmv::BandMatrixView<T> > b;
  std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
  MakeBandList(b,cb);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-2*j;
  a1.diag().AddToAll(T(10)*N);
  a1 /= T(10);

  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

  static tmv::Vector<T> v1(N);
  static tmv::Vector<T> v2(N-1);
  for (int i=0; i<N; ++i) v1(i) = 16.-3*i; 
  for (int i=0; i<N-1; ++i) v2(i) = -6.+i; 
  static tmv::Vector<std::complex<T> > cv1(N);
  static tmv::Vector<std::complex<T> > cv2(N-1);
  for (int i=0; i<N; ++i) cv1(i) = std::complex<T>(16.-3*i,i+4.); 
  for (int i=0; i<N-1; ++i) cv2(i) = std::complex<T>(2*i-3.,-6.+i); 

  tmv::Matrix<T> a3 = a1.Cols(0,N/2);
  tmv::Matrix<std::complex<T> > ca3 = ca1.Cols(0,N/2);
  tmv::Matrix<T> a4 = a1.Rows(0,N/2);
  tmv::Matrix<std::complex<T> > ca4 = ca1.Rows(0,N/2);
  tmv::Matrix<T> a5 = a1.Cols(0,0);
  tmv::Matrix<std::complex<T> > ca5 = ca1.Cols(0,0);
  tmv::Matrix<T> a6 = a1.Rows(0,0);
  tmv::Matrix<std::complex<T> > ca6 = ca1.Rows(0,0);
  tmv::Matrix<T> a7 = a1;
  tmv::Matrix<std::complex<T> > ca7 = ca1;
  a7.diag().AddToAll(T(10)*N);
  ca7.diag().AddToAll(T(10)*N);

  for(size_t i=START;i<b.size();i++) {
    if (showstartdone) 
      std::cout<<"Start loop: i = "<<i<<"\nbi = "<<tmv::Type(b[i])<<"  "<<b[i]<<std::endl;
    const tmv::BandMatrixView<T>& bi = b[i];
    const tmv::BandMatrixView<std::complex<T> >& cbi = cb[i];
    if (dt == tmv::LU && !bi.IsSquare()) continue;

    bi.SaveDiv();
    cbi.SaveDiv();

    tmv::Matrix<T> m(bi);
    m.SaveDiv();
    bi.DivideUsing(dt);
    bi.SetDiv();
    m.DivideUsing(dt);
    m.SetDiv();

    std::ostream* divout = showdiv ? &std::cout : 0;
    Assert(bi.CheckDecomp(divout),"CheckDecomp");
    T eps = m.rowsize()*EPS*Norm(m)*Norm(m.Inverse());

    if (bi.colsize() == size_t(N)) {
      tmv::Vector<T> x1 = v1/bi;
      tmv::Vector<T> x2 = v1/m;
      if (showacc) {
	std::cout<<"v/b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(x1)<<std::endl;
      }
      Assert(Norm(x1-x2) < eps*Norm(x1),"Band v/b");
    }

    if (bi.rowsize() == size_t(N)) {
      tmv::Vector<T> x1 = v1%bi;
      tmv::Vector<T> x2 = v1%m;
      if (showacc) {
	std::cout<<"v%b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(x1)<<std::endl;
      }
      Assert(Norm(x1-x2) < eps*Norm(x1),"Band v%b");
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
      Assert(std::abs(m.Det()-bi.Det()) < eps*std::abs(m.Det()+m.Norm()),"Band Det");
      T msign, bsign;
      Assert(std::abs(m.LogDet(&msign)-bi.LogDet(&bsign)) < N*eps,"Band LogDet");
      Assert(std::abs(msign-bsign) < N*eps,"Band LogDet - sign");
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
      Assert(std::abs(cbi.Det()-cm.Det()) < ceps*std::abs(cm.Det()+cm.Norm()),
	  "Band CDet");
      std::complex<T> cmsign, cbsign;
      Assert(std::abs(cm.LogDet(&cmsign)-cbi.LogDet(&cbsign)) < N*eps,
	  "Band CLogDet");
      Assert(std::abs(cmsign-cbsign) < N*eps,"Band CLogDet - sign");
    }

    tmv::Vector<std::complex<T> > cv(v1 * std::complex<T>(1,1));
    cv(1) += std::complex<T>(-1,5);
    cv(2) -= std::complex<T>(-1,5);

    if (m.colsize() == size_t(N)) {
      // test real / complex
      tmv::Vector<std::complex<T> > y1 = v1/cbi;
      tmv::Vector<std::complex<T> > y2 = v1/cm;
      if (showacc) {
	std::cout<<"v/cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(y1)<<std::endl;
      }
      Assert(Norm(y1-y2) < ceps*Norm(y1),"Band v/cb");

      // test complex / real
      y1 = cv/bi;
      y2 = cv/m;
      if (showacc) {
	std::cout<<"cv/b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(y1)<<std::endl;
      }
      Assert(Norm(y1-y2) < eps*Norm(y1),"Band cv/b");

      // test complex / complex
      y1 = cv/cbi;
      y2 = cv/cm;
      if (showacc) {
	std::cout<<"cv/cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(y1)<<std::endl;
      }
      Assert(Norm(y1-y2) < ceps*Norm(y1),"Band cv/cb");
    }

    if (bi.rowsize() == size_t(N)) {
      tmv::Vector<std::complex<T> > y1 = v1%cbi;
      tmv::Vector<std::complex<T> > y2 = v1%cm;
      if (showacc) {
	std::cout<<"v%cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(y1)<<std::endl;
      }
      Assert(Norm(y1-y2) < ceps*Norm(y1),"Band v%cb");

      y1 = cv%bi;
      y2 = cv%m;
      if (showacc) {
	std::cout<<"cv%b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(y1)<<std::endl;
      }
      Assert(Norm(y1-y2) < eps*Norm(y1),"Band cv%b");
      y1 = cv%cbi;
      y2 = cv%cm;
      if (showacc) {
	std::cout<<"cv%cb: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(y1)<<std::endl;
      }
      Assert(Norm(y1-y2) < ceps*Norm(y1),"Band cv%cb");
    }

  }

  TestBandDiv_A<T>(dt);
  TestBandDiv_B<T>(dt);
  TestBandDiv_C<T>(dt);
  TestBandDiv_D<T>(dt);
  std::cout<<"BandMatrix<"<<tmv::Type(T())<<"> Division using ";
  std::cout<<tmv::Text(dt)<<" passed all tests\n";
}

template <class T, tmv::StorageType stor> static void TestBandDecomp()
{
  for (int mattype = 0; mattype < 5; mattype++) {
    // mattype = 0  is Square
    // mattype = 1  is NonSquare short
    // mattype = 2  is NonSquare tall
    // mattype = 3  is TriDiag
    // mattype = 4  is Singular

    int M = 8;
    const int N = 8;
    int nlo = 2;
    int nhi = 3;
    if (mattype == 1) M = 10;
    else if (mattype == 2) M = 15;
    else if (mattype == 3) nlo = nhi = 1;

    tmv::BandMatrix<T,stor> m(M,N,nlo,nhi);
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
      if (i<=j+nlo && j<=i+nhi) m(i,j) = T(2.+4*i-5*j);
    if (mattype != 4) {
      m(0,0) = T(14);
      m(1,0) = T(-2);
      m(4,5) = T(7);
      m(2,2) = T(-10);
      m.diag() *= T(30);
    }

    tmv::BandMatrix<std::complex<T>,stor> c(M,N,nlo,nhi);
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
      if (i<=j+nlo && j<=i+nhi) c(i,j) = std::complex<T>(2.+4*i-5*j,3.-i);
    if (mattype != 4) {
      c(0,0) = T(14);
      c(1,0) = T(-2);
      c(4,5) = T(7);
      c(2,2) = T(-10);
      c.diag() *= T(30);
    }

    T eps = EPS;
    T ceps = EPS;
    if (mattype != 3) {
      eps *= Norm(m) * Norm(m.Inverse());
      ceps *= Norm(c) * Norm(c.Inverse());
    } else {
      eps *= T(10);
      ceps *= T(10);
    }

    // LU Decomposition
    if (mattype == 0) {
      m.DivideUsing(tmv::LU);
      m.SetDiv();
      tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.LUD().GetL();
      tmv::UpperTriMatrix<T> U = m.LUD().GetU();
      const int* p = m.LUD().GetP();
      tmv::Matrix<T> PLU = L*U;
      PLU.ReversePermuteRows(p);
      if (m.LUD().IsTrans()) PLU.TransposeSelf();
      Assert(Norm(m-PLU) < eps*Norm(m),"Band LU");

      tmv::LowerTriMatrix<T,tmv::UnitDiag> L2(M);
      tmv::BandMatrix<T,stor> U2(M,M,0,nlo+nhi);
      int p2[N];
      LU_Decompose(m,L2.View(),U2.View(),p2);
      PLU = L2*U2;
      PLU.ReversePermuteRows(p2);
      Assert(Norm(m-PLU) < eps*Norm(m),"Band LU2");

      c.DivideUsing(tmv::LU);
      c.SetDiv();
      tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = c.LUD().GetL();
      tmv::UpperTriMatrix<std::complex<T> > cU = c.LUD().GetU();
      p = c.LUD().GetP();
      tmv::Matrix<std::complex<T> > cPLU = cL*cU;
      cPLU.ReversePermuteRows(p);
      if (c.LUD().IsTrans()) cPLU.TransposeSelf();
      Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU");

      tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL2(M);
      tmv::BandMatrix<std::complex<T>,stor> cU2(M,M,0,nlo+nhi);
      LU_Decompose(c,cL2.View(),cU2.View(),p2);
      cPLU = cL2*cU2;
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU2");

      LU_Decompose(c,cL2.Conjugate(),cU2.View(),p2);
      cPLU = cL2.Conjugate()*cU2;
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU3");

      LU_Decompose(c,cL2.View(),cU2.Conjugate(),p2);
      cPLU = cL2*cU2.Conjugate();
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU4");

      LU_Decompose(c,cL2.Conjugate(),cU2.Conjugate(),p2);
      cPLU = cL2.Conjugate()*cU2.Conjugate();
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU5");

      LU_Decompose(c.Conjugate(),cL2.View(),cU2.View(),p2);
      cPLU = cL2*cU2;
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c.Conjugate()-cPLU) < ceps*Norm(c),"Band C LU6");

      LU_Decompose(c.Conjugate(),cL2.Conjugate(),cU2.View(),p2);
      cPLU = cL2.Conjugate()*cU2;
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c.Conjugate()-cPLU) < ceps*Norm(c),"Band C LU7");

      LU_Decompose(c.Conjugate(),cL2.View(),cU2.Conjugate(),p2);
      cPLU = cL2*cU2.Conjugate();
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c.Conjugate()-cPLU) < ceps*Norm(c),"Band C LU8");

      LU_Decompose(c.Conjugate(),cL2.Conjugate(),cU2.Conjugate(),p2);
      cPLU = cL2.Conjugate()*cU2.Conjugate();
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c.Conjugate()-cPLU) < ceps*Norm(c),"Band C LU9");
    }
  
    // QR Decomposition
    if (mattype != 4) {
      m.DivideUsing(tmv::QR);
      m.SetDiv();
      tmv::Matrix<T> Q = m.QRD().GetQ();
      tmv::BandMatrix<T> R = m.QRD().GetR();
      tmv::Matrix<T> QR = Q*R;
      if (m.QRD().IsTrans()) QR.TransposeSelf();
      Assert(Norm(m-QR) < eps*Norm(m),"Band QR");

      QR_Decompose(m,Q.View(),R.View());
      QR = Q*R;
      Assert(Norm(m-QR) < eps*Norm(m),"Band QR2");

      tmv::BandMatrix<T,stor> R2(N,N,0,nlo+nhi);
      QR_Decompose(m,R2.View());
      Assert(Norm(R-R2) < eps*Norm(R),"Band QR3");

      c.DivideUsing(tmv::QR);
      c.SetDiv();
      tmv::Matrix<std::complex<T> > cQ = c.QRD().GetQ();
      tmv::BandMatrix<std::complex<T> > cR = c.QRD().GetR();
      tmv::Matrix<std::complex<T> > cQR = cQ*cR;
      if (c.QRD().IsTrans()) cQR.TransposeSelf();
      Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR");

      QR_Decompose(c,cQ.View(),cR.View());
      cQR = cQ*cR;
      Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR2");

      tmv::BandMatrix<std::complex<T>,stor> cR2(N,N,0,nlo+nhi);
      QR_Decompose(c,cR2.View());
      Assert(Norm(cR-cR2) < ceps*Norm(cR),"Band C QR3");

      QR_Decompose(c,cQ.Conjugate(),cR2.View());
      cQR = cQ.Conjugate()*cR2;
      Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR4");

      QR_Decompose(c,cQ.View(),cR2.Conjugate());
      cQR = cQ*cR2.Conjugate();
      Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR5");

      QR_Decompose(c,cQ.Conjugate(),cR2.Conjugate());
      cQR = cQ.Conjugate()*cR2.Conjugate();
      Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR6");

      QR_Decompose(c,cR2.View());
      Assert(Norm(cR-cR2) < ceps*Norm(cR),"Band C QR7");

      QR_Decompose(c,cR2.Conjugate());
      Assert(Norm(cR-cR2.Conjugate()) < ceps*Norm(cR),"Band C QR8");
    
      QR_Decompose(c.Conjugate(),cQ.View(),cR2.View());
      cQR = cQ*cR2;
      Assert(Norm(c.Conjugate()-cQR) < ceps*Norm(c),"Band C QR9");

      QR_Decompose(c.Conjugate(),cR2.View());
      Assert(Norm(cR.Conjugate()-cR2) < ceps*Norm(cR),"Band C QR10");

      QR_Decompose(c.Conjugate(),cQ.Conjugate(),cR2.View());
      cQR = cQ.Conjugate()*cR2;
      Assert(Norm(c.Conjugate()-cQR) < ceps*Norm(c),"Band C QR11");

      QR_Decompose(c.Conjugate(),cQ.View(),cR2.Conjugate());
      cQR = cQ*cR2.Conjugate();
      Assert(Norm(c.Conjugate()-cQR) < ceps*Norm(c),"Band C QR12");

      QR_Decompose(c.Conjugate(),cQ.Conjugate(),cR2.Conjugate());
      cQR = cQ.Conjugate()*cR2.Conjugate();
      Assert(Norm(c.Conjugate()-cQR) < ceps*Norm(c),"Band C QR13");

      QR_Decompose(c.Conjugate(),cR2.View());
      Assert(Norm(cR.Conjugate()-cR2) < ceps*Norm(cR),"Band C QR14");

      QR_Decompose(c.Conjugate(),cR2.Conjugate());
      Assert(Norm(cR.Conjugate()-cR2.Conjugate()) < ceps*Norm(cR),"Band C QR15");
    }
    
    // SV Decomposition
    {
      m.DivideUsing(tmv::SV);
      m.SaveDiv();
      m.SetDiv();
      tmv::Matrix<T> U = m.SVD().GetU();
      tmv::DiagMatrix<T> S = m.SVD().GetS();
      tmv::Matrix<T> V = m.SVD().GetV();
      Assert(Norm(m-U*S*V) < eps*Norm(m),"SV");

      tmv::Matrix<T> U2(M,N);
      tmv::DiagMatrix<T> S2(N);
      tmv::Matrix<T> V2(N,N);
      SV_Decompose(m,U2.View(),S2.View(),V2.View());
      Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"SV2");

      SV_Decompose(m,S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"SV3");
      SV_Decompose(m,U2.View(),S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"SV4 S");
      Assert(Norm(U2*S2*S2*U2.Transpose()-m*m.Transpose()) < 
	  eps*Norm(m)*Norm(m),"SV3 U");
      SV_Decompose(m,S2.View(),V2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"SV5 S");
      Assert(Norm(V2.Transpose()*S2*S2*V2-m.Transpose()*m) < 
	  eps*Norm(m)*Norm(m),"SV5 V");

      c.DivideUsing(tmv::SV);
      c.SaveDiv();
      c.SetDiv();
      tmv::Matrix<std::complex<T> > cU = c.SVD().GetU();
      S = c.SVD().GetS();
      tmv::Matrix<std::complex<T> > cV = c.SVD().GetV();
      Assert(Norm(c-cU*S*cV) < ceps*Norm(c),"C SV");

      tmv::Matrix<std::complex<T> > cU2(M,N);
      tmv::Matrix<std::complex<T> > cV2(N,N);
      SV_Decompose(c,cU2.View(),S2.View(),cV2.View());
      Assert(Norm(c-cU2*S2*cV2) < eps*Norm(c),"C SV2");

      SV_Decompose(c,S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV3");
      SV_Decompose(c,cU2.View(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV4 S");
      Assert(Norm(cU2*S2*S2*cU2.Adjoint()-c*c.Adjoint()) < 
	  ceps*Norm(c)*Norm(c),"C SV4 U");
      SV_Decompose(c,S2.View(),cV2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV5 S");
      Assert(Norm(cV2.Adjoint()*S2*S2*cV2-c.Adjoint()*c) < 
	  ceps*Norm(c)*Norm(c),"C SV5 V");

      SV_Decompose(c,cU2.Conjugate(),S2.View(),cV2.View());
      Assert(Norm(c-cU2.Conjugate()*S2*cV2) < eps*Norm(c),"C SV6");
      SV_Decompose(c,cU2.View(),S2.View(),cV2.Conjugate());
      Assert(Norm(c-cU2*S2*cV2.Conjugate()) < eps*Norm(c),"C SV7");
      SV_Decompose(c,cU2.Conjugate(),S2.View(),cV2.Conjugate());
      Assert(Norm(c-cU2.Conjugate()*S2*cV2.Conjugate()) < eps*Norm(c),"C SV8");

      SV_Decompose(c,cU2.View(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV9 S");
      Assert(Norm(cU2*S2*S2*cU2.Adjoint()-c*c.Adjoint()) < 
	  ceps*Norm(c)*Norm(c),"C SV9 U");
      SV_Decompose(c,S2.View(),cV2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV10 S");
      Assert(Norm(cV2.Adjoint()*S2*S2*cV2-c.Adjoint()*c) < 
	  ceps*Norm(c)*Norm(c),"C SV10 V");

      SV_Decompose(c.Conjugate(),cU2.View(),S2.View(),cV2.View());
      Assert(Norm(c.Conjugate()-cU2*S2*cV2) < eps*Norm(c),"C SV11");
      SV_Decompose(c.Conjugate(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV12");
      SV_Decompose(c.Conjugate(),cU2.View(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV13 S");
      Assert(Norm(cU2*S2*S2*cU2.Adjoint()-c.Conjugate()*c.Transpose()) < 
	  ceps*Norm(c)*Norm(c),"C SV13 U");
      SV_Decompose(c.Conjugate(),S2.View(),cV2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV14 S");
      Assert(Norm(cV2.Adjoint()*S2*S2*cV2-c.Transpose()*c.Conjugate()) < 
	  ceps*Norm(c)*Norm(c),"C SV14 V");

      SV_Decompose(c.Conjugate(),cU2.Conjugate(),S2.View(),cV2.View());
      Assert(Norm(c.Conjugate()-cU2.Conjugate()*S2*cV2) < eps*Norm(c),"C SV15");
      SV_Decompose(c.Conjugate(),cU2.View(),S2.View(),cV2.Conjugate());
      Assert(Norm(c.Conjugate()-cU2*S2*cV2.Conjugate()) < eps*Norm(c),"C SV16");
      SV_Decompose(c.Conjugate(),cU2.Conjugate(),S2.View(),cV2.Conjugate());
      Assert(Norm(c.Conjugate()-cU2.Conjugate()*S2*cV2.Conjugate()) < 
	  eps*Norm(c),"C SV17");

      SV_Decompose(c.Conjugate(),cU2.Conjugate(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV18 S");
      Assert(Norm(cU2.Conjugate()*S2*S2*cU2.Transpose()-
	    c.Conjugate()*c.Transpose())<ceps*Norm(c)*Norm(c),"C SV18 U");
      SV_Decompose(c.Conjugate(),S2.View(),cV2.Conjugate());
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV19 S");
      Assert(Norm(cV2.Transpose()*S2*S2*cV2.Conjugate()-
	    c.Transpose()*c.Conjugate())<ceps*Norm(c)*Norm(c),"C SV19 V");
    }
  }
}

template <class T> void TestAllBandDiv()
{
  TestBandDecomp<T,tmv::ColMajor>();
  TestBandDecomp<T,tmv::RowMajor>();
  TestBandDecomp<T,tmv::DiagMajor>();
  TestBandDiv<T>(tmv::LU);
  TestBandDiv<T>(tmv::QR);
  TestBandDiv<T>(tmv::SV);
}

#ifdef INST_DOUBLE
template void TestAllBandDiv<double>();
#endif
#ifdef INST_FLOAT
template void TestAllBandDiv<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllBandDiv<long double>();
#endif
#ifdef INST_INT
template void TestAllBandDiv<int>();
#endif

