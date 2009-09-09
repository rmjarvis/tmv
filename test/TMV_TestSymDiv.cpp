#ifdef NDEBUG
#undef NDEBUG
#endif

#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"
#include "TMV_Band.h"
#include "TMV_SymBand.h"
#include "TMV_TestSymArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> static bool IsPosDef(const tmv::GenSymMatrix<T>& m)
{
  try {
    tmv::HermMatrix<T> m2 = m;
    CH_Decompose(m2.View());
  }
  catch (tmv::NonPosDef) {
#ifdef XTEST
    if (showacc)
      std::cout<<"caught\n";
    for(size_t i=1;i<=m.size();i++) {
      T d = tmv::Matrix<T>(m.SubSymMatrix(0,i)).Det();
      if (showacc) 
	std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
      if (!(tmv::REAL(d) > 0)) {
	return false;
      }
    }
    std::cout<<"m = "<<Type(m)<<"  "<<m<<std::endl;
    for(size_t i=1;i<=m.size();i++) {
      T d = tmv::Matrix<T>(m.SubSymMatrix(0,i)).Det();
      std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
    }
    Assert(tmv::FALSE,"Caught NonPosDef, but no determinants of sub-blocks are negative");
#endif
    return false;
  }
#ifdef XTEST
  if (showacc)
    std::cout<<"not caught\n";
  for(size_t i=1;i<=m.size();i++) {
    T d = tmv::Matrix<T>(m.SubSymMatrix(0,i)).Det();
    if (showacc) 
      std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
    if (tmv::REAL(d) < 0) {
      std::cout<<"m = "<<Type(m)<<"  "<<m<<std::endl;
      std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
      Assert(tmv::FALSE,"Didn't catch NonPosDef, but determinant of sub-block is negative");
    }
  }
#endif
  return true;
}

template <class T> void TestSymDiv(tmv::DivType dt, PosDefCode pdc)
{
  const int N = 10;

  std::vector<tmv::SymMatrixView<T> > s;
  std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
  MakeSymList(s,cs,pdc);

  tmv::Vector<T> v1(N);
  tmv::Vector<T> v2(N-1);
  for (int i=0; i<N; ++i) v1(i) = 16.-3*i; 
  for (int i=0; i<N-1; ++i) v2(i) = -7.+2*i; 

  std::ostream* checkout = showdiv ? &std::cout : 0;

  for(size_t i=START;i<s.size();i++) {
    tmv::SymMatrixView<T> si = s[i];
    tmv::SymMatrixView<std::complex<T> > csi = cs[i];
    if (dt == tmv::CH && csi.issym()) continue;
    si.SaveDiv();
    csi.SaveDiv();
    if (showstartdone)
      std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::Type(si)<<"  "<<si<<std::endl;

    Assert(IsPosDef(si) == (pdc==PosDef),"IsPosDef");
    if (csi.isherm()) {
      Assert(IsPosDef(csi) == (pdc==PosDef),"IsPosDef");
    }

    tmv::Matrix<T> m(si);
    m.SaveDiv();
    if (dt == tmv::CH) m.DivideUsing(tmv::LU);
    else m.DivideUsing(dt);
    m.SetDiv();
    Assert(m.CheckDecomp(checkout),"CheckDecomp m"); 
    T eps = EPS;
    if (pdc == Sing) eps *= 1000;
    else eps *= Norm(m)*Norm(m.Inverse());
    si.DivideUsing(dt);
    si.SetDiv();
    if (si.isherm()) {
      tmv::HermMatrix<T> six = si;
      six.DivideUsing(dt);
      six.SetDiv();
      Assert(six.CheckDecomp(checkout),"CheckDecomp six(herm)"); 
    } else {
      tmv::SymMatrix<T> six = si;
      six.DivideUsing(dt);
      six.SetDiv();
      Assert(six.CheckDecomp(checkout),"CheckDecomp six(sym)"); 
    }
    Assert(si.CheckDecomp(checkout),"CheckDecomp si"); 

    tmv::Vector<T> x1 = v1/si;
    tmv::Vector<T> x2 = v1/m;
    if (showacc) {
      std::cout<<"v1 = "<<v1<<std::endl;
      std::cout<<"v1/si = "<<x1<<std::endl;
      std::cout<<"v1/m = "<<x2<<std::endl;
      std::cout<<"si*x1 = "<<si*x1<<std::endl;
      std::cout<<"m*x2 = "<<m*x2<<std::endl;
      std::cout<<"v/b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(x1)<<std::endl;
    }
    Assert(Norm(x1-x2) < eps*Norm(x1),"Sym v/b");

    x1 = v1%si;
    x2 = v1%m;
    if (showacc)
      std::cout<<"v%b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(x1)<<std::endl;
    Assert(Norm(x1-x2) < eps*Norm(x1),"Sym v%b");

    tmv::Matrix<T,tmv::ColMajor> sinv = si.Inverse();
    tmv::Matrix<T,tmv::ColMajor> minv = m.Inverse();
    if (showacc) {
      std::cout<<"sinv = "<<sinv<<std::endl;
      std::cout<<"minv = "<<minv<<std::endl;
      std::cout<<"Norm(minv-sinv) = "<<Norm(minv-sinv)<<"  "<<eps*Norm(sinv)<<std::endl;
    }
    Assert(Norm(sinv-minv) < eps*Norm(sinv),"Sym Inverse");

    if (pdc != Sing) {
      tmv::Matrix<T> m2 = m;
      if (showacc) {
	std::cout<<"si.Det = "<<si.Det()<<", m.Det = "<<m2.Det()<<std::endl;
	std::cout<<"abs(sdet-mdet) = "<<std::abs(si.Det()-m2.Det());
	std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(m2.Det())<<std::endl;
	std::cout<<"abs(abs(sdet)-abs(mdet)) = "<<
	  std::abs(std::abs(si.Det())-std::abs(m2.Det()));
	std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(m2.Det())<<std::endl;
	std::cout<<"m.LogDet() = "<<m2.LogDet()<<std::endl;
	std::cout<<"si.LogDet() = "<<si.LogDet()<<std::endl;
	std::cout<<"abs(diff) = "<<std::abs(m2.LogDet()-si.LogDet())<<"   ";
	std::cout<<"eps = "<<eps<<std::endl;
	T msign,ssign;
	m.LogDet(&msign);
	si.LogDet(&ssign);
	std::cout<<"m sign = "<<msign<<std::endl;
	std::cout<<"si sign = "<<ssign<<std::endl;
	std::cout<<"abs(diff) = "<<std::abs(msign-ssign)<<"   ";
	std::cout<<"eps = "<<eps<<std::endl;
      }
      Assert(std::abs(m2.Det()-si.Det()) < eps*std::abs(m2.Det()),"Sym Det");
      T msign,ssign;
      Assert(std::abs(m2.LogDet(&msign)-si.LogDet(&ssign)) < 10*N*eps,"Sym LogDet");
      Assert(std::abs(msign-ssign) < 10*N*eps,"Sym LogDet - sign");
    }

    tmv::Matrix<std::complex<T> > cm(csi);
    cm.SaveDiv();
    if (dt == tmv::CH) cm.DivideUsing(tmv::LU);
    else cm.DivideUsing(dt);
    cm.SetDiv();
    Assert(cm.CheckDecomp(checkout),"CheckDecomp cm"); 
    csi.DivideUsing(dt);
    csi.SetDiv();
    if (csi.isherm()) {
      tmv::HermMatrix<std::complex<T> > csix = csi;
      csix.DivideUsing(dt);
      csix.SetDiv();
      Assert(csix.CheckDecomp(checkout),"CheckDecomp csix(herm)"); 
    } else {
      tmv::SymMatrix<std::complex<T> > csix = csi;
      csix.DivideUsing(dt);
      csix.SetDiv();
      Assert(csix.CheckDecomp(checkout),"CheckDecomp csix(sym)"); 
    }
    Assert(csi.CheckDecomp(checkout),"CheckDecomp csi"); 

    T ceps = EPS;
    if (pdc == Sing) ceps *= 1000;
    else ceps *= Norm(cm)*Norm(cm.Inverse());

    if (pdc != Sing) {
      tmv::Matrix<std::complex<T> > cm2 = cm;
      //std::cout<<"cm2 = "<<cm2<<std::endl;
      if (showacc) {
	std::cout<<"csi.Det = "<<csi.Det()<<", cm.Det = "<<cm2.Det()<<std::endl;
	std::cout<<"abs(csidet-cmdet) = "<<std::abs(csi.Det()-cm2.Det());
	std::cout<<"  csidet/cmdet = "<<csi.Det()/cm2.Det();
	std::cout<<"  EPS*abs(cmdet) = "<<ceps*std::abs(cm2.Det())<<std::endl;
	std::cout<<"abs(abs(csdet)-abs(cmdet)) = "<<
	  std::abs(std::abs(csi.Det())-std::abs(cm2.Det()));
	std::cout<<"  EPS*abs(cmdet) = "<<ceps*std::abs(cm2.Det())<<std::endl;
	std::cout<<"cm.LogDet() = "<<cm2.LogDet()<<std::endl;
	std::cout<<"csi.LogDet() = "<<csi.LogDet()<<std::endl;
	std::cout<<"abs(diff) = "<<std::abs(cm2.LogDet()-csi.LogDet())<<"   ";
	std::cout<<"eps = "<<eps<<std::endl;
	std::complex<T> cmsign,cssign;
	cm2.LogDet(&cmsign);
	csi.LogDet(&cssign);
	std::cout<<"cm sign = "<<cmsign<<std::endl;
	std::cout<<"csi sign = "<<cssign<<std::endl;
	std::cout<<"abs(diff) = "<<std::abs(cmsign-cssign)<<"   ";
	std::cout<<"eps = "<<eps<<std::endl;
      }
      Assert(std::abs(csi.Det()-cm2.Det()) < ceps*std::abs(cm2.Det()),"Sym CDet");
      std::complex<T> cmsign,cssign;
      Assert(std::abs(cm2.LogDet(&cmsign)-csi.LogDet(&cssign)) < 10*N*eps,
	  "Sym CLogDet");
      Assert(std::abs(cmsign-cssign) < 10*N*eps,"Sym CLogDet - sign");
    }

    tmv::Vector<std::complex<T> > cv = v1 * std::complex<T>(1,1);
    cv(1) += std::complex<T>(-1,5);
    cv(2) -= std::complex<T>(-1,5);

    // test real / complex
    tmv::Vector<std::complex<T> > y1 = v1/csi;
    tmv::Vector<std::complex<T> > y2 = v1/cm;
    if (showacc) 
      std::cout<<"v/cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(y1)<<std::endl;
    Assert(Norm(y1-y2) < ceps*Norm(y1),"Sym v/cs");

    // test complex / real
    y1 = cv/si;
    y2 = cv/m;
    if (showacc) 
      std::cout<<"cv/b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(y1)<<std::endl;
    Assert(Norm(y1-y2) < eps*Norm(y1),"Sym cv/b");

    // test complex / complex
    y1 = cv/csi;
    y2 = cv/cm;
    if (showacc) 
      std::cout<<"cv/cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(y1)<<std::endl;
    Assert(Norm(y1-y2) < ceps*Norm(y1),"Sym cv/cs");

    y1 = v1%csi;
    y2 = v1%cm;
    if (showacc) 
      std::cout<<"v%cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(y1)<<std::endl;
    Assert(Norm(y1-y2) < ceps*Norm(y1),"Sym v%cs");

    y1 = cv%si;
    y2 = cv%m;
    if (showacc) 
      std::cout<<"cv%b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(y1)<<std::endl;
    Assert(Norm(y1-y2) < eps*Norm(y1),"Sym cv%b");
    y1 = cv%csi;
    y2 = cv%cm;
    if (showacc) 
      std::cout<<"cv%cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(y1)<<std::endl;
    Assert(Norm(y1-y2) < ceps*Norm(y1),"Sym cv%cs");
  }

  TestSymDiv_A<T>(dt,pdc);
  TestSymDiv_B<T>(dt,pdc);
  TestSymDiv_C<T>(dt,pdc);
  TestSymDiv_D<T>(dt,pdc);
  TestSymDiv_E<T>(dt,pdc);

  std::cout<<PDLabel(pdc)<<" SymMatrix<"<<tmv::Type(T())<<"> Division using ";
  std::cout<<tmv::Text(dt)<<" passed all tests\n";
}

template <class T, tmv::UpLoType uplo, tmv::StorageType stor> 
void TestHermDecomp()
{
  for (int mattype = 0; mattype < 3; mattype++) {
    //std::cout<<"mattype = "<<mattype<<std::endl;
    // mattype = 0  is PosDef
    // mattype = 1  is Indef
    // mattype = 2  is Singular

    const int N = 8;

    tmv::HermMatrix<T,uplo,stor> m(N);
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
      if ((uplo == tmv::Upper && j>=i) || (uplo == tmv::Lower && i>=j)) 
	m(i,j) = T(2.+4*i-5*j);
    if (mattype == 0) {
      m /= T(100);
      m.diag().AddToAll(T(1));
      for(int i=0;i<N;i++) m.diag()(i) += T(i);
    } else if (mattype == 1) {
      m /= T(8);
      m.diag(0,2,4).Zero();
      m(4,0) = T(50);
      m.diag(1,4,7).AddToAll(T(10));
      m.row(7,0,5).AddToAll(T(10));
    } else {
      m.row(2,0,2).Zero();
      m.row(2,2,8).Zero();
      m.row(4,0,4) = m.row(5,0,4);
      m.row(4,5,8) = m.row(5,5,8);
      m(4,5) = m(4,4) = m(5,5);
    }

    tmv::HermMatrix<std::complex<T>,uplo,stor> c(N);
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
      if ((uplo == tmv::Upper && j>=i) || (uplo == tmv::Lower && i>=j)) 
	c(i,j) = std::complex<T>(2.+4*i-5*j,3.-i);
    c.diag().Imag().Zero();
    if (mattype == 0) {
      c /= T(100);
      c.diag().AddToAll(T(1));
      for(int i=0;i<N;i++) c.diag()(i) += T(i);
    } else if (mattype == 1) {
      c /= T(8);
      c.diag(0,2,4).Zero();
      c(4,0) = T(50);
      c.diag(1,4,7).AddToAll(T(10));
      c.row(7,0,5).AddToAll(T(10));
    } else {
      c.row(2,0,2).Zero();
      c.row(2,2,8).Zero();
      c.row(4,0,4) = c.row(5,0,4);
      c.row(4,5,8) = c.row(5,5,8);
      c(4,5) = c(4,4) = c(5,5);
    }

    T eps = EPS;
    T ceps = EPS;
    if (mattype != 2) {
      eps *= Norm(m) * Norm(m.Inverse());
      ceps *= Norm(c) * Norm(c.Inverse());
    } else {
      eps *= T(100); ceps *= T(100); 
    }
    //std::cout<<"m = "<<m<<std::endl;
    //std::cout<<"c = "<<c<<std::endl;
    //std::cout<<"eps = "<<eps<<"  "<<ceps<<std::endl;

    // CH Decomposition
    try {
      m.DivideUsing(tmv::CH);
      m.SetDiv();
      tmv::LowerTriMatrix<T> L = m.CHD().GetL();
      tmv::Matrix<T> LLt = L*L.Adjoint();
      Assert(Norm(m-LLt) < eps*Norm(m),"Herm CH");

      tmv::HermMatrix<T,uplo,stor> m2 = m;
      CH_Decompose(m2.View());
      L = m2.LowerTri();
      LLt = L*L.Adjoint();
      Assert(Norm(m-LLt) < eps*Norm(m),"Herm CH2");

      c.DivideUsing(tmv::CH);
      c.SetDiv();
      tmv::LowerTriMatrix<std::complex<T> > cL = c.CHD().GetL();
      tmv::Matrix<std::complex<T> > cLLt = cL*cL.Adjoint();
      Assert(Norm(c-cLLt) < ceps*Norm(c),"Herm C CH");

      tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
      CH_Decompose(c2.View());
      cL = c2.LowerTri();
      cLLt = cL*cL.Adjoint();
      Assert(Norm(c-cLLt) < ceps*Norm(c),"Herm C CH2");

      c2.Conjugate() = c;
      CH_Decompose(c2.Conjugate());
      cL = c2.Conjugate().LowerTri();
      cLLt = cL*cL.Adjoint();
      Assert(Norm(c-cLLt) < ceps*Norm(c),"Herm C CH2");
      Assert(mattype == 0,"Didn't throw NonPosDef, mattype != pos def");
    }
    catch(tmv::NonPosDef) 
    { Assert(mattype != 0,"Caught NonPosDef, mattype == pos def"); }

    // LDL Decomposition
    {
      m.DivideUsing(tmv::LU);
      m.SetDiv();
      tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.LUD().GetL();
      tmv::BandMatrix<T> D = m.LUD().GetD();
      const int* p = m.LUD().GetP();
      tmv::Matrix<T> LDL = L*D*L.Adjoint();
      LDL.ReversePermuteRows(p);
      LDL.ReversePermuteCols(p);
      Assert(Norm(m-LDL) < eps*Norm(m),"Herm LDL");

      tmv::HermMatrix<T,uplo,stor> m2 = m;
      tmv::HermBandMatrix<T,uplo,stor> D2(N,1);
      int p2[N];
      LDL_Decompose(m2.View(),D2.View(),p2);
      L = m2.LowerTri(tmv::UnitDiag);
      LDL = L*D2*L.Adjoint();
      LDL.ReversePermuteRows(p2);
      LDL.ReversePermuteCols(p2);
      Assert(Norm(m-LDL) < eps*Norm(m),"Herm LDL2");

      c.DivideUsing(tmv::LU);
      c.SetDiv();
      tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = c.LUD().GetL();
      tmv::BandMatrix<std::complex<T> > cD = c.LUD().GetD();
      p = c.LUD().GetP();
      tmv::Matrix<std::complex<T> > cLDL = cL*cD*cL.Adjoint();
      cLDL.ReversePermuteRows(p);
      cLDL.ReversePermuteCols(p);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL");

      tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
      tmv::HermBandMatrix<std::complex<T>,uplo,stor> cD2(N,1);
      LDL_Decompose(c2.View(),cD2.View(),p2);
      cL = c2.LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2*cL.Adjoint();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL2");

      c2.Conjugate() = c;
      LDL_Decompose(c2.Conjugate(),cD2.View(),p2);
      cL = c2.Conjugate().LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2*cL.Adjoint();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL3");

      c2 = c;
      LDL_Decompose(c2.View(),cD2.Conjugate(),p2);
      cL = c2.LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2.Conjugate()*cL.Adjoint();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL4");

      c2.Conjugate() = c;
      LDL_Decompose(c2.Conjugate(),cD2.Conjugate(),p2);
      cL = c2.Conjugate().LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2.Conjugate()*cL.Adjoint();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL5");
    }

    // SV Decomposition
    {
      m.DivideUsing(tmv::SV);
      m.SaveDiv();
      m.SetDiv();
      tmv::Matrix<T> U = m.SVD().GetU();
      tmv::DiagMatrix<T> S = m.SVD().GetS();
      tmv::Matrix<T> V = m.SVD().GetV();
      Assert(Norm(m-U*S*V) < eps*Norm(m),"Herm SV");

      tmv::Matrix<T> U2(N,N);
      tmv::DiagMatrix<T> S2(N);
      tmv::Matrix<T> V2(N,N);
      SV_Decompose(m,U2.View(),S2.View(),V2.View());
      Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"Herm SV2");

      tmv::HermMatrix<T,uplo,stor> m2 = m;
      SV_Decompose(m2.View(),S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"Herm SV3");
      SV_Decompose(m,U2.View(),S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"Herm SV4 S");
      Assert(Norm(U2*S2*S2*U2.Adjoint()-m*m.Adjoint()) < 
	  eps*Norm(m)*Norm(m),"HermBand C SV4 U");
      SV_Decompose(m,S2.View(),V2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"Herm SV5 S");
      Assert(Norm(m.Adjoint()*m-V2.Adjoint()*S2*S2*V2) < 
	  eps*Norm(m)*Norm(m),"Herm SV5 V");

      c.DivideUsing(tmv::SV);
      c.SaveDiv();
      c.SetDiv();
      tmv::Matrix<std::complex<T> > cU = c.SVD().GetU();
      S = c.SVD().GetS();
      tmv::Matrix<std::complex<T> > cV = c.SVD().GetV();
      Assert(Norm(c-cU*S*cV) < ceps*Norm(c),"Herm C SV");

      tmv::Matrix<std::complex<T> > cU2(N,N);
      tmv::Matrix<std::complex<T> > cV2(N,N);
      SV_Decompose(c,cU2.View(),S2.View(),cV2.View());
      Assert(Norm(c-cU2*S2*cV2) < ceps*Norm(c),"Herm C SV2");

      tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
      SV_Decompose(c2.View(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"Herm C SV3");
      SV_Decompose(c,cU2.View(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"Herm C SV4 S");
      Assert(Norm(cU2*S2*S2*cU2.Adjoint()-c*c.Adjoint()) < 
	  ceps*Norm(c)*Norm(c),"Herm C SV4 U");
      SV_Decompose(c,S2.View(),cV2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"Herm C SV5 S");
      Assert(Norm(c.Adjoint()*c-cV2.Adjoint()*S2*S2*cV2) < 
	  ceps*Norm(c)*Norm(c),"Herm C SV5 V");

      c2 = c.Conjugate();
      SV_Decompose(c2.Conjugate(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"Herm C SV6");
      SV_Decompose(c,cU2.Conjugate(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"Herm C SV4 S");
      Assert(Norm(cU2.Conjugate()*S2*S2*cU2.Transpose()-c*c.Adjoint()) < 
	  ceps*Norm(c)*Norm(c),"Herm C SV4 U");
    }
    
    // Eigen
    {
      tmv::Matrix<T> V(N,N);
      tmv::Vector<T> L(N);
      Eigen(m,V.View(),L.View());
      Assert(Norm(m*V-V*DiagMatrixViewOf(L)) < eps*Norm(m),"Herm Eigen");

      tmv::Vector<T> L2(N);
      tmv::HermMatrix<T,uplo,stor> m2 = m;
      Eigen(m2.View(),L2.View());
      Assert(Norm(L2-L) < eps*Norm(L),"Herm Eigen2");

      tmv::Matrix<std::complex<T> > cV(N,N);
      Eigen(c,cV.View(),L.View());
      Assert(Norm(c*cV-cV*DiagMatrixViewOf(L)) < eps*Norm(c),"Herm C Eigen");

      tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
      Eigen(c2.View(),L2.View());
      Assert(Norm(L2-L) < eps*Norm(L),"Herm C Eigen2");

      Eigen(c,cV.Conjugate(),L.View());
      Assert(Norm(c*cV.Conjugate()-cV.Conjugate()*DiagMatrixViewOf(L)) < 
	  eps*Norm(c),"Herm C Eigen3");

      Eigen(c.Conjugate(),cV.View(),L.View());
      Assert(Norm(c.Conjugate()*cV-cV*DiagMatrixViewOf(L)) < eps*Norm(c),
	  "Herm C Eigen4");

      Eigen(c.Conjugate(),cV.Conjugate(),L.View());
      Assert(Norm(c*cV-cV*DiagMatrixViewOf(L)) < eps*Norm(c),"Herm C Eigen5");

      c2.Conjugate() = c;
      Eigen(c2.Conjugate(),L2.View());
      Assert(Norm(L2-L) < eps*Norm(L),"Herm C Eigen6");
    }

    // Square Root
    try {
      tmv::HermMatrix<T,uplo,stor> S = m;
      SquareRoot(S.View());
      Assert(Norm(m-S*S) < eps*Norm(m),"Herm Square Root");

      tmv::HermMatrix<std::complex<T>,uplo,stor> cS = c;
      SquareRoot(cS.View());

      cS.Conjugate() = c;
      SquareRoot(cS.Conjugate());
      Assert(Norm(c-cS.Conjugate()*cS.Conjugate()) < eps*Norm(c),
	  "Herm C Square Root 2");
      Assert(mattype == 0,"Didn't throw NonPosDef, mattype != pos def");
    }
    catch(tmv::NonPosDef) 
    { Assert(mattype != 0,"Caught NonPosDef, mattype == pos def"); }
  }
}

template <class T, tmv::UpLoType uplo, tmv::StorageType stor> 
void TestSymDecomp()
{
  for (int mattype = 0; mattype < 2; mattype++) {
    // mattype = 0  is Normal
    // mattype = 1  is Singular

    const int N = 8;

    tmv::SymMatrix<T,uplo,stor> m(N);
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
      if ((uplo == tmv::Upper && j>=i) || (uplo == tmv::Lower && i>=j)) 
	m(i,j) = T(2.+4*i-5*j);
    if (mattype == 0) {
      m /= T(100);
      m.diag().AddToAll(T(1));
      for(int i=0;i<N;i++) m.diag()(i) += T(i);
    } else if (mattype == 1) {
      m /= T(8);
      m.diag(0,2,4).Zero();
      m(4,0) = T(50);
      m.diag(1,4,7).AddToAll(T(10));
      m.row(7,0,5).AddToAll(T(10));
    } else {
      m.row(2,0,2).Zero();
      m.row(2,2,8).Zero();
      m.row(4,0,4) = m.row(5,0,4);
      m.row(4,5,8) = m.row(5,5,8);
      m(4,5) = m(4,4) = m(5,5);
    }

    tmv::SymMatrix<std::complex<T>,uplo,stor> c(N);
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
      c(i,j) = std::complex<T>(2.+4*i-5*j,3.-i);
    if (mattype == 0) {
      c /= T(100);
      c.diag().AddToAll(T(1));
      for(int i=0;i<N;i++) c.diag()(i) += T(i);
    } else if (mattype == 1) {
      c /= T(8);
      c.diag(0,2,4).Zero();
      c(4,0) = T(50);
      c.diag(1,4,7).AddToAll(T(10));
      c.row(7,0,5).AddToAll(T(10));
    } else {
      c.row(2,0,2).Zero();
      c.row(2,2,8).Zero();
      c.row(4,0,4) = c.row(5,0,4);
      c.row(4,5,8) = c.row(5,5,8);
      c(4,5) = c(4,4) = c(5,5);
    }

    T eps = EPS;
    T ceps = EPS;
    if (mattype != 1) {
      eps *= Norm(m) * Norm(m.Inverse());
      ceps *= Norm(c) * Norm(c.Inverse());
    } else {
      eps *= T(100);
      ceps *= T(100);
    }

    // LDL Decomposition
    {
      m.DivideUsing(tmv::LU);
      m.SetDiv();
      tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.LUD().GetL();
      tmv::BandMatrix<T> D = m.LUD().GetD();
      const int* p = m.LUD().GetP();
      tmv::Matrix<T> LDL = L*D*L.Transpose();
      LDL.ReversePermuteRows(p);
      LDL.ReversePermuteCols(p);
      Assert(Norm(m-LDL) < eps*Norm(m),"Sym LDL");

      tmv::SymMatrix<T,uplo,stor> m2 = m;
      tmv::SymBandMatrix<T,uplo,stor> D2(N,1);
      int p2[N];
      LDL_Decompose(m2.View(),D2.View(),p2);
      L = m2.LowerTri(tmv::UnitDiag);
      LDL = L*D2*L.Transpose();
      LDL.ReversePermuteRows(p2);
      LDL.ReversePermuteCols(p2);
      Assert(Norm(m-LDL) < eps*Norm(m),"Sym LDL2");

      c.DivideUsing(tmv::LU);
      c.SetDiv();
      tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = c.LUD().GetL();
      tmv::BandMatrix<std::complex<T> > cD = c.LUD().GetD();
      tmv::Matrix<std::complex<T> > cLDL = cL*cD*cL.Transpose();
      cLDL.ReversePermuteRows(p);
      cLDL.ReversePermuteCols(p);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL");

      tmv::SymMatrix<std::complex<T>,uplo,stor> c2 = c;
      tmv::SymBandMatrix<std::complex<T>,uplo,stor> cD2(N,1);
      LDL_Decompose(c2.View(),cD2.View(),p2);
      cL = c2.LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2*cL.Transpose();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL2");

      c2.Conjugate() = c;
      LDL_Decompose(c2.Conjugate(),cD2.View(),p2);
      cL = c2.Conjugate().LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2*cL.Transpose();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL3");

      c2= c;
      LDL_Decompose(c2.View(),cD2.Conjugate(),p2);
      cL = c2.LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2.Conjugate()*cL.Transpose();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL4");

      c2.Conjugate() = c;
      LDL_Decompose(c2.Conjugate(),cD2.Conjugate(),p2);
      cL = c2.Conjugate().LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2.Conjugate()*cL.Transpose();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL5");
    }

    // SV Decomposition
    {
      m.DivideUsing(tmv::SV);
      m.SaveDiv();
      m.SetDiv();
      tmv::Matrix<T> U = m.SVD().GetU();
      tmv::DiagMatrix<T> S = m.SVD().GetS();
      tmv::Matrix<T> V = m.SVD().GetV();
      Assert(Norm(m-U*S*V) < eps*Norm(m),"Sym SV");

      tmv::Matrix<T> U2(N,N);
      tmv::DiagMatrix<T> S2(N);
      tmv::Matrix<T> V2(N,N);
      SV_Decompose(m,U2.View(),S2.View(),V2.View());
      Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"Sym SV2");

      tmv::SymMatrix<T,uplo,stor> m2 = m;
      SV_Decompose(m2.View(),S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"Sym SV3");
      SV_Decompose(m,U2.View(),S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"Sym SV4 S");
      Assert(Norm(U2-U) < eps*Norm(U),"Sym SV4 U");
      SV_Decompose(m,S2.View(),V2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"Sym SV5 S");
      Assert(Norm(m.Adjoint()*m-V2.Adjoint()*S2*S2*V2) < eps*Norm(m)*Norm(m),
	  "Sym SV5 V");

      c.DivideUsing(tmv::SV);
      c.SaveDiv();
      c.SetDiv();
      tmv::Matrix<std::complex<T> > cU = c.SymSVD().GetU();
      S = c.SymSVD().GetS();
      tmv::Matrix<std::complex<T> > cV = c.SymSVD().GetV();
      Assert(Norm(c-cU*S*cV) < ceps*Norm(c),"Sym C SV");

      tmv::Matrix<std::complex<T> > cU2(N,N);
      tmv::Matrix<std::complex<T> > cV2(N,N);
      SV_Decompose(c,cU2.View(),S2.View(),cV2.View());
      Assert(Norm(c-cU2*S2*cV2) < ceps*Norm(c),"Sym C SV2");

      tmv::SymMatrix<std::complex<T>,uplo,stor> c2 = c;
      SV_Decompose(c2.View(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"Sym C SV3");
      SV_Decompose(c,cU2.View(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"Sym C SV4 S");
      Assert(Norm(cU2-cU) < ceps*Norm(cU),"Sym C SV4 U");
      SV_Decompose(c,S2.View(),cV2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"Sym C SV5 S");
      Assert(Norm(c.Adjoint()*c-cV2.Adjoint()*S2*S2*cV2) < ceps*Norm(c)*Norm(c),
	  "Sym C SV5 V");

      c2 = c.Conjugate();
      SV_Decompose(c2.Conjugate(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"Sym C SV6");
      SV_Decompose(c,cU2.Conjugate(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"Sym C SV7 S");
      Assert(Norm(cU2.Conjugate()-cU) < ceps*Norm(cU),"Sym C SV7 U");
    }
  }
}

template <class T, tmv::StorageType stor> 
void TestPolar()
{
  for (int mattype = 0; mattype < 4; mattype++) {
    // mattype = 0  is Square
    // mattype = 1  is NonSquare slightly tall
    // mattype = 2  is NonSquare very tall
    // mattype = 3  is Singular

    int M = 8;
    int N = 8;
    if (mattype == 1) M = 11;
    else if (mattype == 2) M = 45;

    tmv::Matrix<T,stor> m(M,N);
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) m(i,j) = T(2.+4*i-5*j);
    m(0,0) = T(14);
    m(1,0) = T(-2);
    m(2,0) = T(7);
    m(3,4) = T(-10);
    if (mattype != 3) m.diag() *= T(30);

    tmv::Matrix<std::complex<T>,stor> c(M,N);
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
      c(i,j) = std::complex<T>(2.+4*i-5*j,3.-i);
    c(0,0) *= T(14);
    c(1,0) *= T(-2);
    c(2,0) *= T(7);
    c(7,6) *= T(-10);
    if (mattype != 3) c.diag() *= T(30);

    T eps = EPS;
    T ceps = EPS;
    if (mattype != 3) {
      eps *= Norm(m) * Norm(m.Inverse());
      ceps *= Norm(c) * Norm(c.Inverse());
    } else {
      eps *= T(100);
      ceps *= T(100);
    }

    // Matrix Polar Decomposition
    {
      tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor> P(N);
      tmv::Matrix<T,stor> U = m;
      Polar_Decompose(U.View(),P.View());
      Assert(Norm(m-U*P) < eps*Norm(m),"Polar");
      Assert(Norm(U.Adjoint()*U-T(1)) < eps,"Polar UtU");

      U = m;
      Polar_Decompose(U.View(),P.Transpose());
      Assert(Norm(m-U*P) < eps*Norm(m),"Polar2");
      Assert(Norm(U.Adjoint()*U-T(1)) < eps,"Polar2 UtU");

      tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> cP(N);
      tmv::Matrix<std::complex<T>,stor> cU = c;
      Polar_Decompose(cU.View(),cP.View());
      Assert(Norm(c-cU*cP) < ceps*Norm(c),"C Polar");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar UtU");

      cU = c;
      Polar_Decompose(cU.View(),cP.Adjoint());
      Assert(Norm(c-cU*cP) < ceps*Norm(c),"C Polar2");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar2 UtU");

      cU = c;
      Polar_Decompose(cU.View(),cP.Transpose());
      Assert(Norm(c-cU*cP.Transpose()) < ceps*Norm(c),"C Polar3");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar3 UtU");

      cU = c;
      Polar_Decompose(cU.View(),cP.Conjugate());
      Assert(Norm(c-cU*cP.Transpose()) < ceps*Norm(c),"C Polar4");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar4 UtU");

      cU.Conjugate() = c;
      Polar_Decompose(cU.Conjugate(),cP.View());
      Assert(Norm(c-cU.Conjugate()*cP) < ceps*Norm(c),"C Polar5");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar5 UtU");

      cU.Conjugate() = c;
      Polar_Decompose(cU.Conjugate(),cP.Adjoint());
      Assert(Norm(c-cU.Conjugate()*cP) < ceps*Norm(c),"C Polar6");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar6 UtU");

      cU.Conjugate() = c;
      Polar_Decompose(cU.Conjugate(),cP.Transpose());
      Assert(Norm(c-cU.Conjugate()*cP.Transpose()) < ceps*Norm(c),"C Polar7");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar7 UtU");

      cU.Conjugate() = c;
      Polar_Decompose(cU.Conjugate(),cP.Conjugate());
      Assert(Norm(c-cU.Conjugate()*cP.Transpose()) < ceps*Norm(c),"C Polar8");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar8 UtU");
    }

    // BandMatrix Polar Decomposition
    {
      tmv::BandMatrixView<T> b(m.View(),1,3);
      tmv::BandMatrixView<std::complex<T> > cb(c.View(),1,3);

      tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor> P(N);
      tmv::Matrix<T,stor> U(b.colsize(),N);
      Polar_Decompose(b,U.View(),P.View());
      Assert(Norm(b-U*P) < eps*Norm(b),"Band Polar");
      Assert(Norm(U.Adjoint()*U-T(1)) < eps,"Band Polar UtU");

      Polar_Decompose(b,U.View(),P.Transpose());
      Assert(Norm(b-U*P) < eps*Norm(b),"Band Polar2");
      Assert(Norm(U.Adjoint()*U-T(1)) < eps,"Band Polar2 UtU");

      tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> cP(N);
      tmv::Matrix<std::complex<T>,stor> cU(cb.colsize(),N);
      Polar_Decompose(cb,cU.View(),cP.View());
      Assert(Norm(cb-cU*cP) < ceps*Norm(cb),"C Band Polar");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar UtU");

      Polar_Decompose(cb,cU.View(),cP.Adjoint());
      Assert(Norm(cb-cU*cP) < ceps*Norm(cb),"C Band Polar2");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar2 UtU");

      Polar_Decompose(cb,cU.View(),cP.Transpose());
      Assert(Norm(cb-cU*cP.Transpose()) < ceps*Norm(cb),"C Band Polar3");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar3 UtU");

      Polar_Decompose(cb,cU.View(),cP.Conjugate());
      Assert(Norm(cb-cU*cP.Transpose()) < ceps*Norm(cb),"C Band Polar4");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar4 UtU");

      Polar_Decompose(cb,cU.Conjugate(),cP.View());
      Assert(Norm(cb-cU.Conjugate()*cP) < ceps*Norm(cb),"C Band Polar5");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar5 UtU");

      Polar_Decompose(cb,cU.Conjugate(),cP.Adjoint());
      Assert(Norm(cb-cU.Conjugate()*cP) < ceps*Norm(cb),"C Band Polar6");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar6 UtU");

      Polar_Decompose(cb,cU.Conjugate(),cP.Transpose());
      Assert(Norm(cb-cU.Conjugate()*cP.Transpose()) < ceps*Norm(cb),"C Band Polar7");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar7 UtU");

      Polar_Decompose(cb,cU.Conjugate(),cP.Conjugate());
      Assert(Norm(cb-cU.Conjugate()*cP.Transpose()) < ceps*Norm(cb),"C Band Polar8");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar8 UtU");
    }
  }
}

template <class T> void TestAllSymDiv()
{
  TestHermDecomp<T,tmv::Upper,tmv::ColMajor>();
  TestHermDecomp<T,tmv::Upper,tmv::RowMajor>();
  TestHermDecomp<T,tmv::Lower,tmv::ColMajor>();
  TestHermDecomp<T,tmv::Lower,tmv::RowMajor>();
  TestSymDecomp<T,tmv::Upper,tmv::ColMajor>();
  TestSymDecomp<T,tmv::Upper,tmv::RowMajor>();
  TestSymDecomp<T,tmv::Lower,tmv::ColMajor>();
  TestSymDecomp<T,tmv::Lower,tmv::RowMajor>();
  TestPolar<T,tmv::RowMajor>();
  TestPolar<T,tmv::ColMajor>();
  TestSymDiv<T>(tmv::CH,PosDef);
  TestSymDiv<T>(tmv::LU,PosDef);
  TestSymDiv<T>(tmv::LU,InDef);
  TestSymDiv<T>(tmv::SV,PosDef);
  TestSymDiv<T>(tmv::SV,InDef);
  TestSymDiv<T>(tmv::SV,Sing);
}

#ifdef INST_DOUBLE
template void TestAllSymDiv<double>();
#endif
#ifdef INST_FLOAT
template void TestAllSymDiv<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllSymDiv<long double>();
#endif
#ifdef INST_INT
template void TestAllSymDiv<int>();
#endif



