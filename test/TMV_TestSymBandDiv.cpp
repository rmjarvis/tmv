#ifdef NDEBUG
#undef NDEBUG
#endif

#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Band.h"
#include "TMV_Diag.h"
#include "TMV_Sym.h"
#include "TMV_TestSymBandArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> static bool IsPosDef(const tmv::GenSymBandMatrix<T>& m)
{
  try {
    tmv::HermBandMatrix<T> m2 = m;
    CH_Decompose(m2.View());
  }
  catch (tmv::NonPosDef) {
#ifdef XTEST
    if (showacc)
      std::cout<<"caught\n";
    for(size_t i=1;i<=m.size();i++) {
      T d = tmv::Matrix<T>(m.SubSymBandMatrix(0,i)).Det();
      if (showacc) 
	std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
      if (!(tmv::REAL(d) > 0)) {
	return false;
      }
    }
    std::cout<<"m = "<<Type(m)<<"  "<<m<<std::endl;
    for(size_t i=1;i<=m.size();i++) {
      T d = tmv::Matrix<T>(m.SubSymBandMatrix(0,i)).Det();
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
    T d = tmv::Matrix<T>(m.SubSymBandMatrix(0,i)).Det();
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

template <class T> void TestSymBandDiv(tmv::DivType dt, PosDefCode pdc)
{
  const int N = 10;

  std::vector<tmv::SymBandMatrixView<T> > sb;
  std::vector<tmv::SymBandMatrixView<std::complex<T> > > csb;
  MakeSymBandList(sb,csb,pdc);

  size_t ntot = sb.size();

  tmv::Vector<T> v1(N);
  tmv::Vector<T> v2(N-1);
  for (int i=0; i<N; ++i) v1(i) = 16.-3*i; 
  for (int i=0; i<N-1; ++i) v2(i) = -7.+2*i; 

  std::ostream* checkout = showdiv ? &std::cout : 0;

  for(size_t i=START;i<ntot;i++) {
    if (i>START) break;
    tmv::SymBandMatrixView<T> si = sb[i];
    tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];
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
      tmv::SymBandMatrix<T> six = si;
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
    Assert(Norm(x1-x2) < eps*Norm(x1),"SymBand v/b");

    x1 = v1%si;
    x2 = v1%m;
    if (showacc)
      std::cout<<"v%b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(x1)<<std::endl;
    Assert(Norm(x1-x2) < eps*Norm(x1),"SymBand v%b");

    tmv::Matrix<T,tmv::ColMajor> sinv = si.Inverse();
    tmv::Matrix<T,tmv::ColMajor> minv = m.Inverse();
    if (showacc) {
      std::cout<<"sinv = "<<sinv<<std::endl;
      std::cout<<"minv = "<<minv<<std::endl;
      std::cout<<"Norm(minv-sinv) = "<<Norm(minv-sinv)<<"  "<<eps*Norm(sinv)<<std::endl;
    }
    Assert(Norm(sinv-minv) < eps*Norm(sinv),"SymBand Inverse");

    if (showacc) {
      std::cout<<"si.Det = "<<si.Det()<<", m.Det = "<<m.Det()<<std::endl;
      std::cout<<"abs(sdet-mdet) = "<<std::abs(si.Det()-m.Det());
      std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(m.Det())<<std::endl;
      std::cout<<"abs(abs(sdet)-abs(mdet)) = "<<std::abs(std::abs(si.Det())-std::abs(m.Det()));
      std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(m.Det())<<std::endl;
    }
    if (pdc != Sing) {
      Assert(std::abs(m.Det()-si.Det()) < eps*std::abs(m.Det()),"SymBand Det");
      T msign,ssign;
      Assert(std::abs(m.LogDet(&msign)-si.LogDet(&ssign)) < 10*N*eps,
	  "SymBand LogDet");
      Assert(std::abs(msign-ssign) < 10*N*eps, "SymBand LogDet - sign");
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
      tmv::SymBandMatrix<std::complex<T> > csix = csi;
      csix.DivideUsing(dt);
      csix.SetDiv();
      Assert(csix.CheckDecomp(checkout),"CheckDecomp csix(sym)"); 
    }
    Assert(csi.CheckDecomp(checkout),"CheckDecomp csi"); 

    T ceps = EPS;
    if (pdc == Sing) ceps *= 1000;
    else ceps *= Norm(cm)*Norm(cm.Inverse());

    if (showacc) {
      std::cout<<"csi.Det = "<<csi.Det()<<", cm.Det = "<<cm.Det()<<std::endl;
      std::cout<<"abs(csidet-cmdet) = "<<std::abs(csi.Det()-cm.Det());
      std::cout<<"  csidet/cmdet = "<<csi.Det()/cm.Det();
      std::cout<<"  EPS*abs(cmdet) = "<<ceps*std::abs(cm.Det())<<std::endl;
      std::cout<<"abs(abs(csdet)-abs(cmdet)) = "<<std::abs(std::abs(csi.Det())-std::abs(cm.Det()));
      std::cout<<"  EPS*abs(cmdet) = "<<ceps*std::abs(cm.Det())<<std::endl;
    }
    if (pdc != Sing) {
      Assert(std::abs(csi.Det()-cm.Det()) < ceps*std::abs(cm.Det()),"SymBand CDet");
      std::complex<T> cmsign,cssign;
      Assert(std::abs(cm.LogDet(&cmsign)-csi.LogDet(&cssign)) < 10*N*eps,
	  "SymBand CLogDet");
      Assert(std::abs(cmsign-cssign) < 10*N*eps, "SymBand CLogDet - sign");
    }

    tmv::Vector<std::complex<T> > cv = v1 * std::complex<T>(1,1);
    cv(1) += std::complex<T>(-1,5);
    cv(2) -= std::complex<T>(-1,5);

    // test real / complex
    tmv::Vector<std::complex<T> > y1 = v1/csi;
    tmv::Vector<std::complex<T> > y2 = v1/cm;
    if (showacc) 
      std::cout<<"v/cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(y1)<<std::endl;
    Assert(Norm(y1-y2) < ceps*Norm(y1),"SymBand v/cs");

    // test complex / real
    y1 = cv/si;
    y2 = cv/m;
    if (showacc) 
      std::cout<<"cv/b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(y1)<<std::endl;
    Assert(Norm(y1-y2) < eps*Norm(y1),"SymBand cv/b");

    // test complex / complex
    y1 = cv/csi;
    y2 = cv/cm;
    if (showacc) 
      std::cout<<"cv/cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(y1)<<std::endl;
    Assert(Norm(y1-y2) < ceps*Norm(y1),"SymBand cv/cs");

    y1 = v1%csi;
    y2 = v1%cm;
    if (showacc) 
      std::cout<<"v%cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(y1)<<std::endl;
    Assert(Norm(y1-y2) < ceps*Norm(y1),"SymBand v%cs");

    y1 = cv%si;
    y2 = cv%m;
    if (showacc) 
      std::cout<<"cv%b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(y1)<<std::endl;
    Assert(Norm(y1-y2) < eps*Norm(y1),"SymBand cv%b");
    y1 = cv%csi;
    y2 = cv%cm;
    if (showacc) 
      std::cout<<"cv%cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(y1)<<std::endl;
    Assert(Norm(y1-y2) < ceps*Norm(y1),"SymBand cv%cs");
  }

  TestSymBandDiv_A<T>(dt,pdc);
  TestSymBandDiv_B<T>(dt,pdc);
  TestSymBandDiv_C<T>(dt,pdc);
  TestSymBandDiv_D<T>(dt,pdc);
  TestSymBandDiv_E<T>(dt,pdc);
  TestSymBandDiv_F<T>(dt,pdc);

  std::cout<<PDLabel(pdc)<<" SymBandMatrix<"<<tmv::Type(T())<<"> Division using ";
  std::cout<<tmv::Text(dt)<<" passed all tests\n";
}

template <class T, tmv::UpLoType uplo, tmv::StorageType stor> 
void TestHermBandDecomp()
{
  for (int mattype = 0; mattype < 9; mattype++) {
    //std::cout<<"mattype = "<<mattype<<std::endl;
    //std::cout<<"uplo, stor = "<<Text(uplo)<<"  "<<Text(stor)<<std::endl;
    // mattype = 0  is PosDef, nlo = 3
    // mattype = 1  is Indef, nlo = 3
    // mattype = 2  is Singular, nlo = 3
    // mattype = 3  is PosDef, nlo = 1
    // mattype = 4  is Indef, nlo = 1
    // mattype = 5  is Singular, nlo = 1
    // mattype = 6  is PosDef, nlo = 0
    // mattype = 7  is Indef, nlo = 0
    // mattype = 8  is Singular, nlo = 0

    const int N = 8;
    int nlo = 3;
    if (mattype / 3 == 1) nlo = 1;
    else if (mattype / 3 == 2) nlo = 0;

    tmv::HermBandMatrix<T,uplo,stor> m(N,nlo);
    //std::cout<<"m = "<<Type(m)<<std::endl;
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
      if ((uplo == tmv::Upper && j>=i && j<=i+nlo) || 
	  (uplo == tmv::Lower && i>=j && i<=j+nlo)) 
	m(i,j) = T(2.5+4*i-5*j);
    //std::cout<<"m = "<<m<<std::endl;
    if (mattype%3 == 0) {
      m /= T(100);
      m.diag().AddToAll(T(1));
      for(int i=0;i<N;i++) m.diag()(i) += T(i);
    } else if (mattype%3 == 1) {
      m /= T(8);
      if (nlo > 0) {
	m(2,1) = T(50);
	m.diag(1,4,7).AddToAll(T(10));
	m.row(7,7-nlo,7).AddToAll(T(10));
      }
    } else {
      if (nlo > 0) {
	if (nlo > 2) m.row(2,0,2).Zero();
	else m.row(2,2-nlo,2).Zero();
	m.row(2,2,2+nlo+1).Zero();
	m.row(4,5-nlo,4) = m.row(5,5-nlo,4);
	m.row(4,5,4+nlo+1) = m.row(5,5,4+nlo+1);
	if (4-nlo >= 0) m(4,4-nlo) = T(0);
	if (5+nlo < N) m(5,5+nlo) = T(0);
	m(4,5) = m(4,4) = m(5,5);
      } else {
	m.diag(0,2,4).Zero();
      }
    }
    //std::cout<<"m = "<<m<<std::endl;

    tmv::HermBandMatrix<std::complex<T>,uplo,stor> c(N,nlo);
    //std::cout<<"c = "<<Type(c)<<std::endl;
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
      if ((uplo == tmv::Upper && j>=i && j<=i+nlo) || 
	  (uplo == tmv::Lower && i>=j && i<=j+nlo)) 
	c(i,j) = std::complex<T>(2.5+4*i-5*j,3.-i);
    c.diag().Imag().Zero();
    if (mattype%3 == 0) {
      c /= T(100);
      c.diag().AddToAll(T(1));
      for(int i=0;i<N;i++) c.diag()(i) += T(i);
    } else if (mattype%3 == 1) {
      c /= T(8);
      if (nlo > 0) {
	c(2,1) = T(50);
	c.diag(1,4,7).AddToAll(T(10));
	c.row(7,7-nlo,7).AddToAll(T(10));
      }
    } else {
      if (nlo > 0) {
	if (nlo > 2) c.row(2,0,2).Zero();
	else c.row(2,2-nlo,2).Zero();
	c.row(2,3,2+nlo+1).Zero();
	c.row(4,5-nlo,4) = c.row(5,5-nlo,4);
	c.row(4,5,4+nlo+1) = c.row(5,5,4+nlo+1);
	if (4-nlo >= 0) c(4,4-nlo) = T(0);
	if (5+nlo < N) c(5,5+nlo) = T(0);
	c(4,5) = c(4,4) = c(5,5);
      } else {
	c.diag(0,2,4).Zero();
      }
    }
    //std::cout<<"c = "<<c<<std::endl;

    T eps = EPS;
    T ceps = EPS;
    if (mattype % 3 != 2) {
      eps *= Norm(m) * Norm(m.Inverse());
      ceps *= Norm(c) * Norm(c.Inverse());
    } else {
      eps *= T(100);
      ceps *= T(100);
    }
    //std::cout<<"eps = "<<eps<<"  "<<ceps<<std::endl;

    // CH Decomposition
    //std::cout<<"CHD"<<std::endl;
    if (mattype % 3 == 0) {
      tmv::BandMatrix<T> L(N,N,nlo,0);
      tmv::BandMatrix<T> LLt(N,N,nlo,nlo);
      if (mattype == 0) {
	m.DivideUsing(tmv::CH);
	m.SetDiv();
	L = m.CHD().GetL();
	LLt = L*L.Adjoint();
	Assert(Norm(m-LLt) < eps*Norm(m),"HermBand CH");
      }

      tmv::HermBandMatrix<T,uplo,stor> m2 = m;
      CH_Decompose(m2.View());
      L = m2.LowerBand();
      LLt = L*L.Adjoint();
      Assert(Norm(m-LLt) < eps*Norm(m),"HermBand CH2");

      tmv::BandMatrix<std::complex<T> > cL(N,N,nlo,0);
      tmv::BandMatrix<std::complex<T> > cLLt(N,N,nlo,nlo);
      if (mattype == 0) {
	c.DivideUsing(tmv::CH);
	c.SetDiv();
	cL = c.CHD().GetL();
	cLLt = cL*cL.Adjoint();
	Assert(Norm(c-cLLt) < ceps*Norm(c),"HermBand C CH");
      }

      tmv::HermBandMatrix<std::complex<T>,uplo,stor> c2 = c;
      CH_Decompose(c2.View());
      cL = c2.LowerBand();
      cLLt = cL*cL.Adjoint();
      Assert(Norm(c-cLLt) < ceps*Norm(c),"HermBand C CH2");

      c2.Conjugate() = c;
      CH_Decompose(c2.Conjugate());
      cL = c2.Conjugate().LowerBand();
      cLLt = cL*cL.Adjoint();
      Assert(Norm(c-cLLt) < ceps*Norm(c),"HermBand C CH2");
    }

    // LDL Decomposition
    //std::cout<<"LDL"<<std::endl;
    if (mattype /3 == 1) try {
      tmv::BandMatrix<T> L(N,N,1,0);
      tmv::DiagMatrix<T> D(N);
      tmv::BandMatrix<T> LDL(N,N,1,1);
      if (mattype == 3) {
	m.DivideUsing(tmv::CH);
	m.SetDiv();
	L = m.CHD().GetL();
	D = m.CHD().GetD();
	LDL = L*D*L.Adjoint();
	Assert(Norm(m-LDL) < eps*Norm(m),"HermBand LDL");
      }

      tmv::HermBandMatrix<T,uplo,stor> m2 = m;
      LDL_Decompose(m2.View());
      (L = m2.LowerBand()).diag().SetAllTo(T(1));
      D = DiagMatrixViewOf(m2.diag());
      LDL = L*D*L.Adjoint();
      Assert(Norm(m-LDL) < eps*Norm(m),"HermBand LDL2");

      // Alt calculation:
      // (I + L) D (I + Lt)
      // D + LD + DLt = LDLt
      tmv::DiagMatrix<T> E(m2.diag(-1));
      LDL.Zero();
      LDL.diag(0,1,N) = (E*D.SubDiagMatrix(0,N-1)*E).diag();
      LDL.diag(-1) = (E*D.SubDiagMatrix(0,N-1)).diag();
      LDL += D;
      LDL.diag(1) = LDL.diag(-1);
      Assert(Norm(m-LDL) < eps*Norm(m),"HermBand LDL2 Alt");

      tmv::BandMatrix<std::complex<T> > cL(N,N,nlo,0);
      tmv::DiagMatrix<std::complex<T> > cD(N);
      tmv::BandMatrix<std::complex<T> > cLDL(N,N,nlo,nlo);
      if (mattype == 3) {
	c.DivideUsing(tmv::CH);
	c.SetDiv();
	cL = c.CHD().GetL();
	cD = c.CHD().GetD();
	cLDL = cL*cD*cL.Adjoint();
	Assert(Norm(c-cLDL) < ceps*Norm(c),"HermBand C LDL");
      }

      tmv::HermBandMatrix<std::complex<T>,uplo,stor> c2 = c;
      LDL_Decompose(c2.View());
      (cL = c2.LowerBand()).diag().SetAllTo(T(1));
      cD = DiagMatrixViewOf(c2.diag());
      cLDL = cL*cD*cL.Adjoint();
      Assert(Norm(c-cLDL) < ceps*Norm(c),"HermBand C LDL2");
      Assert(Norm(cD.Imag()) < ceps*Norm(c),"HermBand C LDL2 - D");

      // Alt calculation:
      tmv::DiagMatrix<std::complex<T> > cE(c2.diag(-1));
      cLDL.Zero();
      cLDL.diag(0,1,N) = (cE*cD.SubDiagMatrix(0,N-1)*cE.Conjugate()).diag();
      cLDL.diag(-1) = (cE*cD.SubDiagMatrix(0,N-1)).diag();
      cLDL += cD;
      cLDL.diag(1) = cLDL.diag(-1).Conjugate();
      Assert(Norm(c-cLDL) < ceps*Norm(c),"HermBand C LDL2 Alt");

      c2.Conjugate() = c;
      LDL_Decompose(c2.Conjugate());
      (cL = c2.Conjugate().LowerBand()).diag().SetAllTo(T(1));
      cD = DiagMatrixViewOf(c2.Conjugate().diag());
      cLDL = cL*cD*cL.Adjoint();
      Assert(Norm(c-cLDL) < ceps*Norm(c),"HermBand C LDL3");
      Assert(Norm(cD.Imag()) < ceps*Norm(c),"HermBand C LDL2 - D");
      // The Lapack version throws whenever mattype is not posdef, 
      // but native algorithm succeeds when mattype is indefinite.
      Assert(mattype%3 != 2,"didn't throw NonPosDef, mattype != singular"); 
    }
    catch (tmv::NonPosDef)
    { Assert(mattype%3 != 0,"caught NonPosDef, mattype != posdef"); }

    // SV Decomposition
    //std::cout<<"SVD"<<std::endl;
    {
      m.DivideUsing(tmv::SV);
      m.SaveDiv();
      m.SetDiv();
      tmv::Matrix<T> U = m.SVD().GetU();
      tmv::DiagMatrix<T> S = m.SVD().GetS();
      tmv::Matrix<T> V = m.SVD().GetV();
      Assert(Norm(m-U*S*V) < eps*Norm(m),"HermBand SV");

      tmv::Matrix<T> U2(N,N);
      tmv::DiagMatrix<T> S2(N);
      tmv::Matrix<T> V2(N,N);
      SV_Decompose(m,U2.View(),S2.View(),V2.View());
      Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"HermBand SV2");

      SV_Decompose(m,S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"HermBand SV3");
      SV_Decompose(m,U2.View(),S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"HermBand SV4 S");
      Assert(Norm(U2*S2*S2*U2.Adjoint()-m*m.Adjoint()) < 
	  eps*Norm(m)*Norm(m),"HermBand SV4 U");
      SV_Decompose(m,S2.View(),V2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"HermBand SV5 S");
      Assert(Norm(m.Adjoint()*m-V2.Adjoint()*S2*S2*V2) < 
	  eps*Norm(m)*Norm(m),"HermBand SV5 V");

      c.DivideUsing(tmv::SV);
      c.SaveDiv();
      c.SetDiv();
      tmv::Matrix<std::complex<T> > cU = c.SVD().GetU();
      S = c.SVD().GetS();
      tmv::Matrix<std::complex<T> > cV = c.SVD().GetV();
      Assert(Norm(c-cU*S*cV) < ceps*Norm(c),"HermBand C SV");

      tmv::Matrix<std::complex<T> > cU2(N,N);
      tmv::Matrix<std::complex<T> > cV2(N,N);
      SV_Decompose(c,cU2.View(),S2.View(),cV2.View());
      Assert(Norm(c-cU2*S2*cV2) < ceps*Norm(c),"HermBand C SV2");

      SV_Decompose(c,S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"HermBand C SV3");
      SV_Decompose(c,cU2.View(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"HermBand C SV4 S");
      Assert(Norm(cU2*S2*S2*cU2.Adjoint()-c*c.Adjoint()) < 
	  ceps*Norm(c)*Norm(c),"HermBand C SV4 U");
      SV_Decompose(c,S2.View(),cV2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"HermBand C SV5 S");
      Assert(Norm(c.Adjoint()*c-cV2.Adjoint()*S2*S2*cV2) < 
	  ceps*Norm(c)*Norm(c),"Herm C SV5 V");

      SV_Decompose(c.Conjugate(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"HermBand C SV6");
      SV_Decompose(c.Conjugate(),cU2.View(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"HermBand C SV7 S");
      Assert(Norm(cU2*S2*S2*cU2.Adjoint()-c.Conjugate()*c.Transpose()) < 
	  ceps*Norm(c)*Norm(c),"HermBand C SV7 U");
      SV_Decompose(c.Conjugate(),cU2.Conjugate(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"HermBand C SV8 S");
      Assert(Norm(cU2.Conjugate()*S2*S2*cU2.Transpose()-
	    c.Conjugate()*c.Transpose()) < ceps*Norm(c)*Norm(c),
	  "HermBand C SV8 U");
    }
    
    // Eigen
    //std::cout<<"Eigen"<<std::endl;
    {
      tmv::Matrix<T> V(N,N);
      tmv::Vector<T> L(N);
      Eigen(m,V.View(),L.View());
      Assert(Norm(m*V-V*DiagMatrixViewOf(L)) < eps*Norm(m),"HermBand Eigen");

      tmv::Vector<T> L2(N);
      Eigen(m,L2.View());
      Assert(Norm(L2-L) < eps*Norm(L),"HermBand Eigen2");

      tmv::Matrix<std::complex<T> > cV(N,N);
      Eigen(c,cV.View(),L.View());
      Assert(Norm(c*cV-cV*DiagMatrixViewOf(L)) < ceps*Norm(c),
	  "HermBand C Eigen");

      Eigen(c,L2.View());
      Assert(Norm(L2-L) < ceps*Norm(L),"HermBand C Eigen2");

      Eigen(c,cV.Conjugate(),L.View());
      Assert(Norm(c*cV.Conjugate()-cV.Conjugate()*DiagMatrixViewOf(L)) < 
	  ceps*Norm(c),"HermBand C Eigen3");

      Eigen(c.Conjugate(),L2.View());
      Assert(Norm(L2-L) < ceps*Norm(L),"HermBand C Eigen4");
    }

    // Square Root
    //std::cout<<"Square root:"<<std::endl;
    try {
      tmv::HermMatrix<T,uplo,stor> S(N);
      SquareRoot(m,S.View());
      //std::cout<<"SquareRoot of "<<m<<" is "<<S<<std::endl;
      Assert(Norm(m-S*S) < eps*Norm(m),"HermBand Square Root");

      tmv::HermMatrix<std::complex<T>,uplo,stor> cS(N);
      SquareRoot(c,cS.View());
      //std::cout<<"SquareRoot of "<<c<<" is "<<cS<<std::endl;
      Assert(Norm(c-cS*cS) < eps*Norm(c),"HermBand C Square Root");

      SquareRoot(c,cS.Conjugate());
      Assert(Norm(c-cS.Conjugate()*cS.Conjugate()) < ceps*Norm(c),
	  "HermBand C Square Root2");

      SquareRoot(c,cS.Transpose());
      Assert(Norm(c-cS.Transpose()*cS.Transpose()) < ceps*Norm(c),
	  "HermBand C Square Root3");

      SquareRoot(c,cS.Adjoint());
      Assert(Norm(c-cS.Adjoint()*cS.Adjoint()) < ceps*Norm(c),
	  "HermBand C Square Root4");

      SquareRoot(c.Conjugate(),cS.View());
      Assert(Norm(c.Conjugate()-cS*cS) < ceps*Norm(c),"HermBand C Square Root5");

      SquareRoot(c.Conjugate(),cS.Conjugate());
      Assert(Norm(c.Conjugate()-cS.Conjugate()*cS.Conjugate()) < ceps*Norm(c),
	  "HermBand C Square Root6");

      SquareRoot(c.Conjugate(),cS.Transpose());
      Assert(Norm(c.Conjugate()-cS.Transpose()*cS.Transpose()) < ceps*Norm(c),
	  "HermBand C Square Root7");

      SquareRoot(c.Conjugate(),cS.Adjoint());
      Assert(Norm(c.Conjugate()-cS.Adjoint()*cS.Adjoint()) < ceps*Norm(c),
	  "HermBand C Square Root8");
      Assert(mattype%3 == 0,"didn't throw NonPosDef, mattype == pos def"); 
    }
    catch (tmv::NonPosDef)
    { Assert(mattype%3 != 0,"caught NonPosDef, mattype != posdef"); }
  }
}

template <class T, tmv::UpLoType uplo, tmv::StorageType stor> 
void TestSymBandDecomp()
{
  for (int mattype = 0; mattype < 6; mattype++) {
    //std::cout<<"mattype = "<<mattype<<std::endl;
    // mattype = 0  is PosDef, nlo = 3
    // mattype = 1  is Indef, nlo = 3
    // mattype = 2  is Singular, nlo = 3
    // mattype = 3  is PosDef, nlo = 1
    // mattype = 4  is Indef, nlo = 1
    // mattype = 5  is Singular, nlo = 1

    const int N = 8;
    int nlo = 3;
    if (mattype >= 3) nlo = 1;

    tmv::SymBandMatrix<T,uplo,stor> m(N,nlo);
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
      if ((uplo == tmv::Upper && j>=i && j<=i+nlo) || 
	  (uplo == tmv::Lower && i>=j && i<=j+nlo)) 
	m(i,j) = T(2.+4*i-5*j);
    if (mattype%3 == 0) {
      m /= T(100);
      m.diag().AddToAll(T(1));
      for(int i=0;i<N;i++) m.diag()(i) += T(i);
    } else if (mattype%3 == 1) {
      m /= T(8);
      m.diag(0,2,4).Zero();
      m(2,1) = T(50);
      m.diag(1,4,7).AddToAll(T(10));
      m.row(7,7-nlo,7).AddToAll(T(10));
    } else {
      if (nlo > 2) m.row(2,0,2).Zero();
      else m.row(2,2-nlo,2).Zero();
      m.row(2,2,2+nlo+1).Zero();
      m.row(4,5-nlo,4) = m.row(5,5-nlo,4);
      m.row(4,5,4+nlo+1) = m.row(5,5,4+nlo+1);
      if (4-nlo >= 0) m(4,4-nlo) = T(0);
      if (5+nlo < N) m(5,5+nlo) = T(0);
      m(4,5) = m(4,4) = m(5,5);
    }
    //std::cout<<"m = "<<m<<std::endl;

    tmv::SymBandMatrix<std::complex<T>,uplo,stor> c(N,nlo);
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
      if ((uplo == tmv::Upper && j>=i && j<=i+nlo) || 
	  (uplo == tmv::Lower && i>=j && i<=j+nlo)) 
	c(i,j) = std::complex<T>(2.+4*i-5*j,3.-i);
    if (mattype%3 == 0) {
      c /= T(100);
      c.diag().AddToAll(T(1));
      for(int i=0;i<N;i++) c.diag()(i) += T(i);
    } else if (mattype%3 == 1) {
      c /= T(8);
      c.diag(0,2,4).Zero();
      c(2,1) = T(50);
      c.diag(1,4,7).AddToAll(T(10));
      c.row(7,7-nlo,7).AddToAll(T(10));
    } else {
      if (nlo > 2) c.row(2,0,2).Zero();
      else c.row(2,2-nlo,2).Zero();
      c.row(2,2,2+nlo+1).Zero();
      c.row(4,5-nlo,4) = c.row(5,5-nlo,4);
      c.row(4,5,4+nlo+1) = c.row(5,5,4+nlo+1);
      if (4-nlo >= 0) c(4,4-nlo) = T(0);
      if (5+nlo < N) c(5,5+nlo) = T(0);
      c(4,5) = c(4,4) = c(5,5);
    }
    //std::cout<<"c = "<<c<<std::endl;

    T eps = EPS;
    T ceps = EPS;
    if (mattype %3 != 2) {
      eps *= Norm(m) * Norm(m.Inverse());
      ceps *= Norm(c) * Norm(c.Inverse());
    } else {
      eps *= T(100);
      ceps *= T(100);
    }

    // LDL Decomposition
    if (mattype /3 == 1) try {
      tmv::BandMatrix<T> L(N,N,1,0);
      tmv::DiagMatrix<T> D(N);
      tmv::BandMatrix<T> LDL(N,N,1,1);

      tmv::SymBandMatrix<T,uplo,stor> m2 = m;
      LDL_Decompose(m2.View());
      (L = m2.LowerBand()).diag().SetAllTo(T(1));
      D = DiagMatrixViewOf(m2.diag());
      LDL = L*D*L.Transpose();
      Assert(Norm(m-LDL) < eps*Norm(m),"SymBand LDL2");

      // Alt calculation:
      // (I + L) D (I + Lt)
      // D + LD + DLt = LDLt
      tmv::DiagMatrix<T> E(m2.diag(-1));
      LDL.Zero();
      LDL.diag(0,1,N) = (E*D.SubDiagMatrix(0,N-1)*E).diag();
      LDL.diag(-1) = (E*D.SubDiagMatrix(0,N-1)).diag();
      LDL += D;
      LDL.diag(1) = LDL.diag(-1);
      Assert(Norm(m-LDL) < eps*Norm(m),"SymBand LDL2 Alt");

      tmv::BandMatrix<std::complex<T> > cL(N,N,nlo,0);
      tmv::DiagMatrix<std::complex<T> > cD(N);
      tmv::BandMatrix<std::complex<T> > cLDL(N,N,nlo,nlo);

      tmv::SymBandMatrix<std::complex<T>,uplo,stor> c2 = c;
      LDL_Decompose(c2.View());
      (cL = c2.LowerBand()).diag().SetAllTo(T(1));
      cD = DiagMatrixViewOf(c2.diag());
      cLDL = cL*cD*cL.Transpose();
      Assert(Norm(c-cLDL) < ceps*Norm(c),"SymBand C LDL2");

      // Alt calculation:
      tmv::DiagMatrix<std::complex<T> > cE(c2.diag(-1));
      cLDL.Zero();
      cLDL.diag(0,1,N) = (cE*cD.SubDiagMatrix(0,N-1)*cE).diag();
      cLDL.diag(-1) = (cE*cD.SubDiagMatrix(0,N-1)).diag();
      cLDL += cD;
      cLDL.diag(1) = cLDL.diag(-1);
      Assert(Norm(c-cLDL) < eps*Norm(m),"SymBand C LDL2 Alt");

      c2.Conjugate() = c;
      LDL_Decompose(c2.Conjugate());
      (cL = c2.Conjugate().LowerBand()).diag().SetAllTo(T(1));
      cD = DiagMatrixViewOf(c2.Conjugate().diag());
      cLDL = cL*cD*cL.Transpose();
      Assert(Norm(c-cLDL) < ceps*Norm(c),"SymBand C LDL3");
      Assert(mattype%3 != 2,"didn't throw NonPosDef, mattype != singular"); 
    }
    catch (tmv::NonPosDef)
    { Assert(mattype%3 != 0,"caught NonPosDef, mattype != posdef"); }

    // SV Decomposition
    {
      m.DivideUsing(tmv::SV);
      m.SaveDiv();
      m.SetDiv();
      tmv::Matrix<T> U = m.SVD().GetU();
      tmv::DiagMatrix<T> S = m.SVD().GetS();
      tmv::Matrix<T> V = m.SVD().GetV();
      Assert(Norm(m-U*S*V) < eps*Norm(m),"SymBand SV");

      tmv::Matrix<T> U2(N,N);
      tmv::DiagMatrix<T> S2(N);
      tmv::Matrix<T> V2(N,N);
      SV_Decompose(m,U2.View(),S2.View(),V2.View());
      Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"SymBand SV2");

      SV_Decompose(m,S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"SymBand SV3");
      SV_Decompose(m,U2.View(),S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"SymBand SV4 S");
      Assert(Norm(U2-U) < eps*Norm(U),"SymBand SV4 U");
      SV_Decompose(m,S2.View(),V2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"SymBand SV5 S");
      Assert(Norm(m.Adjoint()*m-V2.Adjoint()*S2*S2*V2) < eps*Norm(m)*Norm(m),
	  "SymBand SV5 V");

      c.DivideUsing(tmv::SV);
      c.SaveDiv();
      c.SetDiv();
      tmv::Matrix<std::complex<T> > cU = c.SymSVD().GetU();
      S = c.SymSVD().GetS();
      tmv::Matrix<std::complex<T> > cV = c.SymSVD().GetV();
      Assert(Norm(c-cU*S*cV) < ceps*Norm(c),"SymBand C SV");

      tmv::Matrix<std::complex<T> > cU2(N,N);
      tmv::Matrix<std::complex<T> > cV2(N,N);
      SV_Decompose(c,cU2.View(),S2.View(),cV2.View());
      Assert(Norm(c-cU2*S2*cV2) < ceps*Norm(c),"SymBand C SV2");

      SV_Decompose(c,S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"SymBand C SV3");
      SV_Decompose(c,cU2.View(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"SymBand C SV4 S");
      Assert(Norm(cU2-cU) < ceps*Norm(cU),"SymBand C SV4 U");
      SV_Decompose(c,S2.View(),cV2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"SymBand C SV5 S");
      Assert(Norm(c.Adjoint()*c-cV2.Adjoint()*S2*S2*cV2) < ceps*Norm(c)*Norm(c),
	  "SymBand C SV5 V");

      SV_Decompose(c.Conjugate(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"SymBand C SV6");
      SV_Decompose(c.Conjugate(),cU2.View(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"SymBand C SV7 S");
      Assert(Norm(cU2-cU.Conjugate()) < ceps*Norm(cU),"SymBand C SV7 U");
      SV_Decompose(c.Conjugate(),cU2.Conjugate(),S2.View());
      Assert(Norm(S2-S) < ceps*Norm(S),"SymBand C SV8 S");
      Assert(Norm(cU2.Conjugate()-cU.Conjugate()) < ceps*Norm(cU),
	  "SymBand C SV8 U");
    }
  }
}

template <class T> void TestAllSymBandDiv()
{
  TestHermBandDecomp<T,tmv::Upper,tmv::ColMajor>();
  TestHermBandDecomp<T,tmv::Upper,tmv::RowMajor>();
  TestHermBandDecomp<T,tmv::Lower,tmv::ColMajor>();
  TestHermBandDecomp<T,tmv::Lower,tmv::RowMajor>();
  TestSymBandDecomp<T,tmv::Upper,tmv::ColMajor>();
  TestSymBandDecomp<T,tmv::Upper,tmv::RowMajor>();
  TestSymBandDecomp<T,tmv::Lower,tmv::ColMajor>();
  TestSymBandDecomp<T,tmv::Lower,tmv::RowMajor>();
  TestSymBandDiv<T>(tmv::CH,PosDef);
  TestSymBandDiv<T>(tmv::LU,PosDef);
  TestSymBandDiv<T>(tmv::LU,InDef);
  TestSymBandDiv<T>(tmv::SV,PosDef);
  TestSymBandDiv<T>(tmv::SV,InDef);
  TestSymBandDiv<T>(tmv::SV,Sing);
}

#ifdef INST_DOUBLE
template void TestAllSymBandDiv<double>();
#endif
#ifdef INST_FLOAT
template void TestAllSymBandDiv<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllSymBandDiv<long double>();
#endif
#ifdef INST_INT
template void TestAllSymBandDiv<int>();
#endif
