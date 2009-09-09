// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_TestSymArith.h"

template <class T> static bool IsPosDef(const tmv::GenSymMatrix<T>& m)
{
#ifdef NOTHROW
  for(size_t i=1;i<=m.size();i++) {
    T d = m.SubSymMatrix(0,i).Det();
    if (tmv::REAL(d) < 0) return false;
  }
  return true;
#else
  try {
    tmv::HermMatrix<T> m2 = m;
    CH_Decompose(m2.View());
  }
  catch (tmv::NonPosDef) {
#ifdef XTEST
    if (showacc) std::cout<<"caught nonposdef\n";
    //std::cout<<"m.size ="<<m.size()<<std::endl;
    for(size_t i=1;i<=m.size();i++) {
      //std::cout<<"i = "<<i<<std::endl;
      T d = tmv::Matrix<T>(m.SubSymMatrix(0,i)).Det();
      if (showacc) std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
      if (!(tmv::REAL(d) > 0)) {
        //std::cout<<"!>0 "<<tmv::REAL(d)<<std::endl;
        return false;
      } else {
        //std::cout<<">0 "<<tmv::REAL(d)<<std::endl;
      }
    }
    //std::cout<<"m = "<<TypeText(m)<<"  "<<m<<std::endl;
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
      std::cout<<"m = "<<TypeText(m)<<"  "<<m<<std::endl;
      std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
      Assert(tmv::FALSE,"Didn't catch NonPosDef, but determinant of sub-block is negative");
    }
  }
#endif
  return true;
#endif
}

template <class T> void TestSymDiv(tmv::DivType dt, PosDefCode pdc)
{
  const int N = 10;

  std::vector<tmv::SymMatrixView<T> > s;
  std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
  std::vector<tmv::BaseMatrix<T>*> B;
  std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
  MakeSymList(s,cs,B,CB,pdc);

  tmv::Vector<T> v1(N);
  tmv::Vector<T> v2(N-1);
  for (int i=0; i<N; ++i) v1(i) = T(16-3*i); 
  for (int i=0; i<N-1; ++i) v2(i) = T(-7+2*i); 

  std::ostream* checkout = showdiv ? &std::cout : 0;

  for(size_t i=START;i<s.size();i++) {
    tmv::SymMatrixView<T> si = s[i];
    tmv::SymMatrixView<std::complex<T> > csi = cs[i];
    if (dt == tmv::CH && csi.issym()) continue;
    si.SaveDiv();
    csi.SaveDiv();
    if (showstartdone)
      std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TypeText(si)<<"  "<<si<<std::endl;

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

  if (pdc != Sing) {
    TestSymDiv_A<T>(dt,pdc);
    TestSymDiv_B1<T>(dt,pdc);
    TestSymDiv_C1<T>(dt,pdc);
    TestSymDiv_D1<T>(dt,pdc);
    TestSymDiv_E1<T>(dt,pdc);
  }
  if (pdc == PosDef) {
    if (dt != tmv::CH) TestSymDiv_B2<T>(dt,pdc);
    if (dt == tmv::LU) TestSymDiv_C2<T>(dt,pdc);
    if (dt == tmv::LU) TestSymDiv_D2<T>(dt,pdc);
    if (dt != tmv::CH) TestSymDiv_E2<T>(dt,pdc);
  }
  for(size_t i=0;i<B.size();++i) delete B[i];
  for(size_t i=0;i<CB.size();++i) delete CB[i];

  std::cout<<PDLabel(pdc)<<" SymMatrix<"<<tmv::TypeText(T())<<"> Division using ";
  std::cout<<tmv::Text(dt)<<" passed all tests\n";
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
  std::cout<<"SymMatrix<"<<tmv::TypeText(T())<<"> passed all ";
  std::cout<<"decomposition tests.\n";
  TestSymDiv<T>(tmv::CH,PosDef);
  TestSymDiv<T>(tmv::LU,PosDef);
  TestSymDiv<T>(tmv::LU,InDef);
  TestSymDiv<T>(tmv::SV,PosDef);
  TestSymDiv<T>(tmv::SV,InDef);
  TestSymDiv<T>(tmv::SV,Sing);
}

#ifdef TEST_DOUBLE
template void TestAllSymDiv<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllSymDiv<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllSymDiv<long double>();
#endif
