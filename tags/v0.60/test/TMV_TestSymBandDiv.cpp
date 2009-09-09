
#define START 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestSymBandArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> inline bool IsPosDef(const tmv::GenSymBandMatrix<T>& m)
{
  try {
    tmv::ConstSymBandMatrixView<T> mview = m;
    mview.DivideUsing(tmv::CH);
    mview.SetDiv();
    std::ostream* checkout = showdiv ? &std::cout : 0;
    Assert(mview.CheckDecomp(checkout),"CheckDecomp mview in IsPosDef"); 
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

template <class T> inline void TestSymBandDiv(tmv::DivType dt, PosDefCode pdc)
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
    T eps = EPS*Norm(m)*Norm(m.Inverse());
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

    if (pdc != Sing) {
      if (showacc) {
	std::cout<<"si.Det = "<<si.Det()<<", m.Det = "<<m.Det()<<std::endl;
	std::cout<<"abs(sdet-mdet) = "<<std::abs(si.Det()-m.Det());
	std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(m.Det())<<std::endl;
	std::cout<<"abs(abs(sdet)-abs(mdet)) = "<<std::abs(std::abs(si.Det())-std::abs(m.Det()));
	std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(m.Det())<<std::endl;
      }
      Assert(std::abs(m.Det()-si.Det()) < eps*std::abs(m.Det()),"SymBand Det");
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

    T ceps = EPS*Norm(cm)*Norm(cm.Inverse());

    if (pdc != Sing) {
      if (showacc) {
	std::cout<<"csi.Det = "<<csi.Det()<<", cm.Det = "<<cm.Det()<<std::endl;
	std::cout<<"abs(csidet-cmdet) = "<<std::abs(csi.Det()-cm.Det());
	std::cout<<"  csidet/cmdet = "<<csi.Det()/cm.Det();
	std::cout<<"  EPS*abs(cmdet) = "<<ceps*std::abs(cm.Det())<<std::endl;
	std::cout<<"abs(abs(csdet)-abs(cmdet)) = "<<std::abs(std::abs(csi.Det())-std::abs(cm.Det()));
	std::cout<<"  EPS*abs(cmdet) = "<<ceps*std::abs(cm.Det())<<std::endl;
      }
      Assert(std::abs(csi.Det()-cm.Det()) < ceps*std::abs(cm.Det()),"SymBand CDet");
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

template <class T> void TestAllSymBandDiv()
{
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
