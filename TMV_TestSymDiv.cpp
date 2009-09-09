
#define START 0

#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_TestMatrixDivArith.h"

enum PosDefCode { PosDef, InDef, Sing };
inline std::string PDLabel(PosDefCode pdc)
{
  if (pdc == PosDef) return "Positive Definite";
  else if (pdc == InDef) return "Indefinite";
  else return "Singular";
}

template <class T> inline bool IsPosDef(const tmv::GenSymMatrix<T>& m)
{
  try {
    tmv::ConstSymMatrixView<T> mview = m;
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

template <class T> inline void TestSymDiv(tmv::DivType dt, PosDefCode pdc)
{
  const int N = 10;

  std::vector<tmv::SymMatrixView<T> > s;
  std::vector<tmv::SymMatrixView<std::complex<T> > > cs;

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+0.1*i-0.2*j;
  tmv::Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-0.3*i+0.1*j;
  if (pdc == PosDef) {
    a1 /= T(N*N);
    a2 /= T(N*N);
    a1.diag().AddToAll(T(1));
    a2.diag().AddToAll(T(1));
    for(int i=0;i<N;++i) a1.diag()(i) += T(i);
    for(int i=0;i<2*N;++i) a2.diag()(i) += T(i);
  } else if (pdc == InDef) {
    a1 /= T(N);
    a2 /= T(N);
    a1.diag(0,2,4).Zero();
    a2.diag(0,2,4).Zero();
    a1(4,0) = a1(0,4) = T(50);
    a1.diag(-1,6,9).AddToAll(T(10));
    a2.diag(-1,6,9).AddToAll(T(10));
    a1.diag(1,6,9).AddToAll(T(10));
    a2.diag(1,6,9).AddToAll(T(10));
    a1.row(9,0,5).AddToAll(T(10));
    a2.row(9,0,5).AddToAll(T(10));
    a1(3,3) = T(0);
    a2(3,3) = T(0);
    if (N > 10) {
      a1.diag(0,10,N) *= T(0.0001);
      a2.diag(0,10,N) *= T(0.0001);
    }
  } else {
    a1.row(2).Zero();
    a1.col(2).Zero();
    a1.row(4) = a1.row(5);
    a1.col(4) = a1.col(5);
    a1(4,5) = a1(5,4) = a1(5,5) = a1(4,4);
    a2.row(2).Zero();
    a2.col(2).Zero();
    a2.row(4) = a2.row(5);
    a2.col(4) = a2.col(5);
    a2(4,5) = a2(5,4) = a2(5,5) = a2(4,4);
  }
  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);
  //ca1.diag().Imag().Zero();
  tmv::Matrix<std::complex<T> > ca2 = a2 * std::complex<T>(3,-4);
  //ca2.diag().Imag().Zero();

  tmv::Vector<T> v1(N);
  tmv::Vector<T> v2(N-1);
  for (int i=0; i<N; ++i) v1(i) = 16.-3*i; 
  for (int i=0; i<N-1; ++i) v2(i) = -7.+2*i; 

  tmv::HermMatrix<T,tmv::Upper,tmv::RowMajor> H1(a1);
  tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> CH1(ca1);
  s.push_back(H1.View());
  cs.push_back(CH1.View());
#ifdef XTEST
  tmv::HermMatrix<T,tmv::Upper,tmv::ColMajor> H2(a1);
  tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> CH2(ca1);
  s.push_back(H2.View());
  cs.push_back(CH2.View());
  tmv::HermMatrix<T,tmv::Lower,tmv::RowMajor> H3(a1);
  tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> CH3(ca1);
  s.push_back(H3.View());
  cs.push_back(CH3.View());
  tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor> H4(a1);
  tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> CH4(ca1);
  s.push_back(H4.View());
  cs.push_back(CH4.View());
#endif

  // These need to be outside of (if dt!=CH loop), otherwise they will
  // go out of scope at the end of it and the Views will refer to
  // defunct memory.
  tmv::SymMatrix<T,tmv::Upper,tmv::RowMajor> S1(a1);
  tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> CS1(ca1);
#ifdef XTEST
  tmv::SymMatrix<T,tmv::Upper,tmv::ColMajor> S2(a1);
  tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> CS2(ca1);
  tmv::SymMatrix<T,tmv::Lower,tmv::RowMajor> S3(a1);
  tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> CS3(ca1);
  tmv::SymMatrix<T,tmv::Lower,tmv::ColMajor> S4(a1);
  tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> CS4(ca1);
#endif
  if (dt != tmv::CH) {
    s.push_back(S1.View());
    cs.push_back(CS1.View());
#ifdef XTEST
    s.push_back(S2.View());
    cs.push_back(CS2.View());
    s.push_back(S3.View());
    cs.push_back(CS3.View());
    s.push_back(S4.View());
    cs.push_back(CS4.View());
#endif
  }

#ifdef XTEST
  tmv::HermMatrix<T,tmv::Upper> H5(a2);
  tmv::HermMatrix<std::complex<T>,tmv::Upper> CH5(ca2);
  tmv::HermMatrix<T,tmv::Lower> H6(a2);
  tmv::HermMatrix<std::complex<T>,tmv::Lower> CH6(ca2);
  tmv::SymMatrix<T,tmv::Upper> S5(a2);
  tmv::SymMatrix<std::complex<T>,tmv::Upper> CS5(ca2);
  tmv::SymMatrix<T,tmv::Lower> S6(a2);
  tmv::SymMatrix<std::complex<T>,tmv::Lower> CS6(ca2);
  if (doallarith) {
    s.push_back(H5.SubSymMatrix(0,2*N,2));
    cs.push_back(CH5.SubSymMatrix(0,2*N,2));
    s.push_back(H6.SubSymMatrix(0,2*N,2));
    cs.push_back(CH6.SubSymMatrix(0,2*N,2));

    if (dt != tmv::CH) {
      s.push_back(S5.SubSymMatrix(0,2*N,2));
      cs.push_back(CS5.SubSymMatrix(0,2*N,2));
      s.push_back(S6.SubSymMatrix(0,2*N,2));
      cs.push_back(CS6.SubSymMatrix(0,2*N,2));
    }
  }
#endif

  size_t ntot = s.size();

  tmv::Matrix<T> a3 = a1.Cols(0,N/2);
  tmv::Matrix<std::complex<T> > ca3 = ca1.Cols(0,N/2);
  tmv::Matrix<T> a4 = a1.Rows(0,N/2);
  tmv::Matrix<std::complex<T> > ca4 = ca1.Rows(0,N/2);
  tmv::Matrix<T> a5 = a1.Cols(0,0);
  tmv::Matrix<std::complex<T> > ca5 = ca1.Cols(0,0);
  tmv::Matrix<T> a6 = a1.Rows(0,0);
  tmv::Matrix<std::complex<T> > ca6 = ca1.Rows(0,0);
  tmv::HermMatrix<T> s1 = H1;
  tmv::HermMatrix<std::complex<T> > cs1 = s1;
  cs1.LowerTri().OffDiag() *= std::complex<T>(2,-1);

  std::ostream* checkout = showdiv ? &std::cout : 0;

  for(size_t i=START;i<ntot;i++) {
    tmv::SymMatrixView<T> si = s[i];
    tmv::SymMatrixView<std::complex<T> > csi = cs[i];
    si.SaveDiv();
    csi.SaveDiv();
    if (showstartdone)
      std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::Type(si)<<"  "<<si<<std::endl;

    Assert(IsPosDef(si) == (pdc==PosDef),"IsPosDef");
    if (csi.isherm()) {
      Assert(IsPosDef(csi) == (pdc==PosDef),"IsPosDef");
    }

    tmv::Matrix<T> m = si;
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
      std::cout<<"v/b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(v1)<<std::endl;
    }
    Assert(Norm(x1-x2) < eps*Norm(v1),"Sym v/b");

    x1 = v1%si;
    x2 = v1%m;
    if (showacc)
      std::cout<<"v%b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(v1)<<std::endl;
    Assert(Norm(x1-x2) < eps*Norm(v1),"Sym v%b");

    tmv::Matrix<T,tmv::ColMajor> sinv = si.Inverse();
    tmv::Matrix<T,tmv::ColMajor> minv = m.Inverse();
    if (showacc) {
      std::cout<<"sinv = "<<sinv<<std::endl;
      std::cout<<"minv = "<<minv<<std::endl;
      std::cout<<"Norm(minv-sinv) = "<<Norm(minv-sinv)<<"  "<<eps*Norm(sinv)<<std::endl;
    }
    Assert(Norm(sinv-minv) < eps*Norm(sinv),"Sym Inverse");

    if (pdc != Sing) {
      if (showacc) {
	std::cout<<"si.Det = "<<si.Det()<<", m.Det = "<<m.Det()<<std::endl;
	std::cout<<"abs(sdet-mdet) = "<<std::abs(si.Det()-m.Det());
	std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(m.Det())<<std::endl;
	std::cout<<"abs(abs(sdet)-abs(mdet)) = "<<std::abs(std::abs(si.Det())-std::abs(m.Det()));
	std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(m.Det())<<std::endl;
      }
      Assert(std::abs(m.Det()-si.Det()) < eps*std::abs(m.Det()),"Sym Det");
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
      Assert(std::abs(csi.Det()-cm.Det()) < ceps*std::abs(cm.Det()),"Sym CDet");
    }

    tmv::Vector<std::complex<T> > cv(v1 * std::complex<T>(1,1));
    cv(1) += std::complex<T>(-1,5);
    cv(2) -= std::complex<T>(-1,5);

    // test real / complex
    tmv::Vector<std::complex<T> > y1 = v1/csi;
    tmv::Vector<std::complex<T> > y2 = v1/cm;
    if (showacc) 
      std::cout<<"v/cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(v1)<<std::endl;
    Assert(Norm(y1-y2) < ceps*Norm(v1),"Sym v/cs");

    // test complex / real
    y1 = cv/si;
    y2 = cv/m;
    if (showacc) 
      std::cout<<"cv/b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(cv)<<std::endl;
    Assert(Norm(y1-y2) < eps*Norm(cv),"Sym cv/b");

    // test complex / complex
    y1 = cv/csi;
    y2 = cv/cm;
    if (showacc) 
      std::cout<<"cv/cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(cv)<<std::endl;
    Assert(Norm(y1-y2) < ceps*Norm(cv),"Sym cv/cs");

    y1 = v1%csi;
    y2 = v1%cm;
    if (showacc) 
      std::cout<<"v%cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(v1)<<std::endl;
    Assert(Norm(y1-y2) < ceps*Norm(v1),"Sym v%cs");

    y1 = cv%si;
    y2 = cv%m;
    if (showacc) 
      std::cout<<"cv%b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(cv)<<std::endl;
    Assert(Norm(y1-y2) < eps*Norm(cv),"Sym cv%b");
    y1 = cv%csi;
    y2 = cv%cm;
    if (showacc) 
      std::cout<<"cv%cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(cv)<<std::endl;
    Assert(Norm(y1-y2) < ceps*Norm(cv),"Sym cv%cs");


#ifdef XTEST
    if (pdc == PosDef)
      TestMatrixDivArith<T>(dt==tmv::CH?tmv::LU:dt,a1.View(),si,ca1.View(),csi,
	  "Sym/SquareMatrix");
#endif
    if (pdc != Sing) {
      TestMatrixDivArith<T>(dt,si,a1.View(),csi,ca1.View(),"SquareMatrix/Sym");
#ifdef XTEST
      TestMatrixDivArith<T>(dt,si,a3.View(),csi,ca3.View(),
	  "NonSquareMatrix/Sym");
      TestMatrixDivArith<T>(dt,si,a4.View(),csi,ca4.View(),
	  "NonSquareMatrix/Sym");
      TestMatrixDivArith<T>(dt,si,a5.View(),csi,ca5.View(),
	  "DegenerateMatrix/Sym");
      TestMatrixDivArith<T>(dt,si,a6.View(),csi,ca6.View(),
	  "DegenerateMatrix/Sym");
      TestMatrixDivArith<T>(dt,si,s1.View(),csi,cs1.View(),"Sym/Sym");
#endif
    }
  }

  std::cout<<PDLabel(pdc)<<" SymMatrix<"<<tmv::Type(T())<<"> Division using ";
  std::cout<<tmv::Text(dt)<<" passed all tests\n";
}

template <class T> void TestAllSymDiv()
{
  TestSymDiv<T>(tmv::CH,PosDef);
  TestSymDiv<T>(tmv::LU,PosDef);
  TestSymDiv<T>(tmv::LU,InDef);
  TestSymDiv<T>(tmv::SV,PosDef);
  TestSymDiv<T>(tmv::SV,InDef);
  TestSymDiv<T>(tmv::SV,Sing);
}

template void TestAllSymDiv<double>();
#ifndef NOFLOAT
template void TestAllSymDiv<float>();
#endif
#ifdef LONGDOUBLE
template void TestAllSymDiv<long double>();
#endif
