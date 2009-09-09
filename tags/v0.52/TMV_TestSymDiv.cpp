//#define SHOWACC
//#define SHOWTESTS
//#define SHOWCHECK
//#define DOALLARITH

#define STARTAT 0

#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_TestMatrixDiv.h"

using tmv::Matrix;
using tmv::SymMatrix;
using tmv::HermMatrix;
using tmv::SymMatrixView;
using tmv::RowMajor;
using tmv::ColMajor;
using tmv::Upper;
using tmv::Lower;

enum PosDefCode { PosDef, InDef, Sing };
inline string PDLabel(PosDefCode pdc)
{
  if (pdc == PosDef) return "Positive Definite";
  else if (pdc == InDef) return "Indefinite";
  else return "Singular";
}

template <class T> void TestSymDiv(tmv::DivType dt, PosDefCode pdc)
{
  const int N = 10;

  vector<SymMatrixView<T> > s;
  vector<SymMatrixView<complex<T> > > cs;

  Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+0.1*i-0.2*j;
  Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-0.3*i+0.1*j;
  if (pdc == PosDef) {
    a1 /= T(N*N);
    a2 /= T(N*N);
    a1.diag().AddToAll(T(1));
    a2.diag().AddToAll(T(1));
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
  Matrix<complex<T> > ca1 = a1 * complex<T>(3,-4);
  ca1.diag().Imag().Zero();
  Matrix<complex<T> > ca2 = a2 * complex<T>(3,-4);
  ca2.diag().Imag().Zero();

  Vector<T> v1(N);
  Vector<T> v2(N-1);
  for (int i=0; i<N; ++i) v1(i) = 16.-3*i; 
  for (int i=0; i<N-1; ++i) v2(i) = -7.+2*i; 

  HermMatrix<T,Upper,RowMajor> H1(a1);
  HermMatrix<complex<T>,Upper,RowMajor> CH1(ca1);
  s.push_back(H1.View());
  cs.push_back(CH1.View());
  HermMatrix<T,Upper,ColMajor> H2(a1);
  HermMatrix<complex<T>,Upper,ColMajor> CH2(ca1);
  s.push_back(H2.View());
  cs.push_back(CH2.View());
  HermMatrix<T,Lower,RowMajor> H3(a1);
  HermMatrix<complex<T>,Lower,RowMajor> CH3(ca1);
  s.push_back(H3.View());
  cs.push_back(CH3.View());
  HermMatrix<T,Lower,ColMajor> H4(a1);
  HermMatrix<complex<T>,Lower,ColMajor> CH4(ca1);
  s.push_back(H4.View());
  cs.push_back(CH4.View());

  // These need to be outside of (if dt!=CH loop), otherwise they will
  // go out of scope at the end of it and the Views will refer to
  // defunct memory.
  SymMatrix<T,Upper,RowMajor> S1(a1);
  SymMatrix<complex<T>,Upper,RowMajor> CS1(ca1);
  SymMatrix<T,Upper,ColMajor> S2(a1);
  SymMatrix<complex<T>,Upper,ColMajor> CS2(ca1);
  SymMatrix<T,Lower,RowMajor> S3(a1);
  SymMatrix<complex<T>,Lower,RowMajor> CS3(ca1);
  SymMatrix<T,Lower,ColMajor> S4(a1);
  SymMatrix<complex<T>,Lower,ColMajor> CS4(ca1);
  if (dt != tmv::CH) {
    s.push_back(S1.View());
    cs.push_back(CS1.View());
    s.push_back(S2.View());
    cs.push_back(CS2.View());
    s.push_back(S3.View());
    cs.push_back(CS3.View());
    s.push_back(S4.View());
    cs.push_back(CS4.View());
  }

#ifdef DOALLARITH
  HermMatrix<T,Upper> H5(a2);
  HermMatrix<complex<T>,Upper> CH5(ca2);
  s.push_back(H5.SubSymMatrix(0,2*N,2));
  cs.push_back(CH5.SubSymMatrix(0,2*N,2));
  HermMatrix<T,Lower> H6(a2);
  HermMatrix<complex<T>,Lower> CH6(ca2);
  s.push_back(H6.SubSymMatrix(0,2*N,2));
  cs.push_back(CH6.SubSymMatrix(0,2*N,2));

  SymMatrix<T,Upper> S5(a2);
  SymMatrix<complex<T>,Upper> CS5(ca2);
  SymMatrix<T,Lower> S6(a2);
  SymMatrix<complex<T>,Lower> CS6(ca2);
  if (dt != tmv::CH) {
    s.push_back(S5.SubSymMatrix(0,2*N,2));
    cs.push_back(CS5.SubSymMatrix(0,2*N,2));
    s.push_back(S6.SubSymMatrix(0,2*N,2));
    cs.push_back(CS6.SubSymMatrix(0,2*N,2));
  }
#endif

  size_t ntot = s.size();

  Matrix<T> a3 = a1.Cols(0,N/2);
  Matrix<complex<T> > ca3 = ca1.Cols(0,N/2);
  Matrix<T> a4 = a1.Rows(0,N/2);
  Matrix<complex<T> > ca4 = ca1.Rows(0,N/2);
  Matrix<T> a5 = a1.Cols(0,0);
  Matrix<complex<T> > ca5 = ca1.Cols(0,0);
  Matrix<T> a6 = a1.Rows(0,0);
  Matrix<complex<T> > ca6 = ca1.Rows(0,0);
  HermMatrix<T> s1 = H1;
  HermMatrix<complex<T> > cs1 = s1;
  cs1.LowerTri().OffDiag() *= complex<T>(2,-1);

#ifdef SHOWCHECK
  ostream* checkout = &cout;
#else
  ostream* checkout = 0;
#endif

  for(size_t i=STARTAT;i<ntot;i++) {
    SymMatrixView<T> si = s[i];
    SymMatrixView<complex<T> > csi = cs[i];
    si.SaveDiv();
    csi.SaveDiv();
#ifdef SHOWACC
    cout<<"Start loop: i = "<<i<<", si = "<<si<<endl;
#endif

    Matrix<T> m = si;
    m.SaveDiv();
    if (dt == tmv::CH) m.DivideUsing(tmv::LU);
    else m.DivideUsing(dt);
    m.SetDiv();
    Assert(m.CheckDecomp(checkout),"CheckDecomp m"); 
    T eps = EPS*Norm(m)*Norm(m.Inverse());
    si.DivideUsing(dt);
    si.SetDiv();
    if (si.isherm()) {
      HermMatrix<T> six = si;
      six.DivideUsing(dt);
      six.SetDiv();
      Assert(six.CheckDecomp(checkout),"CheckDecomp si(herm)"); 
    } else {
      SymMatrix<T> six = si;
      six.DivideUsing(dt);
      six.SetDiv();
      Assert(six.CheckDecomp(checkout),"CheckDecomp si(sym)"); 
    }

    Vector<T> x1 = v1/si;
    Vector<T> x2 = v1/m;
#ifdef SHOWACC
    //cerr<<"v1 = "<<v1<<endl;
    //cerr<<"v1/si = "<<x1<<endl;
    //cerr<<"v1/m = "<<x2<<endl;
    //cerr<<"si*x1 = "<<si*x1<<endl;
    //cerr<<"m*x2 = "<<m*x2<<endl;
    cout<<"v/b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(v1)<<endl;
#endif
    Assert(Norm(x1-x2) < eps*Norm(v1),"Sym v/b");

    x1 = v1%si;
    x2 = v1%m;
#ifdef SHOWACC
    cout<<"v%b: Norm(x1-x2) = "<<Norm(x1-x2)<<"  "<<eps*Norm(v1)<<endl;
#endif
    Assert(Norm(x1-x2) < eps*Norm(v1),"Sym v%b");

    Matrix<T,ColMajor> sinv = si.Inverse();
    Matrix<T,ColMajor> minv = m.Inverse();
#ifdef SHOWACC
    //cerr<<"sinv = "<<sinv<<endl;
    //cerr<<"minv = "<<minv<<endl;
    cout<<"Norm(minv-sinv) = "<<Norm(minv-sinv)<<"  "<<eps*Norm(sinv)<<endl;
#endif
    Assert(Norm(sinv-minv) < eps*Norm(sinv),"Sym Inverse");

    if (pdc != Sing) {
#ifdef SHOWACC
      //cout<<"si.Det = "<<si.Det()<<", m.Det = "<<m.Det()<<endl;
      cout<<"abs(sdet-mdet) = "<<abs(si.Det()-m.Det());
      cout<<"  EPS*abs(mdet) = "<<eps*abs(m.Det())<<endl;
      //cout<<"abs(abs(sdet)-abs(mdet)) = "<<abs(abs(si.Det())-abs(m.Det()));
      //cout<<"  EPS*abs(mdet) = "<<eps*abs(m.Det())<<endl;
#endif
      Assert(abs(m.Det()-si.Det()) < eps*abs(m.Det()),"Sym Det");
    }

    Matrix<complex<T> > cm(csi);
    cm.SaveDiv();
    if (dt == tmv::CH) cm.DivideUsing(tmv::LU);
    else cm.DivideUsing(dt);
    cm.SetDiv();
    Assert(cm.CheckDecomp(checkout),"CheckDecomp cm"); 
    csi.DivideUsing(dt);
    csi.SetDiv();
    if (csi.isherm()) {
      HermMatrix<complex<T> > csix = csi;
      csix.DivideUsing(dt);
      csix.SetDiv();
      Assert(csix.CheckDecomp(checkout),"CheckDecomp csi(herm)"); 
    } else {
      SymMatrix<complex<T> > csix = csi;
      csix.DivideUsing(dt);
      csix.SetDiv();
      Assert(csix.CheckDecomp(checkout),"CheckDecomp csi(sym)"); 
    }

    T ceps = EPS*Norm(cm)*Norm(cm.Inverse());

    if (pdc != Sing) {
#ifdef SHOWACC
      //cout<<"csi.Det = "<<csi.Det()<<", cm.Det = "<<cm.Det()<<endl;
      cout<<"abs(csidet-cmdet) = "<<abs(csi.Det()-cm.Det());
      cout<<"  csidet/cmdet = "<<csi.Det()/cm.Det();
      cout<<"  EPS*abs(cmdet) = "<<ceps*abs(cm.Det())<<endl;
      //cout<<"abs(abs(csdet)-abs(cmdet)) = "<<abs(abs(csi.Det())-abs(cm.Det()));
      //cout<<"  EPS*abs(cmdet) = "<<ceps*abs(cm.Det())<<endl;
#endif
      Assert(abs(csi.Det()-cm.Det()) < ceps*abs(cm.Det()),"Sym CDet");
    }

    Vector<complex<T> > cv(v1 * complex<T>(1,1));
    cv(1) += complex<T>(-1,5);
    cv(2) -= complex<T>(-1,5);

    // test real / complex
    Vector<complex<T> > y1 = v1/csi;
    Vector<complex<T> > y2 = v1/cm;
#ifdef SHOWACC
    cout<<"v/cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(v1)<<endl;
#endif
    Assert(Norm(y1-y2) < ceps*Norm(v1),"Sym v/cs");

    // test complex / real
    y1 = cv/si;
    y2 = cv/m;
#ifdef SHOWACC
    cout<<"cv/b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(cv)<<endl;
#endif
    Assert(Norm(y1-y2) < eps*Norm(cv),"Sym cv/b");

    // test complex / complex
    y1 = cv/csi;
    y2 = cv/cm;
#ifdef SHOWACC
    cout<<"cv/cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(cv)<<endl;
#endif
    Assert(Norm(y1-y2) < ceps*Norm(cv),"Sym cv/cs");

    y1 = v1%csi;
    y2 = v1%cm;
#ifdef SHOWACC
    cout<<"v%cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(v1)<<endl;
#endif
    Assert(Norm(y1-y2) < ceps*Norm(v1),"Sym v%cs");

    y1 = cv%si;
    y2 = cv%m;
#ifdef SHOWACC
    cout<<"cv%b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<eps*Norm(cv)<<endl;
#endif
    Assert(Norm(y1-y2) < eps*Norm(cv),"Sym cv%b");
    y1 = cv%csi;
    y2 = cv%cm;
#ifdef SHOWACC
    cout<<"cv%cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<ceps*Norm(cv)<<endl;
#endif
    Assert(Norm(y1-y2) < ceps*Norm(cv),"Sym cv%cs");


    if (pdc == PosDef)
      TestMatrixDivArith<T>(dt==tmv::CH?tmv::LU:dt,a1,si,ca1,csi,
	  "Sym/SquareMatrix");
    if (pdc != Sing) {
      TestMatrixDivArith<T>(dt,si,a1,csi,ca1,"SquareMatrix/Sym");
      TestMatrixDivArith<T>(dt,si,a3,csi,ca3,"NonSquareMatrix/Sym");
      TestMatrixDivArith<T>(dt,si,a4,csi,ca4,"NonSquareMatrix/Sym");
      TestMatrixDivArith<T>(dt,si,a5,csi,ca5,"DegenerateMatrix/Sym");
      TestMatrixDivArith<T>(dt,si,a6,csi,ca6,"DegenerateMatrix/Sym");
      TestMatrixDivArith<T>(dt,si,s1,csi,cs1,"Sym/Sym");
    }
  }

  cout<<PDLabel(pdc)<<" SymMatrix<"<<tmv::Type(T())<<"> Division using ";
  cout<<tmv::Text(dt)<<" passed all tests\n";
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
