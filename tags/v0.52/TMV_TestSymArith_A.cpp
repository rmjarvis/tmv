//#define SHOWTESTS
//#define DOALLARITH

#include "TMV_Test.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_TestMatrixArith.h"

using tmv::Matrix;
using tmv::SymMatrix;
using tmv::HermMatrix;
using tmv::SymMatrixView;
using tmv::ConjItType;
using tmv::RowMajor;
using tmv::ColMajor;
using tmv::Upper;
using tmv::Lower;

template <class T1, class T2> bool CanAdd(
    const SymMatrixView<T1>& a, const SymMatrixView<T2>& b)
{ return a.size() == b.size(); }

template <class T> bool CanAdd(
    const SymMatrixView<complex<T> >& a, const SymMatrixView<complex<T> >& b)
{ return a.size() == b.size() && a.sym() == b.sym(); }

template <class T1, class T2> bool CanMultEq(
    const SymMatrixView<T1>& a, const SymMatrixView<T2>& b)
{ return false; }

template <class T1, class T2> bool CanMultEq2(
    const SymMatrixView<T1>& a, const SymMatrixView<T2>& b)
{ return false; }

template <class T1, class T2> bool CanAddX(
    const SymMatrixView<T1>& a, const T2 x)
{ return tmv::IsReal(x) || !a.isherm(); }

template <class T1, class T2> bool CanMultX(
    const SymMatrixView<T1>& a, const T2 x)
{ return tmv::IsReal(x) || !a.isherm(); }

template <class T> void TestSymMatrixArith_A()
{
  const int N = 10;

  vector<SymMatrixView<T> > s;
  vector<SymMatrixView<complex<T> > > cs;

  Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-5*j;
  Matrix<complex<T> > ca1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) = 
    complex<T>(3.+i-5*j,2.-3.*i);
  Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+6*j;
  Matrix<complex<T> > ca2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) = 
    complex<T>(1.-3.*i+6*j,-4.+2.*j);

  SymMatrix<T,Upper,RowMajor> S1(a1);
  HermMatrix<T,Upper,RowMajor> H1(a1);
  SymMatrix<complex<T>,Upper,RowMajor> CS1(ca1);
  HermMatrix<complex<T>,Upper,RowMajor> CH1(ca1);
  CH1.diag().Imag().Zero();
  s.push_back(S1.View());
  s.push_back(H1.View());
  cs.push_back(CS1.View());
  cs.push_back(CH1.View());
  SymMatrix<T,Upper,ColMajor> S2(a1);
  HermMatrix<T,Upper,ColMajor> H2(a1);
  SymMatrix<complex<T>,Upper,ColMajor> CS2(ca1);
  HermMatrix<complex<T>,Upper,ColMajor> CH2(ca1);
  CH2.diag().Imag().Zero();
  s.push_back(S2.View());
  s.push_back(H2.View());
  cs.push_back(CS2.View());
  cs.push_back(CH2.View());
#ifdef DOALLARITH
  SymMatrix<T,Lower,RowMajor> S3(a1);
  HermMatrix<T,Lower,RowMajor> H3(a1);
  SymMatrix<complex<T>,Lower,RowMajor> CS3(ca1);
  HermMatrix<complex<T>,Lower,RowMajor> CH3(ca1);
  CH3.diag().Imag().Zero();
  s.push_back(S3.View());
  s.push_back(H3.View());
  cs.push_back(CS3.View());
  cs.push_back(CH3.View());
  SymMatrix<T,Lower,ColMajor> S4(a1);
  HermMatrix<T,Lower,ColMajor> H4(a1);
  SymMatrix<complex<T>,Lower,ColMajor> CS4(ca1);
  HermMatrix<complex<T>,Lower,ColMajor> CH4(ca1);
  CH4.diag().Imag().Zero();
  s.push_back(S4.View());
  s.push_back(H4.View());
  cs.push_back(CS4.View());
  cs.push_back(CH4.View());
  Matrix<T> a2a = a2;
  Matrix<T> a2b = a2;
  Matrix<complex<T> > ca2a = ca2;
  Matrix<complex<T> > ca2b = ca2;
  ca2b.diag().Imag().Zero();
  s.push_back(SymMatrixViewOf(a2a,Upper).SubSymMatrix(0,2*N,2));
  s.push_back(HermMatrixViewOf(a2b,Upper).SubSymMatrix(0,2*N,2));
  cs.push_back(SymMatrixViewOf(ca2a,Upper).SubSymMatrix(0,2*N,2));
  cs.push_back(HermMatrixViewOf(ca2b,Upper).SubSymMatrix(0,2*N,2));
#endif

  for(size_t i=0;i<s.size();i++) for(size_t j=0;j<s.size();j++) if (i!=j) {
    //cerr<<"i,j = "<<i<<','<<j<<endl;
    //cerr<<"s[i] = "<<s[i]<<endl;
    //cerr<<"s[j] = "<<s[j]<<endl;
    if (cs[i].isherm())
      TestMatrixArith<T,HermMatrix<T>,HermMatrix<complex<T> > >(
	  s[i],s[j],cs[i],cs[j],"Sym/Sym");
    else
      TestMatrixArith<T,SymMatrix<T>,SymMatrix<complex<T> > >(
	  s[i],s[j],cs[i],cs[j],"Sym/Sym");
  }
}

template void TestSymMatrixArith_A<double>();
#ifndef NOFLOAT
template void TestSymMatrixArith_A<float>();
#endif
#ifdef LONGDOUBLE
template void TestSymMatrixArith_A<long double>();
#endif
