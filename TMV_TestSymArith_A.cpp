#define START1 8
#define START2 0

#include "TMV_Test.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_TestMatrixArith.h"

template <class T1, class T2> inline bool CanAdd(
    const tmv::SymMatrixView<T1>& a, const tmv::SymMatrixView<T2>& b)
{ return a.size() == b.size(); }

template <class T> inline bool CanAdd(
    const tmv::SymMatrixView<std::complex<T> >& a,
    const tmv::SymMatrixView<std::complex<T> >& b)
{ return a.size() == b.size() && a.sym() == b.sym(); }

template <class T1, class T2> inline bool CanMultEq(
    const tmv::SymMatrixView<T1>& , const tmv::SymMatrixView<T2>& )
{ return false; }

template <class T1, class T2> inline bool CanMultEq2(
    const tmv::SymMatrixView<T1>& , const tmv::SymMatrixView<T2>& )
{ return false; }

template <class T1, class T2> inline bool CanAddX(
    const tmv::SymMatrixView<T1>& a, const T2 x)
{ return tmv::IsReal(x) || !a.isherm(); }

template <class T1, class T2> inline bool CanMultX(
    const tmv::SymMatrixView<T1>& a, const T2 x)
{ return tmv::IsReal(x) || !a.isherm(); }

template <class T> void TestSymMatrixArith_A()
{
  const int N = 10;

  std::vector<tmv::SymMatrixView<T> > s;
  std::vector<tmv::SymMatrixView<std::complex<T> > > cs;

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-5*j;
  tmv::Matrix<std::complex<T> > ca1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) = 
    std::complex<T>(3.+i-5*j,2.-3.*i);
  tmv::Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+6*j;
  tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) = 
    std::complex<T>(1.-3.*i+6*j,-4.+2.*j);

  tmv::SymMatrix<T,tmv::Upper,tmv::RowMajor> S1(a1);
  tmv::HermMatrix<T,tmv::Upper,tmv::RowMajor> H1(a1);
  tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> CS1(ca1);
  tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> CH1(ca1);
  CH1.diag().Imag().Zero();
  s.push_back(S1.View());
  s.push_back(H1.View());
  cs.push_back(CS1.View());
  cs.push_back(CH1.View());

#ifdef XTEST
  tmv::SymMatrix<T,tmv::Upper,tmv::ColMajor> S2(a1);
  tmv::HermMatrix<T,tmv::Upper,tmv::ColMajor> H2(a1);
  tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> CS2(ca1);
  tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> CH2(ca1);
  CH2.diag().Imag().Zero();
  s.push_back(S2.View());
  s.push_back(H2.View());
  cs.push_back(CS2.View());
  cs.push_back(CH2.View());

  tmv::SymMatrix<T,tmv::Lower,tmv::RowMajor> S3(a1);
  tmv::HermMatrix<T,tmv::Lower,tmv::RowMajor> H3(a1);
  tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> CS3(ca1);
  tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> CH3(ca1);
  CH3.diag().Imag().Zero();
  tmv::SymMatrix<T,tmv::Lower,tmv::ColMajor> S4(a1);
  tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor> H4(a1);
  tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> CS4(ca1);
  tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> CH4(ca1);
  CH4.diag().Imag().Zero();
  tmv::Matrix<T> a2a = a2;
  tmv::Matrix<T> a2b = a2;
  tmv::Matrix<std::complex<T> > ca2a = ca2;
  tmv::Matrix<std::complex<T> > ca2b = ca2;
  ca2b.diag().Imag().Zero();
  if (doallarith) {
    s.push_back(S3.View());
    s.push_back(H3.View());
    cs.push_back(CS3.View());
    cs.push_back(CH3.View());
    s.push_back(S4.View());
    s.push_back(H4.View());
    cs.push_back(CS4.View());
    cs.push_back(CH4.View());
    s.push_back(SymMatrixViewOf(a2a,tmv::Upper).SubSymMatrix(0,2*N,2));
    s.push_back(HermMatrixViewOf(a2b,tmv::Upper).SubSymMatrix(0,2*N,2));
    cs.push_back(SymMatrixViewOf(ca2a,tmv::Upper).SubSymMatrix(0,2*N,2));
    cs.push_back(HermMatrixViewOf(ca2b,tmv::Upper).SubSymMatrix(0,2*N,2));
  }
#endif

  for(size_t i=START1;i<s.size();i++) 
    for(size_t j=START2;j<s.size();j++) if (i!=j) {
      if (showstartdone) {
	std::cout<<"i,j = "<<i<<','<<j<<std::endl;
	std::cout<<"s[i] = "<<s[i]<<std::endl;
	std::cout<<"s[j] = "<<s[j]<<std::endl;
      }
      if (cs[i].isherm())
	TestMatrixArith<T,tmv::HermMatrix<T>,tmv::HermMatrix<std::complex<T> > >(
	    s[i],s[j],cs[i],cs[j],"Sym/Sym");
      else {
	TestMatrixArith<T,tmv::SymMatrix<T>,tmv::SymMatrix<std::complex<T> > >(
	    s[i],s[j],cs[i],cs[j],"Sym/Sym");
      }
    }
}

template void TestSymMatrixArith_A<double>();
#ifndef NOFLOAT
template void TestSymMatrixArith_A<float>();
#endif
#ifdef LONGDOUBLE
template void TestSymMatrixArith_A<long double>();
#endif
