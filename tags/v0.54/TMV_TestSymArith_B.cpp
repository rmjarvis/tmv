#define START 0

#include "TMV_Test.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_TestMatrixArith.h"
#include "TMV_TestOProd.h"

template <class T1, class T2> inline bool CanAddEq(
    const tmv::SymMatrixView<T1>& , const tmv::MatrixView<T2>& )
{ return false; }

template <class T1, class T2> inline bool CanMultEq(
    const tmv::SymMatrixView<T1>& , const tmv::MatrixView<T2>& )
{ return false; }

template <class T1, class T2> inline bool CanMultEq2(
    const tmv::MatrixView<T1>& , const tmv::SymMatrixView<T2>& )
{ return false; }

template <class T1, class T2> inline bool CanAddEqX(
    const tmv::SymMatrixView<T1>& a, const T2 x)
{ return tmv::IsReal(x) || !a.isherm(); }

template <class T1, class T2> inline bool CanMultEqX(
    const tmv::SymMatrixView<T1>& a, const T2 x)
{ return tmv::IsReal(x) || !a.isherm(); }

template <class T> void TestSymMatrixArith_B()
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

  tmv::Matrix<T,tmv::RowMajor> a3 = a2.Rows(0,N);
  tmv::Matrix<std::complex<T> > ca3 = a3 * std::complex<T>(-3,4.);
  tmv::Matrix<T,tmv::RowMajor> a4 = a1.Cols(0,0);
  tmv::Matrix<std::complex<T> > ca4 = a4;

  tmv::Vector<T> v = a1.col(0);
  tmv::Matrix<T> vv = v^v;
  tmv::Matrix<T> rvrv = v.Reverse()^v.Reverse();
  tmv::Vector<std::complex<T> > cv = ca1.col(0);
  tmv::Matrix<std::complex<T> > cvcv = cv^cv;
  tmv::Matrix<std::complex<T> > cvcvt = cv^cv.Conjugate();
  tmv::Matrix<std::complex<T> > crvcrv = cv.Reverse()^cv.Reverse();
  tmv::Matrix<std::complex<T> > crvcrvt = cv.Reverse()^cv.Reverse().Conjugate();

  for(size_t i=START;i<s.size();i++) {
    if (cs[i].isherm()) {
      TestMatrixArith<T,tmv::HermMatrix<T>,tmv::HermMatrix<std::complex<T> > >(
	  s[i],a1.View(),cs[i].View(),ca1.View(),"Sym/SquareM");
#ifdef XTEST
      TestMatrixArith<T,tmv::HermMatrix<T>,tmv::HermMatrix<std::complex<T> > >(
	  s[i],a3.View(),cs[i].View(),ca3.View(),"Sym/NonSquareM");
      TestMatrixArith<T,tmv::HermMatrix<T>,tmv::HermMatrix<std::complex<T> > >(
	  s[i],a4.View(),cs[i].View(),ca4.View(),"Sym/DegenerateM");
#endif
    } else {
      TestMatrixArith<T,tmv::SymMatrix<T>,tmv::SymMatrix<std::complex<T> > >(
	  s[i],a1.View(),cs[i].View(),ca1.View(),"Sym/SquareM");
#ifdef XTEST
      TestMatrixArith<T,tmv::SymMatrix<T>,tmv::SymMatrix<std::complex<T> > >(
	  s[i],a3.View(),cs[i].View(),ca3.View(),"Sym/NonSquareM");
      TestMatrixArith<T,tmv::SymMatrix<T>,tmv::SymMatrix<std::complex<T> > >(
	  s[i],a4.View(),cs[i].View(),ca4.View(),"Sym/DegenerateM");
#endif
    }
#ifdef XTESt
    if (doallarith) {
      TestMatrixArith<T,Matrix<T>,Matrix<std::complex<T> > >(
	  a1.View(),s[i],ca1.View(),cs[i].View(),"SquareM/Sym");
      TestMatrixArith<T,Matrix<T>,Matrix<std::complex<T> > >(
	  a3.View(),s[i],ca3.View(),cs[i].View(),"NonSquareM/Sym");
      TestMatrixArith<T,Matrix<T>,Matrix<std::complex<T> > >(
	  a4.View(),s[i],ca4.View(),cs[i].View(),"DegenerateM/Sym");
    }
#endif

    tmv::Matrix<T> s0 = s[i];
    TestOProd(s[i],v.View(),v.View(),s0,vv,"Sym RR");
    TestOProd(s[i],a1.col(0),a1.col(0).View(),s0,vv,"Sym Step RR");
    TestOProd(s[i],v.Reverse(),v.Reverse(),s0,rvrv,"Sym Rev RR");

    tmv::Matrix<std::complex<T> > cs0 = cs[i];
    TestOProd(cs[i],v.View(),v.View(),cs0,vv,"Sym CR");
    TestOProd(cs[i],a1.col(0),a1.col(0),cs0,vv,"Sym Step CR");
    TestOProd(cs[i],v.Reverse(),v.Reverse(),cs0,rvrv,"Sym Rev CR");

    if (cs[i].isherm()) {
      TestOProd(cs[i],cv.View(),cv.Conjugate(),cs0,cvcvt,"Herm CC");
      TestOProd(cs[i],ca1.col(0),ca1.col(0).Conjugate(),cs0,cvcvt,
	  "Herm Step CC");
      TestOProd(cs[i],cv.Reverse(),cv.Conjugate().Reverse(),cs0,crvcrvt,
	  "Herm Rev CC");
    } else {
      TestOProd(cs[i],cv.View(),cv.View(),cs0,cvcv,"Sym CC");
      TestOProd(cs[i],ca1.col(0),ca1.col(0),cs0,cvcv,"Sym Step CC");
      TestOProd(cs[i],cv.Reverse(),cv.Reverse(),cs0,crvcrv,"Sym Rev CC");
    }
  }
}

template void TestSymMatrixArith_B<double>();
#ifndef NOFLOAT
template void TestSymMatrixArith_B<float>();
#endif
#ifdef LONGDOUBLE
template void TestSymMatrixArith_B<long double>();
#endif
