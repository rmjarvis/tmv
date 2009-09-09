//#define SHOWTESTS
//#define DOALLARITH

#include "TMV_Test.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_TestMatrixArith.h"
#include "TMV_TestOProd.h"

using tmv::Matrix;
using tmv::SymMatrix;
using tmv::HermMatrix;
using tmv::SymMatrixView;
using tmv::MatrixView;
using tmv::RowMajor;
using tmv::ColMajor;
using tmv::Upper;
using tmv::Lower;

template <class T1, class T2> bool CanAddEq(
    const SymMatrixView<T1>& a, const MatrixView<T2>& b)
{ return false; }

template <class T1, class T2> bool CanMultEq(
    const SymMatrixView<T1>& a, const MatrixView<T2>& b)
{ return false; }

template <class T1, class T2> bool CanMultEq2(
    const MatrixView<T1>& a, const SymMatrixView<T2>& b)
{ return false; }

template <class T1, class T2> bool CanAddEqX(
    const SymMatrixView<T1>& a, const T2 x)
{ return tmv::IsReal(x) || !a.isherm(); }

template <class T1, class T2> bool CanMultEqX(
    const SymMatrixView<T1>& a, const T2 x)
{ return tmv::IsReal(x) || !a.isherm(); }

template <class T> void TestSymMatrixArith_B()
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

  Matrix<T,RowMajor> a3 = a2.Rows(0,N);
  Matrix<complex<T> > ca3 = a3 * complex<T>(-3,4.);
  Matrix<T,RowMajor> a4 = a1.Cols(0,0);
  Matrix<complex<T> > ca4 = a4;

  Vector<T> v = a1.col(0);
  Matrix<T> vv = v^v;
  Matrix<T> rvrv = v.Reverse()^v.Reverse();
  Vector<complex<T> > cv = ca1.col(0);
  Matrix<complex<T> > cvcv = cv^cv;
  Matrix<complex<T> > cvcvt = cv^cv.Conjugate();
  Matrix<complex<T> > crvcrv = cv.Reverse()^cv.Reverse();
  Matrix<complex<T> > crvcrvt = cv.Reverse()^cv.Reverse().Conjugate();

  Matrix<T,ColMajor> w = a1.Cols(0,3);
  Matrix<T> ww = w*w.Transpose();
  Matrix<complex<T>,ColMajor> cw = ca1.Cols(0,3);
  Matrix<complex<T> > cwcw = cw*cw.Transpose();
  Matrix<complex<T> > cwcwt = cw*cw.Adjoint();

  for(size_t i=0;i<s.size();i++) {
    if (cs[i].isherm()) {
      TestMatrixArith<T,HermMatrix<T>,HermMatrix<complex<T> > >(
	  s[i],a1.View(),cs[i].View(),ca1.View(),"Sym/SquareM");
      TestMatrixArith<T,HermMatrix<T>,HermMatrix<complex<T> > >(
	  s[i],a3.View(),cs[i].View(),ca3.View(),"Sym/NonSquareM");
      TestMatrixArith<T,HermMatrix<T>,HermMatrix<complex<T> > >(
	  s[i],a4.View(),cs[i].View(),ca4.View(),"Sym/DegenerateM");
    } else {
      TestMatrixArith<T,SymMatrix<T>,SymMatrix<complex<T> > >(
	  s[i],a1.View(),cs[i].View(),ca1.View(),"Sym/SquareM");
      TestMatrixArith<T,SymMatrix<T>,SymMatrix<complex<T> > >(
	  s[i],a3.View(),cs[i].View(),ca3.View(),"Sym/NonSquareM");
      TestMatrixArith<T,SymMatrix<T>,SymMatrix<complex<T> > >(
	  s[i],a4.View(),cs[i].View(),ca4.View(),"Sym/DegenerateM");
    }
#ifdef DOALLARITH
    TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
	a1.View(),s[i],ca1.View(),cs[i].View(),"SquareM/Sym");
    TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
	a3.View(),s[i],ca3.View(),cs[i].View(),"NonSquareM/Sym");
    TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
	a4.View(),s[i],ca4.View(),cs[i].View(),"DegenerateM/Sym");
#endif

    Matrix<T> s0 = s[i];
    TestOProd(s[i],v.View(),v.View(),s0,vv,"Sym RR");
    TestOProd(s[i],a1.col(0),a1.col(0).View(),s0,vv,"Sym Step RR");
    TestOProd(s[i],v.Reverse(),v.Reverse(),s0,rvrv,"Sym Rev RR");
    TestOProd(s[i],w.View(),w.View(),s0,ww,"Sym CM RR");
    TestOProd(s[i],a1.Cols(0,3),a1.Cols(0,3),s0,ww,"Sym RM RR");

    Matrix<complex<T> > cs0 = cs[i];
    TestOProd(cs[i],v.View(),v.View(),cs0,vv,"Sym CR");
    TestOProd(cs[i],a1.col(0),a1.col(0),cs0,vv,"Sym Step CR");
    TestOProd(cs[i],v.Reverse(),v.Reverse(),cs0,rvrv,"Sym Rev CR");
    TestOProd(cs[i],w.View(),w.View(),cs0,ww,"Sym CM CR");
    TestOProd(cs[i],a1.Cols(0,3),a1.Cols(0,3),cs0,ww,"Sym RM CR");

    if (cs[i].isherm()) {
      TestOProd(cs[i],cv.View(),cv.Conjugate(),cs0,cvcvt,"Herm CC");
      TestOProd(cs[i],ca1.col(0),ca1.col(0).Conjugate(),cs0,cvcvt,
	  "Herm Step CC");
      TestOProd(cs[i],cv.Reverse(),cv.Conjugate().Reverse(),cs0,crvcrvt,
	  "Herm Rev CC");
      TestOProd(cs[i],cw.View(),cw.Conjugate(),cs0,cwcwt,"Herm CM CC");
      TestOProd(cs[i],ca1.Cols(0,3),ca1.Cols(0,3).Conjugate(),cs0,cwcwt,
	  "Herm RM CC");
    } else {
      TestOProd(cs[i],cv.View(),cv.View(),cs0,cvcv,"Sym CC");
      TestOProd(cs[i],ca1.col(0),ca1.col(0),cs0,cvcv,"Sym Step CC");
      TestOProd(cs[i],cv.Reverse(),cv.Reverse(),cs0,crvcrv,"Sym Rev CC");
      TestOProd(cs[i],cw.View(),cw.View(),cs0,cwcw,"Sym CM CC");
      TestOProd(cs[i],ca1.Cols(0,3),ca1.Cols(0,3),cs0,cwcw,"Sym RM CC");
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
