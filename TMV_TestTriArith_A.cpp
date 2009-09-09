#include "TMV_Test.h"
#include "TMV_Tri.h"
#include "TMV.h"

using tmv::Matrix;
using tmv::UpperTriMatrix;
using tmv::LowerTriMatrix;
using tmv::UpperTriMatrixView;
using tmv::LowerTriMatrixView;
using tmv::NonUnitDiag;
using tmv::UnitDiag;
using tmv::RowMajor;
using tmv::IndexStyle;

template <class T1, class T2, IndexStyle I1, IndexStyle I2> bool CanAddEq(
    const UpperTriMatrixView<T1,I1>& a, const UpperTriMatrixView<T2,I2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2, IndexStyle I1, IndexStyle I2> bool CanAddEq(
    const LowerTriMatrixView<T1,I1>& a, const LowerTriMatrixView<T2,I2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2, IndexStyle I1, IndexStyle I2> bool CanAddEq(
    const UpperTriMatrixView<T1,I1>& a, const LowerTriMatrixView<T2,I2>& b)
{ return false; }

template <class T1, class T2, IndexStyle I1, IndexStyle I2> bool CanAddEq(
    const LowerTriMatrixView<T1,I1>& a, const UpperTriMatrixView<T2,I2>& b)
{ return false; }

template <class T1, class T2, IndexStyle I1> bool CanAddEqX(
    const UpperTriMatrixView<T1,I1>& a, const T2 x)
{ return !a.isunit(); }

template <class T1, class T2, IndexStyle I1> bool CanAddEqX(
    const LowerTriMatrixView<T1,I1>& a, const T2 x)
{ return !a.isunit(); }

template <class T1, class T2, IndexStyle I1> bool CanMultEqX(
    const UpperTriMatrixView<T1,I1>& a, const T2 x)
{ return !a.isunit(); }

template <class T1, class T2, IndexStyle I1> bool CanMultEqX(
    const LowerTriMatrixView<T1,I1>& a, const T2 x)
{ return !a.isunit(); }

template <class T1, class T2, IndexStyle I1, IndexStyle I2> bool CanMultEq(
    const UpperTriMatrixView<T1,I1>& a, const UpperTriMatrixView<T2,I2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2, IndexStyle I1, IndexStyle I2> bool CanMultEq(
    const LowerTriMatrixView<T1,I1>& a, const LowerTriMatrixView<T2,I2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2, IndexStyle I1, IndexStyle I2> bool CanMultEq(
    const UpperTriMatrixView<T1,I1>& a, const LowerTriMatrixView<T2,I2>& b)
{ return false; }

template <class T1, class T2, IndexStyle I1, IndexStyle I2> bool CanMultEq(
    const LowerTriMatrixView<T1,I1>& a, const UpperTriMatrixView<T2,I2>& b)
{ return false; }

template <class T1, class T2, IndexStyle I1, IndexStyle I2> bool CanMultEq2(
    const UpperTriMatrixView<T1,I1>& a, const UpperTriMatrixView<T2,I2>& b)
{ return a.size() == b.size() && !b.isunit(); }

template <class T1, class T2, IndexStyle I1, IndexStyle I2> bool CanMultEq2(
    const LowerTriMatrixView<T1,I1>& a, const LowerTriMatrixView<T2,I2>& b)
{ return a.size() == b.size() && !b.isunit(); }

template <class T1, class T2, IndexStyle I1, IndexStyle I2> bool CanMultEq2(
    const UpperTriMatrixView<T1,I1>& a, const LowerTriMatrixView<T2,I2>& b)
{ return false; }

template <class T1, class T2, IndexStyle I1, IndexStyle I2> bool CanMultEq2(
    const LowerTriMatrixView<T1,I1>& a, const UpperTriMatrixView<T2,I2>& b)
{ return false; }

template <class T, IndexStyle I> bool CanDoSV(const UpperTriMatrixView<T,I>& a)
{ return false; }
template <class T, IndexStyle I> bool CanDoSV(const LowerTriMatrixView<T,I>& a)
{ return false; }

#include "TMV_TestMatrixArith.h"

template <class T> void TestTriMatrixArith_A()
{
  const int N = 10;

  Matrix<T> a1(N,N);
  for(int i=0;i<N;i++) for(int j=0;j<N;j++) a1(i,j) = 12+3*i-5*j;

  UpperTriMatrix<T,NonUnitDiag> u1(a1);
  UpperTriMatrix<T,UnitDiag> u2(a1);
  LowerTriMatrix<T,NonUnitDiag> l1(a1);
  LowerTriMatrix<T,UnitDiag> l2(a1);

  Matrix<complex<T> > c1 = a1 * complex<T>(1,2);
  UpperTriMatrix<complex<T>,NonUnitDiag> cu1(c1);
  UpperTriMatrix<complex<T>,UnitDiag> cu2(c1);
  LowerTriMatrix<complex<T>,NonUnitDiag> cl1(c1);
  LowerTriMatrix<complex<T>,UnitDiag> cl2(c1);

  TestMatrixArith<T,LowerTriMatrix<T>,LowerTriMatrix<complex<T> > >(
	l1.View(),u1.View(),cl1.View(),cu1.View(),"Tri/Tri");
#ifdef XTEST
  TestMatrixArith<T,LowerTriMatrix<T>,LowerTriMatrix<complex<T> > >(
	l1.View(),u2.View(),cl1.View(),cu2.View(),"Tri/Tri");
  TestMatrixArith<T,LowerTriMatrix<T,UnitDiag>,
    LowerTriMatrix<complex<T>,UnitDiag> >(
	l2.View(),u1.View(),cl2.View(),cu1.View(),"Tri/Tri");
  TestMatrixArith<T,LowerTriMatrix<T,UnitDiag>,
    LowerTriMatrix<complex<T>,UnitDiag> >(
	l2.View(),u2.View(),cl2.View(),cu2.View(),"Tri/Tri");
  if (doallarith) {
    TestMatrixArith<T,LowerTriMatrix<T>,LowerTriMatrix<complex<T> > >(
	l1.View(),l2.View(),cl1.View(),cl2.View(),"Tri/Tri");
    TestMatrixArith<T,LowerTriMatrix<T,UnitDiag>,LowerTriMatrix<complex<T>,UnitDiag> >(
	  l2.View(),l1.View(),cl2.View(),cl1.View(),"Tri/Tri");
    TestMatrixArith<T,UpperTriMatrix<T>,UpperTriMatrix<complex<T> > >(
	u1.View(),u2.View(),cu1.View(),cu2.View(),"Tri/Tri");
    TestMatrixArith<T,UpperTriMatrix<T,UnitDiag>,UpperTriMatrix<complex<T>,UnitDiag> >(
	  u2.View(),u1.View(),cu2.View(),cu1.View(),"Tri/Tri");
    TestMatrixArith<T,UpperTriMatrix<T>,UpperTriMatrix<complex<T> > >(
	u1.View(),l1.View(),cu1.View(),cl1.View(),"Tri/Tri");
    TestMatrixArith<T,UpperTriMatrix<T>,UpperTriMatrix<complex<T> > >(
	u1.View(),l2.View(),cu1.View(),cl2.View(),"Tri/Tri");
    TestMatrixArith<T,UpperTriMatrix<T,UnitDiag>,UpperTriMatrix<complex<T>,UnitDiag> >(
	  u2.View(),l1.View(),cu2.View(),cl1.View(),"Tri/Tri");
    TestMatrixArith<T,UpperTriMatrix<T,UnitDiag>,UpperTriMatrix<complex<T>,UnitDiag> >(
	  u2.View(),l2.View(),cu2.View(),cl2.View(),"Tri/Tri");
  }
#endif

  cout<<"TriMatrix<"<<tmv::Type(T())<<"> (Tri/Tri) Arithmetic passed all tests\n";
}

template void TestTriMatrixArith_A<double>();
#ifndef NOFLOAT
template void TestTriMatrixArith_A<float>();
#endif
#ifdef LONGDOUBLE
template void TestTriMatrixArith_A<long double>();
#endif

