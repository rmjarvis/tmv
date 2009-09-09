//#define SHOWTESTS
//#define DOALLARITH

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

template <class T1, class T2> bool CanAddEq(
    const UpperTriMatrixView<T1>& a, const tmv::MatrixView<T2>& b)
{ return false; }

template <class T1, class T2> bool CanAddEq(
    const LowerTriMatrixView<T1>& a, const tmv::MatrixView<T2>& b)
{ return false; }

template <class T1, class T2> bool CanAddEqX(const UpperTriMatrixView<T1>& a,
    const T2 x)
{ return !a.isunit(); }

template <class T1, class T2> bool CanAddEqX(const LowerTriMatrixView<T1>& a,
    const T2 x)
{ return !a.isunit(); }

template <class T1, class T2> bool CanMultEqX(const UpperTriMatrixView<T1>& a,
    const T2 x)
{ return !a.isunit(); }

template <class T1, class T2> bool CanMultEqX(const LowerTriMatrixView<T1>& a,
    const T2 x)
{ return !a.isunit(); }

template <class T1, class T2> bool CanMultEq(
    const UpperTriMatrixView<T1>& a, const tmv::MatrixView<T2>& b)
{ return false; }

template <class T1, class T2> bool CanMultEq(
    const LowerTriMatrixView<T1>& a, const tmv::MatrixView<T2>& b)
{ return false; }

template <class T1, class T2> bool CanMultEq2(
    const tmv::MatrixView<T1>& a, const UpperTriMatrixView<T2>& b)
{ return false; }

template <class T1, class T2> bool CanMultEq2(
    const tmv::MatrixView<T1>& a, const LowerTriMatrixView<T2>& b)
{ return false; }

template <class T> bool CanDoSV(const UpperTriMatrixView<T>& m)
{ return false; }
template <class T> bool CanDoSV(const LowerTriMatrixView<T>& m)
{ return false; }

#include "TMV_TestMatrixArith.h"

template <class T> void TestTriMatrixArith_B()
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

  TestMatrixArith<T,UpperTriMatrix<T>,UpperTriMatrix<complex<T> > >(
	u1.View(),a1.View(),cu1.View(),c1.View(),"Tri/Square");
  TestMatrixArith<T,UpperTriMatrix<T,UnitDiag>,
    UpperTriMatrix<complex<T>,UnitDiag> >(
	u2.View(),a1.View(),cu2.View(),c1.View(),"Tri/Square");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a1.View(),u1.View(),c1.View(),cu1.View(),"Square/Tri");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a1.View(),u2.View(),c1.View(),cu2.View(),"Square/Tri");
#ifdef DOALLARITh
  TestMatrixArith<T,LowerTriMatrix<T>,LowerTriMatrix<complex<T> > >(
      l1.View(),a1.View(),cl1.View(),c1.View(),"Tri/Square");
  TestMatrixArith<T,LowerTriMatrix<T,UnitDiag>,
    LowerTriMatrix<complex<T>,UnitDiag> >(
      l2.View(),a1.View(),cl2.View(),c1.View(),"Tri/Square");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a1.View(),l1.View(),c1.View(),cl1.View(),"Square/Tri");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a1.View(),l2.View(),c1.View(),cl2.View(),"Square/Tri");
#endif

  Matrix<T> a2(15,10);
  for(int i=0;i<15;++i) for(int j=0;j<10;++j) a2(i,j) = 1-3*i+2*j;
  Matrix<complex<T> > c2 = a2*complex<T>(1,2);

  TestMatrixArith<T,UpperTriMatrix<T>,UpperTriMatrix<complex<T> > >(
      u1.View(),a2.View(),cu1.View(),c2.View(),"Tri/NonSquare");
  TestMatrixArith<T,UpperTriMatrix<T,UnitDiag>,
    UpperTriMatrix<complex<T>,UnitDiag> >(
      u2.View(),a2.View(),cu2.View(),c2.View(),"Tri/NonSquare");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a2.View(),u1.View(),c2.View(),cu1.View(),"NonSquare/Tri");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a2.View(),u2.View(),c2.View(),cu2.View(),"NonSquare/Tri");
#ifdef DOALLARITH
  TestMatrixArith<T,LowerTriMatrix<T>,LowerTriMatrix<complex<T> > >(
      l1.View(),a2.View(),cl1.View(),c2.View(),"Tri/NonSquare");
  TestMatrixArith<T,LowerTriMatrix<T,UnitDiag>,
    LowerTriMatrix<complex<T>,UnitDiag> >(
      l2.View(),a2.View(),cl2.View(),c2.View(),"Tri/NonSquare");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a2.View(),l1.View(),c2.View(),cl1.View(),"NonSquare/Tri");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a2.View(),l2.View(),c2.View(),cl2.View(),"NonSquare/Tri");
#endif

  Matrix<T> a3(10,0,1);
  Matrix<complex<T> > c3 = a3;

  TestMatrixArith<T,UpperTriMatrix<T>,UpperTriMatrix<complex<T> > >(
      u1.View(),a3.View(),cu1.View(),c3.View(),"Tri/Degenerate");
  TestMatrixArith<T,UpperTriMatrix<T,UnitDiag>,
    UpperTriMatrix<complex<T>,UnitDiag> >(
      u2.View(),a3.View(),cu2.View(),c3.View(),"Tri/Degenerate");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a3.View(),u1.View(),c3.View(),cu1.View(),"Degenerate/Tri");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a3.View(),u2.View(),c3.View(),cu2.View(),"Degenerate/Tri");
#ifdef DOALLARITH
  TestMatrixArith<T,LowerTriMatrix<T>,LowerTriMatrix<complex<T> > >(
      l1.View(),a3.View(),cl1.View(),c3.View(),"Tri/Degenerate");
  TestMatrixArith<T,LowerTriMatrix<T,UnitDiag>,
    LowerTriMatrix<complex<T>,UnitDiag> >(
      l2.View(),a3.View(),cl2.View(),c3.View(),"Tri/Degenerate");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a3.View(),l1.View(),c3.View(),cl1.View(),"Degenerate/Tri");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a3.View(),l2.View(),c3.View(),cl2.View(),"Degenerate/Tri");
#endif

  cout<<"TriMatrix<"<<tmv::Type(T())<<"> (Matrix/Tri) Arithmetic passed all tests\n";
}

template void TestTriMatrixArith_B<double>();
#ifndef NOFLOAT
template void TestTriMatrixArith_B<float>();
#endif
#ifdef LONGDOUBLE
template void TestTriMatrixArith_B<long double>();
#endif

