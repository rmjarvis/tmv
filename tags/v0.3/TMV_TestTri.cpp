//#define SHOWTESTS

#include "TMV_Test.h"
#include "TMV_Tri.h"
#include "TMV.h"

using tmv::Matrix;
using tmv::UpperTriMatrix;
using tmv::LowerTriMatrix;
using tmv::UpperTriMatrixView;
using tmv::LowerTriMatrixView;
using tmv::ConstUpperTriMatrixView;
using tmv::ConstLowerTriMatrixView;
using tmv::NonUnitDiag;
using tmv::UnitDiag;

template <class T> extern void TestTriMatrixArith_A();
template <class T> extern void TestTriMatrixArith_B();

template <class T> void TestTriMatrix()
{
  const int N = 10;

  Matrix<T> a(N,N);
  for(int i=0;i<N;i++) for(int j=0;j<N;j++) a(i,j) = 12+3*i-5*j;

  UpperTriMatrix<T,NonUnitDiag> u1(N);
  UpperTriMatrix<T,UnitDiag> u2(N);
  LowerTriMatrix<T,NonUnitDiag> l1(N);
  LowerTriMatrix<T,UnitDiag> l2(N);

  Assert(u1.colsize() == size_t(N) && u1.rowsize() == size_t(N),
      "Creating UpperTriMatrix(N)");
  Assert(u2.colsize() == size_t(N) && u2.rowsize() == size_t(N),
      "Creating UpperTriMatrix(N)");
  Assert(l1.colsize() == size_t(N) && l1.rowsize() == size_t(N),
      "Creating LowerTriMatrix(N)");
  Assert(l2.colsize() == size_t(N) && l2.rowsize() == size_t(N),
      "Creating LowerTriMatrix(N)");

  u1 = UpperTriMatrixViewOf(a,NonUnitDiag);
  u2 = UpperTriMatrixViewOf(a,UnitDiag);
  l1 = LowerTriMatrixViewOf(a,NonUnitDiag);
  l2 = LowerTriMatrixViewOf(a,UnitDiag);

  const UpperTriMatrix<T,NonUnitDiag>& xu2 = u2;
  const LowerTriMatrix<T,NonUnitDiag>& xl2 = l2;

  for(int i=0;i<N;i++) for(int j=0;j<N;j++) {
    if (i > j) {
      Assert(l1(i,j) == a(i,j),"Read/Write TriMatrix");
      Assert(l2(i,j) == a(i,j),"Read/Write TriMatrix");
    } else if (i==j) {
      Assert(l1(i,j) == a(i,j),"Read/Write TriMatrix");
      Assert(u1(i,j) == a(i,j),"Read/Write TriMatrix");
      Assert(xl2(i,j) == 1,"Read/Write TriMatrix");
      Assert(xu2(i,j) == 1,"Read/Write TriMatrix");
    } else {
      Assert(u1(i,j) == a(i,j),"Read/Write TriMatrix");
      Assert(u2(i,j) == a(i,j),"Read/Write TriMatrix");
    }
  }

  Assert(u1==UpperTriMatrixView<T>(a.View(),NonUnitDiag),"TriMatrix ==");
  Assert(u2==UpperTriMatrixView<T>(a.View(),UnitDiag),"TriMatrix ==");
  Assert(l1==LowerTriMatrixView<T>(a.View(),NonUnitDiag),"TriMatrix ==");
  Assert(l2==LowerTriMatrixView<T>(a.View(),UnitDiag),"TriMatrix ==");

  Matrix<T> b1(u1);
  Matrix<T> b2(u2);
  Matrix<T> b3(l1);
  Matrix<T> b4(l2);

  for(int i=0;i<N;i++) for(int j=0;j<N;j++) {
    if (i > j) {
      Assert(b1(i,j) == 0,"b1 TriMatrix -> Matrix");
      Assert(b2(i,j) == 0,"b2 TriMatrix -> Matrix");
      Assert(b3(i,j) == a(i,j),"b3 TriMatrix -> Matrix");
      Assert(b4(i,j) == a(i,j),"b4 TriMatrix -> Matrix");
    } else if (i==j) {
      Assert(b1(i,j) == a(i,j),"b1 TriMatrix -> Matrix");
      Assert(b2(i,j) == 1,"b2 TriMatrix -> Matrix");
      Assert(b3(i,j) == a(i,j),"b3 TriMatrix -> Matrix");
      Assert(b4(i,j) == 1,"b4 TriMatrix -> Matrix");
    } else {
      Assert(b1(i,j) == a(i,j),"b1 TriMatrix -> Matrix");
      Assert(b2(i,j) == a(i,j),"b2 TriMatrix -> Matrix");
      Assert(b3(i,j) == 0,"b3 TriMatrix -> Matrix");
      Assert(b4(i,j) == 0,"b4 TriMatrix -> Matrix");
    }
  }

  UpperTriMatrix<T> u4 = u1+u1;
  UpperTriMatrix<T> u5 = u1+u2;
  UpperTriMatrix<T> u6 = u2+u2;
  LowerTriMatrix<T> l4 = l1+l1;
  LowerTriMatrix<T> l5 = l1+l2;
  LowerTriMatrix<T> l6 = l2+l2;
  Matrix<T> m1 = l1+u1;
  Matrix<T> m2 = l1+u2;
  Matrix<T> m3 = l2+u2;

  for(int i=0;i<N;i++) for(int j=0;j<N;j++) {
    if (i > j) {
      Assert(l4(i,j) == 2*a(i,j),"Add TriMatrices l4");
      Assert(l5(i,j) == 2*a(i,j),"Add TriMatrices l5");
      Assert(l6(i,j) == 2*a(i,j),"Add TriMatrices l6");
      Assert(m1(i,j) == a(i,j),"Add TriMatrices m1");
      Assert(m2(i,j) == a(i,j),"Add TriMatrices m2");
      Assert(m3(i,j) == a(i,j),"Add TriMatrices m3");
    } else if (i==j) {
      Assert(l4(i,j) == 2*a(i,j),"Add TriMatrices l4");
      Assert(l5(i,j) == 1+a(i,j),"Add TriMatrices l5");
      Assert(l6(i,j) == 2,"Add TriMatrices l6");
      Assert(u4(i,j) == 2*a(i,j),"Add TriMatrices u4");
      Assert(u5(i,j) == 1+a(i,j),"Add TriMatrices u5");
      Assert(u6(i,j) == 2,"Add TriMatrices u6");
      Assert(m1(i,j) == 2*a(i,j),"Add TriMatrices m1");
      Assert(m2(i,j) == 1+a(i,j),"Add TriMatrices m2");
      Assert(m3(i,j) == 2,"Add TriMatrices m3");
    } else {
      Assert(u4(i,j) == 2*a(i,j),"Add TriMatrices u4");
      Assert(u5(i,j) == 2*a(i,j),"Add TriMatrices u5");
      Assert(u6(i,j) == 2*a(i,j),"Add TriMatrices u6");
      Assert(m1(i,j) == a(i,j),"Add TriMatrices m1");
      Assert(m2(i,j) == a(i,j),"Add TriMatrices m2");
      Assert(m3(i,j) == a(i,j),"Add TriMatrices m3");
    }
  }

  u4 = u1-u1;
  u5 = u1-u2;
  u6 = u2-u2;
  l4 = l1-l1;
  l5 = l1-l2;
  l6 = l2-l2;
  m1 = l1-u1;
  m2 = l1-u2;
  m3 = l2-u2;

  for(int i=0;i<N;i++) for(int j=0;j<N;j++) {
    if (i > j) {
      Assert(l4(i,j) == 0,"Subtract TriMatrices l4");
      Assert(l5(i,j) == 0,"Subtract TriMatrices l5");
      Assert(l6(i,j) == 0,"Subtract TriMatrices l6");
      Assert(m1(i,j) == a(i,j),"Subtract TriMatrices m1");
      Assert(m2(i,j) == a(i,j),"Subtract TriMatrices m2");
      Assert(m3(i,j) == a(i,j),"Subtract TriMatrices m3");
    } else if (i==j) {
      Assert(l4(i,j) == 0,"Subtract TriMatrices l4");
      Assert(l5(i,j) == a(i,j)-1,"Subtract TriMatrices l5");
      Assert(l6(i,j) == 0,"Subtract TriMatrices l6");
      Assert(u4(i,j) == 0,"Subtract TriMatrices u4");
      Assert(u5(i,j) == a(i,j)-1,"Subtract TriMatrices u5");
      Assert(u6(i,j) == 0,"Subtract TriMatrices u6");
      Assert(m1(i,j) == 0,"Subtract TriMatrices m1");
      Assert(m2(i,j) == a(i,j)-1,"Subtract TriMatrices m2");
      Assert(m3(i,j) == 0,"Subtract TriMatrices m3");
    } else {
      Assert(u4(i,j) == 0,"Subtract TriMatrices u4");
      Assert(u5(i,j) == 0,"Subtract TriMatrices u5");
      Assert(u6(i,j) == 0,"Subtract TriMatrices u6");
      Assert(m1(i,j) == -a(i,j),"Subtract TriMatrices m1");
      Assert(m2(i,j) == -a(i,j),"Subtract TriMatrices m2");
      Assert(m3(i,j) == -a(i,j),"Subtract TriMatrices m3");
    }
  }

  cout<<"TriMatrix<"<<tmv::Type(T())<<"> passed all tests\n";

  TestTriMatrixArith_A<T>();
  TestTriMatrixArith_B<T>();
}

template void TestTriMatrix<double>();
#ifndef NOFLOAT
template void TestTriMatrix<float>();
#endif
#ifdef LONGDOUBLE
template void TestTriMatrix<long double>();
#endif
