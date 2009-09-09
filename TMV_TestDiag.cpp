
#include "TMV.h"
#include "TMV_Diag.h"
#include "TMV_Test.h"
#include "TMV_TestMatrixArith.h"

using tmv::Matrix;
using tmv::DiagMatrix;
using tmv::DiagMatrixView;
using tmv::MatrixView;
using tmv::RowMajor;
using tmv::IndexStyle;

template <class T1, IndexStyle I1, class T2, IndexStyle I2> bool CanAddEq(
    const DiagMatrixView<T1,I1>& a, const MatrixView<T2,I2>& b)
{ return false; }

template <class T1, IndexStyle I1, class T2, IndexStyle I2> bool CanMultEq(
    const DiagMatrixView<T1,I1>& a, const MatrixView<T2,I2>& b)
{ return false; }

template <class T1, IndexStyle I1, class T2, IndexStyle I2> bool CanMultEq2(
    const MatrixView<T2,I2>& b, const DiagMatrixView<T1,I1>& a)
{ return false; }

template <class T, IndexStyle I> bool CanDoSV(const DiagMatrixView<T,I>& a)
{ return false; }

template <class T> void TestDiagMatrix()
{
  const int N = 10;

  DiagMatrix<T> a(N);
  DiagMatrix<T,tmv::FortranStyle> af(N);
  Assert(a.colsize() == size_t(N) && a.rowsize() == size_t(N),
      "Creating DiagMatrix(N)");
  Assert(af.colsize() == size_t(N) && af.rowsize() == size_t(N),
      "Creating DiagMatrix(N)");

  for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
    if (i == j) a(i,j) = k;
    if (i == j) af(i+1,j+1) = k;
  }
  tmv::ConstDiagMatrixView<T> acv = a.View();
  tmv::DiagMatrixView<T> av = a.View();
  tmv::ConstDiagMatrixView<T,tmv::FortranStyle> afcv = af.View();
  tmv::DiagMatrixView<T,tmv::FortranStyle> afv = af.View();

  for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k)
    if (i == j) {
      Assert(a(i,j) == k,"Read/Write DiagMatrix");
      Assert(acv(i,j) == k,"Access DiagMatrix CV");
      Assert(av(i,j) == k,"Access DiagMatrix V");
      Assert(af(i+1,j+1) == k,"Read/Write DiagMatrixF");
      Assert(afcv(i+1,i+1) == k,"Access DiagMatrixF CV");
      Assert(afv(i+1,i+1) == k,"Access DiagMatrixF V");
      Assert(a(i) == k,"Single argument access for DiagMatrix");
      Assert(acv(i) == k,"Single argument access for DiagMatrix CV");
      Assert(av(i) == k,"Single argument access for DiagMatrix V");
      Assert(af(i+1) == k,"Single argument access for DiagMatrixF");
      Assert(afcv(i+1) == k,"Single argument access for DiagMatrixF CV");
      Assert(afv(i+1) == k,"Single argument access for DiagMatrixF V");
    }

  Assert(a==af,"CStyle Matrix == FortranStyle Matrix");
  Assert(a==acv,"Matrix == ConstMatrixView");
  Assert(a==av,"Matrix == MatrixView");
  Assert(a==afcv,"Matrix == FortranStyle ConstMatrixView");
  Assert(a==afv,"Matrix == FortranStyle MatrixView");

  DiagMatrix<T> b(N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
    if (i == j) {
      a(i,j) = 3.+i+5*j;
      b(i,j) = 5.+2*i+4*j;
    }
  af = a;
  Assert(a==af,"Copy CStyle DiagMatrix to FotranStyle");

  DiagMatrix<T> c(N);
  c = a+b;
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
    if (i == j)
      Assert(c(i,j) == 8.+3*i+9*j,"Add DiagMatrices");

  c = a-b;
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
    if (i == j)
      Assert(c(i,j) == -2.-i+j,"Subtract DiagMatrices");

  Matrix<T> m(a);
  for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k)
    if (i == j)
      Assert(a(i,j) == m(i,j),"DiagMatrix -> Matrix");
  Assert(a == DiagMatrix<T>(m),"Matrix -> DiagMatrix");

  DiagMatrix<complex<T> > ca = a*complex<T>(1,2);
  DiagMatrix<complex<T> > cb = b*complex<T>(-5,-1);

  Matrix<T> p(N,N,5);
  p.diag().AddToAll(N*10);
  Matrix<complex<T> > cp = p*complex<T>(2,3);
  cp += ca;

  Matrix<T> q(2*N,N,-2);
  q.Rows(0,N) += p;
  q.Rows(N,2*N) -= p;
  q.Rows(N/2,3*N/2) += T(4)*p;
  Matrix<complex<T> > cq = q*complex<T>(-1,4);
  cq.Rows(0,N) -= cp;
  cq.Rows(N,2*N) -= cp;
  cq.Rows(N/2,3*N/2) -= cp;

  Matrix<T> r(N,0,1);
  Matrix<complex<T> > cr(N,0,1);

  TestMatrixArith<T,DiagMatrix<T>,DiagMatrix<complex<T> > >(
      a.View(),p.View(),ca.View(),cp.View(), "Diag/SquareM");
  TestMatrixArith<T,DiagMatrix<T>,DiagMatrix<complex<T> > >(
      a.View(),q.View(),ca.View(),cq.View(), "Diag/NonSquareM");
#ifdef XTEST
  TestMatrixArith<T,DiagMatrix<T>,DiagMatrix<complex<T> > >(
      a.View(),r.View(),ca.View(),cr.View(), "Diag/DegenM");
  TestMatrixArith<T,DiagMatrix<T>,DiagMatrix<complex<T> > >(
      a.View(),b.View(),ca.View(),cb.View(),"Diag/Diag");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      p.View(),a.View(),cp.View(),ca.View(), "SquareM/Diag");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      q.View(),a.View(),cq.View(),ca.View(), "NonSquareM/Diag");
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      r.View(),a.View(),cr.View(),ca.View(), "DegenM/Diag");
#endif

  cout<<"DiagMatrix<"<<tmv::Type(T())<<"> passed all tests\n";
}

template void TestDiagMatrix<double>();
#ifndef NOFLOAT
template void TestDiagMatrix<float>();
#endif
#ifdef LONGDOUBLE
template void TestDiagMatrix<long double>();
#endif
