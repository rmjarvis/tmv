//#define SHOWTESTS
//#define DOALLARITH

#include "TMV_Test.h"
#include "TMV.h"
#include "TMV_TestMatrixArith.h"

using tmv::Matrix;
using tmv::Vector;

template <class T, StorageType S> void TestMatrixReal()
{
  const int M = 15;
  const int N = 10;

  Matrix<T,S> m(M,N);
  Assert(m.colsize() == size_t(M) && m.rowsize() == size_t(N),
      "Creating Matrix(M,N)");

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    m(i,j) = k;

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    Assert(m(i,j) == k,"Read/Write Matrix");

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    Assert(m.row(i)(j) == k,"Matrix.row");

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    Assert(m.col(j)(i) == k,"Matrix.col");

  Matrix<T,S> a(M,N);
  Matrix<T,S> b(M,N);
  Matrix<T,S> c(M,N);
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) {
    a(i,j) = 3.+i+5*j;
    b(i,j) = 5.+2*i+4*j;
  }

  c = a+b;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(c(i,j) == 8.+3*i+9*j,"Add Matrices");

  c = a-b;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(c(i,j) == -2.-i+j,"Subtract Matrices");

  cout<<"MatrixReal<"<<tmv::Type(T())<<","<<tmv::Text(S)<<"> passed all tests\n";
}

template <class T, StorageType S> void TestMatrixComplex()
{
  const int M = 15;
  const int N = 10;

  Matrix<complex<T>,S> m(M,N);
  Assert(m.colsize() == size_t(M) && m.rowsize() && size_t(N),
      "Creating CMatrix(M,N)");

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    m(i,j) = complex<T>(k,k+1000);

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    Assert(m(i,j) == complex<T>(k,k+1000),"Read/Write CMatrix");

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    Assert(m.row(i)(j) == complex<T>(k,k+1000),"CMatrix.row");

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    Assert(m.col(j)(i) == complex<T>(k,k+1000),"CMatrix.col");

  Matrix<complex<T>,S> a(M,N);
  Matrix<complex<T>,S> b(M,N);
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) {
    a(i,j) = complex<T>(3.+i+5*j,0.+i-j);
    b(i,j) = complex<T>(3.+2*i+4*j,4.-10*i);
  }

  m = a+b;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(m(i,j) == complex<T>(6.+3*i+9*j,4.-9*i-j),"Add CMatrix");

  m = a-b;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(m(i,j) == complex<T>(0.-i+j,-4.+11*i-j),"Subtract CMatrix");

  m = a;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(m(i,j) == a(i,j),"Copy CMatrix");

  cout<<"MatrixComplex<"<<tmv::Type(T())<<","<<tmv::Text(S)<<"> passed all tests\n";
}

template <class T> void TestAllMatrixArith()
{
  Matrix<T> a1(4,4);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) a1(i,j) = 2+4*i-5*j;
  a1(0,0) = 14.;
  a1(1,0) = -2.;
  a1(2,0) = 7.;
  a1(3,0) = -10.;
  a1(2,2) = 30.;

  Matrix<complex<T> > c1 = a1;
  c1(2,3) += complex<T>(2,3);
  c1(1,0) *= complex<T>(0,2);
  c1.col(1) *= complex<T>(-1,3);
  c1.row(3) += Vector<complex<T> >(4,complex<T>(1,9));

  Matrix<T,ColMajor> a2 = a1.Transpose();
  a2.row(1) *= T(3);
  a2.col(2) -= Vector<T>(4,4.);
  Matrix<complex<T>,ColMajor> c2 = c1;
  c2 -= a2;
  c2 *= complex<T>(1,-2);

  Matrix<T> a3x(12,16);
  for(int i=0;i<12;++i) for(int j=0;j<16;++j) a3x(i,j) = 1-2*i+3*j;
  a3x.diag().AddToAll(30);
  Matrix<complex<T> > c3x = a3x*complex<T>(1,-2);
  c3x.diag().AddToAll(complex<T>(-22,15));
  tmv::MatrixView<T> a3 = a3x.SubMatrix(0,12,0,16,3,4);
  tmv::MatrixView<complex<T> > c3 = c3x.SubMatrix(0,12,0,16,3,4);

  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a1.View(),a2.View(),c1.View(),c2.View(),"Square"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a2.View(),a1.View(),c2.View(),c1.View(),"Square"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a3.View(),a1.View(),c3.View(),c1.View(),"Square"); 
#ifdef DOALLARITH
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a3.View(),a2.View(),c3.View(),c2.View(),"Square"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a1.View(),a3.View(),c1.View(),c3.View(),"Square"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a2.View(),a3.View(),c2.View(),c3.View(),"Square"); 
#endif

  cout<<"SquareMatrix<"<<tmv::Type(T())<<"> Arithmetic passed all tests\n";

  Matrix<T,RowMajor> a4(7,4);
  for(int i=0;i<7;++i) for(int j=0;j<4;++j) a4(i,j) = 1-3*i+2*j;
  Matrix<T,ColMajor> a5 = a4.Transpose();
  a4.SubMatrix(2,6,0,4) += a1;
  a5.SubMatrix(0,4,1,5) -= a2;

  Matrix<complex<T>,RowMajor> c4 = a4*complex<T>(1,2);
  Matrix<complex<T>,ColMajor> c5 = c4.Adjoint();
  c4.SubMatrix(2,6,0,4) += c1;
  c5.SubMatrix(0,4,1,5) -= c2;
  c4.col(1) *= complex<T>(2,1);
  c4.row(6).AddToAll(complex<T>(-7,2));
  c5.col(3) *= complex<T>(-1,3);
  c5.row(0).AddToAll(complex<T>(1,9));

  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a1.View(),a4.View(),c1.View(),c4.View(),"NonSquare"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a3.View(),a4.View(),c3.View(),c4.View(),"NonSquare"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a4.View(),a5.View(),c4.View(),c5.View(),"NonSquare"); 
#ifdef DOALLARITH
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a1.View(),a5.View(),c1.View(),c5.View(),"NonSquare"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a2.View(),a4.View(),c2.View(),c4.View(),"NonSquare"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a2.View(),a5.View(),c2.View(),c5.View(),"NonSquare"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a3.View(),a5.View(),c3.View(),c5.View(),"NonSquare"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a4.View(),a1.View(),c4.View(),c1.View(),"NonSquare"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a4.View(),a2.View(),c4.View(),c2.View(),"NonSquare"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a4.View(),a3.View(),c4.View(),c3.View(),"NonSquare"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a5.View(),a1.View(),c5.View(),c1.View(),"NonSquare"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a5.View(),a2.View(),c5.View(),c2.View(),"NonSquare"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a5.View(),a3.View(),c5.View(),c3.View(),"NonSquare"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a5.View(),a4.View(),c5.View(),c4.View(),"NonSquare"); 
#endif

  cout<<"NonSquare Matrix<"<<tmv::Type(T())<<"> Arithmetic passed all tests\n";

  Matrix<T> a6(4,0,1);
  Matrix<T> a7(0,4,1);
  Matrix<complex<T> > c6 = a6;
  Matrix<complex<T> > c7 = a7;

  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a1.View(),a6.View(),c1.View(),c6.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a3.View(),a6.View(),c3.View(),c6.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a4.View(),a6.View(),c4.View(),c6.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a6.View(),a7.View(),c6.View(),c7.View(),"Degenerate"); 
#ifdef DOALLARITH
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a1.View(),a7.View(),c1.View(),c7.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a2.View(),a6.View(),c2.View(),c6.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a2.View(),a7.View(),c2.View(),c7.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a3.View(),a7.View(),c3.View(),c7.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a4.View(),a7.View(),c4.View(),c7.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a5.View(),a6.View(),c5.View(),c6.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a5.View(),a7.View(),c5.View(),c7.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a6.View(),a1.View(),c6.View(),c1.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a6.View(),a2.View(),c6.View(),c2.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a6.View(),a3.View(),c6.View(),c3.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a6.View(),a4.View(),c6.View(),c4.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a6.View(),a5.View(),c6.View(),c5.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a7.View(),a1.View(),c7.View(),c1.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a7.View(),a2.View(),c7.View(),c2.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a7.View(),a3.View(),c7.View(),c3.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a7.View(),a4.View(),c7.View(),c4.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a7.View(),a5.View(),c7.View(),c5.View(),"Degenerate"); 
  TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
      a7.View(),a6.View().View(),c7.View(),c6.View(),"Degenerate"); 
#endif

  cout<<"Degenerate Matrix<"<<tmv::Type(T())<<"> Arithmetic passed all tests\n";
}

template <class T> void TestAllMatrix()
{
  TestMatrixReal<T,RowMajor>();
  TestMatrixReal<T,ColMajor>();
  TestMatrixComplex<T,RowMajor>();
  TestMatrixComplex<T,ColMajor>();
  TestAllMatrixArith<T>();
}

template void TestAllMatrix<double>();
#ifndef NOFLOAT
template void TestAllMatrix<float>();
#endif
#ifdef LONGDOUBLE
template void TestAllMatrix<long double>();
#endif
