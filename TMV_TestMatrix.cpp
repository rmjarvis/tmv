//#define SHOWTESTS
//#define DOALLARITH

#include "TMV_Test.h"
#include "TMV.h"
#include "TMV_TestMatrixArith.h"
#include "TMV_TestOProd.h"

using tmv::Matrix;
using tmv::Vector;

template <class T, StorageType S> void TestBasicMatrix()
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

  Matrix<complex<T>,S> cm(M,N);
  Matrix<complex<T>,S> ca(M,N);
  Matrix<complex<T>,S> cb(M,N);
  Assert(cm.colsize() == size_t(M) && cm.rowsize() && size_t(N),
      "Creating CMatrix(M,N)");

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    cm(i,j) = complex<T>(k,k+1000);

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    Assert(cm(i,j) == complex<T>(k,k+1000),"Read/Write CMatrix");

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    Assert(cm.row(i)(j) == complex<T>(k,k+1000),"CMatrix.row");

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    Assert(cm.col(j)(i) == complex<T>(k,k+1000),"CMatrix.col");

  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) {
    ca(i,j) = complex<T>(3.+i+5*j,0.+i-j);
    cb(i,j) = complex<T>(3.+2*i+4*j,4.-10*i);
  }

  cm = ca+cb;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(cm(i,j) == complex<T>(6.+3*i+9*j,4.-9*i-j),"Add CMatrix");

  cm = ca-cb;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(cm(i,j) == complex<T>(0.-i+j,-4.+11*i-j),"Subtract CMatrix");

  cm = ca;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(cm(i,j) == ca(i,j),"Copy CMatrix");

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

  Matrix<T,tmv::ColMajor> a2 = a1.Transpose();
  a2.row(1) *= T(3);
  a2.col(2) -= Vector<T>(4,4.);
  Matrix<complex<T>,tmv::ColMajor> c2 = c1;
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

  Vector<T> v1 = a1.col(0);
  Vector<T> v2 = a1.row(2);
  Matrix<T> v1v2 = v1 ^ v2;
  Matrix<T> a0 = a1;
  TestOProd(a1.View(),v1,v2,a0,v1v2,"Square RRR");
  Vector<complex<T> > cv1 = c1.col(1);
  Vector<complex<T> > cv2 = c1.row(3);
  Matrix<complex<T> > cv1cv2 = cv1 ^ cv2;
  Matrix<complex<T> > cv1v2 = cv1 ^ v2;
  Matrix<complex<T> > v1cv2 = v1 ^ cv2;
  Matrix<complex<T> > c0 = c1;
  TestOProd(c1.View(),cv1,cv2,c0,cv1cv2,"Square CCC");
  TestOProd(c1.View(),v1,cv2,c0,v1cv2,"Square CRC");
  TestOProd(c1.View(),cv1,v2,c0,cv1v2,"Square CCR");
  TestOProd(c1.View(),v1,v2,c0,v1v2,"Square CRR");

  Matrix<T,tmv::RowMajor> a4(7,4);
  for(int i=0;i<7;++i) for(int j=0;j<4;++j) a4(i,j) = 1-3*i+2*j;
  Matrix<T,tmv::ColMajor> a5 = a4.Transpose();
  a4.SubMatrix(2,6,0,4) += a1;
  a5.SubMatrix(0,4,1,5) -= a2;

  Matrix<complex<T>,tmv::RowMajor> c4 = a4*complex<T>(1,2);
  Matrix<complex<T>,tmv::ColMajor> c5 = c4.Adjoint();
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

  Vector<T> v3 = a4.col(2);
  Matrix<T> v3v2 = v3 ^ v2;
  Matrix<T> v2v3 = v2 ^ v3;
  Matrix<T> a40 = a4;
  Matrix<T> a50 = a5;
  TestOProd(a4.View(),v3,v2,a40,v3v2,"NonSqaure RRR");
  TestOProd(a5.View(),v2,v3,a50,v2v3,"NonSqaure RRR");
  Vector<complex<T> > cv3 = c4.col(3);
  Matrix<complex<T> > cv3cv2 = cv3 ^ cv2;
  Matrix<complex<T> > cv2cv3 = cv2 ^ cv3;
  Matrix<complex<T> > cv3v2 = cv3 ^ v2;
  Matrix<complex<T> > cv2v3 = cv2 ^ v3;
  Matrix<complex<T> > v3cv2 = v3 ^ cv2;
  Matrix<complex<T> > v2cv3 = v2 ^ cv3;
  Matrix<complex<T> > c40 = c4;
  Matrix<complex<T> > c50 = c5;
  TestOProd(c4.View(),cv3,cv2,c40,cv3cv2,"NonSqaure CCC");
  TestOProd(c5.View(),cv2,cv3,c50,cv2cv3,"NonSqaure CCC");
  TestOProd(c4.View(),v3,cv2,c40,v3cv2,"NonSqaure CRC");
  TestOProd(c5.View(),v2,cv3,c50,v2cv3,"NonSqaure CRC");
  TestOProd(c4.View(),cv3,v2,c40,cv3v2,"NonSqaure CCR");
  TestOProd(c5.View(),cv2,v3,c50,cv2v3,"NonSqaure CCR");
  TestOProd(c4.View(),v3,v2,c40,v3v2,"NonSqaure CRR");
  TestOProd(c5.View(),v2,v3,c50,v2v3,"NonSqaure CRR");

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
      a7.View(),a6.View(),c7.View(),c6.View(),"Degenerate"); 
#endif

  Vector<T> v4 = a6.row(2);
  Matrix<T> v4v2 = v4 ^ v2;
  Matrix<T> v2v4 = v2 ^ v4;
  Matrix<T> a60 = a6;
  Matrix<T> a70 = a7;
  TestOProd(a6.View(),v2,v4,a60,v2v4,"Degenerate RRR");
  TestOProd(a7.View(),v4,v2,a70,v4v2,"Degenerate RRR");
  Vector<complex<T> > cv4 = c6.row(2);
  Matrix<complex<T> > cv4cv2 = cv4 ^ cv2;
  Matrix<complex<T> > cv2cv4 = cv2 ^ cv4;
  Matrix<complex<T> > cv4v2 = cv4 ^ v2;
  Matrix<complex<T> > cv2v4 = cv2 ^ v4;
  Matrix<complex<T> > v4cv2 = v4 ^ cv2;
  Matrix<complex<T> > v2cv4 = v2 ^ cv4;
  Matrix<complex<T> > c60 = c6;
  Matrix<complex<T> > c70 = c7;
  TestOProd(c6.View(),cv2,cv4,c60,cv2cv4,"Degenerate CCC");
  TestOProd(c7.View(),cv4,cv2,c70,cv4cv2,"Degenerate CCC");
  TestOProd(c6.View(),v2,cv4,c60,v2cv4,"Degenerate CRC");
  TestOProd(c7.View(),v4,cv2,c70,v4cv2,"Degenerate CRC");
  TestOProd(c6.View(),cv2,v4,c60,cv2v4,"Degenerate CCR");
  TestOProd(c7.View(),cv4,v2,c70,cv4v2,"Degenerate CCR");
  TestOProd(c6.View(),v2,v4,c60,v2v4,"Degenerate CRR");
  TestOProd(c7.View(),v4,v2,c70,v4v2,"Degenerate CRR");

}

template <class T> void TestAllMatrix()
{
  TestBasicMatrix<T,tmv::RowMajor>();
  TestBasicMatrix<T,tmv::ColMajor>();
  cout<<"Matrix"<<tmv::Type(T())<<"> passed all basic tests\n";

  TestAllMatrixArith<T>();
  cout<<"Matrix<"<<tmv::Type(T())<<"> Arithmetic passed all tests\n";
}

template void TestAllMatrix<double>();
#ifndef NOFLOAT
template void TestAllMatrix<float>();
#endif
#ifdef LONGDOUBLE
template void TestAllMatrix<long double>();
#endif
