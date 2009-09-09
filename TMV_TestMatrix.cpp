#include "TMV_Test.h"
#include "TMV.h"
#include "TMV_TestMatrixArith.h"
#include "TMV_TestOProd.h"
#include <fstream>

template <class T, tmv::StorageType S> inline void TestBasicMatrix()
{
  const int M = 15;
  const int N = 10;

  tmv::Matrix<T,S> m(M,N);
  tmv::Matrix<T,S,tmv::FortranStyle> mf(M,N);
  Assert(m.colsize() == size_t(M) && m.rowsize() == size_t(N),
      "Creating Matrix(M,N)");
  Assert(m.colsize() == size_t(M) && m.rowsize() == size_t(N),
      "Creating MatrixF(M,N)");

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
    m(i,j) = k;
    mf(i+1,j+1) = k;
  }
  tmv::ConstMatrixView<T> mcv = m.View();
  tmv::MatrixView<T> mv = m.View();
  tmv::ConstMatrixView<T,tmv::FortranStyle> mfcv = mf.View();
  tmv::MatrixView<T,tmv::FortranStyle> mfv = mf.View();

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
    Assert(m(i,j) == k,"Read/Write Matrix");
    Assert(mcv(i,j) == k,"Access Matrix CV");
    Assert(mv(i,j) == k,"Access Matrix V");
    Assert(mf(i+1,j+1) == k,"Read/Write MatrixF");
    Assert(mfcv(i+1,j+1) == k,"Access MatrixF CV");
    Assert(mfv(i+1,j+1) == k,"Access MatrixF V");
    Assert(m[i][j] == k,"[] style access of Matrix");
    Assert(mcv[i][j] == k,"[] style access of Matrix CV");
    Assert(mv[i][j] == k,"[] style access of Matrix V");
    Assert(mf[i+1][j+1] == k,"[] style access of MatrixF");
    Assert(mfcv[i+1][j+1] == k,"[] style access of MatrixF CV");
    Assert(mfv[i+1][j+1] == k,"[] style access of MatrixF V");
    Assert(m.row(i)(j) == k,"Matrix.row");
    Assert(mcv.row(i)(j) == k,"Matrix.row CV");
    Assert(mv.row(i)(j) == k,"Matrix.row V");
    Assert(mf.row(i+1)(j+1) == k,"MatrixF.row");
    Assert(mfcv.row(i+1)(j+1) == k,"MatrixF.row CV");
    Assert(mfv.row(i+1)(j+1) == k,"MatrixF.row V");
    Assert(m.row(i,j,N)(0) == k,"Matrix.row2");
    Assert(mcv.row(i,j,N)(0) == k,"Matrix.row2 CV");
    Assert(mv.row(i,j,N)(0) == k,"Matrix.row2 V");
    Assert(mf.row(i+1,j+1,N)(1) == k,"MatrixF.row2");
    Assert(mfcv.row(i+1,j+1,N)(1) == k,"MatrixF.row2 CV");
    Assert(mfv.row(i+1,j+1,N)(1) == k,"MatrixF.row2 V");
    Assert(m.col(j)(i) == k,"Matrix.col");
    Assert(mcv.col(j)(i) == k,"Matrix.col CV");
    Assert(mv.col(j)(i) == k,"Matrix.col V");
    Assert(mf.col(j+1)(i+1) == k,"MatrixF.col");
    Assert(mfcv.col(j+1)(i+1) == k,"MatrixF.col CV");
    Assert(mfv.col(j+1)(i+1) == k,"MatrixF.col V");
    Assert(m.col(j,i,M)(0) == k,"Matrix.col2");
    Assert(mcv.col(j,i,M)(0) == k,"Matrix.col2 CV");
    Assert(mv.col(j,i,M)(0) == k,"Matrix.col2 V");
    Assert(mf.col(j+1,i+1,M)(1) == k,"MatrixF.col2");
    Assert(mfcv.col(j+1,i+1,M)(1) == k,"MatrixF.col2 CV");
    Assert(mfv.col(j+1,i+1,M)(1) == k,"MatrixF.col2 V");
    if (i<j) {
      Assert(m.diag(j-i)(i) == k,"Matrix.diag");
      Assert(mcv.diag(j-i)(i) == k,"Matrix.diag CV");
      Assert(mv.diag(j-i)(i) == k,"Matrix.diag V");
      Assert(mf.diag(j-i)(i+1) == k,"MatrixF.diag");
      Assert(mfcv.diag(j-i)(i+1) == k,"MatrixF.diag CV");
      Assert(mfv.diag(j-i)(i+1) == k,"MatrixF.diag V");
      Assert(m.diag(j-i,i,N-j+i)(0) == k,"Matrix.diag2");
      Assert(mcv.diag(j-i,i,N-j+i)(0) == k,"Matrix.diag2 CV");
      Assert(mv.diag(j-i,i,N-j+i)(0) == k,"Matrix.diag2 V");
      Assert(mf.diag(j-i,i+1,N-j+i)(1) == k,"Matrix.diag2");
      Assert(mfcv.diag(j-i,i+1,N-j+i)(1) == k,"Matrix.diag2 CV");
      Assert(mfv.diag(j-i,i+1,N-j+i)(1) == k,"Matrix.diag2 V");
    } else {
      if (i==j) {
	Assert(m.diag()(i) == k,"Matrix.diag");
	Assert(mcv.diag()(i) == k,"Matrix.diag CV");
	Assert(mv.diag()(i) == k,"Matrix.diag V");
	Assert(mf.diag()(i+1) == k,"MatrixF.diag");
	Assert(mfcv.diag()(i+1) == k,"MatrixF.diag CV");
	Assert(mfv.diag()(i+1) == k,"MatrixF.diag V");
      }
      Assert(m.diag(j-i)(j) == k,"Matrix.diag1");
      Assert(mcv.diag(j-i)(j) == k,"Matrix.diag1 CV");
      Assert(mv.diag(j-i)(j) == k,"Matrix.diag1 V");
      Assert(mf.diag(j-i)(j+1) == k,"MatrixF.diag1");
      Assert(mfcv.diag(j-i)(j+1) == k,"MatrixF.diag1 CV");
      Assert(mfv.diag(j-i)(j+1) == k,"MatrixF.diag1 V");
      if (N+i-j > M) {
	Assert(m.diag(j-i,j,M+j-i)(0) == k,"Matrix.diag2");
	Assert(mcv.diag(j-i,j,M+j-i)(0) == k,"Matrix.diag2 CV");
	Assert(mv.diag(j-i,j,M+j-i)(0) == k,"Matrix.diag2 V");
	Assert(mf.diag(j-i,j+1,M+j-i)(1) == k,"Matrix.diag2");
	Assert(mfcv.diag(j-i,j+1,M+j-i)(1) == k,"Matrix.diag2 CV");
	Assert(mfv.diag(j-i,j+1,M+j-i)(1) == k,"Matrix.diag2 V");
      } else {
	Assert(m.diag(j-i,j,N)(0) == k,"Matrix.diag2");
	Assert(mcv.diag(j-i,j,N)(0) == k,"Matrix.diag2 CV");
	Assert(mv.diag(j-i,j,N)(0) == k,"Matrix.diag2 V");
	Assert(mf.diag(j-i,j+1,N)(1) == k,"Matrix.diag2");
	Assert(mfcv.diag(j-i,j+1,N)(1) == k,"Matrix.diag2 CV");
	Assert(mfv.diag(j-i,j+1,N)(1) == k,"Matrix.diag2 V");
      }
    }
  }
  Assert(m == mf,"CStyle Matrix == FortranStyle Matrix");
  Assert(m == mcv,"Matrix == ConstMatrixView");
  Assert(m == mv,"Matrix == MatrixView");
  Assert(m == mfcv,"Matrix == FortranStyle ConstMatrixView");
  Assert(m == mfv,"Matrix == FortranStyle MatrixView");

  tmv::Matrix<T,S> a(M,N);
  tmv::Matrix<T,S> b(M,N);
  tmv::Matrix<T,S> c(M,N);
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) {
    a(i,j) = 3.+i+5*j;
    b(i,j) = 5.+2*i+4*j;
  }
  mf = a;
  Assert(a == mf,"Copy CStyle Matrix to FortranStyle");

  c = a+b;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(c(i,j) == 8.+3*i+9*j,"Add Matrices");

  c = a-b;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(c(i,j) == -2.-i+j,"Subtract Matrices");

  tmv::Matrix<std::complex<T>,S> cm(M,N);
  tmv::Matrix<std::complex<T>,S> ca(M,N);
  tmv::Matrix<std::complex<T>,S> cb(M,N);
  Assert(cm.colsize() == size_t(M) && cm.rowsize() && size_t(N),
      "Creating CMatrix(M,N)");

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    cm(i,j) = std::complex<T>(k,k+1000);

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
    Assert(cm(i,j) == std::complex<T>(k,k+1000),"Read/Write CMatrix");
    Assert(cm.row(i)(j) == std::complex<T>(k,k+1000),"CMatrix.row");
    Assert(cm.col(j)(i) == std::complex<T>(k,k+1000),"CMatrix.col");
  }

  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) {
    ca(i,j) = std::complex<T>(3.+i+5*j,0.+i-j);
    cb(i,j) = std::complex<T>(3.+2*i+4*j,4.-10*i);
  }

  cm = ca+cb;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(cm(i,j) == std::complex<T>(6.+3*i+9*j,4.-9*i-j),"Add CMatrix");

  cm = ca-cb;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(cm(i,j) == std::complex<T>(0.-i+j,-4.+11*i-j),"Subtract CMatrix");

  cm = ca;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(cm(i,j) == ca(i,j),"Copy CMatrix");

  // Test I/O

  std::ofstream fout("tmvtest_matrix_io.dat");
  if (!fout) throw std::runtime_error(
      "Couldn't open tmvtest_matrix_io.dat for output");
  fout << m << std::endl << cm << std::endl;
  fout.close();

  tmv::Matrix<T,tmv::RowMajor> xm1(M,N);
  tmv::Matrix<std::complex<T>,tmv::RowMajor> xcm1(M,N);
  std::ifstream fin("tmvtest_matrix_io.dat");
  if (!fin) throw std::runtime_error(
      "Couldn't open tmvtest_matrix_io.dat for input");
  fin >> xm1 >> xcm1;
  fin.close();
  Assert(m == xm1,"Matrix I/O check #1");
  Assert(cm == xcm1,"CMatrix I/O check #1");

  tmv::Matrix<T,tmv::ColMajor> xm2(M,N);
  tmv::Matrix<std::complex<T>,tmv::ColMajor> xcm2(M,N);
  fin.open("tmvtest_matrix_io.dat");
  if (!fin) throw std::runtime_error(
      "Couldn't open tmvtest_matrix_io.dat for input");
  fin >> xm2 >> xcm2;
  fin.close();
  Assert(m == xm2,"Matrix I/O check #2");
  Assert(cm == xcm2,"CMatrix I/O check #2");

  std::auto_ptr<tmv::Matrix<T> > xm3;
  std::auto_ptr<tmv::Matrix<std::complex<T> > > xcm3;
  fin.open("tmvtest_matrix_io.dat");
  if (!fin) throw std::runtime_error(
      "Couldn't open tmvtest_matrix_io.dat for input");
  fin >> xm3 >> xcm3;
  fin.close();
  Assert(m == *xm3,"Matrix I/O check #3");
  Assert(cm == *xcm3,"CMatrix I/O check #3");

}

template <class T> inline void TestAllMatrixArith()
{
  tmv::Matrix<T> a1(4,4);
  tmv::Matrix<T,tmv::RowMajor,tmv::FortranStyle> a1f(4,4);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) {
    a1(i,j) = 2+4*i-5*j;
    a1f(i+1,j+1) = 2+4*i-5*j;
  }
  a1(0,0) = 14.; a1f(1,1) = 14.;
  a1(1,0) = -2.; a1f(2,1) = -2.;
  a1(2,0) = 7.; a1f(3,1) = 7.;
  a1(3,0) = -10.; a1f(4,1) = -10.;
  a1(2,2) = 30.; a1f(3,3) = 30.;

  Assert(a1 == a1f, "Define a1f");

  tmv::Matrix<std::complex<T> > c1 = a1;
  c1(2,3) += std::complex<T>(2,3);
  c1(1,0) *= std::complex<T>(0,2);
  c1.col(1) *= std::complex<T>(-1,3);
  c1.row(3) += tmv::Vector<std::complex<T> >(4,std::complex<T>(1,9));
  tmv::Matrix<std::complex<T>,tmv::RowMajor,tmv::FortranStyle> c1f = a1;
  c1f(3,4) += std::complex<T>(2,3);
  c1f(2,1) *= std::complex<T>(0,2);
  c1f.col(2) *= std::complex<T>(-1,3);
  c1f.row(4) += tmv::Vector<std::complex<T> >(4,std::complex<T>(1,9));
  Assert(c1 == c1f, "Define c1f");

  tmv::Matrix<T,tmv::ColMajor> a2 = a1.Transpose();
  a2.row(1) *= T(3);
  a2.col(2) -= tmv::Vector<T>(4,4.);
  tmv::Matrix<std::complex<T>,tmv::ColMajor> c2 = c1;
  c2 -= a2;
  c2 *= std::complex<T>(1,-2);
  tmv::Matrix<T,tmv::ColMajor,tmv::FortranStyle> a2f = a2;
  tmv::Matrix<std::complex<T>,tmv::ColMajor,tmv::FortranStyle> c2f = c2;
  Assert(a2 == a2f, "Define a2f");
  Assert(c2 == c2f, "Define c2f");

  tmv::Matrix<T> a3x(12,16);
  for(int i=0;i<12;++i) for(int j=0;j<16;++j) a3x(i,j) = 1-2*i+3*j;
  a3x.diag().AddToAll(30);
  tmv::Matrix<std::complex<T> > c3x = a3x*std::complex<T>(1,-2);
  c3x.diag().AddToAll(std::complex<T>(-22,15));
  tmv::Matrix<T,tmv::RowMajor,tmv::FortranStyle> a3xf = a3x;
  tmv::MatrixView<T> a3 = a3x.SubMatrix(0,12,0,16,3,4);
  tmv::MatrixView<std::complex<T> > c3 = c3x.SubMatrix(0,12,0,16,3,4);

  TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
      a1.View(),a2.View(),c1.View(),c2.View(),"Square"); 
  TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
      a2.View(),a1.View(),c2.View(),c1.View(),"Square"); 
  TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
      a3.View(),a1.View(),c3.View(),c1.View(),"Square"); 
#ifdef XTEST
  if (doallarith) {
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a3.View(),a2.View(),c3.View(),c2.View(),"Square"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a1.View(),a3.View(),c1.View(),c3.View(),"Square"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a2.View(),a3.View(),c2.View(),c3.View(),"Square"); 
  }
#endif

  tmv::Vector<T> v1 = a1.col(0);
  tmv::Vector<T> v2 = a1.row(2);
  tmv::Matrix<T> v1v2 = v1 ^ v2;
  tmv::Matrix<T> a0 = a1;
  TestOProd(a1.View(),v1,v2,a0,v1v2,"Square RRR");
  tmv::Vector<std::complex<T> > cv1 = c1.col(1);
  tmv::Vector<std::complex<T> > cv2 = c1.row(3);
  tmv::Matrix<std::complex<T> > cv1cv2 = cv1 ^ cv2;
  tmv::Matrix<std::complex<T> > cv1v2 = cv1 ^ v2;
  tmv::Matrix<std::complex<T> > v1cv2 = v1 ^ cv2;
  tmv::Matrix<std::complex<T> > c0 = c1;
  TestOProd(c1.View(),cv1,cv2,c0,cv1cv2,"Square CCC");
  TestOProd(c1.View(),v1,cv2,c0,v1cv2,"Square CRC");
  TestOProd(c1.View(),cv1,v2,c0,cv1v2,"Square CCR");
  TestOProd(c1.View(),v1,v2,c0,v1v2,"Square CRR");

  tmv::Matrix<T,tmv::RowMajor> a4(7,4);
  for(int i=0;i<7;++i) for(int j=0;j<4;++j) a4(i,j) = 1-3*i+2*j;
  tmv::Matrix<T,tmv::ColMajor> a5 = a4.Transpose();
  a4.SubMatrix(2,6,0,4) += a1;
  a5.SubMatrix(0,4,1,5) -= a2;

  tmv::Matrix<std::complex<T>,tmv::RowMajor> c4 = a4*std::complex<T>(1,2);
  tmv::Matrix<std::complex<T>,tmv::ColMajor> c5 = c4.Adjoint();
  c4.SubMatrix(2,6,0,4) += c1;
  c5.SubMatrix(0,4,1,5) -= c2;
  c4.col(1) *= std::complex<T>(2,1);
  c4.row(6).AddToAll(std::complex<T>(-7,2));
  c5.col(3) *= std::complex<T>(-1,3);
  c5.row(0).AddToAll(std::complex<T>(1,9));

  TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
      a1.View(),a4.View(),c1.View(),c4.View(),"NonSquare"); 
  TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
      a3.View(),a4.View(),c3.View(),c4.View(),"NonSquare"); 
  TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
      a4.View(),a5.View(),c4.View(),c5.View(),"NonSquare"); 
#ifdef XTEST
  if (doallarith) {
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a1.View(),a5.View(),c1.View(),c5.View(),"NonSquare"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a2.View(),a4.View(),c2.View(),c4.View(),"NonSquare"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a2.View(),a5.View(),c2.View(),c5.View(),"NonSquare"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a3.View(),a5.View(),c3.View(),c5.View(),"NonSquare"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a4.View(),a1.View(),c4.View(),c1.View(),"NonSquare"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a4.View(),a2.View(),c4.View(),c2.View(),"NonSquare"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a4.View(),a3.View(),c4.View(),c3.View(),"NonSquare"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a5.View(),a1.View(),c5.View(),c1.View(),"NonSquare"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a5.View(),a2.View(),c5.View(),c2.View(),"NonSquare"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a5.View(),a3.View(),c5.View(),c3.View(),"NonSquare"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a5.View(),a4.View(),c5.View(),c4.View(),"NonSquare"); 
  }
#endif

  tmv::Vector<T> v3 = a4.col(2);
  tmv::Matrix<T> v3v2 = v3 ^ v2;
  tmv::Matrix<T> v2v3 = v2 ^ v3;
  tmv::Matrix<T> a40 = a4;
  tmv::Matrix<T> a50 = a5;
  TestOProd(a4.View(),v3,v2,a40,v3v2,"NonSqaure RRR");
  TestOProd(a5.View(),v2,v3,a50,v2v3,"NonSqaure RRR");
  tmv::Vector<std::complex<T> > cv3 = c4.col(3);
  tmv::Matrix<std::complex<T> > cv3cv2 = cv3 ^ cv2;
  tmv::Matrix<std::complex<T> > cv2cv3 = cv2 ^ cv3;
  tmv::Matrix<std::complex<T> > cv3v2 = cv3 ^ v2;
  tmv::Matrix<std::complex<T> > cv2v3 = cv2 ^ v3;
  tmv::Matrix<std::complex<T> > v3cv2 = v3 ^ cv2;
  tmv::Matrix<std::complex<T> > v2cv3 = v2 ^ cv3;
  tmv::Matrix<std::complex<T> > c40 = c4;
  tmv::Matrix<std::complex<T> > c50 = c5;
  TestOProd(c4.View(),cv3,cv2,c40,cv3cv2,"NonSqaure CCC");
  TestOProd(c5.View(),cv2,cv3,c50,cv2cv3,"NonSqaure CCC");
  TestOProd(c4.View(),v3,cv2,c40,v3cv2,"NonSqaure CRC");
  TestOProd(c5.View(),v2,cv3,c50,v2cv3,"NonSqaure CRC");
  TestOProd(c4.View(),cv3,v2,c40,cv3v2,"NonSqaure CCR");
  TestOProd(c5.View(),cv2,v3,c50,cv2v3,"NonSqaure CCR");
  TestOProd(c4.View(),v3,v2,c40,v3v2,"NonSqaure CRR");
  TestOProd(c5.View(),v2,v3,c50,v2v3,"NonSqaure CRR");

#ifdef XTEST
  tmv::Matrix<T> a6(4,0,1);
  tmv::Matrix<T> a7(0,4,1);
  tmv::Matrix<std::complex<T> > c6 = a6;
  tmv::Matrix<std::complex<T> > c7 = a7;

  TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
      a1.View(),a6.View(),c1.View(),c6.View(),"Degenerate"); 
  TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
      a3.View(),a6.View(),c3.View(),c6.View(),"Degenerate"); 
  TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
      a4.View(),a6.View(),c4.View(),c6.View(),"Degenerate"); 
  TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
      a6.View(),a7.View(),c6.View(),c7.View(),"Degenerate"); 
  if (doallarith) {
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a1.View(),a7.View(),c1.View(),c7.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a2.View(),a6.View(),c2.View(),c6.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a2.View(),a7.View(),c2.View(),c7.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a3.View(),a7.View(),c3.View(),c7.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a4.View(),a7.View(),c4.View(),c7.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a5.View(),a6.View(),c5.View(),c6.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a5.View(),a7.View(),c5.View(),c7.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a6.View(),a1.View(),c6.View(),c1.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a6.View(),a2.View(),c6.View(),c2.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a6.View(),a3.View(),c6.View(),c3.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a6.View(),a4.View(),c6.View(),c4.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a6.View(),a5.View(),c6.View(),c5.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a7.View(),a1.View(),c7.View(),c1.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a7.View(),a2.View(),c7.View(),c2.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a7.View(),a3.View(),c7.View(),c3.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a7.View(),a4.View(),c7.View(),c4.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a7.View(),a5.View(),c7.View(),c5.View(),"Degenerate"); 
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a7.View(),a6.View(),c7.View(),c6.View(),"Degenerate"); 
  }

  tmv::Vector<T> v4 = a6.row(2);
  tmv::Matrix<T> v4v2 = v4 ^ v2;
  tmv::Matrix<T> v2v4 = v2 ^ v4;
  tmv::Matrix<T> a60 = a6;
  tmv::Matrix<T> a70 = a7;
  TestOProd(a6.View(),v2,v4,a60,v2v4,"Degenerate RRR");
  TestOProd(a7.View(),v4,v2,a70,v4v2,"Degenerate RRR");
  tmv::Vector<std::complex<T> > cv4 = c6.row(2);
  tmv::Matrix<std::complex<T> > cv4cv2 = cv4 ^ cv2;
  tmv::Matrix<std::complex<T> > cv2cv4 = cv2 ^ cv4;
  tmv::Matrix<std::complex<T> > cv4v2 = cv4 ^ v2;
  tmv::Matrix<std::complex<T> > cv2v4 = cv2 ^ v4;
  tmv::Matrix<std::complex<T> > v4cv2 = v4 ^ cv2;
  tmv::Matrix<std::complex<T> > v2cv4 = v2 ^ cv4;
  tmv::Matrix<std::complex<T> > c60 = c6;
  tmv::Matrix<std::complex<T> > c70 = c7;
  TestOProd(c6.View(),cv2,cv4,c60,cv2cv4,"Degenerate CCC");
  TestOProd(c7.View(),cv4,cv2,c70,cv4cv2,"Degenerate CCC");
  TestOProd(c6.View(),v2,cv4,c60,v2cv4,"Degenerate CRC");
  TestOProd(c7.View(),v4,cv2,c70,v4cv2,"Degenerate CRC");
  TestOProd(c6.View(),cv2,v4,c60,cv2v4,"Degenerate CCR");
  TestOProd(c7.View(),cv4,v2,c70,cv4v2,"Degenerate CCR");
  TestOProd(c6.View(),v2,v4,c60,v2v4,"Degenerate CRR");
  TestOProd(c7.View(),v4,v2,c70,v4v2,"Degenerate CRR");
#endif

}

template <class T> void TestAllMatrix()
{
  TestBasicMatrix<T,tmv::RowMajor>();
  TestBasicMatrix<T,tmv::ColMajor>();
  std::cout<<"Matrix<"<<tmv::Type(T())<<"> passed all basic tests\n";

  TestAllMatrixArith<T>();
  std::cout<<"Matrix<"<<tmv::Type(T())<<"> Arithmetic passed all tests\n";
}

template void TestAllMatrix<double>();
#ifndef NOFLOAT
template void TestAllMatrix<float>();
#endif
#ifdef LONGDOUBLE
template void TestAllMatrix<long double>();
#endif
