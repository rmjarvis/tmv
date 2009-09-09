#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"
#include <fstream>

template <class T, size_t M, size_t N, tmv::StorageType S> 
inline void TestBasicSmallMatrix()
{
  tmv::SmallMatrix<T,M,N,S> m;
  tmv::SmallMatrix<T,M,N,S,tmv::FortranStyle> mf;
  Assert(m.colsize() == size_t(M) && m.rowsize() == size_t(N),
      "Creating SmallMatrix(M,N)");
  Assert(m.colsize() == size_t(M) && m.rowsize() == size_t(N),
      "Creating SmallMatrixF(M,N)");

  for (size_t i=0, k=0; i<M; ++i) for (size_t j=0; j<N; ++j, ++k) {
    m(i,j) = k;
    mf(i+1,j+1) = k;
  }
#define Si (S==tmv::RowMajor ? int(N) : 1)
#define Sj (S==tmv::RowMajor ? 1 : int(M))
  tmv::ConstSmallMatrixView<T,M,N,Si,Sj> mcv = m.View();
  tmv::SmallMatrixView<T,M,N,Si,Sj> mv = m.View();
  tmv::ConstSmallMatrixView<T,M,N,Si,Sj,false,tmv::FortranStyle> mfcv = 
    mf.View();
  tmv::SmallMatrixView<T,M,N,Si,Sj,false,tmv::FortranStyle> mfv = mf.View();
#undef Si
#undef Sj

  for (size_t i=0, k=0; i<M; ++i) for (size_t j=0; j<N; ++j, ++k) {
    Assert(m(i,j) == k,"Read/Write SmallMatrix");
    Assert(mcv(i,j) == k,"Access SmallMatrix CV");
    Assert(mv(i,j) == k,"Access SmallMatrix V");
    Assert(mf(i+1,j+1) == k,"Read/Write SmallMatrixF");
    Assert(mfcv(i+1,j+1) == k,"Access SmallMatrixF CV");
    Assert(mfv(i+1,j+1) == k,"Access SmallMatrixF V");
    Assert(m[i][j] == k,"[] style access of SmallMatrix");
    Assert(mcv[i][j] == k,"[] style access of SmallMatrix CV");
    Assert(mv[i][j] == k,"[] style access of SmallMatrix V");
    Assert(mf[i+1][j+1] == k,"[] style access of SmallMatrixF");
    Assert(mfcv[i+1][j+1] == k,"[] style access of SmallMatrixF CV");
    Assert(mfv[i+1][j+1] == k,"[] style access of SmallMatrixF V");
    Assert(m.row(i)(j) == k,"SmallMatrix.row");
    Assert(mcv.row(i)(j) == k,"SmallMatrix.row CV");
    Assert(mv.row(i)(j) == k,"SmallMatrix.row V");
    Assert(mf.row(i+1)(j+1) == k,"SmallMatrixF.row");
    Assert(mfcv.row(i+1)(j+1) == k,"SmallMatrixF.row CV");
    Assert(mfv.row(i+1)(j+1) == k,"SmallMatrixF.row V");
    Assert(m.row(i,j,N)(0) == k,"SmallMatrix.row2");
    Assert(mcv.row(i,j,N)(0) == k,"SmallMatrix.row2 CV");
    Assert(mv.row(i,j,N)(0) == k,"SmallMatrix.row2 V");
    Assert(mf.row(i+1,j+1,N)(1) == k,"SmallMatrixF.row2");
    Assert(mfcv.row(i+1,j+1,N)(1) == k,"SmallMatrixF.row2 CV");
    Assert(mfv.row(i+1,j+1,N)(1) == k,"SmallMatrixF.row2 V");
    Assert(m.col(j)(i) == k,"SmallMatrix.col");
    Assert(mcv.col(j)(i) == k,"SmallMatrix.col CV");
    Assert(mv.col(j)(i) == k,"SmallMatrix.col V");
    Assert(mf.col(j+1)(i+1) == k,"SmallMatrixF.col");
    Assert(mfcv.col(j+1)(i+1) == k,"SmallMatrixF.col CV");
    Assert(mfv.col(j+1)(i+1) == k,"SmallMatrixF.col V");
    Assert(m.col(j,i,M)(0) == k,"SmallMatrix.col2");
    Assert(mcv.col(j,i,M)(0) == k,"SmallMatrix.col2 CV");
    Assert(mv.col(j,i,M)(0) == k,"SmallMatrix.col2 V");
    Assert(mf.col(j+1,i+1,M)(1) == k,"SmallMatrixF.col2");
    Assert(mfcv.col(j+1,i+1,M)(1) == k,"SmallMatrixF.col2 CV");
    Assert(mfv.col(j+1,i+1,M)(1) == k,"SmallMatrixF.col2 V");
    if (i<j) {
      Assert(m.diag(j-i)(i) == k,"SmallMatrix.diag");
      Assert(mcv.diag(j-i)(i) == k,"SmallMatrix.diag CV");
      Assert(mv.diag(j-i)(i) == k,"SmallMatrix.diag V");
      Assert(mf.diag(j-i)(i+1) == k,"SmallMatrixF.diag");
      Assert(mfcv.diag(j-i)(i+1) == k,"SmallMatrixF.diag CV");
      Assert(mfv.diag(j-i)(i+1) == k,"SmallMatrixF.diag V");
      Assert(m.diag(j-i,i,N-j+i)(0) == k,"SmallMatrix.diag2");
      Assert(mcv.diag(j-i,i,N-j+i)(0) == k,"SmallMatrix.diag2 CV");
      Assert(mv.diag(j-i,i,N-j+i)(0) == k,"SmallMatrix.diag2 V");
      Assert(mf.diag(j-i,i+1,N-j+i)(1) == k,"SmallMatrix.diag2");
      Assert(mfcv.diag(j-i,i+1,N-j+i)(1) == k,"SmallMatrix.diag2 CV");
      Assert(mfv.diag(j-i,i+1,N-j+i)(1) == k,"SmallMatrix.diag2 V");
    } else {
      if (i==j) {
	Assert(m.diag()(i) == k,"SmallMatrix.diag");
	Assert(mcv.diag()(i) == k,"SmallMatrix.diag CV");
	Assert(mv.diag()(i) == k,"SmallMatrix.diag V");
	Assert(mf.diag()(i+1) == k,"SmallMatrixF.diag");
	Assert(mfcv.diag()(i+1) == k,"SmallMatrixF.diag CV");
	Assert(mfv.diag()(i+1) == k,"SmallMatrixF.diag V");
      }
      Assert(m.diag(j-i)(j) == k,"SmallMatrix.diag1");
      Assert(mcv.diag(j-i)(j) == k,"SmallMatrix.diag1 CV");
      Assert(mv.diag(j-i)(j) == k,"SmallMatrix.diag1 V");
      Assert(mf.diag(j-i)(j+1) == k,"SmallMatrixF.diag1");
      Assert(mfcv.diag(j-i)(j+1) == k,"SmallMatrixF.diag1 CV");
      Assert(mfv.diag(j-i)(j+1) == k,"SmallMatrixF.diag1 V");
      if (N+i-j > M) {
	Assert(m.diag(j-i,j,M+j-i)(0) == k,"SmallMatrix.diag2");
	Assert(mcv.diag(j-i,j,M+j-i)(0) == k,"SmallMatrix.diag2 CV");
	Assert(mv.diag(j-i,j,M+j-i)(0) == k,"SmallMatrix.diag2 V");
	Assert(mf.diag(j-i,j+1,M+j-i)(1) == k,"SmallMatrix.diag2");
	Assert(mfcv.diag(j-i,j+1,M+j-i)(1) == k,"SmallMatrix.diag2 CV");
	Assert(mfv.diag(j-i,j+1,M+j-i)(1) == k,"SmallMatrix.diag2 V");
      } else {
	Assert(m.diag(j-i,j,N)(0) == k,"SmallMatrix.diag2");
	Assert(mcv.diag(j-i,j,N)(0) == k,"SmallMatrix.diag2 CV");
	Assert(mv.diag(j-i,j,N)(0) == k,"SmallMatrix.diag2 V");
	Assert(mf.diag(j-i,j+1,N)(1) == k,"SmallMatrix.diag2");
	Assert(mfcv.diag(j-i,j+1,N)(1) == k,"SmallMatrix.diag2 CV");
	Assert(mfv.diag(j-i,j+1,N)(1) == k,"SmallMatrix.diag2 V");
      }
    }
  }
  Assert(m == mf,"CStyle SmallMatrix == FortranStyle SmallMatrix");
  Assert(m == mcv,"SmallMatrix == ConstSmallMatrixView");
  Assert(m == mv,"SmallMatrix == SmallMatrixView");
  Assert(m == mfcv,"SmallMatrix == FortranStyle ConstSmallMatrixView");
  Assert(m == mfv,"SmallMatrix == FortranStyle SmallMatrixView");

  // Test Basic Arithmetic 
  
  tmv::SmallMatrix<T,M,N,S> a;
  tmv::SmallMatrix<T,M,N,S> b;
  tmv::SmallMatrix<T,M,N,S> c;
  for (size_t i=0; i<M; ++i) for (size_t j=0; j<N; ++j) {
    a(i,j) = 3.+i+5*j;
    b(i,j) = 5.+2*i+4*j;
  }
  mf = a;
  Assert(a == mf,"Copy CStyle SmallMatrix to FortranStyle");

  c = a+b;
  for (size_t i=0; i<M; ++i) for (size_t j=0; j<N; ++j) 
    Assert(c(i,j) == 8.+3*i+9*j,"Add Matrices");

  c = a-b;
  for (size_t i=0; i<M; ++i) for (size_t j=0; j<N; ++j) 
    Assert(c(i,j) == -2.-i+j,"Subtract Matrices");

  tmv::SmallMatrix<std::complex<T>,M,N,S> cm;
  tmv::SmallMatrix<std::complex<T>,M,N,S> ca;
  tmv::SmallMatrix<std::complex<T>,M,N,S> cb;
  Assert(cm.colsize() == size_t(M) && cm.rowsize() && size_t(N),
      "Creating CSmallMatrix(M,N)");

  for (size_t i=0, k=0; i<M; ++i) for (size_t j=0; j<N; ++j, ++k)
    cm(i,j) = std::complex<T>(k,k+1000);

  for (size_t i=0, k=0; i<M; ++i) for (size_t j=0; j<N; ++j, ++k) {
    Assert(cm(i,j) == std::complex<T>(k,k+1000),"Read/Write CSmallMatrix");
    Assert(cm.row(i)(j) == std::complex<T>(k,k+1000),"CSmallMatrix.row");
    Assert(cm.col(j)(i) == std::complex<T>(k,k+1000),"CSmallMatrix.col");
  }

  for (size_t i=0; i<M; ++i) for (size_t j=0; j<N; ++j) {
    ca(i,j) = std::complex<T>(3.+i+5*j,0.+i-j);
    cb(i,j) = std::complex<T>(3.+2*i+4*j,4.-10*i);
  }

  cm = ca+cb;
  for (size_t i=0; i<M; ++i) for (size_t j=0; j<N; ++j) 
    Assert(cm(i,j) == std::complex<T>(6.+3*i+9*j,4.-9*i-j),"Add CSmallMatrix");

  cm = ca-cb;
  for (size_t i=0; i<M; ++i) for (size_t j=0; j<N; ++j) 
    Assert(cm(i,j) == std::complex<T>(0.-i+j,-4.+11*i-j),"Subtract CSmallMatrix");

  cm = ca;
  for (size_t i=0; i<M; ++i) for (size_t j=0; j<N; ++j) 
    Assert(cm(i,j) == ca(i,j),"Copy CSmallMatrix");

  // Test I/O

  std::ofstream fout("tmvtest_smallmatrix_io.dat");
  if (!fout) throw std::runtime_error(
      "Couldn't open tmvtest_smallmatrix_io.dat for output");
  fout << m << std::endl << cm << std::endl;
  fout.close();

  tmv::SmallMatrix<T,M,N,tmv::RowMajor> xm1;
  tmv::SmallMatrix<std::complex<T>,M,N,tmv::RowMajor> xcm1;
  std::ifstream fin("tmvtest_smallmatrix_io.dat");
  if (!fin) throw std::runtime_error(
      "Couldn't open tmvtest_smallmatrix_io.dat for input");
  fin >> xm1 >> xcm1;
  fin.close();
  Assert(m == xm1,"SmallMatrix I/O check #1");
  Assert(cm == xcm1,"CSmallMatrix I/O check #1");

  tmv::SmallMatrix<T,M,N,tmv::ColMajor> xm2;
  tmv::SmallMatrix<std::complex<T>,M,N,tmv::ColMajor> xcm2;
  fin.open("tmvtest_smallmatrix_io.dat");
  if (!fin) throw std::runtime_error(
      "Couldn't open tmvtest_smallmatrix_io.dat for input");
  fin >> xm2 >> xcm2;
  fin.close();
  Assert(m == xm2,"SmallMatrix I/O check #2");
  Assert(cm == xcm2,"CSmallMatrix I/O check #2");
}

template <class T> void TestAllSmallMatrix()
{
  TestBasicSmallMatrix<T,15,10,tmv::RowMajor>();
  TestBasicSmallMatrix<T,15,10,tmv::ColMajor>();
  std::cout<<"SmallMatrix<"<<tmv::Type(T())<<"> passed all basic tests\n";

  TestSmallMatrix_Sub1<T>();
  TestSmallMatrix_Sub2<T>();
  TestSmallMatrix_Sub3<T>();
  TestSmallMatrix_Sub4<T>();
  TestSmallMatrix_Sub5<T>();
  std::cout<<"SmallMatrix<"<<tmv::Type(T())<<"> passed all subsets tests\n";
}

#ifdef INST_DOUBLE
template void TestAllSmallMatrix<double>();
#endif
#ifdef INST_FLOAT
template void TestAllSmallMatrix<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllSmallMatrix<long double>();
#endif
#ifdef INST_INT
template void TestAllSmallMatrix<int>();
#endif
