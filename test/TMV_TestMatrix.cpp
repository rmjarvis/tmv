// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV_Mat.h"
#include <fstream>
#include <cstdio>

#define CT std::complex<T>

template <class T, tmv::StorageType S> static void TestBasicMatrix_1()
{
  const int M = 15;
  const int N = 10;

  tmv::Matrix<T,S> m(M,N);
  tmv::MatrixF<T,S> mf(M,N);
  Assert(m.colsize() == size_t(M) && m.rowsize() == size_t(N),
      "Creating Matrix(M,N)");
  Assert(m.colsize() == size_t(M) && m.rowsize() == size_t(N),
      "Creating MatrixF(M,N)");

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
    m(i,j) = T(k);
    mf(i+1,j+1) = T(k);
  }
  tmv::ConstMatrixView<T> mcv = m.View();
  tmv::MatrixView<T> mv = m.View();
  tmv::ConstMatrixViewF<T> mfcv = mf.View();
  tmv::MatrixViewF<T> mfv = mf.View();

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
}

template <class T, tmv::StorageType S> static void TestBasicMatrix_2()
{
  const int M = 15;
  const int N = 10;

  tmv::Matrix<T,S> m(M,N);
  tmv::MatrixF<T,S> mf(M,N);

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
    m(i,j) = T(k);
    mf(i+1,j+1) = T(k);
  }
  tmv::ConstMatrixView<T> mcv = m.View();
  tmv::MatrixView<T> mv = m.View();
  tmv::ConstMatrixViewF<T> mfcv = mf.View();
  tmv::MatrixViewF<T> mfv = mf.View();

  Assert(m.SubMatrix(2,5,1,4) == m.SubMatrix(2,5,1,4,1,1),"SubMatrix");
  Assert(m.SubVector(2,5,4,2,3) == m.SubMatrix(2,14,5,11,4,2).diag(),
      "SubVector");
  Assert(m.ColPair(2,5) == m.SubMatrix(0,M,2,8,1,3),"ColPair");
  Assert(m.ColPair(7,2) == m.SubMatrix(0,M,7,-3,1,-5),"ColPair");
  Assert(m.RowPair(3,7) == m.SubMatrix(3,11,0,N,4,1),"RowPair");
  Assert(m.RowPair(2,0) == m.SubMatrix(2,-2,0,N,-2,1),"RowPair");
  Assert(m.Cols(2,5) == m.SubMatrix(0,M,2,5),"Cols");
  Assert(m.Rows(3,7) == m.SubMatrix(3,7,0,N),"Rows");

  Assert(mf.SubMatrix(3,5,2,4) == mf.SubMatrix(3,5,2,4,1,1),"SubMatrixFF");
  Assert(mf.SubVector(3,6,4,2,3) == mf.SubMatrix(3,11,6,10,4,2).diag(),
      "SubVectorFF");
  Assert(mf.ColPair(3,6) == mf.SubMatrix(1,M,3,6,1,3),"ColPairFF");
  Assert(mf.ColPair(8,3) == mf.SubMatrix(1,M,8,3,1,-5),"ColPairFF");
  Assert(mf.RowPair(4,8) == mf.SubMatrix(4,8,1,N,4,1),"RowPairFF");
  Assert(mf.RowPair(3,1) == mf.SubMatrix(3,1,1,N,-2,1),"RowPairFF");
  Assert(mf.Cols(3,5) == mf.SubMatrix(1,M,3,5),"ColsFF");
  Assert(mf.Rows(4,7) == mf.SubMatrix(4,7,1,N),"RowsFF");

  Assert(m.SubMatrix(2,5,1,4) == mf.SubMatrix(3,5,2,4),"SubMatrixF");
  Assert(m.SubMatrix(2,8,1,10,2,3) == mf.SubMatrix(3,7,2,8,2,3),"SubMatrixF");
  Assert(m.SubVector(2,5,4,2,3) == mf.SubVector(3,6,4,2,3),"SubVectorF");
  Assert(m.SubVector(8,1,-1,2,4) == mf.SubVector(9,2,-1,2,4),"SubVector2F");
  Assert(m.SubVector(12,8,-4,-2,2) == mf.SubVector(13,9,-4,-2,2),
      "SubVector3F");
  Assert(m.ColPair(2,5) == mf.ColPair(3,6),"ColPairF");
  Assert(m.ColPair(7,2) == mf.ColPair(8,3),"ColPairF");
  Assert(m.RowPair(3,7) == mf.RowPair(4,8),"RowPairF");
  Assert(m.RowPair(2,0) == mf.RowPair(3,1),"RowPairF");
  Assert(m.Cols(2,5) == mf.Cols(3,5),"ColsF");
  Assert(m.Rows(3,7) == mf.Rows(4,7),"RowsF");

  Assert(m.SubMatrix(2,5,1,4) == mcv.SubMatrix(2,5,1,4),"SubMatrixCV");
  Assert(m.SubMatrix(2,8,1,10,2,3) == mcv.SubMatrix(2,8,1,10,2,3),
      "SubMatrixCV");
  Assert(m.SubVector(2,5,4,2,3) == mcv.SubVector(2,5,4,2,3),"SubVectorCV");
  Assert(m.SubVector(8,1,-1,2,4) == mcv.SubVector(8,1,-1,2,4),"SubVector2CV");
  Assert(m.SubVector(12,8,-4,-2,2) == mcv.SubVector(12,8,-4,-2,2),
      "SubVector3CV");
  Assert(m.ColPair(2,5) == mcv.ColPair(2,5),"ColPairCV");
  Assert(m.ColPair(7,2) == mcv.ColPair(7,2),"ColPairCV");
  Assert(m.RowPair(3,7) == mcv.RowPair(3,7),"RowPairCV");
  Assert(m.RowPair(2,0) == mcv.RowPair(2,0),"RowPairCV");
  Assert(m.Cols(2,5) == mcv.Cols(2,5),"ColsCV");
  Assert(m.Rows(3,7) == mcv.Rows(3,7),"RowsCV");

  Assert(m.SubMatrix(2,5,1,4) == mv.SubMatrix(2,5,1,4),"SubMatrixV");
  Assert(m.SubMatrix(2,8,1,10,2,3) == mv.SubMatrix(2,8,1,10,2,3),"SubMatrixV");
  Assert(m.SubVector(2,5,4,2,3) == mv.SubVector(2,5,4,2,3),"SubVectorV");
  Assert(m.SubVector(8,1,-1,2,4) == mv.SubVector(8,1,-1,2,4),"SubVector2V");
  Assert(m.SubVector(12,8,-4,-2,2) == mv.SubVector(12,8,-4,-2,2),
      "SubVector3V");
  Assert(m.ColPair(2,5) == mv.ColPair(2,5),"ColPairV");
  Assert(m.ColPair(7,2) == mv.ColPair(7,2),"ColPairV");
  Assert(m.RowPair(3,7) == mv.RowPair(3,7),"RowPairV");
  Assert(m.RowPair(2,0) == mv.RowPair(2,0),"RowPairV");
  Assert(m.Cols(2,5) == mv.Cols(2,5),"ColsV");
  Assert(m.Rows(3,7) == mv.Rows(3,7),"RowsV");

  Assert(mf.SubMatrix(3,5,2,4) == mfcv.SubMatrix(3,5,2,4),"SubMatrixFCV");
  Assert(mf.SubMatrix(3,7,2,8,2,3) == mfcv.SubMatrix(3,7,2,8,2,3),
      "SubMatrixFCV");
  Assert(mf.SubVector(3,6,4,2,3) == mfcv.SubVector(3,6,4,2,3),"SubVectorFCV");
  Assert(mf.SubVector(9,2,-1,2,4) == mfcv.SubVector(9,2,-1,2,4),
      "SubVector2FCV");
  Assert(mf.SubVector(13,9,-4,-2,2) == mfcv.SubVector(13,9,-4,-2,2),
      "SubVector3FCV");
  Assert(mf.ColPair(3,6) == mfcv.ColPair(3,6),"ColPairFCV");
  Assert(mf.ColPair(8,3) == mfcv.ColPair(8,3),"ColPairFCV");
  Assert(mf.RowPair(4,8) == mfcv.RowPair(4,8),"RowPairFCV");
  Assert(mf.RowPair(3,1) == mfcv.RowPair(3,1),"RowPairFCV");
  Assert(mf.Cols(3,5) == mfcv.Cols(3,5),"ColsFCV");
  Assert(mf.Rows(4,7) == mfcv.Rows(4,7),"RowsFCV");

  Assert(mf.SubMatrix(3,5,2,4) == mfv.SubMatrix(3,5,2,4),"SubMatrixFV");
  Assert(mf.SubMatrix(3,7,2,8,2,3) == mfv.SubMatrix(3,7,2,8,2,3),"SubMatrixFV");
  Assert(mf.SubVector(3,6,4,2,3) == mfv.SubVector(3,6,4,2,3),"SubVectorFV");
  Assert(mf.SubVector(9,2,-1,2,4) == mfv.SubVector(9,2,-1,2,4),"SubVector2FV");
  Assert(mf.SubVector(13,9,-4,-2,2) == mfv.SubVector(13,9,-4,-2,2),
      "SubVector3FV");
  Assert(mf.ColPair(3,6) == mfv.ColPair(3,6),"ColPairFV");
  Assert(mf.ColPair(8,3) == mfv.ColPair(8,3),"ColPairFV");
  Assert(mf.RowPair(4,8) == mfv.RowPair(4,8),"RowPairFV");
  Assert(mf.RowPair(3,1) == mfv.RowPair(3,1),"RowPairFV");
  Assert(mf.Cols(3,5) == mfv.Cols(3,5),"ColsFV");
  Assert(mf.Rows(4,7) == mfv.Rows(4,7),"RowsFV");

  tmv::Matrix<T,S> a(M,N);
  tmv::Matrix<T,S> b(M,N);
  tmv::Matrix<T,S> c(M,N);
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) {
    a(i,j) = T(3+i+5*j);
    b(i,j) = T(5+2*i+4*j);
  }
  mf = a;
  Assert(a == mf,"Copy CStyle Matrix to FortranStyle");

  std::vector<T> qv(12);
  tmv::Matrix<T,S> q4(3,4);
  tmv::Matrix<T,S> q5t(4,3);
  tmv::MatrixView<T> q5 = q5t.Transpose();
  if (S == tmv::RowMajor) {
    T qvar[] = { 
      T(0), T(-1), T(-2), T(-3),
      T(2), T(1), T(0), T(-1),
      T(4), T(3), T(2), T(1) 
    };
    for(int i=0;i<12;i++) qv[i] = qvar[i];
    q4 <<
       0, -1, -2, -3,
       2, 1, 0, -1,
       4, 3, 2, 1;
    q5 <<
       0, 2, 4,
       -1, 1, 3,
       -2, 0, 2,
       -3, -1, 1;
  } else {
    T qvar[] = {
      T(0), T(2), T(4),
      T(-1), T(1), T(3),
      T(-2), T(0), T(2),
      T(-3), T(-1), T(1) 
    };
    for(int i=0;i<12;i++) qv[i] = qvar[i];
    q4 <<
       0, 2, 4,
       -1, 1, 3,
       -2, 0, 2,
       -3, -1, 1;
    q5 <<
       0, -1, -2, -3,
       2, 1, 0, -1,
       4, 3, 2, 1;
  }
  const int Si = (S == tmv::RowMajor ? 4 : 1);
  const int Sj = (S == tmv::RowMajor ? 1 : 3);
  T qar[12];
  for(int i=0;i<12;i++) qar[i] = qv[i];
  tmv::Matrix<T,S> q1(3,4,qar);
  tmv::Matrix<T,S> q2(3,4,qv);
  
  tmv::ConstMatrixView<T> q3 = tmv::MatrixViewOf(qar,3,4,S);
  tmv::ConstMatrixView<T,Si,Sj> q6 = tmv::MatrixViewOf(qar,3,4,Si,Sj);

  if (showacc) {
    std::cout<<"q1 = "<<q1<<std::endl;
    std::cout<<"q2 = "<<q2<<std::endl;
    std::cout<<"q3 = "<<q3<<std::endl;
    std::cout<<"q4 = "<<q4<<std::endl;
    std::cout<<"q5 = "<<q5<<std::endl;
    std::cout<<"q6 = "<<q6<<std::endl;
  }

  for(int i=0;i<3;i++) for(int j=0;j<4;j++) {
    Assert(q1(i,j) == T(2*i-j),"Create Matrix from T*");
    Assert(q2(i,j) == T(2*i-j),"Create Matrix from vector");
    Assert(q3(i,j) == T(2*i-j),"Create MatrixView of T* (S)");
    Assert(q4(i,j) == T(2*i-j),"Create Matrix from << list");
    Assert(q5(i,j) == T(2*i-j),"Create MatrixView from << list");
    Assert(q6(i,j) == T(2*i-j),"Create MatrixView of T* (Si,Sj)");
  }

  c = a+b;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(c(i,j) == T(8+3*i+9*j),"Add Matrices");

  c = a-b;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(c(i,j) == T(-2-i+j),"Subtract Matrices");

  tmv::Matrix<CT,S> cm(M,N);
  tmv::Matrix<CT,S> ca(M,N);
  tmv::Matrix<CT,S> cb(M,N);
  Assert(cm.colsize() == size_t(M) && cm.rowsize() && size_t(N),
      "Creating CMatrix(M,N)");

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
    cm(i,j) = CT(T(k),T(k+1000));

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
    Assert(cm(i,j) == CT(T(k),T(k+1000)),"Read/Write CMatrix");
    Assert(cm.row(i)(j) == CT(T(k),T(k+1000)),"CMatrix.row");
    Assert(cm.col(j)(i) == CT(T(k),T(k+1000)),"CMatrix.col");
  }

  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) {
    ca(i,j) = CT(T(3+i+5*j),T(i-j));
    cb(i,j) = CT(T(3+2*i+4*j),T(4-10*i));
  }

  cm = ca+cb;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(cm(i,j) == CT(T(6+3*i+9*j),T(4-9*i-j)),"Add CMatrix");

  cm = ca-cb;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(cm(i,j) == CT(T(-i+j),T(-4+11*i-j)),"Subtract CMatrix");

  cm = ca;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(cm(i,j) == ca(i,j),"Copy CMatrix");
}

template <class T, tmv::StorageType S> static void TestBasicMatrix_IO()
{
  const int M = 15;
  const int N = 10;

  tmv::Matrix<T,S> m(M,N);
  tmv::Matrix<CT,S> cm(M,N);

  for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
    m(i,j) = T(k);
    cm(i,j) = CT(T(k),T(k+1000));
  }

  std::ofstream fout("tmvtest_matrix_io.dat");
  if (!fout) 
#ifdef NOTHROW
  { std::cerr<<"Couldn't open tmvtest_matrix_io.dat for output\n"; exit(1); }
#else
  throw std::runtime_error(
      "Couldn't open tmvtest_matrix_io.dat for output");
#endif
  fout << m << std::endl << cm << std::endl;
  fout.close();

  tmv::Matrix<T,tmv::RowMajor> xm1(M,N);
  tmv::Matrix<CT,tmv::RowMajor> xcm1(M,N);
  std::ifstream fin("tmvtest_matrix_io.dat");
  if (!fin) 
#ifdef NOTHROW
  { std::cerr<<"Couldn't open tmvtest_matrix_io.dat for input\n"; exit(1); }
#else
  throw std::runtime_error(
      "Couldn't open tmvtest_matrix_io.dat for input");
#endif
  fin >> xm1 >> xcm1;
  fin.close();
  Assert(m == xm1,"Matrix I/O check #1");
  Assert(cm == xcm1,"CMatrix I/O check #1");

  tmv::Matrix<T,tmv::ColMajor> xm2(M,N);
  tmv::Matrix<CT,tmv::ColMajor> xcm2(M,N);
  fin.open("tmvtest_matrix_io.dat");
  if (!fin) 
#ifdef NOTHROW
  { std::cerr<<"Couldn't open tmvtest_matrix_io.dat for input\n"; exit(1); }
#else
  throw std::runtime_error(
      "Couldn't open tmvtest_matrix_io.dat for input");
#endif
  fin >> xm2 >> xcm2;
  fin.close();
  Assert(m == xm2,"Matrix I/O check #2");
  Assert(cm == xcm2,"CMatrix I/O check #2");

  std::auto_ptr<tmv::Matrix<T> > xm3;
  std::auto_ptr<tmv::Matrix<CT> > xcm3;
  fin.open("tmvtest_matrix_io.dat");
  if (!fin) 
#ifdef NOTHROW
  { std::cerr<<"Couldn't open tmvtest_matrix_io.dat for input\n"; exit(1); }
#else
  throw std::runtime_error(
      "Couldn't open tmvtest_matrix_io.dat for input");
#endif
  fin >> xm3 >> xcm3;
  fin.close();
  Assert(m == *xm3,"Matrix I/O check #3");
  Assert(cm == *xcm3,"CMatrix I/O check #3");

#ifndef XTEST
  std::remove("tmvtest_matrix_io.dat");
#endif

}

template <class T> void TestAllMatrix()
{
#if 1
  TestBasicMatrix_1<T,tmv::RowMajor>();
  TestBasicMatrix_1<T,tmv::ColMajor>();
  TestBasicMatrix_2<T,tmv::RowMajor>();
  TestBasicMatrix_2<T,tmv::ColMajor>();
  TestBasicMatrix_IO<T,tmv::RowMajor>();
  TestBasicMatrix_IO<T,tmv::ColMajor>();
  std::cout<<"Matrix<"<<tmv::TypeText(T())<<"> passed all basic tests\n";
#endif

#if 1
  TestMatrixArith_1<T>();
  TestMatrixArith_2<T>();
  TestMatrixArith_3<T>();
  TestMatrixArith_4<T>();
  TestMatrixArith_5<T>();
  TestMatrixArith_6<T>();
  TestMatrixArith_7<T>();
  TestMatrixArith_8<T>();
  std::cout<<"Matrix<"<<tmv::TypeText(T())<<"> Arithmetic passed all tests\n";
#endif
}

#ifdef TEST_DOUBLE
template void TestAllMatrix<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllMatrix<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllMatrix<long double>();
#endif
#ifdef TEST_INT
template void TestAllMatrix<int>();
#endif
