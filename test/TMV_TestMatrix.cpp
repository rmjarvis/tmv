#ifdef NDEBUG
#undef NDEBUG
#endif
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"
#include <fstream>
#include <cstdio>

//#define XTEST

#include "TMV_TestMatrixArith.h"
#define CT std::complex<T>

template <class T, tmv::StorageType S> static void TestBasicMatrix()
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
    m(i,j) = T(k);
    mf(i+1,j+1) = T(k);
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
    a(i,j) = 3.+i+5*j;
    b(i,j) = 5.+2*i+4*j;
  }
  mf = a;
  Assert(a == mf,"Copy CStyle Matrix to FortranStyle");

  std::vector<T> qv(6);
  if (S == tmv::RowMajor) {
    T qvar[] = { T(0), T(-1), T(-2),
                 T(2), T(1), T(0) };
    for(int i=0;i<6;i++) qv[i] = qvar[i];
  } else {
    T qvar[] = { T(0), T(2),
                 T(-1), T(1),
	         T(-2), T(0) };
    for(int i=0;i<6;i++) qv[i] = qvar[i];
  }
  T qar[6];
  for(int i=0;i<6;i++) qar[i] = qv[i];
  tmv::Matrix<T,S> q1(2,3,qar);
  tmv::Matrix<T,S> q2(2,3,qv);
  tmv::ConstMatrixView<T> q3 = tmv::MatrixViewOf(qar,2,3,S);
  for(int i=0;i<2;i++) for(int j=0;j<3;j++) {
    Assert(q1(i,j) == T(2*i-j),"Create Matrix from T*");
    Assert(q2(i,j) == T(2*i-j),"Create Matrix from vector");
    Assert(q3(i,j) == T(2*i-j),"Create MatrixView of T*");
  }

  c = a+b;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(c(i,j) == 8.+3*i+9*j,"Add Matrices");

  c = a-b;
  for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
    Assert(c(i,j) == -2.-i+j,"Subtract Matrices");

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

  // Test I/O

  std::ofstream fout("tmvtest_matrix_io.dat");
  if (!fout) throw std::runtime_error(
      "Couldn't open tmvtest_matrix_io.dat for output");
  fout << m << std::endl << cm << std::endl;
  fout.close();

  tmv::Matrix<T,tmv::RowMajor> xm1(M,N);
  tmv::Matrix<CT,tmv::RowMajor> xcm1(M,N);
  std::ifstream fin("tmvtest_matrix_io.dat");
  if (!fin) throw std::runtime_error(
      "Couldn't open tmvtest_matrix_io.dat for input");
  fin >> xm1 >> xcm1;
  fin.close();
  Assert(m == xm1,"Matrix I/O check #1");
  Assert(cm == xcm1,"CMatrix I/O check #1");

  tmv::Matrix<T,tmv::ColMajor> xm2(M,N);
  tmv::Matrix<CT,tmv::ColMajor> xcm2(M,N);
  fin.open("tmvtest_matrix_io.dat");
  if (!fin) throw std::runtime_error(
      "Couldn't open tmvtest_matrix_io.dat for input");
  fin >> xm2 >> xcm2;
  fin.close();
  Assert(m == xm2,"Matrix I/O check #2");
  Assert(cm == xcm2,"CMatrix I/O check #2");

  std::auto_ptr<tmv::Matrix<T> > xm3;
  std::auto_ptr<tmv::Matrix<CT> > xcm3;
  fin.open("tmvtest_matrix_io.dat");
  if (!fin) throw std::runtime_error(
      "Couldn't open tmvtest_matrix_io.dat for input");
  fin >> xm3 >> xcm3;
  fin.close();
  Assert(m == *xm3,"Matrix I/O check #3");
  Assert(cm == *xcm3,"CMatrix I/O check #3");

#ifndef XTEST
  std::remove("tmvtest_matrix_io.dat");
  //system("rm tmvtest_matrix_io.dat");
#endif

}

template <class T> static void TestAllMatrixArith()
{
  // First test that algorithms work for some relevant different sizes
#ifdef XTEST
#if 0
  const int NSIZE = 16;
  const int sizear[NSIZE] = {1,2,3,4,5,6,7,8,9,30,31,32,63,111,128,137};
#else
  const int NSIZE = 4;
  const int sizear[NSIZE] = {1,2,63,137};
#endif
  for(int m1=0;m1<NSIZE;m1++) {
    for(int n1=0;n1<NSIZE;n1++) {
      int m = sizear[m1];
      int n = sizear[n1];
      std::cout<<"m,n = "<<m<<','<<n<<std::endl;
      tmv::Matrix<CT,tmv::ColMajor> c(m,n);
      tmv::Matrix<CT,tmv::ColMajor> c1(m,n);

      // Test various adds:
      {
	tmv::Matrix<CT,tmv::ColMajor> ac(m,n);
	tmv::Matrix<CT,tmv::ColMajor> bc(m,n);
	for(int i=0;i<m;i++) for(int j=0;j<n;j++) 
	  ac(i,j) = CT(T(i+5),T(j-2));
	for(int i=0;i<m;i++) for(int j=0;j<n;j++) 
	  bc(i,j) = CT(T(2*i-13),T(-3*j+8));
	tmv::Matrix<CT,tmv::RowMajor> ar = ac;
	tmv::Matrix<CT,tmv::RowMajor> br = bc;
	T eps = EPS*(T(1)+Norm(ac)+Norm(bc));

	c1 = ac + bc;
	Assert(Norm((c=ar+br)-c1) < eps,"ar+br");
	Assert(Norm((c=ar+bc)-c1) < eps,"ar+bc");
	Assert(Norm((c=ac+br)-c1) < eps,"ac+br");
	c1 = ac + tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm((c=ac+bc.Conjugate())-c1) < eps,"ac+cbc");
	Assert(Norm((c=ar+br.Conjugate())-c1) < eps,"ar+cbr");
	Assert(Norm((c=ar+bc.Conjugate())-c1) < eps,"ar+cbc");
	Assert(Norm((c=ac+br.Conjugate())-c1) < eps,"ac+cbr");
	c1 = tmv::Matrix<CT>(ac.Conjugate()) + bc;
	Assert(Norm((c=ac.Conjugate()+bc)-c1) < eps,"cac+bc");
	Assert(Norm((c=ar.Conjugate()+br)-c1) < eps,"car+br");
	Assert(Norm((c=ar.Conjugate()+bc)-c1) < eps,"car+bc");
	Assert(Norm((c=ac.Conjugate()+br)-c1) < eps,"cac+br");
	c1 = tmv::Matrix<CT>(ac.Conjugate()) + tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm((c=ac.Conjugate()+bc.Conjugate())-c1) < eps,"cac+cbc");
	Assert(Norm((c=ar.Conjugate()+br.Conjugate())-c1) < eps,"car+cbr");
	Assert(Norm((c=ar.Conjugate()+bc.Conjugate())-c1) < eps,"car+cbc");
	Assert(Norm((c=ac.Conjugate()+br.Conjugate())-c1) < eps,"cac+cbr");

	T x1(7);
	T x2(-3);
	c1 = x1*ac + x2*bc;
	Assert(Norm((c=x1*ar+x2*br)-c1) < eps,"xar+xbr");
	Assert(Norm((c=x1*ar+x2*bc)-c1) < eps,"xar+xbc");
	Assert(Norm((c=x1*ac+x2*br)-c1) < eps,"xac+xbr");
	c1 = x1*ac + x2*tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm((c=x1*ac+x2*bc.Conjugate())-c1) < eps,"xac+xcbc");
	Assert(Norm((c=x1*ar+x2*br.Conjugate())-c1) < eps,"xar+xcbr");
	Assert(Norm((c=x1*ar+x2*bc.Conjugate())-c1) < eps,"xar+xcbc");
	Assert(Norm((c=x1*ac+x2*br.Conjugate())-c1) < eps,"xac+xcbr");
	c1 = x1*tmv::Matrix<CT>(ac.Conjugate()) + x2*bc;
	Assert(Norm((c=x1*ac.Conjugate()+x2*bc)-c1) < eps,"xcac+xbc");
	Assert(Norm((c=x1*ar.Conjugate()+x2*br)-c1) < eps,"xcar+xbr");
	Assert(Norm((c=x1*ar.Conjugate()+x2*bc)-c1) < eps,"xcar+xbc");
	Assert(Norm((c=x1*ac.Conjugate()+x2*br)-c1) < eps,"xcac+xbr");
	c1 = x1*tmv::Matrix<CT>(ac.Conjugate()) + x2*tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm((c=x1*ac.Conjugate()+x2*bc.Conjugate())-c1) < eps,
	    "xcac+xcbc");
	Assert(Norm((c=x1*ar.Conjugate()+x2*br.Conjugate())-c1) < eps,
	    "xcar+xcbr");
	Assert(Norm((c=x1*ar.Conjugate()+x2*bc.Conjugate())-c1) < eps,
	    "xcar+xcbc");
	Assert(Norm((c=x1*ac.Conjugate()+x2*br.Conjugate())-c1) < eps,
	    "xcac+xcbr");

	CT z1(7,-9);
	CT z2(-3,-8);
	c1 = z1*ac + z2*bc;
	Assert(Norm((c=z1*ar+z2*br)-c1) < eps,"zar+zbr");
	Assert(Norm((c=z1*ar+z2*bc)-c1) < eps,"zar+zbc");
	Assert(Norm((c=z1*ac+z2*br)-c1) < eps,"zac+zbr");
	c1 = z1*ac + z2*tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm((c=z1*ac+z2*bc.Conjugate())-c1) < eps,"zac+zcbc");
	Assert(Norm((c=z1*ar+z2*br.Conjugate())-c1) < eps,"zar+zcbr");
	Assert(Norm((c=z1*ar+z2*bc.Conjugate())-c1) < eps,"zar+zcbc");
	Assert(Norm((c=z1*ac+z2*br.Conjugate())-c1) < eps,"zac+zcbr");
	c1 = z1*tmv::Matrix<CT>(ac.Conjugate()) + z2*bc;
	Assert(Norm((c=z1*ac.Conjugate()+z2*bc)-c1) < eps,"zcac+zbc");
	Assert(Norm((c=z1*ar.Conjugate()+z2*br)-c1) < eps,"zcar+zbr");
	Assert(Norm((c=z1*ar.Conjugate()+z2*bc)-c1) < eps,"zcar+zbc");
	Assert(Norm((c=z1*ac.Conjugate()+z2*br)-c1) < eps,"zcac+zbr");
	c1 = z1*tmv::Matrix<CT>(ac.Conjugate()) + z2*tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm((c=z1*ac.Conjugate()+z2*bc.Conjugate())-c1) < eps,
	    "zcac+zcbc");
	Assert(Norm((c=z1*ar.Conjugate()+z2*br.Conjugate())-c1) < eps,
	    "zcar+zcbr");
	Assert(Norm((c=z1*ar.Conjugate()+z2*bc.Conjugate())-c1) < eps,
	    "zcar+zcbc");
	Assert(Norm((c=z1*ac.Conjugate()+z2*br.Conjugate())-c1) < eps,
	    "zcac+zcbr");
      }

      // Now test multiplies
      for(int k1=0;k1<NSIZE;k1++) {
	int k = sizear[k1];
	tmv::Matrix<CT,tmv::ColMajor> ac(m,k);
	tmv::Matrix<CT,tmv::ColMajor> bc(k,n);
	for(int i=0;i<m;i++) for(int j=0;j<k;j++) 
	  ac(i,j) = CT(T(i+5),T(j-2));
	for(int i=0;i<k;i++) for(int j=0;j<n;j++) 
	  bc(i,j) = CT(T(2*i-13),T(-3*j+8));
	tmv::Matrix<CT,tmv::RowMajor> ar = ac;
	tmv::Matrix<CT,tmv::RowMajor> br = bc;
	T eps = T(10) * EPS * (T(1) + Norm(ac)*Norm(bc));

	c1 = ac * bc;
	Assert(Norm((c=ar*br)-c1) < eps,"ar*br");
	Assert(Norm((c=ar*bc)-c1) < eps,"ar*bc");
	Assert(Norm((c=ac*br)-c1) < eps,"ac*br");
	c1 = ac * tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm((c=ac*bc.Conjugate())-c1) < eps,"ac*cbc");
	Assert(Norm((c=ar*br.Conjugate())-c1) < eps,"ar*cbr");
	Assert(Norm((c=ar*bc.Conjugate())-c1) < eps,"ar*cbc");
	Assert(Norm((c=ac*br.Conjugate())-c1) < eps,"ac*cbr");
	c1 = tmv::Matrix<CT>(ac.Conjugate()) * bc;
	Assert(Norm((c=ac.Conjugate()*bc)-c1) < eps,"cac*bc");
	Assert(Norm((c=ar.Conjugate()*br)-c1) < eps,"car*br");
	Assert(Norm((c=ar.Conjugate()*bc)-c1) < eps,"car*bc");
	Assert(Norm((c=ac.Conjugate()*br)-c1) < eps,"cac*br");
	c1 = tmv::Matrix<CT>(ac.Conjugate()) * tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm((c=ac.Conjugate()*bc.Conjugate())-c1) < eps,"cac*cbc");
	Assert(Norm((c=ar.Conjugate()*br.Conjugate())-c1) < eps,"car*cbr");
	Assert(Norm((c=ar.Conjugate()*bc.Conjugate())-c1) < eps,"car*cbc");
	Assert(Norm((c=ac.Conjugate()*br.Conjugate())-c1) < eps,"cac*cbr");

	tmv::Matrix<CT> c0 = c1;
	c1 = c0 + ac * bc;
	Assert(Norm(((c=c0)+=ar*br)-c1) < eps,"+ar*br");
	Assert(Norm(((c=c0)+=ar*bc)-c1) < eps,"+ar*bc");
	Assert(Norm(((c=c0)+=ac*br)-c1) < eps,"+ac*br");
	c1 = c0 + ac * tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm(((c=c0)+=ac*bc.Conjugate())-c1) < eps,"+ac*cbc");
	Assert(Norm(((c=c0)+=ar*br.Conjugate())-c1) < eps,"+ar*cbr");
	Assert(Norm(((c=c0)+=ar*bc.Conjugate())-c1) < eps,"+ar*cbc");
	Assert(Norm(((c=c0)+=ac*br.Conjugate())-c1) < eps,"+ac*cbr");
	c1 = c0 + tmv::Matrix<CT>(ac.Conjugate()) * bc;
	Assert(Norm(((c=c0)+=ac.Conjugate()*bc)-c1) < eps,"+cac*bc");
	Assert(Norm(((c=c0)+=ar.Conjugate()*br)-c1) < eps,"+car*br");
	Assert(Norm(((c=c0)+=ar.Conjugate()*bc)-c1) < eps,"+car*bc");
	Assert(Norm(((c=c0)+=ac.Conjugate()*br)-c1) < eps,"+cac*br");
	c1 = c0 + tmv::Matrix<CT>(ac.Conjugate()) * tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm(((c=c0)+=ac.Conjugate()*bc.Conjugate())-c1) < eps,"+cac*cbc");
	Assert(Norm(((c=c0)+=ar.Conjugate()*br.Conjugate())-c1) < eps,"+car*cbr");
	Assert(Norm(((c=c0)+=ar.Conjugate()*bc.Conjugate())-c1) < eps,"+car*cbc");
	Assert(Norm(((c=c0)+=ac.Conjugate()*br.Conjugate())-c1) < eps,"+cac*cbr");

	T x1(7);
	eps *= x1;
	c1 = x1*ac * bc;
	Assert(Norm((c=x1*ar*br)-c1) < eps,"xar*br");
	Assert(Norm((c=x1*ar*bc)-c1) < eps,"xar*bc");
	Assert(Norm((c=x1*ac*br)-c1) < eps,"xac*br");
	c1 = x1*ac * tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm((c=x1*ac*bc.Conjugate())-c1) < eps,"xac*cbc");
	Assert(Norm((c=x1*ar*br.Conjugate())-c1) < eps,"xar*cbr");
	Assert(Norm((c=x1*ar*bc.Conjugate())-c1) < eps,"xar*cbc");
	Assert(Norm((c=x1*ac*br.Conjugate())-c1) < eps,"xac*cbr");
	c1 = x1*tmv::Matrix<CT>(ac.Conjugate()) * bc;
	Assert(Norm((c=x1*ac.Conjugate()*bc)-c1) < eps,"xcac*bc");
	Assert(Norm((c=x1*ar.Conjugate()*br)-c1) < eps,"xcar*br");
	Assert(Norm((c=x1*ar.Conjugate()*bc)-c1) < eps,"xcar*bc");
	Assert(Norm((c=x1*ac.Conjugate()*br)-c1) < eps,"xcac*br");
	c1 = x1*tmv::Matrix<CT>(ac.Conjugate()) * tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm((c=x1*ac.Conjugate()*bc.Conjugate())-c1) < eps,
	    "xcac*cbc");
	Assert(Norm((c=x1*ar.Conjugate()*br.Conjugate())-c1) < eps,
	    "xcar*cbr");
	Assert(Norm((c=x1*ar.Conjugate()*bc.Conjugate())-c1) < eps,
	    "xcar*cbc");
	Assert(Norm((c=x1*ac.Conjugate()*br.Conjugate())-c1) < eps,
	    "xcac*cbr");

	c1 = c0 + x1*ac * bc;
	Assert(Norm(((c=c0)+=x1*ar*br)-c1) < eps,"+xar*br");
	Assert(Norm(((c=c0)+=x1*ar*bc)-c1) < eps,"+xar*bc");
	Assert(Norm(((c=c0)+=x1*ac*br)-c1) < eps,"+xac*br");
	c1 = c0 + x1*ac * tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm(((c=c0)+=x1*ac*bc.Conjugate())-c1) < eps,"+xac*cbc");
	Assert(Norm(((c=c0)+=x1*ar*br.Conjugate())-c1) < eps,"+xar*cbr");
	Assert(Norm(((c=c0)+=x1*ar*bc.Conjugate())-c1) < eps,"+xar*cbc");
	Assert(Norm(((c=c0)+=x1*ac*br.Conjugate())-c1) < eps,"+xac*cbr");
	c1 = c0 + x1*tmv::Matrix<CT>(ac.Conjugate()) * bc;
	Assert(Norm(((c=c0)+=x1*ac.Conjugate()*bc)-c1) < eps,"+xcac*bc");
	Assert(Norm(((c=c0)+=x1*ar.Conjugate()*br)-c1) < eps,"+xcar*br");
	Assert(Norm(((c=c0)+=x1*ar.Conjugate()*bc)-c1) < eps,"+xcar*bc");
	Assert(Norm(((c=c0)+=x1*ac.Conjugate()*br)-c1) < eps,"+xcac*br");
	c1 = c0 + x1*tmv::Matrix<CT>(ac.Conjugate()) * tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm(((c=c0)+=x1*ac.Conjugate()*bc.Conjugate())-c1) < eps,
	    "+xcac*cbc");
	Assert(Norm(((c=c0)+=x1*ar.Conjugate()*br.Conjugate())-c1) < eps,
	    "+xcar*cbr");
	Assert(Norm(((c=c0)+=x1*ar.Conjugate()*bc.Conjugate())-c1) < eps,
	    "+xcar*cbc");
	Assert(Norm(((c=c0)+=x1*ac.Conjugate()*br.Conjugate())-c1) < eps,
	    "+xcac*cbr");

	CT z1(7,-9);
	eps *= std::norm(z1)/x1;
	c1 = z1*ac * bc;
	Assert(Norm((c=z1*ar*br)-c1) < eps,"zar*br");
	Assert(Norm((c=z1*ar*bc)-c1) < eps,"zar*bc");
	Assert(Norm((c=z1*ac*br)-c1) < eps,"zac*br");
	c1 = z1*ac * tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm((c=z1*ac*bc.Conjugate())-c1) < eps,"zac*cbc");
	Assert(Norm((c=z1*ar*br.Conjugate())-c1) < eps,"zar*cbr");
	Assert(Norm((c=z1*ar*bc.Conjugate())-c1) < eps,"zar*cbc");
	Assert(Norm((c=z1*ac*br.Conjugate())-c1) < eps,"zac*cbr");
	c1 = z1*tmv::Matrix<CT>(ac.Conjugate()) * bc;
	Assert(Norm((c=z1*ac.Conjugate()*bc)-c1) < eps,"zcac*bc");
	Assert(Norm((c=z1*ar.Conjugate()*br)-c1) < eps,"zcar*br");
	Assert(Norm((c=z1*ar.Conjugate()*bc)-c1) < eps,"zcar*bc");
	Assert(Norm((c=z1*ac.Conjugate()*br)-c1) < eps,"zcac*br");
	c1 = z1*tmv::Matrix<CT>(ac.Conjugate()) * tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm((c=z1*ac.Conjugate()*bc.Conjugate())-c1) < eps,
	    "zcac*cbc");
	Assert(Norm((c=z1*ar.Conjugate()*br.Conjugate())-c1) < eps,
	    "zcar*cbr");
	Assert(Norm((c=z1*ar.Conjugate()*bc.Conjugate())-c1) < eps,
	    "zcar*cbc");
	Assert(Norm((c=z1*ac.Conjugate()*br.Conjugate())-c1) < eps,
	    "zcac*cbr");

	c1 = c0 + z1*ac * bc;
	Assert(Norm(((c=c0)+=z1*ar*br)-c1) < eps,"+zar*br");
	Assert(Norm(((c=c0)+=z1*ar*bc)-c1) < eps,"+zar*bc");
	Assert(Norm(((c=c0)+=z1*ac*br)-c1) < eps,"+zac*br");
	c1 = c0 + z1*ac * tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm(((c=c0)+=z1*ac*bc.Conjugate())-c1) < eps,"+zac*cbc");
	Assert(Norm(((c=c0)+=z1*ar*br.Conjugate())-c1) < eps,"+zar*cbr");
	Assert(Norm(((c=c0)+=z1*ar*bc.Conjugate())-c1) < eps,"+zar*cbc");
	Assert(Norm(((c=c0)+=z1*ac*br.Conjugate())-c1) < eps,"+zac*cbr");
	c1 = c0 + z1*tmv::Matrix<CT>(ac.Conjugate()) * bc;
	Assert(Norm(((c=c0)+=z1*ac.Conjugate()*bc)-c1) < eps,"+zcac*bc");
	Assert(Norm(((c=c0)+=z1*ar.Conjugate()*br)-c1) < eps,"+zcar*br");
	Assert(Norm(((c=c0)+=z1*ar.Conjugate()*bc)-c1) < eps,"+zcar*bc");
	Assert(Norm(((c=c0)+=z1*ac.Conjugate()*br)-c1) < eps,"+zcac*br");
	c1 = c0 + z1*tmv::Matrix<CT>(ac.Conjugate()) * tmv::Matrix<CT>(bc.Conjugate());
	Assert(Norm(((c=c0)+=z1*ac.Conjugate()*bc.Conjugate())-c1) < eps,
	    "+zcac*cbc");
	Assert(Norm(((c=c0)+=z1*ar.Conjugate()*br.Conjugate())-c1) < eps,
	    "+zcar*cbr");
	Assert(Norm(((c=c0)+=z1*ar.Conjugate()*bc.Conjugate())-c1) < eps,
	    "+zcar*cbc");
	Assert(Norm(((c=c0)+=z1*ac.Conjugate()*br.Conjugate())-c1) < eps,
	    "+zcac*cbr");
      }
    }
  }
#endif

  // Now we use the TestMatrixArith.h file to test lots of different
  // syntaxes for calling matrix arithmetic.  This tests the inline
  // parser more than the algorithms.
  tmv::Matrix<T,tmv::RowMajor> a1(4,4);
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) {
    a1(i,j) = T(2+4*i-5*j);
  }
  a1(0,0) = 14.; 
  a1(1,0) = -2.; 
  a1(2,0) = 7.; 
  a1(3,0) = -10.;
  a1(2,2) = 30.;

  tmv::Matrix<CT,tmv::RowMajor> ca1 = a1;
  ca1(2,3) += CT(2,3);
  ca1(1,0) *= CT(0,2);
  ca1.col(1) *= CT(-1,3);
  ca1.row(3) += tmv::Vector<CT>(4,CT(1,9));
  tmv::MatrixView<T> a1v = a1.View();
  tmv::MatrixView<CT> ca1v = ca1.View();

  tmv::Matrix<T,tmv::ColMajor> a2 = a1.Transpose();
  a2.row(1) *= T(3);
  a2.col(2) -= tmv::Vector<T>(4,4.);
  tmv::Matrix<CT,tmv::ColMajor> ca2 = ca1;
  ca2 -= a2;
  ca2 *= CT(1,-2);
  tmv::MatrixView<T> a2v = a2.View();
  tmv::MatrixView<CT> ca2v = ca2.View();

  tmv::Matrix<T> a3x(12,16);
  for(int i=0;i<12;++i) for(int j=0;j<16;++j) a3x(i,j) = T(1-2*i+3*j);
  a3x.diag().AddToAll(30);
  tmv::Matrix<CT> ca3x = a3x*CT(1,-2);
  ca3x.diag().AddToAll(CT(-22,15));
  tmv::MatrixView<T> a3v = a3x.SubMatrix(0,12,0,16,3,4);
  tmv::MatrixView<CT> ca3v = ca3x.SubMatrix(0,12,0,16,3,4);

  tmv::Matrix<T> a1x(4,4);
  tmv::Matrix<CT> ca1x(4,4);

  TestMatrixArith123<T>(a1x,ca1x,a1v,ca1v,"Square");
  TestMatrixArith123<T>(a1x,ca1x,a2v,ca2v,"Square");
  TestMatrixArith123<T>(a1x,ca1x,a3v,ca3v,"Square");
  TestMatrixArith45<T>(a1x,ca1x,a1v,ca1v,a2v,ca2v,"Square");
  TestMatrixArith45<T>(a1x,ca1x,a2v,ca2v,a1v,ca1v,"Square");
  TestMatrixArith45<T>(a1x,ca1x,a3v,ca3v,a1v,ca1v,"Square");
#ifdef XTEST
  TestMatrixArith45<T>(a1x,ca1x,a1v,ca1v,a3v,ca3v,"Square");
  TestMatrixArith45<T>(a1x,ca1x,a3v,ca3v,a2v,ca2v,"Square");
  TestMatrixArith45<T>(a1x,ca1x,a2v,ca2v,a3v,ca3v,"Square");
#endif

  tmv::Vector<T> v1 = a1.col(0);
  tmv::VectorView<T> v1v = v1.View();
  tmv::Vector<T> v15(20);
  tmv::VectorView<T> v1s = v15.SubVector(0,20,5);
  v1s = v1v;

  tmv::Vector<T> v2 = a1.row(2);
  tmv::VectorView<T> v2v = v2.View();
  tmv::Vector<T> v25(20);
  tmv::VectorView<T> v2s = v25.SubVector(0,20,5);
  v2s = v2v;

  tmv::Vector<CT> cv1 = ca1.col(0);
  tmv::VectorView<CT> cv1v = cv1.View();
  tmv::Vector<CT> cv15(20);
  tmv::VectorView<CT> cv1s = cv15.SubVector(0,20,5);
  cv1s = cv1v;

  tmv::Vector<CT> cv2 = ca1.row(2);
  tmv::VectorView<CT> cv2v = cv2.View();
  tmv::Vector<CT> cv25(20);
  tmv::VectorView<CT> cv2s = cv25.SubVector(0,20,5);
  cv2s = cv2v;

  TestMatrixArith6<T>(a1x,ca1x,a1v,ca1v,v1v,cv1v,v2v,cv2v,"Square");
  TestMatrixArith6<T>(a1x,ca1x,a1v,ca1v,v1s,cv1s,v2v,cv2v,"Square");
  TestMatrixArith6<T>(a1x,ca1x,a1v,ca1v,v1v,cv1v,v2s,cv2s,"Square");
  TestMatrixArith6<T>(a1x,ca1x,a1v,ca1v,v1s,cv1s,v2s,cv2s,"Square");
  TestMatrixArith6<T>(a1x,ca1x,a2v,ca2v,v1v,cv1v,v2v,cv2v,"Square");
  TestMatrixArith6<T>(a1x,ca1x,a2v,ca2v,v1s,cv1s,v2v,cv2v,"Square");
  TestMatrixArith6<T>(a1x,ca1x,a2v,ca2v,v1v,cv1v,v2s,cv2s,"Square");
  TestMatrixArith6<T>(a1x,ca1x,a2v,ca2v,v1s,cv1s,v2s,cv2s,"Square");
  TestMatrixArith6<T>(a1x,ca1x,a3v,ca3v,v1v,cv1v,v2v,cv2v,"Square");
  TestMatrixArith6<T>(a1x,ca1x,a3v,ca3v,v1s,cv1s,v2v,cv2v,"Square");
  TestMatrixArith6<T>(a1x,ca1x,a3v,ca3v,v1v,cv1v,v2s,cv2s,"Square");
  TestMatrixArith6<T>(a1x,ca1x,a3v,ca3v,v1s,cv1s,v2s,cv2s,"Square");

  tmv::Matrix<T,tmv::RowMajor> a4(7,4);
  for(int i=0;i<7;++i) for(int j=0;j<4;++j) a4(i,j) = T(1-3*i+2*j);
  tmv::Matrix<T,tmv::ColMajor> a5 = a4.Transpose();
  a4.SubMatrix(2,6,0,4) += a1;
  a5.SubMatrix(0,4,1,5) -= a2;
  tmv::MatrixView<T> a4v = a4.View();
  tmv::MatrixView<T> a5v = a5.View();

  tmv::Matrix<CT,tmv::RowMajor> ca4 = a4*CT(1,2);
  tmv::Matrix<CT,tmv::ColMajor> ca5 = ca4.Adjoint();
  ca4.SubMatrix(2,6,0,4) += ca1;
  ca5.SubMatrix(0,4,1,5) -= ca2;
  ca4.col(1) *= CT(2,1);
  ca4.row(6).AddToAll(CT(-7,2));
  ca5.col(3) *= CT(-1,3);
  ca5.row(0).AddToAll(CT(1,9));
  tmv::MatrixView<CT> ca4v = ca4.View();
  tmv::MatrixView<CT> ca5v = ca5.View();

  tmv::Matrix<T> a4x(7,4);
  tmv::Matrix<CT> ca4x(7,4);
  tmv::Matrix<T> a5x(4,7);
  tmv::Matrix<CT> ca5x(4,7);
  TestMatrixArith123<T>(a4x,ca4x,a4v,ca4v,"NonSquare");
  TestMatrixArith123<T>(a5x,ca5x,a5v,ca5v,"NonSquare");
  TestMatrixArith45<T>(a1x,ca1x,a1v,ca1v,a4v,ca4v,"NonSquare");
  TestMatrixArith45<T>(a4x,ca4x,a4v,ca4v,a1v,ca1v,"NonSquare");
  TestMatrixArith45<T>(a4x,ca4x,a4v,ca4v,a5v,ca5v,"NonSquare");
  TestMatrixArith45<T>(a5x,ca5x,a5v,ca5v,a4v,ca4v,"NonSquare");
#ifdef XTEST
  TestMatrixArith45<T>(a1x,ca1x,a2v,ca2v,a4v,ca4v,"NonSquare");
  TestMatrixArith45<T>(a1x,ca1x,a3v,ca3v,a4v,ca4v,"NonSquare");
  TestMatrixArith45<T>(a1x,ca1x,a1v,ca1v,a5v,ca4v,"NonSquare");
  TestMatrixArith45<T>(a1x,ca1x,a2v,ca2v,a5v,ca4v,"NonSquare");
  TestMatrixArith45<T>(a1x,ca1x,a3v,ca3v,a5v,ca4v,"NonSquare");
  TestMatrixArith45<T>(a4x,ca4x,a4v,ca4v,a2v,ca2v,"NonSquare");
  TestMatrixArith45<T>(a4x,ca4x,a4v,ca4v,a3v,ca3v,"NonSquare");
  TestMatrixArith45<T>(a5x,ca5x,a5v,ca5v,a1v,ca1v,"NonSquare");
  TestMatrixArith45<T>(a5x,ca5x,a5v,ca5v,a2v,ca2v,"NonSquare");
  TestMatrixArith45<T>(a5x,ca5x,a5v,ca5v,a3v,ca3v,"NonSquare");
#endif

  tmv::Vector<T> v3 = a4.col(2);
  tmv::VectorView<T> v3v = v3.View();
  tmv::Vector<T> v35(35);
  tmv::VectorView<T> v3s = v35.SubVector(0,35,5);
  v3s = v3v;

  tmv::Vector<CT> cv3 = ca4.col(2);
  tmv::VectorView<CT> cv3v = cv3.View();
  tmv::Vector<CT> cv35(35);
  tmv::VectorView<CT> cv3s = cv35.SubVector(0,35,5);
  cv3s = cv3v;

  TestMatrixArith6<T>(a4x,ca4x,a4v,ca4v,v3v,cv3v,v2v,cv2v,"NonSquare");
  TestMatrixArith6<T>(a4x,ca4x,a4v,ca4v,v3s,cv3s,v2v,cv2v,"NonSquare");
  TestMatrixArith6<T>(a4x,ca4x,a4v,ca4v,v3v,cv3v,v2s,cv2s,"NonSquare");
  TestMatrixArith6<T>(a4x,ca4x,a4v,ca4v,v3s,cv3s,v2s,cv2s,"NonSquare");
  TestMatrixArith6<T>(a5x,ca5x,a5v,ca5v,v1v,cv1v,v3v,cv3v,"NonSquare");
  TestMatrixArith6<T>(a5x,ca5x,a5v,ca5v,v1s,cv1s,v3v,cv3v,"NonSquare");
  TestMatrixArith6<T>(a5x,ca5x,a5v,ca5v,v1v,cv1v,v3s,cv3s,"NonSquare");
  TestMatrixArith6<T>(a5x,ca5x,a5v,ca5v,v1s,cv1s,v3s,cv3s,"NonSquare");

#ifdef XTEST
  tmv::Matrix<T> a6(4,0,1);
  tmv::Matrix<T> a7(0,4,1);
  tmv::Matrix<CT> ca6 = a6;
  tmv::Matrix<CT> ca7 = a7;
  tmv::MatrixView<T> a6v = a6.View();
  tmv::MatrixView<T> a7v = a7.View();
  tmv::MatrixView<CT> ca6v = ca6.View();
  tmv::MatrixView<CT> ca7v = ca7.View();

  tmv::Matrix<T> a6x(4,0);
  tmv::Matrix<T> a7x(0,4);
  tmv::Matrix<CT> ca6x(4,0);
  tmv::Matrix<CT> ca7x(0,4);
  TestMatrixArith123<T>(a6x,ca6x,a6v,ca6v,"Degenerate");
  TestMatrixArith123<T>(a7x,ca7x,a7v,ca7v,"Degenerate");
  TestMatrixArith45<T>(a1x,ca1x,a1v,ca1v,a6v,ca6v,"Degenerate");
  TestMatrixArith45<T>(a6x,ca6x,a6v,ca6v,a1v,ca1v,"Degenerate");
  TestMatrixArith45<T>(a6x,ca6x,a6v,ca6v,a7v,ca7v,"Degenerate");
  TestMatrixArith45<T>(a7x,ca7x,a7v,ca7v,a6v,ca6v,"Degenerate");
  TestMatrixArith45<T>(a1x,ca1x,a2v,ca2v,a6v,ca6v,"Degenerate");
  TestMatrixArith45<T>(a1x,ca1x,a3v,ca3v,a6v,ca6v,"Degenerate");
  TestMatrixArith45<T>(a4x,ca4x,a4v,ca4v,a6v,ca6v,"Degenerate");
  TestMatrixArith45<T>(a5x,ca5x,a5v,ca5v,a6v,ca6v,"Degenerate");
  TestMatrixArith45<T>(a1x,ca1x,a1v,ca1v,a7v,ca6v,"Degenerate");
  TestMatrixArith45<T>(a1x,ca1x,a2v,ca2v,a7v,ca6v,"Degenerate");
  TestMatrixArith45<T>(a1x,ca1x,a3v,ca3v,a7v,ca6v,"Degenerate");
  TestMatrixArith45<T>(a4x,ca4x,a4v,ca4v,a7v,ca6v,"Degenerate");
  TestMatrixArith45<T>(a5x,ca5x,a5v,ca5v,a7v,ca6v,"Degenerate");
  TestMatrixArith45<T>(a6x,ca6x,a6v,ca6v,a2v,ca2v,"Degenerate");
  TestMatrixArith45<T>(a6x,ca6x,a6v,ca6v,a3v,ca3v,"Degenerate");
  TestMatrixArith45<T>(a6x,ca6x,a6v,ca6v,a4v,ca4v,"Degenerate");
  TestMatrixArith45<T>(a6x,ca6x,a6v,ca6v,a5v,ca5v,"Degenerate");
  TestMatrixArith45<T>(a7x,ca7x,a7v,ca7v,a1v,ca1v,"Degenerate");
  TestMatrixArith45<T>(a7x,ca7x,a7v,ca7v,a2v,ca2v,"Degenerate");
  TestMatrixArith45<T>(a7x,ca7x,a7v,ca7v,a3v,ca3v,"Degenerate");
  TestMatrixArith45<T>(a7x,ca7x,a7v,ca7v,a4v,ca4v,"Degenerate");
  TestMatrixArith45<T>(a7x,ca7x,a7v,ca7v,a5v,ca5v,"Degenerate");

  tmv::Vector<T> v4 = a6.row(2);
  tmv::VectorView<T> v4v = v4.View();
  tmv::Vector<T> v45(0);
  tmv::VectorView<T> v4s = v45.SubVector(0,0,5);

  tmv::Vector<CT> cv4 = ca6.row(2);
  tmv::VectorView<CT> cv4v = cv4.View();
  tmv::Vector<CT> cv45(0);
  tmv::VectorView<CT> cv4s = cv45.SubVector(0,0,5);

  TestMatrixArith6<T>(a6x,ca6x,a6v,ca6v,v1v,cv1v,v4v,cv4v,"Degenerate");
  TestMatrixArith6<T>(a6x,ca6x,a6v,ca6v,v1s,cv1s,v4v,cv4v,"Degenerate");
  TestMatrixArith6<T>(a6x,ca6x,a6v,ca6v,v1v,cv1v,v4s,cv4s,"Degenerate");
  TestMatrixArith6<T>(a6x,ca6x,a6v,ca6v,v1s,cv1s,v4s,cv4s,"Degenerate");
  TestMatrixArith6<T>(a7x,ca7x,a7v,ca7v,v4v,cv4v,v2v,cv2v,"Degenerate");
  TestMatrixArith6<T>(a7x,ca7x,a7v,ca7v,v4s,cv4s,v2v,cv2v,"Degenerate");
  TestMatrixArith6<T>(a7x,ca7x,a7v,ca7v,v4v,cv4v,v2s,cv2s,"Degenerate");
  TestMatrixArith6<T>(a7x,ca7x,a7v,ca7v,v4s,cv4s,v2s,cv2s,"Degenerate");

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

#ifdef INST_DOUBLE
template void TestAllMatrix<double>();
#endif
#ifdef INST_FLOAT
template void TestAllMatrix<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllMatrix<long double>();
#endif
#ifdef INST_INT
template void TestAllMatrix<int>();
#endif
