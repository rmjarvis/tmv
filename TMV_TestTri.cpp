
#include "TMV_Test.h"
#include "TMV.h"
#include "TMV_Tri.h"

using tmv::Matrix;
using tmv::UpperTriMatrix;
using tmv::LowerTriMatrix;
using tmv::UpperTriMatrixView;
using tmv::LowerTriMatrixView;
using tmv::ConstUpperTriMatrixView;
using tmv::ConstLowerTriMatrixView;
using tmv::NonUnitDiag;
using tmv::UnitDiag;
using tmv::FortranStyle;

template <class T> extern void TestTriMatrixArith_A();
template <class T> extern void TestTriMatrixArith_B();

template <class T, tmv::DiagType D, tmv::StorageType S> void TestBasicTriMatrix()
{
  const int N = 10;

  UpperTriMatrix<T,D,S> u(N);
  LowerTriMatrix<T,D,S> l(N);
  UpperTriMatrix<T,D,S,FortranStyle> uf(N);
  LowerTriMatrix<T,D,S,FortranStyle> lf(N);

  Assert(u.colsize() == size_t(N) && u.rowsize() == size_t(N),
      "Creating UpperTriMatrix(N)");
  Assert(l.colsize() == size_t(N) && l.rowsize() == size_t(N),
      "Creating LowerTriMatrix(N)");

  for (int i=0,k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
    if (i < j || (D==NonUnitDiag && i==j)) {
      u(i,j) = T(k);
      uf(i+1,j+1) = T(k);
    }
    if (j < i || (D==NonUnitDiag && i==j)) {
      l(i,j) = T(k);
      lf(i+1,j+1) = T(k);
    }
  }

  UpperTriMatrixView<T> uv = u.View();
  ConstUpperTriMatrixView<T> ucv = u.View();
  UpperTriMatrixView<T,FortranStyle> ufv = uf.View();
  ConstUpperTriMatrixView<T,FortranStyle> ufcv = uf.View();
  LowerTriMatrixView<T> lv = l.View();
  ConstLowerTriMatrixView<T> lcv = l.View();
  LowerTriMatrixView<T,FortranStyle> lfv = lf.View();
  ConstLowerTriMatrixView<T,FortranStyle> lfcv = lf.View();

  const UpperTriMatrix<T,D,S>& ux = u;
  const LowerTriMatrix<T,D,S>& lx = l;
  const UpperTriMatrix<T,D,S,FortranStyle>& ufx = uf;
  const LowerTriMatrix<T,D,S,FortranStyle>& lfx = lf;

  for(int i=0,k=1;i<N;++i) for(int j=0;j<N;++j,++k) {
    if (i < j) {
      Assert(u(i,j) == T(k),"Read/Write TriMatrix");
      Assert(ux(i,j) == T(k),"Access const TriMatrix");
      Assert(ucv(i,j) == T(k),"Access TriMatrix CV");
      Assert(uv(i,j) == T(k),"Access TriMatrix V");
      Assert(uf(i+1,j+1) == T(k),"Read/Write TriMatrixF");
      Assert(ufx(i+1,j+1) == T(k),"Access const TriMatrixF");
      Assert(ufcv(i+1,j+1) == T(k),"Access TriMatrixF CV");
      Assert(ufv(i+1,j+1) == T(k),"Access TriMatrixF V");
      Assert(lx(i,j) == T(0),"Access const TriMatrix");
      Assert(lcv(i,j) == T(0),"Access TriMatrix CV");
      Assert(lfx(i+1,j+1) == T(0),"Access const TriMatrixF");
      Assert(lfcv(i+1,j+1) == T(0),"Access TriMatrixF CV");
      Assert(u.row(i,i+1,N)(j-i-1) == T(k),"TriMatrix.row1");
      Assert(ux.row(i,i+1,N)(j-i-1) == T(k),"const TriMatrix.row1");
      Assert(ucv.row(i,i+1,N)(j-i-1) == T(k),"TriMatrix.row1 CV");
      Assert(uv.row(i,i+1,N)(j-i-1) == T(k),"TriMatrix.row1 V");
      Assert(uf.row(i+1,i+2,N)(j-i) == T(k),"TriMatrixF.row1");
      Assert(ufx.row(i+1,i+2,N)(j-i) == T(k),"const TriMatrixF.row1");
      Assert(ufcv.row(i+1,i+2,N)(j-i) == T(k),"TriMatrixF.row1 CV");
      Assert(ufv.row(i+1,i+2,N)(j-i) == T(k),"TriMatrixF.row1 V");
      Assert(u.row(i,j,N)(0) == T(k),"TriMatrix.row2");
      Assert(ux.row(i,j,N)(0) == T(k),"const TriMatrix.row2");
      Assert(ucv.row(i,j,N)(0) == T(k),"TriMatrix.row2 CV");
      Assert(uv.row(i,j,N)(0) == T(k),"TriMatrix.row2 V");
      Assert(uf.row(i+1,j+1,N)(1) == T(k),"TriMatrixF.row2");
      Assert(ufx.row(i+1,j+1,N)(1) == T(k),"const TriMatrixF.row2");
      Assert(ufcv.row(i+1,j+1,N)(1) == T(k),"TriMatrixF.row2 CV");
      Assert(ufv.row(i+1,j+1,N)(1) == T(k),"TriMatrixF.row2 V");
      Assert(u.col(j,0,j)(i) == T(k),"TriMatrix.col1");
      Assert(ux.col(j,0,j)(i) == T(k),"const TriMatrix.col1");
      Assert(ucv.col(j,0,j)(i) == T(k),"TriMatrix.col1 CV");
      Assert(uv.col(j,0,j)(i) == T(k),"TriMatrix.col1 V");
      Assert(uf.col(j+1,1,j)(i+1) == T(k),"TriMatrixF.col1");
      Assert(ufx.col(j+1,1,j)(i+1) == T(k),"const TriMatrixF.col1");
      Assert(ufcv.col(j+1,1,j)(i+1) == T(k),"TriMatrixF.col1 CV");
      Assert(ufv.col(j+1,1,j)(i+1) == T(k),"TriMatrixF.col1 V");
      Assert(u.col(j,i,j)(0) == T(k),"TriMatrix.col2");
      Assert(ux.col(j,i,j)(0) == T(k),"const TriMatrix.col2");
      Assert(ucv.col(j,i,j)(0) == T(k),"TriMatrix.col2 CV");
      Assert(uv.col(j,i,j)(0) == T(k),"TriMatrix.col2 V");
      Assert(uf.col(j+1,i+1,j)(1) == T(k),"TriMatrixF.col2");
      Assert(ufx.col(j+1,i+1,j)(1) == T(k),"const TriMatrixF.col2");
      Assert(ufcv.col(j+1,i+1,j)(1) == T(k),"TriMatrixF.col2 CV");
      Assert(ufv.col(j+1,i+1,j)(1) == T(k),"TriMatrixF.col2 V");
      Assert(u.diag(j-i)(i) == T(k),"TriMatrix.diag1");
      Assert(ux.diag(j-i)(i) == T(k),"const TriMatrix.diag1");
      Assert(ucv.diag(j-i)(i) == T(k),"TriMatrix.diag1 CV ");
      Assert(uv.diag(j-i)(i) == T(k),"TriMatrix.diag1 V");
      Assert(uf.diag(j-i)(i+1) == T(k),"TriMatrixF.diag1");
      Assert(ufx.diag(j-i)(i+1) == T(k),"const TriMatrixF.diag1");
      Assert(ufcv.diag(j-i)(i+1) == T(k),"TriMatrixF.diag1 CV ");
      Assert(ufv.diag(j-i)(i+1) == T(k),"TriMatrixF.diag1 V");
      Assert(u.diag(j-i,i,N-j+i)(0) == T(k),"TriMatrix.diag2");
      Assert(ux.diag(j-i,i,N-j+i)(0) == T(k),"const TriMatrix.diag2");
      Assert(ucv.diag(j-i,i,N-j+i)(0) == T(k),"TriMatrix.diag2 CV ");
      Assert(uv.diag(j-i,i,N-j+i)(0) == T(k),"TriMatrix.diag2 V");
      Assert(uf.diag(j-i,i+1,N-j+i)(1) == T(k),"TriMatrixF.diag2");
      Assert(ufx.diag(j-i,i+1,N-j+i)(1) == T(k),"const TriMatrixF.diag2");
      Assert(ufcv.diag(j-i,i+1,N-j+i)(1) == T(k),"TriMatrixF.diag2 CV ");
      Assert(ufv.diag(j-i,i+1,N-j+i)(1) == T(k),"TriMatrixF.diag2 V");
    } else if (i==j) {
      if (D == UnitDiag) {
	Assert(ux(i,i) == T(1),"Access const TriMatrix");
	Assert(ucv(i,i) == T(1),"Access TriMatrix CV");
	Assert(ufx(i+1,i+1) == T(1),"Access const TriMatrixF");
	Assert(ufcv(i+1,i+1) == T(1),"Access TriMatrixF CV");
	Assert(lx(i,i) == T(1),"Access const TriMatrix");
	Assert(lcv(i,i) == T(1),"Access TriMatrix CV");
	Assert(lfx(i+1,i+1) == T(1),"Access const TriMatrixF");
	Assert(lfcv(i+1,i+1) == T(1),"Access TriMatrixF CV");
      } else {
	Assert(u(i,i) == T(k),"Read/Write TriMatrix");
	Assert(ux(i,i) == T(k),"Access const TriMatrix");
	Assert(ucv(i,i) == T(k),"Access TriMatrix CV");
	Assert(uv(i,i) == T(k),"Access TriMatrix V");
	Assert(uf(i+1,i+1) == T(k),"Read/Write TriMatrixF");
	Assert(ufx(i+1,i+1) == T(k),"Access const TriMatrixF");
	Assert(ufcv(i+1,i+1) == T(k),"Access TriMatrixF CV");
	Assert(ufv(i+1,i+1) == T(k),"Access TriMatrixF V");
	Assert(l(i,i) == T(k),"Read/Write TriMatrix");
	Assert(lx(i,i) == T(k),"Access const TriMatrix");
	Assert(lcv(i,i) == T(k),"Access TriMatrix CV");
	Assert(lv(i,i) == T(k),"Access TriMatrix V");
	Assert(lf(i+1,i+1) == T(k),"Read/Write TriMatrixF");
	Assert(lfx(i+1,i+1) == T(k),"Access const TriMatrixF");
	Assert(lfcv(i+1,i+1) == T(k),"Access TriMatrixF CV");
	Assert(lfv(i+1,i+1) == T(k),"Access TriMatrixF V");
	Assert(u.row(i,i,N)(0) == T(k),"TriMatrix.row1");
	Assert(ux.row(i,i,N)(0) == T(k),"const TriMatrix.row1");
	Assert(ucv.row(i,i,N)(0) == T(k),"TriMatrix.row1 CV");
	Assert(uv.row(i,i,N)(0) == T(k),"TriMatrix.row1 V");
	Assert(uf.row(i+1,i+1,N)(1) == T(k),"TriMatrixF.row1");
	Assert(ufx.row(i+1,i+1,N)(1) == T(k),"const TriMatrixF.row1");
	Assert(ufcv.row(i+1,i+1,N)(1) == T(k),"TriMatrixF.row1 CV");
	Assert(ufv.row(i+1,i+1,N)(1) == T(k),"TriMatrixF.row1 V");
	Assert(u.col(i,0,i+1)(i) == T(k),"TriMatrix.col1");
	Assert(ux.col(i,0,i+1)(i) == T(k),"const TriMatrix.col1");
	Assert(ucv.col(i,0,i+1)(i) == T(k),"TriMatrix.col1 CV");
	Assert(uv.col(i,0,i+1)(i) == T(k),"TriMatrix.col1 V");
	Assert(uf.col(i+1,1,i+1)(i+1) == T(k),"TriMatrixF.col1");
	Assert(ufx.col(i+1,1,i+1)(i+1) == T(k),"const TriMatrixF.col1");
	Assert(ufcv.col(i+1,1,i+1)(i+1) == T(k),"TriMatrixF.col1 CV");
	Assert(ufv.col(i+1,1,i+1)(i+1) == T(k),"TriMatrixF.col1 V");
	Assert(u.col(i,i,i+1)(0) == T(k),"TriMatrix.col2");
	Assert(ux.col(i,i,i+1)(0) == T(k),"const TriMatrix.col2");
	Assert(ucv.col(i,i,i+1)(0) == T(k),"TriMatrix.col2 CV");
	Assert(uv.col(i,i,i+1)(0) == T(k),"TriMatrix.col2 V");
	Assert(uf.col(i+1,i+1,i+1)(1) == T(k),"TriMatrixF.col2");
	Assert(ufx.col(i+1,i+1,i+1)(1) == T(k),"const TriMatrixF.col2");
	Assert(ufcv.col(i+1,i+1,i+1)(1) == T(k),"TriMatrixF.col2 CV");
	Assert(ufv.col(i+1,i+1,i+1)(1) == T(k),"TriMatrixF.col2 V");
	Assert(u.diag()(i) == T(k),"TriMatrix.diag");
	Assert(ux.diag()(i) == T(k),"const TriMatrix.diag");
	Assert(ucv.diag()(i) == T(k),"TriMatrix.diag CV ");
	Assert(uv.diag()(i) == T(k),"TriMatrix.diag V");
	Assert(uf.diag()(i+1) == T(k),"TriMatrixF.diag");
	Assert(ufx.diag()(i+1) == T(k),"const TriMatrixF.diag");
	Assert(ufcv.diag()(i+1) == T(k),"TriMatrixF.diag CV ");
	Assert(ufv.diag()(i+1) == T(k),"TriMatrixF.diag V");
	Assert(u.diag(0)(i) == T(k),"TriMatrix.diag1");
	Assert(ux.diag(0)(i) == T(k),"const TriMatrix.diag1");
	Assert(ucv.diag(0)(i) == T(k),"TriMatrix.diag1 CV ");
	Assert(uv.diag(0)(i) == T(k),"TriMatrix.diag1 V");
	Assert(uf.diag(0)(i+1) == T(k),"TriMatrixF.diag1");
	Assert(ufx.diag(0)(i+1) == T(k),"const TriMatrixF.diag1");
	Assert(ufcv.diag(0)(i+1) == T(k),"TriMatrixF.diag1 CV ");
	Assert(ufv.diag(0)(i+1) == T(k),"TriMatrixF.diag1 V");
	Assert(u.diag(0,i,N-i+i)(0) == T(k),"TriMatrix.diag2");
	Assert(ux.diag(0,i,N-i+i)(0) == T(k),"const TriMatrix.diag2");
	Assert(ucv.diag(0,i,N-i+i)(0) == T(k),"TriMatrix.diag2 CV ");
	Assert(uv.diag(0,i,N-i+i)(0) == T(k),"TriMatrix.diag2 V");
	Assert(uf.diag(0,i+1,N-i+i)(1) == T(k),"TriMatrixF.diag2");
	Assert(ufx.diag(0,i+1,N-i+i)(1) == T(k),"const TriMatrixF.diag2");
	Assert(ufcv.diag(0,i+1,N-i+i)(1) == T(k),"TriMatrixF.diag2 CV ");
	Assert(ufv.diag(0,i+1,N-i+i)(1) == T(k),"TriMatrixF.diag2 V");
	Assert(l.col(i,i,N)(0) == T(k),"TriMatrix.col2");
	Assert(lx.col(i,i,N)(0) == T(k),"const TriMatrix.col2");
	Assert(lcv.col(i,i,N)(0) == T(k),"TriMatrix.col2 CV");
	Assert(lv.col(i,i,N)(0) == T(k),"TriMatrix.col2 V");
	Assert(lf.col(i+1,i+1,N)(1) == T(k),"TriMatrixF.col2");
	Assert(lfx.col(i+1,i+1,N)(1) == T(k),"const TriMatrixF.col2");
	Assert(lfcv.col(i+1,i+1,N)(1) == T(k),"TriMatrixF.col2 CV");
	Assert(lfv.col(i+1,i+1,N)(1) == T(k),"TriMatrixF.col2 V");
	Assert(l.row(i,0,i+1)(i) == T(k),"TriMatrix.row1");
	Assert(lx.row(i,0,i+1)(i) == T(k),"const TriMatrix.row1");
	Assert(lcv.row(i,0,i+1)(i) == T(k),"TriMatrix.row1 CV");
	Assert(lv.row(i,0,i+1)(i) == T(k),"TriMatrix.row1 V");
	Assert(lf.row(i+1,1,i+1)(i+1) == T(k),"TriMatrixF.row1");
	Assert(lfx.row(i+1,1,i+1)(i+1) == T(k),"const TriMatrixF.row1");
	Assert(lfcv.row(i+1,1,i+1)(i+1) == T(k),"TriMatrixF.row1 CV");
	Assert(lfv.row(i+1,1,i+1)(i+1) == T(k),"TriMatrixF.row1 V");
	Assert(l.row(i,i,i+1)(0) == T(k),"TriMatrix.row2");
	Assert(lx.row(i,i,i+1)(0) == T(k),"const TriMatrix.row2");
	Assert(lcv.row(i,i,i+1)(0) == T(k),"TriMatrix.row2 CV");
	Assert(lv.row(i,i,i+1)(0) == T(k),"TriMatrix.row2 V");
	Assert(lf.row(i+1,i+1,i+1)(1) == T(k),"TriMatrixF.row2");
	Assert(lfx.row(i+1,i+1,i+1)(1) == T(k),"const TriMatrixF.row2");
	Assert(lfcv.row(i+1,i+1,i+1)(1) == T(k),"TriMatrixF.row2 CV");
	Assert(lfv.row(i+1,i+1,i+1)(1) == T(k),"TriMatrixF.row2 V");
	Assert(l.diag()(i) == T(k),"TriMatrix.diag");
	Assert(lx.diag()(i) == T(k),"const TriMatrix.diag");
	Assert(lcv.diag()(i) == T(k),"TriMatrix.diag CV ");
	Assert(lv.diag()(i) == T(k),"TriMatrix.diag V");
	Assert(lf.diag()(i+1) == T(k),"TriMatrixF.diag");
	Assert(lfx.diag()(i+1) == T(k),"const TriMatrixF.diag");
	Assert(lfcv.diag()(i+1) == T(k),"TriMatrixF.diag CV ");
	Assert(lfv.diag()(i+1) == T(k),"TriMatrixF.diag V");
	Assert(l.diag(0)(i) == T(k),"TriMatrix.diag1");
	Assert(lx.diag(0)(i) == T(k),"const TriMatrix.diag1");
	Assert(lcv.diag(0)(i) == T(k),"TriMatrix.diag1 CV ");
	Assert(lv.diag(0)(i) == T(k),"TriMatrix.diag1 V");
	Assert(lf.diag(0)(i+1) == T(k),"TriMatrixF.diag1");
	Assert(lfx.diag(0)(i+1) == T(k),"const TriMatrixF.diag1");
	Assert(lfcv.diag(0)(i+1) == T(k),"TriMatrixF.diag1 CV ");
	Assert(lfv.diag(0)(i+1) == T(k),"TriMatrixF.diag1 V");
	Assert(l.diag(0,i,N-i+i)(0) == T(k),"TriMatrix.diag2");
	Assert(lx.diag(0,i,N-i+i)(0) == T(k),"const TriMatrix.diag2");
	Assert(lcv.diag(0,i,N-i+i)(0) == T(k),"TriMatrix.diag2 CV ");
	Assert(lv.diag(0,i,N-i+i)(0) == T(k),"TriMatrix.diag2 V");
	Assert(lf.diag(0,i+1,N-i+i)(1) == T(k),"TriMatrixF.diag2");
	Assert(lfx.diag(0,i+1,N-i+i)(1) == T(k),"const TriMatrixF.diag2");
	Assert(lfcv.diag(0,i+1,N-i+i)(1) == T(k),"TriMatrixF.diag2 CV ");
	Assert(lfv.diag(0,i+1,N-i+i)(1) == T(k),"TriMatrixF.diag2 V");
      }
    } else {
      Assert(ux(i,j) == T(0),"Access const TriMatrix");
      Assert(ucv(i,j) == T(0),"Access TriMatrix CV");
      Assert(ufx(i+1,j+1) == T(0),"Access const TriMatrixF");
      Assert(ufcv(i+1,j+1) == T(0),"Access TriMatrixF CV");
      Assert(l(i,j) == T(k),"Read/Write TriMatrix");
      Assert(lx(i,j) == T(k),"Access const TriMatrix");
      Assert(lcv(i,j) == T(k),"Access TriMatrix CV");
      Assert(lv(i,j) == T(k),"Access TriMatrix V");
      Assert(lf(i+1,j+1) == T(k),"Read/Write TriMatrixF");
      Assert(lfx(i+1,j+1) == T(k),"Access const TriMatrixF");
      Assert(lfcv(i+1,j+1) == T(k),"Access TriMatrixF CV");
      Assert(lfv(i+1,j+1) == T(k),"Access TriMatrixF V");
      Assert(l.col(j,j+1,N)(i-j-1) == T(k),"TriMatrix.col1");
      Assert(lx.col(j,j+1,N)(i-j-1) == T(k),"const TriMatrix.col1");
      Assert(lcv.col(j,j+1,N)(i-j-1) == T(k),"TriMatrix.col1 CV");
      Assert(lv.col(j,j+1,N)(i-j-1) == T(k),"TriMatrix.col1 V");
      Assert(lf.col(j+1,j+2,N)(i-j) == T(k),"TriMatrixF.col1");
      Assert(lfx.col(j+1,j+2,N)(i-j) == T(k),"const TriMatrixF.col1");
      Assert(lfcv.col(j+1,j+2,N)(i-j) == T(k),"TriMatrixF.col1 CV");
      Assert(lfv.col(j+1,j+2,N)(i-j) == T(k),"TriMatrixF.col1 V");
      Assert(l.col(j,i,N)(0) == T(k),"TriMatrix.col2");
      Assert(lx.col(j,i,N)(0) == T(k),"const TriMatrix.col2");
      Assert(lcv.col(j,i,N)(0) == T(k),"TriMatrix.col2 CV");
      Assert(lv.col(j,i,N)(0) == T(k),"TriMatrix.col2 V");
      Assert(lf.col(j+1,i+1,N)(1) == T(k),"TriMatrixF.col2");
      Assert(lfx.col(j+1,i+1,N)(1) == T(k),"const TriMatrixF.col2");
      Assert(lfcv.col(j+1,i+1,N)(1) == T(k),"TriMatrixF.col2 CV");
      Assert(lfv.col(j+1,i+1,N)(1) == T(k),"TriMatrixF.col2 V");
      Assert(l.row(i,0,i)(j) == T(k),"TriMatrix.row1");
      Assert(lx.row(i,0,i)(j) == T(k),"const TriMatrix.row1");
      Assert(lcv.row(i,0,i)(j) == T(k),"TriMatrix.row1 CV");
      Assert(lv.row(i,0,i)(j) == T(k),"TriMatrix.row1 V");
      Assert(lf.row(i+1,1,i)(j+1) == T(k),"TriMatrixF.row1");
      Assert(lfx.row(i+1,1,i)(j+1) == T(k),"const TriMatrixF.row1");
      Assert(lfcv.row(i+1,1,i)(j+1) == T(k),"TriMatrixF.row1 CV");
      Assert(lfv.row(i+1,1,i)(j+1) == T(k),"TriMatrixF.row1 V");
      Assert(l.row(i,j,i)(0) == T(k),"TriMatrix.row2");
      Assert(lx.row(i,j,i)(0) == T(k),"const TriMatrix.row2");
      Assert(lcv.row(i,j,i)(0) == T(k),"TriMatrix.row2 CV");
      Assert(lv.row(i,j,i)(0) == T(k),"TriMatrix.row2 V");
      Assert(lf.row(i+1,j+1,i)(1) == T(k),"TriMatrixF.row2");
      Assert(lfx.row(i+1,j+1,i)(1) == T(k),"const TriMatrixF.row2");
      Assert(lfcv.row(i+1,j+1,i)(1) == T(k),"TriMatrixF.row2 CV");
      Assert(lfv.row(i+1,j+1,i)(1) == T(k),"TriMatrixF.row2 V");
      Assert(l.diag(-int(i-j))(j) == T(k),"TriMatrix.diag1");
      Assert(lx.diag(-int(i-j))(j) == T(k),"const TriMatrix.diag1");
      Assert(lcv.diag(-int(i-j))(j) == T(k),"TriMatrix.diag1 CV ");
      Assert(lv.diag(-int(i-j))(j) == T(k),"TriMatrix.diag1 V");
      Assert(lf.diag(-int(i-j))(j+1) == T(k),"TriMatrixF.diag1");
      Assert(lfx.diag(-int(i-j))(j+1) == T(k),"const TriMatrixF.diag1");
      Assert(lfcv.diag(-int(i-j))(j+1) == T(k),"TriMatrixF.diag1 CV ");
      Assert(lfv.diag(-int(i-j))(j+1) == T(k),"TriMatrixF.diag1 V");
      Assert(l.diag(-int(i-j),j,N+j-i)(0) == T(k),"TriMatrix.diag2");
      Assert(lx.diag(-int(i-j),j,N+j-i)(0) == T(k),"const TriMatrix.diag2");
      Assert(lcv.diag(-int(i-j),j,N+j-i)(0) == T(k),"TriMatrix.diag2 CV ");
      Assert(lv.diag(-int(i-j),j,N+j-i)(0) == T(k),"TriMatrix.diag2 V");
      Assert(lf.diag(-int(i-j),j+1,N+j-i)(1) == T(k),"TriMatrixF.diag2");
      Assert(lfx.diag(-int(i-j),j+1,N+j-i)(1) == T(k),"const TriMatrixF.diag2");
      Assert(lfcv.diag(-int(i-j),j+1,N+j-i)(1) == T(k),"TriMatrixF.diag2 CV ");
      Assert(lfv.diag(-int(i-j),j+1,N+j-i)(1) == T(k),"TriMatrixF.diag2 V");
    }
  }

  Matrix<T> mu(u);
  Matrix<T> ml(l);

  for(int i=0,k=0;i<N;++i) for(int j=0;j<N;++j,++k) {
    Assert(mu(i,j) == ux(i,j),"TriMatrix -> Matrix");
    Assert(ml(i,j) == lx(i,j),"TriMatrix -> Matrix");
  }
  Assert(u == UpperTriMatrixViewOf(mu),"TriMatrix == ");
  Assert(u == UpperTriMatrix<T,D,S>(mu),"TriMatrix == ");
  Assert(l == LowerTriMatrixViewOf(ml),"TriMatrix == ");
  Assert(l == LowerTriMatrix<T,D,S>(ml),"TriMatrix == ");

  Matrix<T> a(N,N);
  for(int i=0;i<N;++i) for(int j=0;j<N;++j) a(i,j) = 12+3*i-5*j;

  UpperTriMatrix<T,D,S> u1(a);
  UpperTriMatrix<T,D,S> u2(a.Transpose());
  LowerTriMatrix<T,D,S> l1(a);
  LowerTriMatrix<T,D,S> l2(a.Transpose());
  const UpperTriMatrix<T,D,S>& u1x = u1;
  const UpperTriMatrix<T,D,S>& u2x = u2;
  const LowerTriMatrix<T,D,S>& l1x = l1;
  const LowerTriMatrix<T,D,S>& l2x = l2;

  UpperTriMatrix<T> u4 = u1+u1;
  UpperTriMatrix<T> u5 = u1+u2;
  UpperTriMatrix<T> u6 = u2+u2;
  LowerTriMatrix<T> l4 = l1+l1;
  LowerTriMatrix<T> l5 = l1+l2;
  LowerTriMatrix<T> l6 = l2+l2;
  Matrix<T> m1 = l1+u1;
  Matrix<T> m2 = l1+u2;
  Matrix<T> m3 = l2+u1;
  Matrix<T> m4 = l2+u2;

  const UpperTriMatrix<T>& u4x = u4;
  const UpperTriMatrix<T>& u5x = u5;
  const UpperTriMatrix<T>& u6x = u6;
  const LowerTriMatrix<T>& l4x = l4;
  const LowerTriMatrix<T>& l5x = l5;
  const LowerTriMatrix<T>& l6x = l6;

  for(int i=0;i<N;i++) for(int j=0;j<N;j++) {
    Assert(u4x(i,j) == 2*u1x(i,j),"Add triMatrices");
    Assert(u5x(i,j) == u1x(i,j) + u2x(i,j),"Add triMatrices");
    Assert(u6x(i,j) == 2*u2x(i,j),"Add triMatrices");
    Assert(l4x(i,j) == 2*l1x(i,j),"Add triMatrices");
    Assert(l5x(i,j) == l1x(i,j) + l2x(i,j),"Add triMatrices");
    Assert(l6x(i,j) == 2*l2x(i,j),"Add triMatrices");
    Assert(m1(i,j) == l1x(i,j) + u1x(i,j),"Add triMatrices");
    Assert(m2(i,j) == l1x(i,j) + u2x(i,j),"Add triMatrices");
    Assert(m3(i,j) == l2x(i,j) + u1x(i,j),"Add triMatrices");
    Assert(m4(i,j) == l2x(i,j) + u2x(i,j),"Add triMatrices");
  }

  u4 = u1-u1;
  u5 = u1-u2;
  u6 = u2-u2;
  l4 = l1-l1;
  l5 = l1-l2;
  l6 = l2-l2;
  m1 = l1-u1;
  m2 = l1-u2;
  m3 = l2-u1;
  m4 = l2-u2;

  for(int i=0;i<N;i++) for(int j=0;j<N;j++) {
    Assert(u4x(i,j) == T(0),"Subtract TriMatrices");
    Assert(u5x(i,j) == u1x(i,j) - u2x(i,j),"Subtract TriMatrices");
    Assert(u6x(i,j) == T(0),"Subtract TriMatrices");
    Assert(l4x(i,j) == T(0),"Subtract TriMatrices");
    Assert(l5x(i,j) == l1x(i,j) - l2x(i,j),"Subtract TriMatrices");
    Assert(l6x(i,j) == T(0),"Subtract TriMatrices");
    Assert(m1(i,j) == l1x(i,j) - u1x(i,j),"Subtract TriMatrices");
    Assert(m2(i,j) == l1x(i,j) - u2x(i,j),"Subtract TriMatrices");
    Assert(m3(i,j) == l2x(i,j) - u1x(i,j),"Subtract TriMatrices");
    Assert(m4(i,j) == l2x(i,j) - u2x(i,j),"Subtract TriMatrices");
  }

  UpperTriMatrix<complex<T>,tmv::NonUnitDiag,tmv::RowMajor> cunr = u;
  cunr *= complex<T>(1,2);
  UpperTriMatrix<complex<T>,D,S> cu = cunr;
  Assert(cunr == cu,"TriMatrix == TriMatrix<N,R>");
  Assert(cunr.View() == cu.View(),"TriMatrix View");
  Assert(cunr.Transpose() == cu.Transpose(),"TriMatrix Transpose");
  Assert(cunr.Conjugate() == cu.Conjugate(),"TriMatrix Conjugate");
  Assert(cunr.Adjoint() == cu.Adjoint(),"TriMatrix Adjoint");
  Assert(cunr.OffDiag() == cu.OffDiag(),"TriMatrix OffDiag");
  Assert(cunr.Real() == cu.Real(),"TriMatrix Real");
  Assert(cunr.Imag() == cu.Imag(),"TriMatrix Imag");
  Assert(cunr.SubMatrix(0,N/2,N/2,N) == cu.SubMatrix(0,N/2,N/2,N),"TriMatrix SubMatrix");
  Assert(cunr.SubTriMatrix(0,N/2) == cu.SubTriMatrix(0,N/2),"TriMatrix SubTriMatrix");
  Assert(cunr.SubTriMatrix(N/2,N) == cu.SubTriMatrix(N/2,N),"TriMatrix SubTriMatrix");
  
  LowerTriMatrix<complex<T>,tmv::NonUnitDiag,tmv::RowMajor> clnr = l;
  clnr *= complex<T>(1,2);
  LowerTriMatrix<complex<T>,D,S> cl = clnr;
  Assert(clnr == cl,"TriMatrix == TriMatrix<N,R>");
  Assert(clnr.View() == cl.View(),"TriMatrix View");
  Assert(clnr.Transpose() == cl.Transpose(),"TriMatrix Transpose");
  Assert(clnr.Conjugate() == cl.Conjugate(),"TriMatrix Conjugate");
  Assert(clnr.Adjoint() == cl.Adjoint(),"TriMatrix Adjoint");
  Assert(clnr.OffDiag() == cl.OffDiag(),"TriMatrix OffDiag");
  Assert(clnr.Real() == cl.Real(),"TriMatrix Real");
  Assert(clnr.Imag() == cl.Imag(),"TriMatrix Imag");
  Assert(clnr.SubMatrix(N/2,N,0,N/2) == cl.SubMatrix(N/2,N,0,N/2),"TriMatrix SubMatrix");
  Assert(clnr.SubTriMatrix(0,N/2) == cl.SubTriMatrix(0,N/2),"TriMatrix SubTriMatrix");
  Assert(clnr.SubTriMatrix(N/2,N) == cl.SubTriMatrix(N/2,N),"TriMatrix SubTriMatrix");
  
}

template <class T> void TestTriMatrix()
{
  TestBasicTriMatrix<T,tmv::UnitDiag,tmv::RowMajor>();
  TestBasicTriMatrix<T,tmv::NonUnitDiag,tmv::RowMajor>();
  TestBasicTriMatrix<T,tmv::UnitDiag,tmv::ColMajor>();
  TestBasicTriMatrix<T,tmv::NonUnitDiag,tmv::ColMajor>();

  cout<<"TriMatrix<"<<tmv::Type(T())<<"> passed all basic tests\n";

  TestTriMatrixArith_A<T>();
  cout<<"TriMatrix<"<<tmv::Type(T())<<"> (Tri/Tri) Arithmetic passed all tests\n";
  TestTriMatrixArith_B<T>();
  cout<<"TriMatrix<"<<tmv::Type(T())<<"> (Matrix/Tri) Arithmetic passed all tests\n";
}

template void TestTriMatrix<double>();
#ifndef NOFLOAT
template void TestTriMatrix<float>();
#endif
#ifdef LONGDOUBLE
template void TestTriMatrix<long double>();
#endif
