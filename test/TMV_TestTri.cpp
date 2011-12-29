
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"
#include <fstream>
#include <cstdio>

template <class T, tmv::DiagType D, tmv::StorageType S> 
static void TestBasicUpperTriMatrix_1()
{
    const int N = 10;

    tmv::UpperTriMatrix<T,D,S> u(N);

    Assert(u.colsize() == size_t(N) && u.rowsize() == size_t(N),
           "Creating UpperTriMatrix(N)");

    for (int i=0,k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if (i < j || (D==tmv::NonUnitDiag && i==j)) {
            u(i,j) = T(k);
        }
    }

    tmv::UpperTriMatrixView<T> uv = u.view();
    tmv::ConstUpperTriMatrixView<T> ucv = u.view();

    const tmv::UpperTriMatrix<T,D,S>& ux = u;

    for(int i=0,k=1;i<N;++i) for(int j=0;j<N;++j,++k) {
        if (i < j) {
            Assert(u(i,j) == T(k),"Read/Write TriMatrix");
            Assert(ux(i,j) == T(k),"Access const TriMatrix");
            Assert(ucv(i,j) == T(k),"Access TriMatrix CV");
            Assert(uv(i,j) == T(k),"Access TriMatrix V");
            Assert(u.row(i,i+1,N)(j-i-1) == T(k),"TriMatrix.row1");
            Assert(ux.row(i,i+1,N)(j-i-1) == T(k),"const TriMatrix.row1");
            Assert(ucv.row(i,i+1,N)(j-i-1) == T(k),"TriMatrix.row1 CV");
            Assert(uv.row(i,i+1,N)(j-i-1) == T(k),"TriMatrix.row1 V");
            Assert(u.row(i,j,N)(0) == T(k),"TriMatrix.row2");
            Assert(ux.row(i,j,N)(0) == T(k),"const TriMatrix.row2");
            Assert(ucv.row(i,j,N)(0) == T(k),"TriMatrix.row2 CV");
            Assert(uv.row(i,j,N)(0) == T(k),"TriMatrix.row2 V");
            Assert(u.col(j,0,j)(i) == T(k),"TriMatrix.col1");
            Assert(ux.col(j,0,j)(i) == T(k),"const TriMatrix.col1");
            Assert(ucv.col(j,0,j)(i) == T(k),"TriMatrix.col1 CV");
            Assert(uv.col(j,0,j)(i) == T(k),"TriMatrix.col1 V");
            Assert(u.col(j,i,j)(0) == T(k),"TriMatrix.col2");
            Assert(ux.col(j,i,j)(0) == T(k),"const TriMatrix.col2");
            Assert(ucv.col(j,i,j)(0) == T(k),"TriMatrix.col2 CV");
            Assert(uv.col(j,i,j)(0) == T(k),"TriMatrix.col2 V");
            Assert(u.diag(j-i)(i) == T(k),"TriMatrix.diag1");
            Assert(ux.diag(j-i)(i) == T(k),"const TriMatrix.diag1");
            Assert(ucv.diag(j-i)(i) == T(k),"TriMatrix.diag1 CV ");
            Assert(uv.diag(j-i)(i) == T(k),"TriMatrix.diag1 V");
            Assert(u.diag(j-i,i,N-j+i)(0) == T(k),"TriMatrix.diag2");
            Assert(ux.diag(j-i,i,N-j+i)(0) == T(k),"const TriMatrix.diag2");
            Assert(ucv.diag(j-i,i,N-j+i)(0) == T(k),"TriMatrix.diag2 CV ");
            Assert(uv.diag(j-i,i,N-j+i)(0) == T(k),"TriMatrix.diag2 V");
        } else if (i==j) {
            if (D == tmv::UnitDiag) {
                Assert(ux(i,i) == T(1),"Access const TriMatrix");
                Assert(ucv(i,i) == T(1),"Access TriMatrix CV");
            } else {
                Assert(u(i,i) == T(k),"Read/Write TriMatrix");
                Assert(ux(i,i) == T(k),"Access const TriMatrix");
                Assert(ucv(i,i) == T(k),"Access TriMatrix CV");
                Assert(uv(i,i) == T(k),"Access TriMatrix V");
                Assert(u.row(i,i,N)(0) == T(k),"TriMatrix.row1");
                Assert(ux.row(i,i,N)(0) == T(k),"const TriMatrix.row1");
                Assert(ucv.row(i,i,N)(0) == T(k),"TriMatrix.row1 CV");
                Assert(uv.row(i,i,N)(0) == T(k),"TriMatrix.row1 V");
                Assert(u.col(i,0,i+1)(i) == T(k),"TriMatrix.col1");
                Assert(ux.col(i,0,i+1)(i) == T(k),"const TriMatrix.col1");
                Assert(ucv.col(i,0,i+1)(i) == T(k),"TriMatrix.col1 CV");
                Assert(uv.col(i,0,i+1)(i) == T(k),"TriMatrix.col1 V");
                Assert(u.col(i,i,i+1)(0) == T(k),"TriMatrix.col2");
                Assert(ux.col(i,i,i+1)(0) == T(k),"const TriMatrix.col2");
                Assert(ucv.col(i,i,i+1)(0) == T(k),"TriMatrix.col2 CV");
                Assert(uv.col(i,i,i+1)(0) == T(k),"TriMatrix.col2 V");
                Assert(u.diag()(i) == T(k),"TriMatrix.diag");
                Assert(ux.diag()(i) == T(k),"const TriMatrix.diag");
                Assert(ucv.diag()(i) == T(k),"TriMatrix.diag CV ");
                Assert(uv.diag()(i) == T(k),"TriMatrix.diag V");
                Assert(u.diag(0)(i) == T(k),"TriMatrix.diag1");
                Assert(ux.diag(0)(i) == T(k),"const TriMatrix.diag1");
                Assert(ucv.diag(0)(i) == T(k),"TriMatrix.diag1 CV ");
                Assert(uv.diag(0)(i) == T(k),"TriMatrix.diag1 V");
                Assert(u.diag(0,i,N-i+i)(0) == T(k),"TriMatrix.diag2");
                Assert(ux.diag(0,i,N-i+i)(0) == T(k),"const TriMatrix.diag2");
                Assert(ucv.diag(0,i,N-i+i)(0) == T(k),"TriMatrix.diag2 CV ");
                Assert(uv.diag(0,i,N-i+i)(0) == T(k),"TriMatrix.diag2 V");
            }
        } else {
            Assert(ux(i,j) == T(0),"Access const TriMatrix");
            Assert(ucv(i,j) == T(0),"Access TriMatrix CV");
        }
    }

#if (XTEST & 32)
    tmv::UpperTriMatrix<T,D,S,tmv::FortranStyle> uf(N);
    for (int i=0,k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if (i < j || (D==tmv::NonUnitDiag && i==j)) {
            uf(i+1,j+1) = T(k);
        }
    }

    tmv::UpperTriMatrixView<T,tmv::FortranStyle> ufv = uf.view();
    tmv::ConstUpperTriMatrixView<T,tmv::FortranStyle> ufcv = uf.view();

    const tmv::UpperTriMatrix<T,D,S,tmv::FortranStyle>& ufx = uf;

    for(int i=0,k=1;i<N;++i) for(int j=0;j<N;++j,++k) {
        if (i < j) {
            Assert(uf(i+1,j+1) == T(k),"Read/Write TriMatrixF");
            Assert(ufx(i+1,j+1) == T(k),"Access const TriMatrixF");
            Assert(ufcv(i+1,j+1) == T(k),"Access TriMatrixF CV");
            Assert(ufv(i+1,j+1) == T(k),"Access TriMatrixF V");
            Assert(uf.row(i+1,i+2,N)(j-i) == T(k),"TriMatrixF.row1");
            Assert(ufx.row(i+1,i+2,N)(j-i) == T(k),"const TriMatrixF.row1");
            Assert(ufcv.row(i+1,i+2,N)(j-i) == T(k),"TriMatrixF.row1 CV");
            Assert(ufv.row(i+1,i+2,N)(j-i) == T(k),"TriMatrixF.row1 V");
            Assert(uf.row(i+1,j+1,N)(1) == T(k),"TriMatrixF.row2");
            Assert(ufx.row(i+1,j+1,N)(1) == T(k),"const TriMatrixF.row2");
            Assert(ufcv.row(i+1,j+1,N)(1) == T(k),"TriMatrixF.row2 CV");
            Assert(ufv.row(i+1,j+1,N)(1) == T(k),"TriMatrixF.row2 V");
            Assert(uf.col(j+1,1,j)(i+1) == T(k),"TriMatrixF.col1");
            Assert(ufx.col(j+1,1,j)(i+1) == T(k),"const TriMatrixF.col1");
            Assert(ufcv.col(j+1,1,j)(i+1) == T(k),"TriMatrixF.col1 CV");
            Assert(ufv.col(j+1,1,j)(i+1) == T(k),"TriMatrixF.col1 V");
            Assert(uf.col(j+1,i+1,j)(1) == T(k),"TriMatrixF.col2");
            Assert(ufx.col(j+1,i+1,j)(1) == T(k),"const TriMatrixF.col2");
            Assert(ufcv.col(j+1,i+1,j)(1) == T(k),"TriMatrixF.col2 CV");
            Assert(ufv.col(j+1,i+1,j)(1) == T(k),"TriMatrixF.col2 V");
            Assert(uf.diag(j-i)(i+1) == T(k),"TriMatrixF.diag1");
            Assert(ufx.diag(j-i)(i+1) == T(k),"const TriMatrixF.diag1");
            Assert(ufcv.diag(j-i)(i+1) == T(k),"TriMatrixF.diag1 CV ");
            Assert(ufv.diag(j-i)(i+1) == T(k),"TriMatrixF.diag1 V");
            Assert(uf.diag(j-i,i+1,N-j+i)(1) == T(k),"TriMatrixF.diag2");
            Assert(ufx.diag(j-i,i+1,N-j+i)(1) == T(k),"const TriMatrixF.diag2");
            Assert(ufcv.diag(j-i,i+1,N-j+i)(1) == T(k),"TriMatrixF.diag2 CV ");
            Assert(ufv.diag(j-i,i+1,N-j+i)(1) == T(k),"TriMatrixF.diag2 V");
        } else if (i==j) {
            if (D == tmv::UnitDiag) {
                Assert(ufx(i+1,i+1) == T(1),"Access const TriMatrixF");
                Assert(ufcv(i+1,i+1) == T(1),"Access TriMatrixF CV");
            } else {
                Assert(uf(i+1,i+1) == T(k),"Read/Write TriMatrixF");
                Assert(ufx(i+1,i+1) == T(k),"Access const TriMatrixF");
                Assert(ufcv(i+1,i+1) == T(k),"Access TriMatrixF CV");
                Assert(ufv(i+1,i+1) == T(k),"Access TriMatrixF V");
                Assert(uf.row(i+1,i+1,N)(1) == T(k),"TriMatrixF.row1");
                Assert(ufx.row(i+1,i+1,N)(1) == T(k),"const TriMatrixF.row1");
                Assert(ufcv.row(i+1,i+1,N)(1) == T(k),"TriMatrixF.row1 CV");
                Assert(ufv.row(i+1,i+1,N)(1) == T(k),"TriMatrixF.row1 V");
                Assert(uf.col(i+1,1,i+1)(i+1) == T(k),"TriMatrixF.col1");
                Assert(ufx.col(i+1,1,i+1)(i+1) == T(k),"const TriMatrixF.col1");
                Assert(ufcv.col(i+1,1,i+1)(i+1) == T(k),"TriMatrixF.col1 CV");
                Assert(ufv.col(i+1,1,i+1)(i+1) == T(k),"TriMatrixF.col1 V");
                Assert(uf.col(i+1,i+1,i+1)(1) == T(k),"TriMatrixF.col2");
                Assert(ufx.col(i+1,i+1,i+1)(1) == T(k),"const TriMatrixF.col2");
                Assert(ufcv.col(i+1,i+1,i+1)(1) == T(k),"TriMatrixF.col2 CV");
                Assert(ufv.col(i+1,i+1,i+1)(1) == T(k),"TriMatrixF.col2 V");
                Assert(uf.diag()(i+1) == T(k),"TriMatrixF.diag");
                Assert(ufx.diag()(i+1) == T(k),"const TriMatrixF.diag");
                Assert(ufcv.diag()(i+1) == T(k),"TriMatrixF.diag CV ");
                Assert(ufv.diag()(i+1) == T(k),"TriMatrixF.diag V");
                Assert(uf.diag(0)(i+1) == T(k),"TriMatrixF.diag1");
                Assert(ufx.diag(0)(i+1) == T(k),"const TriMatrixF.diag1");
                Assert(ufcv.diag(0)(i+1) == T(k),"TriMatrixF.diag1 CV ");
                Assert(ufv.diag(0)(i+1) == T(k),"TriMatrixF.diag1 V");
                Assert(uf.diag(0,i+1,N-i+i)(1) == T(k),"TriMatrixF.diag2");
                Assert(ufx.diag(0,i+1,N-i+i)(1) == T(k),"const TriMatrixF.diag2");
                Assert(ufcv.diag(0,i+1,N-i+i)(1) == T(k),"TriMatrixF.diag2 CV ");
                Assert(ufv.diag(0,i+1,N-i+i)(1) == T(k),"TriMatrixF.diag2 V");
            }
        } else {
            Assert(ufx(i+1,j+1) == T(0),"Access const TriMatrixF");
            Assert(ufcv(i+1,j+1) == T(0),"Access TriMatrixF CV");
        }
    }
#endif

    u.resize(2);
    Assert(u.colsize() == 2 && u.rowsize() == 2,"u.resize(2)");
    for (int i=0, k=0; i<2; ++i) for (int j=0; j<2; ++j, ++k) 
        if (i < j || (D==tmv::NonUnitDiag && i==j)) u(i,j) = T(k);
    for (int i=0, k=0; i<2; ++i) for (int j=0; j<2; ++j, ++k) 
        if (i < j || (D==tmv::NonUnitDiag && i==j)) 
            Assert(u(i,j) == k,"Read/Write resized UpperTriMatrix");

    u.resize(2*N);
    Assert(u.colsize() == 2*N && u.rowsize() == 2*N,"m.resize(2*N)");
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<2*N; ++j, ++k) 
        if (i < j || (D==tmv::NonUnitDiag && i==j)) u(i,j) = T(k);
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<2*N; ++j, ++k) 
        if (i < j || (D==tmv::NonUnitDiag && i==j)) 
            Assert(u(i,j) == k,"Read/Write resized UpperTriMatrix");

}

template <class T, tmv::DiagType D, tmv::StorageType S> 
static void TestBasicLowerTriMatrix_1()
{
    const int N = 10;

    tmv::LowerTriMatrix<T,D,S> l(N);

    Assert(l.colsize() == size_t(N) && l.rowsize() == size_t(N),
           "Creating LowerTriMatrix(N)");

    for (int i=0,k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if (j < i || (D==tmv::NonUnitDiag && i==j)) {
            l(i,j) = T(k);
        }
    }

    tmv::LowerTriMatrixView<T> lv = l.view();
    tmv::ConstLowerTriMatrixView<T> lcv = l.view();

    const tmv::LowerTriMatrix<T,D,S>& lx = l;

    for(int i=0,k=1;i<N;++i) for(int j=0;j<N;++j,++k) {
        if (i < j) {
            Assert(lx(i,j) == T(0),"Access const TriMatrix");
            Assert(lcv(i,j) == T(0),"Access TriMatrix CV");
        } else if (i==j) {
            if (D == tmv::UnitDiag) {
                Assert(lx(i,i) == T(1),"Access const TriMatrix");
                Assert(lcv(i,i) == T(1),"Access TriMatrix CV");
            } else {
                Assert(l(i,i) == T(k),"Read/Write TriMatrix");
                Assert(lx(i,i) == T(k),"Access const TriMatrix");
                Assert(lcv(i,i) == T(k),"Access TriMatrix CV");
                Assert(lv(i,i) == T(k),"Access TriMatrix V");
                Assert(l.col(i,i,N)(0) == T(k),"TriMatrix.col2");
                Assert(lx.col(i,i,N)(0) == T(k),"const TriMatrix.col2");
                Assert(lcv.col(i,i,N)(0) == T(k),"TriMatrix.col2 CV");
                Assert(lv.col(i,i,N)(0) == T(k),"TriMatrix.col2 V");
                Assert(l.row(i,0,i+1)(i) == T(k),"TriMatrix.row1");
                Assert(lx.row(i,0,i+1)(i) == T(k),"const TriMatrix.row1");
                Assert(lcv.row(i,0,i+1)(i) == T(k),"TriMatrix.row1 CV");
                Assert(lv.row(i,0,i+1)(i) == T(k),"TriMatrix.row1 V");
                Assert(l.row(i,i,i+1)(0) == T(k),"TriMatrix.row2");
                Assert(lx.row(i,i,i+1)(0) == T(k),"const TriMatrix.row2");
                Assert(lcv.row(i,i,i+1)(0) == T(k),"TriMatrix.row2 CV");
                Assert(lv.row(i,i,i+1)(0) == T(k),"TriMatrix.row2 V");
                Assert(l.diag()(i) == T(k),"TriMatrix.diag");
                Assert(lx.diag()(i) == T(k),"const TriMatrix.diag");
                Assert(lcv.diag()(i) == T(k),"TriMatrix.diag CV ");
                Assert(lv.diag()(i) == T(k),"TriMatrix.diag V");
                Assert(l.diag(0)(i) == T(k),"TriMatrix.diag1");
                Assert(lx.diag(0)(i) == T(k),"const TriMatrix.diag1");
                Assert(lcv.diag(0)(i) == T(k),"TriMatrix.diag1 CV ");
                Assert(lv.diag(0)(i) == T(k),"TriMatrix.diag1 V");
                Assert(l.diag(0,i,N-i+i)(0) == T(k),"TriMatrix.diag2");
                Assert(lx.diag(0,i,N-i+i)(0) == T(k),"const TriMatrix.diag2");
                Assert(lcv.diag(0,i,N-i+i)(0) == T(k),"TriMatrix.diag2 CV ");
                Assert(lv.diag(0,i,N-i+i)(0) == T(k),"TriMatrix.diag2 V");
            }
        } else {
            Assert(l(i,j) == T(k),"Read/Write TriMatrix");
            Assert(lx(i,j) == T(k),"Access const TriMatrix");
            Assert(lcv(i,j) == T(k),"Access TriMatrix CV");
            Assert(lv(i,j) == T(k),"Access TriMatrix V");
            Assert(l.col(j,j+1,N)(i-j-1) == T(k),"TriMatrix.col1");
            Assert(lx.col(j,j+1,N)(i-j-1) == T(k),"const TriMatrix.col1");
            Assert(lcv.col(j,j+1,N)(i-j-1) == T(k),"TriMatrix.col1 CV");
            Assert(lv.col(j,j+1,N)(i-j-1) == T(k),"TriMatrix.col1 V");
            Assert(l.col(j,i,N)(0) == T(k),"TriMatrix.col2");
            Assert(lx.col(j,i,N)(0) == T(k),"const TriMatrix.col2");
            Assert(lcv.col(j,i,N)(0) == T(k),"TriMatrix.col2 CV");
            Assert(lv.col(j,i,N)(0) == T(k),"TriMatrix.col2 V");
            Assert(l.row(i,0,i)(j) == T(k),"TriMatrix.row1");
            Assert(lx.row(i,0,i)(j) == T(k),"const TriMatrix.row1");
            Assert(lcv.row(i,0,i)(j) == T(k),"TriMatrix.row1 CV");
            Assert(lv.row(i,0,i)(j) == T(k),"TriMatrix.row1 V");
            Assert(l.row(i,j,i)(0) == T(k),"TriMatrix.row2");
            Assert(lx.row(i,j,i)(0) == T(k),"const TriMatrix.row2");
            Assert(lcv.row(i,j,i)(0) == T(k),"TriMatrix.row2 CV");
            Assert(lv.row(i,j,i)(0) == T(k),"TriMatrix.row2 V");
            Assert(l.diag(-(i-j))(j) == T(k),"TriMatrix.diag1");
            Assert(lx.diag(-(i-j))(j) == T(k),"const TriMatrix.diag1");
            Assert(lcv.diag(-(i-j))(j) == T(k),"TriMatrix.diag1 CV ");
            Assert(lv.diag(-(i-j))(j) == T(k),"TriMatrix.diag1 V");
            Assert(l.diag(-(i-j),j,N+j-i)(0) == T(k),"TriMatrix.diag2");
            Assert(lx.diag(-(i-j),j,N+j-i)(0) == T(k),"const TriMatrix.diag2");
            Assert(lcv.diag(-(i-j),j,N+j-i)(0) == T(k),"TriMatrix.diag2 CV ");
            Assert(lv.diag(-(i-j),j,N+j-i)(0) == T(k),"TriMatrix.diag2 V");
        }
    }

#if (XTEST & 32)
    tmv::LowerTriMatrix<T,D,S,tmv::FortranStyle> lf(N);

    for (int i=0,k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if (j < i || (D==tmv::NonUnitDiag && i==j)) {
            lf(i+1,j+1) = T(k);
        }
    }

    tmv::LowerTriMatrixView<T,tmv::FortranStyle> lfv = lf.view();
    tmv::ConstLowerTriMatrixView<T,tmv::FortranStyle> lfcv = lf.view();

    const tmv::LowerTriMatrix<T,D,S,tmv::FortranStyle>& lfx = lf;

    for(int i=0,k=1;i<N;++i) for(int j=0;j<N;++j,++k) {
        if (i < j) {
            Assert(lfx(i+1,j+1) == T(0),"Access const TriMatrixF");
            Assert(lfcv(i+1,j+1) == T(0),"Access TriMatrixF CV");
        } else if (i==j) {
            if (D == tmv::UnitDiag) {
                Assert(lfx(i+1,i+1) == T(1),"Access const TriMatrixF");
                Assert(lfcv(i+1,i+1) == T(1),"Access TriMatrixF CV");
            } else {
                Assert(lf(i+1,i+1) == T(k),"Read/Write TriMatrixF");
                Assert(lfx(i+1,i+1) == T(k),"Access const TriMatrixF");
                Assert(lfcv(i+1,i+1) == T(k),"Access TriMatrixF CV");
                Assert(lfv(i+1,i+1) == T(k),"Access TriMatrixF V");
                Assert(lf.col(i+1,i+1,N)(1) == T(k),"TriMatrixF.col2");
                Assert(lfx.col(i+1,i+1,N)(1) == T(k),"const TriMatrixF.col2");
                Assert(lfcv.col(i+1,i+1,N)(1) == T(k),"TriMatrixF.col2 CV");
                Assert(lfv.col(i+1,i+1,N)(1) == T(k),"TriMatrixF.col2 V");
                Assert(lf.row(i+1,1,i+1)(i+1) == T(k),"TriMatrixF.row1");
                Assert(lfx.row(i+1,1,i+1)(i+1) == T(k),"const TriMatrixF.row1");
                Assert(lfcv.row(i+1,1,i+1)(i+1) == T(k),"TriMatrixF.row1 CV");
                Assert(lfv.row(i+1,1,i+1)(i+1) == T(k),"TriMatrixF.row1 V");
                Assert(lf.row(i+1,i+1,i+1)(1) == T(k),"TriMatrixF.row2");
                Assert(lfx.row(i+1,i+1,i+1)(1) == T(k),"const TriMatrixF.row2");
                Assert(lfcv.row(i+1,i+1,i+1)(1) == T(k),"TriMatrixF.row2 CV");
                Assert(lfv.row(i+1,i+1,i+1)(1) == T(k),"TriMatrixF.row2 V");
                Assert(lf.diag()(i+1) == T(k),"TriMatrixF.diag");
                Assert(lfx.diag()(i+1) == T(k),"const TriMatrixF.diag");
                Assert(lfcv.diag()(i+1) == T(k),"TriMatrixF.diag CV ");
                Assert(lfv.diag()(i+1) == T(k),"TriMatrixF.diag V");
                Assert(lf.diag(0)(i+1) == T(k),"TriMatrixF.diag1");
                Assert(lfx.diag(0)(i+1) == T(k),"const TriMatrixF.diag1");
                Assert(lfcv.diag(0)(i+1) == T(k),"TriMatrixF.diag1 CV ");
                Assert(lfv.diag(0)(i+1) == T(k),"TriMatrixF.diag1 V");
                Assert(lf.diag(0,i+1,N-i+i)(1) == T(k),"TriMatrixF.diag2");
                Assert(lfx.diag(0,i+1,N-i+i)(1) == T(k),"const TriMatrixF.diag2");
                Assert(lfcv.diag(0,i+1,N-i+i)(1) == T(k),"TriMatrixF.diag2 CV ");
                Assert(lfv.diag(0,i+1,N-i+i)(1) == T(k),"TriMatrixF.diag2 V");
            }
        } else {
            Assert(lf(i+1,j+1) == T(k),"Read/Write TriMatrixF");
            Assert(lfx(i+1,j+1) == T(k),"Access const TriMatrixF");
            Assert(lfcv(i+1,j+1) == T(k),"Access TriMatrixF CV");
            Assert(lfv(i+1,j+1) == T(k),"Access TriMatrixF V");
            Assert(lf.col(j+1,j+2,N)(i-j) == T(k),"TriMatrixF.col1");
            Assert(lfx.col(j+1,j+2,N)(i-j) == T(k),"const TriMatrixF.col1");
            Assert(lfcv.col(j+1,j+2,N)(i-j) == T(k),"TriMatrixF.col1 CV");
            Assert(lfv.col(j+1,j+2,N)(i-j) == T(k),"TriMatrixF.col1 V");
            Assert(lf.col(j+1,i+1,N)(1) == T(k),"TriMatrixF.col2");
            Assert(lfx.col(j+1,i+1,N)(1) == T(k),"const TriMatrixF.col2");
            Assert(lfcv.col(j+1,i+1,N)(1) == T(k),"TriMatrixF.col2 CV");
            Assert(lfv.col(j+1,i+1,N)(1) == T(k),"TriMatrixF.col2 V");
            Assert(lf.row(i+1,1,i)(j+1) == T(k),"TriMatrixF.row1");
            Assert(lfx.row(i+1,1,i)(j+1) == T(k),"const TriMatrixF.row1");
            Assert(lfcv.row(i+1,1,i)(j+1) == T(k),"TriMatrixF.row1 CV");
            Assert(lfv.row(i+1,1,i)(j+1) == T(k),"TriMatrixF.row1 V");
            Assert(lf.row(i+1,j+1,i)(1) == T(k),"TriMatrixF.row2");
            Assert(lfx.row(i+1,j+1,i)(1) == T(k),"const TriMatrixF.row2");
            Assert(lfcv.row(i+1,j+1,i)(1) == T(k),"TriMatrixF.row2 CV");
            Assert(lfv.row(i+1,j+1,i)(1) == T(k),"TriMatrixF.row2 V");
            Assert(lf.diag(-(i-j))(j+1) == T(k),"TriMatrixF.diag1");
            Assert(lfx.diag(-(i-j))(j+1) == T(k),"const TriMatrixF.diag1");
            Assert(lfcv.diag(-(i-j))(j+1) == T(k),"TriMatrixF.diag1 CV ");
            Assert(lfv.diag(-(i-j))(j+1) == T(k),"TriMatrixF.diag1 V");
            Assert(lf.diag(-(i-j),j+1,N+j-i)(1) == T(k),"TriMatrixF.diag2");
            Assert(lfx.diag(-(i-j),j+1,N+j-i)(1) == T(k),"const TriMatrixF.diag2");
            Assert(lfcv.diag(-(i-j),j+1,N+j-i)(1) == T(k),"TriMatrixF.diag2 CV ");
            Assert(lfv.diag(-(i-j),j+1,N+j-i)(1) == T(k),"TriMatrixF.diag2 V");
        }
    }
#endif

    l.resize(2);
    Assert(l.colsize() == 2 && l.rowsize() == 2,"l.resize(2)");
    for (int i=0, k=0; i<2; ++i) for (int j=0; j<2; ++j, ++k) 
        if (j < i || (D==tmv::NonUnitDiag && i==j)) l(i,j) = T(k);
    for (int i=0, k=0; i<2; ++i) for (int j=0; j<2; ++j, ++k) 
        if (j < i || (D==tmv::NonUnitDiag && i==j)) 
            Assert(l(i,j) == k,"Read/Write resized UpperTriMatrix");

    l.resize(2*N);
    Assert(l.colsize() == 2*N && l.rowsize() == 2*N,"m.resize(2*N)");
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<2*N; ++j, ++k) 
        if (j < i || (D==tmv::NonUnitDiag && i==j)) l(i,j) = T(k);
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<2*N; ++j, ++k) 
        if (j < i || (D==tmv::NonUnitDiag && i==j)) 
            Assert(l(i,j) == k,"Read/Write resized UpperTriMatrix");

}

template <class T, tmv::DiagType D, tmv::StorageType S> 
static void TestBasicTriMatrix_2()
{
    typedef std::complex<T> CT;
    const int N = 10;

    tmv::UpperTriMatrix<T,D,S> u(N);
    tmv::LowerTriMatrix<T,D,S> l(N);
    const tmv::UpperTriMatrix<T,D,S>& ux = u;
    const tmv::LowerTriMatrix<T,D,S>& lx = l;

    Assert(u.colsize() == size_t(N) && u.rowsize() == size_t(N),
           "Creating UpperTriMatrix(N)");
    Assert(l.colsize() == size_t(N) && l.rowsize() == size_t(N),
           "Creating LowerTriMatrix(N)");

    for (int i=0,k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if (i < j || (D==tmv::NonUnitDiag && i==j)) {
            u(i,j) = T(k);
        }
        if (j < i || (D==tmv::NonUnitDiag && i==j)) {
            l(i,j) = T(k);
        }
    }

    tmv::Matrix<T> mu = u;
    tmv::Matrix<T> ml = l;

    for(int i=0,k=0;i<N;++i) for(int j=0;j<N;++j,++k) {
        Assert(mu(i,j) == ux(i,j),"TriMatrix -> Matrix");
        Assert(ml(i,j) == lx(i,j),"TriMatrix -> Matrix");
    }
    Assert(u == mu.upperTri(),"TriMatrix == ");
    tmv::UpperTriMatrix<T,D,S> umu(mu);
    Assert(u == umu,"TriMatrix == ");
    Assert(l == ml.lowerTri(),"TriMatrix == ");
    tmv::LowerTriMatrix<T,D,S> lml(ml);
    Assert(l == lml,"TriMatrix == ");


    // Test assignments and constructors from arrays
    const T quarrmnonunit[] = {
        T(0),  T(-1), T(-2),
               T(1),  T(0),
                      T(2) 
    };
    const T quarrmunit[] = {
        T(1),  T(-1), T(-2),
               T(1),  T(0),
                      T(1) 
    };
    const T* quarrm = (D == tmv::UnitDiag) ? quarrmunit : quarrmnonunit;

    const T qlarrmnonunit[] = {
        T(0),
        T(2),  T(1),
        T(4),  T(3),  T(2) 
    };
    const T qlarrmunit[] = {
        T(1), 
        T(2),  T(1),
        T(4),  T(3),  T(1) 
    };
    const T* qlarrm = (D == tmv::UnitDiag) ? qlarrmunit : qlarrmnonunit;

    const T quarcmnonunit[] = {
        T(0),
        T(-1), T(1),
        T(-2), T(0),  T(2) 
    };
    const T quarcmunit[] = {
        T(1), 
        T(-1), T(1),
        T(-2), T(0),  T(1) 
    };
    const T* quarcm = (D == tmv::UnitDiag) ? quarcmunit : quarcmnonunit;

    const T qlarcmnonunit[] = {
        T(0),  T(2),  T(4),
               T(1),  T(3),
                      T(2) 
    };
    const T qlarcmunit[] = {
        T(1),  T(2),  T(4),
               T(1),  T(3),
                      T(1) 
    };
    const T* qlarcm = (D == tmv::UnitDiag) ? qlarcmunit : qlarcmnonunit;

    std::vector<T> quvecrm(6);
    for(int i=0;i<6;i++) quvecrm[i] = quarrm[i];
    std::vector<T> qlvecrm(6);
    for(int i=0;i<6;i++) qlvecrm[i] = qlarrm[i];
    std::vector<T> quveccm(6);
    for(int i=0;i<6;i++) quveccm[i] = quarcm[i];
    std::vector<T> qlveccm(6);
    for(int i=0;i<6;i++) qlveccm[i] = qlarcm[i];

    tmv::UpperTriMatrix<T,D,S> qu1(3);
    std::copy(quarrm, quarrm+6, qu1.rowmajor_begin());
    tmv::UpperTriMatrix<T,D,S> qu2(3);
    std::copy(quarcm, quarcm+6, qu2.colmajor_begin());

    tmv::LowerTriMatrix<T,D,S> ql1(3);
    std::copy(qlarrm, qlarrm+6, ql1.rowmajor_begin());
    tmv::LowerTriMatrix<T,D,S> ql2(3);
    std::copy(qlarcm, qlarcm+6, ql2.colmajor_begin());

    tmv::UpperTriMatrix<T,D,S> qu3(3);
    std::copy(quvecrm.begin(), quvecrm.end(), qu3.rowmajor_begin());
    tmv::UpperTriMatrix<T,D,S> qu4(3);
    std::copy(quveccm.begin(), quveccm.end(), qu4.colmajor_begin());

    tmv::LowerTriMatrix<T,D,S> ql3(3);
    std::copy(qlvecrm.begin(), qlvecrm.end(), ql3.rowmajor_begin());
    tmv::LowerTriMatrix<T,D,S> ql4(3);
    std::copy(qlveccm.begin(), qlveccm.end(), ql4.colmajor_begin());

    tmv::UpperTriMatrix<T,D,S> qu5x(30);
    tmv::UpperTriMatrixView<T> qu5 = qu5x.subTriMatrix(3,18,5);
    std::copy(quvecrm.begin(), quvecrm.end(), qu5.rowmajor_begin());
    tmv::UpperTriMatrix<T,D,S> qu6x(30);
    tmv::UpperTriMatrixView<T> qu6 = qu6x.subTriMatrix(3,18,5);
    std::copy(quveccm.begin(), quveccm.end(), qu6.colmajor_begin());

    tmv::LowerTriMatrix<T,D,S> ql5x(30);
    tmv::LowerTriMatrixView<T> ql5 = ql5x.subTriMatrix(3,18,5);
    std::copy(qlvecrm.begin(), qlvecrm.end(), ql5.rowmajor_begin());
    tmv::LowerTriMatrix<T,D,S> ql6x(30);
    tmv::LowerTriMatrixView<T> ql6 = ql6x.subTriMatrix(3,18,5);
    std::copy(qlveccm.begin(), qlveccm.end(), ql6.colmajor_begin());

    // Assignment using op<< is always in rowmajor order.
    tmv::UpperTriMatrix<T,D,S> qu7(3);
    tmv::LowerTriMatrix<T,D,S> qu8t(3);
    tmv::UpperTriMatrixView<T> qu8 = qu8t.transpose();

    tmv::LowerTriMatrix<T,D,S> ql7(3);
    tmv::UpperTriMatrix<T,D,S> ql8t(3);
    tmv::LowerTriMatrixView<T> ql8 = ql8t.transpose();

    if (D == tmv::UnitDiag) {
        qu7 <<
            1, -1, -2, 
                1,  0,
                    1;
        qu8 <<
            1, -1, -2, 
                1,  0,
                    1;

        ql7 <<
            1,
            2,  1,
            4,  3,  1;
        ql8 <<
            1,
            2,  1,
            4,  3,  1;
    } else {
        qu7 <<
            0, -1, -2, 
                1,  0,
                    2;
        qu8 <<
            0, -1, -2, 
                1,  0,
                    2;

        ql7 <<
            0,
            2,  1,
            4,  3,  2;
        ql8 <<
            0,
            2,  1,
            4,  3,  2;
    }

    // Can also view memory directly 
    // (Diagonal elements in memory do not have to be 1 if D == UnitDiag.)
    T qarrmfull[] = {
        T(0), T(-1), T(-2),
        T(2), T(1), T(0),
        T(4), T(3), T(2),
    };
    T qarcmfull[] = {
        T(0), T(2),  T(4),
        T(-1), T(1), T(3),
        T(-2), T(0), T(2),
    };
    T* qarfull = (S == tmv::RowMajor) ? qarrmfull : qarcmfull;
    const int Si = (S == tmv::RowMajor) ? 3 : 1;
    const int Sj = (S == tmv::RowMajor) ? 1 : 3;
    const tmv::ConstUpperTriMatrixView<T> qu9 = 
        tmv::UpperTriMatrixViewOf(qarfull,3,S,D);
    const tmv::ConstUpperTriMatrixView<T> qu10 = 
        tmv::UpperTriMatrixViewOf(qarfull,3,Si,Sj,D);
    const tmv::ConstLowerTriMatrixView<T> ql9 = 
        tmv::LowerTriMatrixViewOf(qarfull,3,S,D);
    const tmv::ConstLowerTriMatrixView<T> ql10 = 
        tmv::LowerTriMatrixViewOf(qarfull,3,Si,Sj,D);

    if (showacc) {
        std::cout<<"qu1 = "<<qu1<<std::endl;
        std::cout<<"qu2 = "<<qu2<<std::endl;
        std::cout<<"qu3 = "<<qu3<<std::endl;
        std::cout<<"qu4 = "<<qu4<<std::endl;
        std::cout<<"qu5 = "<<qu5<<std::endl;
        std::cout<<"qu6 = "<<qu6<<std::endl;
        std::cout<<"qu7 = "<<qu7<<std::endl;
        std::cout<<"qu8 = "<<qu8<<std::endl;
        std::cout<<"qu9 = "<<qu9<<std::endl;
        std::cout<<"qu10 = "<<qu10<<std::endl;

        std::cout<<"ql1 = "<<ql1<<std::endl;
        std::cout<<"ql2 = "<<ql2<<std::endl;
        std::cout<<"ql3 = "<<ql3<<std::endl;
        std::cout<<"ql4 = "<<ql4<<std::endl;
        std::cout<<"ql5 = "<<ql5<<std::endl;
        std::cout<<"ql6 = "<<ql6<<std::endl;
        std::cout<<"ql7 = "<<ql7<<std::endl;
        std::cout<<"ql8 = "<<ql8<<std::endl;
        std::cout<<"ql9 = "<<ql9<<std::endl;
        std::cout<<"ql10 = "<<ql10<<std::endl;
    }

    for(int i=0;i<3;i++) for(int j=0;j<3;j++) {
        T val = (D == tmv::UnitDiag && i==j) ? T(1) : T(2*i-j);
        if (j>=i) {
            Assert(qu1(i,j) == val,"Create UpperTriMatrix from T* rm");
            Assert(qu2(i,j) == val,"Create UpperTriMatrix from T* cm");
            Assert(qu3(i,j) == val,"Create UpperTriMatrix from vector rm");
            Assert(qu4(i,j) == val,"Create UpperTriMatrix from vector cm");
            Assert(qu5(i,j) == val,"Create UpperTriMatrixView from vector rm");
            Assert(qu6(i,j) == val,"Create UpperTriMatrixView from vector cm");
            Assert(qu7(i,j) == val,"Create UpperTriMatrix from << list");
            Assert(qu8(i,j) == val,"Create UpperTriMatrixView from << list");
            Assert(qu9(i,j) == val,"Create UpperTriMatrixView of T* (S)");
            Assert(qu10(i,j) == val,"Create UpperTriMatrixView of T* (Si,Sj)");
        }
        if (i>=j) {
            Assert(ql1(i,j) == val,"Create UpperTriMatrix from T* rm");
            Assert(ql2(i,j) == val,"Create UpperTriMatrix from T* cm");
            Assert(ql3(i,j) == val,"Create UpperTriMatrix from vector rm");
            Assert(ql4(i,j) == val,"Create UpperTriMatrix from vector cm");
            Assert(ql5(i,j) == val,"Create UpperTriMatrixView from vector rm");
            Assert(ql6(i,j) == val,"Create UpperTriMatrixView from vector cm");
            Assert(ql7(i,j) == val,"Create UpperTriMatrix from << list");
            Assert(ql8(i,j) == val,"Create UpperTriMatrixView from << list");
            Assert(ql9(i,j) == val,"Create UpperTriMatrixView of T* (S)");
            Assert(ql10(i,j) == val,"Create UpperTriMatrixView of T* (Si,Sj)");
        }
    }

    // Test the span of the iteration (i.e. the validity of begin(), end())
    const tmv::UpperTriMatrix<T,D,S>& qu1_const = qu1;
    tmv::UpperTriMatrixView<T> qu1_view = qu1.view();
    tmv::ConstUpperTriMatrixView<T> qu1_constview = qu1_const.view();
    tmv::ConstUpperTriMatrixView<T> qu5_const = qu5;

    typename tmv::UpperTriMatrix<T,D,S>::rowmajor_iterator rmitu1 = 
        qu1.rowmajor_begin();
    typename tmv::UpperTriMatrix<T,D,S>::const_rowmajor_iterator rmitu2 =
        qu1_const.rowmajor_begin();
    typename tmv::UpperTriMatrixView<T>::rowmajor_iterator rmitu3 =
        qu1_view.rowmajor_begin();
    typename tmv::ConstUpperTriMatrixView<T>::const_rowmajor_iterator rmitu4 =
        qu1_constview.rowmajor_begin();
    typename tmv::UpperTriMatrixView<T>::rowmajor_iterator rmitu5 =
        qu5.rowmajor_begin();
    typename tmv::ConstUpperTriMatrixView<T>::const_rowmajor_iterator rmitu6 =
        qu5_const.rowmajor_begin();
    int i = 0;
    while (rmitu1 != qu1.rowmajor_end()) {
        Assert(*rmitu1++ == quarrm[i], "RowMajor iteration 1");
        Assert(*rmitu2++ == quarrm[i], "RowMajor iteration 2");
        Assert(*rmitu3++ == quarrm[i], "RowMajor iteration 3");
        Assert(*rmitu4++ == quarrm[i], "RowMajor iteration 4");
        Assert(*rmitu5++ == quarrm[i], "RowMajor iteration 5");
        Assert(*rmitu6++ == quarrm[i], "RowMajor iteration 6");
        ++i;
    }
    Assert(i == 6, "RowMajor iteration number of elements");
    Assert(rmitu2 == qu1_const.rowmajor_end(), "rmitu2 reaching end");
    Assert(rmitu3 == qu1_view.rowmajor_end(), "rmitu3 reaching end");
    Assert(rmitu4 == qu1_constview.rowmajor_end(), "rmitu4 reaching end");
    Assert(rmitu5 == qu5.rowmajor_end(), "rmitu5 reaching end");
    Assert(rmitu6 == qu5_const.rowmajor_end(), "rmitu6 reaching end");

    const tmv::LowerTriMatrix<T,D,S>& ql1_const = ql1;
    tmv::LowerTriMatrixView<T> ql1_view = ql1.view();
    tmv::ConstLowerTriMatrixView<T> ql1_constview = ql1_const.view();
    tmv::ConstLowerTriMatrixView<T> ql5_const = ql5;

    typename tmv::LowerTriMatrix<T,D,S>::rowmajor_iterator rmitl1 = 
        ql1.rowmajor_begin();
    typename tmv::LowerTriMatrix<T,D,S>::const_rowmajor_iterator rmitl2 =
        ql1_const.rowmajor_begin();
    typename tmv::LowerTriMatrixView<T>::rowmajor_iterator rmitl3 =
        ql1_view.rowmajor_begin();
    typename tmv::ConstLowerTriMatrixView<T>::const_rowmajor_iterator rmitl4 =
        ql1_constview.rowmajor_begin();
    typename tmv::LowerTriMatrixView<T>::rowmajor_iterator rmitl5 =
        ql5.rowmajor_begin();
    typename tmv::ConstLowerTriMatrixView<T>::const_rowmajor_iterator rmitl6 =
        ql5_const.rowmajor_begin();
    i = 0;
    while (rmitl1 != ql1.rowmajor_end()) {
        Assert(*rmitl1++ == qlarrm[i], "RowMajor iteration 1");
        Assert(*rmitl2++ == qlarrm[i], "RowMajor iteration 2");
        Assert(*rmitl3++ == qlarrm[i], "RowMajor iteration 3");
        Assert(*rmitl4++ == qlarrm[i], "RowMajor iteration 4");
        Assert(*rmitl5++ == qlarrm[i], "RowMajor iteration 5");
        Assert(*rmitl6++ == qlarrm[i], "RowMajor iteration 6");
        ++i;
    }
    Assert(i == 6, "RowMajor iteration number of elements");
    Assert(rmitl2 == ql1_const.rowmajor_end(), "rmitl2 reaching end");
    Assert(rmitl3 == ql1_view.rowmajor_end(), "rmitl3 reaching end");
    Assert(rmitl4 == ql1_constview.rowmajor_end(), "rmitl4 reaching end");
    Assert(rmitl5 == ql5.rowmajor_end(), "rmitl5 reaching end");
    Assert(rmitl6 == ql5_const.rowmajor_end(), "rmitl6 reaching end");

    typename tmv::UpperTriMatrix<T,D,S>::colmajor_iterator cmitu1 = 
        qu1.colmajor_begin();
    typename tmv::UpperTriMatrix<T,D,S>::const_colmajor_iterator cmitu2 =
        qu1_const.colmajor_begin();
    typename tmv::UpperTriMatrixView<T>::colmajor_iterator cmitu3 =
        qu1_view.colmajor_begin();
    typename tmv::ConstUpperTriMatrixView<T>::const_colmajor_iterator cmitu4 =
        qu1_constview.colmajor_begin();
    typename tmv::UpperTriMatrixView<T>::colmajor_iterator cmitu5 =
        qu5.colmajor_begin();
    typename tmv::ConstUpperTriMatrixView<T>::const_colmajor_iterator cmitu6 =
        qu5_const.colmajor_begin();
    i = 0;
    while (cmitu1 != qu1.colmajor_end()) {
        Assert(*cmitu1++ == quarcm[i], "ColMajor iteration 1");
        Assert(*cmitu2++ == quarcm[i], "ColMajor iteration 2");
        Assert(*cmitu3++ == quarcm[i], "ColMajor iteration 3");
        Assert(*cmitu4++ == quarcm[i], "ColMajor iteration 4");
        Assert(*cmitu5++ == quarcm[i], "ColMajor iteration 5");
        Assert(*cmitu6++ == quarcm[i], "ColMajor iteration 6");
        ++i;
    }
    Assert(i == 6, "ColMajor iteration number of elements");
    Assert(cmitu2 == qu1_const.colmajor_end(), "cmitu2 reaching end");
    Assert(cmitu3 == qu1_view.colmajor_end(), "cmitu3 reaching end");
    Assert(cmitu4 == qu1_constview.colmajor_end(), "cmitu4 reaching end");
    Assert(cmitu5 == qu5.colmajor_end(), "cmitu5 reaching end");
    Assert(cmitu6 == qu5_const.colmajor_end(), "cmitu6 reaching end");

    typename tmv::LowerTriMatrix<T,D,S>::colmajor_iterator cmitl1 = 
        ql1.colmajor_begin();
    typename tmv::LowerTriMatrix<T,D,S>::const_colmajor_iterator cmitl2 =
        ql1_const.colmajor_begin();
    typename tmv::LowerTriMatrixView<T>::colmajor_iterator cmitl3 =
        ql1_view.colmajor_begin();
    typename tmv::ConstLowerTriMatrixView<T>::const_colmajor_iterator cmitl4 =
        ql1_constview.colmajor_begin();
    typename tmv::LowerTriMatrixView<T>::colmajor_iterator cmitl5 =
        ql5.colmajor_begin();
    typename tmv::ConstLowerTriMatrixView<T>::const_colmajor_iterator cmitl6 =
        ql5_const.colmajor_begin();
    i = 0;
    while (cmitl1 != ql1.colmajor_end()) {
        Assert(*cmitl1++ == qlarcm[i], "ColMajor iteration 1");
        Assert(*cmitl2++ == qlarcm[i], "ColMajor iteration 2");
        Assert(*cmitl3++ == qlarcm[i], "ColMajor iteration 3");
        Assert(*cmitl4++ == qlarcm[i], "ColMajor iteration 4");
        Assert(*cmitl5++ == qlarcm[i], "ColMajor iteration 5");
        Assert(*cmitl6++ == qlarcm[i], "ColMajor iteration 6");
        ++i;
    }
    Assert(i == 6, "ColMajor iteration number of elements");
    Assert(cmitl2 == ql1_const.colmajor_end(), "cmitl2 reaching end");
    Assert(cmitl3 == ql1_view.colmajor_end(), "cmitl3 reaching end");
    Assert(cmitl4 == ql1_constview.colmajor_end(), "cmitl4 reaching end");
    Assert(cmitl5 == ql5.colmajor_end(), "cmitl5 reaching end");
    Assert(cmitl6 == ql5_const.colmajor_end(), "cmitl6 reaching end");


    // Test Basic Arithmetic
    tmv::Matrix<T> a(N,N);
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) a(i,j) = T(12+3*i-5*j);

    tmv::UpperTriMatrix<T,D,S> u1(a);
    tmv::UpperTriMatrix<T,D,S> u2(a.transpose());
    tmv::LowerTriMatrix<T,D,S> l1(a);
    tmv::LowerTriMatrix<T,D,S> l2(a.transpose());
    const tmv::UpperTriMatrix<T,D,S>& u1x = u1;
    const tmv::UpperTriMatrix<T,D,S>& u2x = u2;
    const tmv::LowerTriMatrix<T,D,S>& l1x = l1;
    const tmv::LowerTriMatrix<T,D,S>& l2x = l2;

    tmv::UpperTriMatrix<T> u4 = u1+u1;
    tmv::UpperTriMatrix<T> u5 = u1+u2;
    tmv::UpperTriMatrix<T> u6 = u2+u2;
    tmv::LowerTriMatrix<T> l4 = l1+l1;
    tmv::LowerTriMatrix<T> l5 = l1+l2;
    tmv::LowerTriMatrix<T> l6 = l2+l2;
    tmv::Matrix<T> m1 = l1+u1;
    tmv::Matrix<T> m2 = l1+u2;
    tmv::Matrix<T> m3 = l2+u1;
    tmv::Matrix<T> m4 = l2+u2;

    const tmv::UpperTriMatrix<T>& u4x = u4;
    const tmv::UpperTriMatrix<T>& u5x = u5;
    const tmv::UpperTriMatrix<T>& u6x = u6;
    const tmv::LowerTriMatrix<T>& l4x = l4;
    const tmv::LowerTriMatrix<T>& l5x = l5;
    const tmv::LowerTriMatrix<T>& l6x = l6;

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

    tmv::UpperTriMatrix<CT,tmv::NonUnitDiag,tmv::RowMajor> cunr = u;
    if (D == tmv::UnitDiag)
        cunr.offDiag() *= CT(1,2);
    else
        cunr *= CT(1,2);
    tmv::UpperTriMatrix<CT,D,S> cu(cunr);
    Assert(cunr == cu,"TriMatrix == TriMatrix<N,R>");
    Assert(cunr.view() == cu.view(),"TriMatrix View");
    Assert(cunr.transpose() == cu.transpose(),"TriMatrix transpose");
    Assert(cunr.conjugate() == cu.conjugate(),"TriMatrix conjugate");
    Assert(cunr.adjoint() == cu.adjoint(),"TriMatrix adjoint");
    Assert(cunr.offDiag() == cu.offDiag(),"TriMatrix offDiag");
    Assert(cunr.realPart() == cu.realPart(),"TriMatrix realPart");
    if (D == tmv::NonUnitDiag)
        Assert(cunr.imagPart() == cu.imagPart(),"TriMatrix imagPart");
    Assert(cunr.subMatrix(0,N/2,N/2,N) == cu.subMatrix(0,N/2,N/2,N),"TriMatrix subMatrix");
    Assert(cunr.subTriMatrix(0,N/2) == cu.subTriMatrix(0,N/2),"TriMatrix subTriMatrix");
    Assert(cunr.subTriMatrix(N/2,N) == cu.subTriMatrix(N/2,N),"TriMatrix subTriMatrix");

    tmv::LowerTriMatrix<CT,tmv::NonUnitDiag,tmv::RowMajor> clnr = l;
    if (D == tmv::UnitDiag)
        clnr.offDiag() *= CT(1,2);
    else
        clnr *= CT(1,2);
    tmv::LowerTriMatrix<CT,D,S> cl = clnr;
    Assert(clnr == cl,"TriMatrix == TriMatrix<N,R>");
    Assert(clnr.view() == cl.view(),"TriMatrix View");
    Assert(clnr.transpose() == cl.transpose(),"TriMatrix transpose");
    Assert(clnr.conjugate() == cl.conjugate(),"TriMatrix conjugate");
    Assert(clnr.adjoint() == cl.adjoint(),"TriMatrix adjoint");
    Assert(clnr.offDiag() == cl.offDiag(),"TriMatrix offDiag");
    Assert(clnr.realPart() == cl.realPart(),"TriMatrix realPart");
    if (D == tmv::NonUnitDiag)
        Assert(clnr.imagPart() == cl.imagPart(),"TriMatrix imagPart");
    Assert(clnr.subMatrix(N/2,N,0,N/2) == cl.subMatrix(N/2,N,0,N/2),"TriMatrix subMatrix");
    Assert(clnr.subTriMatrix(0,N/2) == cl.subTriMatrix(0,N/2),"TriMatrix subTriMatrix");
    Assert(clnr.subTriMatrix(N/2,N) == cl.subTriMatrix(N/2,N),"TriMatrix subTriMatrix");
}

template <class T, tmv::DiagType D, tmv::StorageType S> 
static void TestBasicTriMatrix_IO()
{
    typedef std::complex<T> CT;

    const int N = 10;

    if (showstartdone) {
        std::cout<<"Start TestBasicTriMatrix_IO\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"D = "<<tmv::TMV_Text(D)<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
    }

    tmv::UpperTriMatrix<T,D,S> u(N);
    tmv::UpperTriMatrix<CT,D,S> cu(N);
    tmv::LowerTriMatrix<T,D,S> l(N);
    tmv::LowerTriMatrix<CT,D,S> cl(N);

    for (int i=0,k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if (i < j || (D==tmv::NonUnitDiag && i==j)) {
            u(i,j) = T(k);
            cu(i,j) = CT(k,k+1000);
        }
        if (j < i || (D==tmv::NonUnitDiag && i==j)) {
            l(i,j) = T(k);
            cl(i,j) = CT(k,k+1000);
        }
    }
    u(1,3) = l(3,1) = T(1.e-30);
    cu(1,3) = cl(3,1) = CT(1.e-30,1.e-30);
    u(5,6) = l(6,5) = T(9.e-3);
    cu(5,6) = cl(6,5) = CT(9.e-3,9.e-3);
    u(4,7) = l(7,4) = T(0.123456789);
    cu(4,7) = cl(7,4) = CT(3.123456789,600.987654321);

    // First check clipping function...
    tmv::UpperTriMatrix<T> u2 = u;
    tmv::UpperTriMatrix<CT> cu2 = cu;
    tmv::LowerTriMatrix<T> l2 = l;
    tmv::LowerTriMatrix<CT> cl2 = cl;
    if (!std::numeric_limits<T>::is_integer) {
        u2.clip(1.e-2);
        cu2.clip(1.e-2);
        l2.clip(1.e-2);
        cl2.clip(1.e-2);
    }
    tmv::UpperTriMatrix<T> u3 = u;
    tmv::UpperTriMatrix<CT> cu3 = cu;
    tmv::LowerTriMatrix<T> l3 = l;
    tmv::LowerTriMatrix<CT> cl3 = cl;
    u3(1,3) = l3(3,1) = T(0);
    cu3(1,3) = cl3(3,1) = T(0);
    u3(5,6) = l3(6,5) = T(0); // Others shouldn't get clipped.
    Assert(u2 == u3,"UpperTriMatrix clip");
    Assert(cu2 == cu3,"Complex UpperTriMatrix clip");
    Assert(l2 == l3,"LowerTriMatrix clip");
    Assert(cl2 == cl3,"Complex LowerTriMatrix clip");

    // Write matrices with 4 different styles
    std::ofstream fout("tmvtest_trimatrix_io.dat");
    Assert(fout,"Couldn't open tmvtest_trimatrix_io.dat for output");
    fout << u << std::endl;
    fout << l << std::endl;
    fout << cu << std::endl;
    fout << cl << std::endl;
    fout << tmv::CompactIO() << u << std::endl;
    fout << tmv::CompactIO() << l << std::endl;
    fout << tmv::CompactIO() << cu << std::endl;
    fout << tmv::CompactIO() << cl << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << u << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << l << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << cu << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << cl << std::endl;
    tmv::IOStyle myStyle = 
        tmv::CompactIO().setThresh(1.e-2).setPrecision(4).
        markup("Start","[",",","]","---","Done");
    fout << myStyle << u << std::endl;
    fout << myStyle << l << std::endl;
    fout << myStyle << cu << std::endl;
    fout << myStyle << cl << std::endl;
    fout.close();

    // When using (the default) prec(6), these will be the values read in.
    u(4,7) = l(7,4) = T(0.123457);
    cu(4,7) = cl(7,4) = CT(3.12346,600.988);

    // When using prec(12), the full correct values will be read in.

    // When using prec(4), these will be the values read in.
    u3(4,7) = l3(7,4) = T(0.1235);
    if (std::numeric_limits<T>::is_integer) cu3(4,7) = cl3(7,4) = CT(3,600);
    else cu3(4,7) = cl3(7,4) = CT(3.123,601.0);

    // Read them back in
    tmv::UpperTriMatrix<T,D,tmv::RowMajor> xu1(N);
    tmv::LowerTriMatrix<T,D,tmv::RowMajor> xl1(N);
    tmv::UpperTriMatrix<CT,D,tmv::RowMajor> xcu1(N);
    tmv::LowerTriMatrix<CT,D,tmv::RowMajor> xcl1(N);
    std::ifstream fin("tmvtest_trimatrix_io.dat");
    Assert(fin,"Couldn't open tmvtest_trimatrix_io.dat for input");
    fin >> xu1 >> xl1 >> xcu1 >> xcl1;
    Assert(u == xu1,"UpperTriMatrix I/O check normal");
    Assert(l == xl1,"LowerTriMatrix I/O check normal");
    Assert(cu == xcu1,"CUpperTriMatrix I/O check normal");
    Assert(cl == xcl1,"CLowerTriMatrix I/O check normal");
    fin >> tmv::CompactIO() >> xu1 >> tmv::CompactIO() >> xl1;
    fin >> tmv::CompactIO() >> xcu1 >> tmv::CompactIO() >> xcl1;
    Assert(u == xu1,"UpperTriMatrix I/O check compact");
    Assert(l == xl1,"LowerTriMatrix I/O check compact");
    Assert(cu == xcu1,"CUpperTriMatrix I/O check compact");
    Assert(cl == xcl1,"CLowerTriMatrix I/O check compact");
    fin >> xu1.view() >> xl1.view() >> xcu1.view() >> xcl1.view();
    Assert(u2 == xu1,"UpperTriMatrix I/O check thresh");
    Assert(l2 == xl1,"LowerTriMatrix I/O check thresh");
    Assert(cu2 == xcu1,"CUpperTriMatrix I/O check thresh");
    Assert(cl2 == xcl1,"CLowerTriMatrix I/O check thresh");
    fin >> myStyle >> xu1.view() >> myStyle >> xl1.view();
    fin >> myStyle >> xcu1.view() >> myStyle >> xcl1.view();
    Assert(u3 == xu1,"UpperTriMatrix I/O check compact thresh & prec(4)");
    Assert(l3 == xl1,"LowerTriMatrix I/O check compact thresh & prec(4)");
    Assert(cu3 == xcu1,"CUpperTriMatrix I/O check compact thresh & prec(4)");
    Assert(cl3 == xcl1,"CLowerTriMatrix I/O check compact thresh & prec(4)");
    fin.close();

    // Repeat for column major
    tmv::UpperTriMatrix<T,D,tmv::ColMajor> xu2(N);
    tmv::LowerTriMatrix<T,D,tmv::ColMajor> xl2(N);
    tmv::UpperTriMatrix<CT,D,tmv::ColMajor> xcu2(N);
    tmv::LowerTriMatrix<CT,D,tmv::ColMajor> xcl2(N);
    fin.open("tmvtest_trimatrix_io.dat");
    Assert(fin,"Couldn't open tmvtest_trimatrix_io.dat for input");
    fin >> xu2.view() >> xl2.view() >> xcu2.view() >> xcl2.view();
    Assert(u == xu2,"UpperTriMatrix I/O check normal");
    Assert(l == xl2,"LowerTriMatrix I/O check normal");
    Assert(cu == xcu2,"CUpperTriMatrix I/O check normal");
    Assert(cl == xcl2,"CLowerTriMatrix I/O check normal");
    fin >> tmv::CompactIO() >> xu2.view() >> tmv::CompactIO() >> xl2.view();
    fin >> tmv::CompactIO() >> xcu2.view() >> tmv::CompactIO() >> xcl2.view();
    Assert(u == xu2,"UpperTriMatrix I/O check compact");
    Assert(l == xl2,"LowerTriMatrix I/O check compact");
    Assert(cu == xcu2,"CUpperTriMatrix I/O check compact");
    Assert(cl == xcl2,"CLowerTriMatrix I/O check compact");
    fin >> xu2 >> xl2 >> xcu2 >> xcl2;
    Assert(u2 == xu2,"UpperTriMatrix I/O check thresh");
    Assert(l2 == xl2,"LowerTriMatrix I/O check thresh");
    Assert(cu2 == xcu2,"CUpperTriMatrix I/O check thresh");
    Assert(cl2 == xcl2,"CLowerTriMatrix I/O check thresh");
    fin >> myStyle >> xu2 >> myStyle >> xl2;
    fin >> myStyle >> xcu2 >> myStyle >> xcl2;
    Assert(u3 == xu2,"UpperTriMatrix I/O check compact thresh & prec(4)");
    Assert(l3 == xl2,"LowerTriMatrix I/O check compact thresh & prec(4)");
    Assert(cu3 == xcu2,"CUpperTriMatrix I/O check compact thresh & prec(4)");
    Assert(cl3 == xcl2,"CLowerTriMatrix I/O check compact thresh & prec(4)");
    fin.close();

    // And repeat for matrices that need to be resized.
    // Also check switching the default IOStyle.
    tmv::CompactIO().makeDefault();
    tmv::UpperTriMatrix<T,D> zu1,zu2,zu3,zu4;
    tmv::LowerTriMatrix<T,D> zl1,zl2,zl3,zl4;
    tmv::UpperTriMatrix<CT,D> zcu1,zcu2,zcu3,zcu4;
    tmv::LowerTriMatrix<CT,D> zcl1,zcl2,zcl3,zcl4;
    fin.open("tmvtest_trimatrix_io.dat");
    Assert(fin,"Couldn't open tmvtest_trimatrix_io.dat for input");
    fin >> tmv::NormalIO() >> zu1 >> tmv::NormalIO() >> zl1;
    fin >> tmv::NormalIO() >> zcu1 >> tmv::NormalIO() >> zcl1;
    Assert(u == zu1,"UpperTriMatrix I/O check normal");
    Assert(l == zl1,"LowerTriMatrix I/O check normal");
    Assert(cu == zcu1,"CUpperTriMatrix I/O check normal");
    Assert(cl == zcl1,"CLowerTriMatrix I/O check normal");
    fin >> zu2 >> zl2 >> zcu2 >> zcl2;
    Assert(u == zu2,"UpperTriMatrix I/O check compact");
    Assert(l == zl2,"LowerTriMatrix I/O check compact");
    Assert(cu == zcu2,"CUpperTriMatrix I/O check compact");
    Assert(cl == zcl2,"CLowerTriMatrix I/O check compact");
    fin >> tmv::NormalIO() >> zu3 >> tmv::NormalIO() >> zl3;
    fin >> tmv::NormalIO() >> zcu3 >> tmv::NormalIO() >> zcl3;
    Assert(u2 == zu3,"UpperTriMatrix I/O check thresh");
    Assert(l2 == zl3,"LowerTriMatrix I/O check thresh");
    Assert(cu2 == zcu3,"CUpperTriMatrix I/O check thresh");
    Assert(cl2 == zcl3,"CLowerTriMatrix I/O check thresh");
    fin >> myStyle >> zu4 >> myStyle >> zl4;
    fin >> myStyle >> zcu4 >> myStyle >> zcl4;
    Assert(u3 == zu4,"UpperTriMatrix I/O check compact thresh & prec(4)");
    Assert(l3 == zl4,"LowerTriMatrix I/O check compact thresh & prec(4)");
    Assert(cu3 == zcu4,"CUpperTriMatrix I/O check compact thresh & prec(4)");
    Assert(cl3 == zcl4,"CLowerTriMatrix I/O check compact thresh & prec(4)");
    fin.close();

    // Also try reading with the opposite DiagType.
    // This will succeed if doing Unit -> NonUnit,
    // but it will throw an exception if doing NonUnit -> Unit.
    fin.open("tmvtest_trimatrix_io.dat");
    Assert(fin,"Couldn't open tmvtest_trimatrix_io.dat for input");
    if (D==tmv::UnitDiag) {
        tmv::UpperTriMatrix<T,tmv::NonUnitDiag> zu5,zu6,zu7,zu8;
        tmv::LowerTriMatrix<T,tmv::NonUnitDiag> zl5,zl6,zl7,zl8;
        tmv::UpperTriMatrix<CT,tmv::NonUnitDiag> zcu5,zcu6,zcu7,zcu8;
        tmv::LowerTriMatrix<CT,tmv::NonUnitDiag> zcl5,zcl6,zcl7,zcl8;
        fin >> tmv::NormalIO() >> zu5 >> tmv::NormalIO() >> zl5;
        fin >> tmv::NormalIO() >> zcu5 >> tmv::NormalIO() >> zcl5;
        Assert(u == zu5,"UpperTriMatrix I/O check normal Unit->NonUnit");
        Assert(l == zl5,"LowerTriMatrix I/O check normal Unit->NonUnit");
        Assert(cu == zcu5,"CUpperTriMatrix I/O check normal Unit->NonUnit");
        Assert(cl == zcl5,"CLowerTriMatrix I/O check normal Unit->NonUnit");
        fin >> zu6 >> zl6 >> zcu6 >> zcl6;
        Assert(u == zu6,"UpperTriMatrix I/O check compact Unit->NonUnit");
        Assert(l == zl6,"LowerTriMatrix I/O check compact Unit->NonUnit");
        Assert(cu == zcu6,"CUpperTriMatrix I/O check compact Unit->NonUnit");
        Assert(cl == zcl6,"CLowerTriMatrix I/O check compact Unit->NonUnit");
        fin >> tmv::NormalIO() >> zu7 >> tmv::NormalIO() >> zl7;
        fin >> tmv::NormalIO() >> zcu7 >> tmv::NormalIO() >> zcl7;
        Assert(u2 == zu7,"UpperTriMatrix I/O check thresh Unit->NonUnit");
        Assert(l2 == zl7,"LowerTriMatrix I/O check thresh Unit->NonUnit");
        Assert(cu2 == zcu7,"CUpperTriMatrix I/O check thresh Unit->NonUnit");
        Assert(cl2 == zcl7,"CLowerTriMatrix I/O check thresh Unit->NonUnit");
        fin >> myStyle >> zu8 >> myStyle >> zl8;
        fin >> myStyle >> zcu8 >> myStyle >> zcl8;
        Assert(u3 == zu8,"UpperTriMatrix I/O check compact thresh Unit->NonUnit");
        Assert(l3 == zl8,"LowerTriMatrix I/O check compact thresh Unit->NonUnit");
        Assert(cu3 == zcu8,"CUpperTriMatrix I/O check compact thresh Unit->NonUnit");
        Assert(cl3 == zcl8,"CLowerTriMatrix I/O check compact thresh Unit->NonUnit");
        fin.close();
    } else {
#ifndef NOTHROW
        tmv::UpperTriMatrix<T,tmv::UnitDiag> zu5;
        try {
            fin >> tmv::NormalIO() >> zu5;
            Assert(false,"Throw ReadError for UnitDiag read of NonUnitDiag");
        } catch (tmv::ReadError&) {
            Assert(true,"Catch ReadError for UnitDiag read of NonUnitDiag");
        }
        // stream is already corrupted at this point, so no point 
        // in continuing on with other matrices.
#endif
    }
    fin.close();
    tmv::IOStyle::revertDefault();

    // Finally, check that the NormalIO can be read in as a regular matrix.
    tmv::Matrix<T> zm1,zm2;
    tmv::Matrix<CT> zcm1,zcm2;
    fin.open("tmvtest_trimatrix_io.dat");
    Assert(fin,"Couldn't open tmvtest_trimatrix_io.dat for input");
    fin >> zm1 >> zm2 >> zcm1 >> zcm2;
    Assert(u == zm1,"UpperTriMatrix -> Matrix I/O check");
    Assert(l == zm2,"LowerTriMatrix -> Matrix I/O check");
    Assert(cu == zcm1,"CUpperTriMatrix -> CMatrix I/O check");
    Assert(cl == zcm2,"CLowerTriMatrix -> CMatrix I/O check");
    fin.close();

#if XTEST == 0
    std::remove("tmvtest_trimatrix_io.dat");
#endif
}

template <class T, tmv::DiagType D, tmv::StorageType S> 
static void TestBasicTriMatrix()
{
    TestBasicUpperTriMatrix_1<T,D,S>();
    TestBasicLowerTriMatrix_1<T,D,S>();
    TestBasicTriMatrix_2<T,D,S>();
    TestBasicTriMatrix_IO<T,D,S>();
}

template <class T> void TestTriMatrix()
{
#if 1
    TestBasicTriMatrix<T,tmv::UnitDiag,tmv::RowMajor>();
    TestBasicTriMatrix<T,tmv::NonUnitDiag,tmv::RowMajor>();
    TestBasicTriMatrix<T,tmv::UnitDiag,tmv::ColMajor>();
    TestBasicTriMatrix<T,tmv::NonUnitDiag,tmv::ColMajor>();
    std::cout<<"TriMatrix<"<<Text(T())<<"> passed all basic tests\n";
#endif

#if 1
    TestAllAliasMultUL<T>();
    std::cout<<"TriMatrix<"<<Text(T())<<"> passed all aliased multiplication tests\n";
#endif

#if 1
    TestTriMatrixArith_A1a<T>();
    TestTriMatrixArith_A1b<T>();
    TestTriMatrixArith_A2<T>();
    TestTriMatrixArith_A3<T>();
    TestTriMatrixArith_A4a<T>();
    TestTriMatrixArith_A4b<T>();
    TestTriMatrixArith_A5a<T>();
    TestTriMatrixArith_A5b<T>();
    TestTriMatrixArith_A6a<T>();
    TestTriMatrixArith_A6b<T>();
    TestTriMatrixArith_A6c<T>();
    std::cout<<"TriMatrix<"<<Text(T())<<"> (Tri/Tri) Arithmetic passed all tests\n";
#endif
#if 1
    TestTriMatrixArith_B4a<T>();
    TestTriMatrixArith_B4b<T>();
    TestTriMatrixArith_B5a<T>();
    TestTriMatrixArith_B5b<T>();
    TestTriMatrixArith_B6a<T>();
    TestTriMatrixArith_B6b<T>();
    std::cout<<"TriMatrix<"<<Text(T())<<"> (Matrix/Tri) Arithmetic passed all tests\n";
#endif
#if 1
    TestTriMatrixArith_C4a<T>();
    TestTriMatrixArith_C4b<T>();
    TestTriMatrixArith_C5a<T>();
    TestTriMatrixArith_C5b<T>();
    TestTriMatrixArith_C6a<T>();
    TestTriMatrixArith_C6b<T>();
    std::cout<<"TriMatrix<"<<Text(T())<<"> (Diag/Tri) Arithmetic passed all tests\n";
#endif
}

#ifdef TEST_DOUBLE
template void TestTriMatrix<double>();
#endif
#ifdef TEST_FLOAT
template void TestTriMatrix<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestTriMatrix<long double>();
#endif
#ifdef TEST_INT
template void TestTriMatrix<int>();
#endif
