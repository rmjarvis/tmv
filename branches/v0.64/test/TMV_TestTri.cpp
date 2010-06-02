#define XDEBUG_PRODMM
#define PRINTALGO_UL

#include "TMV_Test.h"
#include "TMV_Test1.h"
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
            Assert(l.diag(-int(i-j))(j) == T(k),"TriMatrix.diag1");
            Assert(lx.diag(-int(i-j))(j) == T(k),"const TriMatrix.diag1");
            Assert(lcv.diag(-int(i-j))(j) == T(k),"TriMatrix.diag1 CV ");
            Assert(lv.diag(-int(i-j))(j) == T(k),"TriMatrix.diag1 V");
            Assert(l.diag(-int(i-j),j,N+j-i)(0) == T(k),"TriMatrix.diag2");
            Assert(lx.diag(-int(i-j),j,N+j-i)(0) == T(k),"const TriMatrix.diag2");
            Assert(lcv.diag(-int(i-j),j,N+j-i)(0) == T(k),"TriMatrix.diag2 CV ");
            Assert(lv.diag(-int(i-j),j,N+j-i)(0) == T(k),"TriMatrix.diag2 V");
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
            Assert(lf.diag(-int(i-j))(j+1) == T(k),"TriMatrixF.diag1");
            Assert(lfx.diag(-int(i-j))(j+1) == T(k),"const TriMatrixF.diag1");
            Assert(lfcv.diag(-int(i-j))(j+1) == T(k),"TriMatrixF.diag1 CV ");
            Assert(lfv.diag(-int(i-j))(j+1) == T(k),"TriMatrixF.diag1 V");
            Assert(lf.diag(-int(i-j),j+1,N+j-i)(1) == T(k),"TriMatrixF.diag2");
            Assert(lfx.diag(-int(i-j),j+1,N+j-i)(1) == T(k),"const TriMatrixF.diag2");
            Assert(lfcv.diag(-int(i-j),j+1,N+j-i)(1) == T(k),"TriMatrixF.diag2 CV ");
            Assert(lfv.diag(-int(i-j),j+1,N+j-i)(1) == T(k),"TriMatrixF.diag2 V");
        }
    }
#endif
}

template <class T, tmv::DiagType D, tmv::StorageType S> 
static void TestBasicTriMatrix_2()
{
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

    std::vector<T> qv(9);
    if (S == tmv::RowMajor) {
        const T qvar[] = { T(0), T(-1), T(-2),
            T(2), T(1),  T(0),
            T(4), T(3),  T(2) };
        for(int i=0;i<9;i++) qv[i] = qvar[i];
    } else {
        const T qvar[] = { T(0),  T(2), T(4),
            T(-1), T(1), T(3),
            T(-2), T(0), T(2) };
        for(int i=0;i<9;i++) qv[i] = qvar[i];
    }
    T qar[9];
    for(int i=0;i<9;i++) qar[i] = qv[i];

    {
        const tmv::UpperTriMatrix<T,D,S> q1(3,qar);
        const tmv::UpperTriMatrix<T,D,S> q2(3,qv);
        for(int i=0;i<3;i++) for(int j=0;j<3;j++) if (j>=i) {
            if (D == tmv::UnitDiag && i==j) {
                Assert(q1(i,j) == T(1),"Create UpperTriMatrix from T*");
                Assert(q2(i,j) == T(1),"Create UpperTriMatrix from vector");
            } else {
                Assert(q1(i,j) == T(2*i-j),"Create UpperTriMatrix from T*");
                Assert(q2(i,j) == T(2*i-j),"Create UpperTriMatrix from vector");
            }
        }
        if (D == tmv::UnitDiag) {
            const tmv::ConstUpperTriMatrixView<T> q3 = (
                tmv::UnitUpperTriMatrixViewOf(qar,3,S));
            for(int i=0;i<3;i++) for(int j=0;j<3;j++) if (j>=i) {
                if (i == j)
                    Assert(q3(i,j) == T(1),"Create UnitUpperTriMatrixView of T*");
                else
                    Assert(q3(i,j) == T(2*i-j),"Create UnitUpperTriMatrixView of T*");
            }
        } else {
            const tmv::ConstUpperTriMatrixView<T> q3 = (
                tmv::UpperTriMatrixViewOf(qar,3,S));
            for(int i=0;i<3;i++) for(int j=0;j<3;j++) if (j>=i) {
                Assert(q3(i,j) == T(2*i-j),"Create UpperTriMatrixView of T*");
            }
        }
    }
    {
        const tmv::LowerTriMatrix<T,D,S> q1(3,qar);
        const tmv::LowerTriMatrix<T,D,S> q2(3,qv);
        for(int i=0;i<3;i++) for(int j=0;j<3;j++) if (i>=j) {
            if (D == tmv::UnitDiag && i==j) {
                Assert(q1(i,j) == T(1),"Create LowerTriMatrix from T*");
                Assert(q2(i,j) == T(1),"Create LowerTriMatrix from vector");
            } else {
                Assert(q1(i,j) == T(2*i-j),"Create LowerTriMatrix from T*");
                Assert(q2(i,j) == T(2*i-j),"Create LowerTriMatrix from vector");
            }
        }
        if (D == tmv::UnitDiag) {
            const tmv::ConstLowerTriMatrixView<T> q3 = (
                tmv::UnitLowerTriMatrixViewOf(qar,3,S));
            for(int i=0;i<3;i++) for(int j=0;j<3;j++) if (i>=j) {
                if (i == j)
                    Assert(q3(i,j) == T(1),"Create UnitLowerTriMatrixView of T*");
                else
                    Assert(q3(i,j) == T(2*i-j),"Create UnitLowerTriMatrixView of T*");
            }
        } else {
            const tmv::ConstLowerTriMatrixView<T> q3 = (
                tmv::LowerTriMatrixViewOf(qar,3,S));
            for(int i=0;i<3;i++) for(int j=0;j<3;j++) if (i>=j) {
                Assert(q3(i,j) == T(2*i-j),"Create LowerTriMatrixView of T*");
            }
        }
    }

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

    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::RowMajor> cunr = u;
    if (D == tmv::UnitDiag)
        cunr.offDiag() *= std::complex<T>(1,2);
    else
        cunr *= std::complex<T>(1,2);
    tmv::UpperTriMatrix<std::complex<T>,D,S> cu = cunr;
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

    tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::RowMajor> clnr = l;
    if (D == tmv::UnitDiag)
        clnr.offDiag() *= std::complex<T>(1,2);
    else
        clnr *= std::complex<T>(1,2);
    tmv::LowerTriMatrix<std::complex<T>,D,S> cl = clnr;
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
    const int N = 10;

    tmv::UpperTriMatrix<T,D,S> u(N);
    tmv::LowerTriMatrix<T,D,S> l(N);

    for (int i=0,k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if (i < j || (D==tmv::NonUnitDiag && i==j)) {
            u(i,j) = T(k);
        }
        if (j < i || (D==tmv::NonUnitDiag && i==j)) {
            l(i,j) = T(k);
        }
    }

    std::ofstream fout("tmvtest_trimatrix_io.dat");
    if (!fout) 
#ifdef NOTHROW
    { std::cerr<<"Couldn't open tmvtest_trimatrix_io.dat for output\n"; exit(1); }
#else
    throw std::runtime_error(
        "Couldn't open tmvtest_trimatrix_io.dat for output");
#endif
    fout << u << std::endl << l << std::endl;
    u.writeCompact(fout);
    l.writeCompact(fout);
    fout.close();

    tmv::Matrix<std::complex<T>,tmv::RowMajor> xum1(N,N);
    tmv::Matrix<std::complex<T>,tmv::RowMajor> xlm1(N,N);
    tmv::UpperTriMatrix<std::complex<T>,D,tmv::RowMajor> xu1(N);
    tmv::LowerTriMatrix<std::complex<T>,D,tmv::RowMajor> xl1(N);
    std::ifstream fin("tmvtest_trimatrix_io.dat");
    if (!fin) 
#ifdef NOTHROW
    { std::cerr<<"Couldn't open tmvtest_trimatrix_io.dat for input\n"; exit(1); }
#else
    throw std::runtime_error(
        "Couldn't open tmvtest_trimatrix_io.dat for input");
#endif
    fin >> xum1 >> xlm1;
    fin >> xu1 >> xl1;
    fin.close();
    Assert(tmv::Matrix<std::complex<T> >(u) == xum1,"UpperTriMatrix I/O check #1");
    Assert(tmv::Matrix<std::complex<T> >(l) == xlm1,"LowerTriMatrix I/O check #1");
    Assert(u == xu1,"UpperTriMatrix Compact I/O check #1");
    Assert(l == xl1,"LowerTriMatrix Compact I/O check #1");

    tmv::Matrix<std::complex<T>,tmv::ColMajor> xum2(N,N);
    tmv::Matrix<std::complex<T>,tmv::ColMajor> xlm2(N,N);
    tmv::UpperTriMatrix<std::complex<T>,D,tmv::ColMajor> xu2(N);
    tmv::LowerTriMatrix<std::complex<T>,D,tmv::ColMajor> xl2(N);
    fin.open("tmvtest_trimatrix_io.dat");
    if (!fin) 
#ifdef NOTHROW
    { std::cerr<<"Couldn't open tmvtest_trimatrix_io.dat\n"; exit(1); }
#else
    throw std::runtime_error(
        "Couldn't open tmvtest_trimatrix_io.dat for input");
#endif
    fin >> xum2 >> xlm2;
    fin >> xu2 >> xl2;
    fin.close();
    Assert(tmv::Matrix<std::complex<T> >(u) == xum2,"UpperTriMatrix I/O check #2");
    Assert(tmv::Matrix<std::complex<T> >(l) == xlm2,"LowerTriMatrix I/O check #2");
    Assert(u == xu2,"UpperTriMatrix Compact I/O check #2");
    Assert(l == xl2,"LowerTriMatrix Compact I/O check #2");

    std::auto_ptr<tmv::Matrix<std::complex<T> > > xum3;
    std::auto_ptr<tmv::Matrix<std::complex<T> > > xlm3;
    std::auto_ptr<tmv::UpperTriMatrix<std::complex<T> > > xu3;
    std::auto_ptr<tmv::LowerTriMatrix<std::complex<T> > > xl3;
    fin.open("tmvtest_trimatrix_io.dat");
    if (!fin) 
#ifdef NOTHROW
    { std::cerr<<"Couldn't open tmvtest_trimatrix_io.dat\n"; exit(1); }
#else
    throw std::runtime_error(
        "Couldn't open tmvtest_trimatrix_io.dat for input");
#endif
    fin >> xum3 >> xlm3;
    fin >> xu3 >> xl3;
    fin.close();
    Assert(tmv::Matrix<std::complex<T> >(u) == *xum3,
           "UpperTriMatrix I/O check #3");
    Assert(tmv::Matrix<std::complex<T> >(l) == *xlm3,
           "LowerTriMatrix I/O check #3");
    Assert(u == *xu3,"UpperTriMatrix Compact I/O check #3");
    Assert(l == *xl3,"LowerTriMatrix Compact I/O check #3");

    if (D==tmv::UnitDiag) {
        tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> xu4(N);
        tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag> xl4(N);
        fin.open("tmvtest_trimatrix_io.dat");
        if (!fin) 
#ifdef NOTHROW
        { std::cerr<<"Couldn't open tmvtest_trimatrix_io.dat\n"; exit(1); }
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_trimatrix_io.dat for input");
#endif
        fin >> xum1 >> xlm1;
        fin >> xu4 >> xl4;
        fin.close();
        Assert(u == xu2,"NonUnitDiag UpperTriMatrix Compact I/O check #4");
        Assert(l == xl2,"NonUnitDiag LowerTriMatrix Compact I/O check #4");
    } else {
        fin.open("tmvtest_trimatrix_io.dat");
        if (!fin) 
#ifdef NOTHROW
        { std::cerr<<"Couldn't open tmvtest_trimatrix_io.dat\n"; exit(1); }
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_trimatrix_io.dat for input");
#endif
        fin >> xum1 >> xlm1;
#ifndef NOTHROW
        try {
            tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::RowMajor> xu4(N);
            fin >> xu4;
            Assert(false,"Throw ReadError for UnitDiag read");
        } catch (tmv::ReadError) {
            Assert(true,"Catch ReadError for UnitDiag read");
        }
        try {
            tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::RowMajor> xl4(N);
            fin >> xl4;
            Assert(false,"Throw ReadError for UnitDiag read");
        } catch (tmv::ReadError) {
            Assert(true,"Catch ReadError for UnitDiag read");
        }
#endif
        fin.close();
    }

#ifndef XTEST
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
    std::cout<<"TriMatrix<"<<tmv::TMV_Text(T())<<"> passed all basic tests\n";
#endif

    if (tmv::TMV_Epsilon<T>() > T(0)) {
#if 1
        TestAllAliasMultUL<T>();
        std::cout<<"TriMatrix<"<<tmv::TMV_Text(T())<<"> passed all aliased multiplication tests\n";
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
        std::cout<<"TriMatrix<"<<tmv::TMV_Text(T())<<"> (Tri/Tri) Arithmetic passed all tests\n";
#endif
#if 1
        TestTriMatrixArith_B4a<T>();
        TestTriMatrixArith_B4b<T>();
        TestTriMatrixArith_B5a<T>();
        TestTriMatrixArith_B5b<T>();
        TestTriMatrixArith_B6a<T>();
        TestTriMatrixArith_B6b<T>();
        std::cout<<"TriMatrix<"<<tmv::TMV_Text(T())<<"> (Matrix/Tri) Arithmetic passed all tests\n";
#endif
#if 1
        TestTriMatrixArith_C4a<T>();
        TestTriMatrixArith_C4b<T>();
        TestTriMatrixArith_C5a<T>();
        TestTriMatrixArith_C5b<T>();
        TestTriMatrixArith_C6a<T>();
        TestTriMatrixArith_C6b<T>();
        std::cout<<"TriMatrix<"<<tmv::TMV_Text(T())<<"> (Diag/Tri) Arithmetic passed all tests\n";
#endif
    }
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
