
#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include <fstream>
#include <cstdio>

template <class T, tmv::UpLoType U, tmv::StorageType S>
static void TestBasicSymBandMatrix_1()
{
    const int N = 10;
    const int noff = 3;

    if (showstartdone) {
        std::cout<<"Start TestBasicSymBandMatrix_1\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"U = "<<tmv::TMV_Text(U)<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
        std::cout<<"noff = "<<noff<<std::endl;
    }

    tmv::SymBandMatrix<std::complex<T>,U,S> s1(N,noff);
    tmv::SymBandMatrix<std::complex<T>,U,S> s2(N,noff);
    tmv::SymBandMatrix<std::complex<T>,U,S,tmv::FortranStyle> s1f(N,noff);
    tmv::SymBandMatrix<std::complex<T>,U,S,tmv::FortranStyle> s2f(N,noff);

    Assert(int(s1.colsize()) == N && int(s1.rowsize()) == N,"Creating SymMatrix(N)");
    Assert(s1.nlo() == noff && s1.nhi() == noff,"Creating SymMatrix(N)");
    Assert(int(s1f.colsize()) == N && int(s1f.rowsize()) == N,"Creating SymMatrix(N)F");
    Assert(s1f.nlo() == noff && s1f.nhi() == noff,"Creating SymMatrix(N)F");

    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        if (i <= j && j <= i + noff) {
            s1(i,j) = value; 
            s1f(i+1,j+1) = value; 
        }
        if (j <= i && i <= j + noff) {
            s2(i,j) = value; 
            s2f(i+1,j+1) = value; 
        }
    }

    tmv::SymBandMatrixView<std::complex<T> > s1v = s1.view();
    tmv::ConstSymBandMatrixView<std::complex<T> > s1cv = s1.view();
    tmv::SymBandMatrixView<std::complex<T> > s2v = s2.view();
    tmv::ConstSymBandMatrixView<std::complex<T> > s2cv = s2.view();
    tmv::SymBandMatrixView<std::complex<T>,tmv::FortranStyle> s1fv = s1f.view();
    tmv::ConstSymBandMatrixView<std::complex<T>,tmv::FortranStyle> s1fcv = 
        s1f.view();
    tmv::SymBandMatrixView<std::complex<T>,tmv::FortranStyle> s2fv = s2f.view();
    tmv::ConstSymBandMatrixView<std::complex<T>,tmv::FortranStyle> s2fcv = 
        s2f.view();
    const tmv::SymBandMatrix<std::complex<T>,U,S>& s1x = s1;
    const tmv::SymBandMatrix<std::complex<T>,U,S>& s2x = s2;
    const tmv::SymBandMatrix<std::complex<T>,U,S,tmv::FortranStyle>& s1fx = s1f;
    const tmv::SymBandMatrix<std::complex<T>,U,S,tmv::FortranStyle>& s2fx = s2f;

    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if ( j <= i + noff && i <= j + noff) {
            std::complex<T> value(T(k),T(2*k));
            if (i<=j) {
                int j2 = i+noff+1;  if (j2>=N) j2 = N;
                Assert(s1(i,j) == value,"Read/Write SymBandMatrix");
                Assert(s1x(i,j) == value,"Access const SymBandMatrix");
                Assert(s1v(i,j) == value,"Access SymBandMatrix V");
                Assert(s1cv(i,j) == value,"Access SymBandMatrix CV");
                Assert(s1(j,i) == value,"Access SymBandMatrix - opposite tri");
                Assert(s1x(j,i) == value,
                       "Access const SymBandMatrix - opposite tri");
                Assert(s1v(j,i) == value,"Access SymBandMatrix V");
                Assert(s1cv(j,i) == value,"Access SymBandMatrix CV");
                Assert(s1f(i+1,j+1) == value,"Read/Write SymBandMatrixF");
                Assert(s1fx(i+1,j+1) == value,"Access const SymBandMatrixF");
                Assert(s1fv(i+1,j+1) == value,"Access SymBandMatrixF V");
                Assert(s1fcv(i+1,j+1) == value,"Access SymBandMatrixF CV");
                Assert(s1f(j+1,i+1) == value,"Access SymBandMatrixF - opposite tri");
                Assert(s1fx(j+1,i+1) == value,"Access const SymBandMatrixF - opposite tri");
                Assert(s1fv(j+1,i+1) == value,"Access SymBandMatrixF V");
                Assert(s1fcv(j+1,i+1) == value,"Access SymBandMatrixF CV");
                Assert(s1.row(i,i,j2)(j-i) == value,"SymBandMatrix.row");
                Assert(s1x.row(i,i,j2)(j-i) == value,"SymBandMatrix.row");
                Assert(s1cv.row(i,i,j2)(j-i) == value,"SymBandMatrix.row CV");
                Assert(s1v.row(i,i,j2)(j-i) == value,"SymBandMatrix.row V");
                Assert(s1.col(i,i,j2)(j-i) == value,"SymBandMatrix.col");
                Assert(s1x.col(i,i,j2)(j-i) == value,"SymBandMatrix.col");
                Assert(s1cv.col(i,i,j2)(j-i) == value,"SymBandMatrix.col CV");
                Assert(s1v.col(i,i,j2)(j-i) == value,"SymBandMatrix.col V");
                Assert(s1.row(i,j,j2)(0) == value,"SymBandMatrix.row2");
                Assert(s1x.row(i,j,j2)(0) == value,"SymBandMatrix.row2");
                Assert(s1cv.row(i,j,j2)(0) == value,"SymBandMatrix.row2 CV");
                Assert(s1v.row(i,j,j2)(0) == value,"SymBandMatrix.row2 V");
                Assert(s1.col(i,j,j2)(0) == value,"SymBandMatrix.col2");
                Assert(s1x.col(i,j,j2)(0) == value,"SymBandMatrix.col2");
                Assert(s1cv.col(i,j,j2)(0) == value,"SymBandMatrix.col2 CV");
                Assert(s1v.col(i,j,j2)(0) == value,"SymBandMatrix.col2 V");
                Assert(s1f.row(i+1,i+1,j2)(j-i+1) == value,
                       "SymBandMatrixF.row");
                Assert(s1fx.row(i+1,i+1,j2)(j-i+1) == value,
                       "const SymBandMatrixF.row");
                Assert(s1fcv.row(i+1,i+1,j2)(j-i+1) == value,
                       "SymBandMatrixF.row CV");
                Assert(s1fv.row(i+1,i+1,j2)(j-i+1) == value,
                       "SymBandMatrixF.row V");
                Assert(s1f.col(i+1,i+1,j2)(j-i+1) == value,
                       "SymBandMatrixF.col");
                Assert(s1fx.col(i+1,i+1,j2)(j-i+1) == value,
                       "const SymBandMatrixF.col");
                Assert(s1fcv.col(i+1,i+1,j2)(j-i+1) == value,
                       "SymBandMatrixF.col CV");
                Assert(s1fv.col(i+1,i+1,j2)(j-i+1) == value,
                       "SymBandMatrixF.col V");
                Assert(s1f.row(i+1,j+1,j2)(1) == value,"SymBandMatrixF.row2");
                Assert(s1fx.row(i+1,j+1,j2)(1) == value,
                       "const SymBandMatrixF.row2");
                Assert(s1fcv.row(i+1,j+1,j2)(1) == value,
                       "SymBandMatrixF.row2 CV");
                Assert(s1fv.row(i+1,j+1,j2)(1) == value,
                       "SymBandMatrixF.row2 V");
                Assert(s1f.col(i+1,j+1,j2)(1) == value,"SymBandMatrixF.col2");
                Assert(s1fx.col(i+1,j+1,j2)(1) == value,
                       "const SymBandMatrixF.col2");
                Assert(s1fcv.col(i+1,j+1,j2)(1) == value,
                       "SymBandMatrixF.col2 CV");
                Assert(s1fv.col(i+1,j+1,j2)(1) == value,
                       "SymBandMatrixF.col2 V");
                int d = int(j)-int(i);
                if (d==0) {
                    Assert(s1.diag()(i) == value,"SymBandMatrix.diag");
                    Assert(s1x.diag()(i) == value,"const SymBandMatrix.diag");
                    Assert(s1cv.diag()(i) == value,"SymBandMatrix.diag CV");
                    Assert(s1v.diag()(i) == value,"SymBandMatrix.diag V");
                    Assert(s1f.diag()(i+1) == value,"SymBandMatrixF.diag");
                    Assert(s1fx.diag()(i+1) == value,
                           "const SymBandMatrixF.diag");
                    Assert(s1fcv.diag()(i+1) == value,"SymBandMatrixF.diag CV");
                    Assert(s1fv.diag()(i+1) == value,"SymBandMatrixF.diag V");
                }
                Assert(s1.diag(d)(i) == value,"SymBandMatrix.diag1");
                Assert(s1x.diag(d)(i) == value,"const SymBandMatrix.diag1");
                Assert(s1cv.diag(d)(i) == value,"SymBandMatrix.diag1 CV");
                Assert(s1v.diag(d)(i) == value,"SymBandMatrix.diag1 V");
                Assert(s1.diag(d,i,N-d)(0) == value,"SymBandMatrix.diag2");
                Assert(s1x.diag(d,i,N-d)(0) == value,
                       "const SymBandMatrix.diag2");
                Assert(s1cv.diag(d,i,N-d)(0) == value,"SymBandMatrix.diag2 CV");
                Assert(s1v.diag(d,i,N-d)(0) == value,"SymBandMatrix.diag2 V");
                Assert(s1f.diag(d)(i+1) == value,"SymBandMatrixF.diag1");
                Assert(s1fx.diag(d)(i+1) == value,"const SymBandMatrixF.diag1");
                Assert(s1fcv.diag(d)(i+1) == value,"SymBandMatrixF.diag1 CV");
                Assert(s1fv.diag(d)(i+1) == value,"SymBandMatrixF.diag1 V");
                Assert(s1f.diag(d,i+1,N-d)(1) == value,"SymBandMatrixF.diag2");
                Assert(s1fx.diag(d,i+1,N-d)(1) == value,
                       "const SymBandMatrixF.diag2");
                Assert(s1fcv.diag(d,i+1,N-d)(1) == value,
                       "SymBandMatrixF.diag2 CV");
                Assert(s1fv.diag(d,i+1,N-d)(1) == value,
                       "SymBandMatrixF.diag2 V");
                Assert(s1.diag(-d)(i) == value,"SymBandMatrix.diag1");
                Assert(s1x.diag(-d)(i) == value,"const SymBandMatrix.diag1");
                Assert(s1cv.diag(-d)(i) == value,"SymBandMatrix.diag1 CV");
                Assert(s1v.diag(-d)(i) == value,"SymBandMatrix.diag1 V");
                Assert(s1.diag(-d,i,N-d)(0) == value,"SymBandMatrix.diag2");
                Assert(s1x.diag(-d,i,N-d)(0) == value,
                       "const SymBandMatrix.diag2");
                Assert(s1cv.diag(-d,i,N-d)(0) == value,
                       "SymBandMatrix.diag2 CV");
                Assert(s1v.diag(-d,i,N-d)(0) == value,"SymBandMatrix.diag2 V");
                Assert(s1f.diag(-d)(i+1) == value,"SymBandMatrixF.diag1");
                Assert(s1fx.diag(-d)(i+1) == value,
                       "const SymBandMatrixF.diag1");
                Assert(s1fcv.diag(-d)(i+1) == value,"SymBandMatrixF.diag1 CV");
                Assert(s1fv.diag(-d)(i+1) == value,"SymBandMatrixF.diag1 V");
                Assert(s1f.diag(-d,i+1,N-d)(1) == value,"SymBandMatrixF.diag2");
                Assert(s1fx.diag(-d,i+1,N-d)(1) == value,
                       "const SymBandMatrixF.diag2");
                Assert(s1fcv.diag(-d,i+1,N-d)(1) == value,
                       "SymBandMatrixF.diag2 CV");
                Assert(s1fv.diag(-d,i+1,N-d)(1) == value,
                       "SymBandMatrixF.diag2 V");
            }
            if (j<=i) {
                int j1 = 0;  if (int(i)>noff) j1 = i-noff;
                Assert(s2(j,i) == value,"Read/Write SymBandMatrix");
                Assert(s2x(j,i) == value,"Access const SymBandMatrix");
                Assert(s2v(j,i) == value,"Access SymBandMatrix V");
                Assert(s2cv(j,i) == value,"Access SymBandMatrix CV");
                Assert(s2(i,j) == value,"Access SymBandMatrix - opposite tri");
                Assert(s2x(i,j) == value,
                       "Access const SymBandMatrix - opposite tri");
                Assert(s2v(i,j) == value,"Access SymBandMatrix V");
                Assert(s2cv(i,j) == value,"Access SymBandMatrix CV");
                Assert(s2f(j+1,i+1) == value,"Read/Write SymBandMatrixF");
                Assert(s2fx(j+1,i+1) == value,"Access const SymBandMatrixF");
                Assert(s2fv(j+1,i+1) == value,"Access SymBandMatrixF V");
                Assert(s2fcv(j+1,i+1) == value,"Access SymBandMatrixF CV");
                Assert(s2f(i+1,j+1) == value,
                       "Access SymBandMatrixF - opposite tri");
                Assert(s2fx(i+1,j+1) == value,
                       "Access const SymBandMatrixF - opposite tri");
                Assert(s2fv(i+1,j+1) == value,"Access SymBandMatrixF V");
                Assert(s2fcv(i+1,j+1) == value,"Access SymBandMatrixF CV");
                Assert(s2.row(i,j1,i+1)(j-j1) == value,"SymBandMatrix.row");
                Assert(s2x.row(i,j1,i+1)(j-j1) == value,
                       "const SymBandMatrix.row");
                Assert(s2cv.row(i,j1,i+1)(j-j1) == value,
                       "SymBandMatrix.row CV");
                Assert(s2v.row(i,j1,i+1)(j-j1) == value,"SymBandMatrix.row V");
                Assert(s2.col(i,j1,i+1)(j-j1) == value,"SymBandMatrix.col");
                Assert(s2x.col(i,j1,i+1)(j-j1) == value,
                       "const SymBandMatrix.col");
                Assert(s2cv.col(i,j1,i+1)(j-j1) == value,
                       "SymBandMatrix.col CV");
                Assert(s2v.col(i,j1,i+1)(j-j1) == value,"SymBandMatrix.col V");
                Assert(s2.row(i,j,i+1)(0) == value,"SymBandMatrix.row2");
                Assert(s2x.row(i,j,i+1)(0) == value,"const SymBandMatrix.row2");
                Assert(s2cv.row(i,j,i+1)(0) == value,"SymBandMatrix.row2 CV");
                Assert(s2v.row(i,j,i+1)(0) == value,"SymBandMatrix.row2 V");
                Assert(s2.col(i,j,i+1)(0) == value,"SymBandMatrix.col2");
                Assert(s2x.col(i,j,i+1)(0) == value,"const SymBandMatrix.col2");
                Assert(s2cv.col(i,j,i+1)(0) == value,"SymBandMatrix.col2 CV");
                Assert(s2v.col(i,j,i+1)(0) == value,"SymBandMatrix.col2 V");
                Assert(s2f.row(i+1,j1+1,i+1)(j+1-j1) == value,
                       "SymBandMatrixF.row");
                Assert(s2fx.row(i+1,j1+1,i+1)(j+1-j1) == value,
                       "const SymBandMatrixF.row");
                Assert(s2fcv.row(i+1,j1+1,i+1)(j+1-j1) == value,
                       "SymBandMatrixF.row CV");
                Assert(s2fv.row(i+1,j1+1,i+1)(j+1-j1) == value,
                       "SymBandMatrixF.row V");
                Assert(s2f.col(i+1,j1+1,i+1)(j+1-j1) == value,
                       "SymBandMatrixF.col");
                Assert(s2fx.col(i+1,j1+1,i+1)(j+1-j1) == value,
                       "const SymBandMatrixF.col");
                Assert(s2fcv.col(i+1,j1+1,i+1)(j+1-j1) == value,
                       "SymBandMatrixF.col CV");
                Assert(s2fv.col(i+1,j1+1,i+1)(j+1-j1) == value,
                       "SymBandMatrixF.col V");
                Assert(s2f.row(i+1,j+1,i+1)(1) == value,"SymBandMatrixF.row2");
                Assert(s2fx.row(i+1,j+1,i+1)(1) == value,
                       "const SymBandMatrixF.row2");
                Assert(s2fcv.row(i+1,j+1,i+1)(1) == value,
                       "SymBandMatrixF.row2 CV");
                Assert(s2fv.row(i+1,j+1,i+1)(1) == value,
                       "SymBandMatrixF.row2 V");
                Assert(s2f.col(i+1,j+1,i+1)(1) == value,"SymBandMatrixF.col2");
                Assert(s2fx.col(i+1,j+1,i+1)(1) == value,
                       "const SymBandMatrixF.col2");
                Assert(s2fcv.col(i+1,j+1,i+1)(1) == value,
                       "SymBandMatrixF.col2 CV");
                Assert(s2fv.col(i+1,j+1,i+1)(1) == value,
                       "SymBandMatrixF.col2 V");
                int d = int(j)-int(i);
                if (d==0) {
                    Assert(s2.diag()(i) == value,"SymBandMatrix.diag");
                    Assert(s2x.diag()(i) == value,"const SymBandMatrix.diag");
                    Assert(s2cv.diag()(i) == value,"SymBandMatrix.diag CV");
                    Assert(s2v.diag()(i) == value,"SymBandMatrix.diag V");
                    Assert(s2f.diag()(i+1) == value,"SymBandMatrixF.diag");
                    Assert(s2fx.diag()(i+1) == value,
                           "const SymBandMatrixF.diag");
                    Assert(s2fcv.diag()(i+1) == value,"SymBandMatrixF.diag CV");
                    Assert(s2fv.diag()(i+1) == value,"SymBandMatrixF.diag V");
                }
                Assert(s2.diag(d)(j) == value,"SymBandMatrix.diag1");
                Assert(s2x.diag(d)(j) == value,"const SymBandMatrix.diag1");
                Assert(s2cv.diag(d)(j) == value,"SymBandMatrix.diag1 CV");
                Assert(s2v.diag(d)(j) == value,"SymBandMatrix.diag1 V");
                Assert(s2.diag(d,j,N+d)(0) == value,"SymBandMatrix.diag2");
                Assert(s2x.diag(d,j,N+d)(0) == value,
                       "const SymBandMatrix.diag2");
                Assert(s2cv.diag(d,j,N+d)(0) == value,"SymBandMatrix.diag2 CV");
                Assert(s2v.diag(d,j,N+d)(0) == value,"SymBandMatrix.diag2 V");
                Assert(s2f.diag(d)(j+1) == value,"SymBandMatrixF.diag1");
                Assert(s2fx.diag(d)(j+1) == value,"const SymBandMatrixF.diag1");
                Assert(s2fcv.diag(d)(j+1) == value,"SymBandMatrixF.diag1 CV");
                Assert(s2fv.diag(d)(j+1) == value,"SymBandMatrixF.diag1 V");
                Assert(s2f.diag(d,j+1,N+d)(1) == value,"SymBandMatrixF.diag2");
                Assert(s2fx.diag(d,j+1,N+d)(1) == value,
                       "const SymBandMatrixF.diag2");
                Assert(s2fcv.diag(d,j+1,N+d)(1) == value,
                       "SymBandMatrixF.diag2 CV");
                Assert(s2fv.diag(d,j+1,N+d)(1) == value,
                       "SymBandMatrixF.diag2 V");
                Assert(s2.diag(-d)(j) == value,"SymBandMatrix.diag1");
                Assert(s2x.diag(-d)(j) == value,"const SymBandMatrix.diag1");
                Assert(s2cv.diag(-d)(j) == value,"SymBandMatrix.diag1 CV");
                Assert(s2v.diag(-d)(j) == value,"SymBandMatrix.diag1 V");
                Assert(s2.diag(-d,j,N+d)(0) == value,"SymBandMatrix.diag2");
                Assert(s2x.diag(-d,j,N+d)(0) == value,
                       "const SymBandMatrix.diag2");
                Assert(s2cv.diag(-d,j,N+d)(0) == value,
                       "SymBandMatrix.diag2 CV");
                Assert(s2v.diag(-d,j,N+d)(0) == value,"SymBandMatrix.diag2 V");
                Assert(s2f.diag(-d)(j+1) == value,"SymBandMatrixF.diag1");
                Assert(s2fx.diag(-d)(j+1) == value,
                       "const SymBandMatrixF.diag1");
                Assert(s2fcv.diag(-d)(j+1) == value,"SymBandMatrixF.diag1 CV");
                Assert(s2fv.diag(-d)(j+1) == value,"SymBandMatrixF.diag1 V");
                Assert(s2f.diag(-d,j+1,N+d)(1) == value,"SymBandMatrixF.diag2");
                Assert(s2fx.diag(-d,j+1,N+d)(1) == value,
                       "const SymBandMatrixF.diag2");
                Assert(s2fcv.diag(-d,j+1,N+d)(1) == value,
                       "SymBandMatrixF.diag2 CV");
                Assert(s2fv.diag(-d,j+1,N+d)(1) == value,
                       "SymBandMatrixF.diag2 V");
            }
        }
    }

    Assert(s1 == s1f,"CStyle SymBandMatrix == FortranStyle SymBandMatrix");
    Assert(s1 == s1cv,"SymBandMatrix == ConstSymBandMatrixView");
    Assert(s1 == s1v,"SymBandMatrix == SymBandMatrixView");
    Assert(s1 == s1fcv,"SymBandMatrix == FortranStyle ConstSymBandMatrixView");
    Assert(s1 == s1fv,"SymBandMatrix == FortranStyle SymBandMatrixView");
    Assert(s2 == s2f,"CStyle SymBandMatrix == FortranStyle SymBandMatrix");
    Assert(s2 == s2cv,"SymBandMatrix == ConstSymBandMatrixView");
    Assert(s2 == s2v,"SymBandMatrix == SymBandMatrixView");
    Assert(s2 == s2fcv,"SymBandMatrix == FortranStyle ConstSymBandMatrixView");
    Assert(s2 == s2fv,"SymBandMatrix == FortranStyle SymBandMatrixView");

    s1.resize(3,1);
    Assert(s1.colsize() == 3 && s1.rowsize() == 3,
           "SymBandMatrix s1.resize(3,1) sizes");
    Assert(s1.colsize() == 3 && s1.rowsize() == 3,
           "SymBandMatrix s1.resize(3,1) sizes");
    s2.resize(3,1);
    Assert(s2.nlo() == 1 && s2.nhi() == 1,
           "SymBandMatrix s2.resize(3,1) nlo,nhi");
    Assert(s2.nlo() == 1 && s2.nhi() == 1,
           "SymBandMatrix s2.resize(3,1) nlo,nhi");
    s2.resize(3,1);
    Assert(s2.colsize() == 3 && s2.rowsize() == 3,
           "SymBandMatrix s2.resize(3,1) sizes");
    Assert(s2.nlo() == 1 && s2.nhi() == 1,
           "SymBandMatrix s2.resize(3,1) nlo,nhi");
    for (int i=0, k=0; i<3; ++i) for (int j=0; j<3; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        if (i <= j && j <= i + 1) s1(i,j) = value; 
        if (j <= i && i <= j + 1) s2(i,j) = value; 
    }
    for (int i=0, k=0; i<3; ++i) for (int j=0; j<3; ++j, ++k) {
        if ( j <= i + 1 && i <= j + 1) {
            std::complex<T> value(T(k),T(2*k));
            if (i<=j) {
                Assert(s1(i,j) == value,"Read/Write resized SymBandMatrix");
                Assert(s1(j,i) == value,
                       "Read/Write resized SymBandMatrix opp tri");
            }
            if (j<=i) {
                Assert(s2(i,j) == value,"Read/Write resized SymBandMatrix");
                Assert(s2(j,i) == value,
                       "Read/Write resized SymBandMatrix opp tri");
            }
        }
    }

    s1.resize(2*N,5);
    Assert(int(s1.colsize()) == 2*N && int(s1.rowsize()) == 2*N,
           "SymBandMatrix s1.resize(2*N,5) sizes");
    Assert(s1.nlo() == 5 && s1.nhi() == 5,
           "SymBandMatrix s1.resize(2*N,5) nlo,nhi");
    s2.resize(2*N,5);
    Assert(int(s2.colsize()) == 2*N && int(s2.rowsize()) == 2*N,
           "SymBandMatrix s2.resize(2*N,5) sizes");
    Assert(s2.nlo() == 5 && s2.nhi() == 5,
           "SymBandMatrix s2.resize(2*N,5) nlo,nhi");
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<2*N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        if (i <= j && j <= i + 5) s1(i,j) = value; 
        if (j <= i && i <= j + 5) s2(i,j) = value; 
    }
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<2*N; ++j, ++k) {
        if ( j <= i + 5 && i <= j + 5) {
            std::complex<T> value(T(k),T(2*k));
            if (i<=j) {
                Assert(s1(i,j) == value,"Read/Write resized SymBandMatrix");
                Assert(s1(j,i) == value,
                       "Read/Write resized SymBandMatrix opp tri");
            }
            if (j<=i) {
                Assert(s2(i,j) == value,"Read/Write resized SymBandMatrix");
                Assert(s2(j,i) == value,
                       "Read/Write resized SymBandMatrix opp tri");
            }
        }
    }
}

template <class T, tmv::UpLoType U, tmv::StorageType S>
static void TestBasicHermBandMatrix_1()
{
    const int N = 10;
    const int noff = 3;

    if (showstartdone) {
        std::cout<<"Start TestBasicHermBandMatrix_1\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"U = "<<tmv::TMV_Text(U)<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
        std::cout<<"noff = "<<noff<<std::endl;
    }

    tmv::HermBandMatrix<std::complex<T>,U,S> h1(N,noff);
    tmv::HermBandMatrix<std::complex<T>,U,S> h2(N,noff);
    tmv::HermBandMatrix<std::complex<T>,U,S,tmv::FortranStyle> h1f(N,noff);
    tmv::HermBandMatrix<std::complex<T>,U,S,tmv::FortranStyle> h2f(N,noff);

    Assert(int(h1.colsize()) == N && int(h1.rowsize()) == N,
           "Creating HermMatrix(N)");
    Assert(h1.nlo() == noff && h1.nhi() == noff,"Creating HermMatrix(N)");
    Assert(int(h1f.colsize()) == N && int(h1f.rowsize()) == N,
           "Creating HermMatrix(N)F");
    Assert(h1f.nlo() == noff && h1f.nhi() == noff,"Creating HermMatrix(N)F");

    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (i <= j && j <= i + noff) {
            h1(i,j) = hvalue; 
            h1f(i+1,j+1) = hvalue;
        }
        if (j <= i && i <= j + noff) {
            h2(i,j) = hvalue; 
            h2f(i+1,j+1) = hvalue; 
        }
    }

    tmv::SymBandMatrixView<std::complex<T> > h1v = h1.view();
    tmv::ConstSymBandMatrixView<std::complex<T> > h1cv = h1.view();
    tmv::SymBandMatrixView<std::complex<T> > h2v = h2.view();
    tmv::ConstSymBandMatrixView<std::complex<T> > h2cv = h2.view();
    tmv::SymBandMatrixView<std::complex<T>,tmv::FortranStyle> h1fv = h1f.view();
    tmv::ConstSymBandMatrixView<std::complex<T>,tmv::FortranStyle> h1fcv = 
        h1f.view();
    tmv::SymBandMatrixView<std::complex<T>,tmv::FortranStyle> h2fv = h2f.view();
    tmv::ConstSymBandMatrixView<std::complex<T>,tmv::FortranStyle> h2fcv = 
        h2f.view();
    const tmv::HermBandMatrix<std::complex<T>,U,S>& h1x = h1;
    const tmv::HermBandMatrix<std::complex<T>,U,S>& h2x = h2;
    const tmv::HermBandMatrix<std::complex<T>,U,S,tmv::FortranStyle>& h1fx = h1f;
    const tmv::HermBandMatrix<std::complex<T>,U,S,tmv::FortranStyle>& h2fx = h2f;

    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if ( j <= i + noff && i <= j + noff) {
            std::complex<T> value(T(k),T(2*k));
            std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
            if (i<=j) {
                int j2 = i+noff+1;  if (j2>=N) j2 = N;
                Assert(h1(i,j) == hvalue,"Read/Write HermBandMatrix");
                Assert(h1x(i,j) == hvalue,"Access const HermBandMatrix");
                Assert(h1v(i,j) == hvalue,"Access HermBandMatrix V");
                Assert(h1cv(i,j) == hvalue,"Access HermBandMatrix CV");
                Assert(h1(j,i) == conj(hvalue),
                       "Access HermBandMatrix - opposite tri");
                Assert(h1x(j,i) == conj(hvalue),
                       "Access const HermBandMatrix - opposite tri");
                Assert(h1v(j,i) == conj(hvalue),"Access HermBandMatrix V");
                Assert(h1cv(j,i) == conj(hvalue),"Access HermBandMatrix CV");
                Assert(h1f(i+1,j+1) == hvalue,"Read/Write HermBandMatrixF");
                Assert(h1fx(i+1,j+1) == hvalue,"Access const HermBandMatrixF");
                Assert(h1fv(i+1,j+1) == hvalue,"Access HermBandMatrixF V");
                Assert(h1fcv(i+1,j+1) == hvalue,"Access HermBandMatrixF CV");
                Assert(h1f(j+1,i+1) == conj(hvalue),
                       "Access HermBandMatrixF - opposite tri");
                Assert(h1fx(j+1,i+1) == conj(hvalue),
                       "Access const HermBandMatrixF - opposite tri");
                Assert(h1fv(j+1,i+1) == conj(hvalue),
                       "Access HermBandMatrixF V");
                Assert(h1fcv(j+1,i+1) == conj(hvalue),
                       "Access HermBandMatrixF CV");
                Assert(h1.row(i,i,j2)(j-i) == hvalue,"HermBandMatrix.row");
                Assert(h1x.row(i,i,j2)(j-i) == hvalue,"HermBandMatrix.row");
                Assert(h1cv.row(i,i,j2)(j-i) == hvalue,"HermBandMatrix.row CV");
                Assert(h1v.row(i,i,j2)(j-i) == hvalue,"HermBandMatrix.row V");
                Assert(h1.col(i,i,j2)(j-i) == conj(hvalue),
                       "HermBandMatrix.col");
                Assert(h1x.col(i,i,j2)(j-i) == conj(hvalue),
                       "HermBandMatrix.col");
                Assert(h1cv.col(i,i,j2)(j-i) == conj(hvalue),
                       "HermBandMatrix.col CV");
                Assert(h1v.col(i,i,j2)(j-i) == conj(hvalue),
                       "HermBandMatrix.col V");
                Assert(h1.row(i,j,j2)(0) == hvalue,"HermBandMatrix.row2");
                Assert(h1x.row(i,j,j2)(0) == hvalue,"HermBandMatrix.row2");
                Assert(h1cv.row(i,j,j2)(0) == hvalue,"HermBandMatrix.row2 CV");
                Assert(h1v.row(i,j,j2)(0) == hvalue,"HermBandMatrix.row2 V");
                Assert(h1.col(i,j,j2)(0) == conj(hvalue),"HermBandMatrix.col2");
                Assert(h1x.col(i,j,j2)(0) == conj(hvalue),
                       "HermBandMatrix.col2");
                Assert(h1cv.col(i,j,j2)(0) == conj(hvalue),
                       "HermBandMatrix.col2 CV");
                Assert(h1v.col(i,j,j2)(0) == conj(hvalue),
                       "HermBandMatrix.col2 V");
                Assert(h1f.row(i+1,i+1,j2)(j-i+1) == hvalue,
                       "HermBandMatrixF.row");
                Assert(h1fx.row(i+1,i+1,j2)(j-i+1) == hvalue,
                       "const HermBandMatrixF.row");
                Assert(h1fcv.row(i+1,i+1,j2)(j-i+1) == hvalue,
                       "HermBandMatrixF.row CV");
                Assert(h1fv.row(i+1,i+1,j2)(j-i+1) == hvalue,
                       "HermBandMatrixF.row V");
                Assert(h1f.col(i+1,i+1,j2)(j-i+1) == conj(hvalue),
                       "HermBandMatrixF.col");
                Assert(h1fx.col(i+1,i+1,j2)(j-i+1) == conj(hvalue),
                       "const HermBandMatrixF.col");
                Assert(h1fcv.col(i+1,i+1,j2)(j-i+1) == conj(hvalue),
                       "HermBandMatrixF.col CV");
                Assert(h1fv.col(i+1,i+1,j2)(j-i+1) == conj(hvalue),
                       "HermBandMatrixF.col V");
                Assert(h1f.row(i+1,j+1,j2)(1) == hvalue,"HermBandMatrixF.row2");
                Assert(h1fx.row(i+1,j+1,j2)(1) == hvalue,
                       "const HermBandMatrixF.row2");
                Assert(h1fcv.row(i+1,j+1,j2)(1) == hvalue,
                       "HermBandMatrixF.row2 CV");
                Assert(h1fv.row(i+1,j+1,j2)(1) == hvalue,
                       "HermBandMatrixF.row2 V");
                Assert(h1f.col(i+1,j+1,j2)(1) == conj(hvalue),
                       "HermBandMatrixF.col2");
                Assert(h1fx.col(i+1,j+1,j2)(1) == conj(hvalue),
                       "const HermBandMatrixF.col2");
                Assert(h1fcv.col(i+1,j+1,j2)(1) == conj(hvalue),
                       "HermBandMatrixF.col2 CV");
                Assert(h1fv.col(i+1,j+1,j2)(1) == conj(hvalue),
                       "HermBandMatrixF.col2 V");
                int d = int(j)-int(i);
                if (d==0) {
                    Assert(h1.diag()(i) == hvalue,"HermBandMatrix.diag");
                    Assert(h1x.diag()(i) == hvalue,"const HermBandMatrix.diag");
                    Assert(h1cv.diag()(i) == hvalue,"HermBandMatrix.diag CV");
                    Assert(h1v.diag()(i) == hvalue,"HermBandMatrix.diag V");
                    Assert(h1f.diag()(i+1) == hvalue,"HermBandMatrixF.diag");
                    Assert(h1fx.diag()(i+1) == hvalue,
                           "const HermBandMatrixF.diag");
                    Assert(h1fcv.diag()(i+1) == hvalue,
                           "HermBandMatrixF.diag CV");
                    Assert(h1fv.diag()(i+1) == hvalue,"HermBandMatrixF.diag V");
                }
                Assert(h1.diag(d)(i) == hvalue,"HermBandMatrix.diag1");
                Assert(h1x.diag(d)(i) == hvalue,"const HermBandMatrix.diag1");
                Assert(h1cv.diag(d)(i) == hvalue,"HermBandMatrix.diag1 CV");
                Assert(h1v.diag(d)(i) == hvalue,"HermBandMatrix.diag1 V");
                Assert(h1.diag(d,i,N-d)(0) == hvalue,"HermBandMatrix.diag2");
                Assert(h1x.diag(d,i,N-d)(0) == hvalue,
                       "const HermBandMatrix.diag2");
                Assert(h1cv.diag(d,i,N-d)(0) == hvalue,
                       "HermBandMatrix.diag2 CV");
                Assert(h1v.diag(d,i,N-d)(0) == hvalue,"HermBandMatrix.diag2 V");
                Assert(h1f.diag(d)(i+1) == hvalue,"HermBandMatrixF.diag1");
                Assert(h1fx.diag(d)(i+1) == hvalue,
                       "const HermBandMatrixF.diag1");
                Assert(h1fcv.diag(d)(i+1) == hvalue,"HermBandMatrixF.diag1 CV");
                Assert(h1fv.diag(d)(i+1) == hvalue,"HermBandMatrixF.diag1 V");
                Assert(h1f.diag(d,i+1,N-d)(1) == hvalue,
                       "HermBandMatrixF.diag2");
                Assert(h1fx.diag(d,i+1,N-d)(1) == hvalue,
                       "const HermBandMatrixF.diag2");
                Assert(h1fcv.diag(d,i+1,N-d)(1) == hvalue,
                       "HermBandMatrixF.diag2 CV");
                Assert(h1fv.diag(d,i+1,N-d)(1) == hvalue,
                       "HermBandMatrixF.diag2 V");
                Assert(h1.diag(-d)(i) == conj(hvalue),"HermBandMatrix.diag1");
                Assert(h1x.diag(-d)(i) == conj(hvalue),
                       "const HermBandMatrix.diag1");
                Assert(h1cv.diag(-d)(i) == conj(hvalue),
                       "HermBandMatrix.diag1 CV");
                Assert(h1v.diag(-d)(i) == conj(hvalue),
                       "HermBandMatrix.diag1 V");
                Assert(h1.diag(-d,i,N-d)(0) == conj(hvalue),
                       "HermBandMatrix.diag2");
                Assert(h1x.diag(-d,i,N-d)(0) == conj(hvalue),
                       "const HermBandMatrix.diag2");
                Assert(h1cv.diag(-d,i,N-d)(0) == conj(hvalue),
                       "HermBandMatrix.diag2 CV");
                Assert(h1v.diag(-d,i,N-d)(0) == conj(hvalue),
                       "HermBandMatrix.diag2 V");
                Assert(h1f.diag(-d)(i+1) == conj(hvalue),
                       "HermBandMatrixF.diag1");
                Assert(h1fx.diag(-d)(i+1) == conj(hvalue),
                       "const HermBandMatrixF.diag1");
                Assert(h1fcv.diag(-d)(i+1) == conj(hvalue),
                       "HermBandMatrixF.diag1 CV");
                Assert(h1fv.diag(-d)(i+1) == conj(hvalue),
                       "HermBandMatrixF.diag1 V");
                Assert(h1f.diag(-d,i+1,N-d)(1) == conj(hvalue),
                       "HermBandMatrixF.diag2");
                Assert(h1fx.diag(-d,i+1,N-d)(1) == conj(hvalue),
                       "const HermBandMatrixF.diag2");
                Assert(h1fcv.diag(-d,i+1,N-d)(1) == conj(hvalue),
                       "HermBandMatrixF.diag2 CV");
                Assert(h1fv.diag(-d,i+1,N-d)(1) == conj(hvalue),
                       "HermBandMatrixF.diag2 V");
            }
            if (j<=i) {
                int j1 = 0;  if (int(i)>noff) j1 = i-noff;
                Assert(h2(j,i) == conj(hvalue),"Read/Write HermBandMatrix");
                Assert(h2x(j,i) == conj(hvalue),"Access const HermBandMatrix");
                Assert(h2v(j,i) == conj(hvalue),"Access HermBandMatrix V");
                Assert(h2cv(j,i) == conj(hvalue),"Access HermBandMatrix CV");
                Assert(h2(i,j) == hvalue,
                       "Access HermBandMatrix - opposite tri");
                Assert(h2x(i,j) == hvalue,
                       "Access const HermBandMatrix - opposite tri");
                Assert(h2v(i,j) == hvalue,"Access HermBandMatrix V");
                Assert(h2cv(i,j) == hvalue,"Access HermBandMatrix CV");
                Assert(h2f(j+1,i+1) == conj(hvalue),
                       "Read/Write HermBandMatrixF");
                Assert(h2fx(j+1,i+1) == conj(hvalue),
                       "Access const HermBandMatrixF");
                Assert(h2fv(j+1,i+1) == conj(hvalue),
                       "Access HermBandMatrixF V");
                Assert(h2fcv(j+1,i+1) == conj(hvalue),
                       "Access HermBandMatrixF CV");
                Assert(h2f(i+1,j+1) == hvalue,
                       "Access HermBandMatrixF - opposite tri");
                Assert(h2fx(i+1,j+1) == hvalue,
                       "Access const HermBandMatrixF - opposite tri");
                Assert(h2fv(i+1,j+1) == hvalue,"Access HermBandMatrixF V");
                Assert(h2fcv(i+1,j+1) == hvalue,"Access HermBandMatrixF CV");
                Assert(h2.row(i,j1,i+1)(j-j1) == hvalue,"HermBandMatrix.row");
                Assert(h2x.row(i,j1,i+1)(j-j1) == hvalue,
                       "const HermBandMatrix.row");
                Assert(h2cv.row(i,j1,i+1)(j-j1) == hvalue,
                       "HermBandMatrix.row CV");
                Assert(h2v.row(i,j1,i+1)(j-j1) == hvalue,
                       "HermBandMatrix.row V");
                Assert(h2.col(i,j1,i+1)(j-j1) == conj(hvalue),
                       "HermBandMatrix.col");
                Assert(h2x.col(i,j1,i+1)(j-j1) == conj(hvalue),
                       "const HermBandMatrix.col");
                Assert(h2cv.col(i,j1,i+1)(j-j1) == conj(hvalue),
                       "HermBandMatrix.col CV");
                Assert(h2v.col(i,j1,i+1)(j-j1) == conj(hvalue),
                       "HermBandMatrix.col V");
                Assert(h2.row(i,j,i+1)(0) == hvalue,"HermBandMatrix.row2");
                Assert(h2x.row(i,j,i+1)(0) == hvalue,
                       "const HermBandMatrix.row2");
                Assert(h2cv.row(i,j,i+1)(0) == hvalue,"HermBandMatrix.row2 CV");
                Assert(h2v.row(i,j,i+1)(0) == hvalue,"HermBandMatrix.row2 V");
                Assert(h2.col(i,j,i+1)(0) == conj(hvalue),
                       "HermBandMatrix.col2");
                Assert(h2x.col(i,j,i+1)(0) == conj(hvalue),
                       "const HermBandMatrix.col2");
                Assert(h2cv.col(i,j,i+1)(0) == conj(hvalue),
                       "HermBandMatrix.col2 CV");
                Assert(h2v.col(i,j,i+1)(0) == conj(hvalue),
                       "HermBandMatrix.col2 V");
                Assert(h2f.row(i+1,j1+1,i+1)(j+1-j1) == hvalue,
                       "HermBandMatrixF.row");
                Assert(h2fx.row(i+1,j1+1,i+1)(j+1-j1) == hvalue,
                       "const HermBandMatrixF.row");
                Assert(h2fcv.row(i+1,j1+1,i+1)(j+1-j1) == hvalue,
                       "HermBandMatrixF.row CV");
                Assert(h2fv.row(i+1,j1+1,i+1)(j+1-j1) == hvalue,
                       "HermBandMatrixF.row V");
                Assert(h2f.col(i+1,j1+1,i+1)(j+1-j1) == conj(hvalue),
                       "HermBandMatrixF.col");
                Assert(h2fx.col(i+1,j1+1,i+1)(j+1-j1) == conj(hvalue),
                       "const HermBandMatrixF.col");
                Assert(h2fcv.col(i+1,j1+1,i+1)(j+1-j1) == conj(hvalue),
                       "HermBandMatrixF.col CV");
                Assert(h2fv.col(i+1,j1+1,i+1)(j+1-j1) == conj(hvalue),
                       "HermBandMatrixF.col V");
                Assert(h2f.row(i+1,j+1,i+1)(1) == hvalue,
                       "HermBandMatrixF.row2");
                Assert(h2fx.row(i+1,j+1,i+1)(1) == hvalue,
                       "const HermBandMatrixF.row2");
                Assert(h2fcv.row(i+1,j+1,i+1)(1) == hvalue,
                       "HermBandMatrixF.row2 CV");
                Assert(h2fv.row(i+1,j+1,i+1)(1) == hvalue,
                       "HermBandMatrixF.row2 V");
                Assert(h2f.col(i+1,j+1,i+1)(1) == conj(hvalue),
                       "HermBandMatrixF.col2");
                Assert(h2fx.col(i+1,j+1,i+1)(1) == conj(hvalue),
                       "const HermBandMatrixF.col2");
                Assert(h2fcv.col(i+1,j+1,i+1)(1) == conj(hvalue),
                       "HermBandMatrixF.col2 CV");
                Assert(h2fv.col(i+1,j+1,i+1)(1) == conj(hvalue),
                       "HermBandMatrixF.col2 V");
                int d = int(j)-int(i);
                if (d==0) {
                    Assert(h2.diag()(i) == hvalue,"HermBandMatrix.diag");
                    Assert(h2x.diag()(i) == hvalue,"const HermBandMatrix.diag");
                    Assert(h2cv.diag()(i) == hvalue,"HermBandMatrix.diag CV");
                    Assert(h2v.diag()(i) == hvalue,"HermBandMatrix.diag V");
                    Assert(h2f.diag()(i+1) == hvalue,"HermBandMatrixF.diag");
                    Assert(h2fx.diag()(i+1) == hvalue,
                           "const HermBandMatrixF.diag");
                    Assert(h2fcv.diag()(i+1) == hvalue,
                           "HermBandMatrixF.diag CV");
                    Assert(h2fv.diag()(i+1) == hvalue,"HermBandMatrixF.diag V");
                }
                Assert(h2.diag(d)(j) == hvalue,"HermBandMatrix.diag1");
                Assert(h2x.diag(d)(j) == hvalue,"const HermBandMatrix.diag1");
                Assert(h2cv.diag(d)(j) == hvalue,"HermBandMatrix.diag1 CV");
                Assert(h2v.diag(d)(j) == hvalue,"HermBandMatrix.diag1 V");
                Assert(h2.diag(d,j,N+d)(0) == hvalue,"HermBandMatrix.diag2");
                Assert(h2x.diag(d,j,N+d)(0) == hvalue,
                       "const HermBandMatrix.diag2");
                Assert(h2cv.diag(d,j,N+d)(0) == hvalue,
                       "HermBandMatrix.diag2 CV");
                Assert(h2v.diag(d,j,N+d)(0) == hvalue,"HermBandMatrix.diag2 V");
                Assert(h2f.diag(d)(j+1) == hvalue,"HermBandMatrixF.diag1");
                Assert(h2fx.diag(d)(j+1) == hvalue,
                       "const HermBandMatrixF.diag1");
                Assert(h2fcv.diag(d)(j+1) == hvalue,"HermBandMatrixF.diag1 CV");
                Assert(h2fv.diag(d)(j+1) == hvalue,"HermBandMatrixF.diag1 V");
                Assert(h2f.diag(d,j+1,N+d)(1) == hvalue,
                       "HermBandMatrixF.diag2");
                Assert(h2fx.diag(d,j+1,N+d)(1) == hvalue,
                       "const HermBandMatrixF.diag2");
                Assert(h2fcv.diag(d,j+1,N+d)(1) == hvalue,
                       "HermBandMatrixF.diag2 CV");
                Assert(h2fv.diag(d,j+1,N+d)(1) == hvalue,
                       "HermBandMatrixF.diag2 V");
                Assert(h2.diag(-d)(j) == conj(hvalue),"HermBandMatrix.diag1");
                Assert(h2x.diag(-d)(j) == conj(hvalue),
                       "const HermBandMatrix.diag1");
                Assert(h2cv.diag(-d)(j) == conj(hvalue),
                       "HermBandMatrix.diag1 CV");
                Assert(h2v.diag(-d)(j) == conj(hvalue),
                       "HermBandMatrix.diag1 V");
                Assert(h2.diag(-d,j,N+d)(0) == conj(hvalue),
                       "HermBandMatrix.diag2");
                Assert(h2x.diag(-d,j,N+d)(0) == conj(hvalue),
                       "const HermBandMatrix.diag2");
                Assert(h2cv.diag(-d,j,N+d)(0) == conj(hvalue),
                       "HermBandMatrix.diag2 CV");
                Assert(h2v.diag(-d,j,N+d)(0) == conj(hvalue),
                       "HermBandMatrix.diag2 V");
                Assert(h2f.diag(-d)(j+1) == conj(hvalue),
                       "HermBandMatrixF.diag1");
                Assert(h2fx.diag(-d)(j+1) == conj(hvalue),
                       "const HermBandMatrixF.diag1");
                Assert(h2fcv.diag(-d)(j+1) == conj(hvalue),
                       "HermBandMatrixF.diag1 CV");
                Assert(h2fv.diag(-d)(j+1) == conj(hvalue),
                       "HermBandMatrixF.diag1 V");
                Assert(h2f.diag(-d,j+1,N+d)(1) == conj(hvalue),
                       "HermBandMatrixF.diag2");
                Assert(h2fx.diag(-d,j+1,N+d)(1) == conj(hvalue),
                       "const HermBandMatrixF.diag2");
                Assert(h2fcv.diag(-d,j+1,N+d)(1) == conj(hvalue),
                       "HermBandMatrixF.diag2 CV");
                Assert(h2fv.diag(-d,j+1,N+d)(1) == conj(hvalue),
                       "HermBandMatrixF.diag2 V");
            }
        }
    }

    Assert(h1 == h1f,"CStyle HermBandMatrix == FortranStyle HermBandMatrix");
    Assert(h1 == h1cv,"HermBandMatrix == ConstSymBandMatrixView");
    Assert(h1 == h1v,"HermBandMatrix == SymBandMatrixView");
    Assert(h1 == h1fcv,"HermBandMatrix == FortranStyle ConstSymBandMatrixView");
    Assert(h1 == h1fv,"HermBandMatrix == FortranStyle SymBandMatrixView");
    Assert(h2 == h2f,"CStyle HermBandMatrix == FortranStyle HermBandMatrix");
    Assert(h2 == h2cv,"HermBandMatrix == ConstSymBandMatrixView");
    Assert(h2 == h2v,"HermBandMatrix == SymBandMatrixView");
    Assert(h2 == h2fcv,"HermBandMatrix == FortranStyle ConstSymBandMatrixView");
    Assert(h2 == h2fv,"HermBandMatrix == FortranStyle SymBandMatrixView");

    h1.resize(3,1);
    Assert(h1.colsize() == 3 && h1.rowsize() == 3,
           "HermBandMatrix h1.resize(3,1) sizes");
    Assert(h1.nlo() == 1 && h1.nhi() == 1,
           "HermBandMatrix h1.resize(3,1) nlo,nhi");
    h2.resize(3,1);
    Assert(h2.colsize() == 3 && h2.rowsize() == 3,
           "HermBandMatrix h2.resize(3,1) sizes");
    Assert(h2.nlo() == 1 && h2.nhi() == 1,
           "HermBandMatrix h2.resize(3,1) nlo,nhi");
    for (int i=0, k=0; i<3; ++i) for (int j=0; j<3; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (i <= j && j <= i + 1) h1(i,j) = hvalue; 
        if (j <= i && i <= j + 1) h2(i,j) = hvalue; 
    }
    for (int i=0, k=0; i<3; ++i) for (int j=0; j<3; ++j, ++k) {
        if ( j <= i + 1 && i <= j + 1) {
            std::complex<T> value(T(k),T(2*k));
            std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
            if (i<=j) {
                Assert(h1(i,j) == hvalue,"Read/Write resized HermBandMatrix");
                Assert(h1(j,i) == std::conj(hvalue),
                       "Read/Write resized HermBandMatrix opp tri");
            }
            if (j<=i) {
                Assert(h2(i,j) == hvalue,"Read/Write resized HermBandMatrix");
                Assert(h2(j,i) == std::conj(hvalue),
                       "Read/Write resized HermBandMatrix opp tri");
            }
        }
    }

    h1.resize(2*N,5);
    Assert(int(h1.colsize()) == 2*N && int(h1.rowsize()) == 2*N,
           "HermBandMatrix h1.resize(2*N,5) sizes");
    Assert(h1.nlo() == 5 && h1.nhi() == 5,
           "HermBandMatrix h1.resize(2*N,5) nlo,nhi");
    h2.resize(2*N,5);
    Assert(int(h2.colsize()) == 2*N && int(h2.rowsize()) == 2*N,
           "HermBandMatrix h2.resize(2*N,5) sizes");
    Assert(h2.nlo() == 5 && h2.nhi() == 5,
           "HermBandMatrix h2.resize(2*N,5) nlo,nhi");
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<2*N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (i <= j && j <= i + 5) h1(i,j) = hvalue; 
        if (j <= i && i <= j + 5) h2(i,j) = hvalue; 
    }
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<2*N; ++j, ++k) {
        if ( j <= i + 5 && i <= j + 5) {
            std::complex<T> value(T(k),T(2*k));
            std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
            if (i<=j) {
                Assert(h1(i,j) == hvalue,"Read/Write resized HermBandMatrix");
                Assert(h1(j,i) == std::conj(hvalue),
                       "Read/Write resized HermBandMatrix opp tri");
            }
            if (j<=i) {
                Assert(h2(i,j) == hvalue,"Read/Write resized HermBandMatrix");
                Assert(h2(j,i) == std::conj(hvalue),
                       "Read/Write resized HermBandMatrix opp tri");
            }
        }
    }
}

template <class T, tmv::UpLoType U, tmv::StorageType S>
static void TestBasicSymBandMatrix_2()
{
    const int N = 10;
    const int noff = 3;

    if (showstartdone) {
        std::cout<<"Start TestBasicSymBandMatrix_2\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"U = "<<tmv::TMV_Text(U)<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
        std::cout<<"noff = "<<noff<<std::endl;
    }

    // Test assignments and constructors from arrays
    const T quarrm[] = {
        T(0), T(2), T(4),
              T(1), T(3), T(5),
                    T(2), T(4), T(6),
                          T(3), T(5),
                                T(4)
    };
    const T qlarrm[] = {
        T(0),
        T(2), T(1),
        T(4), T(3), T(2),
              T(5), T(4), T(3),
                    T(6), T(5), T(4)
    };
    const T quarcm[] = {
        T(0),
        T(2), T(1),
        T(4), T(3), T(2),
              T(5), T(4), T(3),
                    T(6), T(5), T(4)
    };
    const T qlarcm[] = {
        T(0), T(2), T(4),
              T(1), T(3), T(5),
                    T(2), T(4), T(6),
                          T(3), T(5),
                                T(4)
    };
    const T quardm[] = {
        T(0), T(1), T(2), T(3), T(4),
           T(2), T(3), T(4), T(5),
              T(4), T(5), T(6)
    };
    const T qlardm[] = {
              T(4), T(5), T(6),
           T(2), T(3), T(4), T(5),
        T(0), T(1), T(2), T(3), T(4)
    };

    std::vector<T> quvecrm(12);
    for(int i=0;i<12;i++) quvecrm[i] = quarrm[i];
    std::vector<T> qlvecrm(12);
    for(int i=0;i<12;i++) qlvecrm[i] = qlarrm[i];
    std::vector<T> quveccm(12);
    for(int i=0;i<12;i++) quveccm[i] = quarcm[i];
    std::vector<T> qlveccm(12);
    for(int i=0;i<12;i++) qlveccm[i] = qlarcm[i];
    std::vector<T> quvecdm(12);
    for(int i=0;i<12;i++) quvecdm[i] = quardm[i];
    std::vector<T> qlvecdm(12);
    for(int i=0;i<12;i++) qlvecdm[i] = qlardm[i];

    tmv::SymBandMatrix<T,U,S> q1(5,2);
    std::copy(quarrm, quarrm+12, q1.upperBand().rowmajor_begin());
    tmv::SymBandMatrix<T,U,S> q2(5,2);
    std::copy(quarcm, quarcm+12, q2.upperBand().colmajor_begin());
    tmv::SymBandMatrix<T,U,S> q3(5,2);
    std::copy(quardm, quardm+12, q3.upperBand().diagmajor_begin());

    tmv::SymBandMatrix<T,U,S> q4(5,2);
    std::copy(qlarrm, qlarrm+12, q4.lowerBand().rowmajor_begin());
    tmv::SymBandMatrix<T,U,S> q5(5,2);
    std::copy(qlarcm, qlarcm+12, q5.lowerBand().colmajor_begin());
    tmv::SymBandMatrix<T,U,S> q6(5,2);
    std::copy(qlardm, qlardm+12, q6.lowerBand().diagmajor_begin());

    tmv::SymBandMatrix<T,U,S> q7(5,2);
    std::copy(quvecrm.begin(), quvecrm.end(), q7.upperBand().rowmajor_begin());
    tmv::SymBandMatrix<T,U,S> q8(5,2);
    std::copy(quveccm.begin(), quveccm.end(), q8.upperBand().colmajor_begin());
    tmv::SymBandMatrix<T,U,S> q9(5,2);
    std::copy(quvecdm.begin(), quvecdm.end(), q9.upperBand().diagmajor_begin());

    tmv::SymBandMatrix<T,U,S> q10(5,2);
    std::copy(qlvecrm.begin(), qlvecrm.end(), q10.lowerBand().rowmajor_begin());
    tmv::SymBandMatrix<T,U,S> q11(5,2);
    std::copy(qlveccm.begin(), qlveccm.end(), q11.lowerBand().colmajor_begin());
    tmv::SymBandMatrix<T,U,S> q12(5,2);
    std::copy(qlvecdm.begin(), qlvecdm.end(), q12.lowerBand().diagmajor_begin());

    tmv::SymBandMatrix<T,U,S> q13x(50,20);
    tmv::SymBandMatrixView<T> q13 = q13x.subSymBandMatrix(3,28,2,5);
    std::copy(quvecrm.begin(), quvecrm.end(), q13.upperBand().rowmajor_begin());
    tmv::SymBandMatrix<T,U,S> q14x(50,20);
    tmv::SymBandMatrixView<T> q14 = q14x.subSymBandMatrix(3,28,2,5);
    std::copy(quveccm.begin(), quveccm.end(), q14.upperBand().colmajor_begin());
    tmv::SymBandMatrix<T,U,S> q15x(50,20);
    tmv::SymBandMatrixView<T> q15 = q15x.subSymBandMatrix(3,28,2,5);
    std::copy(quvecdm.begin(), quvecdm.end(), q15.upperBand().diagmajor_begin());

    tmv::SymBandMatrix<T,U,S> q16x(50,20);
    tmv::SymBandMatrixView<T> q16 = q16x.subSymBandMatrix(3,28,2,5);
    std::copy(qlvecrm.begin(), qlvecrm.end(), q16.lowerBand().rowmajor_begin());
    tmv::SymBandMatrix<T,U,S> q17x(50,20);
    tmv::SymBandMatrixView<T> q17 = q17x.subSymBandMatrix(3,28,2,5);
    std::copy(qlveccm.begin(), qlveccm.end(), q17.lowerBand().colmajor_begin());
    tmv::SymBandMatrix<T,U,S> q18x(50,20);
    tmv::SymBandMatrixView<T> q18 = q18x.subSymBandMatrix(3,28,2,5);
    std::copy(qlvecdm.begin(), qlvecdm.end(), q18.lowerBand().diagmajor_begin());

    // Assignment using op<< is always in rowmajor order.
    tmv::SymBandMatrix<T,U,S> q19(5,2);
    tmv::SymBandMatrix<T,U,S> q20t(5,2);
    tmv::SymBandMatrixView<T> q20 = q20t.transpose();

    tmv::SymBandMatrix<T,U,S> q21(5,2);
    tmv::SymBandMatrix<T,U,S> q22t(5,2);
    tmv::SymBandMatrixView<T> q22 = q22t.transpose();

    q19.upperBand() <<
        0, 2, 4,
           1, 3, 5,
              2, 4, 6,
                 3, 5, 
                    4;
    q20.upperBand() <<
        0, 2, 4,
           1, 3, 5,
              2, 4, 6,
                 3, 5, 
                    4;
    q21.lowerBand() <<
        0, 
        2, 1,
        4, 3, 2,
           5, 4, 3,
              6, 5, 4;
    q22.lowerBand() <<
        0, 
        2, 1,
        4, 3, 2,
           5, 4, 3,
              6, 5, 4;

    // Can also view memory directly
    T quarrmfull[] = {
        T(0), T(2), T(4),
        T(1), T(3), T(5),
        T(2), T(4), T(6),
        T(3), T(5), T(7),
        T(4)
    };
    T quarcmfull[] = {
                    T(0),
        T(3), T(2), T(1),
        T(4), T(3), T(2),
        T(5), T(4), T(3),
        T(6), T(5), T(4)
    };
    T quardmfull[] = {
        T(0), T(1), T(2), T(3), T(4),
        T(2), T(3), T(4), T(5), T(6),
        T(4), T(5), T(6)
    };
    T qlarrmfull[] = {
                    T(0),
        T(3), T(2), T(1),
        T(4), T(3), T(2),
        T(5), T(4), T(3),
        T(6), T(5), T(4)
    };
    T qlarcmfull[] = {
        T(0), T(2), T(4),
        T(1), T(3), T(5),
        T(2), T(4), T(6),
        T(3), T(5), T(7),
        T(4)
    };
    T qlardmfull[] = {
                    T(4), T(5), T(6),
        T(1), T(2), T(3), T(4), T(5),
        T(0), T(1), T(2), T(3), T(4)
    };
    T* qarfull = (
        (S == tmv::RowMajor) ?
        ( (U == tmv::Upper) ? quarrmfull : qlarrmfull ) :
        (S == tmv::ColMajor) ?
        ( (U == tmv::Upper) ? quarcmfull : qlarcmfull ) :
        ( (U == tmv::Upper) ? quardmfull : qlardmfull ) );
    T* qarfullx = qarfull + (S == tmv::DiagMajor && U == tmv::Lower ? 8 : 0);
    const int Si = 
        (S == tmv::RowMajor) ? 2 :
        (S == tmv::ColMajor) ? 1 :
        -4;
    const int Sj = 
        (S == tmv::RowMajor) ? 1 :
        (S == tmv::ColMajor) ? 2 :
        5;
    const tmv::ConstSymBandMatrixView<T> q23 =
        tmv::SymBandMatrixViewOf(qarfull,5,2,U,S);
    const tmv::ConstSymBandMatrixView<T> q24 =
        tmv::SymBandMatrixViewOf(qarfullx,5,2,U,Si,Sj);

    if (showacc) {
        std::cout<<"q1 = "<<q1<<std::endl;
        std::cout<<"q2 = "<<q2<<std::endl;
        std::cout<<"q3 = "<<q3<<std::endl;
        std::cout<<"q4 = "<<q4<<std::endl;
        std::cout<<"q5 = "<<q5<<std::endl;
        std::cout<<"q6 = "<<q6<<std::endl;
        std::cout<<"q7 = "<<q7<<std::endl;
        std::cout<<"q8 = "<<q8<<std::endl;
        std::cout<<"q9 = "<<q9<<std::endl;
        std::cout<<"q10 = "<<q10<<std::endl;
        std::cout<<"q11 = "<<q11<<std::endl;
        std::cout<<"q12 = "<<q12<<std::endl;
        std::cout<<"q13 = "<<q13<<std::endl;
        std::cout<<"q14 = "<<q14<<std::endl;
        std::cout<<"q15 = "<<q15<<std::endl;
        std::cout<<"q16 = "<<q16<<std::endl;
        std::cout<<"q17 = "<<q17<<std::endl;
        std::cout<<"q18 = "<<q18<<std::endl;
        std::cout<<"q19 = "<<q19<<std::endl;
        std::cout<<"q20 = "<<q20<<std::endl;
        std::cout<<"q21 = "<<q21<<std::endl;
        std::cout<<"q22 = "<<q22<<std::endl;
        std::cout<<"q23 = "<<q23<<std::endl;
        std::cout<<"q24 = "<<q24<<std::endl;
    }

    for(int i=0;i<5;i++) for(int j=0;j<5;j++) if (std::abs(i-j) <= 2) {
        T val = i >= j ? T(2*i-j) : T(2*j-i);
        Assert(q1(i,j) == val,"Create SymBandMatrix from T* rm upper");
        Assert(q2(i,j) == val,"Create SymBandMatrix from T* cm upper");
        Assert(q3(i,j) == val,"Create SymBandMatrix from T* dm upper");
        Assert(q4(i,j) == val,"Create SymBandMatrix from T* rm lower");
        Assert(q5(i,j) == val,"Create SymBandMatrix from T* cm lower");
        Assert(q6(i,j) == val,"Create SymBandMatrix from T* dm lower");
        Assert(q7(i,j) == val,"Create SymBandMatrix from vector rm upper");
        Assert(q8(i,j) == val,"Create SymBandMatrix from vector cm upper");
        Assert(q9(i,j) == val,"Create SymBandMatrix from vector dm upper");
        Assert(q10(i,j) == val,"Create SymBandMatrix from vector rm lower");
        Assert(q11(i,j) == val,"Create SymBandMatrix from vector cm lower");
        Assert(q12(i,j) == val,"Create SymBandMatrix from vector dm lower");
        Assert(q13(i,j) == val,"Create SymBandMatrixView from vector rm upper");
        Assert(q14(i,j) == val,"Create SymBandMatrixView from vector cm upper");
        Assert(q15(i,j) == val,"Create SymBandMatrixView from vector dm upper");
        Assert(q16(i,j) == val,"Create SymBandMatrixView from vector rm lower");
        Assert(q17(i,j) == val,"Create SymBandMatrixView from vector cm lower");
        Assert(q18(i,j) == val,"Create SymBandMatrixView from vector dm lower");
        Assert(q19(i,j) == val,"Create SymBandMatrix from << list upper");
        Assert(q20(i,j) == val,"Create SymBandMatrixView from << list upper");
        Assert(q21(i,j) == val,"Create SymBandMatrix from << list lower");
        Assert(q22(i,j) == val,"Create SymBandMatrixView from << list lower");
        Assert(q23(i,j) == val,"Create SymBandMatrixView of T* (S)");
        Assert(q24(i,j) == val,"Create SymBandMatrixView of T* (Si,Sj)");
    }

    // Test Basic Arithmetic
    tmv::SymBandMatrix<std::complex<T>,U,S> s1(N,noff);
    tmv::SymBandMatrix<std::complex<T>,U,S> s2(N,noff);
    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        if (i <= j && j <= i + noff) {
            s1(i,j) = value; 
        }
        if (j <= i && i <= j + noff) {
            s2(i,j) = value; 
        }
    }

    tmv::SymBandMatrix<std::complex<T>,U,S> s3(N,noff);
    s3 = s1+s2;

    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        if ( j <= i + noff && i <= j + noff) {
            Assert(s3(i,j) == s1(i,j)+s2(i,j),"Add SymBandMatrices1");
            Assert(s3(j,i) == s1(i,j)+s2(i,j),"Add SymBandMatrices2");
        }
    }

    tmv::SymBandMatrix<std::complex<T>,U,S,tmv::FortranStyle> s3f(N,noff);
    s3f = s3;
    Assert(s3f == s3,"Copy CStyle SymBandMatrix to FortranStyle");

    s3 = s1-s2;

    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        if ( j <= i + noff && i <= j + noff) {
            Assert(s3(i,j) == s1(i,j)-s2(i,j),"Subtract SymBandMatrices1");
            Assert(s3(j,i) == s1(i,j)-s2(i,j),"Subtract SymBandMatrices2");
        }
    }

    tmv::BandMatrix<std::complex<T> > bm1 = s1;
    Assert(bm1.nlo() == s1.nlo(),"SymBandMatrix -> BandMatrix (nlo)");
    Assert(bm1.nhi() == s1.nhi(),"SymBandMatrix -> BandMatrix (nhi)");
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        if ( j <= i + noff && i <= j + noff) {
            Assert(s1(i,j) == bm1(i,j),"SymBandMatrix -> BandMatrix");
        }
    }
    tmv::SymBandMatrix<std::complex<T>,U,S> bm1_sb(bm1,noff);
    Assert(s1 == bm1_sb,"BandMatrix -> SymBandMatrix");

    tmv::Matrix<std::complex<T> > m1 = s1;
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        if ( j <= i + noff && i <= j + noff) {
            Assert(s1(i,j) == m1(i,j),"SymBandMatrix -> Matrix");
        } else {
            Assert(m1(i,j) == std::complex<T>(0),"SymBandMatrix -> Matrix (0)");
        }
    }
    tmv::SymBandMatrix<std::complex<T>,U,S> m1_sb(m1,noff);
    Assert(s1 == m1_sb,"Matrix -> SymBandMatrix");

    tmv::SymMatrix<std::complex<T> > sm1 = s1;
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        if ( j <= i + noff && i <= j + noff) {
            Assert(s1(i,j) == sm1(i,j),"SymBandMatrix -> SymMatrix");
        } else {
            Assert(sm1(i,j) == std::complex<T>(0),"SymBandMatrix -> SymMatrix (0)");
        }
    }
    tmv::SymBandMatrix<std::complex<T>,U,S> sm1_sb(sm1,noff);
    Assert(s1 == sm1_sb,"SymMatrix -> SymBandMatrix");

    tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> sur = s1;
    Assert(s1==sur,"SymBandMatrix == SymBandMatrix<U,R>");
    Assert(s1.view()==sur.view(),"SymBandMatrix.view == SymBandMatrix<U,R>.view");
    Assert(s1.transpose()==sur.transpose(),
           "SymBandMatrix.transpose == SymBandMatrix<U,R>.transpose");
    Assert(s1.conjugate()==sur.conjugate(),
           "SymBandMatrix.conjugate == SymBandMatrix<U,R>.conjugate");
    Assert(s1.adjoint()==sur.adjoint(),
           "SymBandMatrix.adjoint == SymBandMatrix<U,R>.adjoint");
    Assert(s1.upperBand()==sur.upperBand(),
           "SymBandMatrix.upperBand == SymBandMatrix<U,R>.upperBand");
    Assert(s1.lowerBand()==sur.lowerBand(),
           "SymBandMatrix.lowerBand == SymBandMatrix<U,R>.lowerBand");
    Assert(s1.realPart()==sur.realPart(),
           "SymBandMatrix.real == SymBandMatrix<U,R>.real");
    Assert(s1.imagPart()==sur.imagPart(),
           "SymBandMatrix.imag == SymBandMatrix<U,R>.imag");
    Assert(s1.subMatrix(N/2,N/2+2,N/2+1,N/2+noff+1)==
           sur.subMatrix(N/2,N/2+2,N/2+1,N/2+noff+1),
           "SymMatrix.subMatrix1 == SymMatrix<U,R>.subMatrix1");
    Assert(s1.subMatrix(3,noff,0,3)==sur.subMatrix(3,noff,0,3),
           "SymMatrix.subMatrix2 == SymMatrix<U,R>.subMatrix2");
    Assert(s1.subSymMatrix(N/4,N/4+noff+1)==sur.subSymMatrix(N/4,N/4+noff+1),
           "SymMatrix.subSymMatrix == SymMatrix<U,R>.subSymMatrix");
    Assert(s1.subBandMatrix(0,N/2,1,N/2+3,1,2)==
           sur.subBandMatrix(0,N/2,1,N/2+3,1,2),
           "SymMatrix.subBandMatrix == SymMatrix<U,R>.subBandMatrix");
    Assert(s1.subSymBandMatrix(N/4,3*N/4)==sur.subSymBandMatrix(N/4,3*N/4),
           "SymMatrix.subSymBandMatrix == SymMatrix<U,R>.subSymBandMatrix");
}

template <class T, tmv::UpLoType U, tmv::StorageType S>
static void TestBasicHermBandMatrix_2()
{
    const int N = 10;
    const int noff = 3;

    if (showstartdone) {
        std::cout<<"Start TestBasicHermBandMatrix_2\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"U = "<<tmv::TMV_Text(U)<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
        std::cout<<"noff = "<<noff<<std::endl;
    }

    // Test assignments and constructors from arrays
    typedef std::complex<T> CT;
    
    const CT quarrm[] = {
        CT(0,0), CT(2,1), CT(4,2),
                 CT(1,0), CT(3,1), CT(5,2),
                          CT(2,0), CT(4,1), CT(6,2),
                                   CT(3,0), CT(5,1),
                                            CT(4,0)
    };
    const CT qlarrm[] = {
        CT(0,0),
        CT(2,-1), CT(1,0),
        CT(4,-2), CT(3,-1), CT(2,0),
                  CT(5,-2), CT(4,-1), CT(3,0),
                            CT(6,-2), CT(5,-1), CT(4,0)
    };
    const CT quarcm[] = {
        CT(0,0),
        CT(2,1), CT(1,0),
        CT(4,2), CT(3,1), CT(2,0),
                 CT(5,2), CT(4,1), CT(3,0),
                          CT(6,2), CT(5,1), CT(4,0)
    };
    const CT qlarcm[] = {
        CT(0,0), CT(2,-1), CT(4,-2),
                 CT(1,0),  CT(3,-1), CT(5,-2),
                           CT(2,0),  CT(4,-1), CT(6,-2),
                                     CT(3,0),  CT(5,-1),
                                               CT(4,0)
    };
    const CT quardm[] = {
        CT(0,0), CT(1,0), CT(2,0), CT(3,0), CT(4,0),
            CT(2,1), CT(3,1), CT(4,1), CT(5,1),
                 CT(4,2), CT(5,2), CT(6,2)
    };
    const CT qlardm[] = {
                  CT(4,-2), CT(5,-2), CT(6,-2),
             CT(2,-1), CT(3,-1), CT(4,-1), CT(5,-1),
        CT(0,0),  CT(1,0),  CT(2,0),  CT(3,0),  CT(4,0)
    };

    std::vector<CT> quvecrm(12);
    for(int i=0;i<12;i++) quvecrm[i] = quarrm[i];
    std::vector<CT> qlvecrm(12);
    for(int i=0;i<12;i++) qlvecrm[i] = qlarrm[i];
    std::vector<CT> quveccm(12);
    for(int i=0;i<12;i++) quveccm[i] = quarcm[i];
    std::vector<CT> qlveccm(12);
    for(int i=0;i<12;i++) qlveccm[i] = qlarcm[i];
    std::vector<CT> quvecdm(12);
    for(int i=0;i<12;i++) quvecdm[i] = quardm[i];
    std::vector<CT> qlvecdm(12);
    for(int i=0;i<12;i++) qlvecdm[i] = qlardm[i];

    tmv::HermBandMatrix<CT,U,S> q1(5,2);
    std::copy(quarrm, quarrm+12, q1.upperBand().rowmajor_begin());
    tmv::HermBandMatrix<CT,U,S> q2(5,2);
    std::copy(quarcm, quarcm+12, q2.upperBand().colmajor_begin());
    tmv::HermBandMatrix<CT,U,S> q3(5,2);
    std::copy(quardm, quardm+12, q3.upperBand().diagmajor_begin());

    tmv::HermBandMatrix<CT,U,S> q4(5,2);
    std::copy(qlarrm, qlarrm+12, q4.lowerBand().rowmajor_begin());
    tmv::HermBandMatrix<CT,U,S> q5(5,2);
    std::copy(qlarcm, qlarcm+12, q5.lowerBand().colmajor_begin());
    tmv::HermBandMatrix<CT,U,S> q6(5,2);
    std::copy(qlardm, qlardm+12, q6.lowerBand().diagmajor_begin());

    tmv::HermBandMatrix<CT,U,S> q7(5,2);
    std::copy(quvecrm.begin(), quvecrm.end(), q7.upperBand().rowmajor_begin());
    tmv::HermBandMatrix<CT,U,S> q8(5,2);
    std::copy(quveccm.begin(), quveccm.end(), q8.upperBand().colmajor_begin());
    tmv::HermBandMatrix<CT,U,S> q9(5,2);
    std::copy(quvecdm.begin(), quvecdm.end(), q9.upperBand().diagmajor_begin());

    tmv::HermBandMatrix<CT,U,S> q10(5,2);
    std::copy(qlvecrm.begin(), qlvecrm.end(), q10.lowerBand().rowmajor_begin());
    tmv::HermBandMatrix<CT,U,S> q11(5,2);
    std::copy(qlveccm.begin(), qlveccm.end(), q11.lowerBand().colmajor_begin());
    tmv::HermBandMatrix<CT,U,S> q12(5,2);
    std::copy(qlvecdm.begin(), qlvecdm.end(), q12.lowerBand().diagmajor_begin());

    tmv::HermBandMatrix<CT,U,S> q13x(50,20);
    tmv::SymBandMatrixView<CT> q13 = q13x.subSymBandMatrix(3,28,2,5);
    std::copy(quvecrm.begin(), quvecrm.end(), q13.upperBand().rowmajor_begin());
    tmv::HermBandMatrix<CT,U,S> q14x(50,20);
    tmv::SymBandMatrixView<CT> q14 = q14x.subSymBandMatrix(3,28,2,5);
    std::copy(quveccm.begin(), quveccm.end(), q14.upperBand().colmajor_begin());
    tmv::HermBandMatrix<CT,U,S> q15x(50,20);
    tmv::SymBandMatrixView<CT> q15 = q15x.subSymBandMatrix(3,28,2,5);
    std::copy(quvecdm.begin(), quvecdm.end(), q15.upperBand().diagmajor_begin());

    tmv::HermBandMatrix<CT,U,S> q16x(50,20);
    tmv::SymBandMatrixView<CT> q16 = q16x.subSymBandMatrix(3,28,2,5);
    std::copy(qlvecrm.begin(), qlvecrm.end(), q16.lowerBand().rowmajor_begin());
    tmv::HermBandMatrix<CT,U,S> q17x(50,20);
    tmv::SymBandMatrixView<CT> q17 = q17x.subSymBandMatrix(3,28,2,5);
    std::copy(qlveccm.begin(), qlveccm.end(), q17.lowerBand().colmajor_begin());
    tmv::HermBandMatrix<CT,U,S> q18x(50,20);
    tmv::SymBandMatrixView<CT> q18 = q18x.subSymBandMatrix(3,28,2,5);
    std::copy(qlvecdm.begin(), qlvecdm.end(), q18.lowerBand().diagmajor_begin());

    // Assignment using op<< is always in rowmajor order.
    tmv::HermBandMatrix<CT,U,S> q19(5,2);
    tmv::HermBandMatrix<CT,U,S> q20t(5,2);
    tmv::SymBandMatrixView<CT> q20 = q20t.transpose();

    tmv::HermBandMatrix<CT,U,S> q21(5,2);
    tmv::HermBandMatrix<CT,U,S> q22t(5,2);
    tmv::SymBandMatrixView<CT> q22 = q22t.transpose();

    q19.upperBand() <<
        CT(0,0), CT(2,1), CT(4,2),
                 CT(1,0), CT(3,1), CT(5,2),
                          CT(2,0), CT(4,1), CT(6,2),
                                   CT(3,0), CT(5,1),
                                            CT(4,0);
    q20.upperBand() <<
        CT(0,0), CT(2,1), CT(4,2),
                 CT(1,0), CT(3,1), CT(5,2),
                          CT(2,0), CT(4,1), CT(6,2),
                                   CT(3,0), CT(5,1),
                                            CT(4,0);
    q21.lowerBand() <<
        CT(0,0),
        CT(2,-1), CT(1,0),
        CT(4,-2), CT(3,-1), CT(2,0),
                  CT(5,-2), CT(4,-1), CT(3,0),
                            CT(6,-2), CT(5,-1), CT(4,0);
    q22.lowerBand() <<
        CT(0,0),
        CT(2,-1), CT(1,0),
        CT(4,-2), CT(3,-1), CT(2,0),
                  CT(5,-2), CT(4,-1), CT(3,0),
                            CT(6,-2), CT(5,-1), CT(4,0);

    // Can also view memory directly
    CT quarrmfull[] = {
        CT(0,0), CT(2,1), CT(4,2),
        CT(1,0), CT(3,1), CT(5,2),
        CT(2,0), CT(4,1), CT(6,2),
        CT(3,0), CT(5,1), CT(7,2),
        CT(4,0)
    };
    CT quarcmfull[] = {
                          CT(0,0),
        CT(3,2), CT(2,1), CT(1,0),
        CT(4,2), CT(3,1), CT(2,0),
        CT(5,2), CT(4,1), CT(3,0),
        CT(6,2), CT(5,1), CT(4,0)
    };
    CT quardmfull[] = {
        CT(0,0), CT(1,0), CT(2,0), CT(3,0), CT(4,0),
        CT(2,1), CT(3,1), CT(4,1), CT(5,1), CT(6,1),
        CT(4,2), CT(5,2), CT(6,2)
    };
    CT qlarrmfull[] = {
                            CT(0,0),
        CT(3,-2), CT(2,-1), CT(1,0),
        CT(4,-2), CT(3,-1), CT(2,0),
        CT(5,-2), CT(4,-1), CT(3,0),
        CT(6,-2), CT(5,-1), CT(4,0)
    };
    CT qlarcmfull[] = {
        CT(0,0), CT(2,-1), CT(4,-2),
        CT(1,0), CT(3,-1), CT(5,-2),
        CT(2,0), CT(4,-1), CT(6,-2),
        CT(3,0), CT(5,-1), CT(7,-2),
        CT(4,0)
    };
    CT qlardmfull[] = {
                            CT(4,-2), CT(5,-2), CT(6,-2),
        CT(1,-1), CT(2,-1), CT(3,-1), CT(4,-1), CT(5,-1),
        CT(0,0),  CT(1,0),  CT(2,0),  CT(3,0),  CT(4,0)
    };
    CT* qarfull = (
        (S == tmv::RowMajor) ?
        ( (U == tmv::Upper) ? quarrmfull : qlarrmfull ) :
        (S == tmv::ColMajor) ?
        ( (U == tmv::Upper) ? quarcmfull : qlarcmfull ) :
        ( (U == tmv::Upper) ? quardmfull : qlardmfull ) );
    CT* qarfullx = qarfull + (S == tmv::DiagMajor && U == tmv::Lower ? 8 : 0);
    const int Si = 
        (S == tmv::RowMajor) ? 2 :
        (S == tmv::ColMajor) ? 1 :
        -4;
    const int Sj = 
        (S == tmv::RowMajor) ? 1 :
        (S == tmv::ColMajor) ? 2 :
        5;
    const tmv::ConstSymBandMatrixView<CT> q23 =
        tmv::HermBandMatrixViewOf(qarfull,5,2,U,S);
    const tmv::ConstSymBandMatrixView<CT> q24 =
        tmv::HermBandMatrixViewOf(qarfullx,5,2,U,Si,Sj);

    if (showacc) {
        std::cout<<"q1 = "<<q1<<std::endl;
        std::cout<<"q2 = "<<q2<<std::endl;
        std::cout<<"q3 = "<<q3<<std::endl;
        std::cout<<"q4 = "<<q4<<std::endl;
        std::cout<<"q5 = "<<q5<<std::endl;
        std::cout<<"q6 = "<<q6<<std::endl;
        std::cout<<"q7 = "<<q7<<std::endl;
        std::cout<<"q8 = "<<q8<<std::endl;
        std::cout<<"q9 = "<<q9<<std::endl;
        std::cout<<"q10 = "<<q10<<std::endl;
        std::cout<<"q11 = "<<q11<<std::endl;
        std::cout<<"q12 = "<<q12<<std::endl;
        std::cout<<"q13 = "<<q13<<std::endl;
        std::cout<<"q14 = "<<q14<<std::endl;
        std::cout<<"q15 = "<<q15<<std::endl;
        std::cout<<"q16 = "<<q16<<std::endl;
        std::cout<<"q17 = "<<q17<<std::endl;
        std::cout<<"q18 = "<<q18<<std::endl;
        std::cout<<"q19 = "<<q19<<std::endl;
        std::cout<<"q20 = "<<q20<<std::endl;
        std::cout<<"q21 = "<<q21<<std::endl;
        std::cout<<"q22 = "<<q22<<std::endl;
        std::cout<<"q23 = "<<q23<<std::endl;
        std::cout<<"q24 = "<<q24<<std::endl;
    }

    for(int i=0;i<5;i++) for(int j=0;j<5;j++) if (std::abs(i-j) <= 2) {
        CT val = i >= j ? CT(2*i-j,j-i) : CT(2*j-i,j-i);
        Assert(q1(i,j) == val,"Create HermBandMatrix from T* rm upper");
        Assert(q2(i,j) == val,"Create HermBandMatrix from T* cm upper");
        Assert(q3(i,j) == val,"Create HermBandMatrix from T* dm upper");
        Assert(q4(i,j) == val,"Create HermBandMatrix from T* rm lower");
        Assert(q5(i,j) == val,"Create HermBandMatrix from T* cm lower");
        Assert(q6(i,j) == val,"Create HermBandMatrix from T* dm lower");
        Assert(q7(i,j) == val,"Create HermBandMatrix from vector rm upper");
        Assert(q8(i,j) == val,"Create HermBandMatrix from vector cm upper");
        Assert(q9(i,j) == val,"Create HermBandMatrix from vector dm upper");
        Assert(q10(i,j) == val,"Create HermBandMatrix from vector rm lower");
        Assert(q11(i,j) == val,"Create HermBandMatrix from vector cm lower");
        Assert(q12(i,j) == val,"Create HermBandMatrix from vector dm lower");
        Assert(q13(i,j) == val,"Create HermBandMatrixView from vector rm upper");
        Assert(q14(i,j) == val,"Create HermBandMatrixView from vector cm upper");
        Assert(q15(i,j) == val,"Create HermBandMatrixView from vector dm upper");
        Assert(q16(i,j) == val,"Create HermBandMatrixView from vector rm lower");
        Assert(q17(i,j) == val,"Create HermBandMatrixView from vector cm lower");
        Assert(q18(i,j) == val,"Create HermBandMatrixView from vector dm lower");
        Assert(q19(i,j) == val,"Create HermBandMatrix from << list upper");
        Assert(q20(i,j) == val,"Create HermBandMatrixView from << list upper");
        Assert(q21(i,j) == val,"Create HermBandMatrix from << list lower");
        Assert(q22(i,j) == val,"Create HermBandMatrixView from << list lower");
        Assert(q23(i,j) == val,"Create HermBandMatrixView of T* (S)");
        Assert(q24(i,j) == val,"Create HermBandMatrixView of T* (Si,Sj)");
    }

    // Test Basic Arithmetic
    tmv::HermBandMatrix<std::complex<T>,U,S> h1(N,noff);
    tmv::HermBandMatrix<std::complex<T>,U,S> h2(N,noff);
    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (i <= j && j <= i + noff) {
            h1(i,j) = hvalue; 
        }
        if (j <= i && i <= j + noff) {
            h2(i,j) = hvalue; 
        }
    }

    tmv::HermBandMatrix<std::complex<T>,U,S> h3(N,noff);
    h3 = h1+h2;

    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        if ( j <= i + noff && i <= j + noff) {
            Assert(h3(i,j) == h1(i,j)+h2(i,j),"Add HermBandMatrices1");
            Assert(h3(j,i) == conj(h1(i,j)+h2(i,j)),"Add HermBandMatrices2");
        }
    }

    tmv::HermBandMatrix<std::complex<T>,U,S,tmv::FortranStyle> h3f(N,noff);
    h3f = h3;
    Assert(h3f == h3,"Copy CStyle HermBandMatrix to FortranStyle");

    h3 = h1-h2;

    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        if ( j <= i + noff && i <= j + noff) {
            Assert(h3(i,j) == h1(i,j)-h2(i,j),"Subtract HermBandMatrices1");
            Assert(h3(j,i) == conj(h1(i,j)-h2(i,j)),
                   "Subtract HermBandMatrices2");
        }
    }

    tmv::BandMatrix<std::complex<T> > bn1 = h1;
    Assert(bn1.nlo() == h1.nlo(),"SymBandMatrix -> BandMatrix (nlo)");
    Assert(bn1.nhi() == h1.nhi(),"SymBandMatrix -> BandMatrix (nhi)");
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        if ( j <= i + noff && i <= j + noff) {
            Assert(h1(i,j) == bn1(i,j),"HermBandMatrix -> BandMatrix");
        }
    }
    tmv::HermBandMatrix<std::complex<T>,U,S> bn1_hb(bn1,noff);
    Assert(h1 == bn1_hb,"BandMatrix -> HermBandMatrix");

    tmv::Matrix<std::complex<T> > n1 = h1;
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        if ( j <= i + noff && i <= j + noff) {
            Assert(h1(i,j) == n1(i,j),"HermBandMatrix -> Matrix");
        } else {
            Assert(n1(i,j) == std::complex<T>(0),"HermBandMatrix -> Matrix (0)");
        }
    }
    tmv::HermBandMatrix<std::complex<T>,U,S> n1_hb(n1,noff);
    Assert(h1 == n1_hb,"Matrix -> HermBandMatrix");

    tmv::HermMatrix<std::complex<T> > hn1 = h1;
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        if ( j <= i + noff && i <= j + noff) {
            Assert(h1(i,j) == hn1(i,j),"HermBandMatrix -> HermMatrix");
        } else {
            Assert(hn1(i,j) == std::complex<T>(0),"HermBandMatrix -> HermMatrix (0)");
        }
    }
    tmv::HermBandMatrix<std::complex<T>,U,S> hn1_hb(hn1,noff);
    Assert(h1 == hn1_hb,"HermMatrix -> HermBandMatrix");

    tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> hur = h1;
    Assert(h1==hur,"HermMatrix == HermMatrix<U,R>");
    Assert(h1.view()==hur.view(),"HermMatrix.view == HermMatrix<U,R>.view");
    Assert(h1.transpose()==hur.transpose(),
           "HermMatrix.transpose == HermMatrix<U,R>.transpose");
    Assert(h1.conjugate()==hur.conjugate(),
           "HermMatrix.conjugate == HermMatrix<U,R>.conjugate");
    Assert(h1.adjoint()==hur.adjoint(),
           "HermMatrix.adjoint == HermMatrix<U,R>.adjoint");
    Assert(h1.upperBand()==hur.upperBand(),
           "HermMatrix.upperBand == HermMatrix<U,R>.upperBand");
    Assert(h1.lowerBand()==hur.lowerBand(),
           "HermMatrix.lowerBand == HermMatrix<U,R>.lowerBand");
    Assert(h1.realPart()==hur.realPart(),
           "HermMatrix.real == HermMatrix<U,R>.real");
    Assert(h1.subMatrix(N/2,N/2+2,N/2+1,N/2+noff+1)==
           hur.subMatrix(N/2,N/2+2,N/2+1,N/2+noff+1),
           "SymMatrix.subMatrix1 == SymMatrix<U,R>.subMatrix1");
    Assert(h1.subMatrix(3,noff,0,3)==hur.subMatrix(3,noff,0,3),
           "SymMatrix.subMatrix2 == SymMatrix<U,R>.subMatrix2");
    Assert(h1.subSymMatrix(N/4,N/4+noff+1)==hur.subSymMatrix(N/4,N/4+noff+1),
           "SymMatrix.subSymMatrix == SymMatrix<U,R>.subSymMatrix");
    Assert(h1.subBandMatrix(0,N/2,1,N/2+3,1,2)==
           hur.subBandMatrix(0,N/2,1,N/2+3,1,2),
           "SymMatrix.subBandMatrix == SymMatrix<U,R>.subBandMatrix");
    Assert(h1.subSymBandMatrix(N/4,3*N/4)==hur.subSymBandMatrix(N/4,3*N/4),
           "SymMatrix.subSymBandMatrix == SymMatrix<U,R>.subSymBandMatrix");
}

template <class T, tmv::UpLoType U, tmv::StorageType S>
static void TestBasicSymBandMatrix_IO()
{
    const int N = 10;
    const int noff = 3;

    if (showstartdone) {
        std::cout<<"Start TestBasicSymBandMatrix_IO\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"U = "<<tmv::TMV_Text(U)<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
        std::cout<<"noff = "<<noff<<std::endl;
    }

    tmv::SymBandMatrix<std::complex<T>,U,S> s1(N,noff);
    tmv::SymBandMatrix<std::complex<T>,U,S> s2(N,noff);
    tmv::HermBandMatrix<std::complex<T>,U,S> h1(N,noff);
    tmv::HermBandMatrix<std::complex<T>,U,S> h2(N,noff);
    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (i <= j && j <= i + noff) {
            s1(i,j) = value; 
            h1(i,j) = hvalue; 
        }
        if (j <= i && i <= j + noff) {
            s2(i,j) = value; 
            h2(i,j) = hvalue; 
        }
    }

    std::ofstream fout("tmvtest_symbandmatrix_io.dat");
    if (!fout) 
#ifdef NOTHROW
    { std::cerr<<"Couldn't open tmvtest_symbandmatrix_io.dat for output\n"; exit(1); }
#else
    throw std::runtime_error(
        "Couldn't open tmvtest_symbandmatrix_io.dat for output");
#endif
    fout << s1 << std::endl << h1 << std::endl;
    s1.writeCompact(fout);
    h1.writeCompact(fout);
    fout.close();

    tmv::Matrix<std::complex<T>,tmv::RowMajor> xsm1(N,N);
    tmv::Matrix<std::complex<T>,tmv::RowMajor> xhm1(N,N);
    tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> xs1(N,noff);
    tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> xh1(N,noff);
    std::ifstream fin("tmvtest_symbandmatrix_io.dat");
    if (!fin) 
#ifdef NOTHROW
    { std::cerr<<"Couldn't open tmvtest_symbandmatrix_io.dat for input\n"; exit(1); }
#else
    throw std::runtime_error(
        "Couldn't open tmvtest_symbandmatrix_io.dat for input");
#endif
    fin >> xsm1 >> xhm1;
    fin >> xs1 >> xh1;
    fin.close();
    Assert(tmv::Matrix<std::complex<T> >(s1) == xsm1,"SymBandMatrix I/O check #1");
    Assert(tmv::Matrix<std::complex<T> >(h1) == xhm1,"HermBandMatrix I/O check #1");
    Assert(s1 == xs1,"SymBandMatrix Compact I/O check #1");
    Assert(h1 == xh1,"HermBandMatrix Compact I/O check #1");

    tmv::Matrix<std::complex<T>,tmv::ColMajor> xsm2(N,N);
    tmv::Matrix<std::complex<T>,tmv::ColMajor> xhm2(N,N);
    tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> xs2(N,noff);
    tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> xh2(N,noff);
    fin.open("tmvtest_symbandmatrix_io.dat");
    if (!fin) 
#ifdef NOTHROW
    { std::cerr<<"Couldn't open tmvtest_symbandmatrix_io.dat for input\n"; exit(1); }
#else
    throw std::runtime_error(
        "Couldn't open tmvtest_symbandmatrix_io.dat for input");
#endif
    fin >> xsm2 >> xhm2;
    fin >> xs2 >> xh2;
    fin.close();
    Assert(tmv::Matrix<std::complex<T> >(s1) == xsm2,"SymBandMatrix I/O check #2");
    Assert(tmv::Matrix<std::complex<T> >(h1) == xhm2,"HermBandMatrix I/O check #2");
    Assert(s1 == xs2,"SymBandMatrix Compact I/O check #2");
    Assert(h1 == xh2,"HermBandMatrix Compact I/O check #2");

    tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> xs3(N,noff);
    tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> xh3(N,noff);
    fin.open("tmvtest_symbandmatrix_io.dat");
    if (!fin) 
#ifdef NOTHROW
    { std::cerr<<"Couldn't open tmvtest_symbandmatrix_io.dat for input\n"; exit(1); }
#else
    throw std::runtime_error(
        "Couldn't open tmvtest_symbandmatrix_io.dat for input");
#endif
    fin >> xsm1 >> xhm1;
    fin >> xs3 >> xh3;
    fin.close();
    Assert(s1 == xs3,"SymBandMatrix Compact I/O check #3");
    Assert(h1 == xh3,"HermBandMatrix Compact I/O check #3");

    tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> xs4(N,noff);
    tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> xh4(N,noff);
    fin.open("tmvtest_symbandmatrix_io.dat");
    if (!fin) 
#ifdef NOTHROW
    { std::cerr<<"Couldn't open tmvtest_symbandmatrix_io.dat for input\n"; exit(1); }
#else
    throw std::runtime_error(
        "Couldn't open tmvtest_symbandmatrix_io.dat for input");
#endif
    fin >> xsm1 >> xhm1;
    fin >> xs4 >> xh4;
    fin.close();
    Assert(s1 == xs4,"SymBandMatrix Compact I/O check #4");
    Assert(h1 == xh4,"HermBandMatrix Compact I/O check #4");

    std::auto_ptr<tmv::Matrix<std::complex<T> > > xsm5;
    std::auto_ptr<tmv::Matrix<std::complex<T> > > xhm5;
    std::auto_ptr<tmv::SymBandMatrix<std::complex<T> > > xs5;
    std::auto_ptr<tmv::HermBandMatrix<std::complex<T> > > xh5;
    fin.open("tmvtest_symbandmatrix_io.dat");
    if (!fin) 
#ifdef NOTHROW
    { std::cerr<<"Couldn't open tmvtest_symbandmatrix_io.dat for input\n"; exit(1); }
#else
    throw std::runtime_error(
        "Couldn't open tmvtest_symbandmatrix_io.dat for input");
#endif
    fin >> xsm5 >> xhm5;
    fin >> xs5 >> xh5;
    fin.close();
    Assert(tmv::Matrix<std::complex<T> >(s1) == *xsm5,
           "SymBandMatrix I/O check #5");
    Assert(tmv::Matrix<std::complex<T> >(h1) == *xhm5,
           "HermBandMatrix I/O check #5");
    Assert(s1 == *xs5,"SymBandMatrix Compact I/O check #5");
    Assert(h1 == *xh5,"HermBandMatrix Compact I/O check #5");

    std::remove("tmvtest_symbandmatrix_io.dat");
}

template <class T, tmv::UpLoType U, tmv::StorageType S>
static void TestBasicSymBandMatrix()
{
    TestBasicSymBandMatrix_1<T,U,S>();
    TestBasicSymBandMatrix_2<T,U,S>();
    TestBasicHermBandMatrix_1<T,U,S>();
    TestBasicHermBandMatrix_2<T,U,S>();
    TestBasicSymBandMatrix_IO<T,U,S>();
}

template <class T> 
void TestSymBandMatrix() 
{
    TestBasicSymBandMatrix<T,tmv::Upper,tmv::ColMajor>();
    TestBasicSymBandMatrix<T,tmv::Lower,tmv::ColMajor>();
#if (XTEST & 2)
    TestBasicSymBandMatrix<T,tmv::Upper,tmv::RowMajor>();
    TestBasicSymBandMatrix<T,tmv::Lower,tmv::RowMajor>();
#endif

    std::cout<<"SymBandMatrix<"<<tmv::TMV_Text(T())<<
        "> passed all basic tests\n";

    TestSymBandMatrixArith_A<T>();
    std::cout<<"SymBandMatrix<"<<tmv::TMV_Text(T())<<
        "> (SymBand/SymBand) Arithmetic passed all tests\n";
    TestSymBandMatrixArith_B1<T>();
    TestSymBandMatrixArith_B2<T>();
    std::cout<<"SymBandMatrix<"<<tmv::TMV_Text(T())<<
        "> (Matrix/SymBand) Arithmetic passed all tests\n";
    TestSymBandMatrixArith_C1<T>();
    TestSymBandMatrixArith_C2<T>();
    std::cout<<"SymBandMatrix<"<<tmv::TMV_Text(T())<<
        "> (Diag/SymBand) Arithmetic passed all tests\n";
    TestSymBandMatrixArith_D1<T>();
    TestSymBandMatrixArith_D2<T>();
    std::cout<<"SymBandMatrix<"<<tmv::TMV_Text(T())<<
        "> (Tri/SymBand) Arithmetic passed all tests\n";
    TestSymBandMatrixArith_E1<T>();
    TestSymBandMatrixArith_E2<T>();
    std::cout<<"SymBandMatrix<"<<tmv::TMV_Text(T())<<
        "> (Band/SymBand) Arithmetic passed all tests\n";
    TestSymBandMatrixArith_F1<T>();
    TestSymBandMatrixArith_F2<T>();
    std::cout<<"SymBandMatrix<"<<tmv::TMV_Text(T())<<
        "> (Sym/SymBand) Arithmetic passed all tests\n";
}

#ifdef TEST_DOUBLE
template void TestSymBandMatrix<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymBandMatrix<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandMatrix<long double>();
#endif
#ifdef TEST_INT
template void TestSymBandMatrix<int>();
#endif
