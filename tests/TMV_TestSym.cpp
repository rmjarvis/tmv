
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include <fstream>
#include <cstdio>

#define CT std::complex<T>

template <class T, tmv::UpLoType U, tmv::StorageType S> 
static void TestBasicSymMatrix_1()
{
    const int N = 6;

    if (showstartdone) {
        std::cout<<"Start TestBasicSymMatrix_1\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"U = "<<tmv::TMV_Text(U)<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
    }

    tmv::SymMatrix<std::complex<T>,U|S> s1(N);
    tmv::SymMatrix<std::complex<T>,U|S> s2(N);
    tmv::SymMatrix<std::complex<T>,U|S|tmv::FortranStyle> s1f(N);
    tmv::SymMatrix<std::complex<T>,U|S|tmv::FortranStyle> s2f(N);

    Assert(int(s1.colsize()) == N && int(s1.rowsize()) == N,
           "Creating SymMatrix(N)");

    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        if (i<=j) {
            s1(i,j) = value; 
            s1f(i+1,j+1) = value; 
        }
        if (j<=i) {
            s2(i,j) = value; 
            s2f(i+1,j+1) = value; 
        }
    }

    tmv::SymMatrixView<std::complex<T> > s1v = s1.view();
    tmv::ConstSymMatrixView<std::complex<T> > s1cv = s1.view();
    tmv::SymMatrixView<std::complex<T> > s2v = s2.view();
    tmv::ConstSymMatrixView<std::complex<T> > s2cv = s2.view();
    tmv::SymMatrixView<std::complex<T>,tmv::FortranStyle> s1fv = s1f.view();
    tmv::ConstSymMatrixView<std::complex<T>,tmv::FortranStyle> s1fcv =
        s1f.view();
    tmv::SymMatrixView<std::complex<T>,tmv::FortranStyle> s2fv = s2f.view();
    tmv::ConstSymMatrixView<std::complex<T>,tmv::FortranStyle> s2fcv =
        s2f.view();
    const tmv::SymMatrix<std::complex<T>,U|S>& s1x = s1;
    const tmv::SymMatrix<std::complex<T>,U|S>& s2x = s2;
    const tmv::SymMatrix<std::complex<T>,U|S|tmv::FortranStyle>& s1fx = s1f;
    const tmv::SymMatrix<std::complex<T>,U|S|tmv::FortranStyle>& s2fx = s2f;

    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        if (showtests) std::cout<<"i,j = "<<i<<','<<j<<std::endl;
        if (i<=j) {
            Assert(s1(i,j) == value,"Read/Write SymMatrix");
            Assert(s1x(i,j) == value,"Access const SymMatrix");
            Assert(s1v(i,j) == value,"Access SymMatrix V");
            Assert(s1cv(i,j) == value,"Access SymMatrix CV");
            Assert(s1(j,i) == value,"Access SymMatrix - opposite tri");
            Assert(s1x(j,i) == value,"Access const SymMatrix - opposite tri");
            Assert(s1v(j,i) == value,"Access SymMatrix V");
            Assert(s1cv(j,i) == value,"Access SymMatrix CV");
            Assert(s1f(i+1,j+1) == value,"Read/Write SymMatrixF");
            Assert(s1fx(i+1,j+1) == value,"Access const SymMatrixF");
            Assert(s1fv(i+1,j+1) == value,"Access SymMatrixF V");
            Assert(s1fcv(i+1,j+1) == value,"Access SymMatrixF CV");
            Assert(s1f(j+1,i+1) == value,"Access SymMatrixF - opposite tri");
            Assert(s1fx(j+1,i+1) == value,
                   "Access const SymMatrixF - opposite tri");
            Assert(s1fv(j+1,i+1) == value,"Access SymMatrixF V");
            Assert(s1fcv(j+1,i+1) == value,"Access SymMatrixF CV");
            Assert(s1.row(i,i,N)(j-i) == value,"SymMatrix.row");
            Assert(s1x.row(i,i,N)(j-i) == value,"SymMatrix.row");
            Assert(s1cv.row(i,i,N)(j-i) == value,"SymMatrix.row CV");
            Assert(s1v.row(i,i,N)(j-i) == value,"SymMatrix.row V");
            Assert(s1.col(i,i,N)(j-i) == value,"SymMatrix.col");
            Assert(s1x.col(i,i,N)(j-i) == value,"SymMatrix.col");
            Assert(s1cv.col(i,i,N)(j-i) == value,"SymMatrix.col CV");
            Assert(s1v.col(i,i,N)(j-i) == value,"SymMatrix.col V");
            Assert(s1.row(i,j,N)(0) == value,"SymMatrix.row2");
            Assert(s1x.row(i,j,N)(0) == value,"SymMatrix.row2");
            Assert(s1cv.row(i,j,N)(0) == value,"SymMatrix.row2 CV");
            Assert(s1v.row(i,j,N)(0) == value,"SymMatrix.row2 V");
            Assert(s1.col(i,j,N)(0) == value,"SymMatrix.col2");
            Assert(s1x.col(i,j,N)(0) == value,"SymMatrix.col2");
            Assert(s1cv.col(i,j,N)(0) == value,"SymMatrix.col2 CV");
            Assert(s1v.col(i,j,N)(0) == value,"SymMatrix.col2 V");
            Assert(s1f.row(i+1,i+1,N)(j-i+1) == value,"SymMatrixF.row");
            Assert(s1fx.row(i+1,i+1,N)(j-i+1) == value,"const SymMatrixF.row");
            Assert(s1fcv.row(i+1,i+1,N)(j-i+1) == value,"SymMatrixF.row CV");
            Assert(s1fv.row(i+1,i+1,N)(j-i+1) == value,"SymMatrixF.row V");
            Assert(s1f.col(i+1,i+1,N)(j-i+1) == value,"SymMatrixF.col");
            Assert(s1fx.col(i+1,i+1,N)(j-i+1) == value,"const SymMatrixF.col");
            Assert(s1fcv.col(i+1,i+1,N)(j-i+1) == value,"SymMatrixF.col CV");
            Assert(s1fv.col(i+1,i+1,N)(j-i+1) == value,"SymMatrixF.col V");
            Assert(s1f.row(i+1,j+1,N)(1) == value,"SymMatrixF.row2");
            Assert(s1fx.row(i+1,j+1,N)(1) == value,"const SymMatrixF.row2");
            Assert(s1fcv.row(i+1,j+1,N)(1) == value,"SymMatrixF.row2 CV");
            Assert(s1fv.row(i+1,j+1,N)(1) == value,"SymMatrixF.row2 V");
            Assert(s1f.col(i+1,j+1,N)(1) == value,"SymMatrixF.col2");
            Assert(s1fx.col(i+1,j+1,N)(1) == value,"const SymMatrixF.col2");
            Assert(s1fcv.col(i+1,j+1,N)(1) == value,"SymMatrixF.col2 CV");
            Assert(s1fv.col(i+1,j+1,N)(1) == value,"SymMatrixF.col2 V");
            int d = int(j)-int(i);
            if (d==0) {
                Assert(s1.diag()(i) == value,"SymMatrix.diag");
                Assert(s1x.diag()(i) == value,"const SymMatrix.diag");
                Assert(s1cv.diag()(i) == value,"SymMatrix.diag CV");
                Assert(s1v.diag()(i) == value,"SymMatrix.diag V");
                Assert(s1f.diag()(i+1) == value,"SymMatrixF.diag");
                Assert(s1fx.diag()(i+1) == value,"const SymMatrixF.diag");
                Assert(s1fcv.diag()(i+1) == value,"SymMatrixF.diag CV");
                Assert(s1fv.diag()(i+1) == value,"SymMatrixF.diag V");
            }
            Assert(s1.diag(d)(i) == value,"SymMatrix.diag1");
            Assert(s1x.diag(d)(i) == value,"const SymMatrix.diag1");
            Assert(s1cv.diag(d)(i) == value,"SymMatrix.diag1 CV");
            Assert(s1v.diag(d)(i) == value,"SymMatrix.diag1 V");
            Assert(s1.diag(d,i,N-d)(0) == value,"SymMatrix.diag2");
            Assert(s1x.diag(d,i,N-d)(0) == value,"const SymMatrix.diag2");
            Assert(s1cv.diag(d,i,N-d)(0) == value,"SymMatrix.diag2 CV");
            Assert(s1v.diag(d,i,N-d)(0) == value,"SymMatrix.diag2 V");
            Assert(s1f.diag(d)(i+1) == value,"SymMatrixF.diag1");
            Assert(s1fx.diag(d)(i+1) == value,"const SymMatrixF.diag1");
            Assert(s1fcv.diag(d)(i+1) == value,"SymMatrixF.diag1 CV");
            Assert(s1fv.diag(d)(i+1) == value,"SymMatrixF.diag1 V");
            Assert(s1f.diag(d,i+1,N-d)(1) == value,"SymMatrixF.diag2");
            Assert(s1fx.diag(d,i+1,N-d)(1) == value,"const SymMatrixF.diag2");
            Assert(s1fcv.diag(d,i+1,N-d)(1) == value,"SymMatrixF.diag2 CV");
            Assert(s1fv.diag(d,i+1,N-d)(1) == value,"SymMatrixF.diag2 V");
            Assert(s1.diag(-d)(i) == value,"SymMatrix.diag1");
            Assert(s1x.diag(-d)(i) == value,"const SymMatrix.diag1");
            Assert(s1cv.diag(-d)(i) == value,"SymMatrix.diag1 CV");
            Assert(s1v.diag(-d)(i) == value,"SymMatrix.diag1 V");
            Assert(s1.diag(-d,i,N-d)(0) == value,"SymMatrix.diag2");
            Assert(s1x.diag(-d,i,N-d)(0) == value,"const SymMatrix.diag2");
            Assert(s1cv.diag(-d,i,N-d)(0) == value,"SymMatrix.diag2 CV");
            Assert(s1v.diag(-d,i,N-d)(0) == value,"SymMatrix.diag2 V");
            Assert(s1f.diag(-d)(i+1) == value,"SymMatrixF.diag1");
            Assert(s1fx.diag(-d)(i+1) == value,"const SymMatrixF.diag1");
            Assert(s1fcv.diag(-d)(i+1) == value,"SymMatrixF.diag1 CV");
            Assert(s1fv.diag(-d)(i+1) == value,"SymMatrixF.diag1 V");
            Assert(s1f.diag(-d,i+1,N-d)(1) == value,"SymMatrixF.diag2");
            Assert(s1fx.diag(-d,i+1,N-d)(1) == value,"const SymMatrixF.diag2");
            Assert(s1fcv.diag(-d,i+1,N-d)(1) == value,"SymMatrixF.diag2 CV");
            Assert(s1fv.diag(-d,i+1,N-d)(1) == value,"SymMatrixF.diag2 V");
        }
        if (j<=i) {
            Assert(s2(j,i) == value,"Read/Write SymMatrix");
            Assert(s2x(j,i) == value,"Access const SymMatrix");
            Assert(s2v(j,i) == value,"Access SymMatrix V");
            Assert(s2cv(j,i) == value,"Access SymMatrix CV");
            Assert(s2(i,j) == value,"Access SymMatrix - opposite tri");
            Assert(s2x(i,j) == value,"Access const SymMatrix - opposite tri");
            Assert(s2v(i,j) == value,"Access SymMatrix V");
            Assert(s2cv(i,j) == value,"Access SymMatrix CV");
            Assert(s2f(j+1,i+1) == value,"Read/Write SymMatrixF");
            Assert(s2fx(j+1,i+1) == value,"Access const SymMatrixF");
            Assert(s2fv(j+1,i+1) == value,"Access SymMatrixF V");
            Assert(s2fcv(j+1,i+1) == value,"Access SymMatrixF CV");
            Assert(s2f(i+1,j+1) == value,"Access SymMatrixF - opposite tri");
            Assert(s2fx(i+1,j+1) == value,
                   "Access const SymMatrixF - opposite tri");
            Assert(s2fv(i+1,j+1) == value,"Access SymMatrixF V");
            Assert(s2fcv(i+1,j+1) == value,"Access SymMatrixF CV");
            Assert(s2.row(i,0,i+1)(j) == value,"SymMatrix.row");
            Assert(s2x.row(i,0,i+1)(j) == value,"const SymMatrix.row");
            Assert(s2cv.row(i,0,i+1)(j) == value,"SymMatrix.row CV");
            Assert(s2v.row(i,0,i+1)(j) == value,"SymMatrix.row V");
            Assert(s2.col(i,0,i+1)(j) == value,"SymMatrix.col");
            Assert(s2x.col(i,0,i+1)(j) == value,"const SymMatrix.col");
            Assert(s2cv.col(i,0,i+1)(j) == value,"SymMatrix.col CV");
            Assert(s2v.col(i,0,i+1)(j) == value,"SymMatrix.col V");
            Assert(s2.row(i,j,i+1)(0) == value,"SymMatrix.row2");
            Assert(s2x.row(i,j,i+1)(0) == value,"const SymMatrix.row2");
            Assert(s2cv.row(i,j,i+1)(0) == value,"SymMatrix.row2 CV");
            Assert(s2v.row(i,j,i+1)(0) == value,"SymMatrix.row2 V");
            Assert(s2.col(i,j,i+1)(0) == value,"SymMatrix.col2");
            Assert(s2x.col(i,j,i+1)(0) == value,"const SymMatrix.col2");
            Assert(s2cv.col(i,j,i+1)(0) == value,"SymMatrix.col2 CV");
            Assert(s2v.col(i,j,i+1)(0) == value,"SymMatrix.col2 V");
            Assert(s2f.row(i+1,1,i+1)(j+1) == value,"SymMatrixF.row");
            Assert(s2fx.row(i+1,1,i+1)(j+1) == value,"const SymMatrixF.row");
            Assert(s2fcv.row(i+1,1,i+1)(j+1) == value,"SymMatrixF.row CV");
            Assert(s2fv.row(i+1,1,i+1)(j+1) == value,"SymMatrixF.row V");
            Assert(s2f.col(i+1,1,i+1)(j+1) == value,"SymMatrixF.col");
            Assert(s2fx.col(i+1,1,i+1)(j+1) == value,"const SymMatrixF.col");
            Assert(s2fcv.col(i+1,1,i+1)(j+1) == value,"SymMatrixF.col CV");
            Assert(s2fv.col(i+1,1,i+1)(j+1) == value,"SymMatrixF.col V");
            Assert(s2f.row(i+1,j+1,i+1)(1) == value,"SymMatrixF.row2");
            Assert(s2fx.row(i+1,j+1,i+1)(1) == value,"const SymMatrixF.row2");
            Assert(s2fcv.row(i+1,j+1,i+1)(1) == value,"SymMatrixF.row2 CV");
            Assert(s2fv.row(i+1,j+1,i+1)(1) == value,"SymMatrixF.row2 V");
            Assert(s2f.col(i+1,j+1,i+1)(1) == value,"SymMatrixF.col2");
            Assert(s2fx.col(i+1,j+1,i+1)(1) == value,"const SymMatrixF.col2");
            Assert(s2fcv.col(i+1,j+1,i+1)(1) == value,"SymMatrixF.col2 CV");
            Assert(s2fv.col(i+1,j+1,i+1)(1) == value,"SymMatrixF.col2 V");
            int d = int(j)-int(i);
            if (d==0) {
                Assert(s2.diag()(i) == value,"SymMatrix.diag");
                Assert(s2x.diag()(i) == value,"const SymMatrix.diag");
                Assert(s2cv.diag()(i) == value,"SymMatrix.diag CV");
                Assert(s2v.diag()(i) == value,"SymMatrix.diag V");
                Assert(s2f.diag()(i+1) == value,"SymMatrixF.diag");
                Assert(s2fx.diag()(i+1) == value,"const SymMatrixF.diag");
                Assert(s2fcv.diag()(i+1) == value,"SymMatrixF.diag CV");
                Assert(s2fv.diag()(i+1) == value,"SymMatrixF.diag V");
            }
            Assert(s2.diag(d)(j) == value,"SymMatrix.diag1");
            Assert(s2x.diag(d)(j) == value,"const SymMatrix.diag1");
            Assert(s2cv.diag(d)(j) == value,"SymMatrix.diag1 CV");
            Assert(s2v.diag(d)(j) == value,"SymMatrix.diag1 V");
            Assert(s2.diag(d,j,N+d)(0) == value,"SymMatrix.diag2");
            Assert(s2x.diag(d,j,N+d)(0) == value,"const SymMatrix.diag2");
            Assert(s2cv.diag(d,j,N+d)(0) == value,"SymMatrix.diag2 CV");
            Assert(s2v.diag(d,j,N+d)(0) == value,"SymMatrix.diag2 V");
            Assert(s2f.diag(d)(j+1) == value,"SymMatrixF.diag1");
            Assert(s2fx.diag(d)(j+1) == value,"const SymMatrixF.diag1");
            Assert(s2fcv.diag(d)(j+1) == value,"SymMatrixF.diag1 CV");
            Assert(s2fv.diag(d)(j+1) == value,"SymMatrixF.diag1 V");
            Assert(s2f.diag(d,j+1,N+d)(1) == value,"SymMatrixF.diag2");
            Assert(s2fx.diag(d,j+1,N+d)(1) == value,"const SymMatrixF.diag2");
            Assert(s2fcv.diag(d,j+1,N+d)(1) == value,"SymMatrixF.diag2 CV");
            Assert(s2fv.diag(d,j+1,N+d)(1) == value,"SymMatrixF.diag2 V");
            Assert(s2.diag(-d)(j) == value,"SymMatrix.diag1");
            Assert(s2x.diag(-d)(j) == value,"const SymMatrix.diag1");
            Assert(s2cv.diag(-d)(j) == value,"SymMatrix.diag1 CV");
            Assert(s2v.diag(-d)(j) == value,"SymMatrix.diag1 V");
            Assert(s2.diag(-d,j,N+d)(0) == value,"SymMatrix.diag2");
            Assert(s2x.diag(-d,j,N+d)(0) == value,"const SymMatrix.diag2");
            Assert(s2cv.diag(-d,j,N+d)(0) == value,"SymMatrix.diag2 CV");
            Assert(s2v.diag(-d,j,N+d)(0) == value,"SymMatrix.diag2 V");
            Assert(s2f.diag(-d)(j+1) == value,"SymMatrixF.diag1");
            Assert(s2fx.diag(-d)(j+1) == value,"const SymMatrixF.diag1");
            Assert(s2fcv.diag(-d)(j+1) == value,"SymMatrixF.diag1 CV");
            Assert(s2fv.diag(-d)(j+1) == value,"SymMatrixF.diag1 V");
            Assert(s2f.diag(-d,j+1,N+d)(1) == value,"SymMatrixF.diag2");
            Assert(s2fx.diag(-d,j+1,N+d)(1) == value,"const SymMatrixF.diag2");
            Assert(s2fcv.diag(-d,j+1,N+d)(1) == value,"SymMatrixF.diag2 CV");
            Assert(s2fv.diag(-d,j+1,N+d)(1) == value,"SymMatrixF.diag2 V");
        }
    }

    Assert(s1 == s1f,"CStyle SymMatrix == FortranStyle SymMatrix");
    Assert(s1 == s1cv,"SymMatrix == ConstSymMatrixView");
    Assert(s1 == s1v,"SymMatrix == SymMatrixView");
    Assert(s1 == s1fcv,"SymMatrix == FortranStyle ConstSymMatrixView");
    Assert(s1 == s1fv,"SymMatrix == FortranStyle SymMatrixView");
    Assert(s2 == s2f,"CStyle SymMatrix == FortranStyle SymMatrix");
    Assert(s2 == s2cv,"SymMatrix == ConstSymMatrixView");
    Assert(s2 == s2v,"SymMatrix == SymMatrixView");
    Assert(s2 == s2fcv,"SymMatrix == FortranStyle ConstSymMatrixView");
    Assert(s2 == s2fv,"SymMatrix == FortranStyle SymMatrixView");

    s1.resize(3);
    Assert(s1.colsize() == 3 && s1.rowsize() == 3,
           "SymMatrix s1.resize(3)");
    s2.resize(3);
    Assert(s2.colsize() == 3 && s2.rowsize() == 3,
           "SymMatrix s2.resize(3)");
    for (int i=0, k=0; i<3; ++i) for (int j=0; j<3; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (i <= j) s1(i,j) = hvalue; 
        if (j <= i) s2(i,j) = hvalue; 
    }
    for (int i=0, k=0; i<3; ++i) for (int j=0; j<3; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (i<=j) {
            Assert(s1(i,j) == hvalue,"Read/Write resized SymMatrix");
            Assert(s1(j,i) == hvalue,
                   "Read/Write resized SymMatrix opp tri");
        }
        if (j<=i) {
            Assert(s2(j,i) == hvalue,"Read/Write resized SymMatrix");
            Assert(s2(i,j) == hvalue,
                   "Read/Write resized SymMatrix opp tri");
        }
    }

    s1.resize(2*N);
    Assert(int(s1.colsize()) == 2*N && int(s1.rowsize()) == 2*N,
           "SymMatrix s1.resize(2*N) sizes");
    s2.resize(2*N);
    Assert(int(s2.colsize()) == 2*N && int(s2.rowsize()) == 2*N,
           "SymMatrix s2.resize(2*N) sizes");
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<2*N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        if (i <= j) s1(i,j) = value; 
        if (j <= i) s2(i,j) = value; 
    }
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<2*N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        if (i<=j) {
            Assert(s1(i,j) == value,"Read/Write resized SymMatrix");
            Assert(s1(j,i) == value,
                   "Read/Write resized SymMatrix opp tri");
        }
        if (j<=i) {
            Assert(s2(j,i) == value,"Read/Write resized SymMatrix");
            Assert(s2(i,j) == value,
                   "Read/Write resized SymMatrix opp tri");
        }
    }
}

template <class T, tmv::UpLoType U, tmv::StorageType S> 
static void TestBasicHermMatrix_1()
{
    const int N = 6;

    if (showstartdone) {
        std::cout<<"Start TestBasicHermMatrix_1\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"U = "<<tmv::TMV_Text(U)<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
    }

    tmv::HermMatrix<std::complex<T>,U|S> h1(N);
    tmv::HermMatrix<std::complex<T>,U|S> h2(N);
    tmv::HermMatrix<std::complex<T>,U|S|tmv::FortranStyle> h1f(N);
    tmv::HermMatrix<std::complex<T>,U|S|tmv::FortranStyle> h2f(N);

    Assert(int(h1.colsize()) == N && int(h1.rowsize()) == N,
           "Creating HermMatrix(N)");

    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (i<=j) {
            h1(i,j) = hvalue; 
            h1f(i+1,j+1) = hvalue;
        }
        if (j<=i) {
            h2(i,j) = hvalue; 
            h2f(i+1,j+1) = hvalue; 
        }
    }

    tmv::SymMatrixView<std::complex<T> > h1v = h1.view();
    tmv::ConstSymMatrixView<std::complex<T> > h1cv = h1.view();
    tmv::SymMatrixView<std::complex<T> > h2v = h2.view();
    tmv::ConstSymMatrixView<std::complex<T> > h2cv = h2.view();
    tmv::SymMatrixView<std::complex<T>,tmv::FortranStyle> h1fv = h1f.view();
    tmv::ConstSymMatrixView<std::complex<T>,tmv::FortranStyle> h1fcv =
        h1f.view();
    tmv::SymMatrixView<std::complex<T>,tmv::FortranStyle> h2fv = h2f.view();
    tmv::ConstSymMatrixView<std::complex<T>,tmv::FortranStyle> h2fcv =
        h2f.view();
    const tmv::HermMatrix<std::complex<T>,U|S>& h1x = h1;
    const tmv::HermMatrix<std::complex<T>,U|S>& h2x = h2;
    const tmv::HermMatrix<std::complex<T>,U|S|tmv::FortranStyle>& h1fx = h1f;
    const tmv::HermMatrix<std::complex<T>,U|S|tmv::FortranStyle>& h2fx = h2f;

    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (showtests) std::cout<<"i,j = "<<i<<','<<j<<std::endl;
        if (i<=j) {
            Assert(h1(i,j) == hvalue,"Read/Write HermMatrix");
            Assert(h1x(i,j) == hvalue,"Access const HermMatrix");
            Assert(h1v(i,j) == hvalue,"Access HermMatrix V");
            Assert(h1cv(i,j) == hvalue,"Access HermMatrix CV");
            Assert(h1(j,i) == conj(hvalue),"Access HermMatrix - opposite tri");
            Assert(h1x(j,i) == conj(hvalue),
                   "Access const HermMatrix - opposite tri");
            Assert(h1v(j,i) == conj(hvalue),"Access HermMatrix V");
            Assert(h1cv(j,i) == conj(hvalue),"Access HermMatrix CV");
            Assert(h1f(i+1,j+1) == hvalue,"Read/Write HermMatrixF");
            Assert(h1fx(i+1,j+1) == hvalue,"Access const HermMatrixF");
            Assert(h1fv(i+1,j+1) == hvalue,"Access HermMatrixF V");
            Assert(h1fcv(i+1,j+1) == hvalue,"Access HermMatrixF CV");
            Assert(h1f(j+1,i+1) == conj(hvalue),
                   "Access HermMatrixF - opposite tri");
            Assert(h1fx(j+1,i+1) == conj(hvalue),
                   "Access const HermMatrixF - opposite tri");
            Assert(h1fv(j+1,i+1) == conj(hvalue),
                   "Access HermMatrixF V");
            Assert(h1fcv(j+1,i+1) == conj(hvalue),"Access HermMatrixF CV");
            Assert(h1.row(i,i,N)(j-i) == hvalue,"HermMatrix.row");
            Assert(h1x.row(i,i,N)(j-i) == hvalue,"HermMatrix.row");
            Assert(h1cv.row(i,i,N)(j-i) == hvalue,"HermMatrix.row CV");
            Assert(h1v.row(i,i,N)(j-i) == hvalue,"HermMatrix.row V");
            Assert(h1.col(i,i,N)(j-i) == conj(hvalue),"HermMatrix.col");
            Assert(h1x.col(i,i,N)(j-i) == conj(hvalue),"HermMatrix.col");
            Assert(h1cv.col(i,i,N)(j-i) == conj(hvalue),"HermMatrix.col CV");
            Assert(h1v.col(i,i,N)(j-i) == conj(hvalue),"HermMatrix.col V");
            Assert(h1.row(i,j,N)(0) == hvalue,"HermMatrix.row2");
            Assert(h1x.row(i,j,N)(0) == hvalue,"HermMatrix.row2");
            Assert(h1cv.row(i,j,N)(0) == hvalue,"HermMatrix.row2 CV");
            Assert(h1v.row(i,j,N)(0) == hvalue,"HermMatrix.row2 V");
            Assert(h1.col(i,j,N)(0) == conj(hvalue),"HermMatrix.col2");
            Assert(h1x.col(i,j,N)(0) == conj(hvalue),"HermMatrix.col2");
            Assert(h1cv.col(i,j,N)(0) == conj(hvalue),"HermMatrix.col2 CV");
            Assert(h1v.col(i,j,N)(0) == conj(hvalue),"HermMatrix.col2 V");
            Assert(h1f.row(i+1,i+1,N)(j-i+1) == hvalue,"HermMatrixF.row");
            Assert(h1fx.row(i+1,i+1,N)(j-i+1) == hvalue,
                   "const HermMatrixF.row");
            Assert(h1fcv.row(i+1,i+1,N)(j-i+1) == hvalue,"HermMatrixF.row CV");
            Assert(h1fv.row(i+1,i+1,N)(j-i+1) == hvalue,"HermMatrixF.row V");
            Assert(h1f.col(i+1,i+1,N)(j-i+1) == conj(hvalue),"HermMatrixF.col");
            Assert(h1fx.col(i+1,i+1,N)(j-i+1) == conj(hvalue),
                   "const HermMatrixF.col");
            Assert(h1fcv.col(i+1,i+1,N)(j-i+1) == conj(hvalue),
                   "HermMatrixF.col CV");
            Assert(h1fv.col(i+1,i+1,N)(j-i+1) == conj(hvalue),
                   "HermMatrixF.col V");
            Assert(h1f.row(i+1,j+1,N)(1) == hvalue,"HermMatrixF.row2");
            Assert(h1fx.row(i+1,j+1,N)(1) == hvalue,"const HermMatrixF.row2");
            Assert(h1fcv.row(i+1,j+1,N)(1) == hvalue,"HermMatrixF.row2 CV");
            Assert(h1fv.row(i+1,j+1,N)(1) == hvalue,"HermMatrixF.row2 V");
            Assert(h1f.col(i+1,j+1,N)(1) == conj(hvalue),"HermMatrixF.col2");
            Assert(h1fx.col(i+1,j+1,N)(1) == conj(hvalue),
                   "const HermMatrixF.col2");
            Assert(h1fcv.col(i+1,j+1,N)(1) == conj(hvalue),
                   "HermMatrixF.col2 CV");
            Assert(h1fv.col(i+1,j+1,N)(1) == conj(hvalue),"HermMatrixF.col2 V");
            int d = int(j)-int(i);
            if (d==0) {
                Assert(h1.diag()(i) == hvalue,"HermMatrix.diag");
                Assert(h1x.diag()(i) == hvalue,"const HermMatrix.diag");
                Assert(h1cv.diag()(i) == hvalue,"HermMatrix.diag CV");
                Assert(h1v.diag()(i) == hvalue,"HermMatrix.diag V");
                Assert(h1f.diag()(i+1) == hvalue,"HermMatrixF.diag");
                Assert(h1fx.diag()(i+1) == hvalue,"const HermMatrixF.diag");
                Assert(h1fcv.diag()(i+1) == hvalue,"HermMatrixF.diag CV");
                Assert(h1fv.diag()(i+1) == hvalue,"HermMatrixF.diag V");
            }
            Assert(h1.diag(d)(i) == hvalue,"HermMatrix.diag1");
            Assert(h1x.diag(d)(i) == hvalue,"const HermMatrix.diag1");
            Assert(h1cv.diag(d)(i) == hvalue,"HermMatrix.diag1 CV");
            Assert(h1v.diag(d)(i) == hvalue,"HermMatrix.diag1 V");
            Assert(h1.diag(d,i,N-d)(0) == hvalue,"HermMatrix.diag2");
            Assert(h1x.diag(d,i,N-d)(0) == hvalue,"const HermMatrix.diag2");
            Assert(h1cv.diag(d,i,N-d)(0) == hvalue,"HermMatrix.diag2 CV");
            Assert(h1v.diag(d,i,N-d)(0) == hvalue,"HermMatrix.diag2 V");
            Assert(h1f.diag(d)(i+1) == hvalue,"HermMatrixF.diag1");
            Assert(h1fx.diag(d)(i+1) == hvalue,"const HermMatrixF.diag1");
            Assert(h1fcv.diag(d)(i+1) == hvalue,"HermMatrixF.diag1 CV");
            Assert(h1fv.diag(d)(i+1) == hvalue,"HermMatrixF.diag1 V");
            Assert(h1f.diag(d,i+1,N-d)(1) == hvalue,"HermMatrixF.diag2");
            Assert(h1fx.diag(d,i+1,N-d)(1) == hvalue,"const HermMatrixF.diag2");
            Assert(h1fcv.diag(d,i+1,N-d)(1) == hvalue,"HermMatrixF.diag2 CV");
            Assert(h1fv.diag(d,i+1,N-d)(1) == hvalue,"HermMatrixF.diag2 V");
            Assert(h1.diag(-d)(i) == conj(hvalue),"HermMatrix.diag1");
            Assert(h1x.diag(-d)(i) == conj(hvalue),"const HermMatrix.diag1");
            Assert(h1cv.diag(-d)(i) == conj(hvalue),"HermMatrix.diag1 CV");
            Assert(h1v.diag(-d)(i) == conj(hvalue),"HermMatrix.diag1 V");
            Assert(h1.diag(-d,i,N-d)(0) == conj(hvalue),"HermMatrix.diag2");
            Assert(h1x.diag(-d,i,N-d)(0) == conj(hvalue),
                   "const HermMatrix.diag2");
            Assert(h1cv.diag(-d,i,N-d)(0) == conj(hvalue),
                   "HermMatrix.diag2 CV");
            Assert(h1v.diag(-d,i,N-d)(0) == conj(hvalue),"HermMatrix.diag2 V");
            Assert(h1f.diag(-d)(i+1) == conj(hvalue),"HermMatrixF.diag1");
            Assert(h1fx.diag(-d)(i+1) == conj(hvalue),
                   "const HermMatrixF.diag1");
            Assert(h1fcv.diag(-d)(i+1) == conj(hvalue),"HermMatrixF.diag1 CV");
            Assert(h1fv.diag(-d)(i+1) == conj(hvalue),"HermMatrixF.diag1 V");
            Assert(h1f.diag(-d,i+1,N-d)(1) == conj(hvalue),"HermMatrixF.diag2");
            Assert(h1fx.diag(-d,i+1,N-d)(1) == conj(hvalue),
                   "const HermMatrixF.diag2");
            Assert(h1fcv.diag(-d,i+1,N-d)(1) == conj(hvalue),
                   "HermMatrixF.diag2 CV");
            Assert(h1fv.diag(-d,i+1,N-d)(1) == conj(hvalue),
                   "HermMatrixF.diag2 V");
        }
        if (j<=i) {
            Assert(h2(j,i) == conj(hvalue),"Read/Write HermMatrix");
            Assert(h2x(j,i) == conj(hvalue),"Access const HermMatrix");
            Assert(h2v(j,i) == conj(hvalue),"Access HermMatrix V");
            Assert(h2cv(j,i) == conj(hvalue),"Access HermMatrix CV");
            Assert(h2(i,j) == hvalue,"Access HermMatrix - opposite tri");
            Assert(h2x(i,j) == hvalue,"Access const HermMatrix - opposite tri");
            Assert(h2v(i,j) == hvalue,"Access HermMatrix V");
            Assert(h2cv(i,j) == hvalue,"Access HermMatrix CV");
            Assert(h2f(j+1,i+1) == conj(hvalue),"Read/Write HermMatrixF");
            Assert(h2fx(j+1,i+1) == conj(hvalue),"Access const HermMatrixF");
            Assert(h2fv(j+1,i+1) == conj(hvalue),"Access HermMatrixF V");
            Assert(h2fcv(j+1,i+1) == conj(hvalue),"Access HermMatrixF CV");
            Assert(h2f(i+1,j+1) == hvalue,"Access HermMatrixF - opposite tri");
            Assert(h2fx(i+1,j+1) == hvalue,
                   "Access const HermMatrixF - opposite tri");
            Assert(h2fv(i+1,j+1) == hvalue,"Access HermMatrixF V");
            Assert(h2fcv(i+1,j+1) == hvalue,"Access HermMatrixF CV");
            Assert(h2.row(i,0,i+1)(j) == hvalue,"HermMatrix.row");
            Assert(h2x.row(i,0,i+1)(j) == hvalue,"const HermMatrix.row");
            Assert(h2cv.row(i,0,i+1)(j) == hvalue,"HermMatrix.row CV");
            Assert(h2v.row(i,0,i+1)(j) == hvalue,"HermMatrix.row V");
            Assert(h2.col(i,0,i+1)(j) == conj(hvalue),"HermMatrix.col");
            Assert(h2x.col(i,0,i+1)(j) == conj(hvalue),"const HermMatrix.col");
            Assert(h2cv.col(i,0,i+1)(j) == conj(hvalue),"HermMatrix.col CV");
            Assert(h2v.col(i,0,i+1)(j) == conj(hvalue),"HermMatrix.col V");
            Assert(h2.row(i,j,i+1)(0) == hvalue,"HermMatrix.row2");
            Assert(h2x.row(i,j,i+1)(0) == hvalue,"const HermMatrix.row2");
            Assert(h2cv.row(i,j,i+1)(0) == hvalue,"HermMatrix.row2 CV");
            Assert(h2v.row(i,j,i+1)(0) == hvalue,"HermMatrix.row2 V");
            Assert(h2.col(i,j,i+1)(0) == conj(hvalue),"HermMatrix.col2");
            Assert(h2x.col(i,j,i+1)(0) == conj(hvalue),"const HermMatrix.col2");
            Assert(h2cv.col(i,j,i+1)(0) == conj(hvalue),"HermMatrix.col2 CV");
            Assert(h2v.col(i,j,i+1)(0) == conj(hvalue),"HermMatrix.col2 V");
            Assert(h2f.row(i+1,1,i+1)(j+1) == hvalue,"HermMatrixF.row");
            Assert(h2fx.row(i+1,1,i+1)(j+1) == hvalue,"const HermMatrixF.row");
            Assert(h2fcv.row(i+1,1,i+1)(j+1) == hvalue,"HermMatrixF.row CV");
            Assert(h2fv.row(i+1,1,i+1)(j+1) == hvalue,"HermMatrixF.row V");
            Assert(h2f.col(i+1,1,i+1)(j+1) == conj(hvalue),"HermMatrixF.col");
            Assert(h2fx.col(i+1,1,i+1)(j+1) == conj(hvalue),
                   "const HermMatrixF.col");
            Assert(h2fcv.col(i+1,1,i+1)(j+1) == conj(hvalue),
                   "HermMatrixF.col CV");
            Assert(h2fv.col(i+1,1,i+1)(j+1) == conj(hvalue),
                   "HermMatrixF.col V");
            Assert(h2f.row(i+1,j+1,i+1)(1) == hvalue,"HermMatrixF.row2");
            Assert(h2fx.row(i+1,j+1,i+1)(1) == hvalue,"const HermMatrixF.row2");
            Assert(h2fcv.row(i+1,j+1,i+1)(1) == hvalue,"HermMatrixF.row2 CV");
            Assert(h2fv.row(i+1,j+1,i+1)(1) == hvalue,"HermMatrixF.row2 V");
            Assert(h2f.col(i+1,j+1,i+1)(1) == conj(hvalue),"HermMatrixF.col2");
            Assert(h2fx.col(i+1,j+1,i+1)(1) == conj(hvalue),
                   "const HermMatrixF.col2");
            Assert(h2fcv.col(i+1,j+1,i+1)(1) == conj(hvalue),
                   "HermMatrixF.col2 CV");
            Assert(h2fv.col(i+1,j+1,i+1)(1) == conj(hvalue),
                   "HermMatrixF.col2 V");
            int d = int(j)-int(i);
            if (d==0) {
                Assert(h2.diag()(i) == hvalue,"HermMatrix.diag");
                Assert(h2x.diag()(i) == hvalue,"const HermMatrix.diag");
                Assert(h2cv.diag()(i) == hvalue,"HermMatrix.diag CV");
                Assert(h2v.diag()(i) == hvalue,"HermMatrix.diag V");
                Assert(h2f.diag()(i+1) == hvalue,"HermMatrixF.diag");
                Assert(h2fx.diag()(i+1) == hvalue,"const HermMatrixF.diag");
                Assert(h2fcv.diag()(i+1) == hvalue,"HermMatrixF.diag CV");
                Assert(h2fv.diag()(i+1) == hvalue,"HermMatrixF.diag V");
            }
            Assert(h2.diag(d)(j) == hvalue,"HermMatrix.diag1");
            Assert(h2x.diag(d)(j) == hvalue,"const HermMatrix.diag1");
            Assert(h2cv.diag(d)(j) == hvalue,"HermMatrix.diag1 CV");
            Assert(h2v.diag(d)(j) == hvalue,"HermMatrix.diag1 V");
            Assert(h2.diag(d,j,N+d)(0) == hvalue,"HermMatrix.diag2");
            Assert(h2x.diag(d,j,N+d)(0) == hvalue,"const HermMatrix.diag2");
            Assert(h2cv.diag(d,j,N+d)(0) == hvalue,"HermMatrix.diag2 CV");
            Assert(h2v.diag(d,j,N+d)(0) == hvalue,"HermMatrix.diag2 V");
            Assert(h2f.diag(d)(j+1) == hvalue,"HermMatrixF.diag1");
            Assert(h2fx.diag(d)(j+1) == hvalue,"const HermMatrixF.diag1");
            Assert(h2fcv.diag(d)(j+1) == hvalue,"HermMatrixF.diag1 CV");
            Assert(h2fv.diag(d)(j+1) == hvalue,"HermMatrixF.diag1 V");
            Assert(h2f.diag(d,j+1,N+d)(1) == hvalue,"HermMatrixF.diag2");
            Assert(h2fx.diag(d,j+1,N+d)(1) == hvalue,"const HermMatrixF.diag2");
            Assert(h2fcv.diag(d,j+1,N+d)(1) == hvalue,"HermMatrixF.diag2 CV");
            Assert(h2fv.diag(d,j+1,N+d)(1) == hvalue,"HermMatrixF.diag2 V");
            Assert(h2.diag(-d)(j) == conj(hvalue),"HermMatrix.diag1");
            Assert(h2x.diag(-d)(j) == conj(hvalue),"const HermMatrix.diag1");
            Assert(h2cv.diag(-d)(j) == conj(hvalue),"HermMatrix.diag1 CV");
            Assert(h2v.diag(-d)(j) == conj(hvalue),"HermMatrix.diag1 V");
            Assert(h2.diag(-d,j,N+d)(0) == conj(hvalue),"HermMatrix.diag2");
            Assert(h2x.diag(-d,j,N+d)(0) == conj(hvalue),
                   "const HermMatrix.diag2");
            Assert(h2cv.diag(-d,j,N+d)(0) == conj(hvalue),
                   "HermMatrix.diag2 CV");
            Assert(h2v.diag(-d,j,N+d)(0) == conj(hvalue),"HermMatrix.diag2 V");
            Assert(h2f.diag(-d)(j+1) == conj(hvalue),"HermMatrixF.diag1");
            Assert(h2fx.diag(-d)(j+1) == conj(hvalue),
                   "const HermMatrixF.diag1");
            Assert(h2fcv.diag(-d)(j+1) == conj(hvalue),"HermMatrixF.diag1 CV");
            Assert(h2fv.diag(-d)(j+1) == conj(hvalue),"HermMatrixF.diag1 V");
            Assert(h2f.diag(-d,j+1,N+d)(1) == conj(hvalue),"HermMatrixF.diag2");
            Assert(h2fx.diag(-d,j+1,N+d)(1) == conj(hvalue),
                   "const HermMatrixF.diag2");
            Assert(h2fcv.diag(-d,j+1,N+d)(1) == conj(hvalue),
                   "HermMatrixF.diag2 CV");
            Assert(h2fv.diag(-d,j+1,N+d)(1) == conj(hvalue),
                   "HermMatrixF.diag2 V");
        }
    }

    Assert(h1 == h1f,"CStyle HermMatrix == FortranStyle HermMatrix");
    Assert(h1 == h1cv,"HermMatrix == ConstSymMatrixView");
    Assert(h1 == h1v,"HermMatrix == SymMatrixView");
    Assert(h1 == h1fcv,"HermMatrix == FortranStyle ConstSymMatrixView");
    Assert(h1 == h1fv,"HermMatrix == FortranStyle SymMatrixView");
    Assert(h2 == h2f,"CStyle HermMatrix == FortranStyle HermMatrix");
    Assert(h2 == h2cv,"HermMatrix == ConstSymMatrixView");
    Assert(h2 == h2v,"HermMatrix == SymMatrixView");
    Assert(h2 == h2fcv,"HermMatrix == FortranStyle ConstSymMatrixView");
    Assert(h2 == h2fv,"HermMatrix == FortranStyle SymMatrixView");

    h1.resize(3);
    Assert(h1.colsize() == 3 && h1.rowsize() == 3,
           "HermBandMatrix h1.resize(3)");
    h2.resize(3);
    Assert(h2.colsize() == 3 && h2.rowsize() == 3,
           "HermBandMatrix h2.resize(3)");
    for (int i=0, k=0; i<3; ++i) for (int j=0; j<3; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (i <= j) h1(i,j) = hvalue; 
        if (j <= i) h2(i,j) = hvalue; 
    }
    for (int i=0, k=0; i<3; ++i) for (int j=0; j<3; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (i<=j) {
            Assert(h1(i,j) == hvalue,"Read/Write resized HermMatrix");
            Assert(h1(j,i) == std::conj(hvalue),
                   "Read/Write resized HermMatrix opp tri");
        }
        if (j<=i) {
            Assert(h2(i,j) == hvalue,"Read/Write resized HermMatrix");
            Assert(h2(j,i) == std::conj(hvalue),
                   "Read/Write resized HermMatrix opp tri");
        }
    }

    h1.resize(2*N);
    Assert(int(h1.colsize()) == 2*N && int(h1.rowsize()) == 2*N,
           "HermMatrix h1.resize(2*N) sizes");
    h2.resize(2*N);
    Assert(int(h2.colsize()) == 2*N && int(h2.rowsize()) == 2*N,
           "HermMatrix h2.resize(2*N) sizes");
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<2*N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (i <= j) h1(i,j) = hvalue; 
        if (j <= i) h2(i,j) = hvalue; 
    }
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<2*N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (i<=j) {
            Assert(h1(i,j) == hvalue,"Read/Write resized HermMatrix");
            Assert(h1(j,i) == std::conj(hvalue),
                   "Read/Write resized HermMatrix opp tri");
        }
        if (j<=i) {
            Assert(h2(i,j) == hvalue,"Read/Write resized HermMatrix");
            Assert(h2(j,i) == std::conj(hvalue),
                   "Read/Write resized HermMatrix opp tri");
        }
    }
}

template <class T, tmv::UpLoType U, tmv::StorageType S> 
static void TestBasicSymMatrix_2()
{
    const int N = 6;

    if (showstartdone) {
        std::cout<<"Start TestBasicSymMatrix_2\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"U = "<<tmv::TMV_Text(U)<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
    }

    // Test assignments and constructors from arrays
    const T quarrm[] = {
        T(0), T(2), T(4),
              T(1), T(3),
                    T(2)
    };
    const T qlarrm[] = {
        T(0),
        T(2), T(1),
        T(4), T(3), T(2)
    };
    const T quarcm[] = {
        T(0),
        T(2), T(1),
        T(4), T(3), T(2)
    };
    const T qlarcm[] = {
        T(0), T(2), T(4),
              T(1), T(3),
                    T(2)
    };

    std::vector<T> quvecrm(6);
    for(int i=0;i<6;i++) quvecrm[i] = quarrm[i];
    std::vector<T> qlvecrm(6);
    for(int i=0;i<6;i++) qlvecrm[i] = qlarrm[i];
    std::vector<T> quveccm(6);
    for(int i=0;i<6;i++) quveccm[i] = quarcm[i];
    std::vector<T> qlveccm(6);
    for(int i=0;i<6;i++) qlveccm[i] = qlarcm[i];

    tmv::SymMatrix<T,U|S> q1(3);
    std::copy(quarrm, quarrm+6, q1.upperTri().rowmajor_begin());
    tmv::SymMatrix<T,U|S> q2(3);
    std::copy(quarcm, quarcm+6, q2.upperTri().colmajor_begin());

    tmv::SymMatrix<T,U|S> q3(3);
    std::copy(qlarrm, qlarrm+6, q3.lowerTri().rowmajor_begin());
    tmv::SymMatrix<T,U|S> q4(3);
    std::copy(qlarcm, qlarcm+6, q4.lowerTri().colmajor_begin());

    tmv::SymMatrix<T,U|S> q5(3);
    std::copy(quvecrm.begin(), quvecrm.end(), q5.upperTri().rowmajor_begin());
    tmv::SymMatrix<T,U|S> q6(3);
    std::copy(quveccm.begin(), quveccm.end(), q6.upperTri().colmajor_begin());

    tmv::SymMatrix<T,U|S> q7(3);
    std::copy(qlvecrm.begin(), qlvecrm.end(), q7.lowerTri().rowmajor_begin());
    tmv::SymMatrix<T,U|S> q8(3);
    std::copy(qlveccm.begin(), qlveccm.end(), q8.lowerTri().colmajor_begin());

    tmv::SymMatrix<T,U|S> q9x(30);
    tmv::SymMatrixView<T> q9 = q9x.subSymMatrix(3,18,5);
    std::copy(quvecrm.begin(), quvecrm.end(), q9.upperTri().rowmajor_begin());
    tmv::SymMatrix<T,U|S> q10x(30);
    tmv::SymMatrixView<T> q10 = q10x.subSymMatrix(3,18,5);
    std::copy(quveccm.begin(), quveccm.end(), q10.upperTri().colmajor_begin());

    tmv::SymMatrix<T,U|S> q11x(30);
    tmv::SymMatrixView<T> q11 = q11x.subSymMatrix(3,18,5);
    std::copy(qlvecrm.begin(), qlvecrm.end(), q11.lowerTri().rowmajor_begin());
    tmv::SymMatrix<T,U|S> q12x(30);
    tmv::SymMatrixView<T> q12 = q12x.subSymMatrix(3,18,5);
    std::copy(qlveccm.begin(), qlveccm.end(), q12.lowerTri().colmajor_begin());

    // Assignment using op<< is always in rowmajor order.
    tmv::SymMatrix<T,U|S> q13(3);
    tmv::SymMatrix<T,U|S> q14t(3);
    tmv::SymMatrixView<T> q14 = q14t.transpose();

    tmv::SymMatrix<T,U|S> q15(3);
    tmv::SymMatrix<T,U|S> q16t(3);
    tmv::SymMatrixView<T> q16 = q16t.transpose();

    q13.upperTri() <<
        0, 2, 4,
           1, 3,
              2;
    q14.upperTri() <<
        0, 2, 4,
           1, 3,
              2;
    q15.lowerTri() <<
        0,
        2,  1,
        4,  3,  2;
    q16.lowerTri() <<
        0,
        2,  1,
        4,  3,  2;

    // Can also view memory directly
    T quarrmfull[] = {
        T(0), T(2),  T(4),
        T(-1), T(1), T(3),
        T(-2), T(0), T(2)
    };
    T quarcmfull[] = {
        T(0), T(-1), T(-2),
        T(2), T(1), T(0),
        T(4), T(3), T(2)
    };
    T qlarrmfull[] = {
        T(0), T(-1), T(-2),
        T(2), T(1), T(0),
        T(4), T(3), T(2)
    };
    T qlarcmfull[] = {
        T(0), T(2),  T(4),
        T(-1), T(1), T(3),
        T(-2), T(0), T(2)
    };
    T* qarfull = (
        (S == tmv::RowMajor) ? 
        ( (U == tmv::Upper) ? quarrmfull : qlarrmfull ) : 
        ( (U == tmv::Upper) ? quarcmfull : qlarcmfull ) );
    const int Si = (S == tmv::RowMajor) ? 3 : 1;
    const int Sj = (S == tmv::RowMajor) ? 1 : 3;
    const tmv::ConstSymMatrixView<T> q17 =
        tmv::SymMatrixViewOf(qarfull,3,U,S);
    const tmv::ConstSymMatrixView<T> q18 =
        tmv::SymMatrixViewOf(qarfull,3,U,Si,Sj);

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
    }

    for(int i=0;i<3;i++) for(int j=0;j<3;j++) {
        T val = i >= j ? T(2*i-j) : T(2*j-i);
        Assert(q1(i,j) == val,"Create SymMatrix from T* rm upperTri");
        Assert(q2(i,j) == val,"Create SymMatrix from T* cm upperTri");
        Assert(q3(i,j) == val,"Create SymMatrix from T* rm lowerTri");
        Assert(q4(i,j) == val,"Create SymMatrix from T* cm lowerTri");
        Assert(q5(i,j) == val,"Create SymMatrix from vector rm upperTri");
        Assert(q6(i,j) == val,"Create SymMatrix from vector cm upperTri");
        Assert(q7(i,j) == val,"Create SymMatrix from vector rm lowerTri");
        Assert(q8(i,j) == val,"Create SymMatrix from vector cm lowerTri");
        Assert(q9(i,j) == val,"Create SymMatrixView from vector rm upperTri");
        Assert(q10(i,j) == val,"Create SymMatrixView from vector cm upperTri");
        Assert(q11(i,j) == val,"Create SymMatrixView from vector rm lowerTri");
        Assert(q12(i,j) == val,"Create SymMatrixView from vector cm lowerTri");
        Assert(q13(i,j) == val,"Create SymMatrix from << list upperTri");
        Assert(q14(i,j) == val,"Create SymMatrixView from << list upperTri");
        Assert(q15(i,j) == val,"Create SymMatrix from << list lowerTri");
        Assert(q16(i,j) == val,"Create SymMatrixView from << list lowerTri");
        Assert(q17(i,j) == val,"Create SymMatrixView of T* (S)");
        Assert(q18(i,j) == val,"Create SymMatrixView of T* (Si,Sj)");
    }

    // Test Basic Arithmetic
    tmv::SymMatrix<std::complex<T>,U|S> s1(N);
    tmv::SymMatrix<std::complex<T>,U|S> s2(N);

    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        if (i<=j) {
            s1(i,j) = value; 
        }
        if (j<=i) {
            s2(i,j) = value; 
        }
    }

    tmv::SymMatrix<std::complex<T>,U|S> s3(N);
    s3 = s1+s2;

    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        Assert(s3(i,j) == s1(i,j)+s2(i,j),"Add SymMatrices1");
        Assert(s3(j,i) == s1(i,j)+s2(i,j),"Add SymMatrices2");
    }

    s3 = s1-s2;

    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        if (i<=j) {
            Assert(s3(i,j) == s1(i,j)-s2(i,j),"Subtract SymMatrices1");
            Assert(s3(j,i) == s1(i,j)-s2(i,j),"Subtract SymMatrices2");
        }

    tmv::Matrix<std::complex<T> > m1 = s1;
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        Assert(s1(i,j) == m1(i,j),"SymMatrix -> Matrix");
    }
    tmv::SymMatrix<std::complex<T>,U|S> sm1(m1);
    Assert(s1 == sm1,"Matrix -> SymMatrix1");
    tmv::SymMatrix<std::complex<T>,U|S> sm2 = m1;
    Assert(s1 == sm2,"Matrix -> SymMatrix2");

    tmv::SymMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor> sur = s1;
    Assert(s1==sur,"SymMatrix == SymMatrix<U,R>");
    Assert(s1.view()==sur.view(),"SymMatrix.view == SymMatrix<U,R>.view");
    Assert(s1.transpose()==sur.transpose(),
           "SymMatrix.transpose == SymMatrix<U,R>.transpose");
    Assert(s1.conjugate()==sur.conjugate(),
           "SymMatrix.conjugate == SymMatrix<U,R>.conjugate");
    Assert(s1.adjoint()==sur.adjoint(),
           "SymMatrix.adjoint == SymMatrix<U,R>.adjoint");
    Assert(s1.upperTri()==sur.upperTri(),
           "SymMatrix.upperTri == SymMatrix<U,R>.upperTri");
    Assert(s1.lowerTri()==sur.lowerTri(),
           "SymMatrix.lowerTri == SymMatrix<U,R>.lowerTri");
    Assert(s1.realPart()==sur.realPart(),
           "SymMatrix.real == SymMatrix<U,R>.real");
    Assert(s1.imagPart()==sur.imagPart(),
           "SymMatrix.imag == SymMatrix<U,R>.imag");
    Assert(s1.subMatrix(N/2,N,0,N/2)==sur.subMatrix(N/2,N,0,N/2),
           "SymMatrix.subMatrix1 == SymMatrix<U,R>.subMatrix1");
    Assert(s1.subMatrix(0,N/2,N/2,N)==sur.subMatrix(0,N/2,N/2,N),
           "SymMatrix.subMatrix2 == SymMatrix<U,R>.subMatrix2");
    Assert(s1.subSymMatrix(0,N/2)==sur.subSymMatrix(0,N/2),
           "SymMatrix.subSymMatrix == SymMatrix<U,R>.subSymMatrix");
}

template <class T, tmv::UpLoType U, tmv::StorageType S> 
static void TestBasicHermMatrix_2()
{
    const int N = 6;

    if (showstartdone) {
        std::cout<<"Start TestBasicHermMatrix_2\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"U = "<<tmv::TMV_Text(U)<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
    }

    // Test assignments and constructors from arrays
    const CT quarrm[] = {
        CT(0,0), CT(2,1), CT(4,2),
                 CT(1,0), CT(3,1),
                          CT(2,0)
    };
    const CT qlarrm[] = {
        CT(0,0),
        CT(2,-1), CT(1,0),
        CT(4,-2), CT(3,-1), CT(2,0)
    };
    const CT quarcm[] = {
        CT(0,0),
        CT(2,1), CT(1,0),
        CT(4,2), CT(3,1), CT(2,0)
    };
    const CT qlarcm[] = {
        CT(0,0), CT(2,-1), CT(4,-2),
                 CT(1,0),  CT(3,-1),
                           CT(2,0)
    };

    std::vector<CT> quvecrm(6);
    for(int i=0;i<6;i++) quvecrm[i] = quarrm[i];
    std::vector<CT> qlvecrm(6);
    for(int i=0;i<6;i++) qlvecrm[i] = qlarrm[i];
    std::vector<CT> quveccm(6);
    for(int i=0;i<6;i++) quveccm[i] = quarcm[i];
    std::vector<CT> qlveccm(6);
    for(int i=0;i<6;i++) qlveccm[i] = qlarcm[i];

    tmv::HermMatrix<CT,U|S> q1(3);
    std::copy(quarrm, quarrm+6, q1.upperTri().rowmajor_begin());
    tmv::HermMatrix<CT,U|S> q2(3);
    std::copy(quarcm, quarcm+6, q2.upperTri().colmajor_begin());

    tmv::HermMatrix<CT,U|S> q3(3);
    std::copy(qlarrm, qlarrm+6, q3.lowerTri().rowmajor_begin());
    tmv::HermMatrix<CT,U|S> q4(3);
    std::copy(qlarcm, qlarcm+6, q4.lowerTri().colmajor_begin());

    tmv::HermMatrix<CT,U|S> q5(3);
    std::copy(quvecrm.begin(), quvecrm.end(), q5.upperTri().rowmajor_begin());
    tmv::HermMatrix<CT,U|S> q6(3);
    std::copy(quveccm.begin(), quveccm.end(), q6.upperTri().colmajor_begin());

    tmv::HermMatrix<CT,U|S> q7(3);
    std::copy(qlvecrm.begin(), qlvecrm.end(), q7.lowerTri().rowmajor_begin());
    tmv::HermMatrix<CT,U|S> q8(3);
    std::copy(qlveccm.begin(), qlveccm.end(), q8.lowerTri().colmajor_begin());

    tmv::HermMatrix<CT,U|S> q9x(30);
    tmv::SymMatrixView<CT> q9 = q9x.subSymMatrix(3,18,5);
    std::copy(quvecrm.begin(), quvecrm.end(), q9.upperTri().rowmajor_begin());
    tmv::HermMatrix<CT,U|S> q10x(30);
    tmv::SymMatrixView<CT> q10 = q10x.subSymMatrix(3,18,5);
    std::copy(quveccm.begin(), quveccm.end(), q10.upperTri().colmajor_begin());

    tmv::HermMatrix<CT,U|S> q11x(30);
    tmv::SymMatrixView<CT> q11 = q11x.subSymMatrix(3,18,5);
    std::copy(qlvecrm.begin(), qlvecrm.end(), q11.lowerTri().rowmajor_begin());
    tmv::HermMatrix<CT,U|S> q12x(30);
    tmv::SymMatrixView<CT> q12 = q12x.subSymMatrix(3,18,5);
    std::copy(qlveccm.begin(), qlveccm.end(), q12.lowerTri().colmajor_begin());

    // Assignment using op<< is always in rowmajor order.
    tmv::HermMatrix<CT,U|S> q13(3);
    tmv::HermMatrix<CT,U|S> q14t(3);
    tmv::SymMatrixView<CT> q14 = q14t.transpose();

    tmv::HermMatrix<CT,U|S> q15(3);
    tmv::HermMatrix<CT,U|S> q16t(3);
    tmv::SymMatrixView<CT> q16 = q16t.transpose();

    q13.upperTri() <<
        CT(0,0), CT(2,1), CT(4,2),
                 CT(1,0), CT(3,1),
                          CT(2,0);
    q14.upperTri() <<
        CT(0,0), CT(2,1), CT(4,2),
                 CT(1,0), CT(3,1),
                          CT(2,0);
    q15.lowerTri() <<
        CT(0,0),
        CT(2,-1), CT(1,0),
        CT(4,-2), CT(3,-1), CT(2,0);
    q16.lowerTri() <<
        CT(0,0),
        CT(2,-1), CT(1,0),
        CT(4,-2), CT(3,-1), CT(2,0);

    // Can also view memory directly
    CT quarrmfull[] = {
        CT(0,0),   CT(2,1),  CT(4,2),
        CT(-1,-1), CT(1,0),  CT(3,1),
        CT(-2,-2), CT(0,-1), CT(2,0)
    };
    CT quarcmfull[] = {
        CT(0,0), CT(-1,-1), CT(-2,-1),
        CT(2,1), CT(1,0),   CT(0,-1),
        CT(4,2), CT(3,1),   CT(2,0)
    };
    CT qlarrmfull[] = {
        CT(0,0),  CT(-1,1), CT(-2,2),
        CT(2,-1), CT(1,0),  CT(0,1),
        CT(4,-2), CT(3,-1), CT(2,0)
    };
    CT qlarcmfull[] = {
        CT(0,0),  CT(2,-1), CT(4,-2),
        CT(-1,1), CT(1,0),  CT(3,-1),
        CT(-2,2), CT(0,1),  CT(2,0)
    };
    CT* qarfull = (
        (S == tmv::RowMajor) ? 
        ( (U == tmv::Upper) ? quarrmfull : qlarrmfull ) : 
        ( (U == tmv::Upper) ? quarcmfull : qlarcmfull ) );
    const int Si = (S == tmv::RowMajor) ? 3 : 1;
    const int Sj = (S == tmv::RowMajor) ? 1 : 3;
    const tmv::ConstSymMatrixView<CT> q17 =
        tmv::HermMatrixViewOf(qarfull,3,U,S);
    const tmv::ConstSymMatrixView<CT> q18 =
        tmv::HermMatrixViewOf(qarfull,3,U,Si,Sj);

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
    }

    for(int i=0;i<3;i++) for(int j=0;j<3;j++) {
        CT val = i >= j ? CT(2*i-j,j-i) : CT(2*j-i,j-i);
        Assert(q1(i,j) == val,"Create HermMatrix from T* rm upperTri");
        Assert(q2(i,j) == val,"Create HermMatrix from T* cm upperTri");
        Assert(q3(i,j) == val,"Create HermMatrix from T* rm lowerTri");
        Assert(q4(i,j) == val,"Create HermMatrix from T* cm lowerTri");
        Assert(q5(i,j) == val,"Create HermMatrix from vector rm upperTri");
        Assert(q6(i,j) == val,"Create HermMatrix from vector cm upperTri");
        Assert(q7(i,j) == val,"Create HermMatrix from vector rm lowerTri");
        Assert(q8(i,j) == val,"Create HermMatrix from vector cm lowerTri");
        Assert(q9(i,j) == val,"Create HermMatrixView from vector rm upperTri");
        Assert(q10(i,j) == val,"Create HermMatrixView from vector cm upperTri");
        Assert(q11(i,j) == val,"Create HermMatrixView from vector rm lowerTri");
        Assert(q12(i,j) == val,"Create HermMatrixView from vector cm lowerTri");
        Assert(q13(i,j) == val,"Create HermMatrix from << list upperTri");
        Assert(q14(i,j) == val,"Create HermMatrixView from << list upperTri");
        Assert(q15(i,j) == val,"Create HermMatrix from << list lowerTri");
        Assert(q16(i,j) == val,"Create HermMatrixView from << list lowerTri");
        Assert(q17(i,j) == val,"Create HermMatrixView of T* (S)");
        Assert(q18(i,j) == val,"Create HermMatrixView of T* (Si,Sj)");
    }

    // Test Basic Arithmetic
    tmv::HermMatrix<std::complex<T>,U|S> h1(N);
    tmv::HermMatrix<std::complex<T>,U|S> h2(N);

    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        std::complex<T> value(T(k),T(2*k));
        std::complex<T> hvalue = i==j ? std::complex<T>(T(k),0) : value;
        if (i<=j) {
            h1(i,j) = hvalue; 
        }
        if (j<=i) {
            h2(i,j) = hvalue; 
        }
    }

    tmv::HermMatrix<std::complex<T>,U|S> h3(N);
    h3 = h1+h2;

    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        Assert(h3(i,j) == h1(i,j)+h2(i,j),"Add HermMatrices1");
        Assert(h3(j,i) == conj(h1(i,j)+h2(i,j)),"Add HermMatrices2");
    }

    h3 = h1-h2;

    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        if (i<=j) {
            Assert(h3(i,j) == h1(i,j)-h2(i,j),"Subtract HermMatrices1");
            Assert(h3(j,i) == conj(h1(i,j)-h2(i,j)),"Subtract HermMatrices2");
        }

    tmv::Matrix<std::complex<T> > n1 = h1;
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        Assert(h1(i,j) == n1(i,j),"HermMatrix -> Matrix");
    }
    tmv::HermMatrix<std::complex<T>,U|S> hn1(n1);
    Assert(h1 == hn1,"Matrix -> HermMatrix1");
    tmv::HermMatrix<std::complex<T>,U|S> hn2 = n1;
    Assert(h1 == hn2,"Matrix -> HermMatrix2");

    tmv::HermMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor> hur = h1;
    Assert(h1==hur,"HermMatrix == HermMatrix<U,R>");
    Assert(h1.view()==hur.view(),"HermMatrix.view == HermMatrix<U,R>.view");
    Assert(h1.transpose()==hur.transpose(),
           "HermMatrix.transpose == HermMatrix<U,R>.transpose");
    Assert(h1.conjugate()==hur.conjugate(),
           "HermMatrix.conjugate == HermMatrix<U,R>.conjugate");
    Assert(h1.adjoint()==hur.adjoint(),
           "HermMatrix.adjoint == HermMatrix<U,R>.adjoint");
    Assert(h1.upperTri()==hur.upperTri(),
           "HermMatrix.upperTri == HermMatrix<U,R>.upperTri");
    Assert(h1.lowerTri()==hur.lowerTri(),
           "HermMatrix.lowerTri == HermMatrix<U,R>.lowerTri");
    Assert(h1.realPart()==hur.realPart(),
           "HermMatrix.real == HermMatrix<U,R>.real");
    Assert(h1.subMatrix(N/2,N,0,N/2)==hur.subMatrix(N/2,N,0,N/2),
           "HermMatrix.subMatrix1 == HermMatrix<U,R>.subMatrix1");
    Assert(h1.subMatrix(0,N/2,N/2,N)==hur.subMatrix(0,N/2,N/2,N),
           "HermMatrix.subMatrix2 == HermMatrix<U,R>.subMatrix2");
    Assert(h1.subSymMatrix(0,N/2)==hur.subSymMatrix(0,N/2),
           "HermMatrix.subSymMatrix == HermMatrix<U,R>.subSymMatrix");
}

template <class T, tmv::UpLoType U, tmv::StorageType S> 
static void TestBasicSymMatrix_IO()
{
    const int N = 10;

    if (showstartdone) {
        std::cout<<"Start TestBasicSymMatrix_IO\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"U = "<<tmv::TMV_Text(U)<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
    }

    tmv::SymMatrix<T,U|S> s(N);
    tmv::HermMatrix<T,U|S> h(N);
    tmv::SymMatrix<std::complex<T>,U|S> cs(N);
    tmv::HermMatrix<std::complex<T>,U|S> ch(N);

    for (int i=0, k=1; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if ( (U == tmv::Upper && i<=j) ||
             (U == tmv::Lower && i>=j) ) {
            std::complex<T> value(T(k),T(k+1000));
            h(i,j) = s(i,j) = T(k);
            cs(i,j) = value;
            if (i==j) ch(i,j) = T(k);
            else ch(i,j) = value;
        }
    }
    s(1,3) = h(3,1) = T(1.e-30);
    cs(1,3) = ch(3,1) = CT(T(1.e-30),T(1.e-30));
    s(5,6) = h(6,5) = T(9.e-3);
    cs(5,6) = ch(6,5) = CT(T(9.e-3),T(9.e-3));
    cs(5,7) = ch(7,5) = CT(T(9),T(9.e-3));
    s(4,7) = h(7,4) = T(0.123456789);
    cs(4,7) = ch(7,4) = CT(T(3.123456789),T(6.987654321));

    // First check clipping function...
    tmv::SymMatrix<T> s2 = s;
    tmv::SymMatrix<CT> cs2 = cs;
    tmv::HermMatrix<T> h2 = h;
    tmv::HermMatrix<CT> ch2 = ch;
    if (!std::numeric_limits<T>::is_integer) {
        s2.clip(T(1.e-2));
        cs2.clip(T(1.e-2));
        h2.clip(T(1.e-2));
        ch2.clip(T(1.e-2));
    }
    tmv::SymMatrix<T> s3 = s;
    tmv::SymMatrix<CT> cs3 = cs;
    tmv::HermMatrix<T> h3 = h;
    tmv::HermMatrix<CT> ch3 = ch;
    s3(1,3) = h3(3,1) = T(0);
    cs3(1,3) = ch3(3,1) = T(0);
    s3(5,6) = h3(6,5) = T(0); // Others shouldn't get clipped.
    Assert(s2 == s3,"SymMatrix clip");
    Assert(cs2 == cs3,"Complex SymMatrix clip");
    Assert(h2 == h3,"HermMatrix clip");
    Assert(ch2 == ch3,"Complex HermMatrix clip");

    // However, ThreshIO for complex works slightly differently than clip.
    // It clips _either_ the real or imag component, so now cm2(5,6) and
    // cm2(6,6) need to be modified.
    cs2(5,6) = cs3(5,6) = ch2(6,5) = ch3(6,5) = T(0);
    cs2(5,7) = cs3(5,7) = ch2(7,5) = ch3(7,5) = T(9);

    // Write matrices with 4 different style
    std::ofstream fout("tmvtest_symmatrix_io.dat");
    Assert(bool(fout),"Couldn't open tmvtest_symmatrix_io.dat for output");
    fout << s << std::endl;
    fout << h << std::endl;
    fout << cs << std::endl;
    fout << ch << std::endl;
    fout << tmv::CompactIO() << s << std::endl;
    fout << tmv::CompactIO() << h << std::endl;
    fout << tmv::CompactIO() << cs << std::endl;
    fout << tmv::CompactIO() << ch << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << s << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << h << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << cs << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << ch << std::endl;
    tmv::IOStyle myStyle =
        tmv::CompactIO().setThresh(1.e-2).setPrecision(4).
        markup("Start","[",",","]","---","Done");
    fout << myStyle << s << std::endl;
    fout << myStyle << h << std::endl;
    fout << myStyle << cs << std::endl;
    fout << myStyle << ch << std::endl;
    fout.close();

    // When using (the default) prec(6), these will be the values read in.
    s(4,7) = h(7,4) = T(0.123457);
    cs(4,7) = ch(7,4) = CT(T(3.12346),T(6.98765));

    // When using prec(12), the full correct values will be read in.

    // When using prec(4), these will be the values read in.
    s3(4,7) = h3(7,4) = T(0.1235);
    if (std::numeric_limits<T>::is_integer) cs3(4,7) = ch3(7,4) = CT(3,6);
    else cs3(4,7) = ch3(7,4)  = CT(T(3.123),T(6.988));

    // Read them back in
    tmv::SymMatrix<T,tmv::Upper|tmv::RowMajor> xs1(N);
    tmv::HermMatrix<T,tmv::Upper|tmv::RowMajor> xh1(N);
    tmv::SymMatrix<CT,tmv::Upper|tmv::RowMajor> xcs1(N);
    tmv::HermMatrix<CT,tmv::Upper|tmv::RowMajor> xch1(N);
    std::ifstream fin("tmvtest_symmatrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_symmatrix_io.dat for input");
    fin >> xs1 >> xh1 >> xcs1 >> xch1;
    Assert(EqualIO(s,xs1,EPS),"SymMatrix I/O check normal");
    Assert(EqualIO(h,xh1,EPS),"HermMatrix I/O check normal");
    Assert(EqualIO(cs,xcs1,EPS),"CSymMatrix I/O check normal");
    Assert(EqualIO(ch,xch1,EPS),"CHermMatrix I/O check normal");
    fin >> tmv::CompactIO() >> xs1 >> tmv::CompactIO() >> xh1;
    fin >> tmv::CompactIO() >> xcs1 >> tmv::CompactIO() >> xch1;
    Assert(EqualIO(s,xs1,EPS),"SymMatrix I/O check compact");
    Assert(EqualIO(h,xh1,EPS),"HermMatrix I/O check compact");
    Assert(EqualIO(cs,xcs1,EPS),"CSymMatrix I/O check compact");
    Assert(EqualIO(ch,xch1,EPS),"CHermMatrix I/O check compact");
    fin >> xs1.view() >> xh1.view() >> xcs1.view() >> xch1.view();
    Assert(EqualIO(s2,xs1,EPS),"SymMatrix I/O check thresh");
    Assert(EqualIO(h2,xh1,EPS),"HermMatrix I/O check thresh");
    Assert(EqualIO(cs2,xcs1,EPS),"CSymMatrix I/O check thresh");
    Assert(EqualIO(ch2,xch1,EPS),"CHermMatrix I/O check thresh");
    fin >> myStyle >> xs1.view() >> myStyle >> xh1.view();
    fin >> myStyle >> xcs1.view() >> myStyle >> xch1.view();
    Assert(EqualIO(s3,xs1,EPS),"SymMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(h3,xh1,EPS),"HermMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(cs3,xcs1,EPS),"CSymMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(ch3,xch1,EPS),"CHermMatrix I/O check compact thresh & prec(4)");
    fin.close();

    // Repeat for column major
    tmv::SymMatrix<T,tmv::Upper|tmv::ColMajor> xs2(N);
    tmv::HermMatrix<T,tmv::Upper|tmv::ColMajor> xh2(N);
    tmv::SymMatrix<CT,tmv::Upper|tmv::ColMajor> xcs2(N);
    tmv::HermMatrix<CT,tmv::Upper|tmv::ColMajor> xch2(N);
    fin.open("tmvtest_symmatrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_symmatrix_io.dat for input");
    fin >> xs2.view() >> xh2.view() >> xcs2.view() >> xch2.view();
    Assert(EqualIO(s,xs2,EPS),"SymMatrix I/O check normal");
    Assert(EqualIO(h,xh2,EPS),"HermMatrix I/O check normal");
    Assert(EqualIO(cs,xcs2,EPS),"CSymMatrix I/O check normal");
    Assert(EqualIO(ch,xch2,EPS),"CHermMatrix I/O check normal");
    fin >> tmv::CompactIO() >> xs2.view() >> tmv::CompactIO() >> xh2.view();
    fin >> tmv::CompactIO() >> xcs2.view() >> tmv::CompactIO() >> xch2.view();
    Assert(EqualIO(s,xs2,EPS),"SymMatrix I/O check compact");
    Assert(EqualIO(h,xh2,EPS),"HermMatrix I/O check compact");
    Assert(EqualIO(cs,xcs2,EPS),"CSymMatrix I/O check compact");
    Assert(EqualIO(ch,xch2,EPS),"CHermMatrix I/O check compact");
    fin >> xs2 >> xh2 >> xcs2 >> xch2;
    Assert(EqualIO(s2,xs2,EPS),"SymMatrix I/O check thresh");
    Assert(EqualIO(h2,xh2,EPS),"HermMatrix I/O check thresh");
    Assert(EqualIO(cs2,xcs2,EPS),"CSymMatrix I/O check thresh");
    Assert(EqualIO(ch2,xch2,EPS),"CHermMatrix I/O check thresh");
    fin >> myStyle >> xs2 >> myStyle >> xh2;
    fin >> myStyle >> xcs2 >> myStyle >> xch2;
    Assert(EqualIO(s3,xs2,EPS),"SymMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(h3,xh2,EPS),"HermMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(cs3,xcs2,EPS),"CSymMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(ch3,xch2,EPS),"CHermMatrix I/O check compact thresh & prec(4)");
    fin.close();

    // Repeat for Lower Storage
    // Also switch s and h for real matrices.  They should be able to 
    // read the opposite letter for the compact storage.
    tmv::SymMatrix<T,tmv::Lower|tmv::RowMajor> xs3(N);
    tmv::HermMatrix<T,tmv::Lower|tmv::RowMajor> xh3(N);
    tmv::SymMatrix<CT,tmv::Lower|tmv::RowMajor> xcs3(N);
    tmv::HermMatrix<CT,tmv::Lower|tmv::RowMajor> xch3(N);
    fin.open("tmvtest_symmatrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_symmatrix_io.dat for input");
    fin >> xh3.view() >> xs3.view() >> xcs3.view() >> xch3.view();
    Assert(EqualIO(s,xh3,EPS),"SymMatrix I/O check normal");
    Assert(EqualIO(h,xs3,EPS),"HermMatrix I/O check normal");
    Assert(EqualIO(cs,xcs3,EPS),"CSymMatrix I/O check normal");
    Assert(EqualIO(ch,xch3,EPS),"CHermMatrix I/O check normal");
    fin >> tmv::CompactIO() >> xh3.view() >> tmv::CompactIO() >> xs3.view();
    fin >> tmv::CompactIO() >> xcs3.view() >> tmv::CompactIO() >> xch3.view();
    Assert(EqualIO(s,xh3,EPS),"SymMatrix I/O check compact");
    Assert(EqualIO(h,xs3,EPS),"HermMatrix I/O check compact");
    Assert(EqualIO(cs,xcs3,EPS),"CSymMatrix I/O check compact");
    Assert(EqualIO(ch,xch3,EPS),"CHermMatrix I/O check compact");
    fin >> xh3 >> xs3 >> xcs3 >> xch3;
    Assert(EqualIO(s2,xh3,EPS),"SymMatrix I/O check thresh");
    Assert(EqualIO(h2,xs3,EPS),"HermMatrix I/O check thresh");
    Assert(EqualIO(cs2,xcs3,EPS),"CSymMatrix I/O check thresh");
    Assert(EqualIO(ch2,xch3,EPS),"CHermMatrix I/O check thresh");
    fin >> myStyle >> xh3 >> myStyle >> xs3;
    fin >> myStyle >> xcs3 >> myStyle >> xch3;
    Assert(EqualIO(s3,xh3,EPS),"SymMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(h3,xs3,EPS),"HermMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(cs3,xcs3,EPS),"CSymMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(ch3,xch3,EPS),"CHermMatrix I/O check compact thresh & prec(4)");
    fin.close();

    tmv::SymMatrix<T,tmv::Lower|tmv::ColMajor> xs4(N);
    tmv::HermMatrix<T,tmv::Lower|tmv::ColMajor> xh4(N);
    tmv::SymMatrix<CT,tmv::Lower|tmv::ColMajor> xcs4(N);
    tmv::HermMatrix<CT,tmv::Lower|tmv::ColMajor> xch4(N);
    fin.open("tmvtest_symmatrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_symmatrix_io.dat for input");
    fin >> xh4.view() >> xs4.view() >> xcs4.view() >> xch4.view();
    Assert(EqualIO(s,xh4,EPS),"SymMatrix I/O check normal");
    Assert(EqualIO(h,xs4,EPS),"HermMatrix I/O check normal");
    Assert(EqualIO(cs,xcs4,EPS),"CSymMatrix I/O check normal");
    Assert(EqualIO(ch,xch4,EPS),"CHermMatrix I/O check normal");
    fin >> tmv::CompactIO() >> xh4.view() >> tmv::CompactIO() >> xs4.view();
    fin >> tmv::CompactIO() >> xcs4.view() >> tmv::CompactIO() >> xch4.view();
    Assert(EqualIO(s,xh4,EPS),"SymMatrix I/O check compact");
    Assert(EqualIO(h,xs4,EPS),"HermMatrix I/O check compact");
    Assert(EqualIO(cs,xcs4,EPS),"CSymMatrix I/O check compact");
    Assert(EqualIO(ch,xch4,EPS),"CHermMatrix I/O check compact");
    fin >> xh4 >> xs4 >> xcs4 >> xch4;
    Assert(EqualIO(s2,xh4,EPS),"SymMatrix I/O check thresh");
    Assert(EqualIO(h2,xs4,EPS),"HermMatrix I/O check thresh");
    Assert(EqualIO(cs2,xcs4,EPS),"CSymMatrix I/O check thresh");
    Assert(EqualIO(ch2,xch4,EPS),"CHermMatrix I/O check thresh");
    fin >> myStyle >> xh4 >> myStyle >> xs4;
    fin >> myStyle >> xcs4 >> myStyle >> xch4;
    Assert(EqualIO(s3,xh4,EPS),"SymMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(h3,xs4,EPS),"HermMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(cs3,xcs4,EPS),"CSymMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(ch3,xch4,EPS),"CHermMatrix I/O check compact thresh & prec(4)");
    fin.close();

    // And repeat for matrices that need to be resized.
    // Also check switching the default IOStyle.
    tmv::CompactIO().makeDefault();
    tmv::SymMatrix<T> zs1,zs2,zs3,zs4;
    tmv::HermMatrix<T> zh1,zh2,zh3,zh4;
    tmv::SymMatrix<CT> zcs1,zcs2,zcs3,zcs4;
    tmv::HermMatrix<CT> zch1,zch2,zch3,zch4;
    fin.open("tmvtest_symmatrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_symmatrix_io.dat for input");
    fin >> tmv::NormalIO() >> zs1 >> tmv::NormalIO() >> zh1;
    fin >> tmv::NormalIO() >> zcs1 >> tmv::NormalIO() >> zch1;
    Assert(EqualIO(s,zs1,EPS),"SymMatrix I/O check normal");
    Assert(EqualIO(h,zh1,EPS),"HermMatrix I/O check normal");
    Assert(EqualIO(cs,zcs1,EPS),"CSymMatrix I/O check normal");
    Assert(EqualIO(ch,zch1,EPS),"CHermMatrix I/O check normal");
    fin >> zs2 >> zh2 >> zcs2 >> zch2;
    Assert(EqualIO(s,zs2,EPS),"SymMatrix I/O check compact");
    Assert(EqualIO(h,zh2,EPS),"HermMatrix I/O check compact");
    Assert(EqualIO(cs,zcs2,EPS),"CSymMatrix I/O check compact");
    Assert(EqualIO(ch,zch2,EPS),"CHermMatrix I/O check compact");
    fin >> tmv::NormalIO() >> zs3 >> tmv::NormalIO() >> zh3;
    fin >> tmv::NormalIO() >> zcs3 >> tmv::NormalIO() >> zch3;
    Assert(EqualIO(s2,zs3,EPS),"SymMatrix I/O check thresh");
    Assert(EqualIO(h2,zh3,EPS),"HermMatrix I/O check thresh");
    Assert(EqualIO(cs2,zcs3,EPS),"CSymMatrix I/O check thresh");
    Assert(EqualIO(ch2,zch3,EPS),"CHermMatrix I/O check thresh");
    fin >> myStyle >> zs4 >> myStyle >> zh4;
    fin >> myStyle >> zcs4 >> myStyle >> zch4;
    Assert(EqualIO(s3,zs4,EPS),"SymMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(h3,zh4,EPS),"HermMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(cs3,zcs4,EPS),"CSymMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(ch3,zch4,EPS),"CHermMatrix I/O check compact thresh & prec(4)");
    fin.close();
    tmv::IOStyle::revertDefault();

    // Finally, check that the NormalIO can be read in as a regular matrix.
    tmv::Matrix<T> zm1,zm2;
    tmv::Matrix<CT> zcm1,zcm2;
    fin.open("tmvtest_symmatrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_symmatrix_io.dat for input");
    fin >> zm1 >> zm2 >> zcm1 >> zcm2;
    Assert(EqualIO(s,zm1,EPS),"SymMatrix -> Matrix I/O check");
    Assert(EqualIO(h,zm2,EPS),"HermMatrix -> Matrix I/O check");
    Assert(EqualIO(cs,zcm1,EPS),"CSymMatrix -> CMatrix I/O check");
    Assert(EqualIO(ch,zcm2,EPS),"CHermMatrix -> CMatrix I/O check");
    fin.close();

#if XTEST == 0
    std::remove("tmvtest_symmatrix_io.dat");
#endif
}

template <class T, tmv::UpLoType U, tmv::StorageType S> 
static void TestBasicSymMatrix()
{
    TestBasicSymMatrix_1<T,U,S>();
    TestBasicSymMatrix_2<T,U,S>();
    TestBasicHermMatrix_1<T,U,S>();
    TestBasicHermMatrix_2<T,U,S>();
    TestBasicSymMatrix_IO<T,U,S>();
}

template <class T> void TestSymMatrix() 
{
    TestBasicSymMatrix<T,tmv::Upper,tmv::ColMajor>();
    TestBasicSymMatrix<T,tmv::Lower,tmv::ColMajor>();
#if (XTEST & 2)
    TestBasicSymMatrix<T,tmv::Upper,tmv::RowMajor>();
    TestBasicSymMatrix<T,tmv::Lower,tmv::RowMajor>();
#endif

    std::cout<<"SymMatrix<"<<tmv::TMV_Text(T())<<"> passed all basic tests\n";

    TestSymMatrixArith_A<T>();
    std::cout<<"SymMatrix<"<<tmv::TMV_Text(T())<<
        "> (Sym/Sym) Arithmetic passed all tests\n";
    TestSymMatrixArith_B1<T>();
    TestSymMatrixArith_B2<T>();
    std::cout<<"SymMatrix<"<<tmv::TMV_Text(T())<<
        "> (Matrix/Sym) Arithmetic passed all tests\n";
    TestSymMatrixArith_C1<T>();
    TestSymMatrixArith_C2<T>();
    std::cout<<"SymMatrix<"<<tmv::TMV_Text(T())<<
        "> (Diag/Sym) Arithmetic passed all tests\n";
    TestSymMatrixArith_D1<T>();
    TestSymMatrixArith_D2<T>();
    std::cout<<"SymMatrix<"<<tmv::TMV_Text(T())<<
        "> (Tri/Sym) Arithmetic passed all tests\n";
    TestSymMatrixArith_E1<T>();
    TestSymMatrixArith_E2<T>();
    std::cout<<"SymMatrix<"<<tmv::TMV_Text(T())<<
        "> (Band/Sym) Arithmetic passed all tests\n";
}

#ifdef TEST_DOUBLE
template void TestSymMatrix<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymMatrix<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymMatrix<long double>();
#endif
#ifdef TEST_INT
template void TestSymMatrix<int>();
#endif
