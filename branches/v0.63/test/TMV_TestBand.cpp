
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_SymBand.h"
#include <fstream>
#include <cstdio>

template <class T, tmv::StorageType S> 
static void TestBasicBandMatrix_1()
{
    const int N = 10;
    const int nhi = 1;
    const int nlo = 3;

    tmv::BandMatrix<T> a1(N,N,nlo,nhi);
    tmv::BandMatrix<T,tmv::RowMajor,tmv::FortranStyle> a1f(N,N,nlo,nhi);

    Assert(a1.colsize() == size_t(N) && a1.rowsize() == size_t(N),
           "Creating BandMatrix(N)");
    Assert(a1.nlo() == nlo && a1.nhi() == nhi,"Creating BandMatrix(nlo,nhi)");
    Assert(a1f.colsize() == size_t(N) && a1f.rowsize() == size_t(N),
           "Creating BandMatrixF(N)");
    Assert(a1f.nlo() == nlo && a1f.nhi() == nhi,
           "Creating BandMatrixF(nlo,nhi)");

    for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k)
        if ( j <= i + nhi && i <= j + nlo) {
            a1(i,j) = T(k);
            a1f(i+1,j+1) = T(k);
        }

    tmv::ConstBandMatrixView<T> a1cv = a1.View();
    tmv::BandMatrixView<T> a1v = a1.View();
    tmv::ConstBandMatrixView<T,tmv::FortranStyle> a1fcv = a1f.View();
    tmv::BandMatrixView<T,tmv::FortranStyle> a1fv = a1f.View();

    for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if ( j <= i + nhi && i <= j + nlo) {
            int j1 = 0;  if (i > nlo) j1 = i-nlo;
            int j2 = i+nhi+1;  if (j2 >= N) j2 = N;
            int i1 = 0;  if (j > nhi) i1 = j-nhi;
            int i2 = j+nlo+1;  if (i2 >= N) i2 = N;
            Assert(a1(i,j) == k,"Read/Write BandMatrix");
            Assert(a1cv(i,j) == k,"Access BandMatrix CV");
            Assert(a1v(i,j) == k,"Access BandMatrix V");
            Assert(a1f(i+1,j+1) == k,"Read/Write BandMatrixF");
            Assert(a1fcv(i+1,j+1) == k,"Access BandMatrixF CV");
            Assert(a1fv(i+1,j+1) == k,"Access BandMatrixF V");
            Assert(a1.row(i,j1,j2)(j-j1) == k,"BandMatrix.row");
            Assert(a1cv.row(i,j1,j2)(j-j1) == k,"BandMatrix.row CV");
            Assert(a1v.row(i,j1,j2)(j-j1) == k,"BandMatrix.row V");
            Assert(a1f.row(i+1,j1+1,j2)(j+1-j1) == k,"BandMatrixF.row");
            Assert(a1fcv.row(i+1,j1+1,j2)(j+1-j1) == k,"BandMatrixF.row CV");
            Assert(a1fv.row(i+1,j1+1,j2)(j+1-j1) == k,"BandMatrixF.row V");
            Assert(a1.col(j,i1,i2)(i-i1) == k,"BandMatrix.col");
            Assert(a1cv.col(j,i1,i2)(i-i1) == k,"BandMatrix.col CV");
            Assert(a1v.col(j,i1,i2)(i-i1) == k,"BandMatrix.col V");
            Assert(a1f.col(j+1,i1+1,i2)(i+1-i1) == k,"BandMatrixF.col");
            Assert(a1fcv.col(j+1,i1+1,i2)(i+1-i1) == k,"BandMatrixF.col CV");
            Assert(a1fv.col(j+1,i1+1,i2)(i+1-i1) == k,"BandMatrixF.col V");
            int d = j-i;
            if (d>0) {
                Assert(a1.diag(d)(i) == k,"BandMatrix.diag1");
                Assert(a1cv.diag(d)(i) == k,"BandMatrix.diag1 CV");
                Assert(a1v.diag(d)(i) == k,"BandMatrix.diag1 V");
                Assert(a1.diag(d,i,N-d)(0) == k,"BandMatrix.diag2");
                Assert(a1cv.diag(d,i,N-d)(0) == k,"BandMatrix.diag2 CV");
                Assert(a1v.diag(d,i,N-d)(0) == k,"BandMatrix.diag2 V");
                Assert(a1f.diag(d)(i+1) == k,"BandMatrixF.diag1");
                Assert(a1fcv.diag(d)(i+1) == k,"BandMatrixF.diag1 CV");
                Assert(a1fv.diag(d)(i+1) == k,"BandMatrixF.diag1 V");
                Assert(a1f.diag(d,i+1,N-d)(1) == k,"BandMatrixF.diag2");
                Assert(a1fcv.diag(d,i+1,N-d)(1) == k,"BandMatrixF.diag2 CV");
                Assert(a1fv.diag(d,i+1,N-d)(1) == k,"BandMatrixF.diag2 V");
            } else {
                if (d==0) {
                    Assert(a1.diag()(j) == k,"BandMatrix.diag");
                    Assert(a1cv.diag()(j) == k,"BandMatrix.diag CV");
                    Assert(a1v.diag()(j) == k,"BandMatrix.diag V");
                    Assert(a1f.diag()(j+1) == k,"BandMatrixF.diag");
                    Assert(a1fcv.diag()(j+1) == k,"BandMatrixF.diag CV");
                    Assert(a1fv.diag()(j+1) == k,"BandMatrixF.diag V");
                }
                Assert(a1.diag(d)(j) == k,"BandMatrix.diag1");
                Assert(a1cv.diag(d)(j) == k,"BandMatrix.diag1 CV");
                Assert(a1v.diag(d)(j) == k,"BandMatrix.diag1 V");
                Assert(a1.diag(d,j,N+d)(0) == k,"BandMatrix.diag2");
                Assert(a1cv.diag(d,j,N+d)(0) == k,"BandMatrix.diag2 CV");
                Assert(a1v.diag(d,j,N+d)(0) == k,"BandMatrix.diag2 V");
                Assert(a1f.diag(d)(j+1) == k,"BandMatrixF.diag1");
                Assert(a1fcv.diag(d)(j+1) == k,"BandMatrixF.diag1 CV");
                Assert(a1fv.diag(d)(j+1) == k,"BandMatrixF.diag1 V");
                Assert(a1f.diag(d,j+1,N+d)(1) == k,"BandMatrixF.diag2");
                Assert(a1fcv.diag(d,j+1,N+d)(1) == k,"BandMatrixF.diag2 CV");
                Assert(a1fv.diag(d,j+1,N+d)(1) == k,"BandMatrixF.diag2 V");
            }
        }
    }
    Assert(a1 == a1f,"CStyle BandMatrix == FortranStyle BandMatrix");
    Assert(a1 == a1cv,"BandMatrix == ConstBandMatrixView");
    Assert(a1 == a1v,"BandMatrix == BandMatrixView");
    Assert(a1 == a1fcv,"BandMatrix == FortranStyle ConstBandMatrixView");
    Assert(a1 == a1fv,"BandMatrix == FortranStyle BandMatrixView");
}

template <class T, tmv::StorageType S> 
static void TestBasicBandMatrix_2()
{
    const int N = 10;
    const int nhi = 1;
    const int nlo = 3;

    tmv::BandMatrix<T> a1(N,N,nlo,nhi);

    for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k)
        if ( j <= i + nhi && i <= j + nlo) {
            a1(i,j) = T(k);
        }

    std::vector<T> qv;
    if (S == tmv::RowMajor) {
        const T qvar[] = { T(0), T(-1), T(-2),
            T(2), T(1),  T(0),  T(-1),
            T(3),  T(2),  T(1),  T(0) };
        qv.resize(11);
        for(size_t i=0;i<qv.size();i++) qv[i] = qvar[i];
    } else if (S == tmv::ColMajor) {
        const T qvar[] = { T(0),  T(2),
            T(888),  T(-1), T(1),  T(3),
            T(-2), T(0),  T(2), T(888),
            T(-1), T(1), T(888), T(888),
            T(0) };
        qv.resize(15);
        for(size_t i=0;i<qv.size();i++) qv[i] = qvar[i];
    } else {
        const T qvar[] = {        T(2),  T(3),
            T(0),  T(1),  T(2),
            T(-1), T(0),  T(1),
            T(-2), T(-1), T(0) };
        qv.resize(11);
        for(size_t i=0;i<qv.size();i++) qv[i] = qvar[i];
    }
    T qar[15];
    for(size_t i=0;i<qv.size();i++) qar[i] = qv[i];
    tmv::BandMatrix<T,S> q1(3,5,1,2,qar);
    tmv::BandMatrix<T,S> q2(3,5,1,2,qv);
    tmv::BandMatrixView<T> q3 = tmv::BandMatrixViewOf(qar,3,5,1,2,S);
    for(int i=0;i<3;i++) for(int j=0;j<5;j++) if (j<=i+2 && i<=j+1) {
        Assert(q1(i,j) == T(2*i-j),"Create BandMatrix from T*");
        Assert(q2(i,j) == T(2*i-j),"Create BandMatrix from vector");
        Assert(q3(i,j) == T(2*i-j),"Create BandMatrixView of T*");
    }

    tmv::BandMatrix<T> a2(N,N,nlo,nhi);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        if ( j <= i + nhi && i <= j + nlo) {
            a1(i,j) = T(3+i-5*j);
            a2(i,j) = T(5-2*i+4*j);
        }
    tmv::BandMatrix<T,tmv::RowMajor,tmv::FortranStyle> a1f(N,N,nlo,nhi);
    a1f = a1;
    Assert(a1f == a1,"Copy CStyle BandMatrix to FortranStyle");

    tmv::BandMatrix<T> c(N,N,nlo,nhi);
    c = a1+a2;
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        if ( j <= i + nhi && i <= j + nlo) 
            Assert(c(i,j) == T(8-i-j),"Add BandMatrices");

    c = a1-a2;
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        if ( j <= i + nhi && i <= j + nlo) 
            Assert(c(i,j) == T(-2+3*i-9*j),"Subtract BandMatrices");

    tmv::Matrix<T> m1 = a1;
    for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k)
        if ( j <= i + nhi && i <= j + nlo) 
            Assert(a1(i,j) == m1(i,j),"BandMatrix -> Matrix");
    Assert(a1 == tmv::BandMatrix<T>(m1,nlo,nhi),"Matrix -> BandMatrix");
}

template <class T, tmv::StorageType S> 
static void TestBasicBandMatrix_IO()
{
    const int N = 10;
    const int nhi = 1;
    const int nlo = 3;

    tmv::BandMatrix<T> a1(N,N,nlo,nhi);

    for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k)
        if ( j <= i + nhi && i <= j + nlo) {
            a1(i,j) = T(k);
        }

    tmv::BandMatrix<std::complex<T>,S> ca1 = a1*std::complex<T>(1,2);
    std::ofstream fout("tmvtest_bandmatrix_io.dat");
    if (!fout) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_bandmatrix_io.dat for output\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_bandmatrix_io.dat for output");
#endif
    }
    fout << ca1 << std::endl;
    ca1.WriteCompact(fout);
    fout.close();

    tmv::Matrix<std::complex<T>,tmv::RowMajor> xm1(N,N);
    tmv::BandMatrix<std::complex<T>,tmv::RowMajor> xb1(N,N,nlo,nhi);
    std::ifstream fin("tmvtest_bandmatrix_io.dat");
    if (!fin) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_bandmatrix_io.dat for input\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_bandmatrix_io.dat for input");
#endif
    }
    fin >> xm1 >> xb1;
    fin.close();
    Assert(tmv::Matrix<std::complex<T> >(ca1) == xm1,"BandMatrix I/O check #1");
    Assert(ca1 == xb1,"BandMatrix Compact I/O check #1");

    tmv::Matrix<std::complex<T>,tmv::ColMajor> xm2(N,N);
    tmv::BandMatrix<std::complex<T>,tmv::ColMajor> xb2(N,N,nlo,nhi);
    fin.open("tmvtest_bandmatrix_io.dat");
    fin >> xm2 >> xb2;
    fin.close();
    Assert(tmv::Matrix<std::complex<T> >(ca1) == xm2,"BandMatrix I/O check #2");
    Assert(ca1 == xb2,"BandMatrix Compact I/O check #2");

    tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> xb3(N,N,nlo,nhi);
    fin.open("tmvtest_bandmatrix_io.dat");
    fin >> xm1 >> xb3;
    fin.close();
    Assert(ca1 == xb3,"BandMatrix Compact I/O check #3");

    std::auto_ptr<tmv::Matrix<std::complex<T> > > xm4;
    std::auto_ptr<tmv::BandMatrix<std::complex<T> > > xb4;
    fin.open("tmvtest_bandmatrix_io.dat");
    fin >> xm4 >> xb4;
    fin.close();
    Assert(tmv::Matrix<std::complex<T> >(ca1) == *xm4,"BandMatrix I/O check #4");
    Assert(ca1 == *xb4,"BandMatrix Compact I/O check #4");

#ifndef XTEST
    std::remove("tmvtest_bandmatrix_io.dat");
#endif
}

template <class T, tmv::StorageType S> 
static void TestBasicBandMatrix()
{
    TestBasicBandMatrix_1<T,S>();
    TestBasicBandMatrix_2<T,S>();
    TestBasicBandMatrix_IO<T,S>();
}

template <class T> 
void TestBandMatrix()
{
    TestBasicBandMatrix<T,tmv::RowMajor>();
    TestBasicBandMatrix<T,tmv::ColMajor>();
    TestBasicBandMatrix<T,tmv::DiagMajor>();

    std::cout<<"BandMatrix<"<<tmv::TMV_Text(T())<<"> passed all basic tests\n";

    if (tmv::TMV_Epsilon<T>() > T(0)) {
        TestBandMatrixArith_A<T>();
        std::cout<<"BandMatrix<"<<tmv::TMV_Text(T())<<
            "> (Band/Band) Arithmetic passed all tests\n";
        TestBandMatrixArith_B1<T>();
        TestBandMatrixArith_B2<T>();
        std::cout<<"BandMatrix<"<<tmv::TMV_Text(T())<<
            "> (Matrix/Band) Arithmetic passed all tests\n";
        TestBandMatrixArith_C1<T>();
        TestBandMatrixArith_C2<T>();
        std::cout<<"BandMatrix<"<<tmv::TMV_Text(T())<<
            "> (Diag/Band) Arithmetic passed all tests\n";
        TestBandMatrixArith_D1<T>();
        TestBandMatrixArith_D2<T>();
        std::cout<<"BandMatrix<"<<tmv::TMV_Text(T())<<
            "> (Tri/Band) Arithmetic passed all tests\n";
    }
}

#ifdef INST_DOUBLE
template void TestBandMatrix<double>();
#endif
#ifdef INST_FLOAT
template void TestBandMatrix<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestBandMatrix<long double>();
#endif
#ifdef INST_INT
template void TestBandMatrix<int>();
#endif
