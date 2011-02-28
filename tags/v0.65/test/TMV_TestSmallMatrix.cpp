#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include <fstream>
#include <cstdio>

template <class T, int M, int N, tmv::StorageType S> 
inline void TestBasicSmallMatrix()
{
    tmv::SmallMatrix<T,M,N,S> m;
    tmv::SmallMatrix<T,M,N,S,tmv::FortranStyle> mf;
    Assert(m.colsize() == size_t(M) && m.rowsize() == size_t(N),
           "Creating SmallMatrix(M,N)");
    Assert(m.colsize() == size_t(M) && m.rowsize() == size_t(N),
           "Creating SmallMatrixF(M,N)");

    for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
        m(i,j) = T(k);
        mf(i+1,j+1) = T(k);
    }

    for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
        Assert(m(i,j) == k,"Read/Write SmallMatrix");
        Assert(mf(i+1,j+1) == k,"Read/Write SmallMatrixF");
        Assert(m[i][j] == k,"[] style access of SmallMatrix");
        Assert(mf[i+1][j+1] == k,"[] style access of SmallMatrixF");
        Assert(m.row(i)(j) == k,"SmallMatrix.row");
        Assert(mf.row(i+1)(j+1) == k,"SmallMatrixF.row");
        Assert(m.row(i,j,N)(0) == k,"SmallMatrix.row2");
        Assert(mf.row(i+1,j+1,N)(1) == k,"SmallMatrixF.row2");
        Assert(m.col(j)(i) == k,"SmallMatrix.col");
        Assert(mf.col(j+1)(i+1) == k,"SmallMatrixF.col");
        Assert(m.col(j,i,M)(0) == k,"SmallMatrix.col2");
        Assert(mf.col(j+1,i+1,M)(1) == k,"SmallMatrixF.col2");
        if (i<j) {
            Assert(m.diag(j-i)(i) == k,"SmallMatrix.diag");
            Assert(mf.diag(j-i)(i+1) == k,"SmallMatrixF.diag");
            Assert(m.diag(j-i,i,N-j+i)(0) == k,"SmallMatrix.diag2");
            Assert(mf.diag(j-i,i+1,N-j+i)(1) == k,"SmallMatrix.diag2");
        } else {
            if (i==j) {
                Assert(m.diag()(i) == k,"SmallMatrix.diag");
                Assert(mf.diag()(i+1) == k,"SmallMatrixF.diag");
            }
            Assert(m.diag(j-i)(j) == k,"SmallMatrix.diag1");
            Assert(mf.diag(j-i)(j+1) == k,"SmallMatrixF.diag1");
            if (N+i-j > M) {
                Assert(m.diag(j-i,j,M+j-i)(0) == k,"SmallMatrix.diag2");
                Assert(mf.diag(j-i,j+1,M+j-i)(1) == k,"SmallMatrix.diag2");
            } else {
                Assert(m.diag(j-i,j,N)(0) == k,"SmallMatrix.diag2");
                Assert(mf.diag(j-i,j+1,N)(1) == k,"SmallMatrix.diag2");
            }
        }
    }
    Assert(m == mf,"CStyle SmallMatrix == FortranStyle SmallMatrix");

    if (S == tmv::RowMajor) {
        const T qar[] = { 
            T(0), T(-1), T(-2),
            T(2), T(1), T(0) };
        tmv::SmallMatrix<T,2,3,S> q1(qar);
        for(int i=0;i<2;i++) for(int j=0;j<3;j++) {
            Assert(q1(i,j) == T(2*i-j),"Create SmallMatrix from T*");
        }
        std::vector<T> qv(6);
        for(int i=0;i<6;i++) qv[i] = qar[i];
        tmv::SmallMatrix<T,2,3,S> q2(qv);
        for(int i=0;i<2;i++) for(int j=0;j<3;j++) {
            Assert(q2(i,j) == T(2*i-j),"Create SmallMatrix from vector");
        }
        tmv::SmallMatrix<T,2,3,S> q3;
        q3 <<
            0, -1, -2,
            2, 1, 0;
        for(int i=0;i<2;i++) for(int j=0;j<3;j++) {
            Assert(q3(i,j) == T(2*i-j),"Create SmallMatrix from <<");
        }
    } else {
        const T qar[] = { 
            T(0), T(2),
            T(-1), T(1),
            T(-2), T(0) };
        tmv::SmallMatrix<T,2,3,S> q1(qar);
        for(int i=0;i<2;i++) for(int j=0;j<3;j++) {
            Assert(q1(i,j) == T(2*i-j),"Create SmallMatrix from T*");
        }
        std::vector<T> qv(6);
        for(int i=0;i<6;i++) qv[i] = qar[i];
        tmv::SmallMatrix<T,2,3,S> q2(qv);
        for(int i=0;i<2;i++) for(int j=0;j<3;j++) {
            Assert(q2(i,j) == T(2*i-j),"Create SmallMatrix from vector");
        }
        tmv::SmallMatrix<T,2,3,S> q3;
        q3 <<
            0, 2,
            -1, 1,
            -2, 0;
        for(int i=0;i<2;i++) for(int j=0;j<3;j++) {
            Assert(q3(i,j) == T(2*i-j),"Create SmallMatrix from <<");
        }
    }
    // Test Basic Arithmetic 

    tmv::SmallMatrix<T,M,N,S> a;
    tmv::SmallMatrix<T,M,N,S> b;
    tmv::SmallMatrix<T,M,N,S> c;
    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) {
        a(i,j) = T(3+i+5*j);
        b(i,j) = T(5+2*i+4*j);
    }
    mf = a;
    Assert(a == mf,"Copy CStyle SmallMatrix to FortranStyle");

    c = a+b;
    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
        Assert(c(i,j) == T(8+3*i+9*j),"Add Matrices");

    c = a-b;
    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
        Assert(c(i,j) == T(-2-i+j),"Subtract Matrices");

    tmv::SmallMatrix<std::complex<T>,M,N,S> cm;
    tmv::SmallMatrix<std::complex<T>,M,N,S> ca;
    tmv::SmallMatrix<std::complex<T>,M,N,S> cb;
    Assert(cm.colsize() == size_t(M) && cm.rowsize() == size_t(N),
           "Creating CSmallMatrix(M,N)");

    for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
        cm(i,j) = std::complex<T>(T(k),T(k+1000));

    for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
        Assert(cm(i,j) == std::complex<T>(T(k),T(k+1000)),"Read/Write CSmallMatrix");
        Assert(cm.row(i)(j) == std::complex<T>(T(k),T(k+1000)),"CSmallMatrix.row");
        Assert(cm.col(j)(i) == std::complex<T>(T(k),T(k+1000)),"CSmallMatrix.col");
    }

    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) {
        ca(i,j) = std::complex<T>(3+i+5*j,T(0)+i-j);
        cb(i,j) = std::complex<T>(3+2*i+4*j,4-10*i);
    }

    cm = ca+cb;
    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
        Assert(cm(i,j) == std::complex<T>(6+3*i+9*j,4-9*i-j),"Add CSmallMatrix");

    cm = ca-cb;
    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
        Assert(cm(i,j) == std::complex<T>(T(0)-i+j,-4+11*i-j),"Subtract CSmallMatrix");

    cm = ca;
    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
        Assert(cm(i,j) == ca(i,j),"Copy CSmallMatrix");

    // Test I/O

    std::ofstream fout("tmvtest_smallmatrix_io.dat");
    if (!fout) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_smallmatrix_io.dat for output\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_smallmatrix_io.dat for output");
#endif
    }
    fout << m << std::endl << cm << std::endl;
    fout.close();

    tmv::SmallMatrix<T,M,N,tmv::RowMajor> xm1;
    tmv::SmallMatrix<std::complex<T>,M,N,tmv::RowMajor> xcm1;
    std::ifstream fin("tmvtest_smallmatrix_io.dat");
    if (!fin) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_smallmatrix_io.dat for input\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_smallmatrix_io.dat for input");
#endif
    }
    fin >> xm1 >> xcm1;
    fin.close();
    Assert(m == xm1,"SmallMatrix I/O check #1");
    Assert(cm == xcm1,"CSmallMatrix I/O check #1");

    tmv::SmallMatrix<T,M,N,tmv::ColMajor> xm2;
    tmv::SmallMatrix<std::complex<T>,M,N,tmv::ColMajor> xcm2;
    fin.open("tmvtest_smallmatrix_io.dat");
    if (!fin) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_smallmatrix_io.dat for input\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_smallmatrix_io.dat for input");
#endif
    }
    fin >> xm2 >> xcm2;
    fin.close();
    Assert(m == xm2,"SmallMatrix I/O check #2");
    Assert(cm == xcm2,"CSmallMatrix I/O check #2");

    std::remove("tmvtest_smallmatrix_io.dat");
}

template <class T> 
void TestAllSmallMatrix()
{
    TestBasicSmallMatrix<T,15,10,tmv::RowMajor>();
    TestBasicSmallMatrix<T,15,10,tmv::ColMajor>();
    TestSmallMatrix_Sub<T>();
    std::cout<<"SmallMatrix<"<<tmv::TMV_Text(T())<<"> passed all basic tests\n";
}

#ifdef TEST_DOUBLE
template void TestAllSmallMatrix<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllSmallMatrix<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllSmallMatrix<long double>();
#endif
#ifdef TEST_INT
template void TestAllSmallMatrix<int>();
#endif
