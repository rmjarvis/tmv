#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include <fstream>
#include <cstdio>

template <class T, int M, int N, tmv::StorageType S> 
inline void TestBasicSmallMatrix()
{
    if (showstartdone) {
        std::cout<<"Start TestBasicSmallMatrix\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"M,N = "<<M<<','<<N<<std::endl;
    }

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


    // Test assignments and constructors from arrays:
    T qarrm[] = {
        T(0), T(-1), T(-2), T(-3),
        T(2), T(1), T(0), T(-1),
        T(4), T(3), T(2), T(1)
    };
    T qarcm[] = {
        T(0), T(2), T(4),
        T(-1), T(1), T(3),
        T(-2), T(0), T(2),
        T(-3), T(-1), T(1)
    };
    std::vector<T> qvecrm(12);
    for(int i=0;i<12;i++) qvecrm[i] = qarrm[i];
    std::vector<T> qveccm(12);
    for(int i=0;i<12;i++) qveccm[i] = qarcm[i];

    tmv::SmallMatrix<T,3,4,S> q1;
    std::copy(qarrm, qarrm+12, q1.rowmajor_begin());
    tmv::SmallMatrix<T,3,4,S> q2;
    std::copy(qarcm, qarcm+12, q2.colmajor_begin());

    tmv::SmallMatrix<T,3,4,S> q3;
    std::copy(qvecrm.begin(), qvecrm.end(), q3.rowmajor_begin());
    tmv::SmallMatrix<T,3,4,S> q4;
    std::copy(qveccm.begin(), qveccm.end(), q4.colmajor_begin());

    tmv::SmallMatrix<T,30,40,S> q5x;
    tmv::MatrixView<T> q5 = q5x.subMatrix(3,18,5,25,5,5);
    std::copy(qvecrm.begin(), qvecrm.end(), q5.rowmajor_begin());

    tmv::SmallMatrix<T,30,40,S> q6x;
    tmv::MatrixView<T> q6 = q6x.subMatrix(3,18,5,25,5,5);
    std::copy(qveccm.begin(), qveccm.end(), q6.colmajor_begin());

    // Assignment using op<< is always in rowmajor order.
    tmv::SmallMatrix<T,3,4,S> q7;
    tmv::SmallMatrix<T,4,3,S> q8t;
    typename tmv::SmallMatrix<T,4,3,S>::transpose_type q8 = q8t.transpose();
    q7 <<
        0, -1, -2, -3,
        2, 1, 0, -1,
        4, 3, 2, 1;
    q8 <<
        0, -1, -2, -3,
        2, 1, 0, -1,
        4, 3, 2, 1;

    if (showacc) {
        std::cout<<"q1 = "<<q1<<std::endl;
        std::cout<<"q2 = "<<q2<<std::endl;
        std::cout<<"q3 = "<<q3<<std::endl;
        std::cout<<"q4 = "<<q4<<std::endl;
        std::cout<<"q5 = "<<q5<<std::endl;
        std::cout<<"q6 = "<<q6<<std::endl;
        std::cout<<"q7 = "<<q7<<std::endl;
        std::cout<<"q8 = "<<q8<<std::endl;
    }

    for(int i=0;i<3;i++) for(int j=0;j<4;j++) {
        Assert(q1(i,j) == T(2*i-j),"Create SmallMatrix from T* rm");
        Assert(q2(i,j) == T(2*i-j),"Create SmallMatrix from T* cm");
        Assert(q3(i,j) == T(2*i-j),"Create SmallMatrix from vector rm");
        Assert(q4(i,j) == T(2*i-j),"Create SmallMatrix from vector cm");
        Assert(q5(i,j) == T(2*i-j),"Create SmallMatrixView from vector rm");
        Assert(q6(i,j) == T(2*i-j),"Create SmallMatrixView from vector cm");
        Assert(q7(i,j) == T(2*i-j),"Create SmallMatrix from << list");
        Assert(q8(i,j) == T(2*i-j),"Create SmallMatrixView from << list");
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
