
#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include <fstream>
#include <cstdio>

#define CT std::complex<T>

template <class T, tmv::StorageType S> 
static void TestBasicMatrix_1()
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
    tmv::ConstMatrixView<T> mcv = m.view();
    tmv::MatrixView<T> mv = m.view();
    tmv::ConstMatrixView<T,tmv::FortranStyle> mfcv = mf.view();
    tmv::MatrixView<T,tmv::FortranStyle> mfv = mf.view();

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

template <class T, tmv::StorageType S> 
static void TestBasicMatrix_2()
{
    const int M = 15;
    const int N = 10;

    tmv::Matrix<T,S> m(M,N);
    tmv::Matrix<T,S,tmv::FortranStyle> mf(M,N);

    for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
        m(i,j) = T(k);
        mf(i+1,j+1) = T(k);
    }
    tmv::ConstMatrixView<T> mcv = m.view();
    tmv::MatrixView<T> mv = m.view();
    tmv::ConstMatrixView<T,tmv::FortranStyle> mfcv = mf.view();
    tmv::MatrixView<T,tmv::FortranStyle> mfv = mf.view();

    Assert(m.subMatrix(2,5,1,4) == m.subMatrix(2,5,1,4,1,1),"SubMatrix");
    Assert(m.subVector(2,5,4,2,3) == m.subMatrix(2,14,5,11,4,2).diag(),
           "SubVector");
    Assert(m.colPair(2,5) == m.subMatrix(0,M,2,8,1,3),"colPair");
    Assert(m.colPair(7,2) == m.subMatrix(0,M,7,-3,1,-5),"colPair");
    Assert(m.rowPair(3,7) == m.subMatrix(3,11,0,N,4,1),"rowPair");
    Assert(m.rowPair(2,0) == m.subMatrix(2,-2,0,N,-2,1),"rowPair");
    Assert(m.colRange(2,5) == m.subMatrix(0,M,2,5),"colRange");
    Assert(m.rowRange(3,7) == m.subMatrix(3,7,0,N),"rowRange");

    Assert(mf.subMatrix(3,5,2,4) == mf.subMatrix(3,5,2,4,1,1),"SubMatrixFF");
    Assert(mf.subVector(3,6,4,2,3) == mf.subMatrix(3,11,6,10,4,2).diag(),
           "SubVectorFF");
    Assert(mf.colPair(3,6) == mf.subMatrix(1,M,3,6,1,3),"colPairFF");
    Assert(mf.colPair(8,3) == mf.subMatrix(1,M,8,3,1,-5),"colPairFF");
    Assert(mf.rowPair(4,8) == mf.subMatrix(4,8,1,N,4,1),"rowPairFF");
    Assert(mf.rowPair(3,1) == mf.subMatrix(3,1,1,N,-2,1),"rowPairFF");
    Assert(mf.colRange(3,5) == mf.subMatrix(1,M,3,5),"colRangeFF");
    Assert(mf.rowRange(4,7) == mf.subMatrix(4,7,1,N),"rowRangeFF");

    Assert(m.subMatrix(2,5,1,4) == mf.subMatrix(3,5,2,4),"SubMatrixF");
    Assert(m.subMatrix(2,8,1,10,2,3) == mf.subMatrix(3,7,2,8,2,3),"SubMatrixF");
    Assert(m.subVector(2,5,4,2,3) == mf.subVector(3,6,4,2,3),"SubVectorF");
    Assert(m.subVector(8,1,-1,2,4) == mf.subVector(9,2,-1,2,4),"SubVector2F");
    Assert(m.subVector(12,8,-4,-2,2) == mf.subVector(13,9,-4,-2,2),
           "SubVector3F");
    Assert(m.colPair(2,5) == mf.colPair(3,6),"colPairF");
    Assert(m.colPair(7,2) == mf.colPair(8,3),"colPairF");
    Assert(m.rowPair(3,7) == mf.rowPair(4,8),"rowPairF");
    Assert(m.rowPair(2,0) == mf.rowPair(3,1),"rowPairF");
    Assert(m.colRange(2,5) == mf.colRange(3,5),"colRangeF");
    Assert(m.rowRange(3,7) == mf.rowRange(4,7),"rowRangeF");

    Assert(m.subMatrix(2,5,1,4) == mcv.subMatrix(2,5,1,4),"SubMatrixCV");
    Assert(m.subMatrix(2,8,1,10,2,3) == mcv.subMatrix(2,8,1,10,2,3),
           "SubMatrixCV");
    Assert(m.subVector(2,5,4,2,3) == mcv.subVector(2,5,4,2,3),"SubVectorCV");
    Assert(m.subVector(8,1,-1,2,4) == mcv.subVector(8,1,-1,2,4),"SubVector2CV");
    Assert(m.subVector(12,8,-4,-2,2) == mcv.subVector(12,8,-4,-2,2),
           "SubVector3CV");
    Assert(m.colPair(2,5) == mcv.colPair(2,5),"colPairCV");
    Assert(m.colPair(7,2) == mcv.colPair(7,2),"colPairCV");
    Assert(m.rowPair(3,7) == mcv.rowPair(3,7),"rowPairCV");
    Assert(m.rowPair(2,0) == mcv.rowPair(2,0),"rowPairCV");
    Assert(m.colRange(2,5) == mcv.colRange(2,5),"colRangeCV");
    Assert(m.rowRange(3,7) == mcv.rowRange(3,7),"rowRangeCV");

    Assert(m.subMatrix(2,5,1,4) == mv.subMatrix(2,5,1,4),"SubMatrixV");
    Assert(m.subMatrix(2,8,1,10,2,3) == mv.subMatrix(2,8,1,10,2,3),"SubMatrixV");
    Assert(m.subVector(2,5,4,2,3) == mv.subVector(2,5,4,2,3),"SubVectorV");
    Assert(m.subVector(8,1,-1,2,4) == mv.subVector(8,1,-1,2,4),"SubVector2V");
    Assert(m.subVector(12,8,-4,-2,2) == mv.subVector(12,8,-4,-2,2),
           "SubVector3V");
    Assert(m.colPair(2,5) == mv.colPair(2,5),"colPairV");
    Assert(m.colPair(7,2) == mv.colPair(7,2),"colPairV");
    Assert(m.rowPair(3,7) == mv.rowPair(3,7),"rowPairV");
    Assert(m.rowPair(2,0) == mv.rowPair(2,0),"rowPairV");
    Assert(m.colRange(2,5) == mv.colRange(2,5),"colRangeV");
    Assert(m.rowRange(3,7) == mv.rowRange(3,7),"rowRangeV");

    Assert(mf.subMatrix(3,5,2,4) == mfcv.subMatrix(3,5,2,4),"SubMatrixFCV");
    Assert(mf.subMatrix(3,7,2,8,2,3) == mfcv.subMatrix(3,7,2,8,2,3),
           "SubMatrixFCV");
    Assert(mf.subVector(3,6,4,2,3) == mfcv.subVector(3,6,4,2,3),"SubVectorFCV");
    Assert(mf.subVector(9,2,-1,2,4) == mfcv.subVector(9,2,-1,2,4),
           "SubVector2FCV");
    Assert(mf.subVector(13,9,-4,-2,2) == mfcv.subVector(13,9,-4,-2,2),
           "SubVector3FCV");
    Assert(mf.colPair(3,6) == mfcv.colPair(3,6),"colPairFCV");
    Assert(mf.colPair(8,3) == mfcv.colPair(8,3),"colPairFCV");
    Assert(mf.rowPair(4,8) == mfcv.rowPair(4,8),"rowPairFCV");
    Assert(mf.rowPair(3,1) == mfcv.rowPair(3,1),"rowPairFCV");
    Assert(mf.colRange(3,5) == mfcv.colRange(3,5),"colRangeFCV");
    Assert(mf.rowRange(4,7) == mfcv.rowRange(4,7),"rowRangeFCV");

    Assert(mf.subMatrix(3,5,2,4) == mfv.subMatrix(3,5,2,4),"SubMatrixFV");
    Assert(mf.subMatrix(3,7,2,8,2,3) == mfv.subMatrix(3,7,2,8,2,3),
           "SubMatrixFV");
    Assert(mf.subVector(3,6,4,2,3) == mfv.subVector(3,6,4,2,3),"SubVectorFV");
    Assert(mf.subVector(9,2,-1,2,4) == mfv.subVector(9,2,-1,2,4),
           "SubVector2FV");
    Assert(mf.subVector(13,9,-4,-2,2) == mfv.subVector(13,9,-4,-2,2),
           "SubVector3FV");
    Assert(mf.colPair(3,6) == mfv.colPair(3,6),"colPairFV");
    Assert(mf.colPair(8,3) == mfv.colPair(8,3),"colPairFV");
    Assert(mf.rowPair(4,8) == mfv.rowPair(4,8),"rowPairFV");
    Assert(mf.rowPair(3,1) == mfv.rowPair(3,1),"rowPairFV");
    Assert(mf.colRange(3,5) == mfv.colRange(3,5),"colRangeFV");
    Assert(mf.rowRange(4,7) == mfv.rowRange(4,7),"rowRangeFV");

    tmv::Matrix<T,S> a(M,N);
    tmv::Matrix<T,S> b(M,N);
    tmv::Matrix<T,S> c(M,N);
    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) {
        a(i,j) = T(3+i+5*j);
        b(i,j) = T(5+2*i+4*j);
    }
    mf = a;
    Assert(a == mf,"Copy CStyle Matrix to FortranStyle");

    std::vector<T> qv(6);
    tmv::Matrix<T,S> q4(2,3);
    tmv::Matrix<T,S> q5t(3,2);
    tmv::MatrixView<T> q5 = q5t.transpose();
    if (S == tmv::RowMajor) {
        T qvar[] = { T(0), T(-1), T(-2),
            T(2), T(1), T(0) };
        for(int i=0;i<6;i++) qv[i] = qvar[i];
        q4 <<
            0, -1, -2,
            2, 1, 0;
        q5 <<
            0, 2,
            -1, 1,
            -2, 0;
    } else {
        T qvar[] = { T(0), T(2),
            T(-1), T(1),
            T(-2), T(0) };
        for(int i=0;i<6;i++) qv[i] = qvar[i];
        q4 <<
            0, 2,
            -1, 1,
            -2, 0;
        q5 <<
            0, -1, -2,
            2, 1, 0;
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
        Assert(q4(i,j) == T(2*i-j),"Create Matrix from <<");
        Assert(q5(i,j) == T(2*i-j),"Create MatrixView of <<");
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

template <class T, tmv::StorageType S> 
static void TestBasicMatrix_IO()
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
    if (!fout) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_matrix_io.dat for output\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_matrix_io.dat for output");
#endif
    }
    fout << m << std::endl << cm << std::endl;
    fout.close();

    tmv::Matrix<T,tmv::RowMajor> xm1(M,N);
    tmv::Matrix<CT,tmv::RowMajor> xcm1(M,N);
    std::ifstream fin("tmvtest_matrix_io.dat");
    if (!fin) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_matrix_io.dat for input\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_matrix_io.dat for input");
#endif
    }
    fin >> xm1 >> xcm1;
    fin.close();
    Assert(m == xm1,"Matrix I/O check #1");
    Assert(cm == xcm1,"CMatrix I/O check #1");

    tmv::Matrix<T,tmv::ColMajor> xm2(M,N);
    tmv::Matrix<CT,tmv::ColMajor> xcm2(M,N);
    fin.open("tmvtest_matrix_io.dat");
    if (!fin) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_matrix_io.dat for input\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_matrix_io.dat for input");
#endif
    }
    fin >> xm2 >> xcm2;
    fin.close();
    Assert(m == xm2,"Matrix I/O check #2");
    Assert(cm == xcm2,"CMatrix I/O check #2");

    std::auto_ptr<tmv::Matrix<T> > xm3;
    std::auto_ptr<tmv::Matrix<CT> > xcm3;
    fin.open("tmvtest_matrix_io.dat");
    if (!fin) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_matrix_io.dat for input\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_matrix_io.dat for input");
#endif
    }
    fin >> xm3 >> xcm3;
    fin.close();
    Assert(m == *xm3,"Matrix I/O check #3");
    Assert(cm == *xcm3,"CMatrix I/O check #3");

#ifndef XTEST
    std::remove("tmvtest_matrix_io.dat");
    //system("rm tmvtest_matrix_io.dat");
#endif

}

template <class T> 
void TestAllMatrix()
{
    TestBasicMatrix_1<T,tmv::RowMajor>();
    TestBasicMatrix_1<T,tmv::ColMajor>();
    TestBasicMatrix_2<T,tmv::RowMajor>();
    TestBasicMatrix_2<T,tmv::ColMajor>();
    TestBasicMatrix_IO<T,tmv::RowMajor>();
    TestBasicMatrix_IO<T,tmv::ColMajor>();
    std::cout<<"Matrix<"<<tmv::TMV_Text(T())<<"> passed all basic tests\n";

    if (tmv::TMV_Epsilon<T>() > T(0)) {
        TestAllMatrixArith<T>();
        std::cout<<"Matrix<"<<tmv::TMV_Text(T())<<"> Arithmetic passed all tests\n";
    }
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
