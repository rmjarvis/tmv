
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"
#include <fstream>
#include <cstdio>

#define CT std::complex<T>

template <class T, tmv::StorageType S> 
static void TestBasicMatrix_1()
{
    const int M = 15;
    const int N = 10;

    if (showstartdone) {
        std::cout<<"Start TestBasicMatrix_1\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"M,N = "<<M<<','<<N<<std::endl;
    }

    tmv::Matrix<T,S> m(M,N);
    tmv::Matrix<T,S|tmv::FortranStyle> mf(M,N);
    Assert(m.colsize() == M && m.rowsize() == N,
           "Creating Matrix(M,N)");
    Assert(m.colsize() == M && m.rowsize() == N,
           "Creating MatrixF(M,N)");

    for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
        m(i,j) = T(k);
        mf(i+1,j+1) = T(k);
    }
    tmv::ConstMatrixView<T> mcv = m.view();
    tmv::MatrixView<T> mv = m.view();
    tmv::ConstMatrixView<T,tmv::FortranStyle> mfcv = mf.view();
    tmv::MatrixView<T,tmv::FortranStyle> mfv = mf.view();

    tmv::VIt<T,1,tmv::NonConj> vit = m.linearView().begin();
    tmv::CVIt<T,1,tmv::NonConj> cvit = m.linearView().begin();
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
        if (S == tmv::RowMajor) {
            Assert(*vit++ == k,"Matrix.linearView");
            Assert(*cvit++ == k,"Matrix.linearView with CVIt");
        }
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

    m.resize(2,3);
    Assert(m.colsize() == 2 && m.rowsize() == 3,"m.resize(2,3)");
    for (int i=0, k=0; i<2; ++i) for (int j=0; j<3; ++j, ++k) 
        m(i,j) = T(k);
    for (int i=0, k=0; i<2; ++i) for (int j=0; j<3; ++j, ++k) 
        Assert(m(i,j) == k,"Read/Write resized Matrix");

    m.resize(2*M,3*N);
    Assert(m.colsize() == 2*M && m.rowsize() == 3*N,"m.resize(2*M,3*N)");
    for (int i=0, k=0; i<2*M; ++i) for (int j=0; j<3*N; ++j, ++k) 
        m(i,j) = T(k);
    for (int i=0, k=0; i<2*M; ++i) for (int j=0; j<3*N; ++j, ++k) 
        Assert(m(i,j) == k,"Read/Write resized Matrix");

    // Test Basic Arithmetic
    tmv::Matrix<T,S> a(M,N);
    tmv::Matrix<T,S> b(M,N);
    tmv::Matrix<T,S> c(M,N);
    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) {
        a(i,j) = T(3+i+5*j);
        b(i,j) = T(5+2*i+4*j);
    }
    mf = a;
    Assert(a == mf,"Copy CStyle Matrix to FortranStyle");

    c = a+b;
    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
        Assert(c(i,j) == T(8+3*i+9*j),"Add Matrices");

    c = a-b;
    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
        Assert(c(i,j) == T(-2-i+j),"Subtract Matrices");

    tmv::Matrix<CT,S> cm(M,N);
    tmv::Matrix<CT,S> ca(M,N);
    tmv::Matrix<CT,S> cb(M,N);
    Assert(cm.colsize() == M && cm.rowsize() && N,
           "Creating CMatrix(M,N)");

    for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
        cm(i,j) = CT(T(k),T(k+1000));

    tmv::VIt<CT,1,tmv::NonConj> vit2 = cm.linearView().begin();
    tmv::CVIt<CT,1,tmv::NonConj> cvit2 = cm.linearView().begin();
    tmv::VIt<CT,1,tmv::Conj> vit3 = cm.conjugate().linearView().begin();
    tmv::CVIt<CT,1,tmv::Conj> cvit3 = cm.conjugate().linearView().begin();
    for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
        Assert(cm(i,j) == CT(T(k),T(k+1000)),"Read/Write CMatrix");
        Assert(cm.row(i)(j) == CT(T(k),T(k+1000)),"CMatrix.row");
        Assert(cm.col(j)(i) == CT(T(k),T(k+1000)),"CMatrix.col");
        if (S == tmv::RowMajor) {
            Assert(*vit2++ == CT(T(k),T(k+1000)),"CMatrix.linearView");
            Assert(*cvit2++ == CT(T(k),T(k+1000)),"CMatrix.linearView with CVIt");
            Assert(*vit3++ == CT(T(k),-T(k+1000)),"CMatrix.linearView with conj");
            Assert(*cvit3++ == CT(T(k),-T(k+1000)),"CMatrix.linearView with conj, CVIt");
        }
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
static void TestBasicMatrix_2()
{
    const int M = 15;
    const int N = 10;

    if (showstartdone) {
        std::cout<<"Start TestBasicMatrix_2\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"M,N = "<<M<<','<<N<<std::endl;
    }

    tmv::Matrix<T,S> m(M,N);
    tmv::Matrix<T,S|tmv::FortranStyle> mf(M,N);

    for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
        m(i,j) = T(k);
        mf(i+1,j+1) = T(k);
    }
    tmv::ConstMatrixView<T> mcv = m.view();
    tmv::MatrixView<T> mv = m.view();
    tmv::ConstMatrixView<T,tmv::FortranStyle> mfcv = mf.view();
    tmv::MatrixView<T,tmv::FortranStyle> mfv = mf.view();

    Assert(m.subMatrix(2,5,1,4) == m.subMatrix(2,5,1,4,1,1),"subMatrix");
    Assert(m.subVector(2,5,4,2,3) == m.subMatrix(2,14,5,11,4,2).diag(),
           "subVector");
    Assert(m.colPair(2,5) == m.subMatrix(0,M,2,8,1,3),"colPair");
    Assert(m.colPair(7,2) == m.subMatrix(0,M,7,-3,1,-5),"colPair");
    Assert(m.rowPair(3,7) == m.subMatrix(3,11,0,N,4,1),"rowPair");
    Assert(m.rowPair(2,0) == m.subMatrix(2,-2,0,N,-2,1),"rowPair");
    Assert(m.colRange(2,5) == m.subMatrix(0,M,2,5),"colRange");
    Assert(m.rowRange(3,7) == m.subMatrix(3,7,0,N),"rowRange");

    Assert(mf.subMatrix(3,5,2,4) == mf.subMatrix(3,5,2,4,1,1),"subMatrixFF");
    Assert(mf.subVector(3,6,4,2,3) == mf.subMatrix(3,11,6,10,4,2).diag(),
           "subVectorFF");
    Assert(mf.colPair(3,6) == mf.subMatrix(1,M,3,6,1,3),"colPairFF");
    Assert(mf.colPair(8,3) == mf.subMatrix(1,M,8,3,1,-5),"colPairFF");
    Assert(mf.rowPair(4,8) == mf.subMatrix(4,8,1,N,4,1),"rowPairFF");
    Assert(mf.rowPair(3,1) == mf.subMatrix(3,1,1,N,-2,1),"rowPairFF");
    Assert(mf.colRange(3,5) == mf.subMatrix(1,M,3,5),"colRangeFF");
    Assert(mf.rowRange(4,7) == mf.subMatrix(4,7,1,N),"rowRangeFF");

    Assert(m.subMatrix(2,5,1,4) == mf.subMatrix(3,5,2,4),"subMatrixF");
    Assert(m.subMatrix(2,8,1,10,2,3) == mf.subMatrix(3,7,2,8,2,3),"subMatrixF");
    Assert(m.subVector(2,5,4,2,3) == mf.subVector(3,6,4,2,3),"subVectorF");
    Assert(m.subVector(8,1,-1,2,4) == mf.subVector(9,2,-1,2,4),"subVector2F");
    Assert(m.subVector(12,8,-4,-2,2) == mf.subVector(13,9,-4,-2,2),
           "subVector3F");
    Assert(m.colPair(2,5) == mf.colPair(3,6),"colPairF");
    Assert(m.colPair(7,2) == mf.colPair(8,3),"colPairF");
    Assert(m.rowPair(3,7) == mf.rowPair(4,8),"rowPairF");
    Assert(m.rowPair(2,0) == mf.rowPair(3,1),"rowPairF");
    Assert(m.colRange(2,5) == mf.colRange(3,5),"colRangeF");
    Assert(m.rowRange(3,7) == mf.rowRange(4,7),"rowRangeF");

    Assert(m.subMatrix(2,5,1,4) == mcv.subMatrix(2,5,1,4),"subMatrixCV");
    Assert(m.subMatrix(2,8,1,10,2,3) == mcv.subMatrix(2,8,1,10,2,3),
           "subMatrixCV");
    Assert(m.subVector(2,5,4,2,3) == mcv.subVector(2,5,4,2,3),"subVectorCV");
    Assert(m.subVector(8,1,-1,2,4) == mcv.subVector(8,1,-1,2,4),"subVector2CV");
    Assert(m.subVector(12,8,-4,-2,2) == mcv.subVector(12,8,-4,-2,2),
           "subVector3CV");
    Assert(m.colPair(2,5) == mcv.colPair(2,5),"colPairCV");
    Assert(m.colPair(7,2) == mcv.colPair(7,2),"colPairCV");
    Assert(m.rowPair(3,7) == mcv.rowPair(3,7),"rowPairCV");
    Assert(m.rowPair(2,0) == mcv.rowPair(2,0),"rowPairCV");
    Assert(m.colRange(2,5) == mcv.colRange(2,5),"colRangeCV");
    Assert(m.rowRange(3,7) == mcv.rowRange(3,7),"rowRangeCV");

    Assert(m.subMatrix(2,5,1,4) == mv.subMatrix(2,5,1,4),"subMatrixV");
    Assert(m.subMatrix(2,8,1,10,2,3) == mv.subMatrix(2,8,1,10,2,3),"subMatrixV");
    Assert(m.subVector(2,5,4,2,3) == mv.subVector(2,5,4,2,3),"subVectorV");
    Assert(m.subVector(8,1,-1,2,4) == mv.subVector(8,1,-1,2,4),"subVector2V");
    Assert(m.subVector(12,8,-4,-2,2) == mv.subVector(12,8,-4,-2,2),
           "subVector3V");
    Assert(m.colPair(2,5) == mv.colPair(2,5),"colPairV");
    Assert(m.colPair(7,2) == mv.colPair(7,2),"colPairV");
    Assert(m.rowPair(3,7) == mv.rowPair(3,7),"rowPairV");
    Assert(m.rowPair(2,0) == mv.rowPair(2,0),"rowPairV");
    Assert(m.colRange(2,5) == mv.colRange(2,5),"colRangeV");
    Assert(m.rowRange(3,7) == mv.rowRange(3,7),"rowRangeV");

    Assert(mf.subMatrix(3,5,2,4) == mfcv.subMatrix(3,5,2,4),"subMatrixFCV");
    Assert(mf.subMatrix(3,7,2,8,2,3) == mfcv.subMatrix(3,7,2,8,2,3),
           "subMatrixFCV");
    Assert(mf.subVector(3,6,4,2,3) == mfcv.subVector(3,6,4,2,3),"subVectorFCV");
    Assert(mf.subVector(9,2,-1,2,4) == mfcv.subVector(9,2,-1,2,4),
           "subVector2FCV");
    Assert(mf.subVector(13,9,-4,-2,2) == mfcv.subVector(13,9,-4,-2,2),
           "subVector3FCV");
    Assert(mf.colPair(3,6) == mfcv.colPair(3,6),"colPairFCV");
    Assert(mf.colPair(8,3) == mfcv.colPair(8,3),"colPairFCV");
    Assert(mf.rowPair(4,8) == mfcv.rowPair(4,8),"rowPairFCV");
    Assert(mf.rowPair(3,1) == mfcv.rowPair(3,1),"rowPairFCV");
    Assert(mf.colRange(3,5) == mfcv.colRange(3,5),"colRangeFCV");
    Assert(mf.rowRange(4,7) == mfcv.rowRange(4,7),"rowRangeFCV");

    Assert(mf.subMatrix(3,5,2,4) == mfv.subMatrix(3,5,2,4),"subMatrixFV");
    Assert(mf.subMatrix(3,7,2,8,2,3) == mfv.subMatrix(3,7,2,8,2,3),"subMatrixFV");
    Assert(mf.subVector(3,6,4,2,3) == mfv.subVector(3,6,4,2,3),"subVectorFV");
    Assert(mf.subVector(9,2,-1,2,4) == mfv.subVector(9,2,-1,2,4),"subVector2FV");
    Assert(mf.subVector(13,9,-4,-2,2) == mfv.subVector(13,9,-4,-2,2),
           "subVector3FV");
    Assert(mf.colPair(3,6) == mfv.colPair(3,6),"colPairFV");
    Assert(mf.colPair(8,3) == mfv.colPair(8,3),"colPairFV");
    Assert(mf.rowPair(4,8) == mfv.rowPair(4,8),"rowPairFV");
    Assert(mf.rowPair(3,1) == mfv.rowPair(3,1),"rowPairFV");
    Assert(mf.colRange(3,5) == mfv.colRange(3,5),"colRangeFV");
    Assert(mf.rowRange(4,7) == mfv.rowRange(4,7),"rowRangeFV");


    // Test assignments and constructors from arrays
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

    tmv::Matrix<T,S> q1(3,4);
    std::copy(qarrm, qarrm+12, q1.rowmajor_begin());
    tmv::Matrix<T,S> q2(3,4);
    std::copy(qarcm, qarcm+12, q2.colmajor_begin());

    tmv::Matrix<T,S> q3(3,4);
    std::copy(qvecrm.begin(), qvecrm.end(), q3.rowmajor_begin());
    tmv::Matrix<T,S> q4(3,4);
    std::copy(qveccm.begin(), qveccm.end(), q4.colmajor_begin());

    tmv::Matrix<T,S> q5x(30,40);
    tmv::MatrixView<T> q5 = q5x.subMatrix(3,18,5,25,5,5);
    std::copy(qvecrm.begin(), qvecrm.end(), q5.rowmajor_begin());

    tmv::Matrix<T,S> q6x(30,40);
    tmv::MatrixView<T> q6 = q6x.subMatrix(3,18,5,25,5,5);
    std::copy(qveccm.begin(), qveccm.end(), q6.colmajor_begin());

    // Assignment using op<< is always in rowmajor order.
    tmv::Matrix<T,S> q7(3,4);
    tmv::Matrix<T,S> q8t(4,3);
    tmv::MatrixView<T> q8 = q8t.transpose();
    q7 <<
        0, -1, -2, -3,
        2, 1, 0, -1,
        4, 3, 2, 1;
    q8 <<
        0, -1, -2, -3,
        2, 1, 0, -1,
        4, 3, 2, 1;

    // Can also view memory directly
    T* qarS = (S == tmv::RowMajor) ? qarrm : qarcm;
    tmv::ConstMatrixView<T> q9 = tmv::MatrixViewOf(qarS,3,4,S);
    const int Si = (S == tmv::RowMajor ? 4 : 1);
    const int Sj = (S == tmv::RowMajor ? 1 : 3);
    tmv::ConstMatrixView<T> q10 = tmv::MatrixViewOf(qarS,3,4,Si,Sj);

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
    }

    for(int i=0;i<3;i++) for(int j=0;j<4;j++) {
        Assert(q1(i,j) == T(2*i-j),"Create Matrix from T* rm");
        Assert(q2(i,j) == T(2*i-j),"Create Matrix from T* cm");
        Assert(q3(i,j) == T(2*i-j),"Create Matrix from vector rm");
        Assert(q4(i,j) == T(2*i-j),"Create Matrix from vector cm");
        Assert(q5(i,j) == T(2*i-j),"Create MatrixView from vector rm");
        Assert(q6(i,j) == T(2*i-j),"Create MatrixView from vector cm");
        Assert(q7(i,j) == T(2*i-j),"Create Matrix from << list");
        Assert(q8(i,j) == T(2*i-j),"Create MatrixView from << list");
        Assert(q9(i,j) == T(2*i-j),"Create MatrixView of T* (S)");
        Assert(q10(i,j) == T(2*i-j),"Create MatrixView of T* (Si,Sj)");
    }

    // Test the span of the iteration (i.e. the validity of begin(), end())
    const tmv::Matrix<T,S>& q1_const = q1;
    tmv::MatrixView<T> q1_view = q1.view();
    tmv::ConstMatrixView<T> q1_constview = q1_const.view();
    tmv::ConstMatrixView<T> q5_const = q5;

    typename tmv::Matrix<T,S>::rowmajor_iterator rmit1 = q1.rowmajor_begin();
    typename tmv::Matrix<T,S>::const_rowmajor_iterator rmit2 = 
        q1_const.rowmajor_begin();
    typename tmv::MatrixView<T>::rowmajor_iterator rmit3 = 
        q1_view.rowmajor_begin();
    typename tmv::ConstMatrixView<T>::const_rowmajor_iterator rmit4 = 
        q1_constview.rowmajor_begin();
    typename tmv::MatrixView<T>::rowmajor_iterator rmit5 = 
        q5.rowmajor_begin();
    typename tmv::ConstMatrixView<T>::const_rowmajor_iterator rmit6 = 
        q5_const.rowmajor_begin();
    int i = 0;
    while (rmit1 != q1.rowmajor_end()) {
        Assert(*rmit1++ == qarrm[i], "RowMajor iteration 1");
        Assert(*rmit2++ == qarrm[i], "RowMajor iteration 2");
        Assert(*rmit3++ == qarrm[i], "RowMajor iteration 3");
        Assert(*rmit4++ == qarrm[i], "RowMajor iteration 4");
        Assert(*rmit5++ == qarrm[i], "RowMajor iteration 5");
        Assert(*rmit6++ == qarrm[i], "RowMajor iteration 6");
        ++i;
    }
    Assert(i == 12, "RowMajor iteration number of elements");
    Assert(rmit2 == q1_const.rowmajor_end(), "rmit2 reaching end");
    Assert(rmit3 == q1_view.rowmajor_end(), "rmit3 reaching end");
    Assert(rmit4 == q1_constview.rowmajor_end(), "rmit4 reaching end");
    Assert(rmit5 == q5.rowmajor_end(), "rmit5 reaching end");
    Assert(rmit6 == q5_const.rowmajor_end(), "rmit6 reaching end");

    typename tmv::Matrix<T,S>::colmajor_iterator cmit1 = q1.colmajor_begin();
    typename tmv::Matrix<T,S>::const_colmajor_iterator cmit2 = 
        q1_const.colmajor_begin();
    typename tmv::MatrixView<T>::colmajor_iterator cmit3 = 
        q1_view.colmajor_begin();
    typename tmv::ConstMatrixView<T>::const_colmajor_iterator cmit4 = 
        q1_constview.colmajor_begin();
    typename tmv::MatrixView<T>::colmajor_iterator cmit5 = 
        q5.colmajor_begin();
    typename tmv::ConstMatrixView<T>::const_colmajor_iterator cmit6 = 
        q5_const.colmajor_begin();
    i = 0;
    while (cmit1 != q1.colmajor_end()) {
        Assert(*cmit1++ == qarcm[i], "ColMajor iteration 1");
        Assert(*cmit2++ == qarcm[i], "ColMajor iteration 2");
        Assert(*cmit3++ == qarcm[i], "ColMajor iteration 3");
        Assert(*cmit4++ == qarcm[i], "ColMajor iteration 4");
        Assert(*cmit5++ == qarcm[i], "ColMajor iteration 5");
        Assert(*cmit6++ == qarcm[i], "ColMajor iteration 6");
        ++i;
    }
    Assert(i == 12, "ColMajor iteration number of elements");
    Assert(cmit2 == q1_const.colmajor_end(), "cmit2 reaching end");
    Assert(cmit3 == q1_view.colmajor_end(), "cmit3 reaching end");
    Assert(cmit4 == q1_constview.colmajor_end(), "cmit4 reaching end");
    Assert(cmit5 == q5.colmajor_end(), "cmit5 reaching end");
    Assert(cmit6 == q5_const.colmajor_end(), "cmit6 reaching end");

}

template <class T, tmv::StorageType S> 
static void TestBasicMatrix_IO()
{
    const int M = 15;
    const int N = 10;

    if (showstartdone) {
        std::cout<<"Start TestBasicMatrix_IO\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"M,N = "<<M<<','<<N<<std::endl;
    }

    tmv::Matrix<T,S> m(M,N);
    tmv::Matrix<CT,S> cm(M,N);

    for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
        m(i,j) = T(k);
        cm(i,j) = CT(k,k+1000);
    }
    m(3,1) = T(1.e-30);
    cm(3,1) = CT(T(1.e-30),T(1.e-30));
    m(5,6) = T(9.e-3);
    cm(5,6) = CT(T(9.e-3),T(9.e-3));
    cm(6,6) = CT(T(9),T(9.e-3));
    m(7,4) = T(0.123456789);
    cm(7,4) = CT(T(3.123456789),T(6.987654321));

    // First check clipping function...
    tmv::Matrix<T> m2 = m;
    tmv::Matrix<CT> cm2 = cm;
    if (!std::numeric_limits<T>::is_integer) {
        m2.clip(T(1.e-2));
        cm2.clip(T(1.e-2));
    }
    tmv::Matrix<T> m3 = m;
    tmv::Matrix<CT> cm3 = cm;
    m3(3,1) = T(0);
    cm3(3,1) = T(0);
    m3(5,6) = T(0); 
    // Others, esp. cm3(5,6) and cm3(6,6), shouldn't get clipped.
    Assert(m2 == m3,"Matrix clip");
    Assert(cm2 == cm3,"Complex Matrix clip");

    // However, ThreshIO for complex works slightly differently than clip.
    // It clips _either_ the real or imag component, so now cm2(5,6) and
    // cm2(6,6) need to be modified.
    cm2(5,6) = cm3(5,6) = T(0);
    cm2(6,6) = cm3(6,6) = T(9);

    // Write matrices with 4 different styles
    std::ofstream fout("tmvtest_matrix_io.dat");
    Assert(bool(fout),"Couldn't open tmvtest_matrix_io.dat for output");
    fout << m << std::endl;
    fout << cm << std::endl;
    fout << tmv::CompactIO() << m << std::endl;
    fout << tmv::CompactIO() << cm << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << m << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << cm << std::endl;
    // Not a very pretty IO style, but it tests being able to read
    // a style that has no whitespace and has more than one 
    // character for some of the markup elements.
    tmv::IOStyle myStyle = 
        tmv::CompactIO().setThresh(1.e-2).setPrecision(4).
        markup("Start","[",",","]","---","Done");
    fout << myStyle << m << std::endl;
    fout << myStyle << cm << std::endl;
    fout.close();

    // When using (the default) prec(6), these will be the values read in.
    m(7,4) = T(0.123457);
    cm(7,4) = CT(T(3.12346),T(6.98765));

    // When using prec(12), the full correct values will be read in. (m2,cm2)

    // When using prec(4), these will be the values read in.
    m3(7,4) = T(0.1235);
    if (std::numeric_limits<T>::is_integer) cm3(7,4) = CT(3,6);
    else cm3(7,4) = CT(T(3.123),T(6.988));

    // Read them back in
    tmv::Matrix<T,tmv::RowMajor> xm1(M,N);
    tmv::Matrix<CT,tmv::RowMajor> xcm1(M,N);
    std::ifstream fin("tmvtest_matrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_matrix_io.dat for input");
    fin >> xm1 >> xcm1;
    Assert(EqualIO(m,xm1,EPS),"Matrix I/O check normal");
    Assert(EqualIO(cm,xcm1,EPS),"CMatrix I/O check normal");
    fin >> tmv::CompactIO() >> xm1 >> tmv::CompactIO() >> xcm1;
    Assert(EqualIO(m,xm1,EPS),"Matrix I/O check compact");
    Assert(EqualIO(cm,xcm1,EPS),"CMatrix I/O check compact");
    fin >> xm1.view() >> xcm1.view();
    Assert(EqualIO(m2,xm1,EPS),"Matrix I/O check thresh");
    Assert(EqualIO(cm2,xcm1,EPS),"CMatrix I/O check thresh");
    fin >> myStyle >> xm1.view() >> myStyle >> xcm1.view();
    Assert(EqualIO(m3,xm1,EPS),"Matrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(cm3,xcm1,EPS),"CMatrix I/O check compact thresh & prec(4)");
    fin.close();

    // Repeat for column major
    tmv::Matrix<T,tmv::ColMajor> xm2(M,N);
    tmv::Matrix<CT,tmv::ColMajor> xcm2(M,N);
    fin.open("tmvtest_matrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_matrix_io.dat for input");
    fin >> xm2.view() >> xcm2.view();
    Assert(EqualIO(m,xm2,EPS),"Matrix I/O check normal");
    Assert(EqualIO(cm,xcm2,EPS),"CMatrix I/O check normal");
    fin >> tmv::CompactIO() >> xm2.view() >> tmv::CompactIO() >> xcm2.view();
    Assert(EqualIO(m,xm2,EPS),"Matrix I/O check compact");
    Assert(EqualIO(cm,xcm2,EPS),"CMatrix I/O check compact");
    fin >> xm2 >> xcm2;
    Assert(EqualIO(m2,xm2,EPS),"Matrix I/O check thresh");
    Assert(EqualIO(cm2,xcm2,EPS),"CMatrix I/O check thresh");
    fin >> myStyle >> xm2 >> myStyle >> xcm2;
    Assert(EqualIO(m3,xm2,EPS),"Matrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(cm3,xcm2,EPS),"CMatrix I/O check compact thresh & prec(4)");
    fin.close();

    // And repeat for matrices that need to be resized.
    // Also check switching the default IOStyle.
    tmv::CompactIO().makeDefault();
    tmv::Matrix<T> zm1,zm2,zm3,zm4;
    tmv::Matrix<CT> zcm1,zcm2,zcm3,zcm4;
    fin.open("tmvtest_matrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_matrix_io.dat for input");
    fin >> tmv::NormalIO() >> zm1 >> tmv::NormalIO() >> zcm1;
    Assert(EqualIO(m,zm1,EPS),"Matrix I/O check normal with resize");
    Assert(EqualIO(cm,zcm1,EPS),"CMatrix I/O check normal with resize");
    fin >> zm2 >> zcm2;
    Assert(EqualIO(m,zm2,EPS),"Matrix I/O check compact with resize");
    Assert(EqualIO(cm,zcm2,EPS),"CMatrix I/O check compact with resize");
    fin >> tmv::NormalIO() >> zm3 >> tmv::NormalIO() >> zcm3;
    Assert(EqualIO(m2,zm3,EPS),"Matrix I/O check thresh with resize");
    Assert(EqualIO(cm2,zcm3,EPS),"CMatrix I/O check thresh with resize");
    fin >> myStyle >> zm4 >> myStyle >> zcm4;
    Assert(EqualIO(m3,zm4,EPS),"Matrix I/O check compact thresh with resize");
    Assert(EqualIO(cm3,zcm4,EPS),"CMatrix I/O check compact thresh with resize");
    fin.close();
    // Switch it back.
    tmv::IOStyle::revertDefault();

#if XTEST == 0
    std::remove("tmvtest_matrix_io.dat");
#endif

}

template <class T> void TestMatrix()
{
#if 1
    TestBasicMatrix_1<T,tmv::RowMajor>();
    TestBasicMatrix_1<T,tmv::ColMajor>();
    TestBasicMatrix_2<T,tmv::RowMajor>();
    TestBasicMatrix_2<T,tmv::ColMajor>();
    TestBasicMatrix_IO<T,tmv::RowMajor>();
    TestBasicMatrix_IO<T,tmv::ColMajor>();
    std::cout<<"Matrix<"<<tmv::TMV_Text(T())<<"> passed all basic tests\n";
#endif

#if 1
    TestMatrixArith_1<T>();
    TestMatrixArith_2<T>();
    TestMatrixArith_3<T>();
    TestMatrixArith_4<T>();
    TestMatrixArith_5<T>();
    TestMatrixArith_6<T>();
    TestMatrixArith_7<T>();
    TestMatrixArith_8<T>();
    std::cout<<"Matrix<"<<tmv::TMV_Text(T())<<"> Arithmetic passed all tests\n";
#endif
}

#ifdef TEST_DOUBLE
template void TestMatrix<double>();
#endif
#ifdef TEST_FLOAT
template void TestMatrix<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestMatrix<long double>();
#endif
#ifdef TEST_INT
template void TestMatrix<int>();
#endif
