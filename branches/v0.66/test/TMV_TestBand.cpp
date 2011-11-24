
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include <fstream>
#include <cstdio>

template <class T, tmv::StorageType S> 
static void TestBasicBandMatrix_1()
{
    const int N = 10;
    const int nhi = 1;
    const int nlo = 3;

    if (showstartdone) {
        std::cout<<"Start TestBasicBandMatrix_1\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
        std::cout<<"nlo, nhi = "<<nlo<<','<<nhi<<std::endl;
    }

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

    tmv::ConstBandMatrixView<T> a1cv = a1.view();
    tmv::BandMatrixView<T> a1v = a1.view();
    tmv::ConstBandMatrixView<T,tmv::FortranStyle> a1fcv = a1f.view();
    tmv::BandMatrixView<T,tmv::FortranStyle> a1fv = a1f.view();

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

    a1.resize(3,3,1,1);
    Assert(a1.colsize() == 3 && a1.rowsize() == 3,
           "BandMatrix a1.resize(3,3,1,1) sizes");
    Assert(a1.nlo() == 1 && a1.nhi() == 1,
           "BandMatrix a1.resize(3,3,1,1) nlo,nhi");
    for (int i=0, k=0; i<3; ++i) for (int j=0; j<3; ++j, ++k) 
        if ( j <= i + 1 && i <= j + 1) 
            a1(i,j) = T(k);
    for (int i=0, k=0; i<3; ++i) for (int j=0; j<3; ++j, ++k) 
        if ( j <= i + 1 && i <= j + 1) 
            Assert(a1(i,j) == k,"Read/Write resized BandMatrix");

    a1.resize(2*N,3*N,5,3);
    Assert(a1.colsize() == 2*N && a1.rowsize() == 3*N,
           "BandMatrix a1.resize(2*N,3*N,5,3) sizes");
    Assert(a1.nlo() == 5 && a1.nhi() == 3,
           "BandMatrix a1.resize(2*N,3*N,5,3) nlo,nhi");
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<3*N; ++j, ++k) 
        if ( j <= i + 3 && i <= j + 5) 
            a1(i,j) = T(k);
    for (int i=0, k=0; i<2*N; ++i) for (int j=0; j<3*N; ++j, ++k) 
        if ( j <= i + 3 && i <= j + 5) 
            Assert(a1(i,j) == k,"Read/Write resized BandMatrix");

}

template <class T, tmv::StorageType S> 
static void TestBasicBandMatrix_2()
{
    if (showstartdone) {
        std::cout<<"Start TestBasicBandMatrix_2\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
    }

    // Test assignments and constructors from arrays
    T qarrm[] = { 
        T(0), T(-1), T(-2),
        T(3), T(2),  T(1),  T(0),
              T(5),  T(4),  T(3),  T(2) 
    };
    T qarcm[] = {
        T(0),  T(3),
        T(-1), T(2),  T(5),
        T(-2), T(1),  T(4),
               T(0),  T(3),
                      T(2) 
    };
    T qardm[] = {        
           T(3),  T(5),
        T(0),  T(2),  T(4),
           T(-1), T(1),  T(3),
               T(-2), T(0),  T(2) 
    };
    std::vector<T> qvecrm(11);
    for(int i=0;i<11;i++) qvecrm[i] = qarrm[i];
    std::vector<T> qveccm(11);
    for(int i=0;i<11;i++) qveccm[i] = qarcm[i];
    std::vector<T> qvecdm(11);
    for(int i=0;i<11;i++) qvecdm[i] = qardm[i];

    tmv::BandMatrix<T,S> q1(3,5,1,2);
    std::copy(qarrm, qarrm+11, q1.rowmajor_begin());
    tmv::BandMatrix<T,S> q2(3,5,1,2);
    std::copy(qarcm, qarcm+11, q2.colmajor_begin());
    tmv::BandMatrix<T,S> q3(3,5,1,2);
    std::copy(qardm, qardm+11, q3.diagmajor_begin());

    tmv::BandMatrix<T,S> q4(3,5,1,2);
    std::copy(qvecrm.begin(), qvecrm.end(), q4.rowmajor_begin());
    tmv::BandMatrix<T,S> q5(3,5,1,2);
    std::copy(qveccm.begin(), qveccm.end(), q5.colmajor_begin());
    tmv::BandMatrix<T,S> q6(3,5,1,2);
    std::copy(qvecdm.begin(), qvecdm.end(), q6.diagmajor_begin());

    tmv::BandMatrix<T,S> q7x(30,40,10,20);
    tmv::BandMatrixView<T> q7 = q7x.subBandMatrix(3,18,5,30,1,2,5,5);
    std::copy(qvecrm.begin(), qvecrm.end(), q7.rowmajor_begin());

    tmv::BandMatrix<T,S> q8x(30,40,10,20);
    tmv::BandMatrixView<T> q8 = q8x.subBandMatrix(3,18,5,30,1,2,5,5);
    std::copy(qveccm.begin(), qveccm.end(), q8.colmajor_begin());

    tmv::BandMatrix<T,S> q9x(30,40,10,20);
    tmv::BandMatrixView<T> q9 = q9x.subBandMatrix(3,18,5,30,1,2,5,5);
    std::copy(qvecdm.begin(), qvecdm.end(), q9.diagmajor_begin());

    // Assignment using op<< is always in rowmajor order.
    tmv::BandMatrix<T,S> q10(3,5,1,2);
    tmv::BandMatrix<T,S> q11t(5,3,2,1);
    tmv::BandMatrixView<T> q11 = q11t.transpose();
    q10 <<
        0, -1, -2,
        3,  2,  1,  0,
            5,  4,  3,  2;
    q11 <<
        0, -1, -2,
        3,  2,  1,  0,
            5,  4,  3,  2;

    // Can also view memory directly
    T qarrmfull[] = {
               T(0),  T(-1), T(-2),
        T(3),  T(2),  T(1),  T(0),
        T(5),  T(4),  T(3),  T(2)
    };
    T qarcmfull[] = {
        T(0),  T(3),  T(6),
        T(-1), T(2),  T(5),
        T(-2), T(1),  T(4), 
        T(-3), T(0),  T(3),
        T(-4), T(-1), T(2),
    };
    T qardmfull[] = {
               T(3),  T(5),
        T(0),  T(2),  T(4),
        T(-1), T(1),  T(3),
        T(-2), T(0),  T(2) 
    };
    T* qarfull = 
        (S == tmv::RowMajor) ? qarrmfull : 
        (S == tmv::ColMajor) ? qarcmfull : 
        qardmfull;
    T* qarfullx = qarfull + (S == tmv::DiagMajor ? 2 : 0);
    const int Si = 
        (S == tmv::RowMajor) ? 3 :
        (S == tmv::ColMajor) ? 1 : 
        -2;
    const int Sj = 
        (S == tmv::RowMajor) ? 1 :
        (S == tmv::ColMajor) ? 3 : 
        3;

    const tmv::ConstBandMatrixView<T> q12 =
        tmv::BandMatrixViewOf(qarfull,3,5,1,2,S);
    const tmv::ConstBandMatrixView<T> q13 =
        tmv::BandMatrixViewOf(qarfullx,3,5,1,2,Si,Sj);

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
    }

    for(int i=0;i<3;i++) for(int j=0;j<5;j++) {
        if (j <= i+2 && i <= j+1) {
            T val = T(3*i-j);
            Assert(q1(i,j) == val,"Create BandMatrix from T* rm");
            Assert(q2(i,j) == val,"Create BandMatrix from T* cm");
            Assert(q3(i,j) == val,"Create BandMatrix from T* dm");
            Assert(q4(i,j) == val,"Create BandMatrix from vector rm");
            Assert(q5(i,j) == val,"Create BandMatrix from vector cm");
            Assert(q6(i,j) == val,"Create BandMatrix from vector dm");
            Assert(q7(i,j) == val,"Create BandMatrixView from vector rm");
            Assert(q8(i,j) == val,"Create BandMatrixView from vector cm");
            Assert(q9(i,j) == val,"Create BandMatrixView from vector dm");
            Assert(q10(i,j) == val,"Create BandMatrix from << list");
            Assert(q11(i,j) == val,"Create BandMatrixView from << list");
            Assert(q12(i,j) == val,"Create BandMatrixView of T* (S)");
            Assert(q13(i,j) == val,"Create BandMatrixView of T* (Si,Sj)");
        }
    }

    // Test the span of the iteration (i.e. the validity of begin(), end())
    const tmv::BandMatrix<T,S>& q1_const = q1;
    tmv::BandMatrixView<T> q1_view = q1.view();
    tmv::ConstBandMatrixView<T> q1_constview = q1_const.view();
    tmv::ConstBandMatrixView<T> q7_const = q7;

    typename tmv::BandMatrix<T,S>::rowmajor_iterator rmit1 = 
        q1.rowmajor_begin();
    typename tmv::BandMatrix<T,S>::const_rowmajor_iterator rmit2 =
        q1_const.rowmajor_begin();
    typename tmv::BandMatrixView<T>::rowmajor_iterator rmit3 =
        q1_view.rowmajor_begin();
    typename tmv::ConstBandMatrixView<T>::const_rowmajor_iterator rmit4 =
        q1_constview.rowmajor_begin();
    typename tmv::BandMatrixView<T>::rowmajor_iterator rmit5 =
        q7.rowmajor_begin();
    typename tmv::ConstBandMatrixView<T>::const_rowmajor_iterator rmit6 =
        q7_const.rowmajor_begin();
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
    Assert(i == 11, "RowMajor iteration number of elements");
    Assert(rmit2 == q1_const.rowmajor_end(), "rmit2 reaching end");
    Assert(rmit3 == q1_view.rowmajor_end(), "rmit3 reaching end");
    Assert(rmit4 == q1_constview.rowmajor_end(), "rmit4 reaching end");
    Assert(rmit5 == q7.rowmajor_end(), "rmit5 reaching end");
    Assert(rmit6 == q7_const.rowmajor_end(), "rmit6 reaching end");

    typename tmv::BandMatrix<T,S>::colmajor_iterator cmit1 = 
        q1.colmajor_begin();
    typename tmv::BandMatrix<T,S>::const_colmajor_iterator cmit2 =
        q1_const.colmajor_begin();
    typename tmv::BandMatrixView<T>::colmajor_iterator cmit3 =
        q1_view.colmajor_begin();
    typename tmv::ConstBandMatrixView<T>::const_colmajor_iterator cmit4 =
        q1_constview.colmajor_begin();
    typename tmv::BandMatrixView<T>::colmajor_iterator cmit5 =
        q7.colmajor_begin();
    typename tmv::ConstBandMatrixView<T>::const_colmajor_iterator cmit6 =
        q7_const.colmajor_begin();
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
    Assert(i == 11, "ColMajor iteration number of elements");
    Assert(cmit2 == q1_const.colmajor_end(), "cmit2 reaching end");
    Assert(cmit3 == q1_view.colmajor_end(), "cmit3 reaching end");
    Assert(cmit4 == q1_constview.colmajor_end(), "cmit4 reaching end");
    Assert(cmit5 == q7.colmajor_end(), "cmit5 reaching end");
    Assert(cmit6 == q7_const.colmajor_end(), "cmit6 reaching end");

    typename tmv::BandMatrix<T,S>::diagmajor_iterator dmit1 = 
        q1.diagmajor_begin();
    typename tmv::BandMatrix<T,S>::const_diagmajor_iterator dmit2 =
        q1_const.diagmajor_begin();
    typename tmv::BandMatrixView<T>::diagmajor_iterator dmit3 =
        q1_view.diagmajor_begin();
    typename tmv::ConstBandMatrixView<T>::const_diagmajor_iterator dmit4 =
        q1_constview.diagmajor_begin();
    typename tmv::BandMatrixView<T>::diagmajor_iterator dmit5 =
        q7.diagmajor_begin();
    typename tmv::ConstBandMatrixView<T>::const_diagmajor_iterator dmit6 =
        q7_const.diagmajor_begin();
    i = 0;
    while (dmit1 != q1.diagmajor_end()) {
        Assert(*dmit1++ == qardm[i], "DiagMajor iteration 1");
        Assert(*dmit2++ == qardm[i], "DiagMajor iteration 2");
        Assert(*dmit3++ == qardm[i], "DiagMajor iteration 3");
        Assert(*dmit4++ == qardm[i], "DiagMajor iteration 4");
        Assert(*dmit5++ == qardm[i], "DiagMajor iteration 5");
        Assert(*dmit6++ == qardm[i], "DiagMajor iteration 6");
        ++i;
    }
    Assert(i == 11, "DiagMajor iteration number of elements");
    Assert(dmit2 == q1_const.diagmajor_end(), "dmit2 reaching end");
    Assert(dmit3 == q1_view.diagmajor_end(), "dmit3 reaching end");
    Assert(dmit4 == q1_constview.diagmajor_end(), "dmit4 reaching end");
    Assert(dmit5 == q7.diagmajor_end(), "dmit5 reaching end");
    Assert(dmit6 == q7_const.diagmajor_end(), "dmit6 reaching end");

    // Test Basic Arithmetic
    const int N = 10;
    const int nhi = 1;
    const int nlo = 3;

    tmv::BandMatrix<T> a1(N,N,nlo,nhi);
    for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k)
        if ( j <= i + nhi && i <= j + nlo) {
            a1(i,j) = T(k);
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

    if (showstartdone) {
        std::cout<<"Start TestBasicBandMatrix_IO\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
        std::cout<<"nlo, nhi = "<<nlo<<','<<nhi<<std::endl;
    }

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
    ca1.writeCompact(fout);
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

    std::auto_ptr<tmv::Matrix<std::complex<T>,tmv::RowMajor> > xm4;
    std::auto_ptr<tmv::BandMatrix<std::complex<T>,tmv::RowMajor> > xb4;
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

#ifdef TEST_DOUBLE
template void TestBandMatrix<double>();
#endif
#ifdef TEST_FLOAT
template void TestBandMatrix<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandMatrix<long double>();
#endif
#ifdef TEST_INT
template void TestBandMatrix<int>();
#endif
