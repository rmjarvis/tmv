#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include <fstream>
#include <cstdio>
#include <vector>

#define CT std::complex<T>

template <class T, int M, int N, tmv::StorageType S> 
inline void TestBasicSmallMatrix_1()
{
    if (showstartdone) {
        std::cout<<"Start TestBasicSmallMatrix_1\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"M,N = "<<M<<','<<N<<std::endl;
    }

    tmv::SmallMatrix<T,M,N,S> m;
    tmv::SmallMatrix<T,M,N,S|tmv::FortranStyle> mf;
    Assert(m.colsize() == M && m.rowsize() == N,
           "Creating SmallMatrix(M,N)");
    Assert(m.colsize() == M && m.rowsize() == N,
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

    tmv::SmallMatrix<CT,M,N,S> cm;
    tmv::SmallMatrix<CT,M,N,S> ca;
    tmv::SmallMatrix<CT,M,N,S> cb;
    Assert(cm.colsize() == M && cm.rowsize() == N,
           "Creating CSmallMatrix(M,N)");

    for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k)
        cm(i,j) = CT(T(k),T(k+1000));

    for (int i=0, k=0; i<M; ++i) for (int j=0; j<N; ++j, ++k) {
        Assert(cm(i,j) == CT(T(k),T(k+1000)),"Read/Write CSmallMatrix");
        Assert(cm.row(i)(j) == CT(T(k),T(k+1000)),"CSmallMatrix.row");
        Assert(cm.col(j)(i) == CT(T(k),T(k+1000)),"CSmallMatrix.col");
    }

    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) {
        ca(i,j) = CT(3+i+5*j,T(0)+i-j);
        cb(i,j) = CT(3+2*i+4*j,4-10*i);
    }

    cm = ca+cb;
    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
        Assert(cm(i,j) == CT(6+3*i+9*j,4-9*i-j),"Add CSmallMatrix");

    cm = ca-cb;
    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
        Assert(cm(i,j) == CT(T(0)-i+j,-4+11*i-j),"Subtract CSmallMatrix");

    cm = ca;
    for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) 
        Assert(cm(i,j) == ca(i,j),"Copy CSmallMatrix");
}

template <class T, int M, int N, tmv::StorageType S> 
inline void TestBasicSmallMatrix_2()
{
    if (showstartdone) {
        std::cout<<"Start TestBasicSmallMatrix_2\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"M,N = "<<M<<','<<N<<std::endl;
    }

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

    // Test the span of the iteration (i.e. the validity of begin(), end())
    const tmv::SmallMatrix<T,3,4,S>& q1_const = q1;
    typename tmv::SmallMatrix<T,3,4,S>::view_type q1_view = q1.view();
    typename tmv::SmallMatrix<T,3,4,S>::const_view_type q1_constview = 
        q1_const.view();
    tmv::ConstMatrixView<T> q5_const = q5;

    typename tmv::SmallMatrix<T,3,4,S>::rowmajor_iterator rmit1 = 
        q1.rowmajor_begin();
    typename tmv::SmallMatrix<T,3,4,S>::const_rowmajor_iterator rmit2 =
        q1_const.rowmajor_begin();
    typename tmv::SmallMatrix<T,3,4,S>::view_type::rowmajor_iterator rmit3 =
        q1_view.rowmajor_begin();
    typename tmv::SmallMatrix<T,3,4,S>::const_view_type::const_rowmajor_iterator
        rmit4 = q1_constview.rowmajor_begin();
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

    typename tmv::SmallMatrix<T,3,4,S>::colmajor_iterator cmit1 = 
        q1.colmajor_begin();
    typename tmv::SmallMatrix<T,3,4,S>::const_colmajor_iterator cmit2 =
        q1_const.colmajor_begin();
    typename tmv::SmallMatrix<T,3,4,S>::view_type::colmajor_iterator cmit3 =
        q1_view.colmajor_begin();
    typename tmv::SmallMatrix<T,3,4,S>::const_view_type::const_colmajor_iterator
        cmit4 = q1_constview.colmajor_begin();
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

template <class T, int M, int N, tmv::StorageType S> 
inline void TestBasicSmallMatrix_IO()
{
    if (showstartdone) {
        std::cout<<"Start TestBasicMatrix_IO\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"S = "<<tmv::TMV_Text(S)<<std::endl;
        std::cout<<"M,N = "<<M<<','<<N<<std::endl;
    }

    tmv::SmallMatrix<T,M,N,S> m;
    tmv::SmallMatrix<CT,M,N,S> cm;

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
    tmv::SmallMatrix<T,M,N> m2 = m;
    tmv::SmallMatrix<CT,M,N> cm2 = cm;
    if (!std::numeric_limits<T>::is_integer) {
        m2.clip(T(1.e-2));
        cm2.clip(T(1.e-2));
    }
    tmv::SmallMatrix<T,M,N> m3 = m;
    tmv::SmallMatrix<CT,M,N> cm3 = cm;
    m3(3,1) = T(0);
    cm3(3,1) = T(0);
    m3(5,6) = T(0); // Others, esp. cm3(5,6), shouldn't get clipped.
    Assert(m2 == m3,"SmallMatrix clip");
    Assert(cm2 == cm3,"Complex SmallMatrix clip");

    // However, ThreshIO for complex works slightly differently than clip.
    // It clips _either_ the real or imag component, so now cm2(5,6) and 
    // cm2(6,6) need to be modified.
    cm2(5,6) = cm3(5,6) = T(0);
    cm2(6,6) = cm3(6,6) = T(9);

    // Write matrices with 4 different styles
    std::ofstream fout("tmvtest_smallmatrix_io.dat");
    Assert(bool(fout),"Couldn't open tmvtest_smallmatrix_io.dat for output");
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
    tmv::SmallMatrix<T,M,N,tmv::RowMajor> xm1;
    tmv::SmallMatrix<CT,M,N,tmv::RowMajor> xcm1;
    std::ifstream fin("tmvtest_smallmatrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_smallmatrix_io.dat for input");
    fin >> xm1 >> xcm1;
    Assert(EqualIO(m,xm1,EPS),"SmallMatrix I/O check normal");
    Assert(EqualIO(cm,xcm1,EPS),"CSmallMatrix I/O check normal");
    fin >> tmv::CompactIO() >> xm1 >> tmv::CompactIO() >> xcm1;
    Assert(EqualIO(m,xm1,EPS),"SmallMatrix I/O check compact");
    Assert(EqualIO(cm,xcm1,EPS),"CSmallMatrix I/O check compact");
    fin >> xm1.view() >> xcm1.view();
    Assert(EqualIO(m2,xm1,EPS),"SmallMatrix I/O check thresh");
    Assert(EqualIO(cm2,xcm1,EPS),"CSmallMatrix I/O check thresh");
    fin >> myStyle >> xm1.view() >> myStyle >> xcm1.view();
    Assert(EqualIO(m3,xm1,EPS),"SmallMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(cm3,xcm1,EPS),"CSmallMatrix I/O check compact thresh & prec(4)");
    fin.close();

    // Repeat for column major
    tmv::SmallMatrix<T,M,N,tmv::ColMajor> xm2;
    tmv::SmallMatrix<CT,M,N,tmv::ColMajor> xcm2;
    fin.open("tmvtest_smallmatrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_smallmatrix_io.dat for input");
    fin >> xm2.view() >> xcm2.view();
    Assert(EqualIO(m,xm2,EPS),"SmallMatrix I/O check normal");
    Assert(EqualIO(cm,xcm2,EPS),"CSmallMatrix I/O check normal");
    fin >> tmv::CompactIO() >> xm2.view() >> tmv::CompactIO() >> xcm2.view();
    Assert(EqualIO(m,xm2,EPS),"SmallMatrix I/O check compact");
    Assert(EqualIO(cm,xcm2,EPS),"CSmallMatrix I/O check compact");
    fin >> xm2 >> xcm2;
    Assert(EqualIO(m2,xm2,EPS),"SmallMatrix I/O check thresh");
    Assert(EqualIO(cm2,xcm2,EPS),"CSmallMatrix I/O check thresh");
    fin >> myStyle >> xm2 >> myStyle >> xcm2;
    Assert(EqualIO(m3,xm2,EPS),"SmallMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(cm3,xcm2,EPS),"CSmallMatrix I/O check compact thresh & prec(4)");
    fin.close();

    // Read back into regular SmallMatrix
    // Also check switching the default IOStyle.
    tmv::CompactIO().makeDefault();
    tmv::Matrix<T> zm1,zm2,zm3,zm4;
    tmv::Matrix<CT> zcm1,zcm2,zcm3,zcm4;
    fin.open("tmvtest_smallmatrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_smallmatrix_io.dat for input");
    fin >> tmv::NormalIO() >> zm1 >> tmv::NormalIO() >> zcm1;
    Assert(EqualIO(m,zm1,EPS),"SmallMatrix I/O check normal -> Matrix");
    Assert(EqualIO(cm,zcm1,EPS),"CSmallMatrix I/O check normal -> Matrix");
    fin >> zm2 >> zcm2;
    Assert(EqualIO(m,zm2,EPS),"SmallMatrix I/O check compact -> Matrix");
    Assert(EqualIO(cm,zcm2,EPS),"CSmallMatrix I/O check compact -> Matrix");
    fin >> tmv::NormalIO() >> zm3 >> tmv::NormalIO() >> zcm3;
    Assert(EqualIO(m2,zm3,EPS),"SmallMatrix I/O check thresh -> Matrix");
    Assert(EqualIO(cm2,zcm3,EPS),"CSmallMatrix I/O check thresh -> Matrix");
    fin >> myStyle >> zm4 >> myStyle >> zcm4;
    Assert(EqualIO(m3,zm4,EPS),"SmallMatrix I/O check compact thresh -> Matrix");
    Assert(EqualIO(cm3,zcm4,EPS),"CSmallMatrix I/O check compact thresh -> Matrix");
    fin.close();
    // Switch it back.
    tmv::IOStyle::revertDefault();

#if XTEST == 0
    std::remove("tmvtest_smallmatrix_io.dat");
#endif

}

template <class T> 
void TestAllSmallMatrix()
{
    TestBasicSmallMatrix_1<T,15,10,tmv::RowMajor>();
    TestBasicSmallMatrix_1<T,15,10,tmv::ColMajor>();
    TestBasicSmallMatrix_2<T,15,10,tmv::RowMajor>();
    TestBasicSmallMatrix_2<T,15,10,tmv::ColMajor>();
    TestBasicSmallMatrix_IO<T,15,10,tmv::RowMajor>();
    TestBasicSmallMatrix_IO<T,15,10,tmv::ColMajor>();
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
