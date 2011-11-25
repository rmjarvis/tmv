
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"
#include <fstream>
#include <cstdio>

// Break this part out, since we need to skip it for integers
template <class T> void TestDiagMatrixInvert()
{
    const int N = 10;
    tmv::DiagMatrix<T> a(N);
    for (int i=0; i<N; ++i) a(i,i) = T(3+5*i);

    tmv::DiagMatrix<T> ainv = a;
    ainv.invertSelf();
    tmv::DiagMatrix<T> ainv2 = a.inverse();
    for(int i=0;i<N;++i)
        Assert(std::abs(a(i)*ainv(i) - T(1)) < 1.e-6,"DiagMatrix invertSelf");
    for(int i=0;i<N;++i)
        Assert(std::abs(a(i)*ainv2(i) - T(1)) < 1.e-6,"DiagMatrix inverse()");
}

template <> void TestDiagMatrixInvert<int>()
{}

template <class T> void TestDiagMatrix()
{
    const int N = 10;

    tmv::DiagMatrix<T> a(N);
    tmv::DiagMatrix<T,tmv::FortranStyle> af(N);
    Assert(a.colsize() == size_t(N) && a.rowsize() == size_t(N),
           "Creating DiagMatrix(N)");
    Assert(af.colsize() == size_t(N) && af.rowsize() == size_t(N),
           "Creating DiagMatrix(N)");

    for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if (i == j) a(i,j) = T(k);
        if (i == j) af(i+1,j+1) = T(k);
    }
    tmv::ConstDiagMatrixView<T> acv = a.view();
    tmv::DiagMatrixView<T> av = a.view();
    tmv::ConstDiagMatrixView<T,tmv::FortranStyle> afcv = af.view();
    tmv::DiagMatrixView<T,tmv::FortranStyle> afv = af.view();

    for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k)
        if (i == j) {
            Assert(a(i,j) == k,"Read/Write DiagMatrix");
            Assert(acv(i,j) == k,"Access DiagMatrix CV");
            Assert(av(i,j) == k,"Access DiagMatrix V");
            Assert(af(i+1,j+1) == k,"Read/Write DiagMatrixF");
            Assert(afcv(i+1,i+1) == k,"Access DiagMatrixF CV");
            Assert(afv(i+1,i+1) == k,"Access DiagMatrixF V");
            Assert(a(i) == k,"Single argument access for DiagMatrix");
            Assert(acv(i) == k,"Single argument access for DiagMatrix CV");
            Assert(av(i) == k,"Single argument access for DiagMatrix V");
            Assert(af(i+1) == k,"Single argument access for DiagMatrixF");
            Assert(afcv(i+1) == k,"Single argument access for DiagMatrixF CV");
            Assert(afv(i+1) == k,"Single argument access for DiagMatrixF V");
        }

    Assert(a==af,"CStyle Matrix == FortranStyle Matrix");
    Assert(a==acv,"Matrix == ConstMatrixView");
    Assert(a==av,"Matrix == MatrixView");
    Assert(a==afcv,"Matrix == FortranStyle ConstMatrixView");
    Assert(a==afv,"Matrix == FortranStyle MatrixView");


    // Test assignments and constructors from arrays
    T qar[] = { T(0), T(3), T(6) };
    T qar2[] = { T(0), T(1), T(2), T(3), T(4), T(5), T(6) };
    std::vector<T> qvec(3);
    for(int i=0;i<3;i++) qvec[i] = qar[i];

    tmv::DiagMatrix<T> q1(3);
    std::copy(qar, qar+3, q1.begin());

    tmv::DiagMatrix<T> q2(3);
    std::copy(qvec.begin(), qvec.end(), q2.begin());

    tmv::DiagMatrix<T> q3x(30);
    tmv::DiagMatrixView<T> q3 = q3x.subDiagMatrix(3,18,5);
    std::copy(qvec.begin(), qvec.end(), q3.begin());

    tmv::DiagMatrix<T> q4(3);
    tmv::DiagMatrix<T> q5x(30);
    tmv::DiagMatrixView<T> q5 = q5x.subDiagMatrix(3,18,5);
    q4 <<
        0,
           3,
              6;
    q5 <<
        0,
           3,
              6;

    tmv::ConstDiagMatrixView<T> q6 = tmv::DiagMatrixViewOf(qar,3);
    tmv::ConstDiagMatrixView<T> q7 = tmv::DiagMatrixViewOf(qar2,3,3);

    if (showacc) {
        std::cout<<"q1 = "<<q1<<std::endl;
        std::cout<<"q2 = "<<q2<<std::endl;
        std::cout<<"q3 = "<<q3<<std::endl;
        std::cout<<"q4 = "<<q4<<std::endl;
        std::cout<<"q5 = "<<q5<<std::endl;
        std::cout<<"q6 = "<<q6<<std::endl;
        std::cout<<"q7 = "<<q7<<std::endl;
    }

    for(int i=0;i<3;i++) {
        Assert(q1(i,i) == T(3*i),"Create DiagMatrix from T*");
        Assert(q2(i,i) == T(3*i),"Create DiagMatrix from vector");
        Assert(q3(i,i) == T(3*i),"Create DiagMatrixView from vector");
        Assert(q4(i,i) == T(3*i),"Create DiagMatrix from << list");
        Assert(q5(i,i) == T(3*i),"Create DiagMatrixView from << list");
        Assert(q6(i,i) == T(3*i),"Create DiagMatrixView from T*");
        Assert(q7(i,i) == T(3*i),"Create DiagMatrixView from T* with step");
    }

    // Test the span of the iteration (i.e. the validity of begin(), end())
    const tmv::DiagMatrix<T>& q1_const = q1;
    tmv::DiagMatrixView<T> q1_view = q1.view();
    tmv::ConstDiagMatrixView<T> q1_constview = q1_const.view();
    tmv::ConstDiagMatrixView<T> q5_const = q5;

    typename tmv::DiagMatrix<T>::iterator it1 = q1.begin();
    typename tmv::DiagMatrix<T>::const_iterator it2 = q1_const.begin();
    typename tmv::DiagMatrixView<T>::iterator it3 = q1_view.begin();
    typename tmv::ConstDiagMatrixView<T>::const_iterator it4 =
        q1_constview.begin();
    typename tmv::DiagMatrixView<T>::iterator it5 = q5.begin();
    typename tmv::ConstDiagMatrixView<T>::const_iterator it6 =
        q5_const.begin();
    int i = 0;
    while (it1 != q1.end()) {
        Assert(*it1++ == qar[i], "DiagMatrix iteration 1");
        Assert(*it2++ == qar[i], "DiagMatrix iteration 2");
        Assert(*it3++ == qar[i], "DiagMatrix iteration 3");
        Assert(*it4++ == qar[i], "DiagMatrix iteration 4");
        Assert(*it5++ == qar[i], "DiagMatrix iteration 5");
        Assert(*it6++ == qar[i], "DiagMatrix iteration 6");
        ++i;
    }
    Assert(i == 3, "DiagMatrix iteration number of elements");
    Assert(it2 == q1_const.end(), "it2 reaching end");
    Assert(it3 == q1_view.end(), "it3 reaching end");
    Assert(it4 == q1_constview.end(), "it4 reaching end");
    Assert(it5 == q5.end(), "it5 reaching end");
    Assert(it6 == q5_const.end(), "it6 reaching end");


    // Test Basic Arithmetic
    tmv::DiagMatrix<T> b(N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        if (i == j) {
            a(i,j) = T(3+i+5*j);
            b(i,j) = T(5+2*i+4*j);
        }
    af = a;
    Assert(a==af,"Copy CStyle DiagMatrix to FotranStyle");

    tmv::DiagMatrix<T> c(N);
    c = a+b;
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        if (i == j)
            Assert(c(i,j) == T(8+3*i+9*j),"Add DiagMatrices");

    c = a-b;
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        if (i == j)
            Assert(c(i,j) == T(-2-i+j),"Subtract DiagMatrices");

    tmv::Matrix<T> m = a;
    for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k)
        if (i == j)
            Assert(a(i,j) == m(i,j),"DiagMatrix -> Matrix");
    Assert(a == tmv::DiagMatrix<T>(m),"Matrix -> DiagMatrix");


    tmv::DiagMatrix<std::complex<T> > ca = a*std::complex<T>(1,2);
    tmv::DiagMatrix<std::complex<T> > cb = b*std::complex<T>(-5,-1);

    a.resize(2);
    Assert(a.size() == 2,"DiagMatrix a.resize(2)");
    for (int i=0; i<2; ++i) a(i,i) = T(i);
    for (int i=0; i<2; ++i) 
        Assert(a(i,i) == i,"Read/Write resized DiagMatrix");

    a.resize(2*N);
    Assert(a.size() == 2*N,"DiagMatrix a.resize(2*N)");
    for (int i=0; i<2*N; ++i) a(i,i) = T(i);
    for (int i=0; i<2*N; ++i) 
        Assert(a(i,i) == i,"Read/Write resized DiagMatrix");

    // Test I/O

    std::ofstream fout("tmvtest_diagmatrix_io.dat");
    if (!fout) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_diagmatrix_io.dat for output\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_diagmatrix_io.dat for output");
#endif
    }
    fout << ca << std::endl;
    ca.writeCompact(fout);
    fout.close();

    tmv::Matrix<std::complex<T> > xcm1(N,N);
    tmv::DiagMatrix<std::complex<T> > xcd1(N);
    std::ifstream fin("tmvtest_diagmatrix_io.dat");
    if (!fin) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_diagmatrix_io.dat for input\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_diagmatrix_io.dat for input");
#endif
    }
    fin >> xcm1 >> xcd1;
    fin.close();
    Assert(tmv::Matrix<std::complex<T> >(ca) == xcm1,"DiagMatrix I/O check #1");
    Assert(ca == xcd1,"DiagMatrix Compact I/O check #1");

    std::auto_ptr<tmv::Matrix<std::complex<T> > > xcm2;
    std::auto_ptr<tmv::DiagMatrix<std::complex<T> > > xcd2;
    fin.open("tmvtest_diagmatrix_io.dat");
    if (!fin) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_diagmatrix_io.dat for input\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_diagmatrix_io.dat for input");
#endif
    }
    fin >> xcm2 >> xcd2;
    fin.close();
    Assert(tmv::Matrix<std::complex<T> >(ca) == *xcm2,"DiagMatrix I/O check #2");
    Assert(ca == *xcd2,"DiagMatrix Compact I/O check #2");

#ifndef XTEST
    std::remove("tmvtest_diagmatrix_io.dat");
#endif

    TestDiagMatrixInvert<T>();

    TestDiagMatrixArith_A1<T>();
    TestDiagMatrixArith_A2<T>();
    TestDiagMatrixArith_A3<T>();
    TestDiagMatrixArith_A4<T>();
    TestDiagMatrixArith_A5<T>();
    TestDiagMatrixArith_A6<T>();
    TestDiagMatrixArith_B4a<T>();
    TestDiagMatrixArith_B4b<T>();
    TestDiagMatrixArith_B5a<T>();
    TestDiagMatrixArith_B5b<T>();
    TestDiagMatrixArith_B6a<T>();
    TestDiagMatrixArith_B6b<T>();

    std::cout<<"DiagMatrix<"<<Text(T())<<"> passed all tests\n";
}


#ifdef TEST_DOUBLE
template void TestDiagMatrix<double>();
#endif
#ifdef TEST_FLOAT
template void TestDiagMatrix<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestDiagMatrix<long double>();
#endif
#ifdef TEST_INT
template void TestDiagMatrix<int>();
#endif
