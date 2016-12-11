
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"
#include <fstream>
#include <cstdio>

#define CT std::complex<T>

template <class T> 
static void TestBasicDiagMatrix_1()
{
    const int N = 10;

    if (showstartdone) {
        std::cout<<"Start TestDiagMatrix\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
    }

    tmv::DiagMatrix<T> a(N);
    tmv::DiagMatrix<T,tmv::FortranStyle> af(N);
    Assert(a.colsize() == N && a.rowsize() == N,
           "Creating DiagMatrix(N)");
    Assert(af.colsize() == N && af.rowsize() == N,
           "Creating DiagMatrix(N)");

    for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if (i == j) a(i,j) = T(k);
        if (i == j) af(i+1,j+1) = T(k);
    }
    const tmv::DiagMatrix<T>& ac = a;
    tmv::ConstDiagMatrixView<T> acv = a.view();
    tmv::DiagMatrixView<T> av = a.view();
    const tmv::DiagMatrix<T,tmv::FortranStyle>& afc = af;
    tmv::ConstDiagMatrixView<T,tmv::FortranStyle> afcv = af.view();
    tmv::DiagMatrixView<T,tmv::FortranStyle> afv = af.view();

    for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k) {
        if (i == j) {
            Assert(a(i,j) == T(k),"Read/Write DiagMatrix");
            Assert(ac(i,j) == T(k),"Access const DiagMatrix");
            Assert(acv(i,j) == T(k),"Access DiagMatrix CV");
            Assert(av(i,j) == T(k),"Access DiagMatrix V");
            Assert(af(i+1,j+1) == T(k),"Read/Write DiagMatrixF");
            Assert(afc(i+1,j+1) == T(k),"Access const DiagMatrixF");
            Assert(afcv(i+1,j+1) == T(k),"Access DiagMatrixF CV");
            Assert(afv(i+1,j+1) == T(k),"Access DiagMatrixF V");
            Assert(a(i) == T(k),"Single argument access for DiagMatrix");
            Assert(acv(i) == T(k),"Single argument access for DiagMatrix CV");
            Assert(av(i) == T(k),"Single argument access for DiagMatrix V");
            Assert(af(i+1) == T(k),"Single argument access for DiagMatrixF");
            Assert(afcv(i+1) == T(k),"Single argument access for DiagMatrixF CV");
            Assert(afv(i+1) == T(k),"Single argument access for DiagMatrixF V");
        } else {
            Assert(ac(i,j) == T(0),"Access const DiagMatrix");
            Assert(acv(i,j) == T(0),"Access DiagMatrix CV");
            Assert(afc(i+1,j+1) == T(0),"Access const DiagMatrixF");
            Assert(afcv(i+1,j+1) == T(0),"Access DiagMatrixF CV");
        }
    }

    Assert(a==af,"CStyle Matrix == FortranStyle Matrix");
    Assert(a==acv,"Matrix == ConstMatrixView");
    Assert(a==av,"Matrix == MatrixView");
    Assert(a==afcv,"Matrix == FortranStyle ConstMatrixView");
    Assert(a==afv,"Matrix == FortranStyle MatrixView");

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
}

template <class T> 
static void TestBasicDiagMatrix_2()
{
    const int N = 10;

    if (showstartdone) {
        std::cout<<"Start TestBasicDiagMatrix_2\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
    }

    tmv::DiagMatrix<T> m(N);
    tmv::DiagMatrix<T,tmv::FortranStyle> mf(N);

    for (int i=0; i<N; ++i) {
        m(i,i) = T(i+1);
        mf(i+1,i+1) = T(i+1);
    }
    tmv::ConstDiagMatrixView<T> mcv = m.view();
    tmv::DiagMatrixView<T> mv = m.view();
    tmv::ConstDiagMatrixView<T,tmv::FortranStyle> mfcv = mf.view();
    tmv::DiagMatrixView<T,tmv::FortranStyle> mfv = mf.view();
    Assert(m.subDiagMatrix(2,5) == m.subDiagMatrix(2,5,1),"subDiagMatrix");
    Assert(m.subDiagMatrix(2,5) == mf.subDiagMatrix(3,5),"subDiagMatrixF");
    Assert(m.subDiagMatrix(2,8,2) == mf.subDiagMatrix(3,7,2),"subDiagMatrixF");
    Assert(m.subDiagMatrix(2,5) == mcv.subDiagMatrix(2,5),"subDiagMatrixCV");
    Assert(m.subDiagMatrix(2,8,2) == mcv.subDiagMatrix(2,8,2),
           "subDiagMatrixCV");
    Assert(m.subDiagMatrix(2,5) == mv.subDiagMatrix(2,5),"subDiagMatrixV");
    Assert(m.subDiagMatrix(2,8,2) == mv.subDiagMatrix(2,8,2),"subDiagMatrixV");
    Assert(mf.subDiagMatrix(3,5) == mfcv.subDiagMatrix(3,5),"subDiagMatrixFCV");
    Assert(mf.subDiagMatrix(3,7,2) == mfcv.subDiagMatrix(3,7,2),
           "subDiagMatrixFCV");
    Assert(mf.subDiagMatrix(3,5) == mfv.subDiagMatrix(3,5),"subDiagMatrixFV");
    Assert(mf.subDiagMatrix(3,7,2) == mfv.subDiagMatrix(3,7,2),
           "subDiagMatrixFV");

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
    int ii = 0;
    while (it1 != q1.end()) {
        Assert(*it1++ == qar[ii], "DiagMatrix iteration 1");
        Assert(*it2++ == qar[ii], "DiagMatrix iteration 2");
        Assert(*it3++ == qar[ii], "DiagMatrix iteration 3");
        Assert(*it4++ == qar[ii], "DiagMatrix iteration 4");
        Assert(*it5++ == qar[ii], "DiagMatrix iteration 5");
        Assert(*it6++ == qar[ii], "DiagMatrix iteration 6");
        ++ii;
    }
    Assert(ii == 3, "DiagMatrix iteration number of elements");
    Assert(it2 == q1_const.end(), "it2 reaching end");
    Assert(it3 == q1_view.end(), "it3 reaching end");
    Assert(it4 == q1_constview.end(), "it4 reaching end");
    Assert(it5 == q5.end(), "it5 reaching end");
    Assert(it6 == q5_const.end(), "it6 reaching end");

    // Test Basic Arithmetic
    tmv::DiagMatrix<T> a(N);
    tmv::DiagMatrix<T> b(N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
        if (i == j) {
            a(i,j) = T(3+i+5*j);
            b(i,j) = T(5+2*i+4*j);
        }
    }
    tmv::DiagMatrix<T,tmv::FortranStyle> af = a;
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

    tmv::Matrix<T> mm = a;
    for (int i=0, k=0; i<N; ++i) for (int j=0; j<N; ++j, ++k)
        if (i == j)
            Assert(a(i,j) == mm(i,j),"DiagMatrix -> Matrix");
    Assert(a == tmv::DiagMatrix<T>(mm),"Matrix -> DiagMatrix");

    if (!(std::numeric_limits<T>::is_integer)) {
        tmv::DiagMatrix<T> ainv = a;
        ainv.invertSelf();
        tmv::DiagMatrix<T> ainv2 = a.inverse();
        for(int i=0;i<N;++i)
            Assert(std::abs(a(i)*ainv(i) - T(1)) < 1.e-6,
                   "DiagMatrix invertSelf");
        for(int i=0;i<N;++i)
            Assert(std::abs(a(i)*ainv2(i) - T(1)) < 1.e-6,
                   "DiagMatrix inverse()");
    }
}

template <class T> 
static void TestBasicDiagMatrix_IO()
{
    const int N = 10;

    if (showstartdone) {
        std::cout<<"Start TestBasicDiagMatrix_IO\n";
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
    }

    tmv::DiagMatrix<T> m(N);
    tmv::DiagMatrix<CT> cm(N);

    for (int i=0;i<N;++i) {
        m(i,i) = T(i+1);
        cm(i,i) = CT(i+1,i+1001);
    }
    m(3) = T(1.e-30);
    cm(3) = CT(T(1.e-30),T(1.e-30));
    m(5) = T(9.e-3);
    cm(5) = CT(T(9.e-3),T(9.e-3));
    cm(6) = CT(T(9),T(9.e-3));
    m(7) = T(0.123456789);
    cm(7) = CT(T(3.123456789),T(6.987654321));

    // First check clipping function...
    tmv::DiagMatrix<T> m2 = m;
    tmv::DiagMatrix<CT> cm2 = cm;
    if (!std::numeric_limits<T>::is_integer) {
        m2.clip(T(1.e-2));
        cm2.clip(T(1.e-2));
    }
    tmv::DiagMatrix<T> m3 = m;
    tmv::DiagMatrix<CT> cm3 = cm;
    m3(3) = T(0);
    cm3(3) = T(0);
    m3(5) = T(0); 
    // Others, esp. cm3(5) and cm3(6), shouldn't get clipped.
    Assert(m2 == m3,"DiagMatrix clip");
    Assert(cm2 == cm3,"Complex DiagMatrix clip");

    // However, ThreshIO for complex works slightly differently than clip.
    // It clips _either_ the real or imag component, so now cm2(5) and
    // cm2(6) need to be modified.
    cm2(5) = cm3(5) = T(0);
    cm2(6) = cm3(6) = T(9);

    // Write matrices with 4 different styles
    std::ofstream fout("tmvtest_diagmatrix_io.dat");
    Assert(bool(fout),"Couldn't open tmvtest_diagmatrix_io.dat for output");
    fout << m << std::endl;
    fout << cm << std::endl;
    fout << tmv::CompactIO() << m << std::endl;
    fout << tmv::CompactIO() << cm << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << m << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << cm << std::endl;
    tmv::IOStyle myStyle = 
        tmv::CompactIO().setThresh(1.e-2).setPrecision(4).
        markup("Start","[",",","]","---","Done");
    fout << myStyle << m << std::endl;
    fout << myStyle << cm << std::endl;
    fout.close();

    // When using (the default) prec(6), these will be the values read in.
    m(7) = T(0.123457);
    cm(7) = CT(T(3.12346),T(6.98765));

    // When using prec(12), the full correct values will be read in. (m2,cm2)

    // When using prec(4), these will be the values read in.
    m3(7) = T(0.1235);
    if (std::numeric_limits<T>::is_integer) cm3(7) = CT(3,6);
    else cm3(7) = CT(T(3.123),T(6.988));

    // Read them back in
    tmv::DiagMatrix<T> xm1(N);
    tmv::DiagMatrix<CT> xcm1(N);
    std::ifstream fin("tmvtest_diagmatrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_diagmatrix_io.dat for input");
    fin >> xm1 >> xcm1;
    Assert(EqualIO(m,xm1,EPS),"DiagMatrix I/O check normal");
    Assert(EqualIO(cm,xcm1,EPS),"CDiagMatrix I/O check normal");
    fin >> tmv::CompactIO() >> xm1 >> tmv::CompactIO() >> xcm1;
    Assert(EqualIO(m,xm1,EPS),"DiagMatrix I/O check compact");
    Assert(EqualIO(cm,xcm1,EPS),"CDiagMatrix I/O check compact");
    fin >> xm1.view() >> xcm1.view();
    Assert(EqualIO(m2,xm1,EPS),"DiagMatrix I/O check thresh");
    Assert(EqualIO(cm2,xcm1,EPS),"CDiagMatrix I/O check thresh");
    fin >> myStyle >> xm1.view() >> myStyle >> xcm1.view();
    Assert(EqualIO(m3,xm1,EPS),"DiagMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(cm3,xcm1,EPS),"CDiagMatrix I/O check compact thresh & prec(4)");
    fin.close();

    // Repeat for matrices that need to be resized.
    // Also check switching the default IOStyle.
    tmv::CompactIO().makeDefault();
    tmv::DiagMatrix<T> zm1,zm2,zm3,zm4;
    tmv::DiagMatrix<CT> zcm1,zcm2,zcm3,zcm4;
    fin.open("tmvtest_diagmatrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_diagmatrix_io.dat for input");
    fin >> tmv::NormalIO() >> zm1 >> tmv::NormalIO() >> zcm1;
    Assert(EqualIO(m,zm1,EPS),"DiagMatrix I/O check normal");
    Assert(EqualIO(cm,zcm1,EPS),"CDiagMatrix I/O check normal");
    fin >> zm2 >> zcm2;
    Assert(EqualIO(m,zm2,EPS),"DiagMatrix I/O check compact");
    Assert(EqualIO(cm,zcm2,EPS),"CDiagMatrix I/O check compact");
    fin >> tmv::NormalIO() >> zm3 >> tmv::NormalIO() >> zcm3;
    Assert(EqualIO(m2,zm3,EPS),"DiagMatrix I/O check thresh");
    Assert(EqualIO(cm2,zcm3,EPS),"CDiagMatrix I/O check thresh");
    fin >> myStyle >> zm4 >> myStyle >> zcm4;
    Assert(EqualIO(m3,zm4,EPS),"DiagMatrix I/O check compact thresh & prec(4)");
    Assert(EqualIO(cm3,zcm4,EPS),"CDiagMatrix I/O check compact thresh & prec(4)");
    fin.close();
    tmv::IOStyle::revertDefault();

    // Finally, check that the NormalIO can be read in as a regular matrix.
    tmv::Matrix<T> zm5;
    tmv::Matrix<CT> zcm5;
    fin.open("tmvtest_diagmatrix_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_diagmatrix_io.dat for input");
    fin >> zm5 >> zcm5;
    Assert(EqualIO(m,zm5,EPS),"DiagMatrix -> Matrix I/O check");
    Assert(EqualIO(cm,zcm5,EPS),"CDiagMatrix -> CMatrix I/O check");
    fin.close();

#if XTEST == 0
    std::remove("tmvtest_diagmatrix_io.dat");
#endif
}

template <class T> void TestDiagMatrix()
{
    TestBasicDiagMatrix_1<T>();
    TestBasicDiagMatrix_2<T>();
    TestBasicDiagMatrix_IO<T>();

#if 1
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
#endif

    std::cout<<"DiagMatrix<"<<tmv::TMV_Text(T())<<"> passed all tests\n";
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
