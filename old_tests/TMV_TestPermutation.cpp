
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"

#define NO_COMPLEX_ARITH
#define NOASSIGN

#include "TMV_TestMatrixArith.h"
#include <fstream>

template <class T>
void TestPermutation()
{
    T eps = EPS;

    ptrdiff_t pp1[10] = { 0,1,2,3,4,5,6,7,8,9 }; // identity
    ptrdiff_t pp2[10] = { 9,8,7,6,5,5,6,7,8,9 }; // reversal
    ptrdiff_t pp3[10] = { 5,3,5,9,4,5,9,8,9,9 }; // "random"
    ptrdiff_t pp3i[10] = { 2,6,5,6,4,5,8,9,9,9 }; // p3 in inverse order

    int det1 = 1;
    int det2 = -1;
    int det3 = -1;

    tmv::Vector<T> v0(10);
    v0 << 0,1,2,3,4,5,6,7,8,9;
    tmv::Vector<T> v1 = v0;
    tmv::Vector<T> v2 = v0.reverse();
    tmv::Vector<T> v3l(10);
    v3l << 5,3,0,9,4,2,1,8,6,7;
    tmv::Vector<T> v3r(10);
    v3r << 2,6,5,1,4,0,8,9,7,3;

    tmv::Permutation p1(10,pp1,false);
    tmv::Permutation p2(10,pp2,false);
    tmv::Permutation p3(10,pp3,false);
    tmv::Permutation p3i(10,pp3i,true);

    tmv::Matrix<T> m1 = p1;
    tmv::Matrix<T> m2 = p2;
    tmv::Matrix<T> m3 = p3;
    tmv::Matrix<T> m3i = p3i;

    Assert(m3 == m3i,"Inverse permutation definition.");

    //
    // Test construction
    //

    Assert(p1.det() == det1,"Identity permutation determinant");
    Assert(p2.det() == det2,"Reversal permutation determinant");
    Assert(p3.det() == det3,"Random permutation determinant");
    Assert(p3i.det() == det3,"Inverse permutation determinant");

    Assert(p1.det() == m1.det(),"Identity permutation determinant = m.det");
    Assert(p2.det() == m2.det(),"Reversal permutation determinant = m.det");
    Assert(p3.det() == m3.det(),"Random permutation determinant = m.det");
    Assert(p3i.det() == m3.det(),"Inverse permutation determinant = m.det");

    for(int i=0;i<10;++i) for(int j=0;j<10;++j) {
        Assert(p1.cref(i,j) == m1.cref(i,j),"Identity permutation cref");
        Assert(p2.cref(i,j) == m2.cref(i,j),"Reverse permutation cref");
        Assert(p3.cref(i,j) == m3.cref(i,j),"Random permutation cref");
        Assert(p3i.cref(i,j) == m3.cref(i,j),"Inverse permutation cref");
    }

    //
    // Test arithmetic with Vector
    //

    tmv::Vector<T> v(10);
    v = p1 * v0;
    Assert(v == v1,"Identity permutation from left");
    v = v0 * p1;
    Assert(v == v1,"Identity permutation from right");
    v = p2 * v0;
    Assert(v == v2,"Reversal permutation from left");
    v = v0 * p2;
    Assert(v == v2,"Reversal permutation from right");
    v = p3 * v0;
    Assert(v == v3l,"Random permutation from left");
    v = v0 * p3;
    Assert(v == v3r,"Random permutation from right");
    v = p3i * v0;
    Assert(v == v3l,"Inverse permutation from left");
    v = v0 * p3i;
    Assert(v == v3r,"Inverse permutation from right");

    v = m1 * v0;
    Assert(v == v1,"Identity permutation from left -- matrix version");
    v = v0 * m1;
    Assert(v == v1,"Identity permutation from right -- matrix version");
    v = m2 * v0;
    Assert(v == v2,"Reversal permutation from left -- matrix version");
    v = v0 * m2;
    Assert(v == v2,"Reversal permutation from right -- matrix version");
    v = m3 * v0;
    Assert(v == v3l,"Random permutation from left -- matrix version");
    v = v0 * m3;
    Assert(v == v3r,"Random permutation from right -- matrix version");

    v = p1.transpose() * v0;
    Assert(v == v1,"Identity permutation from left -- transpose");
    v = v0 * p1.transpose();
    Assert(v == v1,"Identity permutation from right -- transpose");
    v = p2.transpose() * v0;
    Assert(v == v2,"Reversal permutation from left -- transpose");
    v = v0 * p2.transpose();
    Assert(v == v2,"Reversal permutation from right -- transpose");
    v = p3.transpose() * v0;
    Assert(v == v3r,"Random permutation from left -- transpose");
    v = v0 * p3.transpose();
    Assert(v == v3l,"Random permutation from right -- transpose");
    v = p3i.transpose() * v0;
    Assert(v == v3r,"Inverse permutation from left -- transpose");
    v = v0 * p3i.transpose();
    Assert(v == v3l,"Inverse permutation from right -- transpose");

    v = m1.transpose() * v0;
    Assert(v == v1,"Identity permutation from left -- transpose matrix");
    v = v0 * m1.transpose();
    Assert(v == v1,"Identity permutation from right -- transpose matrix");
    v = m2.transpose() * v0;
    Assert(v == v2,"Reversal permutation from left -- transpose matrix");
    v = v0 * m2.transpose();
    Assert(v == v2,"Reversal permutation from right -- transpose matrix");
    v = m3.transpose() * v0;
    Assert(v == v3r,"Random permutation from left -- transpose matrix");
    v = v0 * m3.transpose();
    Assert(v == v3l,"Random permutation from right -- transpose matrix");

    v = v0 / p1;
    Assert(v == v1,"Identity permutation left division");
    v = v0 % p1;
    Assert(v == v1,"Identity permutation right division");
    v = v0 / p2;
    Assert(v == v2,"Reversal permutation left division");
    v = v0 % p2;
    Assert(v == v2,"Reversal permutation right division");
    v = v0 / p3;
    Assert(v == v3r,"Random permutation left division");
    v = v0 % p3;
    Assert(v == v3l,"Random permutation right division");
    v = v0 / p3i;
    Assert(v == v3r,"Inverse permutation left division");
    v = v0 % p3i;
    Assert(v == v3l,"Inverse permutation right division");

    v = v0;
    v *= p3;
    Assert(v == v3r,"Random permutation *=");
    v %= p3;
    Assert(v == v0,"Random permutation %=");
    v /= p3;
    Assert(v == v3r,"Random permutation /=");

    v = v0;
    v *= p3i;
    Assert(v == v3r,"Inverse permutation *=");
    v %= p3i;
    Assert(v == v0,"Inverse permutation %=");
    v /= p3i;
    Assert(v == v3r,"Inverse permutation /=");

    //
    // Test arithmetic with a Matrix
    //

    tmv::Matrix<T> m0(10,10);
    for(int i=0;i<10;++i) for(int j=0;j<10;++j) 
        m0(i,j) = T(2)+T(3)*i-T(5)*j;
    tmv::Matrix<T> m(10,10);
    tmv::Matrix<std::complex<T> > cm0 = std::complex<T>(3,-2) * m0;
    tmv::Matrix<std::complex<T> > cm(10,10);

    m = p1 * m0;
    Assert(Equal(m,(m1*m0),eps),"Identity permutation p*m");
    m = p2 * m0;
    Assert(Equal(m,(m2*m0),eps),"Reversal permutation p*m");
    m = p3 * m0;
    Assert(Equal(m,(m3*m0),eps),"Random permutation p*m");
    m = p3i * m0;
    Assert(Equal(m,(m3*m0),eps),"Inverse permutation p*m");
    m = m0 * p1;
    Assert(Equal(m,(m0*m1),eps),"Identity permutation m*p");
    m = m0 * p2;
    Assert(Equal(m,(m0*m2),eps),"Reversal permutation m*p");
    m = m0 * p3;
    Assert(Equal(m,(m0*m3),eps),"Random permutation m*p");
    m = m0 * p3i;
    Assert(Equal(m,(m0*m3),eps),"Inverse permutation m*p");

    m = p1.transpose() * m0;
    Assert(Equal(m,(m1.transpose()*m0),eps),"Identity permutation pt*m");
    m = p2.transpose() * m0;
    Assert(Equal(m,(m2.transpose()*m0),eps),"Reversal permutation pt*m");
    m = p3.transpose() * m0;
    Assert(Equal(m,(m3.transpose()*m0),eps),"Random permutation pt*m");
    m = p3i.transpose() * m0;
    Assert(Equal(m,(m3.transpose()*m0),eps),"Inverse permutation pt*m");
    m = m0 * p1.transpose();
    Assert(Equal(m,(m0*m1.transpose()),eps),"Identity permutation m*pt");
    m = m0 * p2.transpose();
    Assert(Equal(m,(m0*m2.transpose()),eps),"Reversal permutation m*pt");
    m = m0 * p3.transpose();
    Assert(Equal(m,(m0*m3.transpose()),eps),"Random permutation m*pt");
    m = m0 * p3i.transpose();
    Assert(Equal(m,(m0*m3.transpose()),eps),"Inverse permutation m*pt");

    // The permutation division works for integers, but the comparison
    // to the regular matrix fails.
    if (!(std::numeric_limits<T>::is_integer)) {
        m = m0 / p1;
        Assert(Equal(m,(m0/m1),eps),"Identity permutation m/p");
        m = m0 / p2;
        Assert(Equal(m,(m0/m2),eps),"Reversal permutation m/p");
        m = m0 / p3;
        Assert(Equal(m,(m0/m3),eps),"Random permutation m/p");
        m = m0 / p3i;
        Assert(Equal(m,(m0/m3),eps),"Inverse permutation m/p");

        m = m0 % p1;
        Assert(Equal(m,(m0%m1),eps),"Identity permutation m%p");
        m = m0 % p2;
        Assert(Equal(m,(m0%m2),eps),"Reversal permutation m%p");
        m = m0 % p3;
        Assert(Equal(m,(m0%m3),eps),"Random permutation m%p");
        m = m0 % p3i;
        Assert(Equal(m,(m0%m3),eps),"Inverse permutation m%p");
    }

    m = m0;
    m *= p3;
    Assert(Equal(m,(m0*m3),eps),"Random permutation m*=p");
    if (!(std::numeric_limits<T>::is_integer)) {
        m = m0;
        m %= p3;
        Assert(Equal(m,(m0%m3),eps),"Random permutation m%=p");
        m = m0;
        m /= p3;
        Assert(Equal(m,(m0/m3),eps),"Random permutation m/=p");
    }

    m = m0;
    m *= p3i;
    Assert(Equal(m,(m0*m3),eps),"Inverse permutation m*=p");
    if (!(std::numeric_limits<T>::is_integer)) {
        m = m0;
        m %= p3i;
        Assert(Equal(m,(m0%m3),eps),"Inverse permutation m%=p");
        m = m0;
        m /= p3i;
        Assert(Equal(m,(m0/m3),eps),"Inverse permutation m/=p");
    }

    // Repeat for complex matrices:
    cm = p1 * cm0;
    Assert(Equal(cm,(m1*cm0),eps),"Identity permutation p*cm");
    cm = p2 * cm0;
    Assert(Equal(cm,(m2*cm0),eps),"Reversal permutation p*cm");
    cm = p3 * cm0;
    Assert(Equal(cm,(m3*cm0),eps),"Random permutation p*cm");
    cm = p3i * cm0;
    Assert(Equal(cm,(m3*cm0),eps),"Inverse permutation p*cm");
    cm = cm0 * p1;
    Assert(Equal(cm,(cm0*m1),eps),"Identity permutation cm*p");
    cm = cm0 * p2;
    Assert(Equal(cm,(cm0*m2),eps),"Reversal permutation cm*p");
    cm = cm0 * p3;
    Assert(Equal(cm,(cm0*m3),eps),"Random permutation cm*p");
    cm = cm0 * p3i;
    Assert(Equal(cm,(cm0*m3),eps),"Inverse permutation cm*p");

    cm = p1.transpose() * cm0;
    Assert(Equal(cm,(m1.transpose()*cm0),eps),"Identity permutation pt*cm");
    cm = p2.transpose() * cm0;
    Assert(Equal(cm,(m2.transpose()*cm0),eps),"Reversal permutation pt*cm");
    cm = p3.transpose() * cm0;
    Assert(Equal(cm,(m3.transpose()*cm0),eps),"Random permutation pt*cm");
    cm = p3i.transpose() * cm0;
    Assert(Equal(cm,(m3.transpose()*cm0),eps),"Inverse permutation pt*cm");
    cm = cm0 * p1.transpose();
    Assert(Equal(cm,(cm0*m1.transpose()),eps),"Identity permutation cm*pt");
    cm = cm0 * p2.transpose();
    Assert(Equal(cm,(cm0*m2.transpose()),eps),"Reversal permutation cm*pt");
    cm = cm0 * p3.transpose();
    Assert(Equal(cm,(cm0*m3.transpose()),eps),"Random permutation cm*pt");
    cm = cm0 * p3i.transpose();
    Assert(Equal(cm,(cm0*m3.transpose()),eps),"Inverse permutation cm*pt");

    if (!(std::numeric_limits<T>::is_integer)) {
        cm = cm0 / p1;
        Assert(Equal(cm,(cm0/m1),eps),"Identity permutation cm/p");
        cm = cm0 / p2;
        Assert(Equal(cm,(cm0/m2),eps),"Reversal permutation cm/p");
        cm = cm0 / p3;
        Assert(Equal(cm,(cm0/m3),eps),"Random permutation cm/p");
        cm = cm0 / p3i;
        Assert(Equal(cm,(cm0/m3),eps),"Inverse permutation cm/p");

        cm = cm0 % p1;
        Assert(Equal(cm,(cm0%m1),eps),"Identity permutation cm%p");
        cm = cm0 % p2;
        Assert(Equal(cm,(cm0%m2),eps),"Reversal permutation cm%p");
        cm = cm0 % p3;
        Assert(Equal(cm,(cm0%m3),eps),"Random permutation cm%p");
        cm = cm0 % p3i;
        Assert(Equal(cm,(cm0%m3),eps),"Inverse permutation cm%p");
    }

    cm = cm0;
    cm *= p3;
    Assert(Equal(cm,(cm0*m3),eps),"Random permutation cm*=p");
    if (!(std::numeric_limits<T>::is_integer)) {
        cm = cm0;
        cm %= p3;
        Assert(Equal(cm,(cm0%m3),eps),"Random permutation cm%=p");
        cm = cm0;
        cm /= p3;
        Assert(Equal(cm,(cm0/m3),eps),"Random permutation cm/=p");
    }

    cm = cm0;
    cm *= p3i;
    Assert(Equal(cm,(cm0*m3),eps),"Inverse permutation cm*=p");
    if (!(std::numeric_limits<T>::is_integer)) {
        cm = cm0;
        cm %= p3i;
        Assert(Equal(cm,(cm0%m3),eps),"Inverse permutation cm%=p");
        cm = cm0;
        cm /= p3i;
        Assert(Equal(cm,(cm0/m3),eps),"Inverse permutation cm/=p");
    }

    // 
    // Test various functions of a Permutation
    //

    p1 = p2;
    Assert(Equal(tmv::Matrix<T>(p1),m2,eps),"Permutation op=");
    Assert(p1 == p2,"Permutation op==");
    Assert(p1 != p3,"Permutation op!=");
    Assert(p3 == p3i,"Permutation op==, different storage order");
    Assert(p1 != p3i,"Permutation op!=, different storage order");
    p1 = p3.transpose();
    Assert(Equal(tmv::Matrix<T>(p1),m3.transpose(),eps),
           "Permutation op= pt");
    Assert(p1 == p3.transpose(),"Permutation op== transpose");
    p1.transposeSelf();
    Assert(Equal(tmv::Matrix<T>(p1),m3,eps),"Permutation transposeSelf");
    p3.setToIdentity();
    Assert(Equal(tmv::Matrix<T>(p3),m1,eps),"Permutation setToIdentity");
    Swap(p1,p3);
    Assert(Equal(tmv::Matrix<T>(p1),m1,eps),"Permutation Swap");
    Assert(Equal(tmv::Matrix<T>(p3),m3,eps),"Permutation Swap");
    tmv::Permutation p1x(10);
    Assert(p1x == p1,"Permutation identity constructor");

    if (!(std::numeric_limits<T>::is_integer)) {
        // Norm(p) returns a double, so don't use a smaller eps than that.
        T deps = T(std::max(double(eps), 10.*tmv::TMV_Epsilon<double>()));
        Assert(Equal2(Norm(p3),Norm(m3),deps),"Permutation Norm");
        Assert(Equal2(NormF(p3),NormF(m3),deps),"Permutation NormF");
    }
    Assert(Equal2(NormSq(p3),NormSq(m3),eps),"Permutation NormSq");
    Assert(Equal2(Norm1(p3),Norm1(m3),eps),"Permutation Norm1");
    if (!(std::numeric_limits<T>::is_integer)) {
        Assert(Equal2(Norm2(p3),Norm2(m3),eps),"Permutation Norm2");
    }
    Assert(Equal2(NormInf(p3),NormInf(m3),eps),"Permutation NormInf");
    Assert(Equal2(MaxAbsElement(p3),MaxAbsElement(m3),eps),
           "Permutation MaxAbsElement");
    Assert(Equal2(Trace(p3),Trace(m3),eps),"Permutation Trace");
    Assert(Equal2(Det(p3),Det(m3),eps),"Permutation Det");
    if (!(std::numeric_limits<T>::is_integer)) {
        Assert(Equal2(LogDet(p3),LogDet(m3),eps),"Permutation LogDet");
    }


    //
    // Test I/O
    //

    std::ofstream fout("tmvtest_permutation_io.dat");
    Assert(bool(fout),"Couldn't open tmvtest_permutation_io.dat for output");
    fout << p3 << std::endl;
    fout << tmv::CompactIO() << p3 << std::endl;
    fout.close();

    tmv::Matrix<T> xm1(10,10);
    tmv::Permutation xp1(10);
    std::ifstream fin("tmvtest_permutation_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_permutation_io.dat for input");
    fin >> xm1 >> tmv::CompactIO() >> xp1;
    fin.close();
    Assert(m3 == xm1,"Permutation I/O check #1");
    Assert(p3 == xp1,"Permutation Compact I/O check #1");

    tmv::Matrix<T> xm2;
    tmv::Permutation xp2;
    fin.open("tmvtest_permutation_io.dat");
    fin >> xm2 >> tmv::CompactIO() >> xp2;
    fin.close();
    Assert(m3 == xm2,"Permutation I/O check #2");
    Assert(p3 == xp2,"Permutation Compact I/O check #2");

#if XTEST == 0
    std::remove("tmvtest_permutation_io.dat");
#endif

#if 0
    //
    // Test other arithmetic
    // These will work in v0.90.  But not yet.
    //

    TestMatrixArith1(p3,p3,"Permutation");
    TestMatrixArith2(p3,p3,"Permutation");
    TestMatrixArith3(p3,p3,"Permutation");
    TestMatrixArith4(p3,p3,m2,cm2,"Permutation");
    TestMatrixArith4(m3,cm3,p2,p2,"Permutation");
    TestMatrixArith5(p3,p3,m2,cm2,"Permutation");
    TestMatrixArith5(m3,cm3,p2,p2,"Permutation");
    TestMatrixArith6x(p3,p3,m2,cm2,"Permutation");
    TestMatrixArith6x(m3,cm3,p2,p2,"Permutation");
#endif

    //
    // Test resize
    //

    tmv::Permutation q3(2);
    tmv::Permutation q3i(2);

    q3.resize(10);
    q3i.resize(10);
    q3 = p3;
    q3i = p3i;

    Assert(p3 == q3,"Permutation assignement after resize");
    Assert(p3i == q3i,"Permutation assignement after resize");

    std::cout<<"Permutation passed all tests\n";

}

#ifdef TEST_DOUBLE
template void TestPermutation<double>();
#endif
#ifdef TEST_FLOAT
template void TestPermutation<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestPermutation<long double>();
#endif
#ifdef TEST_INT
template void TestPermutation<int>();
#endif
