
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"

template <class M, class Ma, class Mb, class M1, class M1a, class M1b>
static void TestAliasMultUL2(
    const M& m, const Ma& ma, const Mb& mb,
    M1& m1, const M1a& m1a, const M1b& m1b, std::string label)
{
    if (showstartdone) {
        std::cout<<"Start AliasMultUL2\n";
        std::cout<<"m = "<<m<<std::endl;
        std::cout<<"ma = "<<ma<<std::endl;
        std::cout<<"mb = "<<mb<<std::endl;
    }

    typedef typename M::value_type T;
    typedef typename M::real_type RT;
    const int N = m.colsize();
    tmv::Matrix<T> m2(N,N);

    RT eps = EPS * RT(N);
    if (!std::numeric_limits<RT>::is_integer) eps *= Norm(m);

    m1 = m;
    m1 = m1a * m1b;
    m2 = ma * mb;
    Assert(Equal(m1,m2,eps),label+" m=a*b");

    m1 = m;
    m1 += m1a * m1b;
    m2 = m + ma * mb;
    Assert(Equal(m1,m2,eps),label+" m+=a*b");

    m1 = m;
    m1 = T(10) * m1a * m1b;
    m2 = T(10) * ma * mb;
    Assert(Equal(m1,m2,eps),label+" m=10*a*b");

    m1 = m;
    m1 += T(10) * m1a * m1b;
    m2 = m + T(10) * ma * mb;
    Assert(Equal(m1,m2,eps),label+" m+=10*a*b");

    m1 = m;
    m1 = m1a * m1b.transpose();
    m2 = ma * mb.transpose();
    Assert(Equal(m1,m2,eps),label+" m=a*bT");

    m1 = m;
    m1 += m1a * m1b.transpose();
    m2 = m + ma * mb.transpose();
    Assert(Equal(m1,m2,eps),label+" m+=a*bT");

    m1 = m;
    m1 = T(10) * m1a * m1b.transpose();
    m2 = T(10) * ma * mb.transpose();
    Assert(Equal(m1,m2,eps),label+" m=10*a*bT");

    m1 = m;
    m1 += T(10) * m1a * m1b.transpose();
    m2 = m + T(10) * ma * mb.transpose();
    Assert(Equal(m1,m2,eps),label+" m+=10*a*bT");

    m1 = m;
    m1 = m1a.transpose() * m1b;
    m2 = ma.transpose() * mb;
    Assert(Equal(m1,m2,eps),label+" m=aT*b");

    m1 = m;
    m1 += m1a.transpose() * m1b;
    m2 = m + ma.transpose() * mb;
    Assert(Equal(m1,m2,eps),label+" m+=aT*b");

    m1 = m;
    m1 = T(10) * m1a.transpose() * m1b;
    m2 = T(10) * ma.transpose() * mb;
    Assert(Equal(m1,m2,eps),label+" m=10*aT*b");

    m1 = m;
    m1 += T(10) * m1a.transpose() * m1b;
    m2 = m + T(10) * ma.transpose() * mb;
    Assert(Equal(m1,m2,eps),label+" m+=10*aT*b");

    m1 = m;
    m1 = m1a.transpose() * m1b.transpose();
    m2 = ma.transpose() * mb.transpose();
    Assert(Equal(m1,m2,eps),label+" m=aT*bT");

    m1 = m;
    m1 += m1a.transpose() * m1b.transpose();
    m2 = m + ma.transpose() * mb.transpose();
    Assert(Equal(m1,m2,eps),label+" m+=aT*bT");

    m1 = m;
    m1 = T(10) * m1a.transpose() * m1b.transpose();
    m2 = T(10) * ma.transpose() * mb.transpose();
    Assert(Equal(m1,m2,eps),label+" m=10*aT*bT");

    m1 = m;
    m1 += T(10) * m1a.transpose() * m1b.transpose();
    m2 = m + T(10) * ma.transpose() * mb.transpose();
    Assert(Equal(m1,m2,eps),label+" m+=10*aT*bT");
}

template <class M, class M1>
static void TestAliasMultUL1(const M& m, M1& m1, std::string label)
{
    TestAliasMultUL2(m,m.upperTri(),m.upperTri(),
                     m1,m1.upperTri(),m1.upperTri(),
                     label+" upperTri,upperTri");
    TestAliasMultUL2(m,m.upperTri(),m.lowerTri(),
                     m1,m1.upperTri(),m1.lowerTri(),
                     label+" upperTri,lowerTri");
    TestAliasMultUL2(m,m.lowerTri(),m.upperTri(),
                     m1,m1.lowerTri(),m1.upperTri(),
                     label+" lowerTri,upperTri");
    TestAliasMultUL2(m,m.lowerTri(),m.lowerTri(),
                     m1,m1.lowerTri(),m1.lowerTri(),
                     label+" lowerTri,lowerTri");

    TestAliasMultUL2(m,m.unitUpperTri(),m.upperTri(),
                     m1,m1.unitUpperTri(),m1.upperTri(),
                     label+" unitUpperTri,upperTri");
    TestAliasMultUL2(m,m.unitUpperTri(),m.lowerTri(),
                     m1,m1.unitUpperTri(),m1.lowerTri(),
                     label+" unitUpperTri,lowerTri");
    TestAliasMultUL2(m,m.unitLowerTri(),m.upperTri(),
                     m1,m1.unitLowerTri(),m1.upperTri(),
                     label+" unitLowerTri,upperTri");
    TestAliasMultUL2(m,m.unitLowerTri(),m.lowerTri(),
                     m1,m1.unitLowerTri(),m1.lowerTri(),
                     label+" unitLowerTri,lowerTri");

    TestAliasMultUL2(m,m.upperTri(),m.unitUpperTri(),
                     m1,m1.upperTri(),m1.unitUpperTri(),
                     label+" upperTri,unitUpperTri");
    TestAliasMultUL2(m,m.upperTri(),m.unitLowerTri(),
                     m1,m1.upperTri(),m1.unitLowerTri(),
                     label+" upperTri,unitLowerTri");
    TestAliasMultUL2(m,m.lowerTri(),m.unitUpperTri(),
                     m1,m1.lowerTri(),m1.unitUpperTri(),
                     label+" lowerTri,unitUpperTri");
    TestAliasMultUL2(m,m.lowerTri(),m.unitLowerTri(),
                     m1,m1.lowerTri(),m1.unitLowerTri(),
                     label+" lowerTri,unitLowerTri");

    TestAliasMultUL2(m,m.unitUpperTri(),m.unitUpperTri(),
                     m1,m1.unitUpperTri(),m1.unitUpperTri(),
                     label+" unitUpperTri,unitUpperTri");
    TestAliasMultUL2(m,m.unitUpperTri(),m.unitLowerTri(),
                     m1,m1.unitUpperTri(),m1.unitLowerTri(),
                     label+" unitUpperTri,unitLowerTri");
    TestAliasMultUL2(m,m.unitLowerTri(),m.unitUpperTri(),
                     m1,m1.unitLowerTri(),m1.unitUpperTri(),
                     label+" unitLowerTri,unitUpperTri");
    TestAliasMultUL2(m,m.unitLowerTri(),m.unitLowerTri(),
                     m1,m1.unitLowerTri(),m1.unitLowerTri(),
                     label+" unitLowerTri,unitLowerTri");
}

template <class T>
void TestAllAliasMultUL()
{
    // We make a special effort to do as many Tri*Tri calculations as 
    // possible in place when they are obvious aliases. 
    // e.g. m = m.upperTri() * m.lowerTri() is a common thing to want to do,
    // for example, as poart of the calculation of m.inverse() using the
    // LU decomposition of m.  (Invert each part in place, then combine.)
    // So this function tests all the combinations of upper and lower for
    // each term in the product and checks to make sure that the in-place
    // multiply works when we use it, and temporaries are used otherwise.
    
    tmv::Matrix<T> m(4,4);
    m << 1, 2, 3, 4,
         5, 6, 7, 8,
         9, 10, 11, 12,
         13, 14, 15, 16;

    tmv::Matrix<T,tmv::ColMajor> m1c = m;
    tmv::Matrix<T,tmv::RowMajor> m1r = m;

    TestAliasMultUL1(m,m1c,"alias CM");
    TestAliasMultUL1(m,m1r,"alias RM");
}

#ifdef TEST_DOUBLE
template void TestAllAliasMultUL<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllAliasMultUL<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllAliasMultUL<long double>();
#endif
#ifdef TEST_INT
template void TestAllAliasMultUL<int>();
#endif
