
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"
#include "tmv/TMV_Permutation.h"

#define NO_COMPLEX_ARITH
#define NODIV
#define NOASSIGN

#include "TMV_TestMatrixArith.h"

void TestPermutation()
{
    int pp1[10] = { 0,1,2,3,4,5,6,7,8,9 }; // identity
    int pp2[10] = { 9,8,7,6,5,5,6,7,8,9 }; // reversal
    int pp3[10] = { 5,3,5,9,4,5,9,8,9,9 }; // "random"

    int det1 = 1;
    int det2 = -1;
    int det3 = -1;

    tmv::Vector<int> v0(10);
    v0 << 0,1,2,3,4,5,6,7,8,9;
    tmv::Vector<int> v1 = v0;
    tmv::Vector<int> v2 = v0.reverse();
    tmv::Vector<int> v3l(10);
    v3l << 5,3,0,9,4,2,1,8,6,7;
    tmv::Vector<int> v3r(10);
    v3r << 2,6,5,1,4,0,8,9,7,3;

    tmv::Permutation p1(10,pp1,false,det1);
    tmv::Permutation p2(10,pp2,false,det2);
    tmv::Permutation p3(10,pp3,false,det3);

    tmv::Vector<int> v(10);
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

    tmv::Matrix<int> m1 = p1;
    tmv::Matrix<int> m2 = p2;
    tmv::Matrix<int> m3 = p3;

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

    v = v0;
    v *= p3;
    Assert(v == v3r,"Random permutation *=");
    v %= p3;
    Assert(v == v0,"Random permutation %=");
    v /= p3;
    Assert(v == v3r,"Random permutation /=");

    Assert(p1.det() == m1.det(),"Identity permutation determinant");
    Assert(p2.det() == m2.det(),"Reversal permutation determinant");
    Assert(p3.det() == m3.det(),"Random permutation determinant");

    TestMatrixArith1<int>(p3,p3,"Permutation");
    TestMatrixArith2<int>(p3,p3,"Permutation");
    TestMatrixArith3<int>(p3,p3,"Permutation");
    TestMatrixArith4<int>(p3,p3,p2,p2,"Permutation");
    TestMatrixArith5<int>(p3,p3,p2,p2,"Permutation");
    TestMatrixArith6x<int>(p3,p3,p2,p2,"Permutation");

    std::cout<<"Permutation passed all tests\n";
}

