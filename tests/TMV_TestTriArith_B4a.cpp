#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"

#define NOELEMMULT

#include "TMV_TestMatrixArith.h"

template <class T> void TestTriMatrixArith_B4a()
{
    typedef std::complex<T> CT;

    tmv::Matrix<T,tmv::RowMajor> a1x(4,4);
    for(int i=0;i<4;++i) for(int j=0;j<4;++j) {
        a1x(i,j) = T(3+4*i-6*j);
    }
    a1x(0,0) = 14; 
    a1x(1,0) = -2; 
    a1x(2,0) = 7; 
    a1x(3,0) = -10;
    a1x(2,2) = 30;

    tmv::Matrix<CT,tmv::RowMajor> ca1x = a1x;
    ca1x(2,3) += CT(2,3);
    ca1x(1,0) *= CT(0,2);
    ca1x.col(1) *= CT(-1,3);
    ca1x.row(3) += tmv::Vector<CT>(4,CT(1,9));

    tmv::Matrix<T,tmv::ColMajor> a2x = a1x.transpose();
    a2x.row(1) *= T(3);
    a2x.col(2) -= tmv::Vector<T>(4,4);
    tmv::Matrix<CT,tmv::ColMajor> ca2x = ca1x;
    ca2x -= a2x;
    ca2x *= CT(1,-2);
    ca2x(0,0) = CT(0,-5);

    tmv::UpperTriMatrixView<T> u1 = a1x.upperTri();
    tmv::UpperTriMatrixView<CT> cu1 = ca1x.upperTri();
    tmv::UpperTriMatrixView<T> u2 = a2x.upperTri();
    tmv::UpperTriMatrixView<CT> cu2 = ca2x.upperTri();
    tmv::UpperTriMatrixView<T> u4 = a1x.unitUpperTri();
    tmv::UpperTriMatrixView<CT> cu4 = ca1x.unitUpperTri();
    tmv::UpperTriMatrixView<T> u5 = a2x.unitUpperTri();
    tmv::UpperTriMatrixView<CT> cu5 = ca2x.unitUpperTri();

    tmv::MatrixView<T> a1 = a1x.view();
    tmv::MatrixView<CT> ca1 = ca1x.view();
    tmv::MatrixView<T> a2 = a2x.view();
    tmv::MatrixView<CT> ca2 = ca2x.view();

    TestMatrixArith4(a1,ca1,u1,cu1,"Square/UpperTri 1");
    TestMatrixArith4(a1,ca1,u2,cu2,"Square/UpperTri 2");
    TestMatrixArith4(a1,ca1,u4,cu4,"Square/UpperTri 3");
    TestMatrixArith4(a1,ca1,u5,cu5,"Square/UpperTri 4");
    TestMatrixArith4(a2,ca2,u1,cu1,"Square/UpperTri 5");
#if (XTEST & 2)
    TestMatrixArith4(a2,ca2,u2,cu2,"Square/UpperTri 6");
    TestMatrixArith4(a2,ca2,u4,cu4,"Square/UpperTri 7");
    TestMatrixArith4(a2,ca2,u5,cu5,"Square/UpperTri 8");
#endif
#if (XTEST & 1)
    tmv::Matrix<T> a3x(12,16);
    for(int i=0;i<12;++i) for(int j=0;j<16;++j) a3x(i,j) = T(1-2*i+3*j);
    a3x.diag().addToAll(30);
    tmv::Matrix<CT> ca3x = a3x*CT(1,-2);
    ca3x.diag().addToAll(CT(-22,15));

    tmv::UpperTriMatrixView<T> u3 = a3x.subMatrix(0,12,0,16,3,4).upperTri();
    tmv::UpperTriMatrixView<CT> cu3 = ca3x.subMatrix(0,12,0,16,3,4).upperTri();
    tmv::UpperTriMatrixView<T> u6 = a3x.subMatrix(0,12,0,16,3,4).unitUpperTri();
    tmv::UpperTriMatrixView<CT> cu6 = ca3x.subMatrix(0,12,0,16,3,4).unitUpperTri();
    tmv::MatrixView<T> a3 = a3x.view();
    tmv::MatrixView<CT> ca3 = ca3x.view();
    TestMatrixArith4(a3,ca3,u1,cu1,"Square/UpperTri 9");
    TestMatrixArith4(a3,ca3,u2,cu2,"Square/UpperTri 10");
    TestMatrixArith4(a3,ca3,u3,cu3,"Square/UpperTri 11");
    TestMatrixArith4(a3,ca3,u4,cu4,"Square/UpperTri 12");
    TestMatrixArith4(a3,ca3,u5,cu5,"Square/UpperTri 13");
    TestMatrixArith4(a3,ca3,u6,cu6,"Square/UpperTri 14");
    TestMatrixArith4(a1,ca1,u3,cu3,"Square/UpperTri 15");
    TestMatrixArith4(a2,ca2,u3,cu3,"Square/UpperTri 16");
    TestMatrixArith4(a1,ca1,u6,cu6,"Square/UpperTri 17");
    TestMatrixArith4(a2,ca2,u6,cu6,"Square/UpperTri 18");
#endif

    tmv::LowerTriMatrixView<T> l1 = a1x.lowerTri();
    tmv::LowerTriMatrixView<CT> cl1 = ca1x.lowerTri();
    tmv::LowerTriMatrixView<T> l2 = a2x.lowerTri();
    tmv::LowerTriMatrixView<CT> cl2 = ca2x.lowerTri();
    tmv::LowerTriMatrixView<T> l4 = a1x.unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl4 = ca1x.unitLowerTri();
    tmv::LowerTriMatrixView<T> l5 = a2x.unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl5 = ca2x.unitLowerTri();

    TestMatrixArith4(a1,ca1,l1,cl1,"Square/LowerTri 1");
    TestMatrixArith4(a1,ca1,l2,cl2,"Square/LowerTri 2");
    TestMatrixArith4(a1,ca1,l4,cl4,"Square/LowerTri 3");
    TestMatrixArith4(a1,ca1,l5,cl5,"Square/LowerTri 4");
#if (XTEST & 2)
    TestMatrixArith4(a2,ca2,l1,cl1,"Square/LowerTri 5");
    TestMatrixArith4(a2,ca2,l2,cl2,"Square/LowerTri 6");
    TestMatrixArith4(a2,ca2,l4,cl4,"Square/LowerTri 7");
    TestMatrixArith4(a2,ca2,l5,cl5,"Square/LowerTri 8");
#endif
#if (XTEST & 1)
    tmv::LowerTriMatrixView<T> l3 = a3x.subMatrix(0,12,0,16,3,4).lowerTri();
    tmv::LowerTriMatrixView<CT> cl3 = ca3x.subMatrix(0,12,0,16,3,4).lowerTri();
    tmv::LowerTriMatrixView<T> l6 = a3x.subMatrix(0,12,0,16,3,4).unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl6 = ca3x.subMatrix(0,12,0,16,3,4).unitLowerTri();
    TestMatrixArith4(a3,ca3,l1,cl1,"Square/LowerTri 9");
    TestMatrixArith4(a3,ca3,l2,cl2,"Square/LowerTri 10");
    TestMatrixArith4(a3,ca3,l3,cl3,"Square/LowerTri 11");
    TestMatrixArith4(a3,ca3,l4,cl4,"Square/LowerTri 12");
    TestMatrixArith4(a3,ca3,l5,cl5,"Square/LowerTri 13");
    TestMatrixArith4(a3,ca3,l6,cl6,"Square/LowerTri 14");
    TestMatrixArith4(a1,ca1,l3,cl3,"Square/LowerTri 15");
    TestMatrixArith4(a2,ca2,l3,cl3,"Square/LowerTri 16");
    TestMatrixArith4(a1,ca1,l6,cl6,"Square/LowerTri 17");
    TestMatrixArith4(a2,ca2,l6,cl6,"Square/LowerTri 18");
#endif
}

#ifdef TEST_DOUBLE
template void TestTriMatrixArith_B4a<double>();
#endif
#ifdef TEST_FLOAT
template void TestTriMatrixArith_B4a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestTriMatrixArith_B4a<long double>();
#endif
#ifdef TEST_INT
template void TestTriMatrixArith_B4a<int>();
#endif
