#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"

#define NOADDEQ
#define NOELEMMULT

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestTriMatrixArith_A4b()
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

    tmv::LowerTriMatrixView<T> l1 = a1x.lowerTri();
    tmv::LowerTriMatrixView<CT> cl1 = ca1x.lowerTri();
    tmv::LowerTriMatrixView<T> l2 = a2x.lowerTri();
    tmv::LowerTriMatrixView<CT> cl2 = ca2x.lowerTri();

    TestMatrixArith4(u1,cu1,l1,cl1,"UpperTri/LowerTri 1");
    TestMatrixArith4(u2,cu2,l2,cl2,"UpperTri/LowerTri 2");
#if (XTEST & 2)
    TestMatrixArith4(u1,cu1,l2,cl2,"UpperTri/LowerTri 3");
    TestMatrixArith4(u2,cu2,l1,cl1,"UpperTri/LowerTri 4");
#endif
#if (XTEST & 1)
    tmv::Matrix<T> a3x(12,16);
    for(int i=0;i<12;++i) for(int j=0;j<16;++j) a3x(i,j) = T(1-2*i+3*j);
    a3x.diag().addToAll(30);
    tmv::Matrix<CT> ca3x = a3x*CT(1,-2);
    ca3x.diag().addToAll(CT(-22,15));

    tmv::UpperTriMatrixView<T> u3 = a3x.subMatrix(0,12,0,16,3,4).upperTri();
    tmv::UpperTriMatrixView<CT> cu3 = ca3x.subMatrix(0,12,0,16,3,4).upperTri();
    tmv::LowerTriMatrixView<T> l3 = a3x.subMatrix(0,12,0,16,3,4).lowerTri();
    tmv::LowerTriMatrixView<CT> cl3 = ca3x.subMatrix(0,12,0,16,3,4).lowerTri();
    TestMatrixArith4(u3,cu3,l1,cl1,"UpperTri/LowerTri 5");
    TestMatrixArith4(u3,cu3,l2,cl2,"UpperTri/LowerTri 6");
    TestMatrixArith4(u1,cu1,l3,cl3,"UpperTri/LowerTri 7");
    TestMatrixArith4(u2,cu2,l3,cl3,"UpperTri/LowerTri 8");
#endif

    TestMatrixArith4(l1,cl1,u1,cu1,"LowerTri/UpperTri 1");
    TestMatrixArith4(l2,cl2,u2,cu2,"LowerTri/UpperTri 2");
#if (XTEST & 2)
    TestMatrixArith4(l1,cl1,u2,cu2,"LowerTri/UpperTri 3");
    TestMatrixArith4(l2,cl2,u1,cu1,"LowerTri/UpperTri 4");
#endif
#if (XTEST & 1)
    TestMatrixArith4(l3,cl3,u1,cu1,"LowerTri/UpperTri 5");
    TestMatrixArith4(l3,cl3,u2,cu2,"LowerTri/UpperTri 6");
    TestMatrixArith4(l1,cl1,u3,cu3,"LowerTri/UpperTri 7");
    TestMatrixArith4(l2,cl2,u3,cu3,"LowerTri/UpperTri 8");
#endif

#if (XTEST & 2)
    tmv::UpperTriMatrixView<T> u4 = a1x.unitUpperTri();
    tmv::UpperTriMatrixView<CT> cu4 = ca1x.unitUpperTri();
    tmv::UpperTriMatrixView<T> u5 = a2x.unitUpperTri();
    tmv::UpperTriMatrixView<CT> cu5 = ca2x.unitUpperTri();

    TestMatrixArith4(u4,cu4,u4,cu4,"UpperTri 18");
    TestMatrixArith4(u5,cu5,u4,cu4,"UpperTri 19");
    TestMatrixArith4(u4,cu4,u5,cu5,"UpperTri 20");
    TestMatrixArith4(u5,cu5,u5,cu5,"UpperTri 21");
#if (XTEST & 1)
    tmv::UpperTriMatrixView<T> u6 = a3x.subMatrix(0,12,0,16,3,4).unitUpperTri();
    tmv::UpperTriMatrixView<CT> cu6 = ca3x.subMatrix(0,12,0,16,3,4).unitUpperTri();
    TestMatrixArith4(u6,cu6,u4,cu4,"UpperTri 22");
    TestMatrixArith4(u6,cu6,u5,cu5,"UpperTri 23");
    TestMatrixArith4(u6,cu6,u6,cu6,"UpperTri 24");
    TestMatrixArith4(u4,cu4,u6,cu6,"UpperTri 25");
    TestMatrixArith4(u5,cu5,u6,cu6,"UpperTri 26");
#endif

    tmv::LowerTriMatrixView<T> l4 = a1x.unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl4 = ca1x.unitLowerTri();
    tmv::LowerTriMatrixView<T> l5 = a2x.unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl5 = ca2x.unitLowerTri();

    TestMatrixArith4(l4,cl4,l4,cl4,"LowerTri 18");
    TestMatrixArith4(l5,cl5,l4,cl4,"LowerTri 19");
    TestMatrixArith4(l4,cl4,l5,cl5,"LowerTri 20");
    TestMatrixArith4(l5,cl5,l5,cl5,"LowerTri 21");
#if (XTEST & 1)
    tmv::LowerTriMatrixView<T> l6 = a3x.subMatrix(0,12,0,16,3,4).unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl6 = ca3x.subMatrix(0,12,0,16,3,4).unitLowerTri();
    TestMatrixArith4(l6,cl6,l4,cl4,"LowerTri 22");
    TestMatrixArith4(l6,cl6,l5,cl5,"LowerTri 23");
    TestMatrixArith4(l6,cl6,l6,cl6,"LowerTri 24");
    TestMatrixArith4(l4,cl4,l6,cl6,"LowerTri 25");
    TestMatrixArith4(l5,cl5,l6,cl6,"LowerTri 26");
#endif

    TestMatrixArith4(u1,cu1,l4,cl4,"UpperTri/LowerTri 9");
    TestMatrixArith4(u2,cu2,l5,cl5,"UpperTri/LowerTri 10");
#if (XTEST & 2)
    TestMatrixArith4(u1,cu1,l5,cl5,"UpperTri/LowerTri 11");
    TestMatrixArith4(u2,cu2,l4,cl4,"UpperTri/LowerTri 12");
#endif
#if (XTEST & 1)
    TestMatrixArith4(u3,cu3,l4,cl4,"UpperTri/LowerTri 13");
    TestMatrixArith4(u3,cu3,l5,cl5,"UpperTri/LowerTri 14");
    TestMatrixArith4(u1,cu1,l6,cl6,"UpperTri/LowerTri 15");
    TestMatrixArith4(u2,cu2,l6,cl6,"UpperTri/LowerTri 16");
#endif

    TestMatrixArith4(u4,cu4,l1,cl1,"UpperTri/LowerTri 17");
    TestMatrixArith4(u5,cu5,l2,cl2,"UpperTri/LowerTri 18");
#if (XTEST & 2)
    TestMatrixArith4(u4,cu4,l2,cl2,"UpperTri/LowerTri 19");
    TestMatrixArith4(u5,cu5,l1,cl1,"UpperTri/LowerTri 20");
#endif
#if (XTEST & 1)
    TestMatrixArith4(u6,cu6,l1,cl1,"UpperTri/LowerTri 21");
    TestMatrixArith4(u6,cu6,l2,cl2,"UpperTri/LowerTri 22");
    TestMatrixArith4(u4,cu4,l3,cl3,"UpperTri/LowerTri 23");
    TestMatrixArith4(u5,cu5,l3,cl3,"UpperTri/LowerTri 24");
#endif

    TestMatrixArith4(u4,cu4,l4,cl4,"UpperTri/LowerTri 25");
    TestMatrixArith4(u5,cu5,l5,cl5,"UpperTri/LowerTri 26");
#if (XTEST & 2)
    TestMatrixArith4(u4,cu4,l5,cl5,"UpperTri/LowerTri 27");
    TestMatrixArith4(u5,cu5,l4,cl4,"UpperTri/LowerTri 28");
#endif
#if (XTEST & 1)
    TestMatrixArith4(u6,cu6,l4,cl4,"UpperTri/LowerTri 29");
    TestMatrixArith4(u6,cu6,l5,cl5,"UpperTri/LowerTri 30");
    TestMatrixArith4(u4,cu4,l6,cl6,"UpperTri/LowerTri 31");
    TestMatrixArith4(u5,cu5,l6,cl6,"UpperTri/LowerTri 32");
#endif

    TestMatrixArith4(l1,cl1,u4,cu4,"LowerTri/UpperTri 9");
    TestMatrixArith4(l2,cl2,u5,cu5,"LowerTri/UpperTri 10");
#if (XTEST & 2)
    TestMatrixArith4(l1,cl1,u5,cu5,"LowerTri/UpperTri 11");
    TestMatrixArith4(l2,cl2,u4,cu4,"LowerTri/UpperTri 12");
#endif
#if (XTEST & 1)
    TestMatrixArith4(l3,cl3,u4,cu4,"LowerTri/UpperTri 13");
    TestMatrixArith4(l3,cl3,u5,cu5,"LowerTri/UpperTri 14");
    TestMatrixArith4(l1,cl1,u6,cu6,"LowerTri/UpperTri 15");
    TestMatrixArith4(l2,cl2,u6,cu6,"LowerTri/UpperTri 16");
#endif
    TestMatrixArith4(l4,cl4,u1,cu1,"LowerTri/UpperTri 17");
    TestMatrixArith4(l5,cl5,u2,cu2,"LowerTri/UpperTri 18");
#if (XTEST & 2)
    TestMatrixArith4(l4,cl4,u2,cu2,"LowerTri/UpperTri 19");
    TestMatrixArith4(l5,cl5,u1,cu1,"LowerTri/UpperTri 20");
#endif
#if (XTEST & 1)
    TestMatrixArith4(l6,cl6,u1,cu1,"LowerTri/UpperTri 21");
    TestMatrixArith4(l6,cl6,u2,cu2,"LowerTri/UpperTri 22");
    TestMatrixArith4(l4,cl4,u3,cu3,"LowerTri/UpperTri 23");
    TestMatrixArith4(l5,cl5,u3,cu3,"LowerTri/UpperTri 24");
#endif

    TestMatrixArith4(l4,cl4,u4,cu4,"LowerTri/UpperTri 25");
    TestMatrixArith4(l5,cl5,u5,cu5,"LowerTri/UpperTri 26");
#if (XTEST & 2)
    TestMatrixArith4(l4,cl4,u5,cu5,"LowerTri/UpperTri 27");
    TestMatrixArith4(l5,cl5,u4,cu4,"LowerTri/UpperTri 28");
#endif
#if (XTEST & 1)
    TestMatrixArith4(l6,cl6,u4,cu4,"LowerTri/UpperTri 29");
    TestMatrixArith4(l6,cl6,u5,cu5,"LowerTri/UpperTri 30");
    TestMatrixArith4(l4,cl4,u6,cu6,"LowerTri/UpperTri 31");
    TestMatrixArith4(l5,cl5,u6,cu6,"LowerTri/UpperTri 32");
#endif

#endif
}

#ifdef TEST_DOUBLE
template void TestTriMatrixArith_A4b<double>();
#endif
#ifdef TEST_FLOAT
template void TestTriMatrixArith_A4b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestTriMatrixArith_A4b<long double>();
#endif
#ifdef TEST_INT
template void TestTriMatrixArith_A4b<int>();
#endif
