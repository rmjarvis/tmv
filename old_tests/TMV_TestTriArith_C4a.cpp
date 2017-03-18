#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"

#define NOADDEQ
#define NOELEMMULT

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestTriMatrixArith_C4a()
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

    tmv::DiagMatrixView<T> d1 = DiagMatrixViewOf(a1x.row(0));
    tmv::DiagMatrixView<CT> cd1 = DiagMatrixViewOf(ca1x.row(0));

    TestMatrixArith4(d1,cd1,u1,cu1,"Diag/UpperTri 1");
    TestMatrixArith4(d1,cd1,u2,cu2,"Diag/UpperTri 2");
    TestMatrixArith4(d1,cd1,u4,cu4,"Diag/UpperTri 3");
    TestMatrixArith4(d1,cd1,u5,cu5,"Diag/UpperTri 4");
#if (XTEST & 1)
    tmv::Matrix<T> a3x(12,16);
    for(int i=0;i<12;++i) for(int j=0;j<16;++j) a3x(i,j) = T(1-2*i+3*j);
    a3x.diag().addToAll(30);
    tmv::Matrix<CT> ca3x = a3x*CT(1,-2);
    ca3x.diag().addToAll(CT(-22,15));

    tmv::DiagMatrixView<T> d3 = DiagMatrixViewOf(a1x.diag());
    tmv::DiagMatrixView<CT> cd3 = DiagMatrixViewOf(ca1x.diag());
    tmv::UpperTriMatrixView<T> u3 = a3x.subMatrix(0,12,0,16,3,4).upperTri();
    tmv::UpperTriMatrixView<CT> cu3 = ca3x.subMatrix(0,12,0,16,3,4).upperTri();
    tmv::UpperTriMatrixView<T> u6 = a3x.subMatrix(0,12,0,16,3,4).unitUpperTri();
    tmv::UpperTriMatrixView<CT> cu6 = ca3x.subMatrix(0,12,0,16,3,4).unitUpperTri();
    TestMatrixArith4(d3,cd3,u1,cu1,"Diag/UpperTri 5");
    TestMatrixArith4(d3,cd3,u2,cu2,"Diag/UpperTri 6");
    TestMatrixArith4(d3,cd3,u3,cu3,"Diag/UpperTri 7");
    TestMatrixArith4(d3,cd3,u4,cu4,"Diag/UpperTri 8");
    TestMatrixArith4(d3,cd3,u5,cu5,"Diag/UpperTri 9");
    TestMatrixArith4(d3,cd3,u6,cu6,"Diag/UpperTri 10");
    TestMatrixArith4(d1,cd1,u3,cu3,"Diag/UpperTri 11");
    TestMatrixArith4(d1,cd1,u6,cu6,"Diag/UpperTri 12");
#endif

#if (XTEST & 2)
    tmv::LowerTriMatrixView<T> l1 = a1x.lowerTri();
    tmv::LowerTriMatrixView<CT> cl1 = ca1x.lowerTri();
    tmv::LowerTriMatrixView<T> l2 = a2x.lowerTri();
    tmv::LowerTriMatrixView<CT> cl2 = ca2x.lowerTri();
    tmv::LowerTriMatrixView<T> l4 = a1x.unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl4 = ca1x.unitLowerTri();
    tmv::LowerTriMatrixView<T> l5 = a2x.unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl5 = ca2x.unitLowerTri();

    TestMatrixArith4(d1,cd1,l1,cl1,"Diag/LowerTri 1");
    TestMatrixArith4(d1,cd1,l2,cl2,"Diag/LowerTri 2");
    TestMatrixArith4(d1,cd1,l4,cl4,"Diag/LowerTri 3");
    TestMatrixArith4(d1,cd1,l5,cl5,"Diag/LowerTri 4");
#if (XTEST & 1)
    tmv::LowerTriMatrixView<T> l3 = a3x.subMatrix(0,12,0,16,3,4).lowerTri();
    tmv::LowerTriMatrixView<CT> cl3 = ca3x.subMatrix(0,12,0,16,3,4).lowerTri();
    tmv::LowerTriMatrixView<T> l6 = a3x.subMatrix(0,12,0,16,3,4).unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl6 = ca3x.subMatrix(0,12,0,16,3,4).unitLowerTri();
    TestMatrixArith4(d3,cd3,l1,cl1,"Diag/LowerTri 5");
    TestMatrixArith4(d3,cd3,l2,cl2,"Diag/LowerTri 6");
    TestMatrixArith4(d3,cd3,l3,cl3,"Diag/LowerTri 7");
    TestMatrixArith4(d3,cd3,l4,cl4,"Diag/LowerTri 8");
    TestMatrixArith4(d3,cd3,l5,cl5,"Diag/LowerTri 9");
    TestMatrixArith4(d3,cd3,l6,cl6,"Diag/LowerTri 10");
    TestMatrixArith4(d1,cd1,l3,cl3,"Diag/LowerTri 11");
    TestMatrixArith4(d1,cd1,l6,cl6,"Diag/LowerTri 12");
#endif
#endif
}

#ifdef TEST_DOUBLE
template void TestTriMatrixArith_C4a<double>();
#endif
#ifdef TEST_FLOAT
template void TestTriMatrixArith_C4a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestTriMatrixArith_C4a<long double>();
#endif
#ifdef TEST_INT
template void TestTriMatrixArith_C4a<int>();
#endif
