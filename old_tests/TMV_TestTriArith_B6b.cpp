#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"

#define NOELEMMULT

#include "TMV_TestMatrixArith.h"

template <class T> void TestTriMatrixArith_B6b()
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

    TestMatrixArith6x(u1,cu1,a1,ca1,"UpperTri/Square 1");
    TestMatrixArith6x(u2,cu2,a1,ca1,"UpperTri/Square 2");
    TestMatrixArith6x(u4,cu4,a1,ca1,"UpperTri/Square 3");
    TestMatrixArith6x(u5,cu5,a1,ca1,"UpperTri/Square 4");
    TestMatrixArith6x(u1,cu1,a2,ca2,"UpperTri/Square 5");
#if (XTEST & 2)
    TestMatrixArith6x(u2,cu2,a2,ca2,"UpperTri/Square 6");
    TestMatrixArith6x(u4,cu4,a2,ca2,"UpperTri/Square 7");
    TestMatrixArith6x(u5,cu5,a2,ca2,"UpperTri/Square 8");
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
    TestMatrixArith6x(u3,cu3,a1,ca1,"UpperTri/Square 9");
    TestMatrixArith6x(u3,cu3,a2,ca2,"UpperTri/Square 10");
    TestMatrixArith6x(u3,cu3,a3,ca3,"UpperTri/Square 11");
    TestMatrixArith6x(u6,cu6,a1,ca1,"UpperTri/Square 12");
    TestMatrixArith6x(u6,cu6,a2,ca2,"UpperTri/Square 13");
    TestMatrixArith6x(u6,cu6,a3,ca3,"UpperTri/Square 14");
    TestMatrixArith6x(u1,cu1,a3,ca3,"UpperTri/Square 15");
    TestMatrixArith6x(u2,cu2,a3,ca3,"UpperTri/Square 16");
    TestMatrixArith6x(u4,cu4,a3,ca3,"UpperTri/Square 17");
    TestMatrixArith6x(u5,cu5,a3,ca3,"UpperTri/Square 18");
#endif

    tmv::LowerTriMatrixView<T> l1 = a1x.lowerTri();
    tmv::LowerTriMatrixView<CT> cl1 = ca1x.lowerTri();
    tmv::LowerTriMatrixView<T> l2 = a2x.lowerTri();
    tmv::LowerTriMatrixView<CT> cl2 = ca2x.lowerTri();
    tmv::LowerTriMatrixView<T> l4 = a1x.unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl4 = ca1x.unitLowerTri();
    tmv::LowerTriMatrixView<T> l5 = a2x.unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl5 = ca2x.unitLowerTri();

    TestMatrixArith6x(l1,cl1,a1,ca1,"LowerTri/Square 1");
    TestMatrixArith6x(l2,cl2,a1,ca1,"LowerTri/Square 2");
    TestMatrixArith6x(l4,cl4,a1,ca1,"LowerTri/Square 3");
    TestMatrixArith6x(l5,cl5,a1,ca1,"LowerTri/Square 4");
#if (XTEST & 2)
    TestMatrixArith6x(l1,cl1,a2,ca2,"LowerTri/Square 5");
    TestMatrixArith6x(l2,cl2,a2,ca2,"LowerTri/Square 6");
    TestMatrixArith6x(l4,cl4,a2,ca2,"LowerTri/Square 7");
    TestMatrixArith6x(l5,cl5,a2,ca2,"LowerTri/Square 8");
#endif
#if (XTEST & 1)
    tmv::LowerTriMatrixView<T> l3 = a3x.subMatrix(0,12,0,16,3,4).lowerTri();
    tmv::LowerTriMatrixView<CT> cl3 = ca3x.subMatrix(0,12,0,16,3,4).lowerTri();
    tmv::LowerTriMatrixView<T> l6 = a3x.subMatrix(0,12,0,16,3,4).unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl6 = ca3x.subMatrix(0,12,0,16,3,4).unitLowerTri();
    TestMatrixArith6x(l3,cl3,a1,ca1,"LowerTri/Square 9");
    TestMatrixArith6x(l3,cl3,a2,ca2,"LowerTri/Square 10");
    TestMatrixArith6x(l3,cl3,a3,ca3,"LowerTri/Square 11");
    TestMatrixArith6x(l6,cl6,a1,ca1,"LowerTri/Square 12");
    TestMatrixArith6x(l6,cl6,a2,ca2,"LowerTri/Square 13");
    TestMatrixArith6x(l6,cl6,a3,ca3,"LowerTri/Square 14");
    TestMatrixArith6x(l1,cl1,a3,ca3,"LowerTri/Square 15");
    TestMatrixArith6x(l2,cl2,a3,ca3,"LowerTri/Square 16");
    TestMatrixArith6x(l4,cl4,a3,ca3,"LowerTri/Square 17");
    TestMatrixArith6x(l5,cl5,a3,ca3,"LowerTri/Square 18");
#endif

    tmv::Matrix<T,tmv::RowMajor> a4x(4,6);
    for(int i=0;i<4;++i) for(int j=0;j<6;++j) {
        a4x(i,j) = T(2+4*i-5*j);
    }
    a4x.colRange(0,4) += a1x;

    tmv::Matrix<CT,tmv::RowMajor> ca4x = CT(1,2) * a4x;
    ca4x.colRange(0,4) += ca1x;

    tmv::Matrix<T,tmv::ColMajor> a5x = a4x;
    tmv::Matrix<CT,tmv::ColMajor> ca5x = ca4x;

    tmv::MatrixView<T> a4 = a4x.view();
    tmv::MatrixView<CT> ca4 = ca4x.view();
    tmv::MatrixView<T> a5 = a5x.view();
    tmv::MatrixView<CT> ca5 = ca5x.view();

    TestMatrixArith6x(u1,cu1,u4,cu4,"UpperTri/NonSquare 1");
    TestMatrixArith6x(u2,cu2,u4,cu4,"UpperTri/NonSquare 2");
    TestMatrixArith6x(u4,cu4,u4,cu4,"UpperTri/NonSquare 3");
    TestMatrixArith6x(u5,cu5,u4,cu4,"UpperTri/NonSquare 4");
#if (XTEST & 2)
    TestMatrixArith6x(u1,cu1,u5,cu5,"UpperTri/NonSquare 5");
    TestMatrixArith6x(u2,cu2,u5,cu5,"UpperTri/NonSquare 6");
    TestMatrixArith6x(u4,cu4,u5,cu5,"UpperTri/NonSquare 7");
    TestMatrixArith6x(u5,cu5,u5,cu5,"UpperTri/NonSquare 8");
#endif
#if (XTEST & 1)
    TestMatrixArith6x(u3,cu3,u4,cu4,"UpperTri/NonSquare 9");
    TestMatrixArith6x(u3,cu3,u5,cu5,"UpperTri/NonSquare 10");
    TestMatrixArith6x(u6,cu6,u4,cu4,"UpperTri/NonSquare 11");
    TestMatrixArith6x(u6,cu6,u5,cu5,"UpperTri/NonSquare 12");
#endif

    TestMatrixArith6x(l1,cl1,a4,ca4,"LowerTri/NonSquare 1");
    TestMatrixArith6x(l2,cl2,a4,ca4,"LowerTri/NonSquare 2");
    TestMatrixArith6x(l4,cl4,a4,ca4,"LowerTri/NonSquare 3");
    TestMatrixArith6x(l5,cl5,a4,ca4,"LowerTri/NonSquare 4");
    TestMatrixArith6x(l1,cl1,a5,ca5,"LowerTri/NonSquare 5");
#if (XTEST & 2)
    TestMatrixArith6x(l2,cl2,a5,ca5,"LowerTri/NonSquare 6");
    TestMatrixArith6x(l4,cl4,a5,ca5,"LowerTri/NonSquare 7");
    TestMatrixArith6x(l5,cl5,a5,ca5,"LowerTri/NonSquare 8");
#endif
#if (XTEST & 1)
    TestMatrixArith6x(l3,cl3,a4,ca4,"LowerTri/NonSquare 9");
    TestMatrixArith6x(l3,cl3,a5,ca5,"LowerTri/NonSquare 10");
    TestMatrixArith6x(l6,cl6,a4,ca4,"LowerTri/NonSquare 11");
    TestMatrixArith6x(l6,cl6,a5,ca5,"LowerTri/NonSquare 12");
#endif

#if (XTEST & 8)
    tmv::Matrix<T> a6x(4,0);
    tmv::Matrix<CT> ca6x(4,0);

    tmv::MatrixView<T> a6 = a6x.view();
    tmv::MatrixView<CT> ca6 = ca6x.view();

    TestMatrixArith6x(u1,cu1,a6,ca6,"UpperTri/Degenerate 1");
    TestMatrixArith6x(u2,cu2,a6,ca6,"UpperTri/Degenerate 2");
    TestMatrixArith6x(u4,cu4,a6,ca6,"UpperTri/Degenerate 3");
    TestMatrixArith6x(u5,cu5,a6,ca6,"UpperTri/Degenerate 4");
#if (XTEST & 1)
    TestMatrixArith6x(u3,cu3,a6,ca6,"UpperTri/Degenerate 5");
    TestMatrixArith6x(u6,cu6,a6,ca6,"UpperTri/Degenerate 6");
#endif

#if (XTEST & 2)
    TestMatrixArith6x(l1,cl1,a6,ca6,"LowerTri/Degenerate 1");
    TestMatrixArith6x(l2,cl2,a6,ca6,"LowerTri/Degenerate 2");
    TestMatrixArith6x(l4,cl4,a6,ca6,"LowerTri/Degenerate 3");
    TestMatrixArith6x(l5,cl5,a6,ca6,"LowerTri/Degenerate 4");
#endif
#if (XTEST & 1)
    TestMatrixArith6x(l4,cl4,a6,ca6,"LowerTri/Degenerate 5");
    TestMatrixArith6x(l6,cl6,a6,ca6,"LowerTri/Degenerate 6");
#endif
#endif
}


#ifdef TEST_DOUBLE
template void TestTriMatrixArith_B6b<double>();
#endif
#ifdef TEST_FLOAT
template void TestTriMatrixArith_B6b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestTriMatrixArith_B6b<long double>();
#endif
#ifdef TEST_INT
template void TestTriMatrixArith_B6b<int>();
#endif
