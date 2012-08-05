
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"
#include <fstream>

#include "TMV_TestMatrixArith.h"
#define CT std::complex<T>

template <class T> void TestMatrixArith_4()
{
    // Now we use the TestMatrixArith.h file to test lots of different
    // syntaxes for calling matrix arithmetic.  This tests the inline
    // parser more than the algorithms.
    tmv::Matrix<T,tmv::RowMajor> a1x(4,4);
    for(int i=0;i<4;++i) for(int j=0;j<4;++j) {
        a1x(i,j) = T(2+4*i-5*j);
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

    tmv::MatrixView<T> a1 = a1x.view();
    tmv::MatrixView<CT> ca1 = ca1x.view();
    tmv::MatrixView<T> a2 = a2x.view();
    tmv::MatrixView<CT> ca2 = ca2x.view();

    TestMatrixArith4(a1,ca1,a2,ca2,"Square 1");
    TestMatrixArith4(a2,ca2,a1,ca1,"Square 2");
#if (XTEST & 1)
    tmv::Matrix<T> a3x(12,16);
    for(int i=0;i<12;++i) for(int j=0;j<16;++j) a3x(i,j) = T(1-2*i+3*j);
    tmv::MatrixView<T> a3 = a3x.subMatrix(0,12,0,16,3,4);
    a3.diag().addToAll(30);
    tmv::Matrix<CT> ca3x = a3x*CT(1,-2);
    tmv::MatrixView<CT> ca3 = ca3x.subMatrix(0,12,0,16,3,4);
    ca3.diag().addToAll(CT(-22,15));

    TestMatrixArith4(a3,ca3,a1,ca1,"Square 3");
    TestMatrixArith4(a1,ca1,a3,ca3,"Square 4");
#if (XTEST & 2)
    TestMatrixArith4(a3,ca3,a2,ca2,"Square 5");
    TestMatrixArith4(a2,ca2,a3,ca3,"Square 6");
#endif
#endif

    tmv::Matrix<T,tmv::RowMajor> a4x(7,4);
    for(int i=0;i<7;++i) for(int j=0;j<4;++j) a4x(i,j) = T(1-3*i+2*j);
    tmv::Matrix<T,tmv::ColMajor> a5x = a4x.transpose();
    a4x.subMatrix(2,6,0,4) += a1x;
    a5x.subMatrix(0,4,1,5) -= a2x;

    tmv::Matrix<CT,tmv::RowMajor> ca4x = a4x*CT(1,2);
    tmv::Matrix<CT,tmv::ColMajor> ca5x = ca4x.adjoint();
    ca4x.subMatrix(2,6,0,4) += ca1x;
    ca5x.subMatrix(0,4,1,5) -= ca2x;
    ca4x.col(1) *= CT(2,1);
    ca4x.row(6).addToAll(CT(-7,2));
    ca5x.col(3) *= CT(-1,3);
    ca5x.row(0).addToAll(CT(1,9));

    tmv::MatrixView<T> a4 = a4x.view();
    tmv::MatrixView<CT> ca4 = ca4x.view();
    tmv::MatrixView<T> a5 = a5x.view();
    tmv::MatrixView<CT> ca5 = ca5x.view();

    TestMatrixArith4(a1,ca1,a4,ca4,"NonSquare 1");
    TestMatrixArith4(a4,ca4,a1,ca1,"NonSquare 2");
    TestMatrixArith4(a4,ca4,a5,ca5,"NonSquare 3");
    TestMatrixArith4(a5,ca5,a4,ca4,"NonSquare 4");
#if (XTEST & 2)
    TestMatrixArith4(a5,ca5,a1,ca1,"NonSquare 5");
    TestMatrixArith4(a1,ca2,a5,ca5,"NonSquare 6");
    TestMatrixArith4(a2,ca2,a4,ca4,"NonSquare 7");
    TestMatrixArith4(a2,ca2,a5,ca5,"NonSquare 8");
    TestMatrixArith4(a4,ca4,a2,ca2,"NonSquare 9");
    TestMatrixArith4(a5,ca5,a2,ca2,"NonSquare 10");
#if (XTEST & 1)
    TestMatrixArith4(a3,ca3,a4,ca4,"NonSquare 11");
    TestMatrixArith4(a3,ca3,a5,ca4,"NonSquare 12");
    TestMatrixArith4(a4,ca4,a3,ca3,"NonSquare 13");
    TestMatrixArith4(a5,ca5,a3,ca3,"NonSquare 14");
#endif
#endif

#if (XTEST & 8)
    tmv::Matrix<T> a6x(4,0,1);
    tmv::Matrix<T> a7x(0,4,1);
    tmv::Matrix<CT> ca6x = a6x;
    tmv::Matrix<CT> ca7x = a7x;

    tmv::MatrixView<T> a6 = a6x.view();
    tmv::MatrixView<CT> ca6 = ca6x.view();
    tmv::MatrixView<T> a7 = a7x.view();
    tmv::MatrixView<CT> ca7 = ca7x.view();

    TestMatrixArith4(a1,ca1,a6,ca6,"Degenerate 1");
    TestMatrixArith4(a2,ca2,a6,ca6,"Degenerate 2");
    TestMatrixArith4(a3,ca3,a6,ca6,"Degenerate 3");
    TestMatrixArith4(a4,ca4,a6,ca6,"Degenerate 4");
    TestMatrixArith4(a5,ca5,a6,ca6,"Degenerate 5");
    TestMatrixArith4(a1,ca1,a7,ca7,"Degenerate 6");
    TestMatrixArith4(a2,ca2,a7,ca7,"Degenerate 7");
    TestMatrixArith4(a3,ca3,a7,ca7,"Degenerate 8");
    TestMatrixArith4(a4,ca4,a7,ca7,"Degenerate 9");
    TestMatrixArith4(a5,ca5,a7,ca7,"Degenerate 10");
    TestMatrixArith4(a6,ca6,a1,ca1,"Degenerate 11");
    TestMatrixArith4(a6,ca6,a2,ca2,"Degenerate 12");
    TestMatrixArith4(a6,ca6,a3,ca3,"Degenerate 13");
    TestMatrixArith4(a6,ca6,a4,ca4,"Degenerate 14");
    TestMatrixArith4(a6,ca6,a5,ca5,"Degenerate 15");
    TestMatrixArith4(a7,ca7,a1,ca1,"Degenerate 16");
    TestMatrixArith4(a7,ca7,a2,ca2,"Degenerate 17");
    TestMatrixArith4(a7,ca7,a3,ca3,"Degenerate 18");
    TestMatrixArith4(a7,ca7,a4,ca4,"Degenerate 19");
    TestMatrixArith4(a7,ca7,a5,ca5,"Degenerate 20");
    TestMatrixArith4(a7,ca7,a6,ca6,"Degenerate 21");
    TestMatrixArith4(a6,ca6,a7,ca7,"Degenerate 22");
#endif
}

#ifdef TEST_DOUBLE
template void TestMatrixArith_4<double>();
#endif
#ifdef TEST_FLOAT
template void TestMatrixArith_4<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestMatrixArith_4<long double>();
#endif
#ifdef TEST_INT
template void TestMatrixArith_4<int>();
#endif
