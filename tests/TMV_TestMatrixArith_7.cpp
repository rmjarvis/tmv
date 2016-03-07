
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"
#include <fstream>

#include "TMV_TestMatrixArith.h"
#define CT std::complex<T>

template <class T> void TestMatrixArith_7()
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

    tmv::Vector<T> v1 = a1.col(0);
    tmv::VectorView<T> v1v = v1.view();
    tmv::Vector<T> v15(20);
    tmv::VectorView<T> v1s = v15.subVector(0,20,5);
    v1s = v1v;

    tmv::Vector<T> v2 = a1.row(2);
    tmv::VectorView<T> v2v = v2.view();
    tmv::Vector<T> v25(20);
    tmv::VectorView<T> v2s = v25.subVector(0,20,5);
    v2s = v2v;

    tmv::Vector<CT> cv1 = ca1.col(0);
    tmv::VectorView<CT> cv1v = cv1.view();
    tmv::Vector<CT> cv15(20);
    tmv::VectorView<CT> cv1s = cv15.subVector(0,20,5);
    cv1s = cv1v;

    tmv::Vector<CT> cv2 = ca1.row(2);
    tmv::VectorView<CT> cv2v = cv2.view();
    tmv::Vector<CT> cv25(20);
    tmv::VectorView<CT> cv2s = cv25.subVector(0,20,5);
    cv2s = cv2v;

    TestMatrixArith7(a1,ca1,v1v,cv1v,v2v,cv2v,"Square 1");
    TestMatrixArith7(a1,ca1,v1s,cv1s,v2v,cv2v,"Square 2");
    TestMatrixArith7(a1,ca1,v1v,cv1v,v2s,cv2s,"Square 3");
    TestMatrixArith7(a1,ca1,v1s,cv1s,v2s,cv2s,"Square 4");
    TestMatrixArith7(a2,ca2,v1v,cv1v,v2v,cv2v,"Square 5");
    TestMatrixArith7(a2,ca2,v1s,cv1s,v2v,cv2v,"Square 6");
    TestMatrixArith7(a2,ca2,v1v,cv1v,v2s,cv2s,"Square 7");
    TestMatrixArith7(a2,ca2,v1s,cv1s,v2s,cv2s,"Square 8");
#if (XTEST & 1)
    tmv::Matrix<T> a3x(12,16);
    for(int i=0;i<12;++i) for(int j=0;j<16;++j) a3x(i,j) = T(1-2*i+3*j);
    tmv::MatrixView<T> a3 = a3x.subMatrix(0,12,0,16,3,4);
    a3.diag().addToAll(30);
    tmv::Matrix<CT> ca3x = a3x*CT(1,-2);
    tmv::MatrixView<CT> ca3 = ca3x.subMatrix(0,12,0,16,3,4);
    ca3.diag().addToAll(CT(-22,15));

    TestMatrixArith7(a3,ca3,v1v,cv1v,v2v,cv2v,"Square 9");
    TestMatrixArith7(a3,ca3,v1s,cv1s,v2v,cv2v,"Square 10");
    TestMatrixArith7(a3,ca3,v1v,cv1v,v2s,cv2s,"Square 11");
    TestMatrixArith7(a3,ca3,v1s,cv1s,v2s,cv2s,"Square 12");
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

    tmv::Vector<T> v3 = a4.col(2);
    tmv::VectorView<T> v3v = v3.view();
    tmv::Vector<T> v35(35);
    tmv::VectorView<T> v3s = v35.subVector(0,35,5);
    v3s = v3v;

    tmv::Vector<CT> cv3 = ca4.col(2);
    tmv::VectorView<CT> cv3v = cv3.view();
    tmv::Vector<CT> cv35(35);
    tmv::VectorView<CT> cv3s = cv35.subVector(0,35,5);
    cv3s = cv3v;

    TestMatrixArith7(a4,ca4,v3v,cv3v,v2v,cv2v,"NonSquare 1");
    TestMatrixArith7(a4,ca4,v3s,cv3s,v2v,cv2v,"NonSquare 2");
    TestMatrixArith7(a4,ca4,v3v,cv3v,v2s,cv2s,"NonSquare 3");
    TestMatrixArith7(a4,ca4,v3s,cv3s,v2s,cv2s,"NonSquare 4");
    TestMatrixArith7(a5,ca5,v1v,cv1v,v3v,cv3v,"NonSquare 5");
    TestMatrixArith7(a5,ca5,v1s,cv1s,v3v,cv3v,"NonSquare 6");
    TestMatrixArith7(a5,ca5,v1v,cv1v,v3s,cv3s,"NonSquare 7");
    TestMatrixArith7(a5,ca5,v1s,cv1s,v3s,cv3s,"NonSquare 8");

#if (XTEST & 8)
    tmv::Matrix<T> a6x(4,0,1);
    tmv::Matrix<T> a7x(0,4,1);
    tmv::Matrix<CT> ca6x = a6x;
    tmv::Matrix<CT> ca7x = a7x;

    tmv::MatrixView<T> a6 = a6x.view();
    tmv::MatrixView<CT> ca6 = ca6x.view();
    tmv::MatrixView<T> a7 = a7x.view();
    tmv::MatrixView<CT> ca7 = ca7x.view();

    tmv::Vector<T> v4 = a6.row(2);
    tmv::VectorView<T> v4v = v4.view();
    tmv::Vector<T> v45(0);
    tmv::VectorView<T> v4s = v45.subVector(0,0,5);

    tmv::Vector<CT> cv4 = ca6.row(2);
    tmv::VectorView<CT> cv4v = cv4.view();
    tmv::Vector<CT> cv45(0);
    tmv::VectorView<CT> cv4s = cv45.subVector(0,0,5);

    TestMatrixArith7(a6,ca6,v1v,cv1v,v4v,cv4v,"Degenerate 1");
    TestMatrixArith7(a6,ca6,v1s,cv1s,v4v,cv4v,"Degenerate 2");
    TestMatrixArith7(a6,ca6,v1v,cv1v,v4s,cv4s,"Degenerate 3");
    TestMatrixArith7(a6,ca6,v1s,cv1s,v4s,cv4s,"Degenerate 4");
    TestMatrixArith7(a7,ca7,v4v,cv4v,v2v,cv2v,"Degenerate 5");
    TestMatrixArith7(a7,ca7,v4s,cv4s,v2v,cv2v,"Degenerate 6");
    TestMatrixArith7(a7,ca7,v4v,cv4v,v2s,cv2s,"Degenerate 7");
    TestMatrixArith7(a7,ca7,v4s,cv4s,v2s,cv2s,"Degenerate 8");
#endif
}

#ifdef TEST_DOUBLE
template void TestMatrixArith_7<double>();
#endif
#ifdef TEST_FLOAT
template void TestMatrixArith_7<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestMatrixArith_7<long double>();
#endif
#ifdef TEST_INT
template void TestMatrixArith_7<int>();
#endif
