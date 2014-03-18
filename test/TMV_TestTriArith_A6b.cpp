#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"

template <class T>
inline void CopyBackM(
    const tmv::UpperTriMatrix<T>& m0,
    tmv::UpperTriMatrixView<T>& m1)
{
    if (m1.isunit()) m1 = m0.viewAsUnitDiag();
    else m1 = m0;
}

template <class T>
inline void CopyBackM(
    const tmv::LowerTriMatrix<T>& m0,
    tmv::LowerTriMatrixView<T>& m1)
{
    if (m1.isunit()) m1 = m0.viewAsUnitDiag();
    else m1 = m0;
}

template <class T1, class T2, class T3>
static inline bool CanAddElemMultMM(
    const tmv::UpperTriMatrixView<T1>& a,
    const tmv::UpperTriMatrixView<T2>& b,
    const tmv::UpperTriMatrixView<T3>& c)
{
    return (a.size() == b.size()) && (a.size() == c.size()) &&
        !c.isunit(); 
}

template <class T1, class T2, class T3>
static inline bool CanAddElemMultMM(
    const tmv::LowerTriMatrixView<T1>& a,
    const tmv::LowerTriMatrixView<T2>& b,
    const tmv::LowerTriMatrixView<T3>& c)
{
    return (a.size() == b.size()) && (a.size() == c.size()) &&
        !c.isunit(); 
}


#define BASIC_MULTMM_ONLY
#include "TMV_TestMatrixArith.h"

template <class T> 
void TestTriMatrixArith_A6b()
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

    tmv::UpperTriMatrixView<T> u4 = a1x.unitUpperTri();
    tmv::UpperTriMatrixView<CT> cu4 = ca1x.unitUpperTri();
    tmv::UpperTriMatrixView<T> u5 = a2x.unitUpperTri();
    tmv::UpperTriMatrixView<CT> cu5 = ca2x.unitUpperTri();

    TestMatrixArith6(u4,cu4,u4,cu4,u4,cu4,"UpperTri 79");
    TestMatrixArith6(u5,cu5,u5,cu5,u4,cu4,"UpperTri 80");
#if (XTEST & 2)
    TestMatrixArith6(u4,cu4,u5,cu5,u4,cu4,"UpperTri 81");
    TestMatrixArith6(u5,cu5,u4,cu4,u4,cu4,"UpperTri 82");
    TestMatrixArith6(u4,cu4,u5,cu5,u5,cu5,"UpperTri 83");
    TestMatrixArith6(u5,cu5,u4,cu4,u5,cu5,"UpperTri 84");
    TestMatrixArith6(u4,cu4,u4,cu4,u5,cu5,"UpperTri 85");
    TestMatrixArith6(u5,cu5,u5,cu5,u5,cu5,"UpperTri 86");
#endif
#if (XTEST & 1)
    tmv::Matrix<T> a3x(12,16);
    for(int i=0;i<12;++i) for(int j=0;j<16;++j) a3x(i,j) = T(1-2*i+3*j);
    a3x.diag().addToAll(30);
    tmv::Matrix<CT> ca3x = a3x*CT(1,-2);
    ca3x.diag().addToAll(CT(-22,15));

    tmv::UpperTriMatrixView<T> u6 = a3x.subMatrix(0,12,0,16,3,4).unitUpperTri();
    tmv::UpperTriMatrixView<CT> cu6 = ca3x.subMatrix(0,12,0,16,3,4).unitUpperTri();
    TestMatrixArith6(u6,cu6,u4,cu4,u4,cu4,"UpperTri 87");
    TestMatrixArith6(u6,cu6,u5,cu5,u4,cu4,"UpperTri 88");
    TestMatrixArith6(u4,cu4,u6,cu6,u4,cu4,"UpperTri 89");
    TestMatrixArith6(u5,cu5,u6,cu6,u4,cu4,"UpperTri 90");
    TestMatrixArith6(u4,cu4,u4,cu4,u6,cu6,"UpperTri 91");
    TestMatrixArith6(u5,cu5,u5,cu5,u6,cu6,"UpperTri 92");
    TestMatrixArith6(u4,cu4,u5,cu5,u6,cu6,"UpperTri 93");
    TestMatrixArith6(u5,cu5,u4,cu4,u6,cu6,"UpperTri 94");
#endif

    tmv::LowerTriMatrixView<T> l4 = a1x.unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl4 = ca1x.unitLowerTri();
    tmv::LowerTriMatrixView<T> l5 = a2x.unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl5 = ca2x.unitLowerTri();

    TestMatrixArith6(l4,cl4,l4,cl4,l4,cl4,"LowerTri 79");
    TestMatrixArith6(l5,cl5,l5,cl5,l4,cl4,"LowerTri 80");
#if (XTEST & 2)
    TestMatrixArith6(l4,cl4,l5,cl5,l4,cl4,"LowerTri 81");
    TestMatrixArith6(l5,cl5,l4,cl4,l4,cl4,"LowerTri 82");
    TestMatrixArith6(l4,cl4,l5,cl5,l5,cl5,"LowerTri 83");
    TestMatrixArith6(l5,cl5,l4,cl4,l5,cl5,"LowerTri 84");
    TestMatrixArith6(l4,cl4,l4,cl4,l5,cl5,"LowerTri 85");
    TestMatrixArith6(l5,cl5,l5,cl5,l5,cl5,"LowerTri 86");
#endif
#if (XTEST & 1)
    tmv::LowerTriMatrixView<T> l6 = a3x.subMatrix(0,12,0,16,3,4).unitLowerTri();
    tmv::LowerTriMatrixView<CT> cl6 = ca3x.subMatrix(0,12,0,16,3,4).unitLowerTri();
    TestMatrixArith6(l6,cl6,l4,cl4,l4,cl4,"LowerTri 87");
    TestMatrixArith6(l6,cl6,l5,cl5,l4,cl4,"LowerTri 88");
    TestMatrixArith6(l4,cl4,l6,cl6,l4,cl4,"LowerTri 89");
    TestMatrixArith6(l5,cl5,l6,cl6,l4,cl4,"LowerTri 90");
    TestMatrixArith6(l4,cl4,l4,cl4,l6,cl6,"LowerTri 91");
    TestMatrixArith6(l5,cl5,l5,cl5,l6,cl6,"LowerTri 92");
    TestMatrixArith6(l4,cl4,l5,cl5,l6,cl6,"LowerTri 93");
    TestMatrixArith6(l5,cl5,l4,cl4,l6,cl6,"LowerTri 94");
#endif
}

#ifdef TEST_DOUBLE
template void TestTriMatrixArith_A6b<double>();
#endif
#ifdef TEST_FLOAT
template void TestTriMatrixArith_A6b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestTriMatrixArith_A6b<long double>();
#endif
#ifdef TEST_INT
template void TestTriMatrixArith_A6b<int>();
#endif
