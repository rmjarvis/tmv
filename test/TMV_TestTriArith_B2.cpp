#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Test1.h"

#define NOSV
#include "TMV_TestMatrixArith.h"

template <class T> void TestTriMatrixArith_B2()
{
#ifdef XTEST
    const int N = 10;

    tmv::Matrix<T> a1(N,N);
    for(int i=0;i<N;i++) for(int j=0;j<N;j++) a1(i,j) = T(12+3*i-5*j);

    tmv::LowerTriMatrix<T,tmv::NonUnitDiag,tmv::RowMajor> l1(a1);

    tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(1,2);
    tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::RowMajor> cl1(ca1);

    tmv::LowerTriMatrixView<T> l1v = l1.view();
    tmv::MatrixView<T> a1v = a1.view();
    tmv::LowerTriMatrixView<std::complex<T> > cl1v = cl1.view();
    tmv::MatrixView<std::complex<T> > ca1v = ca1.view();
    tmv::LowerTriMatrix<T,tmv::NonUnitDiag> l1x = l1v;
    tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag> cl1x = cl1v;

    tmv::UpperTriMatrix<T,tmv::NonUnitDiag,tmv::RowMajor> u1(a1);
    tmv::UpperTriMatrix<T,tmv::UnitDiag,tmv::RowMajor> u2(a1);
    tmv::UpperTriMatrix<T,tmv::NonUnitDiag,tmv::ColMajor> u3(a1);
    tmv::UpperTriMatrix<T,tmv::UnitDiag,tmv::ColMajor> u4(a1);
    tmv::LowerTriMatrix<T,tmv::UnitDiag,tmv::RowMajor> l2(a1);
    tmv::LowerTriMatrix<T,tmv::NonUnitDiag,tmv::ColMajor> l3(a1);
    tmv::LowerTriMatrix<T,tmv::UnitDiag,tmv::ColMajor> l4(a1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::RowMajor> cu1(ca1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::RowMajor> cu2(ca1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::ColMajor> cu3(ca1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::ColMajor> cu4(ca1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::RowMajor> cl2(ca1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::ColMajor> cl3(ca1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::ColMajor> cl4(ca1);
    tmv::UpperTriMatrixView<T> u1v = u1.view();
    tmv::UpperTriMatrixView<T> u2v = u2.view();
    tmv::UpperTriMatrixView<T> u3v = u3.view();
    tmv::UpperTriMatrixView<T> u4v = u4.view();
    tmv::LowerTriMatrixView<T> l2v = l2.view();
    tmv::LowerTriMatrixView<T> l3v = l3.view();
    tmv::LowerTriMatrixView<T> l4v = l4.view();
    tmv::UpperTriMatrixView<std::complex<T> > cu1v = cu1.view();
    tmv::UpperTriMatrixView<std::complex<T> > cu2v = cu2.view();
    tmv::UpperTriMatrixView<std::complex<T> > cu3v = cu3.view();
    tmv::UpperTriMatrixView<std::complex<T> > cu4v = cu4.view();
    tmv::LowerTriMatrixView<std::complex<T> > cl2v = cl2.view();
    tmv::LowerTriMatrixView<std::complex<T> > cl3v = cl3.view();
    tmv::LowerTriMatrixView<std::complex<T> > cl4v = cl4.view();
    tmv::UpperTriMatrix<T,tmv::NonUnitDiag> u1x = u1v;
    tmv::UpperTriMatrix<T,tmv::UnitDiag> u2x = u2v;
    tmv::LowerTriMatrix<T,tmv::UnitDiag> l2x = l2v;
    tmv::Matrix<T> a1x = a1v;
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> cu1x = cu1v;
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> cu2x = cu2v;
    tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cl2x = cl2v;
    tmv::Matrix<std::complex<T> > ca1x = ca1v;

    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,u1v,cu1v,"Square/Tri");
    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,u2v,cu2v,"Square/Tri");
    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,u3v,cu3v,"Square/Tri");
    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,u4v,cu4v,"Square/Tri");
    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,l1v,cl1v,"Square/Tri");
    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,l2v,cl2v,"Square/Tri");
    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,l3v,cl3v,"Square/Tri");
    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,l4v,cl4v,"Square/Tri");

    tmv::Matrix<T> a2(15,10);
    for(int i=0;i<15;++i) for(int j=0;j<10;++j) a2(i,j) = T(1-3*i+2*j);
    tmv::Matrix<std::complex<T> > ca2 = a2*std::complex<T>(1,2);
    tmv::MatrixView<T> a2v = a2.view();
    tmv::MatrixView<std::complex<T> > ca2v = ca2.view();
    tmv::Matrix<T> a2x = a2v;
    tmv::Matrix<std::complex<T> > ca2x = ca2v;

    TestMatrixArith456<T>(a2x,ca2x,a2v,ca2v,u1v,cu1v,"NonSquare/Tri");
    TestMatrixArith456<T>(a2x,ca2x,a2v,ca2v,u2v,cu2v,"NonSquare/Tri");
    TestMatrixArith456<T>(a2x,ca2x,a2v,ca2v,u3v,cu3v,"NonSquare/Tri");
    TestMatrixArith456<T>(a2x,ca2x,a2v,ca2v,u4v,cu4v,"NonSquare/Tri");
    TestMatrixArith456<T>(a2x,ca2x,a2v,ca2v,l1v,cl1v,"NonSquare/Tri");
    TestMatrixArith456<T>(a2x,ca2x,a2v,ca2v,l2v,cl2v,"NonSquare/Tri");
    TestMatrixArith456<T>(a2x,ca2x,a2v,ca2v,l3v,cl3v,"NonSquare/Tri");
    TestMatrixArith456<T>(a2x,ca2x,a2v,ca2v,l4v,cl4v,"NonSquare/Tri");

    tmv::Matrix<T> a3(10,0,1);
    tmv::Matrix<std::complex<T> > ca3 = a3;
    tmv::MatrixView<T> a3v = a3.view();
    tmv::MatrixView<std::complex<T> > ca3v = ca3.view();
    tmv::Matrix<T> a3x = a3v;
    tmv::Matrix<std::complex<T> > ca3x = ca3v;

    TestMatrixArith456<T>(a3x,ca3x,a3v,ca3v,u1v,cu1v,"Degenerate/Tri");
    TestMatrixArith456<T>(a3x,ca3x,a3v,ca3v,u2v,cu2v,"Degenerate/Tri");
    TestMatrixArith456<T>(a3x,ca3x,a3v,ca3v,u3v,cu3v,"Degenerate/Tri");
    TestMatrixArith456<T>(a3x,ca3x,a3v,ca3v,u4v,cu4v,"Degenerate/Tri");
    TestMatrixArith456<T>(a3x,ca3x,a3v,ca3v,l1v,cl1v,"Degenerate/Tri");
    TestMatrixArith456<T>(a3x,ca3x,a3v,ca3v,l2v,cl2v,"Degenerate/Tri");
    TestMatrixArith456<T>(a3x,ca3x,a3v,ca3v,l3v,cl3v,"Degenerate/Tri");
    TestMatrixArith456<T>(a3x,ca3x,a3v,ca3v,l4v,cl4v,"Degenerate/Tri");
#endif
}

#ifdef INST_DOUBLE
template void TestTriMatrixArith_B2<double>();
#endif
#ifdef INST_FLOAT
template void TestTriMatrixArith_B2<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestTriMatrixArith_B2<long double>();
#endif
#ifdef INST_INT
template void TestTriMatrixArith_B2<int>();
#endif

