#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

#define NOADDEQ
#define NOMULTEQ
#define NOSV
#define NOTRANSMM
#include "TMV_TestMatrixArith.h"

template <class T> 
void TestTriMatrixArith_A2()
{
    const int N = 10;

    tmv::Matrix<T> a1(N,N);
    for(int i=0;i<N;i++) for(int j=0;j<N;j++) a1(i,j) = T(13+3*i-5*j);

    tmv::UpperTriMatrix<T,tmv::NonUnitDiag,tmv::RowMajor> u1(a1);
    tmv::LowerTriMatrix<T,tmv::NonUnitDiag,tmv::RowMajor> l1(a1);

    tmv::Matrix<std::complex<T> > c1 = a1 * std::complex<T>(1,2);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::RowMajor> cu1(c1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::RowMajor> cl1(c1);

    tmv::UpperTriMatrixView<T> u1v = u1.View();
    tmv::LowerTriMatrixView<T> l1v = l1.View();
    tmv::UpperTriMatrixView<std::complex<T> > cu1v = cu1.View();
    tmv::LowerTriMatrixView<std::complex<T> > cl1v = cl1.View();
    tmv::LowerTriMatrix<T,tmv::NonUnitDiag> l1x = l1v;
    tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag> cl1x = cl1v;

    TestMatrixArith456<T>(l1x,cl1x,l1v,cl1v,u1v,cu1v,"Tri L/U");

#ifdef XTEST
    tmv::UpperTriMatrix<T,tmv::UnitDiag,tmv::RowMajor> u2(a1);
    tmv::UpperTriMatrix<T,tmv::NonUnitDiag,tmv::ColMajor> u3(a1);
    tmv::UpperTriMatrix<T,tmv::UnitDiag,tmv::ColMajor> u4(a1);
    tmv::LowerTriMatrix<T,tmv::UnitDiag,tmv::RowMajor> l2(a1);
    tmv::LowerTriMatrix<T,tmv::NonUnitDiag,tmv::ColMajor> l3(a1);
    tmv::LowerTriMatrix<T,tmv::UnitDiag,tmv::ColMajor> l4(a1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::RowMajor> cu2(c1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::ColMajor> cu3(c1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::ColMajor> cu4(c1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::RowMajor> cl2(c1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::ColMajor> cl3(c1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::ColMajor> cl4(c1);
    tmv::UpperTriMatrixView<T> u2v = u2.View();
    tmv::UpperTriMatrixView<T> u3v = u3.View();
    tmv::UpperTriMatrixView<T> u4v = u4.View();
    tmv::LowerTriMatrixView<T> l2v = l2.View();
    tmv::LowerTriMatrixView<T> l3v = l3.View();
    tmv::LowerTriMatrixView<T> l4v = l4.View();
    tmv::UpperTriMatrixView<std::complex<T> > cu2v = cu2.View();
    tmv::UpperTriMatrixView<std::complex<T> > cu3v = cu3.View();
    tmv::UpperTriMatrixView<std::complex<T> > cu4v = cu4.View();
    tmv::LowerTriMatrixView<std::complex<T> > cl2v = cl2.View();
    tmv::LowerTriMatrixView<std::complex<T> > cl3v = cl3.View();
    tmv::LowerTriMatrixView<std::complex<T> > cl4v = cl4.View();
    tmv::UpperTriMatrix<T,tmv::NonUnitDiag> u1x = u1v;
    tmv::UpperTriMatrix<T,tmv::UnitDiag> u2x = u2v;
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> cu1x = cu1v;
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> cu2x = cu2v;
    tmv::LowerTriMatrix<T,tmv::UnitDiag> l2x = l2v;
    tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cl2x = cl2v;

    TestMatrixArith456<T>(l1x,cl1x,l1v,cl1v,u2v,cu2v,"Tri L/U");
    TestMatrixArith456<T>(l2x,cl2x,l2v,cl2v,u1v,cu1v,"Tri L/U");
    TestMatrixArith456<T>(l2x,cl2x,l2v,cl2v,u2v,cu2v,"Tri L/U");
    TestMatrixArith456<T>(l1x,cl1x,l1v,cl1v,u3v,cu3v,"Tri L/U");
    TestMatrixArith456<T>(l1x,cl1x,l1v,cl1v,u4v,cu4v,"Tri L/U");
    TestMatrixArith456<T>(l2x,cl2x,l2v,cl2v,u3v,cu3v,"Tri L/U");
    TestMatrixArith456<T>(l2x,cl2x,l2v,cl2v,u4v,cu4v,"Tri L/U");
    TestMatrixArith456<T>(l1x,cl1x,l3v,cl3v,u1v,cu1v,"Tri L/U");
    TestMatrixArith456<T>(l1x,cl1x,l3v,cl3v,u2v,cu2v,"Tri L/U");
    TestMatrixArith456<T>(l1x,cl1x,l3v,cl3v,u3v,cu3v,"Tri L/U");
    TestMatrixArith456<T>(l1x,cl1x,l3v,cl3v,u4v,cu4v,"Tri L/U");
    TestMatrixArith456<T>(l2x,cl2x,l4v,cl4v,u1v,cu1v,"Tri L/U");
    TestMatrixArith456<T>(l2x,cl2x,l4v,cl4v,u2v,cu2v,"Tri L/U");
    TestMatrixArith456<T>(l2x,cl2x,l4v,cl4v,u3v,cu3v,"Tri L/U");
    TestMatrixArith456<T>(l2x,cl2x,l4v,cl4v,u4v,cu4v,"Tri L/U");
    TestMatrixArith456<T>(u1x,cu1x,u1v,cu1v,l1v,cl1v,"Tri U/L");
    TestMatrixArith456<T>(u1x,cu1x,u1v,cu1v,l2v,cl2v,"Tri U/L");
    TestMatrixArith456<T>(u1x,cu1x,u1v,cu1v,l3v,cl3v,"Tri U/L");
    TestMatrixArith456<T>(u1x,cu1x,u1v,cu1v,l4v,cl4v,"Tri U/L");
    TestMatrixArith456<T>(u2x,cu2x,u2v,cu2v,l1v,cl1v,"Tri U/L");
    TestMatrixArith456<T>(u2x,cu2x,u2v,cu2v,l2v,cl2v,"Tri U/L");
    TestMatrixArith456<T>(u2x,cu2x,u2v,cu2v,l3v,cl3v,"Tri U/L");
    TestMatrixArith456<T>(u2x,cu2x,u2v,cu2v,l4v,cl4v,"Tri U/L");
    TestMatrixArith456<T>(u1x,cu1x,u3v,cu3v,l1v,cl1v,"Tri U/L");
    TestMatrixArith456<T>(u1x,cu1x,u3v,cu3v,l2v,cl2v,"Tri U/L");
    TestMatrixArith456<T>(u1x,cu1x,u3v,cu3v,l3v,cl3v,"Tri U/L");
    TestMatrixArith456<T>(u1x,cu1x,u3v,cu3v,l4v,cl4v,"Tri U/L");
    TestMatrixArith456<T>(u2x,cu2x,u4v,cu4v,l1v,cl1v,"Tri U/L");
    TestMatrixArith456<T>(u2x,cu2x,u4v,cu4v,l2v,cl2v,"Tri U/L");
    TestMatrixArith456<T>(u2x,cu2x,u4v,cu4v,l3v,cl3v,"Tri U/L");
    TestMatrixArith456<T>(u2x,cu2x,u4v,cu4v,l4v,cl4v,"Tri U/L");
#endif
}

#ifdef INST_DOUBLE
template void TestTriMatrixArith_A2<double>();
#endif
#ifdef INST_FLOAT
template void TestTriMatrixArith_A2<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestTriMatrixArith_A2<long double>();
#endif
#ifdef INST_INT
template void TestTriMatrixArith_A2<int>();
#endif

