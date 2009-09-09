// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

#define NOADDEQ
#define NOMULTEQ
#define NOSV
#define NOTRANSMM
#include "TMV_TestMatrixArith.h"

template <class T> void TestTriMatrixArith_B1()
{
  const int N = 10;

  tmv::Matrix<T> a1(N,N);
  for(int i=0;i<N;i++) for(int j=0;j<N;j++) a1(i,j) = T(12+3*i-5*j);

  tmv::LowerTriMatrix<T,tmv::NonUnitDiag,tmv::RowMajor> l1(a1);

  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(1,2);
  tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::RowMajor> cl1(ca1);

  tmv::LowerTriMatrixView<T> l1v = l1.View();
  tmv::MatrixView<T> a1v = a1.View();
  tmv::LowerTriMatrixView<std::complex<T> > cl1v = cl1.View();
  tmv::MatrixView<std::complex<T> > ca1v = ca1.View();
  tmv::LowerTriMatrix<T,tmv::NonUnitDiag> l1x = l1v;
  tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag> cl1x = cl1v;

  TestMatrixArith456<T>(l1x,cl1x,l1v,cl1v,a1v,ca1v,"Tri/Square L");

#ifdef XTEST
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
  tmv::UpperTriMatrixView<T> u1v = u1.View();
  tmv::UpperTriMatrixView<T> u2v = u2.View();
  tmv::UpperTriMatrixView<T> u3v = u3.View();
  tmv::UpperTriMatrixView<T> u4v = u4.View();
  tmv::LowerTriMatrixView<T> l2v = l2.View();
  tmv::LowerTriMatrixView<T> l3v = l3.View();
  tmv::LowerTriMatrixView<T> l4v = l4.View();
  tmv::UpperTriMatrixView<std::complex<T> > cu1v = cu1.View();
  tmv::UpperTriMatrixView<std::complex<T> > cu2v = cu2.View();
  tmv::UpperTriMatrixView<std::complex<T> > cu3v = cu3.View();
  tmv::UpperTriMatrixView<std::complex<T> > cu4v = cu4.View();
  tmv::LowerTriMatrixView<std::complex<T> > cl2v = cl2.View();
  tmv::LowerTriMatrixView<std::complex<T> > cl3v = cl3.View();
  tmv::LowerTriMatrixView<std::complex<T> > cl4v = cl4.View();
  tmv::UpperTriMatrix<T,tmv::NonUnitDiag> u1x = u1v;
  tmv::UpperTriMatrix<T,tmv::UnitDiag> u2x = u2v;
  tmv::LowerTriMatrix<T,tmv::UnitDiag> l2x = l2v;
  tmv::Matrix<T> a1x = a1v;
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> cu1x = cu1v;
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> cu2x = cu2v;
  tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cl2x = cl2v;
  tmv::Matrix<std::complex<T> > ca1x = ca1v;

  TestMatrixArith456<T>(l2x,cl2x,l2v,cl2v,a1v,ca1v,"Tri/Square L");
  TestMatrixArith456<T>(l1x,cl1x,l3v,cl3v,a1v,ca1v,"Tri/Square L");
  TestMatrixArith456<T>(l2x,cl2x,l4v,cl4v,a1v,ca1v,"Tri/Square L");
  TestMatrixArith456<T>(u1x,cu1x,u1v,cu1v,a1v,ca1v,"Tri/Square U");
  TestMatrixArith456<T>(u2x,cu2x,u2v,cu2v,a1v,ca1v,"Tri/Square U");
  TestMatrixArith456<T>(u1x,cu1x,u3v,cu3v,a1v,ca1v,"Tri/Square U");
  TestMatrixArith456<T>(u2x,cu2x,u4v,cu4v,a1v,ca1v,"Tri/Square U");

  tmv::Matrix<T> a2(15,10);
  for(int i=0;i<15;++i) for(int j=0;j<10;++j) a2(i,j) = T(1-3*i+2*j);
  tmv::Matrix<std::complex<T> > ca2 = a2*std::complex<T>(1,2);
  tmv::MatrixView<T> a2v = a2.View();
  tmv::MatrixView<std::complex<T> > ca2v = ca2.View();
  tmv::Matrix<T> a2x = a2v;
  tmv::Matrix<std::complex<T> > ca2x = ca2v;

  TestMatrixArith456<T>(l1x,cl1x,l1v,cl1v,a2v,ca2v,"Tri/NonSquare");
  TestMatrixArith456<T>(l2x,cl2x,l2v,cl2v,a2v,ca2v,"Tri/NonSquare");
  TestMatrixArith456<T>(l1x,cl1x,l3v,cl3v,a2v,ca2v,"Tri/NonSquare");
  TestMatrixArith456<T>(l2x,cl2x,l4v,cl4v,a2v,ca2v,"Tri/NonSquare");
  TestMatrixArith456<T>(u1x,cu1x,u1v,cu1v,a2v,ca2v,"Tri/NonSquare");
  TestMatrixArith456<T>(u2x,cu2x,u2v,cu2v,a2v,ca2v,"Tri/NonSquare");
  TestMatrixArith456<T>(u1x,cu1x,u3v,cu3v,a2v,ca2v,"Tri/NonSquare");
  TestMatrixArith456<T>(u2x,cu2x,u4v,cu4v,a2v,ca2v,"Tri/NonSquare");

  tmv::Matrix<T> a3(10,0,1);
  tmv::Matrix<std::complex<T> > ca3 = a3;
  tmv::MatrixView<T> a3v = a3.View();
  tmv::MatrixView<std::complex<T> > ca3v = ca3.View();
  tmv::Matrix<T> a3x = a3v;
  tmv::Matrix<std::complex<T> > ca3x = ca3v;

  TestMatrixArith456<T>(l1x,cl1x,l1v,cl1v,a3v,ca3v,"Tri/Degenerate");
  TestMatrixArith456<T>(l2x,cl2x,l2v,cl2v,a3v,ca3v,"Tri/Degenerate");
  TestMatrixArith456<T>(l1x,cl1x,l3v,cl3v,a3v,ca3v,"Tri/Degenerate");
  TestMatrixArith456<T>(l2x,cl2x,l4v,cl4v,a3v,ca3v,"Tri/Degenerate");
  TestMatrixArith456<T>(u1x,cu1x,u1v,cu1v,a3v,ca3v,"Tri/Degenerate");
  TestMatrixArith456<T>(u2x,cu2x,u2v,cu2v,a3v,ca3v,"Tri/Degenerate");
  TestMatrixArith456<T>(u1x,cu1x,u3v,cu3v,a3v,ca3v,"Tri/Degenerate");
  TestMatrixArith456<T>(u2x,cu2x,u4v,cu4v,a3v,ca3v,"Tri/Degenerate");
#endif
}

#ifdef TEST_DOUBLE
template void TestTriMatrixArith_B1<double>();
#endif
#ifdef TEST_FLOAT
template void TestTriMatrixArith_B1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestTriMatrixArith_B1<long double>();
#endif
#ifdef TEST_INT
template void TestTriMatrixArith_B1<int>();
#endif

