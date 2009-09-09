// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

#define NOADDEQ
#define NOMULTEQ
#define NOTRANSMM
#include "TMV_TestMatrixArith.h"

template <class T> void TestTriMatrixArith_C2()
{
#ifdef XTEST
  const int N = 10;

  tmv::Matrix<T> a1(N,N);
  for(int i=0;i<N;i++) for(int j=0;j<N;j++) a1(i,j) = T(12+3*i-5*j);
  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(1,2);

  tmv::DiagMatrix<T> d1(a1.diag());
  tmv::DiagMatrix<std::complex<T> > cd1(ca1.diag());

  tmv::LowerTriMatrix<T,tmv::NonUnitDiag,tmv::RowMajor> l1(a1);
  tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::RowMajor> cl1(ca1);
  tmv::LowerTriMatrixView<T> l1v = l1.View();
  tmv::LowerTriMatrixView<std::complex<T> > cl1v = cl1.View();
  tmv::LowerTriMatrix<T,tmv::NonUnitDiag> l1x = l1v;
  tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag> cl1x = cl1v;
  tmv::DiagMatrixView<T> d1v = d1.View();
  tmv::DiagMatrixView<std::complex<T> > cd1v = cd1.View();

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
  tmv::DiagMatrix<T> d1x = d1v;
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> cu1x = cu1v;
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> cu2x = cu2v;
  tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cl2x = cl2v;
  tmv::DiagMatrix<std::complex<T> > cd1x = cd1v;

  TestMatrixArith456<T>(u1x,cu1x,d1v,cd1v,u1v,cu1v,"Diag/Tri U");
  TestMatrixArith456<T>(u1x,cu1x,d1v,cd1v,u2v,cu2v,"Diag/Tri U");
  TestMatrixArith456<T>(u1x,cu1x,d1v,cd1v,u3v,cu3v,"Diag/Tri U");
  TestMatrixArith456<T>(u1x,cu1x,d1v,cd1v,u4v,cu4v,"Diag/Tri U");
  TestMatrixArith456<T>(l1x,cl1x,d1v,cd1v,l1v,cl1v,"Diag/Tri L");
  TestMatrixArith456<T>(l1x,cl1x,d1v,cd1v,l2v,cl2v,"Diag/Tri L");
  TestMatrixArith456<T>(l1x,cl1x,d1v,cd1v,l3v,cl3v,"Diag/Tri L");
  TestMatrixArith456<T>(l1x,cl1x,d1v,cd1v,l4v,cl4v,"Diag/Tri L");
#endif
}

#ifdef INST_DOUBLE
template void TestTriMatrixArith_C2<double>();
#endif
#ifdef INST_FLOAT
template void TestTriMatrixArith_C2<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestTriMatrixArith_C2<long double>();
#endif
#ifdef INST_INT
template void TestTriMatrixArith_C2<int>();
#endif

