#ifdef NDEBUG
#undef NDEBUG
#endif
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV_Tri.h"
#include "TMV.h"

template <class T1, class T2> inline bool CanAddEq(
    const tmv::UpperTriMatrixView<T1>& a,
    const tmv::UpperTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanAddEq(
    const tmv::LowerTriMatrixView<T1>& a,
    const tmv::LowerTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> inline bool CanAddEqX(
    const tmv::UpperTriMatrixView<T1>& a, const T2 )
{ return !a.isunit(); }

template <class T1, class T2> inline bool CanAddEqX(
    const tmv::LowerTriMatrixView<T1>& a, const T2 )
{ return !a.isunit(); }

template <class T1, class T2> inline bool CanMultEqX(
    const tmv::UpperTriMatrixView<T1>& a, const T2 )
{ return !a.isunit(); }

template <class T1, class T2> inline bool CanMultEqX(
    const tmv::LowerTriMatrixView<T1>& a, const T2 )
{ return !a.isunit(); }

template <class T1, class T2, class T3> inline bool CanMult(
    const tmv::UpperTriMatrixView<T1>& a,
    const tmv::UpperTriMatrixView<T2>& b,
    const tmv::UpperTriMatrixView<T3>& c)
{
  return a.size() == b.size() && a.size() == b.size() && 
    ((a.isunit() && b.isunit() && c.isunit()) || (!c.isunit()));
}

template <class T1, class T2, class T3> inline bool CanMultXM(
    const tmv::UpperTriMatrixView<T1>& a,
    const tmv::UpperTriMatrixView<T2>& b,
    const tmv::UpperTriMatrixView<T3>& c)
{ return a.size() == b.size() && a.size() == b.size() && !c.isunit(); }

template <class T1, class T2, class T3> inline bool CanMult(
    const tmv::LowerTriMatrixView<T1>& a,
    const tmv::LowerTriMatrixView<T2>& b,
    const tmv::LowerTriMatrixView<T3>& c)
{
  return a.size() == b.size() && a.size() == b.size() && 
    ((a.isunit() && b.isunit() && c.isunit()) || (!c.isunit()));
}

template <class T1, class T2, class T3> inline bool CanMultXM(
    const tmv::LowerTriMatrixView<T1>& a,
    const tmv::LowerTriMatrixView<T2>& b,
    const tmv::LowerTriMatrixView<T3>& c)
{ return a.size() == b.size() && a.size() == b.size() && !c.isunit(); }

#define NOSV
#include "TMV_TestMatrixArith.h"

template <class T> void TestTriMatrixArith_A1()
{
  const int N = 10;

  tmv::Matrix<T> a1(N,N);
  for(int i=0;i<N;i++) for(int j=0;j<N;j++) a1(i,j) = T(12+3*i-5*j);

  tmv::UpperTriMatrix<T,tmv::NonUnitDiag,tmv::RowMajor> u1(a1);
  tmv::LowerTriMatrix<T,tmv::NonUnitDiag,tmv::RowMajor> l1(a1);

  tmv::Matrix<std::complex<T> > c1 = a1 * std::complex<T>(1,2);
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::RowMajor> cu1(c1);
  tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::RowMajor> cl1(c1);

  tmv::LowerTriMatrixView<T> l1v = l1.View();
  tmv::LowerTriMatrixView<std::complex<T> > cl1v = cl1.View();
  tmv::LowerTriMatrix<T,tmv::NonUnitDiag> l1x = l1v;
  tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag> cl1x = cl1v;

  TestMatrixArith123<T>(l1x,cl1x,l1v,cl1v,"Tri L");

#ifdef XTEST
  tmv::UpperTriMatrixView<T> u1v = u1.View();
  tmv::UpperTriMatrixView<std::complex<T> > cu1v = cu1.View();
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

  TestMatrixArith123<T>(l2x,cl2x,l2v,cl2v,"Tri L");
  TestMatrixArith123<T>(u1x,cu1x,u1v,cu1v,"Tri U");
  TestMatrixArith123<T>(u2x,cu2x,u2v,cu2v,"Tri U");
  TestMatrixArith123<T>(l1x,cl1x,l3v,cl3v,"Tri L");
  TestMatrixArith123<T>(l2x,cl2x,l4v,cl4v,"Tri L");
  TestMatrixArith123<T>(u1x,cu1x,u3v,cu3v,"Tri U");
  TestMatrixArith123<T>(u2x,cu2x,u4v,cu4v,"Tri U");
  TestMatrixArith45<T>(l1x,cl1x,l1v,cl1v,l2v,cl2v,"Tri L/L");
  TestMatrixArith45<T>(l1x,cl1x,l1v,cl1v,l3v,cl3v,"Tri L/L");
  TestMatrixArith45<T>(l1x,cl1x,l1v,cl1v,l4v,cl4v,"Tri L/L");
  TestMatrixArith45<T>(l2x,cl2x,l2v,cl2v,l1v,cl1v,"Tri L/L");
  TestMatrixArith45<T>(l2x,cl2x,l2v,cl2v,l3v,cl3v,"Tri L/L");
  TestMatrixArith45<T>(l2x,cl2x,l2v,cl2v,l4v,cl4v,"Tri L/L");
  TestMatrixArith45<T>(l1x,cl1x,l3v,cl3v,l1v,cl1v,"Tri L/L");
  TestMatrixArith45<T>(l1x,cl1x,l3v,cl3v,l2v,cl2v,"Tri L/L");
  TestMatrixArith45<T>(l1x,cl1x,l3v,cl3v,l4v,cl4v,"Tri L/L");
  TestMatrixArith45<T>(l2x,cl2x,l4v,cl4v,l1v,cl1v,"Tri L/L");
  TestMatrixArith45<T>(l2x,cl2x,l4v,cl4v,l2v,cl2v,"Tri L/L");
  TestMatrixArith45<T>(l2x,cl2x,l4v,cl4v,l3v,cl3v,"Tri L/L");
  TestMatrixArith45<T>(u1x,cu1x,u1v,cu1v,u2v,cu2v,"Tri U/U");
  TestMatrixArith45<T>(u1x,cu1x,u1v,cu1v,u3v,cu3v,"Tri U/U");
  TestMatrixArith45<T>(u1x,cu1x,u1v,cu1v,u4v,cu4v,"Tri U/U");
  TestMatrixArith45<T>(u2x,cu2x,u2v,cu2v,u1v,cu1v,"Tri U/U");
  TestMatrixArith45<T>(u2x,cu2x,u2v,cu2v,u3v,cu3v,"Tri U/U");
  TestMatrixArith45<T>(u2x,cu2x,u2v,cu2v,u4v,cu4v,"Tri U/U");
  TestMatrixArith45<T>(u1x,cu1x,u3v,cu3v,u1v,cu1v,"Tri U/U");
  TestMatrixArith45<T>(u1x,cu1x,u3v,cu3v,u2v,cu2v,"Tri U/U");
  TestMatrixArith45<T>(u1x,cu1x,u3v,cu3v,u4v,cu4v,"Tri U/U");
  TestMatrixArith45<T>(u2x,cu2x,u4v,cu4v,u1v,cu1v,"Tri U/U");
  TestMatrixArith45<T>(u2x,cu2x,u4v,cu4v,u2v,cu2v,"Tri U/U");
  TestMatrixArith45<T>(u2x,cu2x,u4v,cu4v,u3v,cu3v,"Tri U/U");
#endif
}

#ifdef INST_DOUBLE
template void TestTriMatrixArith_A1<double>();
#endif
#ifdef INST_FLOAT
template void TestTriMatrixArith_A1<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestTriMatrixArith_A1<long double>();
#endif
#ifdef INST_INT
template void TestTriMatrixArith_A1<int>();
#endif

