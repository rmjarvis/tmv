#include "TMV_Test.h"
#include "TMV_Tri.h"
#include "TMV.h"

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanAddEq(
    const tmv::UpperTriMatrixView<T1,I1>& a,
    const tmv::UpperTriMatrixView<T2,I2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanAddEq(
    const tmv::LowerTriMatrixView<T1,I1>& a,
    const tmv::LowerTriMatrixView<T2,I2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanAddEq(
    const tmv::UpperTriMatrixView<T1,I1>& ,
    const tmv::LowerTriMatrixView<T2,I2>& )
{ return false; }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanAddEq(
    const tmv::LowerTriMatrixView<T1,I1>& ,
    const tmv::UpperTriMatrixView<T2,I2>& )
{ return false; }

template <class T1, class T2, tmv::IndexStyle I1> inline bool CanAddEqX(
    const tmv::UpperTriMatrixView<T1,I1>& a, const T2 )
{ return !a.isunit(); }

template <class T1, class T2, tmv::IndexStyle I1> inline bool CanAddEqX(
    const tmv::LowerTriMatrixView<T1,I1>& a, const T2 )
{ return !a.isunit(); }

template <class T1, class T2, tmv::IndexStyle I1> inline bool CanMultEqX(
    const tmv::UpperTriMatrixView<T1,I1>& a, const T2 )
{ return !a.isunit(); }

template <class T1, class T2, tmv::IndexStyle I1> inline bool CanMultEqX(
    const tmv::LowerTriMatrixView<T1,I1>& a, const T2 )
{ return !a.isunit(); }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanMultEq(
    const tmv::UpperTriMatrixView<T1,I1>& a,
    const tmv::UpperTriMatrixView<T2,I2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanMultEq(
    const tmv::LowerTriMatrixView<T1,I1>& a,
    const tmv::LowerTriMatrixView<T2,I2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanMultEq(
    const tmv::UpperTriMatrixView<T1,I1>& ,
    const tmv::LowerTriMatrixView<T2,I2>& )
{ return false; }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanMultEq(
    const tmv::LowerTriMatrixView<T1,I1>& ,
    const tmv::UpperTriMatrixView<T2,I2>& )
{ return false; }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanMultEq2(
    const tmv::UpperTriMatrixView<T1,I1>& a,
    const tmv::UpperTriMatrixView<T2,I2>& b)
{ return a.size() == b.size() && !b.isunit(); }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanMultEq2(
    const tmv::LowerTriMatrixView<T1,I1>& a,
    const tmv::LowerTriMatrixView<T2,I2>& b)
{ return a.size() == b.size() && !b.isunit(); }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanMultEq2(
    const tmv::UpperTriMatrixView<T1,I1>& ,
    const tmv::LowerTriMatrixView<T2,I2>& )
{ return false; }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanMultEq2(
    const tmv::LowerTriMatrixView<T1,I1>& ,
    const tmv::UpperTriMatrixView<T2,I2>& )
{ return false; }

template <class T, tmv::IndexStyle I> inline bool CanDoSV(
    const tmv::UpperTriMatrixView<T,I>& )
{ return false; }
template <class T, tmv::IndexStyle I> inline bool CanDoSV(
    const tmv::LowerTriMatrixView<T,I>& )
{ return false; }

#include "TMV_TestMatrixArith.h"

template <class T> void TestTriMatrixArith_A()
{
  const int N = 10;

  tmv::Matrix<T> a1(N,N);
  for(int i=0;i<N;i++) for(int j=0;j<N;j++) a1(i,j) = 12+3*i-5*j;

  tmv::UpperTriMatrix<T,tmv::NonUnitDiag> u1(a1);
  tmv::UpperTriMatrix<T,tmv::UnitDiag> u2(a1);
  tmv::LowerTriMatrix<T,tmv::NonUnitDiag> l1(a1);
  tmv::LowerTriMatrix<T,tmv::UnitDiag> l2(a1);

  tmv::Matrix<std::complex<T> > c1 = a1 * std::complex<T>(1,2);
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> cu1(c1);
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> cu2(c1);
  tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag> cl1(c1);
  tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cl2(c1);

  TestMatrixArith<T,tmv::LowerTriMatrix<T>,tmv::LowerTriMatrix<std::complex<T> > >(
	l1.View(),u1.View(),cl1.View(),cu1.View(),"Tri/Tri");
#ifdef XTEST
  TestMatrixArith<T,tmv::LowerTriMatrix<T>,tmv::LowerTriMatrix<std::complex<T> > >(
	l1.View(),u2.View(),cl1.View(),cu2.View(),"Tri/Tri");
  TestMatrixArith<T,tmv::LowerTriMatrix<T,tmv::UnitDiag>,tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> >(
	l2.View(),u1.View(),cl2.View(),cu1.View(),"Tri/Tri");
  TestMatrixArith<T,tmv::LowerTriMatrix<T,tmv::UnitDiag>,tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> >(
	l2.View(),u2.View(),cl2.View(),cu2.View(),"Tri/Tri");
  if (doallarith) {
    TestMatrixArith<T,tmv::LowerTriMatrix<T>,tmv::LowerTriMatrix<std::complex<T> > >(
	l1.View(),l2.View(),cl1.View(),cl2.View(),"Tri/Tri");
    TestMatrixArith<T,tmv::LowerTriMatrix<T,tmv::UnitDiag>,tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> >(
	l2.View(),l1.View(),cl2.View(),cl1.View(),"Tri/Tri");
    TestMatrixArith<T,tmv::UpperTriMatrix<T>,tmv::UpperTriMatrix<std::complex<T> > >(
	u1.View(),u2.View(),cu1.View(),cu2.View(),"Tri/Tri");
    TestMatrixArith<T,tmv::UpperTriMatrix<T,tmv::UnitDiag>,tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> >(
	u2.View(),u1.View(),cu2.View(),cu1.View(),"Tri/Tri");
    TestMatrixArith<T,tmv::UpperTriMatrix<T>,tmv::UpperTriMatrix<std::complex<T> > >(
	u1.View(),l1.View(),cu1.View(),cl1.View(),"Tri/Tri");
    TestMatrixArith<T,tmv::UpperTriMatrix<T>,tmv::UpperTriMatrix<std::complex<T> > >(
	u1.View(),l2.View(),cu1.View(),cl2.View(),"Tri/Tri");
    TestMatrixArith<T,tmv::UpperTriMatrix<T,tmv::UnitDiag>,tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> >(
	u2.View(),l1.View(),cu2.View(),cl1.View(),"Tri/Tri");
    TestMatrixArith<T,tmv::UpperTriMatrix<T,tmv::UnitDiag>,tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> >(
	u2.View(),l2.View(),cu2.View(),cl2.View(),"Tri/Tri");
  }
#endif

  std::cout<<"TriMatrix<"<<tmv::Type(T())<<"> (Tri/Tri) Arithmetic passed all tests\n";
}

template void TestTriMatrixArith_A<double>();
#ifndef NOFLOAT
template void TestTriMatrixArith_A<float>();
#endif
#ifdef LONGDOUBLE
template void TestTriMatrixArith_A<long double>();
#endif

