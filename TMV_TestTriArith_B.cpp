#include "TMV_Test.h"
#include "TMV_Tri.h"
#include "TMV.h"

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanAddEq(
    const tmv::UpperTriMatrixView<T1,I1>& , const tmv::MatrixView<T2,I2>& )
{ return false; }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanAddEq(
    const tmv::LowerTriMatrixView<T1,I1>& , const tmv::MatrixView<T2,I2>& )
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
    const tmv::UpperTriMatrixView<T1,I1>& , const tmv::MatrixView<T2,I2>& )
{ return false; }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanMultEq(
    const tmv::LowerTriMatrixView<T1,I1>& , const tmv::MatrixView<T2,I2>& )
{ return false; }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanMultEq2(
    const tmv::MatrixView<T1,I1>& , const tmv::UpperTriMatrixView<T2,I2>& )
{ return false; }

template <class T1, class T2, tmv::IndexStyle I1, tmv::IndexStyle I2> 
inline bool CanMultEq2(
    const tmv::MatrixView<T1,I1>& , const tmv::LowerTriMatrixView<T2,I2>& )
{ return false; }

template <class T, tmv::IndexStyle I> inline bool CanDoSV(
    const tmv::UpperTriMatrixView<T,I>& )
{ return false; }
template <class T, tmv::IndexStyle I> inline bool CanDoSV(
    const tmv::LowerTriMatrixView<T,I>& )
{ return false; }

#include "TMV_TestMatrixArith.h"

template <class T> void TestTriMatrixArith_B()
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
      l1.View(),a1.View(),cl1.View(),c1.View(),"Tri/Square");
#ifdef XTEST
  TestMatrixArith<T,tmv::LowerTriMatrix<T,tmv::UnitDiag>,tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> >(
      l2.View(),a1.View(),cl2.View(),c1.View(),"Tri/Square");
  if (doallarith) {
    TestMatrixArith<T,tmv::UpperTriMatrix<T>,tmv::UpperTriMatrix<std::complex<T> > >(
	u1.View(),a1.View(),cu1.View(),c1.View(),"Tri/Square");
    TestMatrixArith<T,tmv::UpperTriMatrix<T,tmv::UnitDiag>,tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> >(
	u2.View(),a1.View(),cu2.View(),c1.View(),"Tri/Square");
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a1.View(),u1.View(),c1.View(),cu1.View(),"Square/Tri");
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a1.View(),u2.View(),c1.View(),cu2.View(),"Square/Tri");
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a1.View(),l1.View(),c1.View(),cl1.View(),"Square/Tri");
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a1.View(),l2.View(),c1.View(),cl2.View(),"Square/Tri");
  }

  tmv::Matrix<T> a2(15,10);
  for(int i=0;i<15;++i) for(int j=0;j<10;++j) a2(i,j) = 1-3*i+2*j;
  tmv::Matrix<std::complex<T> > c2 = a2*std::complex<T>(1,2);

  TestMatrixArith<T,tmv::LowerTriMatrix<T>,tmv::LowerTriMatrix<std::complex<T> > >(
      l1.View(),a2.View(),cl1.View(),c2.View(),"Tri/NonSquare");
  TestMatrixArith<T,tmv::LowerTriMatrix<T,tmv::UnitDiag>,tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> >(
      l2.View(),a2.View(),cl2.View(),c2.View(),"Tri/NonSquare");
  if (doallarith) {
    TestMatrixArith<T,tmv::UpperTriMatrix<T>,tmv::UpperTriMatrix<std::complex<T> > >(
	u1.View(),a2.View(),cu1.View(),c2.View(),"Tri/NonSquare");
    TestMatrixArith<T,tmv::UpperTriMatrix<T,tmv::UnitDiag>,tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> >(
	u2.View(),a2.View(),cu2.View(),c2.View(),"Tri/NonSquare");
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a2.View(),u1.View(),c2.View(),cu1.View(),"NonSquare/Tri");
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a2.View(),u2.View(),c2.View(),cu2.View(),"NonSquare/Tri");
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a2.View(),l1.View(),c2.View(),cl1.View(),"NonSquare/Tri");
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a2.View(),l2.View(),c2.View(),cl2.View(),"NonSquare/Tri");
  }

  tmv::Matrix<T> a3(10,0,1);
  tmv::Matrix<std::complex<T> > c3 = a3;

  TestMatrixArith<T,tmv::LowerTriMatrix<T>,tmv::LowerTriMatrix<std::complex<T> > >(
      l1.View(),a3.View(),cl1.View(),c3.View(),"Tri/Degenerate");
  TestMatrixArith<T,tmv::LowerTriMatrix<T,tmv::UnitDiag>,tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> >(
      l2.View(),a3.View(),cl2.View(),c3.View(),"Tri/Degenerate");
  if (doallarith) {
    TestMatrixArith<T,tmv::UpperTriMatrix<T>,tmv::UpperTriMatrix<std::complex<T> > >(
	u1.View(),a3.View(),cu1.View(),c3.View(),"Tri/Degenerate");
    TestMatrixArith<T,tmv::UpperTriMatrix<T,tmv::UnitDiag>,tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> >(
	u2.View(),a3.View(),cu2.View(),c3.View(),"Tri/Degenerate");
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a3.View(),u1.View(),c3.View(),cu1.View(),"Degenerate/Tri");
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a3.View(),u2.View(),c3.View(),cu2.View(),"Degenerate/Tri");
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a3.View(),l1.View(),c3.View(),cl1.View(),"Degenerate/Tri");
    TestMatrixArith<T,tmv::Matrix<T>,tmv::Matrix<std::complex<T> > >(
	a3.View(),l2.View(),c3.View(),cl2.View(),"Degenerate/Tri");
  }
#endif

  std::cout<<"TriMatrix<"<<tmv::Type(T())<<"> (Matrix/Tri) Arithmetic passed all tests\n";
}

template void TestTriMatrixArith_B<double>();
#ifndef NOFLOAT
template void TestTriMatrixArith_B<float>();
#endif
#ifdef LONGDOUBLE
template void TestTriMatrixArith_B<long double>();
#endif

