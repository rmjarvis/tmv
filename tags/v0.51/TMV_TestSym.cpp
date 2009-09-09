//#define SHOWTESTS

#include "TMV_Test.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_TestMatrixArith.h"

using tmv::Matrix;
using tmv::SymMatrix;
using tmv::HermMatrix;

template <class T> extern void TestSymMatrixArith_A();
template <class T> extern void TestSymMatrixArith_B();

template <class T, tmv::UpLoType U, tmv::StorageType S> void TestBasicSymMatrix()
{
  const size_t N = 6;

  SymMatrix<T,U,S> s1(N);
  SymMatrix<T,U,S> s2(N);
  HermMatrix<T,U,S> h1(N);
  HermMatrix<T,U,S> h2(N);

  Assert(s1.colsize() == N && s1.rowsize() == N,"Creating SymMatrix(N)");
  Assert(h1.colsize() == N && h1.rowsize() == N,"Creating HermMatrix(N)");

  for (size_t i=0, k=0; i<N; ++i) for (size_t j=0; j<N; ++j, ++k) {
    if (i<=j) { s1(i,j) = k; h1(i,j) = k; }
    if (j<=i) { s2(i,j) = k; h2(i,j) = k; }
  }
  for (size_t i=0, k=0; i<N; ++i) for (size_t j=0; j<N; ++j, ++k) {
    if (i<=j) {
      Assert(s1(i,j) == k,"Read/Write SymMatrix");
      Assert(h1(i,j) == k,"Read/Write SymMatrix");
      Assert(s1(j,i) == k,"Read/Write SymMatrix");
      Assert(h1(j,i) == k,"Read/Write SymMatrix");
    }
    if (j<=i) {
      Assert(s2(i,j) == k,"Read/Write SymMatrix");
      Assert(h2(i,j) == k,"Read/Write SymMatrix");
      Assert(s2(j,i) == k,"Read/Write SymMatrix");
      Assert(h2(j,i) == k,"Read/Write SymMatrix");
    }
  }

  SymMatrix<complex<T>,U,S> cs1(N);
  SymMatrix<complex<T>,U,S> cs2(N);
  HermMatrix<complex<T>,U,S> ch1(N);
  HermMatrix<complex<T>,U,S> ch2(N);

  for (size_t i=0; i<N; ++i) for (size_t j=0; j<N; ++j) {
    if (i<=j) {
      s1(i,j) = 3.+i-5*j;
      s2(i,j) = 5.-2*i+3*j;
      cs1(i,j) = complex<T>(3.+i-5*j,-2.*i+2.*j);
      cs2(i,j) = complex<T>(5.-2*i+3*j,1.*i-j);
      h1(i,j) = 3.+i-5*j;
      h2(i,j) = 5.-2*i+3*j;
      ch1(i,j) = complex<T>(3.+i-5*j,-2.*i+2.*j);
      ch2(i,j) = complex<T>(5.-2*i+3*j,1.*i-j);
    }
  }

  SymMatrix<T,U,S> s3(N);
  s3 = s1+s2;

  HermMatrix<T,U,S> h3(N);
  h3 = h1+h2;

  SymMatrix<complex<T>,U,S> cs3(N);
  cs3 = cs1+cs2;

  HermMatrix<complex<T>,U,S> ch3(N);
  ch3 = ch1+ch2;

  for (size_t i=0; i<N; ++i) for (size_t j=0; j<N; ++j) 
    if (i<=j) {
      Assert(s3(i,j) == 8.-i-2.*j,"Add SymMatrices1");
      Assert(h3(i,j) == 8.-i-2.*j,"Add HermMatrices1");
      Assert(s3(j,i) == 8.-i-2.*j,"Add SymMatrices2");
      Assert(h3(j,i) == 8.-i-2.*j,"Add HermMatrices2");
      Assert(cs3(i,j) == complex<T>(8.-i-2.*j,-1.*i+j),"Add CSymMatrices1");
      Assert(ch3(i,j) == complex<T>(8.-i-2.*j,-1.*i+j),"Add CHermMatrices1");
      Assert(cs3(j,i) == complex<T>(8.-i-2.*j,-1.*i+j),"Add CSymMatrices2");
      Assert(ch3(j,i) == complex<T>(8.-i-2.*j,1.*i-j),"Add CHermMatrices2");
    }

  s3 = s1-s2;
  h3 = h1-h2;
  cs3 = cs1-cs2;
  ch3 = ch1-ch2;

  for (size_t i=0; i<N; ++i) for (size_t j=0; j<N; ++j) 
    if (i<=j) {
      Assert(s3(i,j) == -2.+3.*i-8.*j,"Subtract SymMatrices1");
      Assert(h3(i,j) == -2.+3.*i-8.*j,"Subtract HermMatrices1");
      Assert(s3(j,i) == -2.+3.*i-8.*j,"Subtract SymMatrices2");
      Assert(h3(j,i) == -2.+3.*i-8.*j,"Subtract HermMatrices2");
      Assert(cs3(i,j) == complex<T>(-2.+3.*i-8.*j,-3.*i+3.*j),
	  "Subtract CSymMatrices1");
      Assert(ch3(i,j) == complex<T>(-2.+3.*i-8.*j,-3.*i+3.*j),
	  "Subtract CHermMatrices1");
      Assert(cs3(j,i) == complex<T>(-2.+3.*i-8.*j,-3.*i+3.*j),
	  "Subtract CSymMatrices2");
      Assert(ch3(j,i) == complex<T>(-2.+3.*i-8.*j,3.*i-3.*j),
	  "Subtract CHermMatrices2");
    }

  Matrix<T> m1 = s1;
  Matrix<T> n1 = h1;
  for (size_t i=0, k=0; i<N; ++i) for (size_t j=0; j<N; ++j, ++k) {
    Assert(s1(i,j) == m1(i,j),"SymMatrix -> Matrix");
    Assert(h1(i,j) == n1(i,j),"HermMatrix -> Matrix");
  }
  Assert(s1 == SymMatrix<T,U,S>(m1),"Matrix -> SymMatrix");

  SymMatrix<complex<T>,tmv::Upper,tmv::RowMajor> csur = cs1;
  Assert(cs1==csur,"SymMatrix == SymMatrix<U,R>");
  Assert(cs1.View()==csur.View(),"SymMatrix.View == SymMatrix<U,R>.View");
  Assert(cs1.Transpose()==csur.Transpose(),
      "SymMatrix.Transpose == SymMatrix<U,R>.Transpose");
  Assert(cs1.Conjugate()==csur.Conjugate(),
      "SymMatrix.Conjugate == SymMatrix<U,R>.Conjugate");
  Assert(cs1.Adjoint()==csur.Adjoint(),
      "SymMatrix.Adjoint == SymMatrix<U,R>.Adjoint");
  Assert(cs1.UpperTri()==csur.UpperTri(),
      "SymMatrix.UpperTri == SymMatrix<U,R>.UpperTri");
  Assert(cs1.LowerTri()==csur.LowerTri(),
      "SymMatrix.LowerTri == SymMatrix<U,R>.LowerTri");
  Assert(cs1.Real()==csur.Real(),"SymMatrix.Real == SymMatrix<U,R>.Real");
  Assert(cs1.Imag()==csur.Imag(),"SymMatrix.Imag == SymMatrix<U,R>.Imag");
  Assert(cs1.SubMatrix(N/2,N,0,N/2)==csur.SubMatrix(N/2,N,0,N/2),
      "SymMatrix.SubMatrix1 == SymMatrix<U,R>.SubMatrix1");
  Assert(cs1.SubMatrix(0,N/2,N/2,N)==csur.SubMatrix(0,N/2,N/2,N),
      "SymMatrix.SubMatrix2 == SymMatrix<U,R>.SubMatrix2");
  Assert(cs1.SubSymMatrix(0,N/2)==csur.SubSymMatrix(0,N/2),
      "SymMatrix.SubSymMatrix == SymMatrix<U,R>.SubSymMatrix");

  HermMatrix<complex<T>,tmv::Upper,tmv::RowMajor> chur = ch1;
  Assert(ch1==chur,"HermMatrix == HermMatrix<U,R>");
  Assert(ch1.View()==chur.View(),"HermMatrix.View == HermMatrix<U,R>.View");
  Assert(ch1.Transpose()==chur.Transpose(),
      "HermMatrix.Transpose == HermMatrix<U,R>.Transpose");
  Assert(ch1.Conjugate()==chur.Conjugate(),
      "HermMatrix.Conjugate == HermMatrix<U,R>.Conjugate");
  Assert(ch1.Adjoint()==chur.Adjoint(),
      "HermMatrix.Adjoint == HermMatrix<U,R>.Adjoint");
  Assert(ch1.UpperTri()==chur.UpperTri(),
      "HermMatrix.UpperTri == HermMatrix<U,R>.UpperTri");
  Assert(ch1.LowerTri()==chur.LowerTri(),
      "HermMatrix.LowerTri == HermMatrix<U,R>.LowerTri");
  Assert(ch1.Real()==chur.Real(),"HermMatrix.Real == HermMatrix<U,R>.Real");
  Assert(ch1.SubMatrix(N/2,N,0,N/2)==chur.SubMatrix(N/2,N,0,N/2),
      "HermMatrix.SubMatrix1 == HermMatrix<U,R>.SubMatrix1");
  Assert(ch1.SubMatrix(0,N/2,N/2,N)==chur.SubMatrix(0,N/2,N/2,N),
      "HermMatrix.SubMatrix2 == HermMatrix<U,R>.SubMatrix2");
  Assert(ch1.SubSymMatrix(0,N/2)==chur.SubSymMatrix(0,N/2),
      "HermMatrix.SubSymMatrix == HermMatrix<U,R>.SubSymMatrix");
}

template <class T> void TestSymMatrix() 
{
  TestBasicSymMatrix<T,tmv::Upper,tmv::RowMajor>();
  TestBasicSymMatrix<T,tmv::Lower,tmv::RowMajor>();
  TestBasicSymMatrix<T,tmv::Upper,tmv::ColMajor>();
  TestBasicSymMatrix<T,tmv::Lower,tmv::ColMajor>();

  cout<<"SymMatrix<"<<tmv::Type(T())<<"> passed all basic tests\n";

  TestSymMatrixArith_A<T>();
  cout<<"SymMatrix<"<<tmv::Type(T())<<"> (Sym/Sym) Arithmetic passed all tests\n";
  TestSymMatrixArith_B<T>();
  cout<<"SymMatrix<"<<tmv::Type(T())<<"> (Matrix/Sym) Arithmetic passed all tests\n";
}

template void TestSymMatrix<double>();
#ifndef NOFLOAT
template void TestSymMatrix<float>();
#endif
#ifdef LONGDOUBLE
template void TestSymMatrix<long double>();
#endif
