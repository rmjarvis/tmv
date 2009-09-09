//#define SHOWTESTS

#include "TMV_Test.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestMatrixArith.h"

using tmv::Matrix;
using tmv::BandMatrix;
using tmv::RowMajor;

template <class T> extern void TestBandMatrixArith_A();
template <class T> extern void TestBandMatrixArith_B();

template <class T> void TestBandMatrix()
{
  const size_t N = 10;
  const int nhi = 1;
  const int nlo = 3;

  BandMatrix<T> a1(N,N,nlo,nhi);

  Assert(a1.colsize() == N && a1.rowsize() == N,"Creating BandMatrix(N)");
  Assert(a1.nlo() == nlo && a1.nhi() == nhi,"Creating BandMatrix(nlo,nhi)");

  for (size_t i=0, k=0; i<N; ++i) for (size_t j=0; j<N; ++j, ++k)
    if ( j <= i + nhi && i <= j + nlo) a1(i,j) = k;

  for (size_t i=0, k=0; i<N; ++i) for (size_t j=0; j<N; ++j, ++k)
    if ( j <= i + nhi && i <= j + nlo) 
      Assert(a1(i,j) == k,"Read/Write BandMatrix");

  BandMatrix<T> a2(N,N,nlo,nhi);
  for (size_t i=0; i<N; ++i) for (size_t j=0; j<N; ++j) 
    if ( j <= i + nhi && i <= j + nlo) {
      a1(i,j) = 3.+i-5*j;
      a2(i,j) = 5.-2*i+4*j;
    }

  BandMatrix<T> c(N,N,nlo,nhi);
  c = a1+a2;
  for (size_t i=0; i<N; ++i) for (size_t j=0; j<N; ++j) 
    if ( j <= i + nhi && i <= j + nlo) 
      Assert(c(i,j) == 8.-i-j,"Add BandMatrices");

  c = a1-a2;
  for (size_t i=0; i<N; ++i) for (size_t j=0; j<N; ++j) 
    if ( j <= i + nhi && i <= j + nlo) 
      Assert(c(i,j) == -2.+3.*i-9.*j,"Subtract BandMatrices");

  Matrix<T> m1 = a1;
  for (size_t i=0, k=0; i<N; ++i) for (size_t j=0; j<N; ++j, ++k)
    if ( j <= i + nhi && i <= j + nlo) 
      Assert(a1(i,j) == m1(i,j),"BandMatrix -> Matrix");
  Assert(a1 == BandMatrix<T>(m1,nlo,nhi),"Matrix -> BandMatrix");

  cout<<"BandMatrix<"<<tmv::Type(T())<<"> passed all tests\n";

  TestBandMatrixArith_A<T>();
  cout<<"BandMatrix<"<<tmv::Type(T())<<"> (Band/Band) Arithmetic passed all tests\n";
  TestBandMatrixArith_B<T>();
  cout<<"BandMatrix<"<<tmv::Type(T())<<"> (Matrix/Band) Arithmetic passed all tests\n";
}

template void TestBandMatrix<double>();
#ifndef NOFLOAT
template void TestBandMatrix<float>();
#endif
#ifdef LONGDOUBLE
template void TestBandMatrix<long double>();
#endif
