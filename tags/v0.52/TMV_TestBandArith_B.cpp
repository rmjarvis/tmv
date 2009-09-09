//#define SHOWTESTS
//#define DOALLARITH

#include "TMV_Test.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestMatrixArith.h"

using tmv::Matrix;
using tmv::BandMatrix;
using tmv::BandMatrixView;
using tmv::MatrixView;
using tmv::ConjItType;
using tmv::RowMajor;
using tmv::ColMajor;
using tmv::DiagMajor;

template <class T1, class T2> bool CanAddEq(
    const BandMatrixView<T1>& a, const MatrixView<T2>& b)
{ return false; }

template <class T1, class T2> bool CanMultEq(
    const BandMatrixView<T1>& a, const MatrixView<T2>& b)
{ return false; }

template <class T1, class T2> bool CanMultEq2(
    const MatrixView<T1>& a, const BandMatrixView<T2>& b)
{ return false; }

template <class T> void TestBandMatrixArith_B()
{
  const int N = 10;

  vector<BandMatrixView<T> > b;

  Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-5*j;
  Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+3*j;
  Vector<T> v1(N);
  Vector<T> v2(N-1);
  for (int i=0; i<N; ++i) v1(i) = 16.-3*i; 
  for (int i=0; i<N-1; ++i) v2(i) = -6.+i; 

  BandMatrix<T,RowMajor> B1(a1,3,1);
  b.push_back(B1.View());
  BandMatrix<T,ColMajor> B2(a1,3,1);
  b.push_back(B2.View());
  BandMatrix<T,DiagMajor> B3(a1,3,1);
  b.push_back(B3.View());
  BandMatrix<T,DiagMajor> B3b(a1,1,3);
  b.push_back(B3b.View());
  BandMatrix<T> B4(a2,6,6);
  b.push_back(B4.SubBandMatrix(0,2*N,0,N,3,3));
  b.push_back(B4.SubBandMatrix(0,N+2,0,N,4,4));
  b.push_back(B4.SubBandMatrix(0,2*N,0,2*N,3,3,2,2));
  BandMatrix<T,DiagMajor> B5(B1,1,1);
  b.push_back(B5.View());
  BandMatrix<T> B6(B1,3,0);
  b.push_back(B6.View());
  BandMatrix<T> B7(a1,0,3);
  b.push_back(B7.View());
#ifdef DOALLARITH
  BandMatrix<T,DiagMajor> B8 = tmv::UpperBiDiagMatrix(v1,v2);
  b.push_back(B8.View());
  BandMatrix<T,DiagMajor> B9 = tmv::LowerBiDiagMatrix(v2,v1);
  b.push_back(B9.View());
  BandMatrix<T,DiagMajor> B10 = tmv::TriDiagMatrix(v2,v1,v2);
  b.push_back(B10.View());
  BandMatrix<T,DiagMajor> B11 = tmv::UpperBiDiagMatrix(v1,v1);
  b.push_back(B11.View());
  BandMatrix<T,DiagMajor> B12 = tmv::LowerBiDiagMatrix(v1,v1);
  b.push_back(B12.View());
  BandMatrix<T,DiagMajor> B13 = tmv::TriDiagMatrix(v1,v1,v2);
  b.push_back(B13.View());
  BandMatrix<T,DiagMajor> B14 = tmv::TriDiagMatrix(v2,v1,v1);
  b.push_back(B14.View());
  b.push_back(B4.SubBandMatrix(0,N,0,2*N,3,3));
  b.push_back(B4.SubBandMatrix(0,N,0,N+2,4,4));
  BandMatrix<T> B15(a1,1,9);
  b.push_back(B15.View());
  BandMatrix<T> B16(a1,3,7);
  b.push_back(B16.View());
  BandMatrix<T> B17(a1,0,0);
  b.push_back(B17.View());
#endif

  Matrix<complex<T> > ca1 = a1 * complex<T>(-3,4.);
  Matrix<T,RowMajor> a3 = a2.Rows(0,N);
  Matrix<complex<T> > ca3 = a3 * complex<T>(-3,4.);
  Matrix<T,RowMajor> a4 = a1.Cols(0,0);
  Matrix<complex<T> > ca4 = a4;

  for(size_t i=0;i<b.size();i++) {
    BandMatrixView<T> bi = b[i];
    BandMatrix<complex<T> > cbi = bi * complex<T>(2.,1.);
    TestMatrixArith<T,BandMatrix<T>,BandMatrix<complex<T> > >(
	bi,a1.View(),cbi.View(),ca1.View(),"Band/SquareM");
    TestMatrixArith<T,BandMatrix<T>,BandMatrix<complex<T> > >(
	bi,a3.View(),cbi.View(),ca3.View(),"Band/NonSquareM");
    TestMatrixArith<T,BandMatrix<T>,BandMatrix<complex<T> > >(
	bi,a4.View(),cbi.View(),ca4.View(),"Band/DegenerateM");
#ifdef DOALLARITH
    TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
	a1.View(),bi,ca1.View(),cbi.View(),"SquareM/Band");
    TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
	a3.View(),bi,ca3.View(),cbi.View(),"NonSquareM/Band");
    TestMatrixArith<T,Matrix<T>,Matrix<complex<T> > >(
	a4.View(),bi,ca4.View(),cbi.View(),"DegenerateM/Band");
#endif
  }
}

template void TestBandMatrixArith_B<double>();
#ifndef NOFLOAT
template void TestBandMatrixArith_B<float>();
#endif
#ifdef LONGDOUBLE
template void TestBandMatrixArith_B<long double>();
#endif
