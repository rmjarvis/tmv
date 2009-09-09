#define START1 0
#define START2 0

#include "TMV_Test.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestMatrixArith.h"

template <class T1, tmv::IndexStyle I1, class T2, tmv::IndexStyle I2> 
inline bool CanAddEq(
    const tmv::BandMatrixView<T1,I1>& a, const tmv::BandMatrixView<T2,I2>& b)
{ 
  return a.colsize() == b.colsize() && a.rowsize() == b.rowsize() &&
    a.nhi() >= b.nhi() && a.nlo() >= b.nlo();
}

template <class T1, tmv::IndexStyle I1, class T2, tmv::IndexStyle I2> 
inline bool CanMultEq(
    const tmv::BandMatrixView<T1,I1>& a, const tmv::BandMatrixView<T2,I2>& b)
{ 
  return a.rowsize() == b.colsize() && b.IsSquare() && 
    (b.nhi()==int(0) || a.nhi() == int(a.rowsize())-1) && 
    (b.nlo()==int(0) || a.nlo() == int(a.colsize())-1);
}

template <class T1, tmv::IndexStyle I1, class T2, tmv::IndexStyle I2> 
inline bool CanMultEq2(
    const tmv::BandMatrixView<T1,I1>& a, const tmv::BandMatrixView<T2,I2>& b)
{ 
  return a.rowsize() == b.colsize() && a.IsSquare() && 
    (a.nhi()==int(0) || b.nhi() == int(b.rowsize())-1) && 
    (a.nlo()==int(0) || b.nlo() == int(b.colsize())-1);
}

template <class T> void TestBandMatrixArith_A()
{
  const int N = 8;

  std::vector<tmv::BandMatrixView<T> > b;

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-5*j;
  tmv::Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+6*j;
  tmv::Vector<T> v1(N);
  tmv::Vector<T> v2(N-1);
  for (int i=0; i<N; ++i) v1(i) = 16.-3*i; 
  for (int i=0; i<N-1; ++i) v2(i) = -6.+i; 

  tmv::BandMatrix<T,tmv::RowMajor> B1(a1,3,1);
  b.push_back(B1.View());
  tmv::BandMatrix<T,tmv::DiagMajor> B3(a1,3,1);
  b.push_back(B3.View());
  tmv::BandMatrix<T,tmv::DiagMajor> B5(B1,1,1);
  b.push_back(B5.View());
  tmv::BandMatrix<T> B6(B1,3,0);
  b.push_back(B6.View());
#ifdef XTEST
  tmv::BandMatrix<T> B4(a2,6,6);
  b.push_back(B4.SubBandMatrix(0,2*N,0,N,3,3));
  tmv::BandMatrix<T> B4a(B4);
  b.push_back(B4a.SubBandMatrix(0,N+2,0,N,4,4));
  tmv::BandMatrix<T> B4b(B4);
  b.push_back(B4b.SubBandMatrix(0,2*N,0,2*N,3,3,2,2));
  tmv::BandMatrix<T,tmv::ColMajor> B2(a1,3,1);
  tmv::BandMatrix<T,tmv::DiagMajor> B3b(a1,1,3);
  tmv::BandMatrix<T> B7(a1,0,3);
  tmv::BandMatrix<T,tmv::DiagMajor> B8 = tmv::UpperBiDiagMatrix(v1,v2);
  tmv::BandMatrix<T,tmv::DiagMajor> B9 = tmv::LowerBiDiagMatrix(v2,v1);
  tmv::BandMatrix<T,tmv::DiagMajor> B10 = tmv::TriDiagMatrix(v2,v1,v2);
  tmv::BandMatrix<T,tmv::DiagMajor> B11 = tmv::UpperBiDiagMatrix(v1,v1);
  tmv::BandMatrix<T,tmv::DiagMajor> B12 = tmv::LowerBiDiagMatrix(v1,v1);
  tmv::BandMatrix<T,tmv::DiagMajor> B13 = tmv::TriDiagMatrix(v1,v1,v2);
  tmv::BandMatrix<T,tmv::DiagMajor> B14 = tmv::TriDiagMatrix(v2,v1,v1);
  tmv::BandMatrix<T> B4c(B4);
  tmv::BandMatrix<T> B4d(B4);
  tmv::BandMatrix<T,tmv::DiagMajor> B15(a1,1,N-1);
  tmv::BandMatrix<T,tmv::DiagMajor> B16(a1,3,N-2);
  tmv::BandMatrix<T> B17(a1,0,0);
  if (doallarith) {
    b.push_back(B2.View());
    b.push_back(B3b.View());
    b.push_back(B7.View());
    b.push_back(B8.View());
    b.push_back(B9.View());
    b.push_back(B10.View());
    b.push_back(B11.View());
    b.push_back(B12.View());
    b.push_back(B13.View());
    b.push_back(B14.View());
    b.push_back(B4c.SubBandMatrix(0,N,0,2*N,3,3));
    b.push_back(B4d.SubBandMatrix(0,N,0,N+2,4,4));
    b.push_back(BandMatrixViewOf(a1,3,N-1));
    b.push_back(B15.View());
    b.push_back(B16.View());
    b.push_back(B17.View());
  }
#endif

  for(size_t i=START1;i<b.size();i++) 
    for(size_t j=START2;j<b.size();j++) if (i!=j) {
      if (showstartdone) {
	std::cerr<<"Start loop "<<i<<','<<j<<std::endl;
	std::cerr<<"bi = "<<b[i]<<std::endl;
	std::cerr<<"bj = "<<b[j]<<std::endl;
      }
      tmv::BandMatrixView<T> bi = b[i];
      tmv::BandMatrixView<T> bj = b[j];
      tmv::BandMatrix<std::complex<T> > cbi = bi * std::complex<T>(2.,1.);
      tmv::BandMatrix<std::complex<T> > cbj = bj * std::complex<T>(-4.,3.);
      TestMatrixArith<T,tmv::BandMatrix<T>,tmv::BandMatrix<std::complex<T> > >(
	  bi,bj,cbi.View(),cbj.View(),"Band/Band");
    }
}

template void TestBandMatrixArith_A<double>();
#ifndef NOFLOAT
template void TestBandMatrixArith_A<float>();
#endif
#ifdef LONGDOUBLE
template void TestBandMatrixArith_A<long double>();
#endif
