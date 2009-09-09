#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"

#include "TMV_TestMatrixArith.h"

template <class T> void TestSmallMatrixArith_A4b()
{
  tmv::SmallMatrix<T,4,4,tmv::RowMajor> a1;
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) {
    a1(i,j) = 2+4*i-5*j;
  }
  a1(0,0) = 14.;
  a1(1,0) = -2.;
  a1(2,0) = 7.;
  a1(3,0) = -10.;
  a1(2,2) = 30.;

  tmv::SmallMatrix<std::complex<T>,4,4> ca1 = a1;
  ca1(2,3) += std::complex<T>(2,3);
  ca1(1,0) *= std::complex<T>(0,2);
  ca1.col(1) *= std::complex<T>(-1,3);
  ca1.row(3) += tmv::SmallVector<std::complex<T>,4>(std::complex<T>(1,9));
  tmv::SmallMatrixView<T,4,4,4,1> a1v = a1.View();
  tmv::SmallMatrixView<std::complex<T>,4,4,4,1> ca1v = ca1.View();

  tmv::SmallMatrix<T,4,4,tmv::ColMajor> a2 = a1.Transpose();
  a2.row(1) *= T(3);
  a2.col(2) -= tmv::SmallVector<T,4>(4.);
  tmv::SmallMatrix<std::complex<T>,4,4,tmv::ColMajor> ca2 = ca1;
  ca2 -= a2;
  ca2 *= std::complex<T>(1,-2);
  tmv::SmallMatrixView<T,4,4,1,4> a2v = a2.View();
  tmv::SmallMatrixView<std::complex<T>,4,4,1,4> ca2v = ca2.View();
  tmv::SmallMatrix<T,4,4> a2x = a2;
  tmv::SmallMatrix<std::complex<T>,4,4> ca2x = ca2;

#ifdef XTEST
  tmv::SmallMatrix<T,12,16> a3x;
  for(int i=0;i<12;++i) for(int j=0;j<16;++j) a3x(i,j) = 1-2*i+3*j;
  a3x.diag().AddToAll(30);
  tmv::SmallMatrix<std::complex<T>,12,16> ca3x = a3x*std::complex<T>(1,-2);
  ca3x.diag().AddToAll(std::complex<T>(-22,15));
  tmv::SmallMatrixView<T,4,4,48,4> a3v = a3x.SubMatrix(0,12,0,16,3,4);
  tmv::SmallMatrixView<std::complex<T>,4,4,48,4> ca3v = 
    ca3x.SubMatrix(0,12,0,16,3,4);
#endif

  if (showstartdone) {
    std::cout<<"A4b\n";
  }
  TestMatrixArith4<T>(a2x,ca2x,a2v,ca2v,a1v,ca1v,"Square");
#ifdef XTEST
  TestMatrixArith4<T>(a2x,ca2x,a2v,ca2v,a3v,ca3v,"Square");
#endif
}

#ifdef INST_DOUBLE
template void TestSmallMatrixArith_A4b<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrixArith_A4b<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrixArith_A4b<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrixArith_A4b<int>();
#endif
