#ifdef NDEBUG
#undef NDEBUG
#endif
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"

#include "TMV_TestMatrixArith.h"

template <class T> void TestSmallMatrixArith_A4c()
{
  tmv::SmallMatrix<T,4,4,tmv::RowMajor> a1;
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) {
    a1(i,j) = T(2+4*i-5*j);
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

#ifdef XTEST
  tmv::SmallMatrix<T,4,4,tmv::ColMajor> a2 = a1.Transpose();
  a2.row(1) *= T(3);
  a2.col(2) -= tmv::SmallVector<T,4>(4.);
  tmv::SmallMatrix<std::complex<T>,4,4,tmv::ColMajor> ca2 = ca1;
  ca2 -= a2;
  ca2 *= std::complex<T>(1,-2);
  tmv::SmallMatrixView<T,4,4,1,4> a2v = a2.View();
  tmv::SmallMatrixView<std::complex<T>,4,4,1,4> ca2v = ca2.View();
#endif

  tmv::SmallMatrix<T,12,16> a3;
  for(int i=0;i<12;++i) for(int j=0;j<16;++j) a3(i,j) = T(1-2*i+3*j);
  a3.diag().AddToAll(T(30));
  tmv::SmallMatrix<std::complex<T>,12,16> ca3 = a3*std::complex<T>(1,-2);
  ca3.diag().AddToAll(std::complex<T>(-22,15));
  tmv::SmallMatrixView<T,4,4,48,4> a3v = a3.SubMatrix(0,12,0,16,3,4);
  tmv::SmallMatrixView<std::complex<T>,4,4,48,4> ca3v = 
    ca3.SubMatrix(0,12,0,16,3,4);
  tmv::SmallMatrix<T,4,4> a3x = a3v;
  tmv::SmallMatrix<std::complex<T>,4,4> ca3x = ca3v;

  if (showstartdone) {
    std::cout<<"A4c\n";
  }
  TestMatrixArith4<T>(a3x,ca3x,a3v,ca3v,a1v,ca1v,"Square");
#ifdef XTEST
  TestMatrixArith4<T>(a3x,ca3x,a3v,ca3v,a2v,ca2v,"Square");
#endif
}

#ifdef INST_DOUBLE
template void TestSmallMatrixArith_A4c<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrixArith_A4c<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrixArith_A4c<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrixArith_A4c<int>();
#endif
