#ifdef NDEBUG
#undef NDEBUG
#endif
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"

#define NONSQUARE
#include "TMV_TestMatrixArith.h"

template <class T> void TestSmallMatrixArith_B1()
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

  tmv::SmallMatrix<T,4,4,tmv::ColMajor> a2 = a1.Transpose();
  a2.row(1) *= T(3);
  a2.col(2) -= tmv::SmallVector<T,4>(4.);
  tmv::SmallMatrix<std::complex<T>,4,4,tmv::ColMajor> ca2 = ca1;
  ca2 -= a2;
  ca2 *= std::complex<T>(1,-2);

  tmv::SmallMatrix<T,7,4> a4;
  for(int i=0;i<7;++i) for(int j=0;j<4;++j) a4(i,j) = T(1-3*i+2*j);
  tmv::SmallMatrix<T,4,7,tmv::ColMajor> a5 = a4.Transpose();
  a4.SubMatrix(2,6,0,4) += a1;
  a5.SubMatrix(0,4,1,5) -= a2;
  tmv::SmallMatrixView<T,7,4,4,1> a4v = a4.View();

  tmv::SmallMatrix<std::complex<T>,7,4> ca4 = a4*std::complex<T>(1,2);
  ca4.SubMatrix(2,6,0,4) += ca1;
  ca4.col(1) *= std::complex<T>(2,1);
  ca4.row(6).AddToAll(std::complex<T>(-7,2));
  tmv::SmallMatrixView<std::complex<T>,7,4,4,1> ca4v = ca4.View();
  tmv::SmallMatrix<T,7,4> a4x;
  tmv::SmallMatrix<std::complex<T>,7,4> ca4x;

#ifdef XTEST
  tmv::SmallMatrix<std::complex<T>,4,7,tmv::ColMajor> ca5 = ca4.Adjoint();
  ca5.SubMatrix(0,4,1,5) -= ca2;
  ca5.col(3) *= std::complex<T>(-1,3);
  ca5.row(0).AddToAll(std::complex<T>(1,9));
  tmv::SmallMatrixView<T,4,7,1,4> a5v = a5.View();
  tmv::SmallMatrixView<std::complex<T>,4,7,1,4 > ca5v = ca5.View();
  tmv::SmallMatrix<T,4,7> a5x;
  tmv::SmallMatrix<std::complex<T>,4,7> ca5x;
#endif

  if (showstartdone) {
    std::cout<<"B1\n";
  }
  TestMatrixArith1<T>(a4x,ca4x,a4v,ca4v,"NonSquare");
#ifdef XTEST
  TestMatrixArith1<T>(a5x,ca5x,a5v,ca5v,"NonSquare");
#endif
}

#ifdef INST_DOUBLE
template void TestSmallMatrixArith_B1<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrixArith_B1<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrixArith_B1<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrixArith_B1<int>();
#endif
