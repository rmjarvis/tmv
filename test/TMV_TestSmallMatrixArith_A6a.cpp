#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"

#include "TMV_TestMatrixArith.h"

template <class T> void TestSmallMatrixArith_A6a()
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

  tmv::SmallVector<T,4> v1 = a1.col(0);
  tmv::SmallVector<T,4> v2 = a1.row(2);
  tmv::SmallVectorView<T,4,1> v1v = v1.View();
  tmv::SmallVectorView<T,4,1> v2v = v2.View();
  tmv::SmallVector<std::complex<T>,4> cv1 = ca1.col(1);
  tmv::SmallVector<std::complex<T>,4> cv2 = ca1.row(3);
  tmv::SmallVectorView<std::complex<T>,4,1> cv1v = cv1.View();
  tmv::SmallVectorView<std::complex<T>,4,1> cv2v = cv2.View();
#ifdef XTEST
  tmv::SmallVector<T,20> v15;
  tmv::SmallVector<T,20> v25;
  tmv::SmallVector<std::complex<T>,20> cv15;
  tmv::SmallVector<std::complex<T>,20> cv25;
  tmv::SmallVectorView<T,4,5> v1s = v15.SubVector(0,20,5);
  tmv::SmallVectorView<T,4,5> v2s = v25.SubVector(0,20,5);
  tmv::SmallVectorView<std::complex<T>,4,5> cv1s = cv15.SubVector(0,20,5);
  tmv::SmallVectorView<std::complex<T>,4,5> cv2s = cv25.SubVector(0,20,5);
#endif

  if (showstartdone) {
    std::cout<<"A6a\n";
  }
  TestMatrixArith6<T>(a1v,ca1v,v1v,cv1v,v2v,cv2v,"Square");
#ifdef XTEST
  TestMatrixArith6<T>(a1v,ca1v,v1s,cv1s,v2v,cv2v,"Square");
  TestMatrixArith6<T>(a1v,ca1v,v1v,cv1v,v2s,cv2s,"Square");
  TestMatrixArith6<T>(a1v,ca1v,v1s,cv1s,v2s,cv2s,"Square");
#endif
}

#ifdef INST_DOUBLE
template void TestSmallMatrixArith_A6a<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrixArith_A6a<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrixArith_A6a<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrixArith_A6a<int>();
#endif
