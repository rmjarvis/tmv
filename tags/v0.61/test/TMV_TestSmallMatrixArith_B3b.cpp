#ifdef NDEBUG
#undef NDEBUG
#endif
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"

#include "TMV_TestMatrixArith.h"

template <class T> void TestSmallMatrixArith_B3b()
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
  a5.SubMatrix(0,4,1,5) -= a2;
  tmv::SmallMatrixView<T,4,7,1,4> a5v = a5.View();

  tmv::SmallMatrix<std::complex<T>,7,4> ca4 = a4*std::complex<T>(1,2);
  tmv::SmallMatrix<std::complex<T>,4,7,tmv::ColMajor> ca5 = ca4.Adjoint();
  ca5.SubMatrix(0,4,1,5) -= ca2;
  ca5.col(3) *= std::complex<T>(-1,3);
  ca5.row(0).AddToAll(std::complex<T>(1,9));
  tmv::SmallMatrixView<std::complex<T>,4,7,1,4 > ca5v = ca5.View();

  tmv::SmallVector<T,4> v1 = a1.row(2);
  tmv::SmallVectorView<T,4,1> v1v = v1.View();
  tmv::SmallVector<std::complex<T>,4> cv1 = ca1.row(2);
  tmv::SmallVectorView<std::complex<T>,4,1> cv1v = cv1.View();
  tmv::SmallVector<T,7> v2 = a5.row(2);
  tmv::SmallVectorView<T,7,1> v2v = v2.View();
  tmv::SmallVector<std::complex<T>,7> cv2 = ca5.row(2);
  tmv::SmallVectorView<std::complex<T>,7,1> cv2v = cv2.View();
#ifdef XTEST
  tmv::SmallVector<T,20> v15;
  tmv::SmallVector<std::complex<T>,20> cv15;
  tmv::SmallVector<T,35> v25;
  tmv::SmallVector<std::complex<T>,35> cv25;
  tmv::SmallVectorView<T,4,5> v1s = v15.SubVector(0,20,5);
  tmv::SmallVectorView<std::complex<T>,4,5> cv1s = cv15.SubVector(0,20,5);
  tmv::SmallVectorView<T,7,5> v2s = v25.SubVector(0,35,5);
  tmv::SmallVectorView<std::complex<T>,7,5> cv2s = cv25.SubVector(0,35,5);
  v1s = v1v;
  cv1s = cv1v;
  v2s = v2v;
  cv2s = cv2v;
#endif

  if (showstartdone) {
    std::cout<<"B3\n";
  }
  TestMatrixArith3<T>(a5v,ca5v,v2v,cv2v,v1v,cv1v,"NonSquare");
#ifdef XTEST
  TestMatrixArith3<T>(a5v,ca5v,v2s,cv2s,v1v,cv1v,"NonSquare");
  TestMatrixArith3<T>(a5v,ca5v,v2v,cv2v,v1v,cv1s,"NonSquare");
  TestMatrixArith3<T>(a5v,ca5v,v2s,cv2s,v1v,cv1s,"NonSquare");
#endif
}

#ifdef INST_DOUBLE
template void TestSmallMatrixArith_B3b<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrixArith_B3b<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrixArith_B3b<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrixArith_B3b<int>();
#endif
