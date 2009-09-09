// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestBandArith.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> void TestBandDiv_B2(tmv::DivType dt)
{
  const int N = 10;

  std::vector<tmv::BandMatrixView<T> > b;
  std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
  std::vector<tmv::BaseMatrix<T>*> B;
  std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
  MakeBandList(b,cb,B,CB);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-2*j);
  a1.diag().AddToAll(T(10)*N);
  a1 /= T(10);
  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);
  a1.diag().AddToAll(T(10)*N);
  ca1.diag().AddToAll(T(10)*N);

  tmv::MatrixView<T> a1v = a1.View();
  tmv::MatrixView<std::complex<T> > ca1v = ca1.View();

#ifdef XTEST
  tmv::Matrix<T> a3 = a1.Cols(0,N/2);
  tmv::Matrix<std::complex<T> > ca3 = ca1.Cols(0,N/2);
  tmv::Matrix<T> a4 = a1.Rows(0,N/2);
  tmv::Matrix<std::complex<T> > ca4 = ca1.Rows(0,N/2);
  tmv::Matrix<T> a5(2*N,N);
  a5.Rows(0,N) = a1;
  a5.Rows(N,2*N) = a1;
  tmv::Matrix<std::complex<T> > ca5(2*N,N);
  ca5.Rows(0,N) = ca1;
  ca5.Rows(N,2*N) = ca1;
  tmv::Matrix<T> a6 = a5.Transpose();
  tmv::Matrix<std::complex<T> > ca6 = ca5.Transpose();

  tmv::MatrixView<T> a3v = a3.View();
  tmv::MatrixView<T> a4v = a4.View();
  tmv::MatrixView<T> a5v = a5.View();
  tmv::MatrixView<T> a6v = a6.View();
  tmv::MatrixView<std::complex<T> > ca3v = ca3.View();
  tmv::MatrixView<std::complex<T> > ca4v = ca4.View();
  tmv::MatrixView<std::complex<T> > ca5v = ca5.View();
  tmv::MatrixView<std::complex<T> > ca6v = ca6.View();
#endif

  for(size_t i=START;i<b.size();i++) {
    if (showstartdone) 
      std::cout<<"Start loop: i = "<<i<<"\nbi = "<<tmv::TypeText(b[i])<<"  "<<b[i]<<std::endl;
    const tmv::BandMatrixView<T>& bi = b[i];
    const tmv::BandMatrixView<std::complex<T> >& cbi = cb[i];

    tmv::BandMatrix<T> bx = bi;
    tmv::BandMatrix<std::complex<T> > cbx = cbi;

    TestMatrixDivArith1<T>(dt,bx,cbx,a1v,bi,ca1v,cbi, "Band/SquareMatrix");
    if (dt == tmv::LU) continue;
#ifdef XTEST
    TestMatrixDivArith1<T>(dt,bx,cbx,a3v,bi,ca3v,cbi, "Band/NonSquareMatrix");
    TestMatrixDivArith1<T>(dt,bx,cbx,a4v,bi,ca4v,cbi, "Band/NonSquareMatrix");
    TestMatrixDivArith1<T>(dt,bx,cbx,a5v,bi,ca5v,cbi, "Band/NonSquareMatrix");
    TestMatrixDivArith1<T>(dt,bx,cbx,a6v,bi,ca6v,cbi, "Band/NonSquareMatrix");
#endif
  }
  for(size_t i=0;i<B.size();++i) delete B[i];
  for(size_t i=0;i<CB.size();++i) delete CB[i];
}

#ifdef TEST_DOUBLE
template void TestBandDiv_B2<double>(tmv::DivType dt);
#endif
#ifdef TEST_FLOAT
template void TestBandDiv_B2<float>(tmv::DivType dt);
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandDiv_B2<long double>(tmv::DivType dt);
#endif
