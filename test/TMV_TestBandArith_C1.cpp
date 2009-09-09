#ifdef NDEBUG
#undef NDEBUG
#endif
#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Diag.h"
#include "TMV_TestBandArith.h"

#include "TMV_TestMatrixArith.h"

template <class T> void TestBandMatrixArith_C1()
{
  const int N = 10;

  std::vector<tmv::BandMatrixView<T> > b;
  std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
  MakeBandList(b,cb);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-5*j;
  tmv::Matrix<std::complex<T> > ca1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j)
    ca1(i,j) = std::complex<T>(3.+i-5*j,4.-8*i-j);

  tmv::DiagMatrix<T> d1(a1);
  tmv::DiagMatrix<std::complex<T> > cd1(ca1);
  tmv::DiagMatrixView<T> d1v = d1.View();
  tmv::DiagMatrixView<std::complex<T> > cd1v = cd1.View();
  tmv::DiagMatrix<T> d1x = d1v;
  tmv::DiagMatrix<std::complex<T> > cd1x = cd1v;

  for(size_t i=START;i<b.size();i++) {
    if (showstartdone) {
      std::cerr<<"Start loop "<<i<<std::endl;
      std::cerr<<"bi = "<<b[i]<<std::endl;
    }
    tmv::BandMatrixView<T> bi = b[i];
    tmv::BandMatrixView<std::complex<T> > cbi = cb[i];
    tmv::BandMatrix<T> bx = bi;
    tmv::BandMatrix<std::complex<T> > cbx = cbi;

    TestMatrixArith45<T>(bx,cbx,bi,cbi,d1v,cd1v,"Band/Diag");
  }
}

#ifdef INST_DOUBLE
template void TestBandMatrixArith_C1<double>();
#endif
#ifdef INST_FLOAT
template void TestBandMatrixArith_C1<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestBandMatrixArith_C1<long double>();
#endif
#ifdef INST_INT
template void TestBandMatrixArith_C1<int>();
#endif
