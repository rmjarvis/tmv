// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestBandArith.h"

#define NOADDEQ

template <class T1, class T2, class T3> inline bool CanMultMM(
    const tmv::DiagMatrixView<T1>& a, const tmv::BandMatrixView<T2>& b,
    const tmv::DiagMatrixView<T3>& c)
{ 
  return b.colsize() == a.size() && b.rowsize() == a.size() && 
  c.size() == a.size() && b.nlo() == 0 && b.nhi() == 0;
}

template <class T1, class T2, class T3> inline bool CanMultMM(
    const tmv::BandMatrixView<T2>& a, const tmv::DiagMatrixView<T1>& b,
    const tmv::DiagMatrixView<T3>& c)
{ 
  return a.colsize() == b.size() && a.rowsize() == b.size() &&
  c.size() == a.size() && a.nlo() == 0 && a.nhi() == 0;
}

template <class T1, class T2, class T3> inline bool CanMultMM(
    const tmv::BandMatrixView<T1>& a, const tmv::BandMatrixView<T2>& b,
    const tmv::DiagMatrixView<T3>& c)
{ 
  return a.colsize() == c.size() && b.rowsize() == c.size() && 
  a.rowsize() == b.colsize() && 
  a.nlo() == 0 && a.nhi() == 0 && b.nlo() == 0 && b.nhi() == 0;
}

#include "TMV_TestMatrixArith.h"

template <class T> void TestBandMatrixArith_C2()
{
#ifdef XTEST
  const int N = 10;

  std::vector<tmv::BandMatrixView<T> > b;
  std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
  std::vector<tmv::BaseMatrix<T>*> B;
  std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
  MakeBandList(b,cb,B,CB);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
  tmv::Matrix<std::complex<T> > ca1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j)
    ca1(i,j) = std::complex<T>(3+i-5*j,4-8*i-j);

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

    TestMatrixArith456<T>(d1x,cd1x,d1v,cd1v,bi,cbi,"Diag/Band");
  }
  for(size_t i=0;i<B.size();i++) delete B[i];
  for(size_t i=0;i<CB.size();i++) delete CB[i];
#endif
}

#ifdef TEST_DOUBLE
template void TestBandMatrixArith_C2<double>();
#endif
#ifdef TEST_FLOAT
template void TestBandMatrixArith_C2<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandMatrixArith_C2<long double>();
#endif
#ifdef TEST_INT
template void TestBandMatrixArith_C2<int>();
#endif
