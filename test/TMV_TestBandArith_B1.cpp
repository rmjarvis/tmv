// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestBandArith.h"

template <class T1, class T2> inline bool CanAddEq(
    const tmv::BandMatrixView<T1>& , const tmv::MatrixView<T2>& )
{ return false; }

template <class T1, class T2, class T3> inline bool CanMultMM(
    const tmv::BandMatrixView<T1>& , const tmv::MatrixView<T2>& ,
    const tmv::BandMatrixView<T3>& )
{ return false; }

template <class T1, class T2, class T3> inline bool CanMultMM(
    const tmv::MatrixView<T2>& , const tmv::BandMatrixView<T1>& ,
    const tmv::BandMatrixView<T3>& )
{ return false; }

template <class T1, class T2, class T3> inline bool CanMultMM(
    const tmv::MatrixView<T2>& , const tmv::MatrixView<T1>& ,
    const tmv::BandMatrixView<T3>& )
{ return false; }

#define NOADDEQ
#define NOMULTEQ
#include "TMV_TestMatrixArith.h"

template <class T> void TestBandMatrixArith_B1()
{
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

  tmv::MatrixView<T> a1v = a1.View();
  tmv::MatrixView<std::complex<T> > ca1v = ca1.View();

#ifdef XTEST
  tmv::Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = T(1-3*i+3*j);
  tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j)
    ca2(i,j) = std::complex<T>(1-3*i+6*j,8+2*i-6*j);
  tmv::Matrix<T,tmv::RowMajor> a3 = a2.Rows(0,N);
  tmv::Matrix<std::complex<T> > ca3 = ca2.Rows(0,N);
  tmv::Matrix<T,tmv::RowMajor> a4 = a1.Cols(0,0);
  tmv::Matrix<std::complex<T> > ca4 = ca1.Cols(0,0);

  tmv::MatrixView<T> a3v = a3.View();
  tmv::MatrixView<T> a4v = a4.View();
  tmv::MatrixView<std::complex<T> > ca3v = ca3.View();
  tmv::MatrixView<std::complex<T> > ca4v = ca4.View();

  tmv::Matrix<T> a1x = a1;
  tmv::Matrix<T> a3x = a3;
  tmv::Matrix<T> a4x = a4;
  tmv::Matrix<std::complex<T> > ca1x = ca1;
  tmv::Matrix<std::complex<T> > ca3x = ca3;
  tmv::Matrix<std::complex<T> > ca4x = ca4;
#endif

  for(size_t i=START;i<b.size();i++) {
    if (showstartdone) {
      std::cerr<<"Start loop "<<i<<std::endl;
      std::cerr<<"bi = "<<b[i]<<std::endl;
    }
    tmv::BandMatrixView<T> bi = b[i];
    tmv::BandMatrixView<std::complex<T> > cbi = cb[i];
    tmv::Matrix<T> mx = bi;
    tmv::Matrix<std::complex<T> > cmx = cbi;

    TestMatrixArith456<T>(mx,cmx,bi,cbi,a1v,ca1v,"Band/SquareM");
#ifdef XTEST
    TestMatrixArith456<T>(mx,cmx,bi,cbi,a3v,ca3v,"Band/NonSquareM");
    TestMatrixArith456<T>(mx,cmx,bi,cbi,a4v,ca4v,"Band/DegenerateM");
#endif
  }
  for(size_t i=0;i<B.size();i++) delete B[i];
  for(size_t i=0;i<CB.size();i++) delete CB[i];
}

#ifdef INST_DOUBLE
template void TestBandMatrixArith_B1<double>();
#endif
#ifdef INST_FLOAT
template void TestBandMatrixArith_B1<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestBandMatrixArith_B1<long double>();
#endif
#ifdef INST_INT
template void TestBandMatrixArith_B1<int>();
#endif
