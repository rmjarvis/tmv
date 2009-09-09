#ifdef NDEBUG
#undef NDEBUG
#endif
#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_TestSymBandArith.h"

#define NOADDEQ
#define NOMULTEQ

#include "TMV_TestMatrixArith.h"

template <class T> void TestSymBandMatrixArith_B1()
{
  const int N = 10;

  std::vector<tmv::SymBandMatrixView<T> > s;
  std::vector<tmv::SymBandMatrixView<std::complex<T> > > cs;
  MakeSymBandList(s,cs,InDef);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-5*j;
  tmv::Matrix<std::complex<T> > ca1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) =
    std::complex<T>(3.+i-5.*j,2.-3.*i);

  tmv::MatrixView<T> a1v = a1.View();
  tmv::MatrixView<std::complex<T> > ca1v = ca1.View();

#ifdef XTEST
  tmv::Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+6*j;
  tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) =
    std::complex<T>(1.-3.*i+6.*j,-4.+2.*j);

  tmv::Matrix<T,tmv::RowMajor> a3 = a2.Rows(0,N);
  tmv::Matrix<std::complex<T> > ca3 = a3 * std::complex<T>(-3,4.);
  tmv::Matrix<T,tmv::RowMajor> a4 = a1.Cols(0,0);
  tmv::Matrix<std::complex<T> > ca4 = a4;

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

  for(size_t i=START;i<s.size();i++) {
    if (showstartdone) {
      std::cout<<"Start loop i = "<<i<<std::endl;
      std::cout<<"si = "<<s[i]<<std::endl;
    }

    tmv::SymBandMatrixView<T> si = s[i];
    tmv::SymBandMatrixView<std::complex<T> > csi = cs[i];

    if (csi.isherm()) {
      tmv::HermBandMatrix<T> sx = si;
      tmv::HermBandMatrix<std::complex<T> > csx = csi;

      TestMatrixArith45<T>(sx,csx,si,csi,a1v,ca1v,"HermBand/SquareM");
#ifdef XTEST
      TestMatrixArith45<T>(sx,csx,si,csi,a3v,ca3v,"HermBand/NonSquareM");
      TestMatrixArith45<T>(sx,csx,si,csi,a4v,ca4v,"HermBand/DegenerateM");
#endif
    } else {
      tmv::SymBandMatrix<T> sx = si;
      tmv::SymBandMatrix<std::complex<T> > csx = csi;

      TestMatrixArith45<T>(sx,csx,si,csi,a1v,ca1v,"SymBand/SquareM");
#ifdef XTEST
      TestMatrixArith45<T>(sx,csx,si,csi,a3v,ca3v,"SymBand/NonSquareM");
      TestMatrixArith45<T>(sx,csx,si,csi,a4v,ca4v,"SymBand/DegenerateM");
#endif
    }
  }
}

#ifdef INST_DOUBLE
template void TestSymBandMatrixArith_B1<double>();
#endif
#ifdef INST_FLOAT
template void TestSymBandMatrixArith_B1<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSymBandMatrixArith_B1<long double>();
#endif
#ifdef INST_INT
template void TestSymBandMatrixArith_B1<int>();
#endif
