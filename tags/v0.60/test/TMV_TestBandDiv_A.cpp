
#define START 0

#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_Band.h"
#include "TMV_TestBandArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestBandDiv_A(tmv::DivType dt)
{
  const int N = 10;

  std::vector<tmv::BandMatrixView<T> > b;
  std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
  MakeBandList(b,cb);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-2*j;
  a1.diag().AddToAll(T(10)*N);
  a1 /= T(10);
  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

  tmv::Matrix<T> a3 = a1.Cols(0,N/2);
  tmv::Matrix<std::complex<T> > ca3 = ca1.Cols(0,N/2);
  tmv::Matrix<T> a4 = a1.Rows(0,N/2);
  tmv::Matrix<std::complex<T> > ca4 = ca1.Rows(0,N/2);
  tmv::Matrix<T> a5 = a1.Cols(0,0);
  tmv::Matrix<std::complex<T> > ca5 = ca1.Cols(0,0);
  tmv::Matrix<T> a6 = a1.Rows(0,0);
  tmv::Matrix<std::complex<T> > ca6 = ca1.Rows(0,0);
  tmv::Matrix<T> a7 = a1;
  tmv::Matrix<std::complex<T> > ca7 = ca1;
  a7.diag().AddToAll(T(10)*N);
  ca7.diag().AddToAll(T(10)*N);

  for(size_t i=START;i<b.size();i++) {
    if (showstartdone) 
      std::cout<<"Start loop: i = "<<i<<"\nbi = "<<tmv::Type(b[i])<<"  "<<b[i]<<std::endl;
    const tmv::BandMatrixView<T>& bi = b[i];
    const tmv::BandMatrixView<std::complex<T> >& cbi = cb[i];
    if (dt == tmv::LU && !bi.IsSquare()) continue;

    bi.SaveDiv();
    cbi.SaveDiv();

    TestMatrixDivArith2<T>(dt,bi,a1.View(),cbi,ca1.View(),
	"SquareMatrix/Band");
    TestMatrixDivArith1<T>(dt,a7.View(),bi,ca7.View(),cbi,
	"Band/SquareMatrix");
#ifdef XTEST
    TestMatrixDivArith1<T>(dt,bi,a3.View(),cbi,ca3.View(),
	"NonSquareMatrix/Band");
    TestMatrixDivArith1<T>(dt,bi,a4.View(),cbi,ca4.View(),
	"NonSquareMatrix/Band");
    TestMatrixDivArith1<T>(dt,bi,a5.View(),cbi,ca5.View(),
	"DegenerateMatrix/Band");
    TestMatrixDivArith1<T>(dt,bi,a6.View(),cbi,ca6.View(),
	"DegenerateMatrix/Band");
#endif
  }
}

#ifdef INST_DOUBLE
template void TestBandDiv_A<double>(tmv::DivType dt);
#endif
#ifdef INST_FLOAT
template void TestBandDiv_A<float>(tmv::DivType dt);
#endif
#ifdef INST_LONGDOUBLE
template void TestBandDiv_A<long double>(tmv::DivType dt);
#endif
#ifdef INST_INT
template void TestBandDiv_A<int>(tmv::DivType dt);
#endif
