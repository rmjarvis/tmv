#ifdef NDEBUG
#undef NDEBUG
#endif

#define START 0

//#define XXDEBUG5

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestBandArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestBandDiv_B(tmv::DivType dt)
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

  tmv::BandMatrix<T> b1(a1,3,1);
  tmv::BandMatrix<std::complex<T> > cb1 = b1 * std::complex<T>(2.,-3.);
  tmv::BandMatrix<T> b2 = b1.SubBandMatrix(0,N/2,0,N,b1.nlo(),b1.nhi());
  tmv::BandMatrix<std::complex<T> > cb2 = cb1.SubBandMatrix(
      0,N/2,0,N,cb1.nlo(),cb1.nhi());
  tmv::BandMatrix<T> b3 = b1.SubBandMatrix(0,N,0,N/2,b1.nlo(),b1.nhi());
  tmv::BandMatrix<std::complex<T> > cb3 = cb1.SubBandMatrix(
      0,N,0,N/2,cb1.nlo(),cb1.nhi());

  for(size_t i=START;i<b.size();i++) {
    if (showstartdone) 
      std::cout<<"Start loop: i = "<<i<<"\nbi = "<<tmv::Type(b[i])<<"  "<<b[i]<<std::endl;
    const tmv::BandMatrixView<T>& bi = b[i];
    const tmv::BandMatrixView<std::complex<T> >& cbi = cb[i];
    if (dt == tmv::LU && !bi.IsSquare()) continue;

    bi.SaveDiv();
    cbi.SaveDiv();

    TestMatrixDivArith1<T>(dt,bi,b1.View(),cbi,cb1.View(),"SquareBand/Band");
#ifdef XTEST
    TestMatrixDivArith1<T>(dt,bi,b2.View(),cbi,cb2.View(),"NonSquareBand/Band");
    TestMatrixDivArith1<T>(dt,bi,b3.View(),cbi,cb3.View(),"NonSquareBand/Band");
#endif
  }
}

#ifdef INST_DOUBLE
template void TestBandDiv_B<double>(tmv::DivType dt);
#endif
#ifdef INST_FLOAT
template void TestBandDiv_B<float>(tmv::DivType dt);
#endif
#ifdef INST_LONGDOUBLE
template void TestBandDiv_B<long double>(tmv::DivType dt);
#endif
#ifdef INST_INT
template void TestBandDiv_B<int>(tmv::DivType dt);
#endif
