
#define START 0

//#define XXDEBUG5

#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_Band.h"
#include "TMV_TestBandArith.h"
#include "TMV_Diag.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestBandDiv_C(tmv::DivType dt)
{
  const int N = 10;

  std::vector<tmv::BandMatrixView<T> > b;
  std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
  MakeBandList(b,cb);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 1.-3*i+j;
  a1.diag().AddToAll(T(10)*N);
  a1 /= T(10);
  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

  tmv::DiagMatrix<T> d(a1);
  tmv::DiagMatrix<std::complex<T> > cd(ca1);

  for(size_t i=START;i<b.size();i++) {
    if (showstartdone) 
      std::cout<<"Start loop: i = "<<i<<"\nbi = "<<tmv::Type(b[i])<<"  "<<b[i]<<std::endl;
    const tmv::BandMatrixView<T>& bi = b[i];
    const tmv::BandMatrixView<std::complex<T> >& cbi = cb[i];
    if (dt == tmv::LU && !bi.IsSquare()) continue;

    bi.SaveDiv();
    cbi.SaveDiv();

    TestMatrixDivArith1<T>(dt,bi,d.View(),cbi,cd.View(),"DiagMatrix/Band");
    if (dt == tmv::LU) 
      TestMatrixDivArith1<T>(dt,d.View(),bi,cd.View(),cbi,"Band/DiagMatrix");
  }
}

#ifdef INST_DOUBLE
template void TestBandDiv_C<double>(tmv::DivType dt);
#endif
#ifdef INST_FLOAT
template void TestBandDiv_C<float>(tmv::DivType dt);
#endif
#ifdef INST_LONGDOUBLE
template void TestBandDiv_C<long double>(tmv::DivType dt);
#endif
#ifdef INST_INT
template void TestBandDiv_C<int>(tmv::DivType dt);
#endif
