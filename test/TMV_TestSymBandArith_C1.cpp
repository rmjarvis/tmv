// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_TestSymBandArith.h"

template <class T1, class T2> inline bool CanAddEq(
    const tmv::SymBandMatrixView<T1>& a, const tmv::DiagMatrixView<T2>& )
{ return a.issym() || tmv::IsReal(T2()); }

#define NOMULTEQ

#include "TMV_TestMatrixArith.h"

template <class T> void TestSymBandMatrixArith_C1()
{
  const int N = 10;

  std::vector<tmv::SymBandMatrixView<T> > s;
  std::vector<tmv::SymBandMatrixView<std::complex<T> > > cs;
  std::vector<tmv::BaseMatrix<T>*> B;
  std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
  MakeSymBandList(s,cs,B,CB,InDef);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
  tmv::Matrix<std::complex<T> > ca1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) =
    std::complex<T>(3+i-5*j,2-3*i);

  tmv::DiagMatrix<T> d1(a1);
  tmv::DiagMatrix<std::complex<T> > cd1(ca1);
  tmv::DiagMatrixView<T> d1v = d1.View();
  tmv::DiagMatrixView<std::complex<T> > cd1v = cd1.View();
  tmv::DiagMatrix<T> d1x = d1v;
  tmv::DiagMatrix<std::complex<T> > cd1x = cd1v;

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
      TestMatrixArith456<T>(sx,csx,si,csi,d1v,cd1v,"HermBand/Diag");
    } else {
      tmv::SymBandMatrix<T> sx = si;
      tmv::SymBandMatrix<std::complex<T> > csx = csi;
      TestMatrixArith456<T>(sx,csx,si,csi,d1v,cd1v,"SymBand/Diag");
    }
  }
  for(size_t i=0;i<B.size();++i) delete B[i];
  for(size_t i=0;i<CB.size();++i) delete CB[i];
}

#ifdef TEST_DOUBLE
template void TestSymBandMatrixArith_C1<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymBandMatrixArith_C1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandMatrixArith_C1<long double>();
#endif
#ifdef TEST_INT
template void TestSymBandMatrixArith_C1<int>();
#endif
