// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_TestSymBandArith.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> void TestSymBandDiv_E2(tmv::DivType dt, PosDefCode pdc)
{
  const int N = 10;

  std::vector<tmv::SymBandMatrixView<T> > sb;
  std::vector<tmv::SymBandMatrixView<std::complex<T> > > csb;
  std::vector<tmv::BaseMatrix<T>*> B;
  std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
  MakeSymBandList(sb,csb,B,CB,pdc);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(1-3*i+j);
  a1.diag().AddToAll(T(10)*N);
  a1 /= T(10);
  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

  tmv::BandMatrix<T> b1(a1,1,3);
  tmv::BandMatrix<std::complex<T> > cb1(ca1,1,3);

  tmv::BandMatrix<T> b1v = b1.View();
  tmv::BandMatrix<std::complex<T> > cb1v = cb1.View();

#ifdef XTEST
  tmv::BandMatrix<T> b3(a1.Cols(0,N-2),1,3);
  tmv::BandMatrix<std::complex<T> > cb3(ca1.Cols(0,N-2),1,3);
  tmv::BandMatrix<T> b4(a1.Rows(0,N-2),1,3);
  tmv::BandMatrix<std::complex<T> > cb4(ca1.Rows(0,N-2),1,3);

  tmv::BandMatrix<T> b3v = b3.View();
  tmv::BandMatrix<std::complex<T> > cb3v = cb3.View();
  tmv::BandMatrix<T> b4v = b4.View();
  tmv::BandMatrix<std::complex<T> > cb4v = cb4.View();
#endif

  for(size_t i=START;i<sb.size();i++) {
    if (showstartdone)
      std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TypeText(sb[i])<<"  "<<sb[i]<<std::endl;
    const tmv::SymBandMatrixView<T>& si = sb[i];
    const tmv::SymBandMatrixView<std::complex<T> >& csi = csb[i];

    if (csi.issym()) {
      tmv::SymBandMatrix<T> sx = si;
      tmv::SymBandMatrix<std::complex<T> > csx = csi;

      TestMatrixDivArith1<T>(dt,sx,csx,b1v,si,cb1v,csi,
          "SymBand/SquareBandMatrix");
      if (dt == tmv::LU) continue;
#ifdef XTEST
      TestMatrixDivArith1<T>(dt,sx,csx,b3v,si,cb3v,csi,
          "SymBand/NonSquareBandMatrix");
      TestMatrixDivArith1<T>(dt,sx,csx,b4v,si,cb4v,csi,
          "SymBand/NonSquareBandMatrix");
#endif
    } else {
      tmv::HermBandMatrix<T> hx = si;
      tmv::HermBandMatrix<std::complex<T> > chx = csi;

      TestMatrixDivArith1<T>(dt,hx,chx,b1v,si,cb1v,csi,
          "SymBand/SquareBandMatrix");
      if (dt == tmv::LU) continue;
#ifdef XTEST
      TestMatrixDivArith1<T>(dt,hx,chx,b3v,si,cb3v,csi,
          "SymBand/NonSquareBandMatrix");
      TestMatrixDivArith1<T>(dt,hx,chx,b4v,si,cb4v,csi,
          "SymBand/NonSquareBandMatrix");
#endif
    }
  }
  for(size_t i=0;i<B.size();++i) delete B[i];
  for(size_t i=0;i<CB.size();++i) delete CB[i];
}

#ifdef TEST_DOUBLE
template void TestSymBandDiv_E2<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef TEST_FLOAT
template void TestSymBandDiv_E2<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandDiv_E2<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
