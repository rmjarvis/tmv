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

template <class T> void TestSymBandDiv_E1(tmv::DivType dt, PosDefCode pdc)
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

  tmv::BandMatrix<T> b1x = b1;
  tmv::BandMatrix<std::complex<T> > cb1x = cb1;

#ifdef XTEST
  tmv::BandMatrix<T> b3(a1.Cols(0,N-2),1,3);
  tmv::BandMatrix<std::complex<T> > cb3(ca1.Cols(0,N-2),1,3);
  tmv::BandMatrix<T> b4(a1.Rows(0,N-2),1,3);
  tmv::BandMatrix<std::complex<T> > cb4(ca1.Rows(0,N-2),1,3);

  tmv::BandMatrix<T> b3x = b3;
  tmv::BandMatrix<T> b4x = b4;  
  tmv::BandMatrix<std::complex<T> > cb3x = cb3;
  tmv::BandMatrix<std::complex<T> > cb4x = cb4;
#endif

  for(size_t i=START;i<sb.size();i++) {
    if (showstartdone)
      std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TypeText(sb[i])<<"  "<<sb[i]<<std::endl;
    const tmv::SymBandMatrixView<T>& si = sb[i];
    const tmv::SymBandMatrixView<std::complex<T> >& csi = csb[i];
    if (dt == tmv::CH && csi.issym()) continue;
    si.SaveDiv();
    csi.SaveDiv();

    TestMatrixDivArith1<T>(dt,b1x,cb1x,si,b1.View(),csi,cb1.View(),
        "SquareBandMatrix/SymBand");
#ifdef XTEST
    TestMatrixDivArith1<T>(dt,b3x,cb3x,si,b3.View(),csi,cb3.View(),
        "SquareBandMatrix/SymBand");
    TestMatrixDivArith1<T>(dt,b4x,cb4x,si,b4.View(),csi,cb4.View(),
        "SquareBandMatrix/SymBand");
#endif
  }
  for(size_t i=0;i<B.size();++i) delete B[i];
  for(size_t i=0;i<CB.size();++i) delete CB[i];
}

#ifdef INST_DOUBLE
template void TestSymBandDiv_E1<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_FLOAT
template void TestSymBandDiv_E1<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymBandDiv_E1<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
