// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#define STARTI 0
#define STARTJ 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_TestSymBandArith.h"
#include "TMV_TestBandArith.h"

template <class T1, class T2> inline bool CanAddEq(
    const tmv::BandMatrixView<T1>& a, const tmv::SymBandMatrixView<T2>& b)
{ 
  return a.colsize() == b.size() && a.rowsize() == b.size() && 
  a.nlo() >= b.nlo() && a.nhi() >= b.nhi(); 
}

template <class T1, class T2, class T3> inline bool CanMultMM(
    const tmv::BandMatrixView<T1>& a, const tmv::SymBandMatrixView<T2>& b,
    const tmv::BandMatrixView<T3>& c)
{ 
  return a.rowsize() == b.colsize() && a.colsize() == c.colsize() &&
  b.rowsize() == c.rowsize() &&
  (c.nlo() >= a.nlo() + b.nlo() || c.nlo() == int(c.colsize())-1) &&
  (c.nhi() >= a.nhi() + b.nhi() || c.nhi() == int(c.rowsize())-1);
}

template <class T1, class T2, class T3> inline bool CanMultMM(
    const tmv::SymBandMatrixView<T1>& a, const tmv::BandMatrixView<T2>& b,
    const tmv::BandMatrixView<T3>& c)
{ 
  return a.rowsize() == b.colsize() && a.colsize() == c.colsize() &&
  b.rowsize() == c.rowsize() &&
  (c.nlo() >= a.nlo() + b.nlo() || c.nlo() == int(c.colsize())-1) &&
  (c.nhi() >= a.nhi() + b.nhi() || c.nhi() == int(c.rowsize())-1);
}

template <class T1, class T2, class T3> inline bool CanMultMM(
    const tmv::SymBandMatrixView<T1>& a, 
    const tmv::SymBandMatrixView<T2>& b, const tmv::BandMatrixView<T3>& c)
{ 
  return a.rowsize() == b.colsize() && a.colsize() == c.colsize() &&
  b.rowsize() == c.rowsize() &&
  (c.nlo() >= a.nlo() + b.nlo() || c.nlo() == int(c.colsize())-1) &&
  (c.nhi() >= a.nhi() + b.nhi() || c.nhi() == int(c.rowsize())-1);
}

#include "TMV_TestMatrixArith.h"

template <class T> void TestSymBandMatrixArith_E2()
{
#ifdef XTEST
  std::vector<tmv::SymBandMatrixView<T> > s;
  std::vector<tmv::SymBandMatrixView<std::complex<T> > > cs;
  std::vector<tmv::BaseMatrix<T>*> B;
  std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
  MakeSymBandList(s,cs,B,CB,InDef);

  std::vector<tmv::BandMatrixView<T> > b;
  std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
  MakeBandList(b,cb,B,CB);

  for(size_t i=STARTI;i<s.size();i++) {
    if (showstartdone) {
      std::cout<<"Start loop i = "<<i<<std::endl;
      std::cout<<"si = "<<s[i]<<std::endl;
    }

    tmv::SymBandMatrixView<T> si = s[i];
    tmv::SymBandMatrixView<std::complex<T> > csi = cs[i];

    for(size_t j=STARTJ;j<b.size();j++) {
      if (showstartdone) {
        std::cerr<<"Start sub-loop "<<j<<std::endl;
        std::cerr<<"bj = "<<b[j]<<std::endl;
      }
      tmv::BandMatrixView<T> bj = b[j];
      tmv::BandMatrixView<std::complex<T> > cbj = cb[j];
      tmv::BandMatrix<T> bx = bj;
      tmv::BandMatrix<std::complex<T> > cbx = cbj;

      if (csi.isherm()) {
        tmv::HermBandMatrix<T> sx = si;
        tmv::HermBandMatrix<std::complex<T> > csx = csi;
        TestMatrixArith456<T>(bx,cbx,bj,cbj,si,csi,"Band/HermBand");
      } else {
        tmv::SymBandMatrix<T> sx = si;
        tmv::SymBandMatrix<std::complex<T> > csx = csi;
        TestMatrixArith456<T>(bx,cbx,bj,cbj,si,csi,"Band/SymBand");
      }
    }
  }
  for(size_t i=0;i<B.size();++i) delete B[i];
  for(size_t i=0;i<CB.size();++i) delete CB[i];
#endif
}

#ifdef INST_DOUBLE
template void TestSymBandMatrixArith_E2<double>();
#endif
#ifdef INST_FLOAT
template void TestSymBandMatrixArith_E2<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSymBandMatrixArith_E2<long double>();
#endif
#ifdef INST_INT
template void TestSymBandMatrixArith_E2<int>();
#endif
