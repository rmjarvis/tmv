// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#define START1 0
#define START2 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_TestSymBandArith.h"

template <class T1, class T2> inline bool CanAddEq(
    const tmv::SymBandMatrixView<T1>& a, const tmv::SymBandMatrixView<T2>& b)
{ return a.size() == b.size() && a.nlo() >= b.nlo(); }

template <class T> inline bool CanAddEq(
    const tmv::SymBandMatrixView<std::complex<T> >& a, 
    const tmv::SymBandMatrixView<std::complex<T> >& b)
{ return a.size() == b.size() && a.nlo() >= b.nlo() && a.sym() == b.sym(); }

template <class T1, class T2> inline bool CanAddEqX(
    const tmv::SymBandMatrixView<T1>& a, const T2 x)
{ return tmv::IsReal(x) || !a.isherm(); }

template <class T1, class T2> inline bool CanMultEqX(
    const tmv::SymBandMatrixView<T1>& a, const T2 x)
{ return tmv::IsReal(x) || !a.isherm(); }

#define NOMULTEQ

#include "TMV_TestMatrixArith.h"

template <class T> void TestSymBandMatrixArith_A()
{
  std::vector<tmv::SymBandMatrixView<T> > s;
  std::vector<tmv::SymBandMatrixView<std::complex<T> > > cs;
  std::vector<tmv::BaseMatrix<T>*> B;
  std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
  MakeSymBandList(s,cs,B,CB,InDef);

  for(size_t i=START1;i<s.size();i++) {
    if (showstartdone) {
      std::cout<<"Start loop i = "<<i<<std::endl;
      std::cout<<"si = "<<s[i]<<std::endl;
    }
    tmv::SymBandMatrixView<T> si = s[i];
    tmv::SymBandMatrixView<std::complex<T> > csi = cs[i];

    if (cs[i].isherm()) {
      tmv::HermBandMatrix<T> sx = si;
      tmv::HermBandMatrix<std::complex<T> > csx = csi;
      TestMatrixArith123<T>(sx,csx,si,csi,"HermBand");

      for(size_t j=START2;j<s.size();j++) if (i!=j) {
        if (showstartdone) {
          std::cout<<"Start sub-loop j = "<<j<<std::endl;
          std::cout<<"sj = "<<s[j]<<std::endl;
        }
        TestMatrixArith456<T>(sx,csx,si,csi,s[j],cs[j],"HermBand/HermBand");
      }
    } else {
      tmv::SymBandMatrix<T> sx = si;
      tmv::SymBandMatrix<std::complex<T> > csx = csi;
      TestMatrixArith123<T>(sx,csx,si,csi,"SymBand");

      for(size_t j=START2;j<s.size();j++) if (i!=j) {
        if (showstartdone) {
          std::cout<<"Start sub-loop j = "<<j<<std::endl;
          std::cout<<"sj = "<<s[j]<<std::endl;
        }
        TestMatrixArith456<T>(sx,csx,si,csi,s[j],cs[j],"SymBand/SymBand");
      }
    }
  }
  for(size_t i=0;i<B.size();++i) delete B[i];
  for(size_t i=0;i<CB.size();++i) delete CB[i];
}

#ifdef TEST_DOUBLE
template void TestSymBandMatrixArith_A<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymBandMatrixArith_A<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandMatrixArith_A<long double>();
#endif
#ifdef TEST_INT
template void TestSymBandMatrixArith_A<int>();
#endif
