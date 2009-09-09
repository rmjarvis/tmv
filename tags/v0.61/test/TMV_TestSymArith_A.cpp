#ifdef NDEBUG
#undef NDEBUG
#endif
#define START1 0
#define START2 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_TestSymArith.h"

template <class T> inline bool CanAddEq(
    const tmv::SymMatrixView<std::complex<T> >& a,
    const tmv::SymMatrixView<std::complex<T> >& b)
{ return a.size() == b.size() && a.sym() == b.sym(); }

template <class T1, class T2> inline bool CanAddEqX(
    const tmv::SymMatrixView<T1>& a, const T2 x)
{ return tmv::IsReal(x) || a.issym(); }

template <class T1, class T2> inline bool CanMultEqX(
    const tmv::SymMatrixView<T1>& a, const T2 x)
{ return tmv::IsReal(x) || a.issym(); }

#define NOMULTEQ

#include "TMV_TestMatrixArith.h"

template <class T> void TestSymMatrixArith_A()
{
  std::vector<tmv::SymMatrixView<T> > s;
  std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
  MakeSymList(s,cs,InDef);

  for(size_t i=START1;i<s.size();i++) {
    if (showstartdone) {
      std::cout<<"Start loop i = "<<i<<std::endl;
      std::cout<<"si = "<<s[i]<<std::endl;
    }
    tmv::SymMatrixView<T> si = s[i];
    tmv::SymMatrixView<std::complex<T> > csi = cs[i];

    if (cs[i].isherm()) {
      tmv::HermMatrix<T> sx = si;
      tmv::HermMatrix<std::complex<T> > csx = csi;
      TestMatrixArith123<T>(sx,csx,si,csi,"Herm");

      for(size_t j=START2;j<s.size();j++) if (i!=j) {
	if (showstartdone) {
	  std::cout<<"Start sub-loop j = "<<j<<std::endl;
	  std::cout<<"sj = "<<s[j]<<std::endl;
	}
	TestMatrixArith45<T>(sx,csx,si,csi,s[j],cs[j],"Herm/Herm");
      }
    } else {
      tmv::SymMatrix<T> sx = si;
      tmv::SymMatrix<std::complex<T> > csx = csi;
      TestMatrixArith123<T>(sx,csx,si,csi,"Sym");

      for(size_t j=START2;j<s.size();j++) if (i!=j) {
	if (showstartdone) {
	  std::cout<<"Start sub-loop j = "<<j<<std::endl;
	  std::cout<<"sj = "<<s[j]<<std::endl;
	}
	TestMatrixArith45<T>(sx,csx,si,csi,s[j],cs[j],"Sym/Sym");
      }
    }
  }
}

#ifdef INST_DOUBLE
template void TestSymMatrixArith_A<double>();
#endif
#ifdef INST_FLOAT
template void TestSymMatrixArith_A<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSymMatrixArith_A<long double>();
#endif
#ifdef INST_INT
template void TestSymMatrixArith_A<int>();
#endif
