
#define START 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestSymBandArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestSymBandDiv_F(tmv::DivType dt, PosDefCode pdc)
{
  const int N = 10;

  std::vector<tmv::SymBandMatrixView<T> > sb;
  std::vector<tmv::SymBandMatrixView<std::complex<T> > > csb;
  MakeSymBandList(sb,csb,pdc);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 1.-3*i+j;
  a1.diag().AddToAll(T(10)*N);
  a1 /= T(10);
  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

  size_t ntot = sb.size();

  tmv::HermMatrix<T> h(a1);
  tmv::HermMatrix<std::complex<T> > ch(ca1);
  tmv::SymMatrix<T> s(a1);
  tmv::SymMatrix<std::complex<T> > cs(ca1);

  for(size_t i=START;i<ntot;i++) {
    tmv::SymBandMatrixView<T> si = sb[i];
    tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];
    if (dt == tmv::CH && csi.issym()) continue;
    si.SaveDiv();
    csi.SaveDiv();
    if (showstartdone)
      std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::Type(si)<<"  "<<si<<std::endl;

    if (pdc == PosDef) {
      if (dt == tmv::LU) {
	TestMatrixDivArith1<T>(dt,h.View(),si,ch.View(),csi,
	    "SymBand/HermMatrix");
	TestMatrixDivArith1<T>(dt,s.View(),si,cs.View(),csi,
	    "SymBand/SymMatrix");
      }
    }
    if (pdc != Sing) {
      TestMatrixDivArith1<T>(dt,si,h.View(),csi,ch.View(),
	  "HermMatrix/SymBand");
      TestMatrixDivArith1<T>(dt,si,s.View(),csi,cs.View(),
	  "SymMatrix/SymBand");
    }
  }
}

#ifdef INST_DOUBLE
template void TestSymBandDiv_F<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_FLOAT
template void TestSymBandDiv_F<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymBandDiv_F<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_INT
template void TestSymBandDiv_F<int>(tmv::DivType dt, PosDefCode pdc);
#endif
