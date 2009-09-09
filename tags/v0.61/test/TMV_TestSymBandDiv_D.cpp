#ifdef NDEBUG
#undef NDEBUG
#endif

#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Tri.h"
#include "TMV_TestSymBandArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestSymBandDiv_D(tmv::DivType dt, PosDefCode pdc)
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

  tmv::UpperTriMatrix<T> u(a1);
  tmv::UpperTriMatrix<std::complex<T> > cu(ca1);
  tmv::LowerTriMatrix<T> l(a1);
  tmv::LowerTriMatrix<std::complex<T> > cl(ca1);

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
	TestMatrixDivArith1<T>(dt,u.View(),si,cu.View(),csi,
	    "SymBand/UpperTriMatrix");
	TestMatrixDivArith1<T>(dt,l.View(),si,cl.View(),csi,
	    "SymBand/LowerTriMatrix");
      }
    }
    if (pdc != Sing) {
      TestMatrixDivArith1<T>(dt,si,u.View(),csi,cu.View(),
	  "UpperTriMatrix/SymBand");
      TestMatrixDivArith1<T>(dt,si,l.View(),csi,cl.View(),
	  "LowerTriMatrix/SymBand");
    }
  }
}

#ifdef INST_DOUBLE
template void TestSymBandDiv_D<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_FLOAT
template void TestSymBandDiv_D<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymBandDiv_D<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_INT
template void TestSymBandDiv_D<int>(tmv::DivType dt, PosDefCode pdc);
#endif
