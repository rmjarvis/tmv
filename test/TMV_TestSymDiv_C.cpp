#ifdef NDEBUG
#undef NDEBUG
#endif

#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Diag.h"
#include "TMV_TestSymArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestSymDiv_C(tmv::DivType dt, PosDefCode pdc)
{
  const int N = 10;

  std::vector<tmv::SymMatrixView<T> > s;
  std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
  MakeSymList(s,cs,pdc);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 1.-3*i+j;
  a1.diag().AddToAll(T(10)*N);
  a1 /= T(10);
  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

  tmv::DiagMatrix<T> d(a1);
  tmv::DiagMatrix<std::complex<T> > cd(ca1);

  for(size_t i=START;i<s.size();i++) {
    tmv::SymMatrixView<T> si = s[i];
    tmv::SymMatrixView<std::complex<T> > csi = cs[i];
    if (dt == tmv::CH && csi.issym()) continue;
    si.SaveDiv();
    csi.SaveDiv();
    if (showstartdone)
      std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::Type(si)<<"  "<<si<<std::endl;

    if (pdc == PosDef) {
      if (dt == tmv::LU) {
	TestMatrixDivArith1<T>(dt,d.View(),si,cd.View(),csi,
	    "Sym/DiagMatrix");
      }
    }
    if (pdc != Sing) {
      TestMatrixDivArith1<T>(dt,si,d.View(),csi,cd.View(),
	  "DiagMatrix/Sym");
    }
  }
}

#ifdef INST_DOUBLE
template void TestSymDiv_C<double>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_FLOAT
template void TestSymDiv_C<float>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymDiv_C<long double>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_INT
template void TestSymDiv_C<int>(tmv::DivType dt, PosDefCode pc);
#endif
