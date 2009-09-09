
#define START 0

#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestSymArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestSymDiv_D(tmv::DivType dt, PosDefCode pdc)
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

  tmv::UpperTriMatrix<T> u(a1);
  tmv::UpperTriMatrix<std::complex<T> > cu(ca1);
  tmv::LowerTriMatrix<T> l(a1);
  tmv::LowerTriMatrix<std::complex<T> > cl(ca1);

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
	TestMatrixDivArith1<T>(dt,u.View(),si,cu.View(),csi,
	    "Sym/UpperTriMatrix");
	TestMatrixDivArith1<T>(dt,l.View(),si,cl.View(),csi,
	    "Sym/LowerTriMatrix");
      }
    }
    if (pdc != Sing) {
      TestMatrixDivArith1<T>(dt,si,u.View(),csi,cu.View(),
	  "UpperTriMatrix/Sym");
      TestMatrixDivArith1<T>(dt,si,l.View(),csi,cl.View(),
	  "LowerTriMatrix/Sym");
    }
  }
}

#ifdef INST_DOUBLE
template void TestSymDiv_D<double>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_FLOAT
template void TestSymDiv_D<float>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymDiv_D<long double>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_INT
template void TestSymDiv_D<int>(tmv::DivType dt, PosDefCode pc);
#endif
