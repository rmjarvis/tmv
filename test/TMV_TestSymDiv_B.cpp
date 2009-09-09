#ifdef NDEBUG
#undef NDEBUG
#endif

#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_TestSymArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestSymDiv_B(tmv::DivType dt, PosDefCode pdc)
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

  tmv::SymMatrix<T> s1(a1);
  tmv::SymMatrix<std::complex<T> > cs1(ca1);
  tmv::HermMatrix<T> h1(a1);
  tmv::HermMatrix<std::complex<T> > ch1(ca1);

  for(size_t i=START;i<s.size();i++) {
    tmv::SymMatrixView<T> si = s[i];
    tmv::SymMatrixView<std::complex<T> > csi = cs[i];
    if (dt == tmv::CH && csi.issym()) continue;
    si.SaveDiv();
    csi.SaveDiv();
    if (showstartdone)
      std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::Type(si)<<"  "<<si<<std::endl;

    if (pdc != Sing) {
      TestMatrixDivArith2<T>(dt,si,s1.View(),csi,cs1.View(),"Sym/Sym");
      TestMatrixDivArith2<T>(dt,si,h1.View(),csi,ch1.View(),"Sym/Sym");
    }
  }
}

#ifdef INST_DOUBLE
template void TestSymDiv_B<double>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_FLOAT
template void TestSymDiv_B<float>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymDiv_B<long double>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_INT
template void TestSymDiv_B<int>(tmv::DivType dt, PosDefCode pc);
#endif
