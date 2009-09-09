// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_TestSymArith.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> void TestSymDiv_D2(tmv::DivType dt, PosDefCode pdc)
{
  const int N = 10;

  std::vector<tmv::SymMatrixView<T> > s;
  std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
  std::vector<tmv::BaseMatrix<T>*> B;
  std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
  MakeSymList(s,cs,B,CB,pdc);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(1-3*i+j);
  a1.diag().AddToAll(T(10)*N);
  a1 /= T(10);
  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

  tmv::UpperTriMatrix<T> u(a1);
  tmv::UpperTriMatrix<std::complex<T> > cu(ca1);
  tmv::LowerTriMatrix<T> l(a1);
  tmv::LowerTriMatrix<std::complex<T> > cl(ca1);

  for(size_t i=START;i<s.size();i++) {
    if (showstartdone)
      std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TypeText(s[i])<<"  "<<s[i]<<std::endl;

    const tmv::SymMatrixView<T>& si = s[i];
    const tmv::SymMatrixView<std::complex<T> >& csi = cs[i];

    if (csi.issym()) {
      tmv::SymMatrix<T> sx = si;
      tmv::SymMatrix<std::complex<T> > csx = csi;

      TestMatrixDivArith1<T>(dt,sx,csx,u.View(),si,cu.View(),csi,
          "Sym/UpperTriMatrix");
      TestMatrixDivArith1<T>(dt,sx,csx,l.View(),si,cl.View(),csi,
          "Sym/LowerTriMatrix");
    } else {
      tmv::HermMatrix<T> hx = si;
      tmv::HermMatrix<std::complex<T> > chx = csi;

      TestMatrixDivArith1<T>(dt,hx,chx,u.View(),si,cu.View(),csi,
          "Herm/UpperTriMatrix");
      TestMatrixDivArith1<T>(dt,hx,chx,l.View(),si,cl.View(),csi,
          "Herm/LowerTriMatrix");
    }
  }
  for(size_t i=0;i<B.size();++i) delete B[i];
  for(size_t i=0;i<CB.size();++i) delete CB[i];
}

#ifdef INST_DOUBLE
template void TestSymDiv_D2<double>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_FLOAT
template void TestSymDiv_D2<float>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymDiv_D2<long double>(tmv::DivType dt, PosDefCode pc);
#endif
