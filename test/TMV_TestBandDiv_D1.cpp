// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestBandArith.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> void TestBandDiv_D1(tmv::DivType dt)
{
  const int N = 10;

  std::vector<tmv::BandMatrixView<T> > b;
  std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
  std::vector<tmv::BaseMatrix<T>*> B;
  std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
  MakeBandList(b,cb,B,CB);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(1-3*i+j);
  a1.diag().AddToAll(T(10)*N);
  a1 /= T(10);
  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

  tmv::UpperTriMatrix<T> u(a1);
  tmv::UpperTriMatrix<std::complex<T> > cu(ca1);
  tmv::LowerTriMatrix<T> l(a1);
  tmv::LowerTriMatrix<std::complex<T> > cl(ca1);

  tmv::UpperTriMatrix<T> ux = u;
  tmv::UpperTriMatrix<std::complex<T> > cux = cu;
  tmv::LowerTriMatrix<T> lx = l;
  tmv::LowerTriMatrix<std::complex<T> > clx = cl;

  for(size_t i=START;i<b.size();i++) {
    if (showstartdone) 
      std::cout<<"Start loop: i = "<<i<<"\nbi = "<<tmv::TypeText(b[i])<<"  "<<b[i]<<std::endl;
    const tmv::BandMatrixView<T>& bi = b[i];
    const tmv::BandMatrixView<std::complex<T> >& cbi = cb[i];
    if (dt == tmv::LU && !bi.IsSquare()) continue;

    bi.SaveDiv();
    cbi.SaveDiv();

    TestMatrixDivArith1<T>(dt,ux,cux,bi,u.View(),cbi,cu.View(),
        "UpperTriMatrix/Band");
    TestMatrixDivArith1<T>(dt,lx,clx,bi,l.View(),cbi,cl.View(),
        "LowerTriMatrix/Band");
  }
  for(size_t i=0;i<B.size();++i) delete B[i];
  for(size_t i=0;i<CB.size();++i) delete CB[i];
}

#ifdef TEST_DOUBLE
template void TestBandDiv_D1<double>(tmv::DivType dt);
#endif
#ifdef TEST_FLOAT
template void TestBandDiv_D1<float>(tmv::DivType dt);
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandDiv_D1<long double>(tmv::DivType dt);
#endif
