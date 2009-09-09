#ifdef NDEBUG
#undef NDEBUG
#endif

#define START 0

//#define XXDEBUG5

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestBandArith.h"
#include "TMV_Tri.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestBandDiv_D(tmv::DivType dt)
{
  const int N = 10;

  std::vector<tmv::BandMatrixView<T> > b;
  std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
  MakeBandList(b,cb);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 1.-3*i+j;
  a1.diag().AddToAll(T(10)*N);
  a1 /= T(10);
  tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

  tmv::UpperTriMatrix<T> u(a1);
  tmv::UpperTriMatrix<std::complex<T> > cu(ca1);

  tmv::LowerTriMatrix<T> l(a1);
  tmv::LowerTriMatrix<std::complex<T> > cl(ca1);

  for(size_t i=START;i<b.size();i++) {
    if (showstartdone) 
      std::cout<<"Start loop: i = "<<i<<"\nbi = "<<tmv::Type(b[i])<<"  "<<b[i]<<std::endl;
    const tmv::BandMatrixView<T>& bi = b[i];
    const tmv::BandMatrixView<std::complex<T> >& cbi = cb[i];
    if (dt == tmv::LU && !bi.IsSquare()) continue;

    bi.SaveDiv();
    cbi.SaveDiv();

    TestMatrixDivArith1<T>(dt,b[i],u.View(),cbi,cu.View(),
	"UpperTriMatrix/Band");
    TestMatrixDivArith1<T>(dt,b[i],l.View(),cbi,cl.View(),
	"LowerTriMatrix/Band");
    if (dt == tmv::LU) {
      TestMatrixDivArith1<T>(dt,u.View(),b[i],cu.View(),cbi,
	  "Band/UpperTriMatrix");
      TestMatrixDivArith1<T>(dt,l.View(),b[i],cl.View(),cbi,
	  "Band/LowerTriMatrix");
    }
  }
}

#ifdef INST_DOUBLE
template void TestBandDiv_D<double>(tmv::DivType dt);
#endif
#ifdef INST_FLOAT
template void TestBandDiv_D<float>(tmv::DivType dt);
#endif
#ifdef INST_LONGDOUBLE
template void TestBandDiv_D<long double>(tmv::DivType dt);
#endif
#ifdef INST_INT
template void TestBandDiv_D<int>(tmv::DivType dt);
#endif
