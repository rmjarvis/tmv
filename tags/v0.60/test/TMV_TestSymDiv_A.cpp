
#define START 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestSymArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> void TestSymDiv_A(tmv::DivType dt, PosDefCode pdc)
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

  tmv::Matrix<T> a3 = a1.Cols(0,N/2);
  tmv::Matrix<std::complex<T> > ca3 = ca1.Cols(0,N/2);
  tmv::Matrix<T> a4 = a1.Rows(0,N/2);
  tmv::Matrix<std::complex<T> > ca4 = ca1.Rows(0,N/2);
  tmv::Matrix<T> a5 = a1.Cols(0,0);
  tmv::Matrix<std::complex<T> > ca5 = ca1.Cols(0,0);
  tmv::Matrix<T> a6 = a1.Rows(0,0);
  tmv::Matrix<std::complex<T> > ca6 = ca1.Rows(0,0);

  for(size_t i=START;i<s.size();i++) {
    tmv::SymMatrixView<T> si = s[i];
    tmv::SymMatrixView<std::complex<T> > csi = cs[i];
    if (dt == tmv::CH && csi.issym()) continue;
    si.SaveDiv();
    csi.SaveDiv();
    if (showstartdone)
      std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::Type(si)<<"  "<<si<<std::endl;

    if (pdc == PosDef) {
      TestMatrixDivArith2<T>(dt==tmv::CH?tmv::LU:dt,a1.View(),si,
	  ca1.View(),csi,"Sym/SquareMatrix");
    }
    if (pdc != Sing) {
      TestMatrixDivArith1<T>(dt,si,a1.View(),csi,ca1.View(),"SquareMatrix/Sym");
#ifdef XTEST
      TestMatrixDivArith1<T>(dt,si,a3.View(),csi,ca3.View(),
	  "NonSquareMatrix/Sym");
      TestMatrixDivArith1<T>(dt,si,a4.View(),csi,ca4.View(),
	  "NonSquareMatrix/Sym");
      TestMatrixDivArith1<T>(dt,si,a5.View(),csi,ca5.View(),
	  "DegenerateMatrix/Sym");
      TestMatrixDivArith1<T>(dt,si,a6.View(),csi,ca6.View(),
	  "DegenerateMatrix/Sym");
#endif
    }
  }
}

#ifdef INST_DOUBLE
template void TestSymDiv_A<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_FLOAT
template void TestSymDiv_A<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymDiv_A<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_INT
template void TestSymDiv_A<int>(tmv::DivType dt, PosDefCode pdc);
#endif
