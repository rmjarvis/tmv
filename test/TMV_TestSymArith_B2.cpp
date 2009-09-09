// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_TestSymArith.h"

#include "TMV_TestMatrixArith.h"

template <class T> void TestSymMatrixArith_B2()
{
#ifdef XTEST
  const int N = 10;

  std::vector<tmv::SymMatrixView<T> > s;
  std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
  std::vector<tmv::BaseMatrix<T>*> B;
  std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
  MakeSymList(s,cs,B,CB,InDef);

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
  tmv::Matrix<std::complex<T> > ca1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) = 
    std::complex<T>(3+i-5*j,2-3*i);

  tmv::Vector<T> v = a1.col(0);
  tmv::Vector<std::complex<T> > cv = ca1.col(0);

  tmv::MatrixView<T> a1v = a1.View();
  tmv::MatrixView<std::complex<T> > ca1v = ca1.View();

  tmv::Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = T(1-3*i+6*j);
  tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) = 
    std::complex<T>(1-3*i+6*j,-4+2*j);

  tmv::Vector<T> v5(5*N);
  tmv::VectorView<T> vs = v5.SubVector(0,5*N,5);
  vs = v;
  tmv::Vector<std::complex<T> > cv5(5*N);
  tmv::VectorView<std::complex<T> > cvs = cv5.SubVector(0,5*N,5);
  cvs = cv;

  tmv::Matrix<T,tmv::RowMajor> a3 = a2.Rows(0,N);
  tmv::Matrix<std::complex<T> > ca3 = a3 * std::complex<T>(-3,4);
  tmv::Matrix<T,tmv::RowMajor> a4 = a1.Cols(0,0);
  tmv::Matrix<std::complex<T> > ca4 = a4;

  tmv::MatrixView<T> a3v = a3.View();
  tmv::MatrixView<T> a4v = a4.View();
  tmv::MatrixView<std::complex<T> > ca3v = ca3.View();
  tmv::MatrixView<std::complex<T> > ca4v = ca4.View();

  tmv::Matrix<T> a1x = a1;
  tmv::Matrix<T> a3x = a3;
  tmv::Matrix<T> a4x = a4;
  tmv::Matrix<std::complex<T> > ca1x = ca1;
  tmv::Matrix<std::complex<T> > ca3x = ca3;
  tmv::Matrix<std::complex<T> > ca4x = ca4;

  symoprod = true;

  for(size_t i=START;i<s.size();i++) {
    if (showstartdone) {
      std::cout<<"Start loop i = "<<i<<std::endl;
      std::cout<<"si = "<<s[i]<<std::endl;
    }

    tmv::SymMatrixView<T> si = s[i];
    tmv::SymMatrixView<std::complex<T> > csi = cs[i];

    if (csi.isherm()) {
      tmv::HermMatrix<T> sx = si;
      tmv::HermMatrix<std::complex<T> > csx = csi;

      TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,si,csi,"SquareM/Herm");
      TestMatrixArith456<T>(a3x,ca3x,a3v,ca3v,si,csi,"NonSquareM/Herm");
      TestMatrixArith456<T>(a4x,ca4x,a4v,ca4v,si,csi,"DegenerateM/Herm");
    } else {
      tmv::SymMatrix<T> sx = si;
      tmv::SymMatrix<std::complex<T> > csx = csi;

      TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,si,csi,"SquareM/Sym");
      TestMatrixArith456<T>(a3x,ca3x,a3v,ca3v,si,csi,"NonSquareM/Sym");
      TestMatrixArith456<T>(a4x,ca4x,a4v,ca4v,si,csi,"DegenerateM/Sym");
    }
  }
  for(size_t i=0;i<B.size();++i) delete B[i];
  for(size_t i=0;i<CB.size();++i) delete CB[i];
#endif
}

#ifdef TEST_DOUBLE
template void TestSymMatrixArith_B2<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymMatrixArith_B2<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymMatrixArith_B2<long double>();
#endif
#ifdef TEST_INT
template void TestSymMatrixArith_B2<int>();
#endif
