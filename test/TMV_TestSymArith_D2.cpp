// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_TestSymArith.h"

#define NOADDEQ
#define NOMULTEQ

#include "TMV_TestMatrixArith.h"

template <class T> void TestSymMatrixArith_D2()
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

  tmv::UpperTriMatrix<T,tmv::NonUnitDiag,tmv::RowMajor> u1(a1);
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::RowMajor> cu1(ca1);
  tmv::UpperTriMatrixView<T> u1v = u1.View();
  tmv::UpperTriMatrixView<std::complex<T> > cu1v = cu1.View();
  tmv::UpperTriMatrix<T,tmv::NonUnitDiag> u1x = u1v;
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> cu1x = cu1v;

  tmv::UpperTriMatrix<T,tmv::UnitDiag,tmv::RowMajor> u2(a1);
  tmv::UpperTriMatrix<T,tmv::NonUnitDiag,tmv::ColMajor> u3(a1);
  tmv::UpperTriMatrix<T,tmv::UnitDiag,tmv::ColMajor> u4(a1);
  tmv::LowerTriMatrix<T,tmv::NonUnitDiag,tmv::RowMajor> l1(a1);
  tmv::LowerTriMatrix<T,tmv::UnitDiag,tmv::RowMajor> l2(a1);
  tmv::LowerTriMatrix<T,tmv::NonUnitDiag,tmv::ColMajor> l3(a1);
  tmv::LowerTriMatrix<T,tmv::UnitDiag,tmv::ColMajor> l4(a1);

  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::RowMajor> cu2(ca1);
  tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::ColMajor> cu3(ca1);
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::ColMajor> cu4(ca1);
  tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::RowMajor> cl1(ca1);
  tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::RowMajor> cl2(ca1);
  tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::ColMajor> cl3(ca1);
  tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag,tmv::ColMajor> cl4(ca1);

  tmv::UpperTriMatrixView<T> u2v = u2.View();
  tmv::UpperTriMatrixView<T> u3v = u3.View();
  tmv::UpperTriMatrixView<T> u4v = u4.View();
  tmv::LowerTriMatrixView<T> l1v = l1.View();
  tmv::LowerTriMatrixView<T> l2v = l2.View();
  tmv::LowerTriMatrixView<T> l3v = l3.View();
  tmv::LowerTriMatrixView<T> l4v = l4.View();
  tmv::UpperTriMatrixView<std::complex<T> > cu2v = cu2.View();
  tmv::UpperTriMatrixView<std::complex<T> > cu3v = cu3.View();
  tmv::UpperTriMatrixView<std::complex<T> > cu4v = cu4.View();
  tmv::LowerTriMatrixView<std::complex<T> > cl1v = cl1.View();
  tmv::LowerTriMatrixView<std::complex<T> > cl2v = cl2.View();
  tmv::LowerTriMatrixView<std::complex<T> > cl3v = cl3.View();
  tmv::LowerTriMatrixView<std::complex<T> > cl4v = cl4.View();

  tmv::UpperTriMatrix<T,tmv::UnitDiag> u2x = u2v;
  tmv::LowerTriMatrix<T,tmv::NonUnitDiag> l1x = l1v;
  tmv::LowerTriMatrix<T,tmv::UnitDiag> l2x = l2v;
  tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> cu2x = cu2v;
  tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag> cl1x = cl1v;
  tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cl2x = cl2v;

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

      TestMatrixArith456<T>(u1x,cu1x,u1v,cu1v,si,csi,"UpperTri/Herm");
      TestMatrixArith456<T>(l1x,cl1x,l1v,cl1v,si,csi,"LowerTri/Herm");
      TestMatrixArith456<T>(u2x,cu2x,u2v,cu2v,si,csi,"UpperTri/Herm");
      TestMatrixArith456<T>(l2x,cl2x,l2v,cl2v,si,csi,"LowerTri/Herm");
      TestMatrixArith456<T>(u1x,cu1x,u3v,cu3v,si,csi,"UpperTri/Herm");
      TestMatrixArith456<T>(l1x,cl1x,l3v,cl3v,si,csi,"LowerTri/Herm");
      TestMatrixArith456<T>(u2x,cu2x,u4v,cu4v,si,csi,"UpperTri/Herm");
      TestMatrixArith456<T>(l2x,cl2x,l4v,cl4v,si,csi,"LowerTri/Herm");
    } else {
      tmv::SymMatrix<T> sx = si;
      tmv::SymMatrix<std::complex<T> > csx = csi;

      TestMatrixArith456<T>(u1x,cu1x,u1v,cu1v,si,csi,"UpperTri/Sym");
      TestMatrixArith456<T>(l1x,cl1x,l1v,cl1v,si,csi,"LowerTri/Sym");
      TestMatrixArith456<T>(u2x,cu2x,u2v,cu2v,si,csi,"UpperTri/Sym");
      TestMatrixArith456<T>(l2x,cl2x,l2v,cl2v,si,csi,"LowerTri/Sym");
      TestMatrixArith456<T>(u1x,cu1x,u3v,cu3v,si,csi,"UpperTri/Sym");
      TestMatrixArith456<T>(l1x,cl1x,l3v,cl3v,si,csi,"LowerTri/Sym");
      TestMatrixArith456<T>(u2x,cu2x,u4v,cu4v,si,csi,"UpperTri/Sym");
      TestMatrixArith456<T>(l2x,cl2x,l4v,cl4v,si,csi,"LowerTri/Sym");
    }
  }
  for(size_t i=0;i<B.size();++i) delete B[i];
  for(size_t i=0;i<CB.size();++i) delete CB[i];
#endif
}

#ifdef INST_DOUBLE
template void TestSymMatrixArith_D2<double>();
#endif
#ifdef INST_FLOAT
template void TestSymMatrixArith_D2<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSymMatrixArith_D2<long double>();
#endif
#ifdef INST_INT
template void TestSymMatrixArith_D2<int>();
#endif
