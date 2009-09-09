#ifdef NDEBUG
#undef NDEBUG
#endif
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"
#include <fstream>

template <class T, size_t M, size_t N, tmv::StorageType S> 
static void DoTestSmallMatrix_Sub2()
{
  tmv::SmallMatrix<T,M,N,S> m;
  tmv::SmallMatrix<T,M,N,S,tmv::FortranStyle> mf;
  Assert(m.colsize() == size_t(M) && m.rowsize() == size_t(N),
      "Creating SmallMatrix(M,N)");
  Assert(m.colsize() == size_t(M) && m.rowsize() == size_t(N),
      "Creating SmallMatrixF(M,N)");

  for (size_t i=0, k=0; i<M; ++i) for (size_t j=0; j<N; ++j, ++k) {
    m(i,j) = T(k);
    mf(i+1,j+1) = T(k);
  }

#define Si (S==tmv::RowMajor ? int(N) : 1)
#define Sj (S==tmv::RowMajor ? 1 : int(M))

  tmv::ConstSmallMatrixView<T,M,N,Si,Sj> mcv = m.View();
  tmv::SmallMatrixView<T,M,N,Si,Sj> mv = m.View();
  tmv::ConstSmallMatrixView<T,M,N,Si,Sj,false,tmv::FortranStyle> mfcv = 
    mf.View();
  tmv::SmallMatrixView<T,M,N,Si,Sj,false,tmv::FortranStyle> mfv = mf.View();

  Assert(m.SubVector(2,5,4,2,3) == m.SubMatrix(2,14,5,11,4,2).diag(),
      "SubVector");

  tmv::ConstSmallVectorView<T,3,4*Si+2*Sj> v1 = mf.SubVector(3,6,4,2,3);
  tmv::ConstSmallVectorView<T,3,4*Si+2*Sj> v2 = mf.SubMatrix(3,11,6,10,4,2).diag();
  Assert(v1 == v2,"SubVectorFF_alt");
  Assert(mf.SubVector(3,6,4,2,3) == mf.SubMatrix(3,11,6,10,4,2).diag(),
      "SubVectorFF");

  Assert(m.SubVector(2,5,4,2,3) == mf.SubVector(3,6,4,2,3),"SubVectorF");
  Assert(m.SubVector(8,1,-1,2,4) == mf.SubVector(9,2,-1,2,4),"SubVector2F");
  Assert(m.SubVector(12,8,-4,-2,2) == mf.SubVector(13,9,-4,-2,2),
      "SubVector3F");

  Assert(m.SubVector(2,5,4,2,3) == mcv.SubVector(2,5,4,2,3),"SubVectorCV");
  Assert(m.SubVector(8,1,-1,2,4) == mcv.SubVector(8,1,-1,2,4),"SubVector2CV");
  Assert(m.SubVector(12,8,-4,-2,2) == mcv.SubVector(12,8,-4,-2,2),
      "SubVector3CV");

  Assert(m.SubVector(2,5,4,2,3) == mv.SubVector(2,5,4,2,3),"SubVectorV");
  Assert(m.SubVector(8,1,-1,2,4) == mv.SubVector(8,1,-1,2,4),"SubVector2V");
  Assert(m.SubVector(12,8,-4,-2,2) == mv.SubVector(12,8,-4,-2,2),
      "SubVector3V");

  Assert(mf.SubVector(3,6,4,2,3) == mfcv.SubVector(3,6,4,2,3),"SubVectorFCV");
  Assert(mf.SubVector(9,2,-1,2,4) == mfcv.SubVector(9,2,-1,2,4),
      "SubVector2FCV");
  Assert(mf.SubVector(13,9,-4,-2,2) == mfcv.SubVector(13,9,-4,-2,2),
      "SubVector3FCV");

  Assert(mf.SubVector(3,6,4,2,3) == mfv.SubVector(3,6,4,2,3),"SubVectorFV");
  Assert(mf.SubVector(9,2,-1,2,4) == mfv.SubVector(9,2,-1,2,4),"SubVector2FV");
  Assert(mf.SubVector(13,9,-4,-2,2) == mfv.SubVector(13,9,-4,-2,2),
      "SubVector3FV");

  tmv::SmallVectorView<T,3,4*Si+2*Sj,false> sub3(m.SubVector(2,5,4,2,3));
  Assert(m.SubVector(2,5,4,2,3) == sub3,"SubVector");
  tmv::SmallVectorView<T,4,-Si+2*Sj,false> sub4(m.SubVector(8,1,-1,2,4));
  Assert(m.SubVector(8,1,-1,2,4) == sub4,"SubVector2");
  tmv::SmallVectorView<T,2,-4*Si-2*Sj,false> sub5(m.SubVector(12,8,-4,-2,2));
  Assert(m.SubVector(12,8,-4,-2,2) == sub5,"SubVector3");

  tmv::ConstSmallVectorView<T,3,4*Si+2*Sj,false> sub14(
      mcv.SubVector(2,5,4,2,3));
  Assert(m.SubVector(2,5,4,2,3) == sub14,"SubVectorCV");
  tmv::ConstSmallVectorView<T,4,-Si+2*Sj,false> sub15(
      mcv.SubVector(8,1,-1,2,4));
  Assert(m.SubVector(8,1,-1,2,4) == sub15,"SubVector2CV");
  tmv::ConstSmallVectorView<T,2,-4*Si-2*Sj,false> sub16(
      mcv.SubVector(12,8,-4,-2,2));
  Assert(m.SubVector(12,8,-4,-2,2) == sub16,"SubVector3CV");

  tmv::SmallVectorView<T,3,4*Si+2*Sj,false> sub25(mv.SubVector(2,5,4,2,3));
  Assert(m.SubVector(2,5,4,2,3) == sub25,"SubVectorV");
  tmv::SmallVectorView<T,4,-Si+2*Sj,false> sub26(mv.SubVector(8,1,-1,2,4));
  Assert(m.SubVector(8,1,-1,2,4) == sub26,"SubVector2V");
  tmv::SmallVectorView<T,2,-4*Si-2*Sj,false> sub27(mv.SubVector(12,8,-4,-2,2));
  Assert(m.SubVector(12,8,-4,-2,2) == sub27,"SubVector3V");

  tmv::SmallVectorView<T,3,4*Si+2*Sj,false,tmv::FortranStyle> sub36(
      mf.SubVector(3,6,4,2,3));
  Assert(m.SubVector(2,5,4,2,3) == sub36,"SubVectorF");
  tmv::SmallVectorView<T,4,-Si+2*Sj,false,tmv::FortranStyle> sub37(
      mf.SubVector(9,2,-1,2,4));
  Assert(m.SubVector(8,1,-1,2,4) == sub37,"SubVector2F");
  tmv::SmallVectorView<T,2,-4*Si-2*Sj,false,tmv::FortranStyle> sub38(
      mf.SubVector(13,9,-4,-2,2));
  Assert(m.SubVector(12,8,-4,-2,2) == sub38,"SubVector3F");

  tmv::ConstSmallVectorView<T,3,4*Si+2*Sj,false,tmv::FortranStyle> sub47(
      mfcv.SubVector(3,6,4,2,3));
  Assert(m.SubVector(2,5,4,2,3) == sub47,"SubVectorFCV");
  tmv::ConstSmallVectorView<T,4,-Si+2*Sj,false,tmv::FortranStyle> sub48(
      mfcv.SubVector(9,2,-1,2,4));
  Assert(m.SubVector(8,1,-1,2,4) == sub48,"SubVector2FCV");
  tmv::ConstSmallVectorView<T,2,-4*Si-2*Sj,false,tmv::FortranStyle> sub49(
      mfcv.SubVector(13,9,-4,-2,2));
  Assert(m.SubVector(12,8,-4,-2,2) == sub49,"SubVector3FCV");

  tmv::SmallVectorView<T,3,4*Si+2*Sj,false,tmv::FortranStyle> sub58(
      mfv.SubVector(3,6,4,2,3));
  Assert(m.SubVector(2,5,4,2,3) == sub58,"SubVectorFV");
  tmv::SmallVectorView<T,4,-Si+2*Sj,false,tmv::FortranStyle> sub59(
      mfv.SubVector(9,2,-1,2,4));
  Assert(m.SubVector(8,1,-1,2,4) == sub59,"SubVector2FV");
  tmv::SmallVectorView<T,2,-4*Si-2*Sj,false,tmv::FortranStyle> sub60(
      mfv.SubVector(13,9,-4,-2,2));
  Assert(m.SubVector(12,8,-4,-2,2) == sub60,"SubVector3FV");

#undef Si
#undef Sj
}

template <class T> void TestSmallMatrix_Sub2()
{
  DoTestSmallMatrix_Sub2<T,15,10,tmv::RowMajor>();
  DoTestSmallMatrix_Sub2<T,15,10,tmv::ColMajor>();
}

#ifdef INST_DOUBLE
template void TestSmallMatrix_Sub2<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrix_Sub2<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrix_Sub2<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrix_Sub2<int>();
#endif
