#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"
#include <fstream>

template <class T, size_t M, size_t N, tmv::StorageType S> 
inline void DoTestSmallMatrix_Sub3()
{
  tmv::SmallMatrix<T,M,N,S> m;
  tmv::SmallMatrix<T,M,N,S,tmv::FortranStyle> mf;
  Assert(m.colsize() == size_t(M) && m.rowsize() == size_t(N),
      "Creating SmallMatrix(M,N)");
  Assert(m.colsize() == size_t(M) && m.rowsize() == size_t(N),
      "Creating SmallMatrixF(M,N)");

  for (size_t i=0, k=0; i<M; ++i) for (size_t j=0; j<N; ++j, ++k) {
    m(i,j) = k;
    mf(i+1,j+1) = k;
  }

#define Si (S==tmv::RowMajor ? int(N) : 1)
#define Sj (S==tmv::RowMajor ? 1 : int(M))

  tmv::ConstSmallMatrixView<T,M,N,Si,Sj> mcv = m.View();
  tmv::SmallMatrixView<T,M,N,Si,Sj> mv = m.View();
  tmv::ConstSmallMatrixView<T,M,N,Si,Sj,false,tmv::FortranStyle> mfcv = 
    mf.View();
  tmv::SmallMatrixView<T,M,N,Si,Sj,false,tmv::FortranStyle> mfv = mf.View();

  Assert(m.ColPair(2,5) == m.SubMatrix(0,M,2,8,1,3),"ColPair");
  Assert(m.ColPair(7,2) == m.SubMatrix(0,M,7,-3,1,-5),"ColPair");

  Assert(mf.ColPair(3,6) == mf.SubMatrix(1,M,3,6,1,3),"ColPairFF");
  Assert(mf.ColPair(8,3) == mf.SubMatrix(1,M,8,3,1,-5),"ColPairFF");

  Assert(m.ColPair(2,5) == mf.ColPair(3,6),"ColPairF");
  Assert(m.ColPair(7,2) == mf.ColPair(8,3),"ColPairF");

  Assert(m.ColPair(2,5) == mcv.ColPair(2,5),"ColPairCV");
  Assert(m.ColPair(7,2) == mcv.ColPair(7,2),"ColPairCV");

  Assert(m.ColPair(2,5) == mv.ColPair(2,5),"ColPairV");
  Assert(m.ColPair(7,2) == mv.ColPair(7,2),"ColPairV");

  Assert(mf.ColPair(3,6) == mfcv.ColPair(3,6),"ColPairFCV");
  Assert(mf.ColPair(8,3) == mfcv.ColPair(8,3),"ColPairFCV");

  Assert(mf.ColPair(3,6) == mfv.ColPair(3,6),"ColPairFV");
  Assert(mf.ColPair(8,3) == mfv.ColPair(8,3),"ColPairFV");

  tmv::SmallMatrixView<T,M,2,Si,3*Sj,false> sub6(m.SubMatrix(0,M,2,8,1,3));
  Assert(m.ColPair(2,5) == sub6,"ColPair");
  tmv::SmallMatrixView<T,M,2,Si,-5*Sj,false> sub7(m.SubMatrix(0,M,7,-3,1,-5));
  Assert(m.ColPair(7,2) == sub7,"ColPair");

  tmv::ConstSmallMatrixView<T,M,2,Si,3*Sj,false> sub17(mcv.ColPair(2,5));
  Assert(m.ColPair(2,5) == sub17,"ColPairCV");
  tmv::ConstSmallMatrixView<T,M,2,Si,-5*Sj,false> sub18(mcv.ColPair(7,2));
  Assert(m.ColPair(7,2) == sub18,"ColPairCV");

  tmv::SmallMatrixView<T,M,2,Si,3*Sj,false> sub28(mv.ColPair(2,5));
  Assert(m.ColPair(2,5) == sub28,"ColPairV");
  tmv::SmallMatrixView<T,M,2,Si,-5*Sj,false> sub29(mv.ColPair(7,2));
  Assert(m.ColPair(7,2) == sub29,"ColPairV");

  tmv::SmallMatrixView<T,M,2,Si,3*Sj,false,tmv::FortranStyle> sub39(
      mf.ColPair(3,6));
  Assert(m.ColPair(2,5) == sub39,"ColPairF");
  tmv::SmallMatrixView<T,M,2,Si,-5*Sj,false,tmv::FortranStyle> sub40(
      mf.ColPair(8,3));
  Assert(m.ColPair(7,2) == sub40,"ColPairF");

  tmv::ConstSmallMatrixView<T,M,2,Si,3*Sj,false,tmv::FortranStyle> sub50(
      mfcv.ColPair(3,6));
  Assert(m.ColPair(2,5) == sub50,"ColPairFCV");
  tmv::ConstSmallMatrixView<T,M,2,Si,-5*Sj,false,tmv::FortranStyle> sub51(
      mfcv.ColPair(8,3));
  Assert(m.ColPair(7,2) == sub51,"ColPairFCV");

  tmv::SmallMatrixView<T,M,2,Si,3*Sj,false,tmv::FortranStyle> sub61(
      mfv.ColPair(3,6));
  Assert(m.ColPair(2,5) == sub61,"ColPairFV");
  tmv::SmallMatrixView<T,M,2,Si,-5*Sj,false,tmv::FortranStyle> sub62(
      mfv.ColPair(8,3));
  Assert(m.ColPair(7,2) == sub62,"ColPairFV");

#undef Si
#undef Sj
}

template <class T> void TestSmallMatrix_Sub3()
{
  DoTestSmallMatrix_Sub3<T,15,10,tmv::RowMajor>();
  DoTestSmallMatrix_Sub3<T,15,10,tmv::ColMajor>();
}

#ifdef INST_DOUBLE
template void TestSmallMatrix_Sub3<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrix_Sub3<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrix_Sub3<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrix_Sub3<int>();
#endif
