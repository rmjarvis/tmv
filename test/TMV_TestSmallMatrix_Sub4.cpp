#ifdef NDEBUG
#undef NDEBUG
#endif
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"
#include <fstream>

template <class T, size_t M, size_t N, tmv::StorageType S> 
static void DoTestSmallMatrix_Sub4()
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

  Assert(m.RowPair(3,7) == m.SubMatrix(3,11,0,N,4,1),"RowPair");
  Assert(m.RowPair(2,0) == m.SubMatrix(2,-2,0,N,-2,1),"RowPair");

  Assert(mf.RowPair(4,8) == mf.SubMatrix(4,8,1,N,4,1),"RowPairFF");
  Assert(mf.RowPair(3,1) == mf.SubMatrix(3,1,1,N,-2,1),"RowPairFF");

  Assert(m.RowPair(3,7) == mf.RowPair(4,8),"RowPairF");
  Assert(m.RowPair(2,0) == mf.RowPair(3,1),"RowPairF");

  Assert(m.RowPair(3,7) == mcv.RowPair(3,7),"RowPairCV");
  Assert(m.RowPair(2,0) == mcv.RowPair(2,0),"RowPairCV");

  Assert(m.RowPair(3,7) == mv.RowPair(3,7),"RowPairV");
  Assert(m.RowPair(2,0) == mv.RowPair(2,0),"RowPairV");

  Assert(mf.RowPair(4,8) == mfcv.RowPair(4,8),"RowPairFCV");
  Assert(mf.RowPair(3,1) == mfcv.RowPair(3,1),"RowPairFCV");

  Assert(mf.RowPair(4,8) == mfv.RowPair(4,8),"RowPairFV");
  Assert(mf.RowPair(3,1) == mfv.RowPair(3,1),"RowPairFV");

  tmv::SmallMatrixView<T,2,N,4*Si,Sj,false> sub8(m.SubMatrix(3,11,0,N,4,1));
  Assert(m.RowPair(3,7) == sub8,"RowPair");
  tmv::SmallMatrixView<T,2,N,-2*Si,Sj,false> sub9(m.SubMatrix(2,-2,0,N,-2,1));
  Assert(m.RowPair(2,0) == sub9,"RowPair");

  tmv::ConstSmallMatrixView<T,2,N,4*Si,Sj,false> sub19(mcv.RowPair(3,7));
  Assert(m.RowPair(3,7) == sub19,"RowPairCV");
  tmv::ConstSmallMatrixView<T,2,N,-2*Si,Sj,false> sub20(mcv.RowPair(2,0));
  Assert(m.RowPair(2,0) == sub20,"RowPairCV");

  tmv::SmallMatrixView<T,2,N,4*Si,Sj,false> sub30(mv.RowPair(3,7));
  Assert(m.RowPair(3,7) == sub30,"RowPairV");
  tmv::SmallMatrixView<T,2,N,-2*Si,Sj,false> sub31(mv.RowPair(2,0));
  Assert(m.RowPair(2,0) == sub31,"RowPairV");

  tmv::SmallMatrixView<T,2,N,4*Si,Sj,false,tmv::FortranStyle> sub41(
      mf.RowPair(4,8));
  Assert(m.RowPair(3,7) == sub41,"RowPairF");
  tmv::SmallMatrixView<T,2,N,-2*Si,Sj,false,tmv::FortranStyle> sub42(
      mf.RowPair(3,1));
  Assert(m.RowPair(2,0) == sub42,"RowPairF");

  tmv::ConstSmallMatrixView<T,2,N,4*Si,Sj,false,tmv::FortranStyle> sub52(
      mfcv.RowPair(4,8));
  Assert(m.RowPair(3,7) == sub52,"RowPairFCV");
  tmv::ConstSmallMatrixView<T,2,N,-2*Si,Sj,false,tmv::FortranStyle> sub53(
      mfcv.RowPair(3,1));
  Assert(m.RowPair(2,0) == sub53,"RowPairFCV");

  tmv::SmallMatrixView<T,2,N,4*Si,Sj,false,tmv::FortranStyle> sub63(
      mfv.RowPair(4,8));
  Assert(m.RowPair(3,7) == sub63,"RowPairFV");
  tmv::SmallMatrixView<T,2,N,-2*Si,Sj,false,tmv::FortranStyle> sub64(
      mfv.RowPair(3,1));
  Assert(m.RowPair(2,0) == sub64,"RowPairFV");

#undef Si
#undef Sj
}

template <class T> void TestSmallMatrix_Sub4()
{
  DoTestSmallMatrix_Sub4<T,15,10,tmv::RowMajor>();
  DoTestSmallMatrix_Sub4<T,15,10,tmv::ColMajor>();
}

#ifdef INST_DOUBLE
template void TestSmallMatrix_Sub4<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrix_Sub4<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrix_Sub4<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrix_Sub4<int>();
#endif
