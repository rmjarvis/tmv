#ifdef NDEBUG
#undef NDEBUG
#endif
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"
#include <fstream>

template <class T, size_t M, size_t N, tmv::StorageType S> 
static void DoTestSmallMatrix_Sub5()
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

  Assert(m.Cols(2,5) == m.SubMatrix(0,M,2,5),"Cols");
  Assert(m.Rows(3,7) == m.SubMatrix(3,7,0,N),"Rows");

  Assert(mf.Cols(3,5) == mf.SubMatrix(1,M,3,5),"ColsFF");
  Assert(mf.Rows(4,7) == mf.SubMatrix(4,7,1,N),"RowsFF");

  Assert(m.Cols(2,5) == mf.Cols(3,5),"ColsF");
  Assert(m.Rows(3,7) == mf.Rows(4,7),"RowsF");

  Assert(m.Cols(2,5) == mcv.Cols(2,5),"ColsCV");
  Assert(m.Rows(3,7) == mcv.Rows(3,7),"RowsCV");

  Assert(m.Cols(2,5) == mv.Cols(2,5),"ColsV");
  Assert(m.Rows(3,7) == mv.Rows(3,7),"RowsV");

  Assert(mf.Cols(3,5) == mfcv.Cols(3,5),"ColsFCV");
  Assert(mf.Rows(4,7) == mfcv.Rows(4,7),"RowsFCV");

  Assert(mf.Cols(3,5) == mfv.Cols(3,5),"ColsFV");
  Assert(mf.Rows(4,7) == mfv.Rows(4,7),"RowsFV");

  tmv::SmallMatrixView<T,M,3,Si,Sj,false> sub10(m.SubMatrix(0,M,2,5));
  Assert(m.Cols(2,5) == sub10,"Cols");
  tmv::SmallMatrixView<T,4,N,Si,Sj,false> sub11(m.SubMatrix(3,7,0,N));
  Assert(m.Rows(3,7) == sub11,"Rows");

  tmv::ConstSmallMatrixView<T,M,3,Si,Sj,false> sub21(mcv.Cols(2,5));
  Assert(m.Cols(2,5) == sub21,"ColsCV");
  tmv::ConstSmallMatrixView<T,4,N,Si,Sj,false> sub22(mcv.Rows(3,7));
  Assert(m.Rows(3,7) == sub22,"RowsCV");

  tmv::SmallMatrixView<T,M,3,Si,Sj,false> sub32(mv.Cols(2,5));
  Assert(m.Cols(2,5) == sub32,"ColsV");
  tmv::SmallMatrixView<T,4,N,Si,Sj,false> sub33(mv.Rows(3,7));
  Assert(m.Rows(3,7) == sub33,"RowsV");

  tmv::SmallMatrixView<T,M,3,Si,Sj,false,tmv::FortranStyle> sub43(
      mf.Cols(3,5));
  Assert(m.Cols(2,5) == sub43,"ColsF");
  tmv::SmallMatrixView<T,4,N,Si,Sj,false,tmv::FortranStyle> sub44(
      mf.Rows(4,7));
  Assert(m.Rows(3,7) == sub44,"RowsF");

  tmv::ConstSmallMatrixView<T,M,3,Si,Sj,false,tmv::FortranStyle> sub54(
      mfcv.Cols(3,5));
  Assert(m.Cols(2,5) == sub54,"ColsFCV");
  tmv::ConstSmallMatrixView<T,4,N,Si,Sj,false,tmv::FortranStyle> sub55(
      mfcv.Rows(4,7));
  Assert(m.Rows(3,7) == sub55,"RowsFCV");

  tmv::SmallMatrixView<T,M,3,Si,Sj,false,tmv::FortranStyle> sub65(
      mfv.Cols(3,5));
  Assert(m.Cols(2,5) == sub65,"ColsFV");
  tmv::SmallMatrixView<T,4,N,Si,Sj,false,tmv::FortranStyle> sub66(
      mfv.Rows(4,7));
  Assert(m.Rows(3,7) == sub66,"RowsFV");

#undef Si
#undef Sj
}

template <class T> void TestSmallMatrix_Sub5()
{
  DoTestSmallMatrix_Sub5<T,15,10,tmv::RowMajor>();
  DoTestSmallMatrix_Sub5<T,15,10,tmv::ColMajor>();
}

#ifdef INST_DOUBLE
template void TestSmallMatrix_Sub5<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrix_Sub5<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrix_Sub5<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrix_Sub5<int>();
#endif
