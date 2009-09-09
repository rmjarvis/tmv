#ifdef NDEBUG
#undef NDEBUG
#endif
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"
#include <fstream>

template <class T, size_t M, size_t N, tmv::StorageType S> 
static void DoTestSmallMatrix_Sub1()
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

  Assert(m.SubMatrix(2,5,1,4) == m.SubMatrix(2,5,1,4,1,1),"SubMatrix");
  Assert(m.SubVector(2,5,4,2,3) == m.SubMatrix(2,14,5,11,4,2).diag(),
      "SubVector");

  Assert(mf.SubMatrix(3,5,2,4) == mf.SubMatrix(3,5,2,4,1,1),"SubMatrixFF");
  Assert(mf.SubVector(3,6,4,2,3) == mf.SubMatrix(3,11,6,10,4,2).diag(),
      "SubVectorFF");

  Assert(m.SubMatrix(2,5,1,4) == mf.SubMatrix(3,5,2,4),"SubMatrixF");
  Assert(m.SubMatrix(2,8,1,10,2,3) == mf.SubMatrix(3,7,2,8,2,3),"SubMatrixF");

  Assert(m.SubMatrix(2,5,1,4) == mcv.SubMatrix(2,5,1,4),"SubMatrixCV");
  Assert(m.SubMatrix(2,8,1,10,2,3) == mcv.SubMatrix(2,8,1,10,2,3),
      "SubMatrixCV");

  Assert(m.SubMatrix(2,5,1,4) == mv.SubMatrix(2,5,1,4),"SubMatrixV");
  Assert(m.SubMatrix(2,8,1,10,2,3) == mv.SubMatrix(2,8,1,10,2,3),"SubMatrixV");

  Assert(mf.SubMatrix(3,5,2,4) == mfcv.SubMatrix(3,5,2,4),"SubMatrixFCV");
  Assert(mf.SubMatrix(3,7,2,8,2,3) == mfcv.SubMatrix(3,7,2,8,2,3),
      "SubMatrixFCV");

  Assert(mf.SubMatrix(3,5,2,4) == mfv.SubMatrix(3,5,2,4),"SubMatrixFV");
  Assert(mf.SubMatrix(3,7,2,8,2,3) == mfv.SubMatrix(3,7,2,8,2,3),"SubMatrixFV");

  tmv::SmallMatrixView<T,3,3,Si,Sj,false> sub1(m.SubMatrix(2,5,1,4));
  Assert(m.SubMatrix(2,5,1,4) == sub1,"SubMatrix");
  tmv::SmallMatrixView<T,3,3,2*Si,3*Sj,false> sub2(m.SubMatrix(2,8,1,10,2,3));
  Assert(m.SubMatrix(2,8,1,10,2,3) == sub2,"SubMatrix");

  tmv::ConstSmallMatrixView<T,3,3,Si,Sj,false> sub12(mcv.SubMatrix(2,5,1,4));
  Assert(m.SubMatrix(2,5,1,4) == sub12,"SubMatrixCV");
  tmv::ConstSmallMatrixView<T,3,3,2*Si,3*Sj,false> sub13(
      mcv.SubMatrix(2,8,1,10,2,3));
  Assert(m.SubMatrix(2,8,1,10,2,3) == sub13,"SubMatrixCV");

  tmv::SmallMatrixView<T,3,3,Si,Sj,false> sub23(mv.SubMatrix(2,5,1,4));
  Assert(m.SubMatrix(2,5,1,4) == sub23,"SubMatrixV");
  tmv::SmallMatrixView<T,3,3,2*Si,3*Sj,false> sub24(mv.SubMatrix(2,8,1,10,2,3));
  Assert(m.SubMatrix(2,8,1,10,2,3) == sub24,"SubMatrixV");

  tmv::SmallMatrixView<T,3,3,Si,Sj,false,tmv::FortranStyle> sub34(
      mf.SubMatrix(3,5,2,4));
  Assert(m.SubMatrix(2,5,1,4) == sub34,"SubMatrixF");
  tmv::SmallMatrixView<T,3,3,2*Si,3*Sj,false,tmv::FortranStyle> sub35(
      mf.SubMatrix(3,7,2,8,2,3));
  Assert(m.SubMatrix(2,8,1,10,2,3) == sub35,"SubMatrixF");

  tmv::ConstSmallMatrixView<T,3,3,Si,Sj,false,tmv::FortranStyle> sub45(
      mfcv.SubMatrix(3,5,2,4));
  Assert(m.SubMatrix(2,5,1,4) == sub45,"SubMatrixFCV");
  tmv::ConstSmallMatrixView<T,3,3,2*Si,3*Sj,false,tmv::FortranStyle> sub46(
      mfcv.SubMatrix(3,7,2,8,2,3));
  Assert(m.SubMatrix(2,8,1,10,2,3) == sub46,"SubMatrixFCV");

  tmv::SmallMatrixView<T,3,3,Si,Sj,false,tmv::FortranStyle> sub56(
      mfv.SubMatrix(3,5,2,4));
  Assert(m.SubMatrix(2,5,1,4) == sub56,"SubMatrixFV");
  tmv::SmallMatrixView<T,3,3,2*Si,3*Sj,false,tmv::FortranStyle> sub57(
      mfv.SubMatrix(3,7,2,8,2,3));
  Assert(m.SubMatrix(2,8,1,10,2,3) == sub57,"SubMatrixFV");

#undef Si
#undef Sj
}

template <class T> void TestSmallMatrix_Sub1()
{
  DoTestSmallMatrix_Sub1<T,15,10,tmv::RowMajor>();
  DoTestSmallMatrix_Sub1<T,15,10,tmv::ColMajor>();
}

#ifdef INST_DOUBLE
template void TestSmallMatrix_Sub1<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrix_Sub1<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrix_Sub1<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrix_Sub1<int>();
#endif
