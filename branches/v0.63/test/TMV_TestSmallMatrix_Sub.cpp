#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_Small.h"
#include <fstream>

template <class T, size_t M, size_t N, tmv::StorageType S> 
static void DoTestSmallMatrix_Sub()
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

    Assert(m.SubMatrix(2,5,1,4) == m.SubMatrix(2,5,1,4,1,1),"SubMatrix");
    Assert(m.SubVector(2,5,4,2,3) == m.SubMatrix(2,14,5,11,4,2).diag(),
           "SubVector");

    Assert(mf.SubMatrix(3,5,2,4) == mf.SubMatrix(3,5,2,4,1,1),"SubMatrixFF");
    Assert(mf.SubVector(3,6,4,2,3) == mf.SubMatrix(3,11,6,10,4,2).diag(),
           "SubVectorFF");

    Assert(m.SubMatrix(2,5,1,4) == mf.SubMatrix(3,5,2,4),"SubMatrixF");
    Assert(m.SubMatrix(2,8,1,10,2,3) == mf.SubMatrix(3,7,2,8,2,3),"SubMatrixF");

    Assert(m.SubVector(2,5,4,2,3) == m.SubMatrix(2,14,5,11,4,2).diag(),
           "SubVector");

    Assert(mf.SubVector(3,6,4,2,3) == mf.SubMatrix(3,11,6,10,4,2).diag(),
           "SubVectorFF");

    Assert(m.SubVector(2,5,4,2,3) == mf.SubVector(3,6,4,2,3),"SubVectorF");
    Assert(m.SubVector(8,1,-1,2,4) == mf.SubVector(9,2,-1,2,4),"SubVector2F");
    Assert(m.SubVector(12,8,-4,-2,2) == mf.SubVector(13,9,-4,-2,2),
           "SubVector3F");

    Assert(m.ColPair(2,5) == m.SubMatrix(0,M,2,8,1,3),"ColPair");
    Assert(m.ColPair(7,2) == m.SubMatrix(0,M,7,-3,1,-5),"ColPair");

    Assert(mf.ColPair(3,6) == mf.SubMatrix(1,M,3,6,1,3),"ColPairFF");
    Assert(mf.ColPair(8,3) == mf.SubMatrix(1,M,8,3,1,-5),"ColPairFF");

    Assert(m.ColPair(2,5) == mf.ColPair(3,6),"ColPairF");
    Assert(m.ColPair(7,2) == mf.ColPair(8,3),"ColPairF");


    Assert(m.RowPair(3,7) == m.SubMatrix(3,11,0,N,4,1),"RowPair");
    Assert(m.RowPair(2,0) == m.SubMatrix(2,-2,0,N,-2,1),"RowPair");

    Assert(mf.RowPair(4,8) == mf.SubMatrix(4,8,1,N,4,1),"RowPairFF");
    Assert(mf.RowPair(3,1) == mf.SubMatrix(3,1,1,N,-2,1),"RowPairFF");

    Assert(m.RowPair(3,7) == mf.RowPair(4,8),"RowPairF");
    Assert(m.RowPair(2,0) == mf.RowPair(3,1),"RowPairF");


    Assert(m.Cols(2,5) == m.SubMatrix(0,M,2,5),"Cols");
    Assert(m.Rows(3,7) == m.SubMatrix(3,7,0,N),"Rows");

    Assert(mf.Cols(3,5) == mf.SubMatrix(1,M,3,5),"ColsFF");
    Assert(mf.Rows(4,7) == mf.SubMatrix(4,7,1,N),"RowsFF");

    Assert(m.Cols(2,5) == mf.Cols(3,5),"ColsF");
    Assert(m.Rows(3,7) == mf.Rows(4,7),"RowsF");

#undef Si
#undef Sj
}

template <class T> 
void TestSmallMatrix_Sub()
{
    DoTestSmallMatrix_Sub<T,15,10,tmv::RowMajor>();
    DoTestSmallMatrix_Sub<T,15,10,tmv::ColMajor>();
}

#ifdef INST_DOUBLE
template void TestSmallMatrix_Sub<double>();
#endif
#ifdef INST_FLOAT
template void TestSmallMatrix_Sub<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSmallMatrix_Sub<long double>();
#endif
#ifdef INST_INT
template void TestSmallMatrix_Sub<int>();
#endif
