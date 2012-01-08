#define START 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymBandArith.h"

#define NOADDEQ
#define NOMULTEQ
#define NOELEMMULT

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymBandMatrixArith_D1()
{
    std::vector<tmv::SymBandMatrixView<T> > sb;
    std::vector<tmv::SymBandMatrixView<std::complex<T> > > csb;
    MakeSymBandList(sb,csb,InDef);

    const int N = sb[0].size();

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) = 
        std::complex<T>(3+i-5*j,2-3*i);

    tmv::UpperTriMatrix<T,tmv::NonUnitDiag|tmv::RowMajor> u1(a1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag|tmv::RowMajor> cu1(ca1);
    tmv::UpperTriMatrixView<T> u1v = u1.view();
    tmv::UpperTriMatrixView<std::complex<T> > cu1v = cu1.view();

#if (XTEST & 2)
    tmv::UpperTriMatrix<T,tmv::UnitDiag|tmv::RowMajor> u2(a1);
    tmv::UpperTriMatrix<T,tmv::NonUnitDiag|tmv::ColMajor> u3(a1);
    tmv::UpperTriMatrix<T,tmv::UnitDiag|tmv::ColMajor> u4(a1);
    tmv::LowerTriMatrix<T,tmv::NonUnitDiag|tmv::RowMajor> l1(a1);
    tmv::LowerTriMatrix<T,tmv::UnitDiag|tmv::RowMajor> l2(a1);
    tmv::LowerTriMatrix<T,tmv::NonUnitDiag|tmv::ColMajor> l3(a1);
    tmv::LowerTriMatrix<T,tmv::UnitDiag|tmv::ColMajor> l4(a1);

    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag|tmv::RowMajor> cu2(ca1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag|tmv::ColMajor> cu3(ca1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag|tmv::ColMajor> cu4(ca1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag|tmv::RowMajor> cl1(ca1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag|tmv::RowMajor> cl2(ca1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag|tmv::ColMajor> cl3(ca1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag|tmv::ColMajor> cl4(ca1);

    tmv::UpperTriMatrixView<T> u2v = u2.view();
    tmv::UpperTriMatrixView<T> u3v = u3.view();
    tmv::UpperTriMatrixView<T> u4v = u4.view();
    tmv::LowerTriMatrixView<T> l1v = l1.view();
    tmv::LowerTriMatrixView<T> l2v = l2.view();
    tmv::LowerTriMatrixView<T> l3v = l3.view();
    tmv::LowerTriMatrixView<T> l4v = l4.view();
    tmv::UpperTriMatrixView<std::complex<T> > cu2v = cu2.view();
    tmv::UpperTriMatrixView<std::complex<T> > cu3v = cu3.view();
    tmv::UpperTriMatrixView<std::complex<T> > cu4v = cu4.view();
    tmv::LowerTriMatrixView<std::complex<T> > cl1v = cl1.view();
    tmv::LowerTriMatrixView<std::complex<T> > cl2v = cl2.view();
    tmv::LowerTriMatrixView<std::complex<T> > cl3v = cl3.view();
    tmv::LowerTriMatrixView<std::complex<T> > cl4v = cl4.view();
#endif

    for(size_t i=START;i<sb.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<sb[i]<<std::endl;
        }

        tmv::SymBandMatrixView<T> si = sb[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];

        TestMatrixArith4(si,csi,u1v,cu1v,"SymBand/UpperTri");
        TestMatrixArith5(si,csi,u1v,cu1v,"SymBand/UpperTri");
        TestMatrixArith6x(si,csi,u1v,cu1v,"SymBand/UpperTri");
#if (XTEST & 2)
        TestMatrixArith4(si,csi,l1v,cl1v,"SymBand/LowerTri");
        TestMatrixArith5(si,csi,l1v,cl1v,"SymBand/LowerTri");
        TestMatrixArith6x(si,csi,l1v,cl1v,"SymBand/LowerTri");
        TestMatrixArith4(si,csi,u2v,cu2v,"SymBand/UpperTri");
        TestMatrixArith5(si,csi,u2v,cu2v,"SymBand/UpperTri");
        TestMatrixArith6x(si,csi,u2v,cu2v,"SymBand/UpperTri");
        TestMatrixArith4(si,csi,l2v,cl2v,"SymBand/LowerTri");
        TestMatrixArith5(si,csi,l2v,cl2v,"SymBand/LowerTri");
        TestMatrixArith6x(si,csi,l2v,cl2v,"SymBand/LowerTri");
        TestMatrixArith4(si,csi,u3v,cu3v,"SymBand/UpperTri");
        TestMatrixArith5(si,csi,u3v,cu3v,"SymBand/UpperTri");
        TestMatrixArith6x(si,csi,u3v,cu3v,"SymBand/UpperTri");
        TestMatrixArith4(si,csi,l3v,cl3v,"SymBand/LowerTri");
        TestMatrixArith5(si,csi,l3v,cl3v,"SymBand/LowerTri");
        TestMatrixArith6x(si,csi,l3v,cl3v,"SymBand/LowerTri");
        TestMatrixArith4(si,csi,u4v,cu4v,"SymBand/UpperTri");
        TestMatrixArith5(si,csi,u4v,cu4v,"SymBand/UpperTri");
        TestMatrixArith6x(si,csi,u4v,cu4v,"SymBand/UpperTri");
        TestMatrixArith4(si,csi,l4v,cl4v,"SymBand/LowerTri");
        TestMatrixArith5(si,csi,l4v,cl4v,"SymBand/LowerTri");
        TestMatrixArith6x(si,csi,l4v,cl4v,"SymBand/LowerTri");
#endif
    }
}

#ifdef TEST_DOUBLE
template void TestSymBandMatrixArith_D1<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymBandMatrixArith_D1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandMatrixArith_D1<long double>();
#endif
#ifdef TEST_INT
template void TestSymBandMatrixArith_D1<int>();
#endif
