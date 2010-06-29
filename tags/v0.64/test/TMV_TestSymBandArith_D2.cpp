#define START 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestSymBandArith.h"

#define NOADDEQ
#define NOMULTEQ

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymBandMatrixArith_D2()
{
#if (XTEST & 2)
    const int N = 10;

    std::vector<tmv::SymBandMatrixView<T> > sb;
    std::vector<tmv::SymBandMatrixView<std::complex<T> > > csb;
    std::vector<tmv::BaseMatrix<T>*> B;
    std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
    MakeSymBandList(sb,csb,B,CB,InDef);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) = 
        std::complex<T>(3+i-5*j,2-3*i);

    tmv::UpperTriMatrix<T,tmv::NonUnitDiag,tmv::RowMajor> u1(a1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,tmv::RowMajor> cu1(ca1);
    tmv::UpperTriMatrixView<T> u1v = u1.view();
    tmv::UpperTriMatrixView<std::complex<T> > cu1v = cu1.view();
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

    for(size_t i=START;i<sb.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<sb[i]<<std::endl;
        }

        tmv::SymBandMatrixView<T> si = sb[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];

        TestMatrixArith4<T>(u1v,cu1v,si,csi,"UpperTri/SymBand");
        TestMatrixArith5<T>(u1v,cu1v,si,csi,"UpperTri/SymBand");
        TestMatrixArith6x<T>(u1v,cu1v,si,csi,"UpperTri/SymBand");
        TestMatrixArith4<T>(l1v,cl1v,si,csi,"LowerTri/SymBand");
        TestMatrixArith5<T>(l1v,cl1v,si,csi,"LowerTri/SymBand");
        TestMatrixArith6x<T>(l1v,cl1v,si,csi,"LowerTri/SymBand");
        TestMatrixArith4<T>(u2v,cu2v,si,csi,"UpperTri/SymBand");
        TestMatrixArith5<T>(u2v,cu2v,si,csi,"UpperTri/SymBand");
        TestMatrixArith6x<T>(u2v,cu2v,si,csi,"UpperTri/SymBand");
        TestMatrixArith4<T>(l2v,cl2v,si,csi,"LowerTri/SymBand");
        TestMatrixArith5<T>(l2v,cl2v,si,csi,"LowerTri/SymBand");
        TestMatrixArith6x<T>(l2v,cl2v,si,csi,"LowerTri/SymBand");
        TestMatrixArith4<T>(u3v,cu3v,si,csi,"UpperTri/SymBand");
        TestMatrixArith5<T>(u3v,cu3v,si,csi,"UpperTri/SymBand");
        TestMatrixArith6x<T>(u3v,cu3v,si,csi,"UpperTri/SymBand");
        TestMatrixArith4<T>(l3v,cl3v,si,csi,"LowerTri/SymBand");
        TestMatrixArith5<T>(l3v,cl3v,si,csi,"LowerTri/SymBand");
        TestMatrixArith6x<T>(l3v,cl3v,si,csi,"LowerTri/SymBand");
        TestMatrixArith4<T>(u4v,cu4v,si,csi,"UpperTri/SymBand");
        TestMatrixArith5<T>(u4v,cu4v,si,csi,"UpperTri/SymBand");
        TestMatrixArith6x<T>(u4v,cu4v,si,csi,"UpperTri/SymBand");
        TestMatrixArith4<T>(l4v,cl4v,si,csi,"LowerTri/SymBand");
        TestMatrixArith5<T>(l4v,cl4v,si,csi,"LowerTri/SymBand");
        TestMatrixArith6x<T>(l4v,cl4v,si,csi,"LowerTri/SymBand");
    }
    for(size_t i=0;i<B.size();++i) delete B[i];
    for(size_t i=0;i<CB.size();++i) delete CB[i];
#endif
}

#ifdef TEST_DOUBLE
template void TestSymBandMatrixArith_D2<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymBandMatrixArith_D2<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandMatrixArith_D2<long double>();
#endif
#ifdef TEST_INT
template void TestSymBandMatrixArith_D2<int>();
#endif