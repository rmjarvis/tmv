#define START 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestSymBandArith.h"

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymBandMatrixArith_B2()
{
#ifdef XTEST
    const int N = 10;

    std::vector<tmv::SymBandMatrixView<T> > s;
    std::vector<tmv::SymBandMatrixView<std::complex<T> > > cs;
    std::vector<tmv::BaseMatrix<T>*> B;
    std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
    MakeSymBandList(s,cs,B,CB,InDef);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) =
        std::complex<T>(3+i-5*j,2-3*i);

    tmv::MatrixView<T> a1v = a1.view();
    tmv::MatrixView<std::complex<T> > ca1v = ca1.view();

    tmv::Matrix<T> a2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = T(1-3*i+6*j);
    tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) =
        std::complex<T>(1-3*i+6*j,-4+2*j);

    tmv::Matrix<T,tmv::RowMajor> a3 = a2.rowRange(0,N);
    tmv::Matrix<std::complex<T> > ca3 = a3 * std::complex<T>(-3,4);
    tmv::Matrix<T,tmv::RowMajor> a4 = a1.colRange(0,0);
    tmv::Matrix<std::complex<T> > ca4 = a4;

    tmv::MatrixView<T> a3v = a3.view();
    tmv::MatrixView<T> a4v = a4.view();
    tmv::MatrixView<std::complex<T> > ca3v = ca3.view();
    tmv::MatrixView<std::complex<T> > ca4v = ca4.view();

    tmv::Matrix<T> a1x = a1;
    tmv::Matrix<T> a3x = a3;
    tmv::Matrix<T> a4x = a4;
    tmv::Matrix<std::complex<T> > ca1x = ca1;
    tmv::Matrix<std::complex<T> > ca3x = ca3;
    tmv::Matrix<std::complex<T> > ca4x = ca4;

    for(size_t i=START;i<s.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<s[i]<<std::endl;
        }

        tmv::SymBandMatrixView<T> si = s[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = cs[i];

        if (csi.isherm()) {
            tmv::HermBandMatrix<T> sx = si;
            tmv::HermBandMatrix<std::complex<T> > csx = csi;

            TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,si,csi,"SquareM/HermBand");
            TestMatrixArith456<T>(a3x,ca3x,a3v,ca3v,si,csi,
                                  "NonSquareM/HermBand");
            TestMatrixArith456<T>(a4x,ca4x,a4v,ca4v,si,csi,
                                  "DegenerateM/HermBand");
        } else {
            tmv::SymBandMatrix<T> sx = si;
            tmv::SymBandMatrix<std::complex<T> > csx = csi;

            TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,si,csi,"SquareM/SymBand");
            TestMatrixArith456<T>(a3x,ca3x,a3v,ca3v,si,csi,
                                  "NonSquareM/SymBand");
            TestMatrixArith456<T>(a4x,ca4x,a4v,ca4v,si,csi,
                                  "DegenerateM/SymBand");
        }
    }
    for(size_t i=0;i<B.size();++i) delete B[i];
    for(size_t i=0;i<CB.size();++i) delete CB[i];
#endif
}

#ifdef INST_DOUBLE
template void TestSymBandMatrixArith_B2<double>();
#endif
#ifdef INST_FLOAT
template void TestSymBandMatrixArith_B2<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSymBandMatrixArith_B2<long double>();
#endif
#ifdef INST_INT
template void TestSymBandMatrixArith_B2<int>();
#endif
