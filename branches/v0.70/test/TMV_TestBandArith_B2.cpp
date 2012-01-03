#define START 0

#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestBandArith.h"

#define NOELEMMULT

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestBandMatrixArith_B2()
{
#if (XTEST & 2)
    std::vector<tmv::BandMatrixView<T> > b;
    std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
    MakeBandList(b,cb);

    const int N = b[0].rowsize();

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);

    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j)
        ca1(i,j) = std::complex<T>(3+i-5*j,4-8*i-j);

    tmv::MatrixView<T> a1v = a1.view();
    tmv::MatrixView<std::complex<T> > ca1v = ca1.view();

    tmv::Matrix<T> a2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = T(1-3*i+3*j);
    tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j)
        ca2(i,j) = std::complex<T>(1-3*i+6*j,8+2*i-6*j);
    tmv::Matrix<T,tmv::RowMajor> a3 = a2.rowRange(0,N);
    tmv::Matrix<std::complex<T> > ca3 = ca2.rowRange(0,N);
    tmv::Matrix<T,tmv::RowMajor> a4 = a1.colRange(0,0);
    tmv::Matrix<std::complex<T> > ca4 = ca1.colRange(0,0);

    tmv::MatrixView<T> a3v = a3.view();
    tmv::MatrixView<T> a4v = a4.view();
    tmv::MatrixView<std::complex<T> > ca3v = ca3.view();
    tmv::MatrixView<std::complex<T> > ca4v = ca4.view();

    for(size_t i=START;i<b.size();i++) {
        if (showstartdone) {
            std::cerr<<"Start loop "<<i<<std::endl;
            std::cerr<<"bi = "<<b[i]<<std::endl;
        }
        tmv::BandMatrixView<T> bi = b[i];
        tmv::BandMatrixView<std::complex<T> > cbi = cb[i];

        TestMatrixArith4(a1v,ca1v,bi,cbi,"SquareM/Band");
        TestMatrixArith5(a1v,ca1v,bi,cbi,"SquareM/Band");
        TestMatrixArith6x(a1v,ca1v,bi,cbi,"SquareM/Band");
        TestMatrixArith4(a3v,ca3v,bi,cbi,"NonSquareM/Band");
        TestMatrixArith5(a3v,ca3v,bi,cbi,"NonSquareM/Band");
        TestMatrixArith6x(a3v,ca3v,bi,cbi,"NonSquareM/Band");
        TestMatrixArith4(a4v,ca4v,bi,cbi,"DegenerateM/Band");
        TestMatrixArith5(a4v,ca4v,bi,cbi,"DegenerateM/Band");
        TestMatrixArith6x(a4v,ca4v,bi,cbi,"DegenerateM/Band");
    }
#endif
}

#ifdef TEST_DOUBLE
template void TestBandMatrixArith_B2<double>();
#endif
#ifdef TEST_FLOAT
template void TestBandMatrixArith_B2<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandMatrixArith_B2<long double>();
#endif
#ifdef TEST_INT
template void TestBandMatrixArith_B2<int>();
#endif
