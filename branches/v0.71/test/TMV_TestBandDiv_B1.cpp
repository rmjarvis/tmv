
#define START 0

#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestBandArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestBandDiv_B1(tmv::DivType dt)
{
    const int N = 10;

    std::vector<tmv::BandMatrixView<T> > b;
    std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
    MakeBandList(b,cb);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-2*j);
    a1.diag().addToAll(T(10)*N);
    a1 /= T(10);
    tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);
    tmv::MatrixView<T> a1v = a1.view();
    tmv::MatrixView<std::complex<T> > ca1v = ca1.view();

#if (XTEST & 2)
    tmv::Matrix<T> a3 = a1.colRange(0,N/2);
    tmv::Matrix<std::complex<T> > ca3 = ca1.colRange(0,N/2);
    tmv::Matrix<T> a4 = a1.rowRange(0,N/2);
    tmv::Matrix<std::complex<T> > ca4 = ca1.rowRange(0,N/2);
    tmv::Matrix<T> a5 = a1.colRange(0,0);
    tmv::Matrix<std::complex<T> > ca5 = ca1.colRange(0,0);
    tmv::Matrix<T> a6 = a1.rowRange(0,0);
    tmv::Matrix<std::complex<T> > ca6 = ca1.rowRange(0,0);

    tmv::MatrixView<T> a3v = a3.view();
    tmv::MatrixView<std::complex<T> > ca3v = ca3.view();
    tmv::MatrixView<T> a4v = a4.view();
    tmv::MatrixView<std::complex<T> > ca4v = ca4.view();
    tmv::MatrixView<T> a5v = a5.view();
    tmv::MatrixView<std::complex<T> > ca5v = ca5.view();
    tmv::MatrixView<T> a6v = a6.view();
    tmv::MatrixView<std::complex<T> > ca6v = ca6.view();

#endif

    for(size_t i=START;i<b.size();i++) {
        if (showstartdone) 
            std::cout<<"Start B1 loop: i = "<<i<<"\nbi = "<<tmv::TMV_Text(b[i])<<
                "  "<<b[i]<<std::endl;
        tmv::BandMatrixView<T> bi = b[i];
        tmv::BandMatrixView<std::complex<T> > cbi = cb[i];
        if (dt == tmv::LU && !bi.isSquare()) continue;

        bi.saveDiv();
        cbi.saveDiv();

        TestMatrixDivArith1(dt,bi,a1v,cbi,ca1v,"SquareMatrix/Band");
#if (XTEST & 2)
        TestMatrixDivArith1(dt,bi,a3v,cbi,ca3v,"NonSquareMatrix/Band");
        TestMatrixDivArith1(dt,bi,a4v,cbi,ca4v,"NonSquareMatrix/Band");
        TestMatrixDivArith1(dt,bi,a5v,cbi,ca5v,"DegenerateMatrix/Band");
        TestMatrixDivArith1(dt,bi,a6v,cbi,ca6v,"DegenerateMatrix/Band");
#endif
    }
}

#ifdef TEST_DOUBLE
template void TestBandDiv_B1<double>(tmv::DivType dt);
#endif
#ifdef TEST_FLOAT
template void TestBandDiv_B1<float>(tmv::DivType dt);
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandDiv_B1<long double>(tmv::DivType dt);
#endif
