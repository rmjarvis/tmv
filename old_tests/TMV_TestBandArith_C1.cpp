#define START 0

#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestBandArith.h"

#define NOELEMMULT

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestBandMatrixArith_C1()
{
    std::vector<tmv::BandMatrixView<T> > b;
    std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
    MakeBandList(b,cb);

    const int N = b[0].rowsize();

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j)
        ca1(i,j) = std::complex<T>(3+i-5*j,4-8*i-j);

    tmv::DiagMatrix<T> d1(a1);
    tmv::DiagMatrix<std::complex<T> > cd1(ca1);
    tmv::DiagMatrixView<T> d1v = d1.view();
    tmv::DiagMatrixView<std::complex<T> > cd1v = cd1.view();
    tmv::DiagMatrix<T> d1x = d1v;
    tmv::DiagMatrix<std::complex<T> > cd1x = cd1v;

    for(size_t i=START;i<b.size();i++) {
        if (showstartdone) {
            std::cerr<<"Start loop "<<i<<std::endl;
            std::cerr<<"bi = "<<b[i]<<std::endl;
        }
        tmv::BandMatrixView<T> bi = b[i];
        tmv::BandMatrixView<std::complex<T> > cbi = cb[i];

        TestMatrixArith4(bi,cbi,d1v,cd1v,"Band/Diag");
        TestMatrixArith5(bi,cbi,d1v,cd1v,"Band/Diag");
        TestMatrixArith6x(bi,cbi,d1v,cd1v,"Band/Diag");
    }
}

#ifdef TEST_DOUBLE
template void TestBandMatrixArith_C1<double>();
#endif
#ifdef TEST_FLOAT
template void TestBandMatrixArith_C1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandMatrixArith_C1<long double>();
#endif
#ifdef TEST_INT
template void TestBandMatrixArith_C1<int>();
#endif
