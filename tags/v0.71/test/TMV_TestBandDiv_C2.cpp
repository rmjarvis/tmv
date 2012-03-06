
#define START 0

#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestBandArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestBandDiv_C2(tmv::DivType dt)
{
    const int N = 10;

    std::vector<tmv::BandMatrixView<T> > b;
    std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
    MakeBandList(b,cb);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(1-3*i+j);
    a1.diag().addToAll(T(10)*N);
    a1 /= T(10);
    tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

    tmv::DiagMatrix<T> d(a1);
    tmv::DiagMatrix<std::complex<T> > cd(ca1);
    tmv::DiagMatrixView<T> dv = d.view();
    tmv::DiagMatrixView<std::complex<T> > cdv = cd.view();

    for(size_t i=START;i<b.size();i++) {
        if (showstartdone) 
            std::cout<<"Start C2 loop: i = "<<i<<"\nbi = "<<tmv::TMV_Text(b[i])<<
                "  "<<b[i]<<std::endl;
        tmv::BandMatrixView<T> bi = b[i];
        tmv::BandMatrixView<std::complex<T> > cbi = cb[i];

        TestMatrixDivArith1(dt,dv,bi,cdv,cbi,"Band/DiagMatrix");
    }
}

#ifdef TEST_DOUBLE
template void TestBandDiv_C2<double>(tmv::DivType dt);
#endif
#ifdef TEST_FLOAT
template void TestBandDiv_C2<float>(tmv::DivType dt);
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandDiv_C2<long double>(tmv::DivType dt);
#endif
