
#define START 0

#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestBandArith.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestBandDiv_A(tmv::DivType dt)
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

    tmv::BandMatrix<T> b1(a1,3,1);
    tmv::BandMatrix<std::complex<T> > cb1 = b1 * std::complex<T>(2,-3);
    tmv::BandMatrix<T> b2 = b1.subBandMatrix(0,N/2,0,N,b1.nlo(),b1.nhi());
    tmv::BandMatrix<std::complex<T> > cb2 = cb1.subBandMatrix(
        0,N/2,0,N,cb1.nlo(),cb1.nhi());
    tmv::BandMatrix<T> b3 = b1.subBandMatrix(0,N,0,N/2,b1.nlo(),b1.nhi());
    tmv::BandMatrix<std::complex<T> > cb3 = cb1.subBandMatrix(
        0,N,0,N/2,cb1.nlo(),cb1.nhi());

    tmv::BandMatrixView<T> b1v = b1.view();
    tmv::BandMatrixView<std::complex<T> > cb1v = cb1.view();
#if (XTEST & 2)
    tmv::BandMatrixView<T> b2v = b2.view();
    tmv::BandMatrixView<std::complex<T> > cb2v = cb2.view();
    tmv::BandMatrixView<T> b3v = b3.view();
    tmv::BandMatrixView<std::complex<T> > cb3v = cb3.view();
#endif

    for(size_t i=START;i<b.size();i++) {
        if (showstartdone) 
            std::cout<<"Start A loop: i = "<<i<<"\nbi = "<<tmv::TMV_Text(b[i])<<
                "  "<<b[i]<<std::endl;
        tmv::BandMatrixView<T> bi = b[i];
        tmv::BandMatrixView<std::complex<T> > cbi = cb[i];
        if (dt == tmv::LU && !bi.isSquare()) continue;

        bi.saveDiv();
        cbi.saveDiv();

        TestMatrixDivArith2(dt,bi,b1v,cbi,cb1v,"SquareBand/Band");
#if (XTEST & 2)
        TestMatrixDivArith1(dt,bi,b2v,cbi,cb2v,"NonSquareBand/Band");
        TestMatrixDivArith1(dt,bi,b3v,cbi,cb3v,"NonSquareBand/Band");
#endif
    }
}

#ifdef TEST_DOUBLE
template void TestBandDiv_A<double>(tmv::DivType dt);
#endif
#ifdef TEST_FLOAT
template void TestBandDiv_A<float>(tmv::DivType dt);
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandDiv_A<long double>(tmv::DivType dt);
#endif
