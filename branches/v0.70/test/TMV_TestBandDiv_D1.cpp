
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
void TestBandDiv_D1(tmv::DivType dt)
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

    tmv::UpperTriMatrix<T> u(a1);
    tmv::UpperTriMatrix<std::complex<T> > cu(ca1);
    tmv::LowerTriMatrix<T> l(a1);
    tmv::LowerTriMatrix<std::complex<T> > cl(ca1);
    tmv::UpperTriMatrixView<T> uv = u.view();
    tmv::UpperTriMatrixView<std::complex<T> > cuv = cu.view();
    tmv::LowerTriMatrixView<T> lv = l.view();
    tmv::LowerTriMatrixView<std::complex<T> > clv = cl.view();

    for(size_t i=START;i<b.size();i++) {
        if (showstartdone) 
            std::cout<<"Start D1 loop: i = "<<i<<"\nbi = "<<tmv::TMV_Text(b[i])<<
                "  "<<b[i]<<std::endl;
        tmv::BandMatrixView<T> bi = b[i];
        tmv::BandMatrixView<std::complex<T> > cbi = cb[i];
        if (dt == tmv::LU && !bi.isSquare()) continue;

        bi.saveDiv();
        cbi.saveDiv();

        TestMatrixDivArith1(dt,bi,uv,cbi,cuv,"UpperTriMatrix/Band");
        TestMatrixDivArith1(dt,bi,lv,cbi,clv,"LowerTriMatrix/Band");
    }
}

#ifdef TEST_DOUBLE
template void TestBandDiv_D1<double>(tmv::DivType dt);
#endif
#ifdef TEST_FLOAT
template void TestBandDiv_D1<float>(tmv::DivType dt);
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandDiv_D1<long double>(tmv::DivType dt);
#endif
