
#define START 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymArith.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestSymDiv_E2(tmv::DivType dt, PosDefCode pdc)
{
    const int N = 10;

    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    MakeSymList(s,cs,pdc);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(1-3*i+j);
    a1.diag().addToAll(T(10)*N);
    a1 /= T(10);
    tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

    tmv::BandMatrix<T> b1(a1,1,3);
    tmv::BandMatrix<std::complex<T> > cb1(ca1,1,3);

    tmv::BandMatrix<T> b1v = b1.view();
    tmv::BandMatrix<std::complex<T> > cb1v = cb1.view();

#if (XTEST & 2)
    tmv::BandMatrix<T> b3(a1.colRange(0,N-2),1,3);
    tmv::BandMatrix<std::complex<T> > cb3(ca1.colRange(0,N-2),1,3);
    tmv::BandMatrix<T> b4(a1.rowRange(0,N-2),1,3);
    tmv::BandMatrix<std::complex<T> > cb4(ca1.rowRange(0,N-2),1,3);

    tmv::BandMatrix<T> b3v = b3.view();
    tmv::BandMatrix<std::complex<T> > cb3v = cb3.view();
    tmv::BandMatrix<T> b4v = b4.view();
    tmv::BandMatrix<std::complex<T> > cb4v = cb4.view();
#endif

    for(size_t i=START;i<s.size();i++) {
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TMV_Text(s[i])<<
                "  "<<s[i]<<std::endl;
        tmv::SymMatrixView<T> si = s[i];
        tmv::SymMatrixView<std::complex<T> > csi = cs[i];

        TestMatrixDivArith1(dt,b1v,si,cb1v,csi,"Sym/SquareBandMatrix");
        if (dt == tmv::LU) continue;
#if (XTEST & 2)
        TestMatrixDivArith1(dt,b3v,si,cb3v,csi,"Sym/NonSquareBandMatrix");
        TestMatrixDivArith1(dt,b4v,si,cb4v,csi,"Sym/NonSquareBandMatrix");
#endif
    }
}

#ifdef TEST_DOUBLE
template void TestSymDiv_E2<double>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef TEST_FLOAT
template void TestSymDiv_E2<float>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymDiv_E2<long double>(tmv::DivType dt, PosDefCode pc);
#endif
