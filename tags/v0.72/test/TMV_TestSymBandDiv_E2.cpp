
#define START 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymBandArith.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestSymBandDiv_E2(tmv::DivType dt, PosDefCode pdc)
{
    const int N = 10;

    std::vector<tmv::SymBandMatrixView<T> > sb;
    std::vector<tmv::SymBandMatrixView<std::complex<T> > > csb;
    MakeSymBandList(sb,csb,pdc);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(1-3*i+j);
    a1.diag().addToAll(T(10)*N);
    a1 /= T(10);
    tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

    tmv::BandMatrix<T> b1(a1,1,3);
    tmv::BandMatrix<std::complex<T> > cb1(ca1,1,3);

    tmv::BandMatrixView<T> b1v = b1.view();
    tmv::BandMatrixView<std::complex<T> > cb1v = cb1.view();

#if (XTEST & 2)
    tmv::BandMatrix<T> b3(a1.colRange(0,N-2),1,3);
    tmv::BandMatrix<std::complex<T> > cb3(ca1.colRange(0,N-2),1,3);
    tmv::BandMatrix<T> b4(a1.rowRange(0,N-2),1,3);
    tmv::BandMatrix<std::complex<T> > cb4(ca1.rowRange(0,N-2),1,3);

    tmv::BandMatrixView<T> b3v = b3.view();
    tmv::BandMatrixView<std::complex<T> > cb3v = cb3.view();
    tmv::BandMatrixView<T> b4v = b4.view();
    tmv::BandMatrixView<std::complex<T> > cb4v = cb4.view();
#endif

    for(size_t i=START;i<sb.size();i++) {
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TMV_Text(sb[i])<<
                "  "<<sb[i]<<std::endl;
        tmv::SymBandMatrixView<T> si = sb[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];

        TestMatrixDivArith1(dt,b1v,si,cb1v,csi,"SymBand/SquareBandMatrix");
        if (dt == tmv::LU) continue;
#if (XTEST & 2)
        TestMatrixDivArith1(dt,b3v,si,cb3v,csi,"SymBand/NonSquareBandMatrix");
        TestMatrixDivArith1(dt,b4v,si,cb4v,csi,"SymBand/NonSquareBandMatrix");
#endif
    }
}

#ifdef TEST_DOUBLE
template void TestSymBandDiv_E2<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef TEST_FLOAT
template void TestSymBandDiv_E2<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandDiv_E2<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
