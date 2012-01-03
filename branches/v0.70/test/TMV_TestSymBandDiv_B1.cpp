
#define START 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymBandArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestSymBandDiv_B1(tmv::DivType dt, PosDefCode pdc)
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

    for(size_t i=START;i<sb.size();i++) {
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TMV_Text(sb[i])<<
                "  "<<sb[i]<<std::endl;

        tmv::SymBandMatrixView<T> si = sb[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];
        if (dt == tmv::CH && csi.issym()) continue;

        si.saveDiv();
        csi.saveDiv();

        TestMatrixDivArith1(dt,si,a1v,csi,ca1v,"SquareMatrix/SymBand");
#if (XTEST & 2)
        TestMatrixDivArith1(dt,si,a3v,csi,ca3v,"NonSquareMatrix/SymBand");
        TestMatrixDivArith1(dt,si,a4v,csi,ca4v,"NonSquareMatrix/SymBand");
        TestMatrixDivArith1(dt,si,a5v,csi,ca5v,"DegenerateMatrix/SymBand");
        TestMatrixDivArith1(dt,si,a6v,csi,ca6v,"DegenerateMatrix/SymBand");
#endif
    }
}

#ifdef TEST_DOUBLE
template void TestSymBandDiv_B1<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef TEST_FLOAT
template void TestSymBandDiv_B1<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandDiv_B1<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
