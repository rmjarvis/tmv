
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
void TestSymBandDiv_A(tmv::DivType dt, PosDefCode pdc)
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

    tmv::SymBandMatrix<T> s1(a1,2);
    tmv::SymBandMatrix<std::complex<T> > cs1(ca1,2);
    tmv::HermBandMatrix<T> h1(a1,2);
    tmv::HermBandMatrix<std::complex<T> > ch1(ca1,2);

    tmv::SymBandMatrixView<T> s1v = s1.view();
    tmv::SymBandMatrixView<std::complex<T> > cs1v = cs1.view();
    tmv::SymBandMatrixView<T> h1v = h1.view();
    tmv::SymBandMatrixView<std::complex<T> > ch1v = ch1.view();

    for(size_t i=START;i<sb.size();i++) {
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TMV_Text(sb[i])<<
                "  "<<sb[i]<<std::endl;

        tmv::SymBandMatrixView<T> si = sb[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];
        if (dt == tmv::CH && csi.issym()) continue;
        si.saveDiv();
        csi.saveDiv();

        TestMatrixDivArith2(dt,si,s1v,csi,cs1v,"SymBand/SymBand");
        TestMatrixDivArith1(dt,si,h1v,csi,ch1v,"HermBand/SymBand");
    }
}

#ifdef TEST_DOUBLE
template void TestSymBandDiv_A<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef TEST_FLOAT
template void TestSymBandDiv_A<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandDiv_A<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
