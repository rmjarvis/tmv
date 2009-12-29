
#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_SymBand.h"
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
    std::vector<tmv::BaseMatrix<T>*> B;
    std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
    MakeSymBandList(sb,csb,B,CB,pdc);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(1-3*i+j);
    a1.diag().AddToAll(T(10)*N);
    a1 /= T(10);
    tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

    tmv::SymBandMatrix<T> s1(a1,2);
    tmv::SymBandMatrix<std::complex<T> > cs1(ca1,2);
    tmv::HermBandMatrix<T> h1(a1,2);
    tmv::HermBandMatrix<std::complex<T> > ch1(ca1,2);

    tmv::SymBandMatrix<T> s1x = s1;
    tmv::SymBandMatrix<std::complex<T> > cs1x = s1;
    tmv::HermBandMatrix<T> h1x = h1;
    tmv::HermBandMatrix<std::complex<T> > ch1x = ch1;

    for(size_t i=START;i<sb.size();i++) {
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TMV_Text(sb[i])<<
                "  "<<sb[i]<<std::endl;

        tmv::SymBandMatrixView<T> si = sb[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];
        if (dt == tmv::CH && csi.issym()) continue;
        si.SaveDiv();
        csi.SaveDiv();

        TestMatrixDivArith2<T>(dt,s1x,cs1x,si,s1.View(),csi,cs1.View(),
                               "SymBand/SymBand");
        TestMatrixDivArith1<T>(dt,h1x,ch1x,si,h1.View(),csi,ch1.View(),
                               "HermBand/SymBand");
    }
    for(size_t i=0;i<B.size();++i) delete B[i];
    for(size_t i=0;i<CB.size();++i) delete CB[i];
}

#ifdef INST_DOUBLE
template void TestSymBandDiv_A<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_FLOAT
template void TestSymBandDiv_A<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymBandDiv_A<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
