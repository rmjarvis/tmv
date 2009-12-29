
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
void TestSymBandDiv_C2(tmv::DivType dt, PosDefCode pdc)
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

    tmv::DiagMatrix<T> d(a1);
    tmv::DiagMatrix<std::complex<T> > cd(ca1);

    for(size_t i=START;i<sb.size();i++) {
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TMV_Text(sb[i])<<
                "  "<<sb[i]<<std::endl;
        const tmv::SymBandMatrixView<T>& si = sb[i];
        const tmv::SymBandMatrixView<std::complex<T> >& csi = csb[i];

        if (csi.issym()) {
            tmv::SymBandMatrix<T> sx = si;
            tmv::SymBandMatrix<std::complex<T> > csx = csi;

            TestMatrixDivArith1<T>(dt,sx,csx,d.View(),si,cd.View(),csi,
                                   "SymBand/DiagMatrix");
        } else {
            tmv::HermBandMatrix<T> hx = si;
            tmv::HermBandMatrix<std::complex<T> > chx = csi;

            TestMatrixDivArith1<T>(dt,hx,chx,d.View(),si,cd.View(),csi,
                                   "HermBand/DiagMatrix");
        }
    }
    for(size_t i=0;i<B.size();++i) delete B[i];
    for(size_t i=0;i<CB.size();++i) delete CB[i];
}

#ifdef INST_DOUBLE
template void TestSymBandDiv_C2<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_FLOAT
template void TestSymBandDiv_C2<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymBandDiv_C2<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
