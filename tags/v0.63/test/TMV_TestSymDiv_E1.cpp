
#define START 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestSymArith.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestSymDiv_E1(tmv::DivType dt, PosDefCode pdc)
{
    const int N = 10;

    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    std::vector<tmv::BaseMatrix<T>*> B;
    std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
    MakeSymList(s,cs,B,CB,pdc);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(1-3*i+j);
    a1.diag().addToAll(T(10)*N);
    a1 /= T(10);
    tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

    tmv::BandMatrix<T> b1(a1,1,3);
    tmv::BandMatrix<std::complex<T> > cb1(ca1,1,3);

    tmv::BandMatrix<T> b1x = b1;
    tmv::BandMatrix<std::complex<T> > cb1x = cb1;

#ifdef XTEST
    tmv::BandMatrix<T> b3(a1.colRange(0,N-2),1,3);
    tmv::BandMatrix<std::complex<T> > cb3(ca1.colRange(0,N-2),1,3);
    tmv::BandMatrix<T> b4(a1.rowRange(0,N-2),1,3);
    tmv::BandMatrix<std::complex<T> > cb4(ca1.rowRange(0,N-2),1,3);

    tmv::BandMatrix<T> b3x = b3;
    tmv::BandMatrix<T> b4x = b4;  
    tmv::BandMatrix<std::complex<T> > cb3x = cb3;
    tmv::BandMatrix<std::complex<T> > cb4x = cb4;  
#endif

    for(size_t i=START;i<s.size();i++) {
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TMV_Text(s[i])<<
                "  "<<s[i]<<std::endl;
        const tmv::SymMatrixView<T>& si = s[i];
        const tmv::SymMatrixView<std::complex<T> >& csi = cs[i];
        if (dt == tmv::CH && csi.issym()) continue;

        si.saveDiv();
        csi.saveDiv();

        TestMatrixDivArith1<T>(dt,b1x,cb1x,si,b1.view(),csi,cb1.view(),
                               "SquareBandMatrix/Sym");
#ifdef XTEST
        TestMatrixDivArith1<T>(dt,b3x,cb3x,si,b3.view(),csi,cb3.view(),
                               "NonSquareBandMatrix/Sym");
        TestMatrixDivArith1<T>(dt,b4x,cb4x,si,b4.view(),csi,cb4.view(),
                               "NonSquareBandMatrix/Sym");
#endif
    }
    for(size_t i=0;i<B.size();++i) delete B[i];
    for(size_t i=0;i<CB.size();++i) delete CB[i];
}

#ifdef INST_DOUBLE
template void TestSymDiv_E1<double>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_FLOAT
template void TestSymDiv_E1<float>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymDiv_E1<long double>(tmv::DivType dt, PosDefCode pc);
#endif
