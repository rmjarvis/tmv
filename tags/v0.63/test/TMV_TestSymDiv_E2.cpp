
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
void TestSymDiv_E2(tmv::DivType dt, PosDefCode pdc)
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

    tmv::BandMatrix<T> b1v = b1.view();
    tmv::BandMatrix<std::complex<T> > cb1v = cb1.view();

#ifdef XTEST
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
        const tmv::SymMatrixView<T>& si = s[i];
        const tmv::SymMatrixView<std::complex<T> >& csi = cs[i];

        if (csi.issym()) {
            tmv::SymMatrix<T> sx = si;
            tmv::SymMatrix<std::complex<T> > csx = csi;

            TestMatrixDivArith1<T>(dt,sx,csx,b1v,si,cb1v,csi,
                                   "Sym/SquareBandMatrix");
            if (dt == tmv::LU) continue;
#ifdef XTEST
            TestMatrixDivArith1<T>(dt,sx,csx,b3v,si,cb3v,csi,
                                   "Sym/NonSquareBandMatrix");
            TestMatrixDivArith1<T>(dt,sx,csx,b4v,si,cb4v,csi,
                                   "Sym/NonSquareBandMatrix");
#endif
        } else {
            tmv::HermMatrix<T> hx = si;
            tmv::HermMatrix<std::complex<T> > chx = csi;

            TestMatrixDivArith1<T>(dt,hx,chx,b1v,si,cb1v,csi,
                                   "Herm/SquareBandMatrix");
            if (dt == tmv::LU) continue;
#ifdef XTEST
            TestMatrixDivArith1<T>(dt,hx,chx,b3v,si,cb3v,csi,
                                   "Herm/NonSquareBandMatrix");
            TestMatrixDivArith1<T>(dt,hx,chx,b4v,si,cb4v,csi,
                                   "Herm/NonSquareBandMatrix");
#endif
        }
    }
    for(size_t i=0;i<B.size();++i) delete B[i];
    for(size_t i=0;i<CB.size();++i) delete CB[i];
}

#ifdef INST_DOUBLE
template void TestSymDiv_E2<double>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_FLOAT
template void TestSymDiv_E2<float>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymDiv_E2<long double>(tmv::DivType dt, PosDefCode pc);
#endif
