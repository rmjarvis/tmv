#define STARTI 0
#define STARTJ 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestSymArith.h"
#include "TMV_TestBandArith.h"

#define NOADDEQ
#define NOMULTEQ

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymMatrixArith_E2()
{
#ifdef XTEST
    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    std::vector<tmv::BaseMatrix<T>*> B;
    std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
    MakeSymList(s,cs,B,CB,InDef);

    std::vector<tmv::BandMatrixView<T> > b;
    std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
    MakeBandList(b,cb,B,CB);

    for(size_t i=STARTI;i<s.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<s[i]<<std::endl;
        }

        tmv::SymMatrixView<T> si = s[i];
        tmv::SymMatrixView<std::complex<T> > csi = cs[i];

        for(size_t j=STARTJ;j<b.size();j++) {
            if (showstartdone) {
                std::cerr<<"Start sub-loop "<<j<<std::endl;
                std::cerr<<"bj = "<<b[j]<<std::endl;
            }
            tmv::BandMatrixView<T> bj = b[j];
            tmv::BandMatrixView<std::complex<T> > cbj = cb[j];
            tmv::BandMatrix<T> bx = bj;
            tmv::BandMatrix<std::complex<T> > cbx = cbj;

            if (csi.isherm()) {
                tmv::HermMatrix<T> sx = si;
                tmv::HermMatrix<std::complex<T> > csx = csi;
                TestMatrixArith456<T>(bx,cbx,bj,cbj,si,csi,"Band/Herm");
            } else {
                tmv::SymMatrix<T> sx = si;
                tmv::SymMatrix<std::complex<T> > csx = csi;
                TestMatrixArith456<T>(bx,cbx,bj,cbj,si,csi,"Band/Sym");
            }
        }
    }
    for(size_t i=0;i<B.size();++i) delete B[i];
    for(size_t i=0;i<CB.size();++i) delete CB[i];
#endif
}

#ifdef INST_DOUBLE
template void TestSymMatrixArith_E2<double>();
#endif
#ifdef INST_FLOAT
template void TestSymMatrixArith_E2<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestSymMatrixArith_E2<long double>();
#endif
#ifdef INST_INT
template void TestSymMatrixArith_E2<int>();
#endif
