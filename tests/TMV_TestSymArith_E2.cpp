#define STARTI 0
#define STARTJ 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymArith.h"
#include "TMV_TestBandArith.h"

#define NOADDEQ
#define NOMULTEQ
#define NOELEMMULT

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymMatrixArith_E2()
{
#if (XTEST & 2)
    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    MakeSymList(s,cs,InDef);

    std::vector<tmv::BandMatrixView<T> > b;
    std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
    MakeBandList(b,cb);

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

            TestMatrixArith4(bj,cbj,si,csi,"Band/Sym");
            TestMatrixArith5(bj,cbj,si,csi,"Band/Sym");
            TestMatrixArith6x(bj,cbj,si,csi,"Band/Sym");
        }
    }
#endif
}

#ifdef TEST_DOUBLE
template void TestSymMatrixArith_E2<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymMatrixArith_E2<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymMatrixArith_E2<long double>();
#endif
#ifdef TEST_INT
template void TestSymMatrixArith_E2<int>();
#endif
