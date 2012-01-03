#define STARTI 0
#define STARTJ 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymBandArith.h"
#include "TMV_TestBandArith.h"

#define NOADDEQ
#define NOMULTEQ
#define NOELEMMULT

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymBandMatrixArith_E1()
{
    std::vector<tmv::SymBandMatrixView<T> > sb;
    std::vector<tmv::SymBandMatrixView<std::complex<T> > > csb;
    MakeSymBandList(sb,csb,InDef);

    std::vector<tmv::BandMatrixView<T> > b;
    std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
    MakeBandList(b,cb);

    for(size_t i=STARTI;i<sb.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<sb[i]<<std::endl;
        }

        tmv::SymBandMatrixView<T> si = sb[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];

        for(size_t j=STARTJ;j<b.size();j++) {
            if (showstartdone) {
                std::cerr<<"Start sub-loop "<<j<<std::endl;
                std::cerr<<"bj = "<<b[j]<<std::endl;
            }
            tmv::BandMatrixView<T> bj = b[j];
            tmv::BandMatrixView<std::complex<T> > cbj = cb[j];

            TestMatrixArith4(si,csi,bj,cbj,"SymBand/Band");
            TestMatrixArith5(si,csi,bj,cbj,"SymBand/Band");
            TestMatrixArith6x(si,csi,bj,cbj,"SymBand/Band");
        }
    }
}

#ifdef TEST_DOUBLE
template void TestSymBandMatrixArith_E1<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymBandMatrixArith_E1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandMatrixArith_E1<long double>();
#endif
#ifdef TEST_INT
template void TestSymBandMatrixArith_E1<int>();
#endif
