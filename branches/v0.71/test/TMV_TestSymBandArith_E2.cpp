#define STARTI 0
#define STARTJ 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymBandArith.h"
#include "TMV_TestBandArith.h"

#define NOELEMMULT

template <class T1, class T2> 
inline bool CanAddEq(
    const tmv::BandMatrixView<T1>& a, const tmv::SymBandMatrixView<T2>& b)
{ 
    return a.colsize() == b.size() && a.rowsize() == b.size() && 
        a.nlo() >= b.nlo() && a.nhi() >= b.nhi(); 
}

template <class T1, class T2, class T3> 
inline bool CanMultMM(
    const tmv::BandMatrixView<T1>& a, const tmv::SymBandMatrixView<T2>& b,
    const tmv::BandMatrixView<T3>& c)
{ 
    return a.rowsize() == b.colsize() && a.colsize() == c.colsize() &&
        b.rowsize() == c.rowsize() &&
        (c.nlo() >= a.nlo() + b.nlo() || c.nlo() == int(c.colsize())-1) &&
        (c.nhi() >= a.nhi() + b.nhi() || c.nhi() == int(c.rowsize())-1);
}

template <class T1, class T2, class T3> 
inline bool CanMultMM(
    const tmv::SymBandMatrixView<T1>& a, const tmv::BandMatrixView<T2>& b,
    const tmv::BandMatrixView<T3>& c)
{ 
    return a.rowsize() == b.colsize() && a.colsize() == c.colsize() &&
        b.rowsize() == c.rowsize() &&
        (c.nlo() >= a.nlo() + b.nlo() || c.nlo() == int(c.colsize())-1) &&
        (c.nhi() >= a.nhi() + b.nhi() || c.nhi() == int(c.rowsize())-1);
}

template <class T1, class T2, class T3> 
inline bool CanMultMM(
    const tmv::SymBandMatrixView<T1>& a, 
    const tmv::SymBandMatrixView<T2>& b, const tmv::BandMatrixView<T3>& c)
{ 
    return a.rowsize() == b.colsize() && a.colsize() == c.colsize() &&
        b.rowsize() == c.rowsize() &&
        (c.nlo() >= a.nlo() + b.nlo() || c.nlo() == int(c.colsize())-1) &&
        (c.nhi() >= a.nhi() + b.nhi() || c.nhi() == int(c.rowsize())-1);
}

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymBandMatrixArith_E2()
{
#if (XTEST & 2)
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

            TestMatrixArith4(bj,cbj,si,csi,"Band/SymBand");
            TestMatrixArith5(bj,cbj,si,csi,"Band/SymBand");
            TestMatrixArith6x(bj,cbj,si,csi,"Band/SymBand");
        }
    }
#endif
}

#ifdef TEST_DOUBLE
template void TestSymBandMatrixArith_E2<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymBandMatrixArith_E2<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandMatrixArith_E2<long double>();
#endif
#ifdef TEST_INT
template void TestSymBandMatrixArith_E2<int>();
#endif
