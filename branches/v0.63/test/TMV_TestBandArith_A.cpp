
#define START1 0
#define START2 7

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestBandArith.h"

template <class T1, class T2> 
inline bool CanAddEq(
    const tmv::BandMatrixView<T1>& a, const tmv::BandMatrixView<T2>& b)
{ 
    return a.colsize() == b.colsize() && a.rowsize() == b.rowsize() &&
        a.nhi() >= b.nhi() && a.nlo() >= b.nlo();
}

template <class T1, class T2, class T3> 
inline bool CanMultMM(
    const tmv::BandMatrixView<T1>& a, const tmv::BandMatrixView<T2>& b,
    const tmv::BandMatrixView<T3>& c)
{ 
    return a.rowsize() == b.colsize() && a.colsize() == c.colsize() &&
        b.rowsize() == c.rowsize() &&
        (c.nlo() >= a.nlo() + b.nlo() || c.nlo() == int(c.colsize())-1) &&
        (c.nhi() >= a.nhi() + b.nhi() || c.nhi() == int(c.rowsize())-1);
}

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestBandMatrixArith_A()
{
    std::vector<tmv::BandMatrixView<T> > b;
    std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
    std::vector<tmv::BaseMatrix<T>*> B;
    std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
    MakeBandList(b,cb,B,CB);

    for(size_t i=START1;i<b.size();i++) {
        if (showstartdone) {
            std::cerr<<"Start loop "<<i<<std::endl;
            std::cerr<<"b[i] = "<<b[i]<<std::endl;
        }
        tmv::BandMatrixView<T> bi = b[i];
        tmv::BandMatrixView<std::complex<T> > cbi = cb[i];
        tmv::BandMatrix<T> bx = bi;
        tmv::BandMatrix<std::complex<T> > cbx = cbi;
        TestMatrixArith123<T>(bx,cbx,bi,cbi,"Band");

        for(size_t j=START2;j<b.size();j++) if (i!=j) {
            if (showstartdone) {
                std::cerr<<"Start sub-loop "<<j<<std::endl;
                std::cerr<<"bj = "<<b[j]<<std::endl;
            }
            TestMatrixArith456<T>(bx,cbx,bi,cbi,b[j],cb[j],"Band/Band");
        }
    }
    for(size_t i=0;i<B.size();i++) delete B[i];
    for(size_t i=0;i<CB.size();i++) delete CB[i];
}

#ifdef INST_DOUBLE
template void TestBandMatrixArith_A<double>();
#endif
#ifdef INST_FLOAT
template void TestBandMatrixArith_A<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestBandMatrixArith_A<long double>();
#endif
#ifdef INST_INT
template void TestBandMatrixArith_A<int>();
#endif
