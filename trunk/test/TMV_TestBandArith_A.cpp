
#define START1 0
#define START2 0

#define NODIV 
#define NOASSIGN
#define PRINTALGO_NormB

#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestBandArith.h"

template <class M1, class M2> 
inline bool CanAddEq(
    const tmv::BaseMatrix_Band_Mutable<M1>& a,
    const tmv::BaseMatrix_Band<M2>& b)
{ 
    return a.colsize() == b.colsize() && a.rowsize() == b.rowsize() &&
        a.nhi() >= b.nhi() && a.nlo() >= b.nlo();
}

template <class M1, class M2, class M3> 
inline bool CanMult(
    const tmv::BaseMatrix_Band<M1>& a, const tmv::BaseMatrix_Band<M2>& b,
    const tmv::BaseMatrix_Band_Mutable<M3>& c)
{ 
    return a.rowsize() == b.colsize() && a.colsize() == c.colsize() &&
        b.rowsize() == c.rowsize() &&
        (c.nlo() >= a.nlo() + b.nlo() || c.nlo() == c.colsize()-1) &&
        (c.nhi() >= a.nhi() + b.nhi() || c.nhi() == c.rowsize()-1);
}

template <class M1, class M2, class M3> 
inline bool CanElemMult(
    const tmv::BaseMatrix_Band<M1>& a, const tmv::BaseMatrix_Band<M2>& b,
    const tmv::BaseMatrix_Band_Mutable<M3>& c)
{ 
    return a.rowsize() == c.rowsize() && a.colsize() == c.colsize() &&
        b.rowsize() == c.rowsize() && b.colsize() == c.colsize() &&
        c.nlo() >= std::min(a.nlo(),b.nlo()) &&
        c.nhi() >= std::min(a.nhi(),b.nhi());
}


#include "TMV_TestMatrixArith.h"

template <class T> 
void TestBandMatrixArith_A()
{
    std::vector<tmv::BandMatrixView<T> > b;
    std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
    MakeBandList(b,cb);

    for(size_t i=START1;i<b.size();i++) {
        if (showstartdone) {
            std::cerr<<"Start loop "<<i<<std::endl;
            std::cerr<<"b[i] = "<<b[i]<<std::endl;
        }
        tmv::BandMatrixView<T> bi = b[i];
        tmv::BandMatrixView<std::complex<T> > cbi = cb[i];
        TestMatrixArith1(bi,cbi,"Band");
        TestMatrixArith2(bi,cbi,"Band");
        TestMatrixArith3(bi,cbi,"Band");

        for(size_t j=START2;j<b.size();j++) if (i!=j) {
            if (showstartdone) {
                std::cerr<<"Start sub-loop "<<j<<std::endl;
                std::cerr<<"bj = "<<b[j]<<std::endl;
            }
            TestMatrixArith4(bi,cbi,b[j],cb[j],"Band/Band");
            TestMatrixArith5(bi,cbi,b[j],cb[j],"Band/Band");
            TestMatrixArith6x(bi,cbi,b[j],cb[j],"Band/Band");
        }
    }
}

#ifdef TEST_DOUBLE
template void TestBandMatrixArith_A<double>();
#endif
#ifdef TEST_FLOAT
template void TestBandMatrixArith_A<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandMatrixArith_A<long double>();
#endif
#ifdef TEST_INT
template void TestBandMatrixArith_A<int>();
#endif
