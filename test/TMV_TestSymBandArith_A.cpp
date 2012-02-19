
#define START1 0
#define START2 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymBandArith.h"

#define NOMULTEQ

template <class T1, class T2> 
inline bool CanAddEq(
    const tmv::SymBandMatrixView<T1>& a, const tmv::SymBandMatrixView<T2>& b)
{ return a.size() == b.size() && a.nlo() >= b.nlo(); }

template <class T> 
inline bool CanAddEq(
    const tmv::SymBandMatrixView<std::complex<T> >& a, 
    const tmv::SymBandMatrixView<std::complex<T> >& b)
{ return a.size() == b.size() && a.nlo() >= b.nlo() && a.sym() == b.sym(); }

template <class T1, class T2> 
inline bool CanAddEqX(
    const tmv::SymBandMatrixView<T1>& a, const T2 x)
{ return tmv::isReal(x) || !a.isherm(); }

template <class T1, class T2> 
inline bool CanMultEqX(
    const tmv::SymBandMatrixView<T1>& a, const T2 x)
{ return tmv::isReal(x) || !a.isherm(); }

template <class T1, class T2>
static inline bool CanElemMultMM(
    const tmv::SymBandMatrixView<T1>& a, const tmv::SymBandMatrixView<T2>& b)
{
    return a.size() == b.size() &&
        a.issym() == b.issym() && a.isherm() == b.isherm();
}

template <class T1, class T2, class T3>
static inline bool CanAddElemMultMM(
    const tmv::SymBandMatrixView<T1>& a, const tmv::SymBandMatrixView<T2>& b,
    const tmv::SymBandMatrixView<T3>& c)
{
    return a.size() == b.size() && a.size() == c.size() &&
        a.issym() == b.issym() && a.issym() == c.issym() &&
        a.isherm() == b.isherm() && a.isherm() == c.isherm() &&
        c.nlo() > std::min(a.nlo(),b.nlo());
}


#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymBandMatrixArith_A()
{
    std::vector<tmv::SymBandMatrixView<T> > sb;
    std::vector<tmv::SymBandMatrixView<std::complex<T> > > csb;
    MakeSymBandList(sb,csb,InDef);

    for(size_t i=START1;i<sb.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<sb[i]<<std::endl;
        }
        tmv::SymBandMatrixView<T> si = sb[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];

        TestMatrixArith1(si,csi,"SymBand");
        TestMatrixArith2(si,csi,"SymBand");
        TestMatrixArith3(si,csi,"SymBand");

        for(size_t j=START2;j<sb.size();j++) if (i!=j) {
            if (showstartdone) {
                std::cout<<"Start sub-loop j = "<<j<<std::endl;
                std::cout<<"sj = "<<sb[j]<<std::endl;
            }
            TestMatrixArith4(si,csi,sb[j],csb[j],"SymBand/SymBand");
            TestMatrixArith5(si,csi,sb[j],csb[j],"SymBand/SymBand");
            TestMatrixArith6x(si,csi,sb[j],csb[j],"SymBand/SymBand");
        }
    }
}

#ifdef TEST_DOUBLE
template void TestSymBandMatrixArith_A<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymBandMatrixArith_A<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandMatrixArith_A<long double>();
#endif
#ifdef TEST_INT
template void TestSymBandMatrixArith_A<int>();
#endif
