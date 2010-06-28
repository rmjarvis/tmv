
#define START1 0
#define START2 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestSymBandArith.h"

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

#define NOMULTEQ

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymBandMatrixArith_A()
{
    std::vector<tmv::SymBandMatrixView<T> > sb;
    std::vector<tmv::SymBandMatrixView<std::complex<T> > > csb;
    std::vector<tmv::BaseMatrix<T>*> B;
    std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
    MakeSymBandList(sb,csb,B,CB,InDef);

    for(size_t i=START1;i<sb.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<sb[i]<<std::endl;
        }
        tmv::SymBandMatrixView<T> si = sb[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];

        TestMatrixArith1<T>(si,csi,"SymBand");
        TestMatrixArith2<T>(si,csi,"SymBand");
        TestMatrixArith3<T>(si,csi,"SymBand");

        for(size_t j=START2;j<sb.size();j++) if (i!=j) {
            if (showstartdone) {
                std::cout<<"Start sub-loop j = "<<j<<std::endl;
                std::cout<<"sj = "<<sb[j]<<std::endl;
            }
            TestMatrixArith4<T>(si,csi,sb[j],csb[j],"SymBand/SymBand");
            TestMatrixArith5<T>(si,csi,sb[j],csb[j],"SymBand/SymBand");
            TestMatrixArith6x<T>(si,csi,sb[j],csb[j],"SymBand/SymBand");
        }
    }
    for(size_t i=0;i<B.size();++i) delete B[i];
    for(size_t i=0;i<CB.size();++i) delete CB[i];
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
