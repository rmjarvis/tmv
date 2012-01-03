
#define START1 0
#define START2 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymArith.h"

#define NOMULTEQ

template <class T> 
inline bool CanAddEq(
    const tmv::SymMatrixView<std::complex<T> >& a,
    const tmv::SymMatrixView<std::complex<T> >& b)
{ return a.size() == b.size() && a.sym() == b.sym(); }

template <class T1, class T2> 
inline bool CanAddEqX(
    const tmv::SymMatrixView<T1>& a, const T2 x)
{ return tmv::isReal(x) || a.issym(); }

template <class T1, class T2> 
inline bool CanMultEqX(
    const tmv::SymMatrixView<T1>& a, const T2 x)
{ return tmv::isReal(x) || a.issym(); }

template <class T1, class T2>
static inline bool CanElemMultMM(
    const tmv::SymMatrixView<T1>& a, const tmv::SymMatrixView<T2>& b)
{
    return a.size() == b.size() &&
        a.issym() == b.issym() && a.isherm() == b.isherm();
}

template <class T1, class T2, class T3>
static inline bool CanAddElemMultMM(
    const tmv::SymMatrixView<T1>& a, const tmv::SymMatrixView<T2>& b,
    const tmv::SymMatrixView<T3>& c)
{
    return a.size() == b.size() && a.size() == c.size() &&
        a.issym() == b.issym() && a.issym() == c.issym() &&
        a.isherm() == b.isherm() && a.isherm() == c.isherm();
}

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymMatrixArith_A()
{
    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    MakeSymList(s,cs,InDef);
 
    for(size_t i=START1;i<s.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<s[i]<<std::endl;
        }
        tmv::SymMatrixView<T> si = s[i];
        tmv::SymMatrixView<std::complex<T> > csi = cs[i];

        TestMatrixArith1(si,csi,"Sym");
        TestMatrixArith2(si,csi,"Sym");
        TestMatrixArith3(si,csi,"Sym");

        for(size_t j=START2;j<s.size();j++) if (i!=j) {
            if (showstartdone) {
                std::cout<<"Start sub-loop j = "<<j<<std::endl;
                std::cout<<"sj = "<<s[j]<<std::endl;
            }
            TestMatrixArith4(si,csi,s[j],cs[j],"Sym/Sym");
            TestMatrixArith5(si,csi,s[j],cs[j],"Sym/Sym");
            TestMatrixArith6x(si,csi,s[j],cs[j],"Sym/Sym");
        }
    }
}

#ifdef TEST_DOUBLE
template void TestSymMatrixArith_A<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymMatrixArith_A<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymMatrixArith_A<long double>();
#endif
#ifdef TEST_INT
template void TestSymMatrixArith_A<int>();
#endif
