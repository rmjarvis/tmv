
#define START1 0
#define START2 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestSymArith.h"

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

#define NOMULTEQ

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymMatrixArith_A()
{
    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    std::vector<tmv::BaseMatrix<T>*> B;
    std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
    MakeSymList(s,cs,B,CB,InDef);

    for(size_t i=START1;i<s.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<s[i]<<std::endl;
        }
        tmv::SymMatrixView<T> si = s[i];
        tmv::SymMatrixView<std::complex<T> > csi = cs[i];

        TestMatrixArith1<T>(si,csi,"Sym");
        TestMatrixArith2<T>(si,csi,"Sym");
        TestMatrixArith3<T>(si,csi,"Sym");

        for(size_t j=START2;j<s.size();j++) if (i!=j) {
            if (showstartdone) {
                std::cout<<"Start sub-loop j = "<<j<<std::endl;
                std::cout<<"sj = "<<s[j]<<std::endl;
            }
            TestMatrixArith4<T>(si,csi,s[j],cs[j],"Sym/Sym");
            TestMatrixArith5<T>(si,csi,s[j],cs[j],"Sym/Sym");
            TestMatrixArith6x<T>(si,csi,s[j],cs[j],"Sym/Sym");
        }
    }
    for(size_t i=0;i<B.size();++i) delete B[i];
    for(size_t i=0;i<CB.size();++i) delete CB[i];
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
