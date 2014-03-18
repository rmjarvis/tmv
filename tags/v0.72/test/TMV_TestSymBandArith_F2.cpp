#define STARTI 0
#define STARTJ 0

//#define XXDEBUG8

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymBandArith.h"
#include "TMV_TestSymArith.h"

#define NOMULTEQ
#define NOELEMMULT

template <class T> 
inline bool CanAddEq(
    const tmv::SymMatrixView<std::complex<T> >& a, 
    const tmv::SymBandMatrixView<std::complex<T> >& b)
{ return a.size() == b.size() && a.sym() == b.sym(); }

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymBandMatrixArith_F2()
{
#if (XTEST & 2)
    std::vector<tmv::SymBandMatrixView<T> > sb;
    std::vector<tmv::SymBandMatrixView<std::complex<T> > > csb;
    MakeSymBandList(sb,csb,InDef);

    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    MakeSymList(s,cs,InDef);

    for(size_t i=STARTI;i<sb.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<sb[i]<<std::endl;
        }

        tmv::SymBandMatrixView<T> si = sb[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];

        for(size_t j=STARTJ;j<s.size();j++) {
            if (showstartdone) {
                std::cerr<<"Start sub-loop "<<j<<std::endl;
                std::cerr<<"sj = "<<s[j]<<std::endl;
            }
            tmv::SymMatrixView<T> sj = s[j];
            tmv::SymMatrixView<std::complex<T> > csj = cs[j];

            TestMatrixArith4(sj,csj,si,csi,"Sym/SymBand");
            TestMatrixArith5(sj,csj,si,csi,"Sym/SymBand");
            TestMatrixArith6x(sj,csj,si,csi,"Sym/SymBand");
        }
    }
#endif
}

#ifdef TEST_DOUBLE
template void TestSymBandMatrixArith_F2<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymBandMatrixArith_F2<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandMatrixArith_F2<long double>();
#endif
#ifdef TEST_INT
template void TestSymBandMatrixArith_F2<int>();
#endif
