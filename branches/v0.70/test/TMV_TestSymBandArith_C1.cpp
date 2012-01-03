#define START 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymBandArith.h"

#define NOMULTEQ
#define NOELEMMULT

template <class T1, class T2> 
inline bool CanAddEq(
    const tmv::SymBandMatrixView<T1>& a, const tmv::DiagMatrixView<T2>& )
{ return a.issym() || tmv::isReal(T2()); }

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymBandMatrixArith_C1()
{
    std::vector<tmv::SymBandMatrixView<T> > sb;
    std::vector<tmv::SymBandMatrixView<std::complex<T> > > csb;
    MakeSymBandList(sb,csb,InDef);

    const int N = sb[0].size();

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) =
        std::complex<T>(3+i-5*j,2-3*i);

    tmv::DiagMatrix<T> d1(a1);
    tmv::DiagMatrix<std::complex<T> > cd1(ca1);
    tmv::DiagMatrixView<T> d1v = d1.view();
    tmv::DiagMatrixView<std::complex<T> > cd1v = cd1.view();

    for(size_t i=START;i<sb.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<sb[i]<<std::endl;
        }

        tmv::SymBandMatrixView<T> si = sb[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];

        TestMatrixArith4(si,csi,d1v,cd1v,"SymBand/Diag");
        TestMatrixArith5(si,csi,d1v,cd1v,"SymBand/Diag");
        TestMatrixArith6x(si,csi,d1v,cd1v,"SymBand/Diag");
    }
}

#ifdef TEST_DOUBLE
template void TestSymBandMatrixArith_C1<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymBandMatrixArith_C1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandMatrixArith_C1<long double>();
#endif
#ifdef TEST_INT
template void TestSymBandMatrixArith_C1<int>();
#endif
