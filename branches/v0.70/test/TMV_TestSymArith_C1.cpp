#define START 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymArith.h"

#define NOMULTEQ
#define NOELEMMULT

template <class T1, class T2> 
inline bool CanAddEq(
    const tmv::SymMatrixView<T1>& a, const tmv::DiagMatrixView<T2>& )
{ return a.issym() || tmv::isReal(T2()); }

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymMatrixArith_C1()
{
    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    MakeSymList(s,cs,InDef);

    const int N = s[0].size();

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) = 
        std::complex<T>(3+i-5*j,2-3*i);

    tmv::DiagMatrix<T> d1(a1);
    tmv::DiagMatrix<std::complex<T> > cd1(ca1);
    tmv::DiagMatrixView<T> d1v = d1.view();
    tmv::DiagMatrixView<std::complex<T> > cd1v = cd1.view();

    for(size_t i=START;i<s.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<s[i]<<std::endl;
        }

        tmv::SymMatrixView<T> si = s[i];
        tmv::SymMatrixView<std::complex<T> > csi = cs[i];

        TestMatrixArith4(si,csi,d1v,cd1v,"Sym/Diag");
        TestMatrixArith5(si,csi,d1v,cd1v,"Sym/Diag");
        TestMatrixArith6x(si,csi,d1v,cd1v,"Sym/Diag");
    }
}

#ifdef TEST_DOUBLE
template void TestSymMatrixArith_C1<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymMatrixArith_C1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymMatrixArith_C1<long double>();
#endif
#ifdef TEST_INT
template void TestSymMatrixArith_C1<int>();
#endif
