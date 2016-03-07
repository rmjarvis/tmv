
#define START 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymArith.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestSymDiv_C1(tmv::DivType dt, PosDefCode pdc)
{
    const int N = 10;

    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    MakeSymList(s,cs,pdc);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(1-3*i+j);
    a1.diag().addToAll(T(10)*N);
    a1 /= T(10);
    tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

    tmv::DiagMatrix<T> d(a1);
    tmv::DiagMatrix<std::complex<T> > cd(ca1);
    tmv::DiagMatrixView<T> dv = d.view();
    tmv::DiagMatrixView<std::complex<T> > cdv = cd.view();

    for(size_t i=START;i<s.size();i++) {
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TMV_Text(s[i])<<
                "  "<<s[i]<<std::endl;
        tmv::SymMatrixView<T> si = s[i];
        tmv::SymMatrixView<std::complex<T> > csi = cs[i];
        if (dt == tmv::CH && csi.issym()) continue;

        si.saveDiv();
        csi.saveDiv();

        TestMatrixDivArith1(dt,si,dv,csi,cdv,"DiagMatrix/Sym");
    }
}

#ifdef TEST_DOUBLE
template void TestSymDiv_C1<double>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef TEST_FLOAT
template void TestSymDiv_C1<float>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymDiv_C1<long double>(tmv::DivType dt, PosDefCode pc);
#endif
