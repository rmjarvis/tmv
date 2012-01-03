
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
void TestSymDiv_D1(tmv::DivType dt, PosDefCode pdc)
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

    tmv::UpperTriMatrix<T> u(a1);
    tmv::UpperTriMatrix<std::complex<T> > cu(ca1);
    tmv::LowerTriMatrix<T> l(a1);
    tmv::LowerTriMatrix<std::complex<T> > cl(ca1);

    tmv::UpperTriMatrixView<T> uv = u.view();
    tmv::UpperTriMatrixView<std::complex<T> > cuv = cu.view();
    tmv::LowerTriMatrixView<T> lv = l.view();
    tmv::LowerTriMatrixView<std::complex<T> > clv = cl.view();

    for(size_t i=START;i<s.size();i++) {
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TMV_Text(s[i])<<
                "  "<<s[i]<<std::endl;
        tmv::SymMatrixView<T> si = s[i];
        tmv::SymMatrixView<std::complex<T> > csi = cs[i];
        if (dt == tmv::CH && csi.issym()) continue;

        si.saveDiv();
        csi.saveDiv();

        TestMatrixDivArith1(dt,si,uv,csi,cuv,"UpperTriMatrix/Sym");
        TestMatrixDivArith1(dt,si,lv,csi,clv,"LowerTriMatrix/Sym");
    }
}

#ifdef TEST_DOUBLE
template void TestSymDiv_D1<double>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef TEST_FLOAT
template void TestSymDiv_D1<float>(tmv::DivType dt, PosDefCode pc);
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymDiv_D1<long double>(tmv::DivType dt, PosDefCode pc);
#endif
