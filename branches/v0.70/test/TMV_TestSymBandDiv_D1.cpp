
#define START 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymBandArith.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestSymBandDiv_D1(tmv::DivType dt, PosDefCode pdc)
{
    const int N = 10;

    std::vector<tmv::SymBandMatrixView<T> > sb;
    std::vector<tmv::SymBandMatrixView<std::complex<T> > > csb;
    MakeSymBandList(sb,csb,pdc);

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

    for(size_t i=START;i<sb.size();i++) {
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TMV_Text(sb[i])<<
                "  "<<sb[i]<<std::endl;
        tmv::SymBandMatrixView<T> si = sb[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];
        if (dt == tmv::CH && csi.issym()) continue;
        si.saveDiv();
        csi.saveDiv();

        TestMatrixDivArith1(dt,si,uv,csi,cuv,"UpperTriMatrix/SymBand");
        TestMatrixDivArith1(dt,si,lv,csi,clv,"LowerTriMatrix/SymBand");
    }
}

#ifdef TEST_DOUBLE
template void TestSymBandDiv_D1<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef TEST_FLOAT
template void TestSymBandDiv_D1<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymBandDiv_D1<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
