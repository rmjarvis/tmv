
#define START 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestSymArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestSymDiv_B1(tmv::DivType dt, PosDefCode pdc)
{
    const int N = 10;

    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    std::vector<tmv::BaseMatrix<T>*> B;
    std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
    MakeSymList(s,cs,B,CB,pdc);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(1-3*i+j);
    a1.diag().addToAll(T(10)*N);
    a1 /= T(10);
    tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

    tmv::Matrix<T> a1x = a1;
    tmv::Matrix<std::complex<T> > ca1x = ca1;
#ifdef XTEST
    tmv::Matrix<T> a3 = a1.colRange(0,N/2);
    tmv::Matrix<std::complex<T> > ca3 = ca1.colRange(0,N/2);
    tmv::Matrix<T> a4 = a1.rowRange(0,N/2);
    tmv::Matrix<std::complex<T> > ca4 = ca1.rowRange(0,N/2);
    tmv::Matrix<T> a5 = a1.colRange(0,0);
    tmv::Matrix<std::complex<T> > ca5 = ca1.colRange(0,0);
    tmv::Matrix<T> a6 = a1.rowRange(0,0);
    tmv::Matrix<std::complex<T> > ca6 = ca1.rowRange(0,0);

    tmv::Matrix<T> a3x = a3;
    tmv::Matrix<T> a4x = a4;
    tmv::Matrix<T> a5x = a5;
    tmv::Matrix<T> a6x = a6;
    tmv::Matrix<std::complex<T> > ca3x = ca3;
    tmv::Matrix<std::complex<T> > ca4x = ca4;
    tmv::Matrix<std::complex<T> > ca5x = ca5;
    tmv::Matrix<std::complex<T> > ca6x = ca6;
#endif

    for(size_t i=START;i<s.size();i++) {
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<
                tmv::TMV_Text(s[i])<<"  "<<s[i]<<std::endl;
        const tmv::SymMatrixView<T>& si = s[i];
        const tmv::SymMatrixView<std::complex<T> >& csi = cs[i];
        if (dt == tmv::CH && csi.issym()) continue;

        si.saveDiv();
        csi.saveDiv();

        TestMatrixDivArith1<T>(dt,a1x,ca1x,si,a1.view(),csi,ca1.view(),
                               "SquareMatrix/Sym");
#ifdef XTEST
        TestMatrixDivArith1<T>(dt,a3x,ca3x,si,a3.view(),csi,ca3.view(),
                               "NonSquareMatrix/Sym");
        TestMatrixDivArith1<T>(dt,a4x,ca4x,si,a4.view(),csi,ca4.view(),
                               "NonSquareMatrix/Sym");
        TestMatrixDivArith1<T>(dt,a5x,ca5x,si,a5.view(),csi,ca5.view(),
                               "DegenerateMatrix/Sym");
        TestMatrixDivArith1<T>(dt,a6x,ca6x,si,a6.view(),csi,ca6.view(),
                               "DegenerateMatrix/Sym");
#endif
    }
    for(size_t i=0;i<B.size();++i) delete B[i];
    for(size_t i=0;i<CB.size();++i) delete CB[i];
}

#ifdef INST_DOUBLE
template void TestSymDiv_B1<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_FLOAT
template void TestSymDiv_B1<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymDiv_B1<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
